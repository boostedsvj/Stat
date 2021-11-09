from __future__ import print_function
import os, sys, traceback
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
from collections import namedtuple
from copy import deepcopy
import uuid, time, inspect, itertools

debugws = False

# representations of info needed to construct RooFit objects on the fly
VarInfo = namedtuple('VarInfo', ['name', 'title', 'val', 'vmin', 'vmax', 'err', 'unit', 'floating'])
PdfInfo = namedtuple('PdfInfo', ['name', 'title', 'formula', 'x', 'pars', 'hist'])

def silence():
    import ROOT as r
    # reduce printouts
    r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.ERROR);

def checkSuff(suff):
    if suff is not None and len(suff)>0 and suff[0] != "_": suff = "_"+suff
    return suff

def makeVar(info, val=None, suff=""):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    suff = checkSuff(suff)
    if info.floating:
        var = r.RooRealVar(info.name+suff, info.title, info.vmin, info.vmax, info.unit)
    else:
        var = r.RooRealVar(info.name+suff, info.title, val if val is not None else info.val, info.vmin, info.vmax)
    return var

def makePdf(info, inits=None, suff=""):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    suff = checkSuff(suff)
    x = makeVar(info.x, suff=suff)
    pars = [makeVar(p, inits[i] if inits is not None else None, suff) for i,p in enumerate(info.pars)]
    allPars = [x]+pars
    pdf_rgp = r.RooGenericPdf(info.name+"_rgp"+suff, info.title, info.formula, r.RooArgList(*allPars))
    pdf = r.RooParametricShapeBinPdf(info.name+suff, info.title, pdf_rgp, x, r.RooArgList(*pars), info.hist)

    # write workspace for debugging
    if debugws:
        wtmp = r.RooWorkspace("w")
        getattr(wtmp,'import')(pdf)
        ofname = "ws{}.root".format(suff)
        if os.path.exists(ofname): os.remove(ofname)
        wtmp.writeToFile(ofname)

    # return all objects to avoid GC
    allObjs = allPars + [pdf_rgp]
    return pdf, allObjs

def makeVarList(varlist):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    iter = varlist.createIterator()
    var = iter.Next()
    newlist = []
    while var:
        newlist.append(var)
        var = iter.Next()

    return newlist

def makeVarInfoList(varlist):
    return [varToInfo(p) for p in makeVarList(varlist)]

def varToInfo(var, floating=False):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    name = var.GetName()
    title = var.GetTitle()
    val = var.getValV()
    vmin = var.getBinning().lowBound()
    vmax = var.getBinning().highBound()
    err = var.getError()
    unit = var.getUnit()
    return VarInfo(name, title, val, vmin, vmax, err, unit, floating)

def pdfToInfo(pdf):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    # access RooGenericPdf formula
    try:
        r.RooGenericPdf2
    except:
        r.gROOT.ProcessLine("""class RooGenericPdf2 : public RooGenericPdf {
public:
TString formExpr() { return _formExpr; }
static TString getFormExpr(RooGenericPdf* pdf){ return ((RooGenericPdf2*)pdf)->formExpr(); }
};""")

    name = pdf.GetName()
    title = pdf.GetTitle()
    formula = str(r.RooGenericPdf2.getFormExpr(pdf.getPdf()))
    x = varToInfo(pdf.getX().arg(), True)
    pars = makeVarInfoList(pdf.getPars())
    bins = pdf.getBins()
    nbins = pdf.getNbins()
    bins.SetSize(nbins+1)
    hist = r.TH1F("hbins","",nbins,bins)

    return PdfInfo(name, title, formula, x, pars, hist)

def remakePdf(pdf, suff="_old"):
    name = pdf.GetName()
    info = pdfToInfo(pdf)
    pdf.SetName(name+suff)
    return makePdf(info)

def remakeData(data, xname, x, suff=""):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    suff = checkSuff(suff)
    hdata = data.createHistogram(xname)
    newdata = r.RooDataHist(data.GetName()+suff,"",r.RooArgList(x),hdata)
    return newdata

# expected args: info (PdfInfo), inits (list of initial values), data (RooAbsData)
def fitOnce(args, tmp=False):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    # suffix for temporary pdfs and vars to keep them separate
    suff = uuid.uuid4().hex if tmp else ""
    pdf, objs = makePdf(args["info"], inits=args["inits"], suff=suff)
    data = remakeData(args["data"], args["info"].x.name, objs[0], suff=suff) if tmp else args["data"]

    ncalls = 100000
    mopt = r.Math.MinimizerOptions()
    mopt.SetMaxFunctionCalls(ncalls)
    mopt.SetMaxIterations(ncalls)

    fitRes = pdf.fitTo(data, r.RooFit.Extended(False), r.RooFit.Save(1), r.RooFit.SumW2Error(True), r.RooFit.Strategy(2), r.RooFit.Minimizer("Minuit2"), r.RooFit.PrintLevel(-1), r.RooFit.Range("Full"))
    
    if tmp:
        # recommended way to compute chi2 in RooFit
        frame = objs[0].frame(r.RooFit.Title("frame_{}".format(suff)))
        data.plotOn(frame, r.RooFit.Name(data.GetName()))
        pdf.plotOn(frame, r.RooFit.Name(pdf.GetName()))
        rchi2 = frame.chiSquare(pdf.GetName(), data.GetName())

        # computes reduced chi2, so get ndf by hand
        dhist = frame.findObject(data.GetName(),r.RooHist.Class())
        ndf = 0
        for i in range(dhist.GetN()):
            x = r.Double(0.)
            y = r.Double(0.)
            dhist.GetPoint(i,x,y)
            if y!=0: ndf += 1

        args["ndf"] = ndf
        args["chi2"] = rchi2*ndf
        args["status"] = fitRes.status()
        pinfo = makeVarInfoList(pdf.getPars())
        args["fitpars"] = [x.val for x in pinfo]
        args["fiterrs"] = [x.err for x in pinfo]
        return args
    else:
        return pdf, objs, fitRes

# one-param version for pool map
def fitOnceTmp(args):
    try:
        return fitOnce(args, True)
    except:
        if args["verbosity"]>=1: print("Crashed combination: {}".format(args["inits"]))
        traceback.print_exc()

# recursively make a list of all combinations
# paramlist: list of lists of allowed values
def varyAll(paramlist, pos=0, val=[], tups=None):
    if tups is None: tups = set()
    vals = paramlist[pos]
    for v in vals:
        tmp = val[:]+[v]
        # check if last param
        if pos+1==len(paramlist):
            tups.add(tuple(tmp))
        else:
            varyAll(paramlist, pos=pos+1, val=tmp, tups=tups)
    if pos==0: return tups

def bruteForce(info, data, initvals, npool, pmax, verbosity=1):
    # use allowed initial values for each parameter of pdf (up to max)
    npars = len(info.pars)
    paramlist = [initvals for p in range(min(npars, 1e10 if pmax is None else pmax[0]))]
    allInits = list(sorted(varyAll(paramlist)))
    if pmax is not None and npars>pmax[0]:
        extras = tuple([pmax[1]]*(npars-pmax[0]))
        allInits = [init+extras for init in allInits]

    # make list of arg combinations for fitOnce
    allArgs = [{"info": deepcopy(info), "inits": inits, "data": data, "verbosity": verbosity} for inits in allInits]

    tstart = time.time()
    resultArgs = []
    if npool==0:
        # run in series
        for a in allArgs:
            resultArgs.append(fitOnceTmp(a))
    else:
        # run in parallel
        npool = min(len(allInits), npool)
        p = Pool(npool)
        resultArgs = p.map(fitOnceTmp, allArgs)
        p.close()
        p.join()
    tstop = time.time()

    # handle any failures, sort by chi2
    total = len(resultArgs)
    passedArgs = sorted([x for x in resultArgs if x is not None and (x["status"]==0 or x["status"]==1)], key = lambda x: x["chi2"])
    passed = len(passedArgs)
    if verbosity>=1: print("bruteForce result: {} out of {} succeeded in {:.2f} sec".format(passed,total,tstop-tstart))

    # debugging info
    if verbosity>=2:
        # determine ndf from non-empty bins
        hist = data.createHistogram("data_hist",makeVar(info.x))
        ndf = sum([hist.GetBinContent(b+1)!=0 and hist.GetBinError(b+1)!=0 for b in range(hist.GetNbinsX())]) - npars
        headers = ["chi2", "chi2/ndf", "inits", "params", "status"]
        nheaders = len(headers)
        headerLengths = [len(x) for x in headers]
        columnFormats = ["{:.2f}", "{:.2f}", "("+', '.join("{:>6.1f}" for i in range(npars))+")", "("+', '.join("{:>8.3f} ({:>8.3f})" for i in range(npars))+")", "{:>3}"]
        def printColumns(args):
            # apply formatting row-wise
            rows = [[
                columnFormats[0].format(arg["chi2"]),
                columnFormats[1].format(arg["chi2"]/float(ndf)),
                columnFormats[2].format(*arg["inits"]),
                # flatten list of tuple pairs
                columnFormats[3].format(*list(itertools.chain(*[(arg["fitpars"][i],arg["fiterrs"][i]) for i in range(npars)]))),
                columnFormats[4].format(arg["status"]),
            ] for arg in args]
            # transpose to find max length for each column
            colLengths = [max(len(row[i]) for row in rows) for i in range(nheaders)] if len(rows)>0 else [0]*nheaders
            colLengths = [max(colLengths[i], headerLengths[i]) for i in range(nheaders)]
            print('  '.join(["{0:<{1}}".format(headers[i], colLengths[i]) for i in range(nheaders)]))
            for row in rows:
                print('  '.join(["{0:>{1}}".format(row[i], colLengths[i]) for i in range(nheaders)]))

        failedArgs = sorted([x for x in resultArgs if x is not None and x["status"]!=0 and x["status"]!=1], key = lambda x: x["chi2"])

        if len(passedArgs)>0:
            print("Passed: (ndf = {})".format(ndf))
            printColumns(passedArgs)
            print("")
        if len(failedArgs)>0:
            print("Failed: (ndf = {})".format(ndf))
            printColumns(failedArgs)
            print("")

    # repeat fit w/ best inits
    pdf, objs, fitRes = fitOnce(passedArgs[0], tmp=False)

    return pdf, objs, fitRes

def main(args):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    # get ROOT objects
    file = r.TFile.Open(args.file)
    ws = file.Get(args.workspace)
    pdf = ws.pdf(args.pdf)

    # convert pdf back to info
    info = pdfToInfo(pdf)

    # convert gen pdf to histogram
    if args.gen:
        gpdf = ws.pdf(args.gen)
        # todo: get this from data?
        nevents = 10000
        hist = gpdf.generateBinned(r.RooArgSet(gpdf.getX().arg()), nevents, True)
    else:
        hist = ws.data(args.data)

    opdf, objs, fitRes = bruteForce(info, hist, args.initvals, args.npool, args.max, args.verbosity)

    # print parameter vals in combine arg format
    opars = makeVarInfoList(opdf.getPars())
    setargs = ["{}={}".format(p.name,p.val) for p in opars]
    print('setParameters {}'.format(','.join(setargs)))

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file", dest="file", type=str, required=True, help="input ROOT file")
    parser.add_argument("-w", "--workspace", dest="workspace", type=str, default="SVJ", help="workspace name")
    parser.add_argument("-p", "--pdf", dest="pdf", type=str, required=True, help="pdf name")
    parser.add_argument("-i", "--initvals", dest="initvals", type=float, default=[-10.0,-1.0,-0.1,0.1,1.0,10.0], nargs='+', help="list of allowed initial values")
    #parser.add_argument("-i", "--initvals", dest="initvals", type=float, default=[-8.0,-0.8,-0.1,0.1,0.8,8.0], nargs='+', help="list of allowed initial values")
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=1, help="number of processes")
    parser.add_argument("-m", "--max", dest="max", type=float, default=None, nargs=2, help="[max # parameters (to generate initvals combinations)] [single initial value (for subsequent parameters)]")
    parser.add_argument("-v", "--verbosity", dest="verbosity", type=int, default=1, help="verbosity level")
    data_group = parser.add_mutually_exclusive_group(required=True)
    data_group.add_argument("-d", "--data", dest="data", type=str, default="", help="dataset name")
    data_group.add_argument("-g", "--gen", dest="gen", type=str, default="", help="pdf name to generate data")
    args = parser.parse_args()

    if args.max is not None: max[0] = int(max[0])
    main(args)

