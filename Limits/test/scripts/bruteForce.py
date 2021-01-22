from __future__ import print_function
import os, sys, traceback
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
from collections import namedtuple
from copy import deepcopy
import uuid

debugws = False

# representations of info needed to construct RooFit objects on the fly
VarInfo = namedtuple('VarInfo', ['name', 'title', 'val', 'vmin', 'vmax', 'unit', 'floating'])
PdfInfo = namedtuple('PdfInfo', ['name', 'title', 'formula', 'x', 'pars', 'hist'])

def silence():
    import ROOT as r
    # reduce printouts
    r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.ERROR);

def checkSuff(suff):
    if len(suff)>0 and suff[0] != "_": suff = "_"+suff
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
    x = makeVar(info.x, suff)
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
    unit = var.getUnit()
    return VarInfo(name, title, val, vmin, vmax, unit, floating)    

def pdfToInfo(pdf):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    # access RooGenericPdf formula
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

# expected args: info (PdfInfo), inits (list of initial values), data (RooAbsData)
def fitOnce(args, tmp=False):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    silence()

    # suffix for temporary pdfs and vars to keep them separate
    suff = uuid.uuid4().hex if tmp else ""
    pdf, objs = makePdf(args["info"], inits=args["inits"], suff=suff)
    data = args["data"]

    ncalls = 100000
    mopt = r.Math.MinimizerOptions()
    mopt.SetMaxFunctionCalls(ncalls)
    mopt.SetMaxIterations(ncalls)

    fitRes = pdf.fitTo(data, r.RooFit.Extended(False), r.RooFit.Save(1), r.RooFit.SumW2Error(True), r.RooFit.Strategy(2), r.RooFit.Minimizer("Minuit2"), r.RooFit.PrintLevel(-1), r.RooFit.Range("Full"))
    
    if tmp:
        args["chi2"] = pdf.createChi2(data).getValV()
        return args
    else:
        return pdf, objs

# one-param version for pool map
def fitOnceTmp(args):
    try:
        return fitOnce(args, True)
    except:
        traceback.print_exc()

# recursively make a list of all combinations
# paramlist: list of lists of allowed values
def varyAll(paramlist, pos=0, val=[], tups=set()):
    vals = paramlist[pos]
    for v in vals:
        tmp = val[:]+[v]
        # check if last param
        if pos+1==len(paramlist):
            tups.add(tuple(tmp))
        else:
            varyAll(paramlist, pos=pos+1, val=tmp, tups=tups)
    if pos==0: return tups

def bruteForce(info, data, initvals, npool):
    # use allowed initial values for each parameter of pdf
    paramlist = [initvals for p in info.pars]
    allInits = list(sorted(varyAll(paramlist)))

    # make list of arg combinations for fitOnce
    allArgs = [{"info": deepcopy(info), "inits": inits, "data": data} for inits in allInits]

    resultArgs = []
    if npool==0:
        # run in series
        for a in allArgs:
            resultArgs.append(fitOnceTmp(a))
    else:
        # run in parallel
        p = Pool(npool)
        resultArgs = p.map(fitOnceTmp, allArgs)
        p.close()
        p.join()

    # sort by chi2
    sortedArgs = sorted(resultArgs, key = lambda x: x["chi2"])

    # repeat fit w/ best inits
    pdf, objs = fitOnce(sortedArgs[0], tmp=False)

    return pdf, objs

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

    opdf, objs = bruteForce(info, hist, args.initvals, args.npool)

    # print parameter vals in combine arg format
    opars = makeVarInfoList(opdf.getPars())
    setargs = ["{}={}".format(p.name,p.val) for p in opars]
    print(','.join(setargs))

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--file", dest="file", type=str, required=True, help="input ROOT file")
    parser.add_argument("-w", "--workspace", dest="workspace", type=str, default="SVJ", help="workspace name")
    parser.add_argument("-p", "--pdf", dest="pdf", type=str, required=True, help="pdf name")
    parser.add_argument("-i", "--initvals", dest="initvals", type=float, default=[-10.0,-1.0,-0.1,0.1,1.0,10.0], nargs='+', help="list of allowed initial values")
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=1, help="number of processes")
    data_group = parser.add_mutually_exclusive_group(required=True)
    data_group.add_argument("-d", "--data", dest="data", type=str, default="", help="dataset name")
    data_group.add_argument("-g", "--gen", dest="gen", type=str, default="", help="pdf name to generate data")
    args = parser.parse_args()
    
    main(args)

