import os,sys,shlex,subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamsTracked, makeSigDict, getSigname, getFname, getWname, getPname, getCombos, getInitFromBF
from collections import OrderedDict
from array import array

def copyVals(tree,index,size):
    tmp = tree.GetVal(index)
    tmp.SetSize(size)
    # deep copy
    return [vv for vv in tmp]

def printCan(can,oname):
    pformats = ["png","pdf"]
    for pformat in pformats:
        if pformat=="pdf":
            pformat = "eps"
        can.Print(oname+"."+pformat,pformat)
        if pformat=="eps":
            os.system("epstopdf {0} && rm {0}".format(oname+".eps"))

def main(sig,name,step,combo,seed,init):
    combos = getCombos()

    # get init vals
    setargs = {}
    errargs = {}
    if init:
        for region in combos[combo]:
            stmp, etmp = getInitFromBF(getWname(sig, region), "SVJ", getPname(region, "Alt" in name), region)
            setargs.update(stmp)
            errargs.update(etmp)

    infname = getFname("Step1"+name, "AsymptoticLimits", combo, sig=sig, seed=seed)
    iparams = getParamsTracked(infname, 0.5)
    ieparams = getParamsTracked(infname, 0.5, includeParam=False, includeErr=True)
    ivals = OrderedDict()
    ierrs = OrderedDict()
    for p in sorted(iparams):
        pname = p.replace("trackedParam_","")
        if not any(x in p for x in combos[combo]): continue
        if init:
            ivals[pname] = setargs[pname]
            ierrs[pname] = errargs[pname]
        else:
            ivals[pname] = iparams[p]
            ierrs[pname] = ieparams[p.replace("trackedParam_","trackedError_")]

    # get likelihood scan vals
    import ROOT as r
    r.gROOT.SetBatch(True)
    mdfname = getFname(step+name, "MultiDimFit", combo, sig=sig, seed=seed)
    mdffile = r.TFile.Open(mdfname)
    try:
        mdftree = mdffile.Get("limit")
    except:
        raise RuntimeError("Could not open {}".format(mdfname))
    allParams = ":".join(["trackedParam_{p}:trackedError_{p}".format(p=p) for p in ivals])
    allQtys = "deltaNLL:r:{}".format(allParams)
    # get bestfit vals first
    mdftree.Draw(allQtys,"quantileExpected==-1","goff")
    brval = mdftree.GetVal(1)[0]
    bvals = OrderedDict()
    berrs = OrderedDict()
    ctr = 2
    for p in ivals:
        bvals[p] = mdftree.GetVal(ctr)[0]
        ctr += 1
        berrs[p] = mdftree.GetVal(ctr)[0]
        ctr += 1

    # get all likelihood vals
    npts = mdftree.Draw(allQtys,"quantileExpected!=-1","goff")
    dnll = copyVals(mdftree,0,npts)
    rval = copyVals(mdftree,1,npts)
    vals = OrderedDict()
    errs = OrderedDict()
    ctr = 2
    for p in ivals:
        vals[p] = copyVals(mdftree,ctr,npts)
        ctr += 1
        errs[p] = copyVals(mdftree,ctr,npts)
        ctr += 1

    # make graphs
    oname = "params_vs_r__{}__{}_{}_{}".format(getSigname(sig),name,step,combo)
    max_params = max(sum(region in k for k in vals) for region in combos[combo])
    pad_size = 500
    can = r.TCanvas(oname,"",max_params*pad_size,len(combos[combo])*pad_size)
    can.Divide(max_params,len(combos[combo]))
    ctr = 1
    graphs = []
    x = array('d', rval)
    for region in combos[combo]:
        for p in sorted(vals):
            if not region in p: continue
            can.cd(ctr)
            y = array('d', vals[p])
            ye = array('d', errs[p])
            # hack to ignore error bars in axis
            gax = r.TGraph(npts,x,y)
            gax.SetMarkerStyle(20)
            gax.GetXaxis().SetTitle("r")
            gax.GetYaxis().SetTitle(p)
            gax.SetTitle("")
            gax.Draw("ap")
            graphs.append(gax)

            gtmp = r.TGraphErrors(npts,x,y,r.nullptr,ye)
            gtmp.SetMarkerStyle(20)
            gtmp.Draw("pz same")
            graphs.append(gtmp)

            # init values
            iy = array('d', [ivals[p]]*npts)
            iye = array('d', [ierrs[p]]*npts)
            igtmp = r.TGraphErrors(npts,x,iy,r.nullptr,iye)
            igtmp.SetLineColor(r.kBlue)
            igtmp.SetLineWidth(2)
            igtmp.SetFillColor(r.kBlue)
            igtmp.SetFillStyle(3444)
            igtmp.Draw("l3 same")
            graphs.append(igtmp)

            # bestfit values
            by = array('d', [bvals[p]]*npts)
            bye = array('d', [berrs[p]]*npts)
            bgtmp = r.TGraphErrors(npts,x,by,r.nullptr,bye)
            bgtmp.SetLineColor(r.kRed)
            bgtmp.SetLineWidth(2)
            bgtmp.SetLineStyle(7)
            bgtmp.SetFillColor(r.kRed)
            bgtmp.SetFillStyle(3444)
            bgtmp.Draw("l3 same")
            graphs.append(bgtmp)

            # draw again on top
            gtmp.Draw("pz same")

            ctr += 1

    printCan(can,oname)

    oname2 = oname.replace("params_","dnll_",1)
    pad_size2 = 700
    can2 = r.TCanvas(oname2,"",pad_size2,pad_size2)
    y = array('d', dnll)
    g = r.TGraph(npts,x,y)
    g.SetMarkerStyle(20)
    g.GetXaxis().SetTitle("r")
    g.GetYaxis().SetTitle("dnll")
    g.SetTitle("")
    g.Draw("ap")
    rline = r.TLine(brval,min(y),brval,max(y))
    rline.SetLineStyle(7)
    rline.SetLineWidth(2)
    rline.SetLineColor(r.kRed)
    rline.Draw("same")
    printCan(can2,oname2)

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-S","--signal", dest="signals", metavar=("mZprime","mDark","rinv","alpha"), type=str, required=True, nargs=4, help="signal parameters")
    parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[step][name])")
    parser.add_argument("-s", "--step", dest="step", type=str, default="Test", help="step name (higgsCombine[step][name])")
    parser.add_argument("-c", "--combo", dest="combo", type=str, required=True, choices=sorted(list(getCombos())), help="combo to plot")
    parser.add_argument("-I", "--init", dest="init", default=False, action='store_true', help="use existing initial values of parameters")
    parser.add_argument("-e", "--seed", dest="seed", type=str, default=None, help="random seed (if toy data)")
    args = parser.parse_args()

    main(makeSigDict(args.sig),args.name,args.step,args.combo,args.seed,args.init)
