import shlex,subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from getParamsTracked import getParamsTracked, getFname
from collections import OrderedDict
from array import array

def copyVals(tree,index,size):
    tmp = tree.GetVal(index)
    tmp.SetSize(size)
    # deep copy
    return [vv for vv in tmp]

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-m", "--mass", dest="mass", type=int, required=True, help="Zprime mass")
parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[step][name])")
parser.add_argument("-s", "--step", dest="step", type=str, default="Test", help="step name (higgsCombine[step][name])")
parser.add_argument("-c", "--combo", dest="combo", type=str, required=True, help="combo to plot")
parser.add_argument("-b", "--batch", dest="batch", default=False, action="store_true", help="batch mode")
args = parser.parse_args()

pformats = ["png"]
combos = {
"cut": ["highCut","lowCut"],
"bdt": ["highSVJ2","lowSVJ2"],
}

# get InitCLs vals
infname = getFname(args.mass, "Step1"+args.name, "AsymptoticLimits", args.combo)
iparams = getParamsTracked(infname, 0.5)
ieparams = getParamsTracked(infname, 0.5, includeParam=False, includeErr=True)
ivals = OrderedDict()
ierrs = OrderedDict()
for p in sorted(iparams):
    pname = p.replace("trackedParam_","")
    if not any(x in p for x in combos[args.combo]): continue
    ivals[pname] = iparams[p]
    ierrs[pname] = ieparams[p.replace("trackedParam_","trackedError_")]

# get likelihood scan vals
import ROOT as r
if args.batch: r.gROOT.SetBatch(True)
mdfname = getFname(args.mass, args.step+args.name, "MultiDimFit", args.combo)
mdffile = r.TFile.Open(mdfname)
mdftree = mdffile.Get("limit")
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
oname = "params_vs_r_{}_{}_{}_{}".format(args.mass,args.name,args.step,args.combo)
max_params = max(sum(region in k for k in vals) for region in combos[args.combo])
pad_size = 500
can = r.TCanvas(oname,"",max_params*pad_size,len(combos[args.combo])*pad_size)
can.Divide(max_params,len(combos[args.combo]))
ctr = 1
graphs = []
x = array('d', rval)
for region in combos[args.combo]:
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
        bgtmp.SetFillColor(r.kRed)
        bgtmp.SetFillStyle(3444)
        bgtmp.Draw("l3 same")
        graphs.append(bgtmp)
        
        # draw again on top
        gtmp.Draw("pz same")
        
        ctr += 1

for pformat in pformats:
    can.Print(oname+"."+pformat,pformat)

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
for pformat in pformats:
    can2.Print(oname2+"."+pformat,pformat)

