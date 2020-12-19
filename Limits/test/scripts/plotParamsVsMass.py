import shlex,subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from getParamsTracked import getParamsTracked, getFname
from collections import defaultdict
from array import array

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[name])")
parser.add_argument("-M", "--method", dest="method", type=str, default="AsymptoticLimits", help="method name (higgsCombineTest.[method])")
parser.add_argument("-q", "--quantile", dest="quantile", type=float, default=-1, help="quantile to plot fits")
parser.add_argument("-c", "--combo", dest="combo", type=str, required=True, help="combo to plot")
parser.add_argument("-b", "--batch", dest="batch", default=False, action="store_true", help="batch mode")
args = parser.parse_args()

masses = [
1500,
1700,
1900,
2100,
2300,
2500,
2700,
2900,
3100,
3300,
3500,
3700,
3900,
4100,
4300,
4500,
4700,
4900,
5100,
]

pformats = ["png"]
combos = {
"cut": ["highCut","lowCut"],
"bdt": ["highSVJ2","lowSVJ2"],
}

x = array('d', masses)

vals = defaultdict(list)
errs = defaultdict(list)
for mass in masses:
    fname = getFname(mass, args.name, args.method, args.combo)
    extraCondition = "abs(trackedParam_mZprime-{})<0.01".format(mass)
    params = getParamsTracked(fname, args.quantile, extraCondition=extraCondition)
    eparams = getParamsTracked(fname, args.quantile, includeParam=False, includeErr=True, extraCondition=extraCondition)
    for p in params:
        pname = p.replace("trackedParam_","")
        if not "Cut" in p and not "SVJ2" in p: continue
        vals[pname].append(params[p])
        errs[pname].append(eparams[p.replace("trackedParam_","trackedError_")])

# make graphs
oname = "params_vs_mass_{}_{}_q{}_{}".format(args.name,args.method,args.quantile,args.combo)
import ROOT as r
if args.batch: r.gROOT.SetBatch(True)
max_params = max(sum(region in k for k in vals) for region in combos[args.combo])
pad_size = 500
can = r.TCanvas(oname,"",max_params*pad_size,len(combos[args.combo])*pad_size)
can.Divide(max_params,len(combos[args.combo]))
ctr = 1
graphs = []
for region in combos[args.combo]:
    for p in sorted(vals):
        if not region in p: continue
        can.cd(ctr)
        y = array('d', vals[p])
        ye = array('d', errs[p])
        # hack to ignore error bars in axis
        gax = r.TGraph(len(masses),x,y)
        gax.SetMarkerStyle(20)
        gax.GetXaxis().SetTitle("m_{Z'} [GeV]")
        gax.GetYaxis().SetTitle(p)
        gax.SetTitle("")
        gax.Draw("ap")
        graphs.append(gax)
        gtmp = r.TGraphErrors(len(masses),x,y,r.nullptr,ye)
        gtmp.SetMarkerStyle(20)
        gtmp.Draw("pz same")
        graphs.append(gtmp)
        ctr += 1
for pformat in pformats:
    can.Print(oname+"."+pformat,pformat)

