import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import getParamNames, paramVal, makeSigDict, getSigname

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-f","--file", dest="file", type=str, required=True, help="file w/ interpolated histograms")
parser.add_argument("-m","--method", dest="method", type=str, default="none", choices=["none","integral","moment"], help="interpolation method")
parser.add_argument("-r","--region", dest="region", type=str, required=True, help="signal region")
args = parser.parse_args()

import ROOT as r
r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

infile = r.TFile.Open(args.file)
rdir = infile.Get("{}_2018".format(args.region))
names = [k.GetName() for k in rdir.GetListOfKeys()]

xvarmin = 1500.
xvarmax = 8000.
mT = r.RooRealVar("mH{}_2018".format(args.region), "m_{T}", xvarmin, xvarmax, "GeV")
w = r.RooWorkspace("SVJ", "workspace")
for name in names:
    htmp = rdir.Get(name)
    name = name.replace("_"+args.method,"")
    rtmp = r.RooDataHist(name, name, r.RooArgList(mT), htmp, 1.)
    getattr(w, "import")(rtmp)

dirname = os.path.dirname(args.file)
outname = "{}/ws_{}_{}_2018_interp.root".format(dirname if len(dirname)>0 else ".",os.path.basename(args.file).replace(".root",""),args.region)
if os.path.exists(outname): os.remove(outname)
w.writeToFile(outname)
