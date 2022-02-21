from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import alphaVal, getParamNames, makeSigDict, getSigname, getSignameShort, paramVal

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-r", "--region", dest="region", type=str, required=True, help="region name")
parser.add_argument("-d", "--dir", dest="dir", type=str, default="root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig7/sigfull", help="datacard directory")
parser.add_argument("-s", "--signals", dest="signals", type=str, required=True, help="file w/ list of signals")
args = parser.parse_args()

import ROOT as r

with open('dict_xsec_Zprime.txt','r') as xfile:
    xsecs = {xline.split('\t')[0]: float(xline.split('\t')[1]) for xline in xfile}
lumi = 35920+41530+59740

param_names = getParamNames()+["xsec"]
param_values = []
with open(args.signals,'r') as sfile:
    for line in sfile:
        line = line.rstrip()
        if len(line)==0: continue
        param_values.append(line.split())
        param_values[-1].append(xsecs[param_values[-1][0]])

signals = [makeSigDict(param_values[i],param_names) for i in range(len(param_values))]

# make a Combine-esque tree in order to reuse plotLimit code
base_qtys = ["quantileExpected","limit"]
qtys = base_qtys + ["trackedParam_{}".format(q) for q in param_names]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")
qobj = r.quantile_t()
qobj.quantileExpected = 0.5
qobj.trackedParam_xsec = 1.0

file = r.TFile.Open("sigAccEff_{}.root".format(args.region),"RECREATE")
tree = r.TTree("limit","limit")
for qty in qtys:
    tree.Branch(qty, r.AddressOf(qobj,qty), '{}/D'.format(qty))

for signal in signals:
    sfile = r.TFile.Open("{}/datacard_{}.root".format(args.dir,getSignameShort(signal)))
    shist = sfile.Get("{}_2018/{}".format(args.region,getSigname(signal)))
    qobj.limit = shist.Integral(-1,-1)/(lumi*signal["xsec"])
    for q in getParamNames():
        setattr(qobj,"trackedParam_{}".format(q),paramVal(q,signal[q]))
    # skip xsec
    tree.Fill()

file.Write()
file.Close()
