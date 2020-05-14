import os,sys

region = sys.argv[1]

import ROOT as r
from array import array
from string import Template

def fill_template(inname, outname, **kwargs):
    if inname==outname:
        raise ValueError("Attempt to overwrite template file: "+inname)
    with open(inname,'r') as temp:
        old_lines = Template(temp.read())
        new_lines = old_lines.substitute(**kwargs)
    with open(outname,'w') as temp:
        temp.write(new_lines)

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

hfile = r.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/dijetbins/datacard.root")

# get histos

hdat = hfile.Get(region+"_2018/data_obs")
hdat.SetDirectory(0)
datrate = hdat.Integral()
hsig = hfile.Get(region+"_2018/SVJ_mZprime3000_mDark20_rinv03_alphapeak")
hsig.SetDirectory(0)
# fix sig bins
#for b in range(1,hsig.GetNbinsX()+1):
#    if hsig.GetBinContent(b)>0 and hsig.GetBinError(b)/hsig.GetBinContent(b)>0.5:
#        hsig.SetBinContent(b,0)
#        hsig.SetBinError(b,0)
sigrate = hsig.Integral()

infile = r.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/gaus/svj-dijet/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","READ")

w = infile.Get("SVJ")

mH = w.var("mH")
roomultipdf_norm = w.var("roomultipdf_norm")

# remove vars to be modified from list
vars = w.allVars()
vars.remove(mH)
vars.remove(roomultipdf_norm)

# modify var range
xmin = hdat.GetXaxis().GetXmin()
xmax = hdat.GetXaxis().GetXmax()
nbins = hdat.GetNbinsX()
mT = r.RooRealVar("mH","mT",xmin,xmax,"GeV")
mT.setBins(nbins)

# make roofit objs
pdat = r.RooDataHist("data_obs","data",r.RooArgList(mT),hdat,1.)
psig = r.RooDataHist("SVJ_mZprime3000_mDark20_rinv03_alphapeak","sig",r.RooArgList(mT),hsig,1.)

new_norm = r.RooRealVar("roomultipdf_norm","newnorm",datrate,0.,1.e6)

# add to list
vars.add(mT)
vars.add(new_norm)

ofname = "ws_svj_"+region+"_dijetbins.root"

# make new workspace, import contents from prev (but take modified mH var)
w2 = r.RooWorkspace("SVJd")
getattr(w2,"import")(vars)
getattr(w2,"import")(w.allPdfs())
getattr(w2,"import")(pdat)
getattr(w2,"import")(psig)

# write out modified workspace
if os.path.exists(ofname): os.remove(ofname)
w2.writeToFile(ofname)

# make datacard
fill_template(
    "SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template_bias.txt",
    "SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_bias.txt",
    sigrate = sigrate,
    datrate = datrate
)

