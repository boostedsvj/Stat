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

infile = r.TFile.Open("ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","READ")

w = infile.Get("SVJ")

mH = w.var("mH")
roomultipdf_norm = w.var("roomultipdf_norm")

# remove vars to be modified from list
vars = w.allVars()
vars.remove(mH)
vars.remove(roomultipdf_norm)

# modify var range
xmin = 1500.
xmax = 3000.
nbins = 30
mT = r.RooRealVar("mH","mT",xmin,xmax,"GeV")
mT.setBins(nbins)

# modify histo range
hdat = w.data("data_obs").createHistogram("data_obs",mT)
hsig = w.data("SVJ_mZprime3000_mDark20_rinv03_alphapeak").createHistogram("SVJ_mZprime3000_mDark20_rinv03_alphapeak",mT)
# remove overflow
hdat.SetBinContent(nbins+1,0)
hdat.SetBinError(nbins+1,0)
datrate = hdat.Integral()
hsig.SetBinContent(nbins+1,0)
hsig.SetBinError(nbins+1,0)
sigrate = hsig.Integral()

# make roofit objs
pdat = r.RooDataHist("data_obs","data",r.RooArgList(mT),hdat,1.)
psig = r.RooDataHist("SVJ_mZprime3000_mDark20_rinv03_alphapeak","sig",r.RooArgList(mT),hsig,1.)

new_norm = r.RooRealVar("roomultipdf_norm","newnorm",datrate,0.,1.e6)

# add to list
vars.add(mT)
vars.add(new_norm)

ofname = "ws_svj_"+region+"_3000.root"

# make new workspace, import contents from prev (but take modified mH var)
w2 = r.RooWorkspace("SVJ3000")
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

