import os,sys

region = sys.argv[1]

import ROOT as r

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

fname = "ws_SVJ_mZprime3100_mDark20_rinv03_alphapeak_"+region+"_2018_template.root"
infile = r.TFile.Open(fname,"READ")
w = infile.Get("SVJ")
mH = w.var("mH")
infile.Close()

hfile = r.TFile.Open("datacard_final_SVJ_3100_20_0.3_peak.root","READ")
hbkg = hfile.Get(region+"_2018/Bkg_norm")
rbkg = r.RooDataHist("Bkg", "MC Bkg",  r.RooArgList(mH), hbkg, 1.)

getattr(w, "import")(rbkg)

# write out modified workspace
os.remove(fname)
w.writeToFile(fname)

