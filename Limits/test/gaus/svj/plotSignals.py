import os, sys

region = sys.argv[1]

import ROOT as r

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

infile = r.TFile.Open("ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template_mod.root")
w = infile.Get("SVJ")
mH = w.var("mH")
data_obs = w.data("data_obs")

ymax = {
"highCut" : 2e3,
"highSVJ2": 50,
"lowCut": 3e4,
"lowSVJ2": 400,
}
means = [2141]
sigmas = [1.0,2.0,3.0,4.0,5.0]
colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
for iMean,mean in enumerate(means):
    can = r.TCanvas()
    frame = mH.frame()
    frame.SetTitle("")
    data_obs.plotOn(frame,r.RooFit.MarkerColor(r.kGray+2),r.RooFit.LineColor(r.kGray+2))
    for iSigma,sigma in enumerate(sigmas):
        hname = "gaus_m{}_s{}".format(iMean,iSigma)
        sig = w.data(hname)
        sig.plotOn(frame,r.RooFit.MarkerColor(colors[iSigma]),r.RooFit.LineColor(colors[iSigma]))
    frame.Draw()
    frame.GetYaxis().SetRangeUser(1,ymax[region])
    can.SetLogy()
    can.Update()
    can.Print("plotSignals_svj_"+region+"_m{}.png".format(iMean),"png")
