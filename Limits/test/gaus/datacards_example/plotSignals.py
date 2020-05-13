import os
import ROOT as r

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

infile = r.TFile.Open("dijet_combine_qq_1900_lumi-137.500_PFDijet2017_mod.root")
w = infile.Get("wPFDijet2017")
th1x = w.var("th1x")
data_obs = w.data("data_obs")

means = [4.21729, 14.987]
sigmas = [1.0,2.0,3.0,4.0,5.0]
colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
for iMean,mean in enumerate(means):
    can = r.TCanvas()
    frame = th1x.frame()
    frame.SetTitle("")
    data_obs.plotOn(frame,r.RooFit.MarkerColor(r.kGray+2),r.RooFit.LineColor(r.kGray+2))
    for iSigma,sigma in enumerate(sigmas):
        hname = "gaus_m{}_s{}".format(iMean,iSigma)
        sig = w.data(hname)
        sig.plotOn(frame,r.RooFit.MarkerColor(colors[iSigma]),r.RooFit.LineColor(colors[iSigma]))
    frame.Draw()
    frame.GetYaxis().SetRangeUser(1,2e6)
    can.SetLogy()
    can.Update()
    can.Print("plotSignals_dijet_m{}.png".format(iMean),"png")
