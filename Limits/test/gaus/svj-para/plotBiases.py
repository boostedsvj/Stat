import os, sys

region = sys.argv[1]

import ROOT as r
import numpy as np

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
#r.TH1.AddDirectory(r.kFALSE)

colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
haxis = r.TH1F("axis","",50,-5,5)
haxis.GetYaxis().SetRangeUser(0,25)
haxis.GetXaxis().SetTitle("b = (r-0)/rErr")
bias_mu = []
bias_sigma = []
can = r.TCanvas()
can.cd()
haxis.Draw()
can.Update()
hists = []
for iopt,ename in enumerate(["","svjOpt_"]):
    sfile = r.TFile.Open("fitDiagnosticsbias_"+ename+region+".root")
    tree = sfile.Get("tree_fit_sb")
    hname = region+ename
    hist = r.TH1F(hname,"",50,-5,5)
    tree.Draw("(r-0)/rErr>>"+hname,"fit_status==0","goff")
    gaus = r.TF1(hname,"gaus",-5,5)
    hist.Fit(gaus,"NQ")
    bias_mu.append(gaus.GetParameter(1))
    bias_sigma.append(gaus.GetParameter(2))
    gaus.SetLineColor(colors[iopt])
    gaus.SetLineWidth(2)
    gaus.Draw("same")
    hist.SetDirectory(0)
    hists.append(hist)
    hist.SetLineColor(colors[iopt])
    hist.SetLineStyle(7)
    hist.Draw("hist same")
can.Update()
can.Print("plotBiases_svj_"+region+".png","png")

for vals in [bias_mu, bias_sigma]:
    print region+"\t"+"\t".join([str(x) for x in vals])
