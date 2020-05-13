import os,sys

extra = len(sys.argv)>1

import ROOT as r
import numpy as np

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
#r.TH1.AddDirectory(r.kFALSE)

if extra: ename = "_svjOpt"
else: ename = ""

means = [4.21729]
sigmas = [1.0,2.0,3.0,4.0,5.0]
colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
haxis = r.TH1F("axis","",50,-5,5)
haxis.GetYaxis().SetRangeUser(0,25)
haxis.GetXaxis().SetTitle("b = (r-0)/rErr")
bias_mu = []
bias_sigma = []
for iMean,mean in enumerate(means):
    can = r.TCanvas()
    can.cd()
    haxis.Draw()
    can.Update()
    for iSigma,sigma in enumerate(sigmas):
        hname = "gaus_m{}_s{}".format(iMean,iSigma)
        sfile = r.TFile.Open("fitDiagnosticsbias"+ename+"_"+hname+".root")
        tree = sfile.Get("tree_fit_sb")
        hist = r.TH1F(hname,"",50,-5,5)
        tree.Draw("(r-0)/rErr>>"+hname,"","goff")
        gaus = r.TF1(hname,"gaus",-5,5)
        hist.Fit(gaus,"N")
        bias_mu.append(gaus.GetParameter(1))
        bias_sigma.append(gaus.GetParameter(2))
        gaus.SetLineColor(colors[iSigma])
        gaus.Draw("same")
    can.Update()
    can.Print("plotBiases_dijet"+ename+"_m{}.png".format(iMean),"png")

sigmas = np.array(sigmas)
bias_mu = np.array(bias_mu)
bias_sigma = np.array(bias_sigma)
cang = r.TCanvas()

g_mu = r.TGraph(len(sigmas),sigmas,bias_mu)
g_mu.SetTitle("")
g_mu.GetXaxis().SetTitle("#sigma_{signal}")
g_mu.GetYaxis().SetTitle("#mu_{b}")
g_mu.Draw("ap")
cang.Print("biasMu_dijet"+ename+"_m{}.png".format(iMean),"png")

g_sigma = r.TGraph(len(sigmas),sigmas,bias_sigma)
g_sigma.SetTitle("")
g_sigma.GetXaxis().SetTitle("#sigma_{signal}")
g_sigma.GetYaxis().SetTitle("#sigma_{b}")
g_sigma.Draw("ap")
cang.Print("biasSigma_dijet"+ename+"_m{}.png".format(iMean),"png")


