import os, sys

region = sys.argv[1]
extra = len(sys.argv)>2 and sys.argv[2]!="-b"

import ROOT as r
import numpy as np

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
#r.TH1.AddDirectory(r.kFALSE)

if extra: ename = "svjOpt_"
else: ename = ""

sregions = ["dijet","lowCut","highCut","lowSVJ2","highSVJ2"]
sigmas = [2.0,4.0]
colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
styles = [1,7]
markers = [20,24]
hists = []
graphs_mu = []
graphs_sigma = []
ymin_mu = 1000
ymax_mu = -1000
ymin_sigma = 1000
ymax_sigma = -1000
haxis = r.TH1F("axis","",50,-5,5)
haxis.GetYaxis().SetRangeUser(0,40)
haxis.GetXaxis().SetTitle("b = (r-0)/rErr")
can = r.TCanvas()
can.cd()
haxis.Draw()
can.Update()
canh = r.TCanvas("canh","canh")
canh.cd()
haxis.Draw()
canh.Update()
for iSigma,sigma in enumerate(sigmas):
    bias_mu = []
    bias_sigma = []
    for iRegion,sregion in enumerate(sregions):
        print sregion
        hname = "{}_gaus_r{}_s{}".format(region,sregion,sigma)
        sfile = r.TFile.Open("fitDiagnosticsbias_"+ename+hname+".root")
        tree = sfile.Get("tree_fit_sb")
        hist = r.TH1F(hname,"",50,-5,5)
        tree.Draw("(r-0)/rErr>>"+hname,"fit_status==0","goff")
        hist.SetLineColor(colors[iRegion])
        hist.SetLineStyle(styles[iSigma])
        canh.cd()
        hist.Draw("hist same")
        hist.SetDirectory(0)
        hists.append(hist)
        gaus = r.TF1(hname,"gaus",-5,5)
        hist.Fit(gaus,"N")
        bias_mu.append(gaus.GetParameter(1))
        bias_sigma.append(gaus.GetParameter(2))
        gaus.SetLineColor(colors[iRegion])
        gaus.SetLineStyle(styles[iSigma])
        can.cd()
        gaus.Draw("same")
    iRegions = np.array(range(len(sregions)),dtype="float")
    bias_mu = np.array(bias_mu)
    bias_sigma = np.array(bias_sigma)
    graph_mu = r.TGraph(len(iRegions),iRegions,bias_mu)
    graph_mu.SetTitle("")
    graph_mu.GetXaxis().SetTitle("signal:data ratio")
    graph_mu.GetYaxis().SetTitle("#mu_{b}")
    graph_mu.SetMarkerStyle(markers[iSigma])
    graphs_mu.append(graph_mu)
    ymin_mu = min(ymin_mu,min(bias_mu))
    ymax_mu = max(ymax_mu,max(bias_mu))
    graph_sigma = r.TGraph(len(iRegions),iRegions,bias_sigma)
    graph_sigma.SetTitle("")
    graph_sigma.GetXaxis().SetTitle("signal:data ratio")
    graph_sigma.GetYaxis().SetTitle("#sigma_{b}")
    graph_sigma.SetMarkerStyle(markers[iSigma])
    graphs_sigma.append(graph_sigma)
    ymin_sigma = min(ymin_sigma,min(bias_sigma))
    ymax_sigma = max(ymax_sigma,max(bias_sigma))

can.Update()
can.Print("plotBiases_dijet-svj_"+ename+region+".png","png")
canh.Update()
canh.Print("plotBiases_hist_dijet-svj_"+ename+region+".png","png")

cang = r.TCanvas()
for ig,g in enumerate(graphs_mu):
    if ig==0:
        g.GetYaxis().SetRangeUser(ymin_mu,ymax_mu)
        g.Draw("ap")
    else: g.Draw("p same")
cang.Print("biasMu_svj_"+ename+region+".png","png")

for ig,g in enumerate(graphs_sigma):
    if ig==0:
        g.GetYaxis().SetRangeUser(ymin_sigma,ymax_sigma)
        g.Draw("ap")
    else: g.Draw("p same")
cang.Print("biasSigma_svj_"+ename+region+".png","png")
