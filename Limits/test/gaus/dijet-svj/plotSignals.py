import os, sys

region = sys.argv[1]

from collections import OrderedDict
import ROOT as r

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

infile = r.TFile.Open("dijet_combine_{}_PFDijet2017.root".format(region))
w = infile.Get("wPFDijet2018")
th1x = w.var("th1x")
data_obs = w.data("data_obs")

ymax = OrderedDict([
("dijet", 2e6),
("lowCut", 3e4),
("highCut" , 4e3),
("lowSVJ2", 400),
("highSVJ2", 50),
])
sigmas = [2.0,4.0]
colors = [r.kBlack, r.kMagenta, r.kRed, r.kBlue, r.kCyan]
styles = [20,24]

can = r.TCanvas()
frame = th1x.frame()
frame.SetTitle("")
data_obs.plotOn(frame,r.RooFit.MarkerColor(r.kGray+2),r.RooFit.LineColor(r.kGray+2))

for iRegion,sregion in enumerate(ymax.keys()):
    for iSigma,sigma in enumerate(sigmas):
        hname = "{}_gaus_r{}_s{}".format(region,sregion,sigma)
        sig = w.data(hname)
        sig.plotOn(frame,r.RooFit.MarkerColor(colors[iRegion]),r.RooFit.LineColor(colors[iRegion]),r.RooFit.MarkerStyle(styles[iSigma]))
frame.Draw()
frame.GetYaxis().SetRangeUser(0.001,ymax[region])
can.SetLogy()
can.Update()
can.Print("plotSignals_dijet-svj_{}.png".format(region),"png")
