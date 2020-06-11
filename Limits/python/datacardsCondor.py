import ROOT
import json

#changed to no longer need the settings.py file
# -syst list is created here
# -rateParams list is created here
# -processes list is created here
# files now saved to base directory and exported to EOS after jobs completion


from ROOT import RooRealVar, RooDataHist, RooArgList, RooGenericPdf, RooBernstein, RooExtendPdf, RooCmdArg, RooWorkspace, RooFit, RooDataSet, RooArgSet, RooCategory, RooFitResult, RooCurve, RooParametricShapeBinPdf
import os, sys
from array import array
import copy, math, pickle
import collections
import numpy as np
from numpy import ndarray

ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetPadColor(0)


syst = collections.OrderedDict()
syst["lumi"] = ("lnN", "all", 1.10)

rateParams = {}
rateParams["lowSVJ1_2018"] = "TMath::Power(TMath::Range(0.01,0.99,@0),1)*TMath::Power(1-TMath::Range(0.01,0.99,@0*%s),1)/(TMath::Power(1-%s,1))"
rateParams["lowSVJ2_2018"] = "TMath::Power(TMath::Range(0.01,0.99,@0),2)*TMath::Power(1-TMath::Range(0.01,0.99,@0*%s),0)/(TMath::Power(1-%s,0))"
rateParams["highSVJ1_2018"] = "TMath::Power(TMath::Range(0.01,0.99,@0),1)*TMath::Power(1-TMath::Range(0.01,0.99,@0*%s),1)/(TMath::Power(1-%s,1))"
rateParams["highSVJ2_2018"] = "TMath::Power(TMath::Range(0.01,0.99,@0),2)*TMath::Power(1-TMath::Range(0.01,0.99,@0*%s),0)/(TMath::Power(1-%s,0))"


processes = ["QCD"]

#ROOT.Math.MinimizerOptions.SetDefaultTolerance(1e-3); 
#ROOT.Math.MinimizerOptions.SetDefaultPrecision(1e-8)
#*****************************************************
#                                                       #
#   getRate(process, ifile)                             #
#                                                       #
#   getCard(sig, ch, ifilename, outdir)                 #
#                                                       #
#*******************************************************#


#*******************************************************#
#                                                       #
#                     Utility Functions                 #
#                                                       #
#*******************************************************#

isData = True

def getRate(ch, process, ifile):
       hName = ch + "/"+ process
       h = ifile.Get(hName)
       #return h.Integral()
       return h.Integral(1,h.GetXaxis().GetNbins()-1)
       #return h.Integral(1,91)

def getHist(ch, process, ifile):
       hName = ch + "/"+ process
       print "Histo Name ", hName
       h = ifile.Get(hName)
       h.SetDirectory(0)
       return h




def getEfficiency(sig, channels, ifilename):

       try:
              ifile = ROOT.TFile.Open(ifilename)
       except IOError:
              print "Cannot open ", ifilename
       else:
              print "Opening file ",  ifilename
              ifile.cd()
              
       cdfilename = os.getcwd()+"/"+ifilename
       
       rates = {}
       for ch in channels:

              rates[ch] = getRate(ch, sig, ifile)

       sumjets = sum(rates.itervalues())*2
       print sumjets

       evts2SVJ = sum(i for ch, i in rates.iteritems() if "SVJ2" in ch)
       print evts2SVJ
       evts1SVJ = sum(i for ch, i in rates.iteritems() if "SVJ1" in ch)

       num = evts1SVJ + (2 * evts2SVJ)
       
       eff = num/sumjets
       return eff




def getRSS(sig, ch, variable, model, dataset, fitRes, carddir,  norm = -1, label = "nom"):
       name = model.GetName()
       order = int(name[-1])
       npar = fitRes[0].floatParsFinal().getSize() if len(fitRes)>0 else 0
       varArg = ROOT.RooArgSet(variable)
      
       #frame = variable.frame()
       frame = variable.frame(ROOT.RooFit.Title(""))
       dataset.plotOn(frame, RooFit.Invisible())
       print(fitRes[0])
       #if len(fitRes) > 0: model.plotOn(frame, RooFit.VisualizeError(fitRes[0], 1, False), RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("FL"), RooFit.Range("Full"))
       if len(fitRes) > 0: graphFit = model.plotOn(frame, RooFit.VisualizeError(fitRes[0], 1, False), RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("FL"), RooFit.Range("Full"))

       model.plotOn(frame, RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("L"), RooFit.Name(model.GetName()),  RooFit.Range("Full"))
       model.paramOn(frame, RooFit.Label(model.GetTitle()), RooFit.Layout(0.45, 0.95, 0.94), RooFit.Format("NEAU"))
       
       graphData = dataset.plotOn(frame, RooFit.DataError(ROOT.RooAbsData.Poisson if isData else ROOT.RooAbsData.SumW2), RooFit.DrawOption("PE0"), RooFit.Name(dataset.GetName()))

       #(ROOT.TVirtualFitter.GetFitter()).GetConfidenceIntervals(model)       

       #f1 = ROOT.TF1(model)
       #f1.SetName("f")
       #f2 = ROOT.TF1(model)
       #f2.SetName("g")
       #nPar = f2.GetNpar()
       #errorpar = np.zeros(nPar)

       #f2.SetParErrors()
       #f3 = ROOT.TF1("error", "f/g")
       
       pulls = frame.pullHist(dataset.GetName(), model.GetName(), True)
       residuals = frame.residHist(dataset.GetName(), model.GetName(), False, True) # this is y_i - f(x_i)
    
       roochi2 = frame.chiSquare(model.GetName(), dataset.GetName(),npar)#dataset.GetName(), model.GetName()) #model.GetName(), dataset.GetName()
       print "forcing bins: 90"
       nbins = 90
       chi = roochi2 * ( nbins - npar)
       print "pls: ", chi,  nbins
       roopro = ROOT.TMath.Prob(chi, nbins - npar)

       frame.SetMaximum(frame.GetMaximum()*10.)
       frame.SetMinimum(0.1)
       print "==========> len(sig): ", len(sig)
       length = 1
       if(length<2):

              c = ROOT.TCanvas("c_"+ch+model.GetName(), ch, 800, 800)

              c.cd()
              pad1 = ROOT.TPad("pad1", "pad1", 0., 0.35, 1., 1.0)
              ROOT.SetOwnership(pad1, False)
              pad1.SetBottomMargin(0.)
              pad1.SetGridx()
              pad1.SetGridy()
              pad1.SetLogy()
              pad1.Draw()
              pad1.cd()

              frame.Draw()
              frame.SetTitle("")
              
              #txt = ROOT.TText(2500, 1000, "chiSquared: " + str(roochi2))
              #txt = ROOT.TText(2500, 6., "chiSquared: " + str(roochi2))
              roochi2_small = format(roochi2, '.2f')
              txt = ROOT.TText(0.65, 0.7, "chiSquared: " + str(roochi2_small))              
              txt.SetNDC()
              txt.SetTextSize(0.04) 
              txt.SetTextColor(ROOT.kRed) 
              txt.Draw();

              #txt_2 = ROOT.TText(2500, 800, "Prob: " + str(roopro))
              #txt_2 = ROOT.TText(2500, 2., "Prob: " + str(roopro))
              roopro_small = format(roopro, '.2f')
              txt_2 = ROOT.TText(0.65, 0.65, "prob: " + str(roopro_small))              
              txt_2.SetNDC()
              txt_2.SetTextSize(0.04) 
              txt_2.SetTextColor(ROOT.kRed) 
              txt_2.Draw();
              
              c.cd()
              c.Update()
              c.cd()
              pad2 = ROOT.TPad("pad2", "pad2", 0, 0.1, 1, 0.3)
              ROOT.SetOwnership(pad2, False)
              pad2.SetTopMargin(0);
              pad2.SetBottomMargin(0.25);
              pad2.SetGridx();
              pad2.SetGridy();
              pad2.Draw();
              pad2.cd()
              pad2.Clear()
       #frame.Draw()
              c.Update()
              c.Modified()

              frame_res = variable.frame(ROOT.RooFit.Title(""))
              frame_res.addPlotable(pulls, "P")
              frame_res.GetYaxis().SetRangeUser(-3.5, 3.5)
              frame_res.SetTitle("")
              frame_res.GetYaxis().SetTitle("(N^{data}-f(i))/#sigma")
              frame_res.GetYaxis().SetTitleOffset(.3)
              frame_res.GetXaxis().SetTitleOffset(1.)
              frame_res.GetYaxis().SetTitleSize(.13)
              frame_res.GetXaxis().SetTitle("m_{T} (GeV)")
              frame_res.GetXaxis().SetTitleSize(.13)
              frame_res.GetXaxis().SetLabelSize(.12)
              frame_res.GetYaxis().SetLabelSize(.12)
              frame_res.Draw("e0SAME")
              frame_res.SetTitle("")

       hist = graphData.getHist()
       xmin_ = graphData.GetXaxis().GetXmin()
       xmax_ = graphData.GetXaxis().GetXmax()
       
       length = 1
       if(length<2):
              c.Update()
              c.Range(xmin_, -3.5, xmax_, 3.5)
              line = ROOT.TLine(xmin_, 0., xmax_, 0.)
              line.SetLineColor(ROOT.kRed)
              line.SetLineWidth(2)
              line.Draw("same")
                            
              c.SaveAs("Residuals_"+ch+"_"+name+"_log.pdf")              
              '''             
              c2 = ROOT.TCanvas("c2_"+ch+model.GetName(), ch, 800, 800)
              c2.cd()
              pad1_2 = ROOT.TPad("pad1_2", "pad1_2", 0., 0.35, 1., 1.0)
              ROOT.SetOwnership(pad1_2, False)
              pad1_2.SetBottomMargin(0.);
              pad1_2.SetGridx();
              pad1_2.SetGridy();
              pad1_2.SetLogy()
              pad1_2.Draw()
              pad1_2.cd()
              frame.Draw()
              
              txt = ROOT.TText(3000, 4., "chiSquared: " + str(roochi2))
              txt.SetTextSize(0.04) 
              txt.SetTextColor(ROOT.kRed) 
              txt.Draw()
              
              c2.cd()
              c2.Update()
              c2.cd()
              pad2_2 = ROOT.TPad("pad2_2", "pad2_2", 0, 0.05, 1, 0.3)
              ROOT.SetOwnership(pad2_2, False)
              pad2_2.SetTopMargin(0);
              pad2_2.SetBottomMargin(0.1);
              pad2_2.SetGridx();
              pad2_2.SetGridy();
              pad2_2.Draw();
              pad2_2.cd()
              pad2_2.Clear()
       #frame.Draw()
              c2.Update()
              c2.Modified()

              '''
       # calculate RSS
       res, rss, chi1, chi2 = 0, 0., 0., 0., 

       xmin, xmax = array('d', [0.]), array('d', [0.])
       dataset.getRange(variable, xmin, xmax)
#       print "N bins: ", hist.GetN()
       #nBins = graphData.GetNbinsX()
       ratioHist = ROOT.TH1F("RatioPlot", "Ratio Plot", hist.GetN(), xmin_, xmax_)
       ROOT.SetOwnership(ratioHist, False)

 #      ratioHist.SetBins(nBins, graphData.GetXaxis().GetXbins().GetArray())

       #hist.Dump()
       #ratioHist.Dump()
       #ratioHist.Print('all')

       sumErrors = 0
       graphFit.Print()
       fitmodel = graphFit.getHist()
       fitmodel_curve = graphFit.getCurve("Bkg_"+str(ch)+str(order)+"_rgp_Norm[mH]_errorband")
       if label == "alt":
              fitmodel_curve = graphFit.getCurve("Bkg_Alt_"+str(ch)+str(order)+"_rgp_Norm[mH]_errorband")
       print "fitmodel_curve: "
       fitmodel_curve.Print()

       c_x = np.ndarray(fitmodel_curve.GetN(), 'd', fitmodel_curve.GetX())
       c_y_all = np.ndarray(fitmodel_curve.GetN(), 'd', fitmodel_curve.GetY())
       c_y_up = c_y_all[0:len(c_y_all)/2]
       c_y_down = list(reversed(c_y_all[len(c_y_all)/2+1:]))
       
       def findClosestLowerBin(x, bins):
              for ib in range(len(bins)):
                     if bins[ib] > x:
                            return ib-1

       
       for i in xrange(0, hist.GetN()):
              
              print "bin x and y: ",  hist.GetX()[i], hist.GetY()[i]
              print hist.GetN(), fitmodel_curve.GetN()
              err_up = c_y_up[findClosestLowerBin(hist.GetX()[i], c_x)]
              err_down = c_y_down[findClosestLowerBin(hist.GetX()[i], c_x)]
              
              print "x, y, eyup, eydo: ", hist.GetX()[i], hist.GetY()[i], err_up, err_down
              if(hist.GetY()[i]>0):
                     diffErrors = abs(err_up - err_down)/hist.GetY()[i]
              else:
                     diffErrors = 0
              sumErrors += diffErrors
              
              hx, hy = hist.GetX()[i], hist.GetY()[i]
              hey = hist.GetErrorY(i)
              hexlo, hexhi = hist.GetErrorXlow(i), hist.GetErrorXhigh(i)
              ry = residuals.GetY()[i]
              pull = pulls.GetY()[i]
              resi = residuals.GetY()[i]

              if (hx - hexlo) > xmax[0] and hx + hexhi > xmax[0]:
                     continue

              if hy <= 0.:
                     continue
              
              res += ry
              rss += ry**2 

              chi1 += abs(pull)
              chi2 += abs((resi**2)/pull)


              ratioHist.SetBinContent(i+1, (ry - hy)/(-1*hy))
              ratioHist.SetBinError(i+1, (hey/ hy**2))
              ratioHist.Print('all')

       print "=======> sumErrors: ", sumErrors
        
       rss = math.sqrt(rss)
       parValList = []
       print(len(fitRes[0].floatParsFinal()))
       for iPar in range(len(fitRes[0].floatParsFinal())):
              print(iPar)
              parValList.append((fitRes[0].floatParsFinal().at(iPar)).getValV())
       out = {"chiSquared":roochi2,"chi2" : chi2, "chi1" : chi1, "rss" : rss, "res" : res, "nbins" : hist.GetN(), "npar" : npar, "parVals": parValList}
       length=1
       if(length<2):

              c2 = ROOT.TCanvas("c2_"+ch+model.GetName(), ch, 800, 800)
              c2.cd()
              pad1_2 = ROOT.TPad("pad1_2", "pad1_2", 0., 0.35, 1., 1.0)
              ROOT.SetOwnership(pad1_2, False)
              pad1_2.SetBottomMargin(0.);
              pad1_2.SetGridx();
              pad1_2.SetGridy();
              pad1_2.SetLogy()
              pad1_2.Draw()
              pad1_2.cd()
              frame.Draw()
              frame.SetTitle("")

              #txt = ROOT.TText(2500, 6., "chiSquared: " + str(roochi2))
              #txt = ROOT.TText(2500, 1000, "chiSquared: " + str(roochi2))
              txt = ROOT.TText(0.65, 0.7, "chiSquared: " + str(roochi2_small))              
              txt.SetNDC()
              txt.SetTextSize(0.04) 
              txt.SetTextColor(ROOT.kRed) 
              txt.Draw()


              #txt_2 = ROOT.TText(2500, 800, "Prob: " + str(roopro))
              txt_2 = ROOT.TText(0.65, 0.65, "prob: " + str(roopro_small))              
              txt_2.SetNDC()
              #txt_2 = ROOT.TText(2500, 2., "Prob: " + str(roopro))
              txt_2.SetTextSize(0.04) 
              txt_2.SetTextColor(ROOT.kRed) 
              txt_2.Draw();


#              txt2 = ROOT.TText(3000, 1., "Error: "+  str(sumErrors))
#              txt2.SetTextSize(0.04) 
#              txt2.SetTextColor(ROOT.kOrange) 
              #txt2.Draw("same");
              
              c2.cd()
              c2.Update()
              c2.cd()
              pad2_2 = ROOT.TPad("pad2_2", "pad2_2", 0, 0.05, 1, 0.3)
              ROOT.SetOwnership(pad2_2, False)
              pad2_2.SetTopMargin(0);
              pad2_2.SetBottomMargin(0.25);
              pad2_2.SetGridx();
              pad2_2.SetGridy();
              pad2_2.Draw();
              pad2_2.cd()
              pad2_2.Clear()
       #frame.Draw()
              c2.Update()
              c2.Modified()
              frame.SetTitle("")

              ratioHist.GetYaxis().SetRangeUser(0.,2.)
              ratioHist.SetTitle("")
              ratioHist.GetYaxis().SetTitle("f(i)/N^{data}")
              ratioHist.GetXaxis().SetTitle("m_{T} (GeV)")
              ratioHist.GetYaxis().SetTitleSize(.13)
              ratioHist.GetYaxis().SetTitleOffset(.3)
              ratioHist.GetXaxis().SetTitleOffset(1.)
              ratioHist.GetXaxis().SetTitleSize(.13)
              ratioHist.GetXaxis().SetLabelSize(.12)
              ratioHist.GetYaxis().SetLabelSize(.12)

       #ratioHist.SetMarkerSize(2)
              ratioHist.SetMarkerStyle(20)
              ratioHist.Draw("PE")
              line2 = ROOT.TLine(xmin_, 1., xmax_, 1.)
              line2.SetLineColor(ROOT.kRed)
              line2.SetLineWidth(2)
              line2.Draw("same")
              frame.SetTitle("")

              c2.SaveAs("Residuals_"+ch+"_"+name + "_ratio_log.pdf")

       return out


def fisherTest(RSS1, RSS2, o1, o2, N):
       print "Testing functions with parameters o1 ", o1, " and o2 ", o2, " with RSS RSS1 ", RSS1, " and RSS2 ", RSS2, " and N ", N
       n1 = o1
       n2 = o2
       F = ((RSS1-RSS2)/(n2-n1)) / (RSS2/(N-n2))
       CL =  1.-ROOT.TMath.FDistI(F, n2-n1, N-n2)
       print "***************** Value of F, CL: ", F, CL
       return CL, F


#*******************************************************#
#                                                       #
#                      Datacard                         #
#                                                       #
#*******************************************************#

def getCard(sig, ch, ifilename, outdir, doModelling, mode = "histo", bias = False, verbose = False):


       ### Cmpute signal efficiencies 
       
       effname = "Efficiencies.txt"
       efile = open(effname, 'r')
       effline = efile.readline()

       effs = json.loads(effline)
       
       eff = effs[sig]

       try:
              ifile = ROOT.TFile.Open(ifilename)
       except IOError:
              print "Cannot open ", ifilename
       else:
              print "Opening file ",  ifilename
              ifile.cd()


       print "BIAS?", bias
       #workdir_ = ifilename.split("/")[:-1]
       #WORKDIR = "/".join(workdir_) + "/"
       carddir = outdir+  sig + "/"
       #if not os.path.isdir(outdir): os.system('mkdir ' +outdir)
       #if not os.path.isdir(outdir + "/" + sig): os.system('mkdir ' +carddir) 
       #if not os.path.isdir(outdir + carddir + "plots/"): os.system('mkdir ' +carddir + "plots/") 
       #if not os.path.isdir(outdir + carddir + "Fisher/"): os.system('mkdir ' +carddir + "Fisher/") 
       #if not os.path.isdir(outdir + carddir + "Residuals/"): os.system('mkdir ' + carddir + "Residuals/") 

       hist_filename = os.getcwd()+"/"+ifilename
       hist = getHist(ch, sig, ifile)


       #*******************************************************#
       #                                                       #
       #                   Generate workspace                  #
       #                                                       #
       #*******************************************************#
    
       if(mode == "template"):

              ###ATT: CHECK BKG AND DATA NORMALIZATION AND DISTRIBUTION
              histBkgData = getHist(ch, "Bkg", ifile)
              histData = getHist(ch, "data_obs", ifile)
              print "channel ", ch             
              print "channel ", sig

              histSig = getHist(ch, sig, ifile)
              print "histSigData: ", histSig.Integral()
              #print "histBkgData: ", histBkgData.Integral()
              xvarmin = 1500.
              xvarmax = 6000.
              mT = RooRealVar(  "mH",    "m_{T}", xvarmin, xvarmax, "GeV")
              binMax = histData.FindBin(xvarmax)
              bkgData = RooDataHist("bkgdata", "MC Bkg",  RooArgList(mT), histBkgData, 1.)
              obsData = RooDataHist("data_obs", "Data",  RooArgList(mT), histData, 1.)
              sigData = RooDataHist("sigdata", "MC Sig",  RooArgList(mT), histSig, 1.)
              print "SUM ENTRIES: ", sigData.sumEntries()
              print "Bkg Integral: ", histData.Integral() 
              #nBkgEvts = histBkgData.Integral(1, histBkgData.GetXaxis().GetNbins()-1) 
              nBkgEvts = histBkgData.Integral(1, binMax)
              nDataEvts = histData.Integral(1, binMax)
              nSigEvts = histSig.Integral(1, binMax)
#              nBkgEvts = histData.Integral(1, histData.GetXaxis().GetNbins()-1)
#              nDataEvts = histData.Integral(1, histData.GetXaxis().GetNbins()-1) 
              #print "Bkg Events: ", nBkgEvts

              print "channel: ", ch
              normBkg = RooRealVar("Bkg_"+ch+"_norm", "Number of background events", nBkgEvts, 0., 2.e4)
              normData = RooRealVar("Data_"+ch+"_norm", "Number of background events", nDataEvts, 0., 2.e4)
              ch_red = ch[:-5]
              modelName = "Bkg_"+ch
              modelAltName =  "Bkg_Alt_"+ch
              
              if(doModelling):
                     #ch_red = ch
                     print "channel: ", ch_red


                     p1_1 = RooRealVar(ch_red + "_p1_1", "p1", 1., -1000., 1000.)
                     p1_2 = RooRealVar(ch_red + "_p1_2", "p1", 1., -1000., 1000.)
                     p1_3 = RooRealVar(ch_red + "_p1_3", "p1", 1., -1000., 1000.)
                     p1_4 = RooRealVar(ch_red + "_p1_4", "p1", 1., -1000., 1000.)

                     p2_2 = RooRealVar(ch_red + "_p2_2", "p2", 1., -1000., 1000.)
                     p2_3 = RooRealVar(ch_red + "_p2_3", "p2", 1., -1000., 1000.)
                     p2_4 = RooRealVar(ch_red + "_p2_4", "p2", 1., -1000., 1000.)

                     p3_3 = RooRealVar(ch_red + "_p3_3", "p3", 1., -1000., 1000.)
                     p3_4 = RooRealVar(ch_red + "_p3_4", "p3", 1., -1000., 1000.)

                     p4_4 = RooRealVar(ch_red + "_p4_4", "p4", 1., -1000., 1000.)


                     #Function from Theorists, combo testing, sequence E, 1, 11, 12, 22
                     # model NM has N params on 1-x and M params on x. exponents are (p_i + p_{i+1} * log(x))
                     # these are the RooGenericPdf verisons, convert to RooParametricShapeBinPdf below
                     modelBkg1_rgp = RooGenericPdf(modelName+"1_rgp", "Thry. fit (01)", "pow(@0/13000, -@1)", RooArgList(mT, p1_1))
                     modelBkg2_rgp = RooGenericPdf(modelName+"2_rgp", "Thry. fit (11)", "pow(1 - @0/13000, @2) *pow(@0/13000, -@1)", RooArgList(mT, p1_2, p2_2))
                     modelBkg3_rgp = RooGenericPdf(modelName+"3_rgp", "Thry. fit (12)", "pow(1 - @0/13000, @2) * pow(@0/13000, -@1-@3*log(@0/13000))", RooArgList(mT, p1_3, p2_3, p3_3))
                     modelBkg4_rgp = RooGenericPdf(modelName+"4_rgp", "Thry. fit (22)", "pow(1 - @0/13000, @2+@4*log(@0/13000)) * pow(@0/13000, -@1-@3*log(@0/13000))", RooArgList(mT, p1_4, p2_4, p3_4, p4_4))

                     modelBkg1 = RooParametricShapeBinPdf(modelName+"1", "Thry. Fit (01)", modelBkg1_rgp, mT, RooArgList(p1_1), histBkgData)
                     modelBkg2 = RooParametricShapeBinPdf(modelName+"2", "Thry. Fit (11)", modelBkg2_rgp, mT, RooArgList(p1_2, p2_2), histBkgData)
                     modelBkg3 = RooParametricShapeBinPdf(modelName+"3", "Thry. Fit (12)", modelBkg3_rgp, mT, RooArgList(p1_3, p2_3, p3_3), histBkgData)
                     modelBkg4 = RooParametricShapeBinPdf(modelName+"4", "Thry. Fit (22)", modelBkg4_rgp, mT, RooArgList(p1_4, p2_4, p3_4, p4_4), histBkgData)
                     RSS = {}
                     fitrange = "Full"
                     #if (ch == "highSVJ1_2016"): fitrange = "Low,High"

                     fitRes1 = modelBkg1.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                     fitRes2 = modelBkg2.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                     fitRes3 = modelBkg3.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                     fitRes4 = modelBkg4.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                     orderBkg = [len(fitRes1.floatParsFinal()),len(fitRes2.floatParsFinal()),len(fitRes3.floatParsFinal()),len(fitRes4.floatParsFinal())]
                     
                     xframe = mT.frame(ROOT.RooFit.Title("extended ML fit example"))

                     c1 = ROOT.TCanvas()
                     c1.cd()
                     obsData.plotOn(xframe, RooFit.Name("obsData"))

                     modelBkg1.plotOn(xframe, RooFit.Name("modelBkg1"), ROOT.RooFit.LineColor(ROOT.kPink + 6), RooFit.Range("Full"))
                     modelBkg2.plotOn(xframe, RooFit.Name("modelBkg2"), ROOT.RooFit.LineColor(ROOT.kBlue -4), RooFit.Range("Full"))
                     modelBkg3.plotOn(xframe, RooFit.Name("modelBkg3"), ROOT.RooFit.LineColor(ROOT.kRed -4), RooFit.Range("Full"))
                     modelBkg4.plotOn(xframe, RooFit.Name("modelBkg4"), ROOT.RooFit.LineColor(ROOT.kGreen +1), RooFit.Range("Full"))
                     #modelBkg1.paramOn(xframe, RooFit.Name("modelBkg1"), RooFit.Layout(0.0,0.5,0.5))
                     #modelBkg2.paramOn(xframe, RooFit.Name("modelBkg2"), RooFit.Layout(0.5,1.0,0.5))
                     #modelBkg3.paramOn(xframe, RooFit.Name("modelBkg3"), RooFit.Layout(0.0,0.5,0.25))
                     #modelBkg4.paramOn(xframe, RooFit.Name("modelBkg4"), RooFit.Layout(0.5,1.0,0.25))

                     chi2s_bkg1 = xframe.chiSquare("modelBkg1", "obsData",1)
                     chi2s_bkg2 = xframe.chiSquare("modelBkg2", "obsData",2)
                     chi2s_bkg3 = xframe.chiSquare("modelBkg3", "obsData",3)
                     chi2s_bkg4 = xframe.chiSquare("modelBkg4", "obsData",4)

                     xframe.SetMinimum(0.0002)
                     xframe.Draw()

                     txt1 = ROOT.TText(2000., 10., "model1, nP {}, chi2: {}".format(orderBkg[0],chi2s_bkg1))
                     txt1.SetTextSize(0.04) 
                     txt1.SetTextColor(ROOT.kPink+6) 
                     xframe.addObject(txt1) 
                     txt2 = ROOT.TText(2000., 1, "model2, nP {}, chi2: {}".format(orderBkg[1],chi2s_bkg2))
                     txt2.SetTextSize(0.04) 
                     txt2.SetTextColor(ROOT.kBlue-4) 
                     xframe.addObject(txt2) 
                     txt3 = ROOT.TText(2000., 0.1, "model3, nP {}, chi2: {}".format(orderBkg[2],chi2s_bkg3))
                     txt3.SetTextSize(0.04) 
                     txt3.SetTextColor(ROOT.kRed-4) 
                     xframe.addObject(txt3) 
                     txt4 = ROOT.TText(2000., 0.01, "model4, nP {}, chi2: {}".format(orderBkg[3],chi2s_bkg4))
                     txt4.SetTextSize(0.04) 
                     txt4.SetTextColor(ROOT.kGreen+1) 
                     xframe.addObject(txt4) 
                     txt1.Draw()
                     txt2.Draw()
                     txt3.Draw()
                     txt4.Draw()

                     c1.SetLogy()
                     c1.SaveAs("TestAfterFit_"+ch+".pdf")

                     RSS[1] = getRSS(sig, ch, mT, modelBkg1, obsData,  [fitRes1], carddir, nDataEvts)
                     RSS[2] = getRSS(sig, ch, mT, modelBkg2, obsData,  [fitRes2], carddir, nDataEvts)
                     RSS[3] = getRSS(sig, ch, mT, modelBkg3, obsData,  [fitRes3], carddir, nDataEvts)
                     RSS[4] = getRSS(sig, ch, mT, modelBkg4, obsData,  [fitRes4], carddir, nDataEvts)

                     print RSS[2]

                     #**********************************************************
                     #                    ALTERNATIVE MODEL                    *
                     #**********************************************************
                     if bias:

                            # alternative function is Silvio's nominal function, but with +1 instead of -1
                            normAlt = RooRealVar("Bkg_"+ch+"alt_norm", "Number of background events", nBkgEvts, 0., 2.e4)
                            normData = RooRealVar("Data_"+ch+"alt_norm", "Number of background events", nDataEvts, 0., 2.e4) 
                            p1_1_alt = RooRealVar(ch_red + "_p1_1_alt", "p1", 1., 0., 50.)
                            p1_2_alt = RooRealVar(ch_red + "_p1_2_alt", "p1", 1., 0., 50.)
                            if(ch == "lowSVJ2_2017" or ch == "highSVJ2_2017"):
                                   #p1_3_alt = RooRealVar(ch_red + "_p1_3_alt", "p1", 1., 0., 10.)
                                   p1_3_alt = RooRealVar(ch_red + "_p1_3_alt", "p1", 1., -1000., 1000.)
                            else:
                                   p1_3_alt = RooRealVar(ch_red + "_p1_3_alt", "p1", 1., -1000., 1000.)
                            p1_4_alt = RooRealVar(ch_red + "_p1_4_alt", "p1", 1., 0., 1000.)
                            p2_2_alt = RooRealVar(ch_red + "_p2_2_alt", "p2", 1., 0., 20.)
                            if(ch == "highSVJ2_2016" or ch == "highSVJ2_2017" or ch == "highSVJ2_2018" or ch == "lowSVJ2_2017" or ch == "highSVJ"):
                                   #p2_3_alt = RooRealVar(ch_red + "_p2_3_alt", "p2", 1., 0., 10.)
                                   p2_3_alt = RooRealVar(ch_red + "_p2_3_alt", "p2", 1., -1000., 1000.)
                            else:
                                   p2_3_alt = RooRealVar(ch_red + "_p2_3_alt", "p2", 1., -1000., 1000.)
                            p2_4_alt = RooRealVar(ch_red + "_p2_4_alt", "p2", 1., 0., 1000.)
                            if(ch == "highSVJ2_2016" or ch == "highSVJ2_2017" or ch == "highSVJ2_2018" or ch == "lowSVJ2_2017" or ch == "highSVJ2"):
                                   #p3_3_alt = RooRealVar(ch_red + "_p3_3_alt", "p3", 1., 0., 10.)
                                   p3_3_alt = RooRealVar(ch_red + "_p3_3_alt", "p3", 1., -1000., 1000.)
                            else:
                                   p3_3_alt = RooRealVar(ch_red + "_p3_3_alt", "p3", 1., 0., 1000.)
                            p3_4_alt = RooRealVar(ch_red + "_p3_4_alt", "p3", 1., 0., 1000.)
                            p4_4_alt = RooRealVar(ch_red + "_p4_4_alt", "p4", 1., 0., 1000.) 



                            # Dijet Function Alt
                            modelAlt1_rgp = RooGenericPdf(modelAltName+"1_rgp", "Dij. fit (1 par.)", "pow(1 - @0/13000, abs(@1)) ", RooArgList(mT, p1_1_alt))
                            modelAlt2_rgp = RooGenericPdf(modelAltName+"2_rgp", "Dij. fit (2 par.)", "pow(1 - @0/13000, abs(@1)) * pow(@0/13000, -abs(@2))", RooArgList(mT, p1_2_alt, p2_2_alt))
                            modelAlt3_rgp = RooGenericPdf(modelAltName+"3_rgp", "Dij. fit (3 par.)", "pow(1 - @0/13000, abs(@1)) * pow(@0/13000, -abs(@2)-abs(@3)*log(@0/13000))", RooArgList(mT, p1_3_alt, p2_3_alt, p3_3_alt))
                            modelAlt4_rgp = RooGenericPdf(modelAltName+"4_rgp", "Dij. fit (4 par.)", "pow(1 - @0/13000, abs(@1)) * pow(@0/13000, -abs(@2)-log(@0/13000)*(abs(@3) + abs(@4)* log(@0/13000)))", RooArgList(mT, p1_4_alt, p2_4_alt, p3_4_alt, p4_4_alt))

                            # New Alt Function
                            #modelAlt1_rgp = RooGenericPdf(modelAltName+"1_rgp", "Alt. fit (1 par.)", "exp(@1*(@0/13000))", RooArgList(mT, p1_1_alt))
                            #modelAlt2_rgp = RooGenericPdf(modelAltName+"2_rgp", "Alt. fit (2 par.)", "exp(@1*(@0/13000)+@2*log(@0/13000))", RooArgList(mT, p1_2_alt, p2_2_alt))
                            #modelAlt3_rgp = RooGenericPdf(modelAltName+"3_rgp", "Alt. fit (3 par.)", "exp(@1*(@0/13000)+@2*log(@0/13000)+@3*pow(log(@0/13000),2))", RooArgList(mT, p1_3_alt, p2_3_alt, p3_3_alt))
                            #modelAlt4_rgp = RooGenericPdf(modelAltName+"4_rgp", "Alt. fit (4 par.)", "exp(@1*(@0/13000)+@2*log(@0/13000)+@3*pow(log(@0/13000),2)+@4*pow(log(@0/13000),3))", RooArgList(mT, p1_4_alt, p2_4_alt, p3_4_alt, p4_4_alt))

                            modelAlt1 = RooParametricShapeBinPdf(modelAltName+"1", "Alt. Fit 1par", modelAlt1_rgp, mT, RooArgList(p1_1_alt), histBkgData)
                            modelAlt2 = RooParametricShapeBinPdf(modelAltName+"2", "Alt. Fit 2par", modelAlt2_rgp, mT, RooArgList(p1_2_alt, p2_2_alt), histBkgData)
                            modelAlt3 = RooParametricShapeBinPdf(modelAltName+"3", "Alt. Fit 3par", modelAlt3_rgp, mT, RooArgList(p1_3_alt, p2_3_alt, p3_3_alt), histBkgData)
                            modelAlt4 = RooParametricShapeBinPdf(modelAltName+"4", "Alt. Fit 4par", modelAlt4_rgp, mT, RooArgList(p1_4_alt, p2_4_alt, p3_4_alt, p4_4_alt), histBkgData)

                            RSS_alt = {}
                            fitrange = "Full"
                            #if (ch == "highSVJ1_2016"): fitrange = "Low,High"
                            fitRes1_alt = modelAlt1.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                            fitRes2_alt = modelAlt2.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                            fitRes3_alt = modelAlt3.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))
                            fitRes4_alt = modelAlt4.fitTo(obsData, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(2), RooFit.Range(fitrange))

                            RSS_alt[1] = getRSS(sig, ch, mT, modelAlt1, obsData,  [fitRes1_alt], carddir,  nDataEvts, label = "alt")
                            RSS_alt[2] = getRSS(sig, ch, mT, modelAlt2, obsData,  [fitRes2_alt], carddir,  nDataEvts, label = "alt")
                            RSS_alt[3] = getRSS(sig, ch, mT, modelAlt3, obsData,  [fitRes3_alt], carddir,  nDataEvts, label = "alt")
                            RSS_alt[4] = getRSS(sig, ch, mT, modelAlt4, obsData,  [fitRes4_alt], carddir,  nDataEvts, label = "alt")
                            length = 1
                            if(length<2):
                                   xframeAlt = mT.frame(ROOT.RooFit.Title("extended ML fit example"))

                                   c2 = ROOT.TCanvas()
                                   c2.cd()
                                   obsData.plotOn(xframeAlt, RooFit.Name("obsData"))

                                   modelAlt1.plotOn(xframeAlt, RooFit.Name("modelAlt1"), ROOT.RooFit.LineColor(ROOT.kPink))
                                   modelAlt2.plotOn(xframeAlt, RooFit.Name("modelAlt2"), ROOT.RooFit.LineColor(ROOT.kBlue))
                                   modelAlt3.plotOn(xframeAlt, RooFit.Name("modelAlt3"), ROOT.RooFit.LineColor(ROOT.kRed))
                                   modelAlt4.plotOn(xframeAlt, RooFit.Name("modelAlt4"), ROOT.RooFit.LineColor(ROOT.kGreen))

                                   chi2s_alt1 = xframeAlt.chiSquare("modelAlt1", "obsData",1)
                                   chi2s_alt2 = xframeAlt.chiSquare("modelAlt2", "obsData",2)
                                   chi2s_alt3 = xframeAlt.chiSquare("modelAlt3", "obsData",3)
                                   chi2s_alt4 = xframeAlt.chiSquare("modelAlt4", "obsData",4)

                                   xframeAlt.SetMinimum(0.002)
                                   xframeAlt.Draw()

                                   c2.SetLogy()
                                   c2.SaveAs("TestAfterFit_"+ch+"_Alt.pdf")



                     #*******************************************************#
                     #                                                       #
                     #                         Fisher                        #
                     #                                                       #
                     #*******************************************************#


                     ofile = open("FisherTest_%s.txt"%(ch),"w")
                     report = "Results from Fisher Test for category %s \n" % (ch)
                     #{"chiSquared":roochi2,"chi2" : chi2, "chi1" : chi1, "rss" : rss, "res" : res, "nbins" : hist.GetN(), "npar" : npar}
                     report += "func\tchi2\trss\tnBins\tnPar \n"
                     for i in range(1,len(RSS)+1):
                            report += "   %i\t%.2f\t%.2f\t%i\t%i\n" % (i, RSS[i]['chiSquared'],RSS[i]['rss'],RSS[i]['nbins'],RSS[i]['npar'])
                     report += "*******************************************************\n"
                     report += "fTest = ((RSS1-RSS2)/(npar2-npar1))/(RSS2/(nBins-npar2))\n"
                     report += "*******************************************************\n"


                     order = 0
                     print "-"*25
                     print "function & $\\chi^2$ & RSS & ndof & F-test & result \\\\"
                     print "\\hline"
                     # BEGIN NEW WAY, tests in a 'single elimination' bracket
                     # begin with 1v2, test the `winner' vs 3, then test that winner vs 4. 
                     # which ever function wins the vs 4 test, is the function we use.
                     aCrit = 0.05
                     RTDict = {"1v2":-1,"1v3":-1,"1v4":-1,"2v3":-1,"2v4":-1,"3v4":-1}
                     FDict = {"1v2":-1,"1v3":-1,"1v4":-1,"2v3":-1,"2v4":-1,"3v4":-1}
                     for o1 in xrange(1, len(RSS)):
                            for o2 in xrange(o1+1,len(RSS)+1):
                                   RTDict[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)] = fisherTest(RSS[o1]['rss'], RSS[o2]['rss'], RSS[o1]['npar'], RSS[o2]['npar'], RSS[o2]["nbins"])
                                   report += "%d par vs %d par: CL=%.5f F_t=%.5f\n" % (RSS[o1]['npar'], RSS[o2]['npar'], RTDict[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)])

                     if RTDict["1v2"] > aCrit:
                            if RTDict["1v3"] > aCrit:
                                   if RTDict["1v4"] > aCrit:
                                          order = 1
                                   else:
                                          order = 4
                            else:
                                   if RTDict["3v4"] > aCrit:
                                          order = 3
                                   else:
                                          order = 4
                     else:
                            if RTDict["2v3"] > aCrit:
                                   if RTDict["2v4"] > aCrit:
                                          order = 2
                                   else:
                                          order = 4
                            else:
                                   if RTDict["3v4"] > aCrit:
                                          order = 3
                                   else:
                                          order = 4
                     #OLD WAY OF DOING IT, only does n vs n+1
                     #for o1 in xrange(1, len(RSS)):
                     #       o2 = o1+1
                     #
                     #       CL = fisherTest(RSS[o1]['rss'], RSS[o2]['rss'], RSS[o1]['npar'], RSS[o2]['npar'], RSS[o2]["nbins"])
                     #
                     #       print "%d par vs %d par CL=%.2f  \n" % (RSS[o1]['npar'], RSS[o2]['npar'], CL),
                     #       report += "%d par vs %d par: CL=%.2f \n" % (RSS[o1]['npar'], RSS[o2]['npar'], CL)
                     #       if CL > 0.10: # The function with less parameters is enough
                     #              if not order:
                     #                     order = o1
                     #                     print "%d par are sufficient" % (RSS[o1]['npar'])
                     #                     #report += "%d par are sufficient \n" % (o1)
                     #       else:
                     #              print "%d par are needed" % (RSS[o2]['npar']),
                     #
                     #       print "\\\\"
                     moreParFlag = 0
                     if order == 0:
                            order = 4
                            moreParFlag = 1
                     print "\\hline"
                     print "-"*25   
                     print "Order is", order, "("+ch+")"
                     report += "Order is %d (%s)\n" % (order, ch)
                     report += ("1 param: " + ", ".join(['%.2f']*len(RSS[1]["parVals"])) + "\n") % tuple(RSS[1]["parVals"])
                     report += ("2 param: " + ", ".join(['%.2f']*len(RSS[2]["parVals"])) + "\n") % tuple(RSS[2]["parVals"])
                     report += ("3 param: " + ", ".join(['%.2f']*len(RSS[3]["parVals"])) + "\n") % tuple(RSS[3]["parVals"])
                     report += ("4 param: " + ", ".join(['%.2f']*len(RSS[4]["parVals"])) + "\n") % tuple(RSS[4]["parVals"])
                     report += "%d par are sufficient\n" % (RSS[order]['npar'])
                     if moreParFlag: report += "BUT really, more Pars are needed!"
                     for i in range(1,len(RSS)+1):
                            report += "model%i chi2: %.2f\n" % (i, RSS[i]['chiSquared'])

                     ofile.write(report)
                     ofile.close()

                     if bias: 
                            print "Running in BIAS mode"

                            ofile_alt = open("FisherTest_alt_%s.txt"%(ch),"w")
                            report = "Results from Fisher Test on the alternative function for category %s \n" % (ch)
                            report += "func\tchi2\trss\tnBins\tnPar \n"
                            for i in range(1,len(RSS)+1):
                                   report += "   %i\t%.2f\t%.2f\t%i\t%i\n" % (i, RSS_alt[i]['chiSquared'],RSS_alt[i]['rss'],RSS_alt[i]['nbins'],RSS_alt[i]['npar'])
                            report += "*******************************************************\n"
                            report += "fTest = ((RSS1-RSS2)/(npar2-npar1))/(RSS2/(nBins-npar2))\n"
                            report += "*******************************************************\n"


                            order_alt = 0
                            print "-"*25
                            print "function & $\\chi^2$ & RSS & ndof & F-test & result \\\\"
                            print "\\hline"
                            # BEGIN NEW WAY, tests in a 'single elimination' bracket
                            # begin with 1v2, test the `winner' vs 3, then test that winner vs 4. 
                            # which ever function wins the vs 4 test, is the function we use.
                            RTDict_alt = {"1v2":-1,"1v3":-1,"1v4":-1,"2v3":-1,"2v4":-1,"3v4":-1}
                            FDict = {"1v2":-1,"1v3":-1,"1v4":-1,"2v3":-1,"2v4":-1,"3v4":-1}
                            for o1 in xrange(1, len(RSS_alt)):
                                   for o2 in xrange(o1+1,len(RSS_alt)+1):
                                          RTDict_alt[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)] = fisherTest(RSS_alt[o1]['rss'], RSS_alt[o2]['rss'], RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], RSS_alt[o2]["nbins"])
                                          report += "%d par vs %d par: CL=%.5f F_t=%.5f\n" % (RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], RTDict_alt[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)])

                            if RTDict_alt["1v2"] > aCrit:
                                   if RTDict_alt["1v3"] > aCrit:
                                          if RTDict_alt["1v4"] > aCrit:
                                                 order_alt = 1
                                          else:
                                                 order_alt = 4
                                   else:
                                          if RTDict_alt["3v4"] > aCrit:
                                                 order_alt = 3
                                          else:
                                                 order_alt = 4
                            else:
                                   if RTDict_alt["2v3"] > aCrit:
                                          if RTDict_alt["2v4"] > aCrit:
                                                 order_alt = 2
                                          else:
                                                 order_alt = 4
                                   else:
                                          if RTDict_alt["3v4"] > aCrit:
                                                 order_alt = 3
                                          else:
                                                 order_alt = 4

                            #OLD WAY OF DOING IT, only does n vs n+1
                            #for o1 in xrange(1, len(RSS)):
                            #       #o2 = o1+1
                            #       for o2 in xrange(o1+1,len(RSS)+1):
                            #
                            #              CL = fisherTest(RSS_alt[o1]['rss'], RSS_alt[o2]['rss'], RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], RSS_alt[o2]["nbins"])
                            #
                            #              print "%d par vs %d par CL=%.5f  \n" % (RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], CL),
                            #              report += "%d par vs %d par: CL=%.5f \n" % (RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], CL)
                            #              if CL > 0.1: # The function with less parameters is enough
                            #                     if not order_alt:
                            #                            order_alt = o1
                            #                            print "%d par are sufficient" % (RSS_alt[o1]['npar'])
                            #                            #report += "%d par are sufficient \n" % (o1)
                            #              else:
                            #                     print "%d par are needed" % (RSS_alt[o2]['npar']),
                            #
                            #              print "\\\\"
                            #END OLD WAY OF DOING IT
                            moreParFlag_alt = 0
                            if order_alt == 0:
                                   order_alt = 4
                                   moreParFlag_alt = 1
                            print "\\hline"
                            print "-"*25   
                            print "Order is", order_alt, "("+ch+")"
                            report += "Order is %d (%s)\n" % (order_alt, ch)
                            report += ("alt 1: " + ", ".join(['%.2f']*len(RSS_alt[1]["parVals"])) + "\n") % tuple(RSS_alt[1]["parVals"])
                            report += ("alt 2: " + ", ".join(['%.2f']*len(RSS_alt[2]["parVals"])) + "\n") % tuple(RSS_alt[2]["parVals"])
                            report += ("alt 3: " + ", ".join(['%.2f']*len(RSS_alt[3]["parVals"])) + "\n") % tuple(RSS_alt[3]["parVals"])
                            report += ("alt 4: " + ", ".join(['%.2f']*len(RSS_alt[4]["parVals"])) + "\n") % tuple(RSS_alt[4]["parVals"])
                            report += "%d par are sufficient\n" % (RSS_alt[order_alt]['npar'])
                            if moreParFlag_alt: report += "BUT really, more Pars are needed!"
                            for i in range(1,min(len(RSS),len(RSS_alt))+1):
                                   report += "model%i chi2, chi2_alt: %.2f %.2f\n" % (i, RSS[i]['chiSquared'],RSS_alt[i]['chiSquared'])

                            ofile_alt.write(report)
                            ofile_alt.close()
                     #temporarily change to force order 2
                     #order = 3
                     if order==1:
                            modelBkg = modelBkg1#.Clone("Bkg")
                            #normzBkg = normzBkg2#.Clone("Bkg_norm")
                            fitRes = fitRes1
                     elif order==2:
                            modelBkg = modelBkg2#.Clone("Bkg")
                            #normzBkg = normzBkg2#.Clone("Bkg_norm")
                            fitRes = fitRes2
                     elif order==3:
                            modelBkg = modelBkg3#.Clone("Bkg")
                            #normzBkg = normzBkg3#.Clone("Bkg_norm")
                            fitRes = fitRes3
                     elif order==4:
                            modelBkg = modelBkg4#.Clone("Bkg")
                            #normzBkg = normzBkg4#.Clone("Bkg_norm")
                            fitRes = fitRes4
                     else:
                            print "Functions with", order+1, "or more parameters are needed to fit the background"
                            exit()

                     modelBkg.SetName(modelName)
                     normData.SetName(modelName+"_norm")
                     wsfilename = "ws.root"
                     wfile_ = ROOT.TFile(wsfilename)
                     w_ = wfile_.FindObjectAny("BackgroundWS")
                     if not w_.__nonzero__() :  w_ = RooWorkspace("BackgroundWS", "workspace")
                     else:  w_ =  wfile_.Get("BackgroundWS")
                     print "Storing ", modelName
                     getattr(w_, "import")(modelBkg, RooFit.Rename(modelBkg.GetName()))

                     getattr(w_, "import")(obsData, RooFit.Rename("data_obs"))


                     if bias:
                            #temporarily change to force order 2
                            #order_alt = 3
                            if order_alt==1:
                                   print "Alt Model is modelAlt1"
                                   modelAlt = modelAlt1#.Clone("Alt")
                                   fitRes = fitRes1_alt
                            elif order_alt==2:
                                   print "Alt Model is modelAlt2"
                                   modelAlt = modelAlt2#.Clone("Alt")
                                   fitRes = fitRes2_alt
                            elif order_alt==3:
                                   print "Alt Model is modelAlt3"
                                   modelAlt = modelAlt3#.Clone("Alt")
                                   fitRes = fitRes3_alt
                            elif order_alt==4:
                                   print "Alt Model is modelAlt4"
                                   modelAlt = modelAlt4#.Clone("Alt")
                                   fitRes = fitRes4_alt
                            else:
                                   print "Functions with", order_alt+1, "or more parameters are needed to fit the background"
                                   exit()
                            modelAlt.SetName(modelAltName)
                            normAlt.SetName(modelAltName+"_norm")
                            getattr(w_, "import")(modelAlt, RooFit.Rename(modelAlt.GetName()))

                     wstatus = w_.writeToFile(wsfilename, False)







              wfile = ROOT.TFile("ws.root")
              wBkg =  wfile.Get("BackgroundWS")
              print "workspace ", wBkg.Print()
              print "Model name: ", modelName
              modelBkg = wBkg.pdf(modelName)
              modelAlt = wBkg.pdf(modelAltName)
              normAlt = wBkg.pdf(modelAltName+"_norm")
              normData = wBkg.pdf(modelName+"_norm")
              parSet = modelBkg.getParameters(obsData)
              parSet.Print()
              argList = ROOT.RooArgList(parSet)
              parNames = [ argList[i].GetName() for i in xrange(0, len(parSet))]                                 


             # modelExt = RooExtendPdf(modelBkg.GetName(), modelBkg.GetTitle(), modelBkg, normzBkg)

              for i in xrange(0, len(parSet)):
                    print  argList[i].GetName(), "   ", argList[i].getVal()


              # BIAS STUDY
              # Generate pseudo data

              #setToys = RooDataSet()
              #setToys.SetName("data_toys")
              #setToys.SetTitle("Data (toys)")

              #if not isData:
              #       print " - Generating", nBkgEvts, "events for toy data"
              #       if bias: setToys = modelAlt.generateBinned(RooArgSet(mT), nBkgEvts, 0,0)
                     #fitRes = modelBkg.fitTo(setToys, RooFit.Extended(True), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(1))

              if bias:
                     ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
                     from ROOT import RooMultiPdf
                     pdf_index_string = "pdf_index_%s" % (ch)
                     cat = RooCategory(pdf_index_string, "Index of Pdf which is active")
                     pdfs = RooArgList()
                     pdfs.add(modelBkg)
                     pdfs.add(modelAlt)
                     roomultipdf = RooMultiPdf("roomultipdf", "All Pdfs", cat, pdfs)
                     normulti = RooRealVar("roomultipdf_norm", "Number of background events", nDataEvts, 0., 1.e6)

              # create workspace
              w = RooWorkspace("SVJ", "workspace")
              # Dataset
              # ATT: include isData

              #getattr(w, "import")(obsData, RooFit.Rename("Bkg"))
              #if isData: getattr(w, "import")(obsData, RooFit.Rename("data_obs"))
              getattr(w, "import")(obsData, RooFit.Rename("data_obs"))
              #else: getattr(w, "import")(setToys, RooFit.Rename("data_obs"))
              #getattr(w, "import")(obsData, RooFit.Rename("data_obs"))
              getattr(w, "import")(sigData, RooFit.Rename(sig))
              
              if bias:
                     getattr(w, "import")(cat)#, RooFit.Rename(cat.GetName()))
                     getattr(w, "import")(normulti)#, RooFit.Rename(normulti.GetName()))
                     getattr(w, "import")(roomultipdf)#, RooFit.Rename(roomultipdf.GetName()))

              for i in xrange(hist.GetNbinsX()):
                     mcstatSysName = "mcstat_%s_%s_MC2018bin%d"  % (ch, sig, i+1)
                     #print mcstatSysName
                     #print sig + "_" + mcstatSysName + "Up"
                     mcstatSigUp = getHist(ch, sig + "_" + mcstatSysName + "Up", ifile)

                     #print "Integral  ", mcstatSigUp.Integral()
                     mcstatSigDown = getHist(ch, sig + "_" + mcstatSysName + "Down", ifile)
                     mcstatSigHistUp = RooDataHist(sig + "_" + mcstatSysName + "Up", "Data (MC sig)",  RooArgList(mT), mcstatSigUp, 1.)
                     mcstatSigHistDown = RooDataHist(sig + "_" + mcstatSysName + "Down", "Data (MC sig)",  RooArgList(mT), mcstatSigDown, 1.)
                     getattr(w, "import")(mcstatSigHistUp, RooFit.Rename(sig + "_" + mcstatSysName + "Up") )
                     getattr(w, "import")(mcstatSigHistDown, RooFit.Rename(sig + "_" + mcstatSysName + "Down") )

              
              for sysName,sysValue  in syst.iteritems():
                     if(sysValue[0]=="shape" and "mcstat" not in sysName):              
                            sysUp =  getHist(ch, sig + "_" + sysName + "Up", ifile)
                            sysDown =  getHist(ch, sig + "_" + sysName + "Down", ifile)
                            #print "==> Trigg sys name: ", sig + "_" + sysName + "Down"
                            sysSigHistUp = RooDataHist(sig + "_" + sysName + "Up", sysName + " uncertainty",  RooArgList(mT), sysUp, 1.)
                            sysSigHistDown = RooDataHist(sig + "_" + sysName + "Down", sysName + " uncertainty",  RooArgList(mT), sysDown, 1.)
                            getattr(w, "import")(sysSigHistUp, RooFit.Rename(sig + "_" + sysName + "Up") )
                            getattr(w, "import")(sysSigHistDown, RooFit.Rename(sig + "_" + sysName + "Down") )
                            


              #else: getattr(w, "import")(setToys, RooFit.Rename("data_obs"))
              getattr(w, "import")(modelBkg, (modelBkg.GetName()))
              if bias:
                     getattr(w, "import")(modelAlt, (modelAltName))
              #getattr(w, "import")(normzBkg, RooFit.Rename(normzBkg.GetName()))
              wstatus = w.writeToFile("ws_%s_%s_%s.root" % (sig, ch, mode), True)

              if wstatus == False : print "Workspace", "ws_%s_%s_%s.root" % (sig, ch, mode) , "saved successfully"
              else: print "Workspace", "ws_%s_%s_%s.root" % (sig, ch, mode) , "not saved successfully"
              workfile = "ws_%s_%s_%s.root" % ( sig, ch, mode)
              # ======   END MODEL GENERATION   ======       
              if(verbose):w.Print()


       rates = {}
       procLine = ""
       procNumbLine = ""
       rateLine = ""
       binString = ""

       if(mode == "template"):       

              processes.append("Bkg")
              processes[:-1] = []
              rates["Bkg"] = nBkgEvts
              procLine += ("%-43s") % (modelBkg.GetName())
              rateLine += ("%-43s") % (rates["Bkg"])
              binString += (("%-43s") % (ch) ) * (2)
              procNumbLine = 1 
              
       else:
              i = 1
              bkgrate = 0
              for p in processes:
                     rates[p] = getRate(ch, p, ifile)
                     bkgrate =  rates[p]
                     procNumbLine += ("%-43s") % (i)
                     procLine += ("%-43s") % (p)
                     rateLine += ("%-43.1f") % (bkgrate)
                     i+=1
              binString += (("%-43s") % (ch) ) * (len(processes)+1)


       rates["data_obs"] = getRate(ch, "data_obs", ifile)
       print "TEST bkg/obs: ", rates["data_obs"], nDataEvts, binMax
       rates[sig] = getRate(ch, sig, ifile)
       print "TEST sig: ", rates[sig], nSigEvts, binMax



       card  = "imax 1 number of channels \n"
       card += "jmax * number of backgrounds \n"
       card += "kmax * number of nuisance parameters\n"
       card += "-----------------------------------------------------------------------------------\n"

       if(mode == "template"):
              #              card += "shapes   %s  %s    %s    %s    %s\n" % (sig, ch, ifilename, "$CHANNEL/$PROCESS", "$CHANNEL/$PROCESS_SYSTEMATIC")
              #              card += "shapes            %-15s  %-5s    %s%s.root    %s\n" % (sig, ch, WORKDIR, ch, "SVJ:$PROCESS")
              card += "shapes   %s  %s    %s    %s\n" % (modelBkg.GetName(), ch, workfile, "SVJ:$PROCESS")
              card += "shapes   %s  %s    %s    %s    %s\n" % (sig, ch, workfile, "SVJ:$PROCESS", "SVJ:$PROCESS_$SYSTEMATIC")
              card += "shapes   %s  %s    %s    %s\n" % ("data_obs", ch, workfile, "SVJ:$PROCESS")

       else:  
              card += "shapes   *      *   %s    %s    %s\n" % (hist_filename, "$CHANNEL/$PROCESS", "$CHANNEL/$PROCESS_$SYSTEMATIC")
              card += "shapes   data_obs      *   %s    %s\n" % (hist_filename, "$CHANNEL/$PROCESS")
       card += "-----------------------------------------------------------------------------------\n"
       card += "bin               %s\n" % ch
       print "===> Observed data: ", rates["data_obs"]
       card += "observation       %.6f\n" % (rates["data_obs"]) 
       card += "-----------------------------------------------------------------------------------\n"
       card += "bin                                     %-43s\n" % (binString)
       card += "process                                 %-43s%-43s\n" % (sig, procLine) #"roomultipdf"
       card += "process                                 %-43s%-43s\n" % ("0", procNumbLine)
       card += "rate                                    %-43.6f%-43s\n" % (rates[sig], rateLine) #signalYield[m].getVal(), nevents
       card += "-----------------------------------------------------------------------------------\n"



       for sysName,sysValue  in syst.iteritems():

              if(sysValue[0]=="lnN"): 
                     card += "%-20s%-20s" % (sysName, sysValue[0])
                     if(sysValue[1]=="all"and len(sysValue)>2):
                            if(mode == "template"): card += "%-20s" % (sysValue[2]) * (2)
                            else: card += "%-20s" % (sysValue[2]) * (len(processes) + 1)

                     else:
                            if (sysValue[1]=="all"):
                                   sysValue[1] = copy.deepcopy(processes)
                                   sysValue[1].append(sig)

                            hsysName =  "_" + sysName  
                            hsysNameUp = "_" + sysName + "UP"  
                            hsysNameDown = "_" + sysName + "DOWN" 
                                                        #print "Applying syst on ", sysValue[1]
                            if("sig" in sysValue[1]):
                                   if(getRate(ch, sig, ifile) != 0.): sigSys = abs((getRate(ch, sig+hsysNameUp, ifile) - getRate(ch, sig+hsysNameDown, ifile))/ (2* getRate(ch, sig, ifile)))
                                   else: sigSys = 1
                                   if(sigSys<1.and sigSys >0.): sigSys = sigSys + 1
                                   card += "%-20s" % (sigSys)
                            else:  card += "%-20s" % ("-")
                            for p in processes:
                                   if (p in sysValue[1]):
                                          if (getRate(ch, p, ifile) != 0.): bkgSys = abs((getRate(ch, p+hsysNameUp, ifile) - getRate(ch, p+hsysNameDown, ifile))/ (2* getRate(ch, p, ifile)) )
                                          else: bkgSys = 1
                                          if(bkgSys<1.and bkgSys >0.): bkgSys = bkgSys + 1
                                          card += "%-20s" % (bkgSys)
                                   else:  card += "%-20s" % ("-")



              elif(sysValue[0]=="lnU"):
                     card += "%-20s%-20s%-20s%-20s" % (sysName, sysValue[0], "-", sysValue[2])
              elif(sysValue[0]=="shape"):
                  if("mcstat" not in sysName):
                            card += "%-20s     shape     " % (sysName)
                            if ("sig" in sysValue[1]): card += "%-20s" % ( "1") 
                            else: card += "%-20s" % ( "-") 
                            for p in processes:
                                   if (p in sysValue[1]): card += "%-20s" % ( "1") 
                                   else: card += "%-20s" % ( "-") 
              

                  else:
                            # CAMBIARE NOME DELLA SYST      , change name of syst?               
                            for samp in sysValue[1]:
                                   sampName = ""
                                   line = ""
                                   if (samp == "sig" or samp == "Sig"): 
            
                                          line = "%-20s" % ( "1") 
                                          line += "%-20s" % ("-") * (len(processes)) 
                                          sampName = sig
                                   elif(mode != "template"):
                                          
                                          line = "%-20s" % ( "-") 
                                          lineProc = ["%-20s" % ( "-") for x in xrange (len(processes))]
                                          if samp in processes: 
                                                 index = processes.index(samp)  
                                                 lineProc[index] = "1"
                                                 
                                          lineProc = "         ".join(lineProc)
                                          line += lineProc
                                          sampName =  samp

                                   else: continue


                                   for i in xrange(hist.GetNbinsX()):

                                          sysName = "mcstat_%s_%s_MC2018bin%d      "  % (ch, sampName, i+1)
                                          card += "%-20s   shape   " % (sysName)
                                          card += line
                                          card += "\n"        

              card += "\n"
       

       for par in parNames: card += "%-20s%-20s\n" % (par, "flatParam")

       #card += "SF            extArg     1 [0.75,1.25]\n"
       #card += "SF param 1 0.25\n"
       #for par, formula in rateParams.iteritems(): 
       #       formula = formula.replace("%s", '%.3f' % (eff))
       #       card += "%-20s %-20s %-20s %-20s %s %-20s\n" % (par + "_rate", "rateParam", ch, sig, formula, "SF")


      # card += "\n"

       outname =  "%s_%s_%s.txt" % (sig, ch, mode)
       cardfile = open(outname, 'w')
       cardfile.write(card)
       cardfile.close()

       if bias:

              card = card.replace(modelBkg.GetName(), "roomultipdf")
              card.replace("rate                                    %-20.6f%-20.6f\n" % (1, 1), "rate                                    %-20.6f%-20.6f\n" % (10, 1))
              card += "%-35s     discrete\n" % (pdf_index_string)
              outname = "%s_%s_%s_bias.txt" % (sig, ch, mode)
              cardfile = open(outname, 'w')
              cardfile.write(card)
              cardfile.close()
    

       #print card
       return card




