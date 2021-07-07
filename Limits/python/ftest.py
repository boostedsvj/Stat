import ROOT
import json

#changed to no longer need the settings.py file
# files now saved to base directory and exported to EOS after jobs completion

#ftest.py - part I of separating ftests from datacard writing
# this part will create a ws.root file containing the two histograms and two funcs as well as funcs' parameters


from ROOT import RooRealVar, RooDataHist, RooArgList, RooGenericPdf, RooBernstein, RooExtendPdf, RooCmdArg, RooWorkspace, RooFit, RooDataSet, RooArgSet, RooCategory, RooFitResult, RooCurve, RooParametricShapeBinPdf
import os, sys
from array import array
import copy, math, pickle
import collections
import numpy as np
from numpy import ndarray
from bruteForce import silence, remakePdf
silence()

ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetPadColor(0)

#*******************************************************#
#                                                       #
#                     Utility Functions                 #
#                                                       #
#*******************************************************#

isData = True

def getHist(ch, process, ifile):
       hName = ch + "/"+ process
       #print "Histo Name ", hName
       h = ifile.Get(hName)
       h.SetDirectory(0)
       return h

def altMerge(l1, l2):
	result = [None]*(len(l1)+len(l2))
	result[::2] = l1
	result[1::2] = l2
	return result


def getRSS(ch, variable, model, dataset, fitRes, carddir, doplots, norm = -1, label = "nom"):
       name = model.GetName()
       order = int(name[-1])
       npar = fitRes[0].floatParsFinal().getSize() if len(fitRes)>0 else 0
       varArg = ROOT.RooArgSet(variable)
      
       frame = variable.frame(ROOT.RooFit.Title(""))
       dataset.plotOn(frame, RooFit.Invisible())
       #print(fitRes[0])
       if doplots and len(fitRes) > 0 and npar <= 3: graphFit = model.plotOn(frame, RooFit.VisualizeError(fitRes[0], 1, False), RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("FL"), RooFit.Range("Full"))

       model.plotOn(frame, RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("L"), RooFit.Name(model.GetName()),  RooFit.Range("Full"))
       model.paramOn(frame, RooFit.Label(model.GetTitle()), RooFit.Layout(0.45, 0.95, 0.94), RooFit.Format("NEAU"))
       
       graphData = dataset.plotOn(frame, RooFit.DataError(ROOT.RooAbsData.Poisson if isData else ROOT.RooAbsData.SumW2), RooFit.DrawOption("PE0"), RooFit.Name(dataset.GetName()))

       pulls = frame.pullHist(dataset.GetName(), model.GetName(), True)
       residuals = frame.residHist(dataset.GetName(), model.GetName(), False, True) # this is y_i - f(x_i)
    
       roochi2 = frame.chiSquare(model.GetName(), dataset.GetName(),npar)#dataset.GetName(), model.GetName()) #model.GetName(), dataset.GetName()
       # get ndf by hand
       dhist = frame.findObject(dataset.GetName(),ROOT.RooHist.Class())
       nbins = 0
       for i in range(dhist.GetN()):
            x = ROOT.Double(0.)
            y = ROOT.Double(0.)
            dhist.GetPoint(i,x,y)
            if y!=0: nbins += 1
       chi2 = roochi2 * ( nbins - npar)
       roopro = ROOT.TMath.Prob(chi2, nbins - npar)

       frame.SetMaximum(frame.GetMaximum()*10.)
       frame.SetMinimum(0.1)
       if doplots and npar <= 3:

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
              
              roochi2_small = format(roochi2, '.2f')
              txt = ROOT.TText(0.65, 0.7, "chiSquared: " + str(roochi2_small))              
              txt.SetNDC()
              txt.SetTextSize(0.04) 
              txt.SetTextColor(ROOT.kRed) 
              txt.Draw()

              roopro_small = format(roopro, '.2f')
              txt_2 = ROOT.TText(0.65, 0.65, "prob: " + str(roopro_small))              
              txt_2.SetNDC()
              txt_2.SetTextSize(0.04) 
              txt_2.SetTextColor(ROOT.kRed) 
              txt_2.Draw()
              
              c.cd()
              c.Update()
              c.cd()
              pad2 = ROOT.TPad("pad2", "pad2", 0, 0.1, 1, 0.34)
              ROOT.SetOwnership(pad2, False)
              pad2.SetTopMargin(0)
              pad2.SetBottomMargin(0.25)
              pad2.SetGridx()
              pad2.SetGridy()
              pad2.Draw()
              pad2.cd()
              pad2.Clear()
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
       
       if doplots and npar <= 3:
              c.Update()
              c.Range(xmin_, -3.5, xmax_, 3.5)
              line = ROOT.TLine(xmin_, 0., xmax_, 0.)
              line.SetLineColor(ROOT.kRed)
              line.SetLineWidth(2)
              line.Draw("same")
                            
              c.SaveAs("Residuals_"+ch+"_"+name+"_log.pdf")              

       # calculate RSS
       rss = 0.

       xmin, xmax = array('d', [0.]), array('d', [0.])
       dataset.getRange(variable, xmin, xmax)
       ratioHist = ROOT.TH1F("RatioPlot", "Ratio Plot", hist.GetN(), xmin_, xmax_)
       ROOT.SetOwnership(ratioHist, False)

       for i in xrange(0, hist.GetN()):
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
              
              rss += ry**2 

              ratioHist.SetBinContent(i+1, (ry - hy)/(-1*hy))
              ratioHist.SetBinError(i+1, (hey/ hy**2))
        
       rss = math.sqrt(rss)
       parValList = []
       parErrList = []
       #print(len(fitRes[0].floatParsFinal()))
       for iPar in range(len(fitRes[0].floatParsFinal())):
              #print(iPar)
              parValList.append((fitRes[0].floatParsFinal().at(iPar)).getValV())
              parErrList.append((fitRes[0].floatParsFinal().at(iPar)).getError())
       out = {"chiSquared":roochi2,"chi2" : chi2, "rss" : rss, "nbins" : nbins, "npar" : npar, "parVals": parValList, "parErr":parErrList}
       if doplots and npar <= 3:

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

              txt = ROOT.TText(0.65, 0.7, "chiSquared: " + str(roochi2_small))              
              txt.SetNDC()
              txt.SetTextSize(0.04) 
              txt.SetTextColor(ROOT.kRed) 
              txt.Draw()

              txt_2 = ROOT.TText(0.65, 0.65, "prob: " + str(roopro_small))              
              txt_2.SetNDC()
              txt_2.SetTextSize(0.04) 
              txt_2.SetTextColor(ROOT.kRed) 
              txt_2.Draw();

              c2.cd()
              c2.Update()
              c2.cd()
              pad2_2 = ROOT.TPad("pad2_2", "pad2_2", 0, 0.05, 1, 0.34)
              ROOT.SetOwnership(pad2_2, False)
              pad2_2.SetTopMargin(0);
              pad2_2.SetBottomMargin(0.25);
              pad2_2.SetGridx();
              pad2_2.SetGridy();
              pad2_2.Draw();
              pad2_2.cd()
              pad2_2.Clear()
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

              ratioHist.SetMarkerStyle(20)
              ratioHist.Draw("PE")
              line2 = ROOT.TLine(xmin_, 1., xmax_, 1.)
              line2.SetLineColor(ROOT.kRed)
              line2.SetLineWidth(2)
              line2.Draw("same")
              frame.SetTitle("")

              c2.SaveAs("Residuals_"+ch+"_"+name + "_ratio_log.pdf")

       return out

def drawTwoFuncs(ch, variable, modelA, modelB, dataset, fitRes, carddir,  norm = -1, label = "nom"):
       name = modelA.GetName()+"_"+modelB.GetName()
       varArg = ROOT.RooArgSet(variable)
      
       frame = variable.frame(ROOT.RooFit.Title(""))
       dataset.plotOn(frame, RooFit.Invisible())
       modelA.plotOn(frame, RooFit.VisualizeError(fitRes[0], 1, False), RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineStyle(2), RooFit.LineColor(ROOT.kBlue), RooFit.FillColor(ROOT.kCyan), RooFit.FillStyle(1001), RooFit.DrawOption("FL"), RooFit.Name(modelA.GetName()),  RooFit.Range("Full"))
       modelB.plotOn(frame, RooFit.VisualizeError(fitRes[1], 1, False), RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineStyle(2), RooFit.LineColor(ROOT.kRed),  RooFit.FillColor(ROOT.kOrange), RooFit.FillStyle(1001), RooFit.DrawOption("FL"), RooFit.Name(modelB.GetName()),  RooFit.Range("Full"))
       dataset.plotOn(frame, RooFit.DataError(ROOT.RooAbsData.Poisson if isData else ROOT.RooAbsData.SumW2), RooFit.DrawOption("PE0"), RooFit.Name(dataset.GetName()))
       modelA.plotOn(frame, RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kBlue), RooFit.DrawOption("L"), RooFit.Name(modelA.GetName()),  RooFit.Range("Full"))
       modelB.plotOn(frame, RooFit.Normalization(norm if norm>0 else dataset.sumEntries(), ROOT.RooAbsReal.NumEvent), RooFit.LineColor(ROOT.kRed), RooFit.DrawOption("L"), RooFit.Name(modelB.GetName()),  RooFit.Range("Full"))
       frame.SetMaximum(frame.GetMaximum()*10.)
       frame.SetMinimum(0.1)

       c = ROOT.TCanvas("c_"+ch+name, ch, 800, 800)
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
              
       txt = ROOT.TText(0.5, 0.7, modelA.GetName())              
       txt.SetNDC()
       txt.SetTextSize(0.04) 
       txt.SetTextColor(ROOT.kBlue) 
       txt.Draw()
       txt_2 = ROOT.TText(0.5, 0.8, modelB.GetName())              
       txt_2.SetNDC()
       txt_2.SetTextSize(0.04) 
       txt_2.SetTextColor(ROOT.kRed) 
       txt_2.Draw()

       c.cd()
       c.Update()
       c.cd()
       pad2 = ROOT.TPad("pad2", "pad2", 0, 0.1, 1, 0.3)
       ROOT.SetOwnership(pad2, False)
       pad2.SetTopMargin(0)
       pad2.SetBottomMargin(0.25)
       pad2.SetGridx()
       pad2.SetGridy()
       pad2.Draw()
       pad2.cd()
       pad2.Clear()

       frame_res = variable.frame(ROOT.RooFit.Title(""))
       frame_res.GetYaxis().SetRangeUser(-3.5, 3.5)
       frame_res.SetTitle("")
       frame_res.GetYaxis().SetTitle("Ratio Main/Alt")
       frame_res.GetYaxis().SetTitleOffset(.3)
       frame_res.GetXaxis().SetTitleOffset(1.)
       frame_res.GetYaxis().SetTitleSize(.13)
       frame_res.GetXaxis().SetTitle("m_{T} (GeV)")
       frame_res.GetXaxis().SetTitleSize(.13)
       frame_res.GetXaxis().SetLabelSize(.12)
       frame_res.GetYaxis().SetLabelSize(.12)
       frame_res.Draw("SAME")
       frame_res.SetTitle("")
       
       c.Update()
       c.Range(1850, -3.5, 8050, 3.5)
       # create two histograms, one for modelA (main) and one for modelB (alt)
       # then divide A by B, and plot in Red
       histA = ROOT.TH1F("histA","histA",65,1850,8050)
       histB = ROOT.TH1F("histB","histB",65,1850,8050)
       histA.Sumw2()
       histB.Sumw2()
       modelA.fillHistogram(histA, RooArgList(varArg))
       modelB.fillHistogram(histB, RooArgList(varArg))
       histA.Divide(histB)
       histA.SetLineColor(ROOT.kMagenta)
       histA.SetLineWidth(2)
       histA.Draw("hist c same")
       frame_res.GetYaxis().SetRangeUser(0.0, 2.0)
       
       c.Update()
       c.Modified()
       c.SaveAs("Residuals_"+ch+"_"+name+"_TwoFuncs_14.pdf")
       c.Delete()

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

def getCard(ch, idir, bias = False, useChi2 = False, verbose = False, doplots = True):
       ch_red = ch[:-5]
       rfilename = idir+"/fitResults_{}.root".format(ch_red)
       print "Opening file ", rfilename
       rfile = ROOT.TFile.Open(rfilename)
       wfilename = idir+"/ws_allFits_{}.root".format(ch_red)
       print "Opening file ", wfilename
       wfile = ROOT.TFile.Open(wfilename)
       ws = wfile.Get("FitWS")

       print "BIAS?", bias
       carddir = ""

       #*******************************************************#
       #                                                       #
       #                   Generate workspace                  #
       #                                                       #
       #*******************************************************#
       mode = "template"
       doModelling = True    
       if(mode == "template"):
              xvarmin = 1850.
              xvarmax = 8050.
              mT = RooRealVar(  "mH"+ch,    "m_{T}", xvarmin, xvarmax, "GeV")

              obsData = ws.data("data_obs")
              nDataEvts = obsData.sum(False)
              nFits = 4
              modelName = "Bkg_"+ch
              nFitsAlt = 4
              modelAltName =  "Bkg_Alt_"+ch
              
              if(doModelling):
                     print "channel: ", ch_red
                     modelBkg = [ws.pdf(modelName+str(i+1)) for i in range(nFits)]
                     fitRes = [rfile.Get("fitresult_{}_data_obs".format(m.GetName())) for m in modelBkg]

                     orderBkg = [len(f.floatParsFinal()) for f in fitRes]

                     if doplots:
                         xframe = mT.frame(ROOT.RooFit.Title("extended ML fit example"))

                         c1 = ROOT.TCanvas()
                         c1.cd()
                         obsData.plotOn(xframe, RooFit.Name("obsData"))

                         modelBkg[0].plotOn(xframe, RooFit.Name("modelBkg1"), ROOT.RooFit.LineColor(ROOT.kPink + 6), RooFit.Range("Full"))
                         modelBkg[1].plotOn(xframe, RooFit.Name("modelBkg2"), ROOT.RooFit.LineColor(ROOT.kBlue -4), RooFit.Range("Full"))
                         modelBkg[2].plotOn(xframe, RooFit.Name("modelBkg3"), ROOT.RooFit.LineColor(ROOT.kRed -4), RooFit.Range("Full"))
                         modelBkg[3].plotOn(xframe, RooFit.Name("modelBkg4"), ROOT.RooFit.LineColor(ROOT.kGreen +1), RooFit.Range("Full"))

                         chi2s_bkg = [xframe.chiSquare("modelBkg{}".format(i+1), "obsData", orderBkg[i]) for i in range(len(modelBkg))]

                         xframe.SetMinimum(0.0002)
                         xframe.Draw()

                         txt1 = ROOT.TText(2000., 10., "model1, nP {}, chi2: {}".format(orderBkg[0],chi2s_bkg[0]))
                         txt1.SetTextSize(0.04) 
                         txt1.SetTextColor(ROOT.kPink+6) 
                         xframe.addObject(txt1) 
                         txt2 = ROOT.TText(2000., 1, "model2, nP {}, chi2: {}".format(orderBkg[1],chi2s_bkg[1]))
                         txt2.SetTextSize(0.04) 
                         txt2.SetTextColor(ROOT.kBlue-4) 
                         xframe.addObject(txt2) 
                         txt3 = ROOT.TText(2000., 0.1, "model3, nP {}, chi2: {}".format(orderBkg[2],chi2s_bkg[2]))
                         txt3.SetTextSize(0.04) 
                         txt3.SetTextColor(ROOT.kRed-4) 
                         xframe.addObject(txt3) 
                         txt4 = ROOT.TText(2000., 0.01, "model4, nP {}, chi2: {}".format(orderBkg[3],chi2s_bkg[3]))
                         txt4.SetTextSize(0.04) 
                         txt4.SetTextColor(ROOT.kGreen+1) 
                         xframe.addObject(txt4) 
                         txt1.Draw()
                         txt2.Draw()
                         txt3.Draw()
                         txt4.Draw()

                         c1.SetLogy()
                         c1.SaveAs("TestAfterFit_"+ch+".pdf")

                     RSS = {orderBkg[i] : getRSS(ch, mT, modelBkg[i], obsData, [fitRes[i]], carddir, doplots, nDataEvts) for i in range(len(modelBkg))}

                     #**********************************************************
                     #                    ALTERNATIVE MODEL                    *
                     #**********************************************************
                     if bias:
                            modelAlt = [ws.pdf(modelAltName+str(i+1)) for i in range(nFitsAlt)]
                            fitResAlt = [rfile.Get("fitresult_{}_data_obs".format(m.GetName())) for m in modelAlt]
                            RSS_alt = {}
                            for i in range(len(modelAlt)):
                                RSS_alt[fitResAlt[i].floatParsFinal().getSize()] = getRSS(ch, mT, modelAlt[i], obsData,  [fitResAlt[i]], carddir, doplots, nDataEvts, label = "alt")

                            if doplots:
                                #for i in range(len(modelBkg)):
                                #    for j in range(len(modelAlt)):
                                #        drawTwoFuncs(ch, mT, modelBkg[i], modelAlt[j], obsData, [fitRes[i], fitResAlt[j]], carddir,  nDataEvts, label = "dual")

                                xframeAlt = mT.frame(ROOT.RooFit.Title("extended ML fit example"))

                                c2 = ROOT.TCanvas()
                                c2.cd()
                                obsData.plotOn(xframeAlt, RooFit.Name("obsData"))

                                modelAlt[0].plotOn(xframeAlt, RooFit.Name("modelAlt1"), ROOT.RooFit.LineColor(ROOT.kPink))
                                modelAlt[1].plotOn(xframeAlt, RooFit.Name("modelAlt2"), ROOT.RooFit.LineColor(ROOT.kBlue))
                                modelAlt[2].plotOn(xframeAlt, RooFit.Name("modelAlt3"), ROOT.RooFit.LineColor(ROOT.kRed))
                                modelAlt[3].plotOn(xframeAlt, RooFit.Name("modelAlt4"), ROOT.RooFit.LineColor(ROOT.kGreen))

                                chi2s_alt = [xframeAlt.chiSquare("modelAlt{}".format(i+1), "obsData", i+1) for i in range(len(modelAlt))]

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
                     report += "func\tchi2\tndof\trchi2\trss\n"
                     nParMin = RSS[fitRes[0].floatParsFinal().getSize()]['npar']
                     for i in xrange(nParMin, nParMin+len(RSS)):
                            report += "   %i\t%.2f\t%i\t%.2f\t%.2f\n" % (i, RSS[i]['chi2'], RSS[i]['nbins']-RSS[i]['npar'], RSS[i]['chiSquared'], RSS[i]['rss'])
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
                     aCrit = 0.07
                     RTDict = {}
                     FDict = {}
                     for o1 in xrange(nParMin, nParMin+len(RSS)-1):
                            for o2 in xrange(o1+1,nParMin+len(RSS)):
                                   print(o1, o2)
                                   if useChi2: qty = 'chi2'
                                   else: qty = 'rss'
                                   RTDict[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)] = fisherTest(RSS[o1][qty], RSS[o2][qty], RSS[o1]['npar'], RSS[o2]['npar'], RSS[o2]["nbins"])
                                   report += "%d par vs %d par: CL=%.5f F_t=%.5f\n" % (RSS[o1]['npar'], RSS[o2]['npar'], RTDict[str(o1)+"v"+str(o2)], FDict[str(o1)+"v"+str(o2)])

                     if RTDict["{}v{}".format(RSS[nParMin]['npar'],RSS[nParMin+1]['npar'])] > aCrit:
                            if RTDict["{}v{}".format(RSS[nParMin]['npar'],RSS[nParMin+2]['npar'])] > aCrit:
                                   if RTDict["{}v{}".format(RSS[nParMin]['npar'],RSS[nParMin+3]['npar'])] > aCrit:
                                          order = RSS[RSS[nParMin]['npar']]['npar']
                                   else:
                                          order = RSS[RSS[nParMin+3]['npar']]['npar']
                            else:
                                   if RTDict["{}v{}".format(RSS[nParMin+2]['npar'],RSS[nParMin+3]['npar'])] > aCrit:
                                          order = RSS[RSS[nParMin+2]['npar']]['npar']
                                   else:
                                          order = RSS[RSS[nParMin+3]['npar']]['npar']
                     else:
                            if RTDict["{}v{}".format(RSS[nParMin+1]['npar'],RSS[nParMin+2]['npar'])] > aCrit:
                                   if RTDict["{}v{}".format(RSS[nParMin+1]['npar'],RSS[nParMin+3]['npar'])] > aCrit:
                                          order = RSS[RSS[nParMin+1]['npar']]['npar']
                                   else:
                                          order = RSS[RSS[nParMin+3]['npar']]['npar']
                            else:
                                   if RTDict["{}v{}".format(RSS[nParMin+2]['npar'],RSS[nParMin+3]['npar'])] > aCrit:
                                          order = RSS[RSS[nParMin+2]['npar']]['npar']
                                   else:
                                          order = RSS[RSS[nParMin+3]['npar']]['npar']
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
                            order = RSS[nParMin+3]["npar"]
                            moreParFlag = 1
                     print "\\hline"
                     print "-"*25   
                     print "Order is", order, "("+ch+")"
                     report+= "aCrit values is {}\n".format(aCrit)
                     report += "Order is %d (%s)\n" % (order, ch)
                     report += ("2 param: " + ", ".join(['%.2f']*len(RSS[nParMin]["parVals"])) + "\n") % tuple(RSS[nParMin]["parVals"])
                     report += ("3 param: " + ", ".join(['%.2f']*len(RSS[nParMin+1]["parVals"])) + "\n") % tuple(RSS[nParMin+1]["parVals"])
                     report += ("4 param: " + ", ".join(['%.2f']*len(RSS[nParMin+2]["parVals"])) + "\n") % tuple(RSS[nParMin+2]["parVals"])
                     report += ("5 param: " + ", ".join(['%.2f']*len(RSS[nParMin+3]["parVals"])) + "\n") % tuple(RSS[nParMin+3]["parVals"])
                     report += "%d (%d) par are sufficient\n" % (RSS[order]['npar'], order)
                     if moreParFlag: report += "BUT really, more Pars are needed!"
                     for i in range(nParMin, nParMin+4):
                            report += "model%i chi2: %.2f\n" % (i, RSS[i]['chiSquared'])

                     ofile.write(report)
                     ofile.close()

                     if bias: 
                            print "Running in BIAS mode"

                            ofile_alt = open("FisherTest_alt_%s.txt"%(ch),"w")
                            report = "Results from Fisher Test on the alternative function for category %s \n" % (ch)
                            report += "func\tchi2\tndof\trchi2\trss\n"
                            nParMin_alt = RSS_alt[fitResAlt[0].floatParsFinal().getSize()]['npar']
                            for i in xrange(nParMin_alt, nParMin_alt+len(RSS_alt)):
                                report += "   %i\t%.2f\t%i\t%.2f\t%.2f\n" % (i, RSS_alt[i]['chi2'], RSS_alt[i]['nbins']-RSS_alt[i]['npar'], RSS_alt[i]['chiSquared'], RSS_alt[i]['rss'])
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
                            RTDict_alt = {}
                            FDict_alt = {}
                            for o1 in xrange(nParMin_alt, nParMin_alt+len(RSS_alt)-1):
                                   for o2 in xrange(o1+1,nParMin_alt+len(RSS_alt)):
                                          if useChi2: qty = 'chi2'
                                          else: qty = 'rss'
                                          RTDict_alt[str(o1)+"v"+str(o2)], FDict_alt[str(o1)+"v"+str(o2)] = fisherTest(RSS_alt[o1][qty], RSS_alt[o2][qty], RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], RSS_alt[o2]["nbins"])
                                          report += "%d par vs %d par: CL=%.5f F_t=%.5f\n" % (RSS_alt[o1]['npar'], RSS_alt[o2]['npar'], RTDict_alt[str(o1)+"v"+str(o2)], FDict_alt[str(o1)+"v"+str(o2)])

                            if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt]['npar'],RSS_alt[nParMin_alt+1]['npar'])] > aCrit:
                                   if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt]['npar'],RSS_alt[nParMin_alt+2]['npar'])] > aCrit:
                                          if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt]['npar'],RSS_alt[nParMin_alt+3]['npar'])] > aCrit:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt]['npar']]['npar']
                                          else:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+3]['npar']]['npar']
                                   else:
                                          if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt+2]['npar'],RSS_alt[nParMin_alt+3]['npar'])] > aCrit:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+2]['npar']]['npar']
                                          else:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+3]['npar']]['npar']
                            else:
                                   if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt+1]['npar'],RSS_alt[nParMin_alt+2]['npar'])] > aCrit:
                                          if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt+1]['npar'],RSS_alt[nParMin_alt+3]['npar'])] > aCrit:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+1]['npar']]['npar']
                                          else:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+3]['npar']]['npar']
                                   else:
                                          if RTDict_alt["{}v{}".format(RSS_alt[nParMin_alt+2]['npar'],RSS_alt[nParMin_alt+3]['npar'])] > aCrit:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+2]['npar']]['npar']
                                          else:
                                                 order_alt = RSS_alt[RSS_alt[nParMin_alt+3]['npar']]['npar']

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
                            report += ("alt 1: " + ", ".join(['%.2f (%.2f)']*len(RSS_alt[nParMin_alt]["parVals"])) + "\n") % tuple(altMerge(RSS_alt[nParMin_alt]["parVals"],RSS_alt[nParMin_alt]["parErr"]))
                            report += ("alt 2: " + ", ".join(['%.2f (%.2f)']*len(RSS_alt[nParMin_alt+1]["parVals"])) + "\n") % tuple(altMerge(RSS_alt[nParMin_alt+1]["parVals"],RSS_alt[nParMin_alt+1]["parErr"]))
                            report += ("alt 3: " + ", ".join(['%.2f (%.2f)']*len(RSS_alt[nParMin_alt+2]["parVals"])) + "\n") % tuple(altMerge(RSS_alt[nParMin_alt+2]["parVals"],RSS_alt[nParMin_alt+2]["parErr"]))
                            report += ("alt 4: " + ", ".join(['%.2f (%.2f)']*len(RSS_alt[nParMin_alt+3]["parVals"])) + "\n") % tuple(altMerge(RSS_alt[nParMin_alt+3]["parVals"],RSS_alt[nParMin_alt+3]["parErr"]))
                            report += "%d par are sufficient\n" % (RSS_alt[order_alt]['npar'])
                            if moreParFlag_alt: report += "BUT really, more Pars are needed!"
                            iPar, iPar_alt = nParMin, nParMin_alt
                            for i in range(1,5):
                                   report += "model%i chi2, chi2_alt: %.2f %.2f\n" % (i, RSS[iPar]['chiSquared'],RSS_alt[iPar_alt]['chiSquared'])
                                   iPar+=1
                                   iPar_alt+=1

                            ofile_alt.write(report)
                            ofile_alt.close()
                     modelBkgF = None
                     for i in range(len(modelBkg)):
                         if order==RSS[nParMin+i]['npar']:
                            modelBkgF = modelBkg[i]
                            print "Model Bkg Saved is", order
                            break
                     if modelBkgF is None:
                            print "Main functions with", RSS[nParMin+3]['npar']+1, "or more parameters are needed to fit the background"
                            #exit()
                            modelBkgF = modelBkg[-1]

                     modelBkgF.SetName(modelName)
                     modelBkgF, objsBkg = remakePdf(modelBkgF)


                     #*******************************************************#
                     #                                                       #
                     #                  Saving RooWorkspace                  #
                     #                                                       #
                     #*******************************************************#

                     wsfilename = "ws_{}.root".format(ch_red)
                     w_ = RooWorkspace("BackgroundWS", "workspace")
                     print "Storing ", modelName
                     getattr(w_, "import")(modelBkgF) #Bkg func with optimal num Paras
                     getattr(w_, "import")(obsData) # data_obs histogram

                     if bias:
                        modelAltF = None
                        for i in range(len(modelAlt)):
                            if order_alt==RSS_alt[nParMin_alt+i]['npar']:
                                modelAltF = modelAlt[i]
                                print "Model Alt Saved is", order_alt
                                break
                        if modelAltF is None:
                               print "Alt functions with", RSS_alt[nParMin_alt+3]['npar']+1, "or more parameters are needed to fit the background"
                               #exit()
                               modelAltF = modelAlt[-1]
                        modelAltF.SetName(modelAltName)
                        modelAltF, objsAlt = remakePdf(modelAltF)
                        getattr(w_, "import")(modelAltF) # Alt func with optimal num Paras

                     wstatus = w_.writeToFile(wsfilename, True)
