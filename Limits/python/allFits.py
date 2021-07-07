import ROOT
import json

#changed to no longer need the settings.py file
# files now saved to base directory and exported to EOS after jobs completion

#allFits.py - part 0 of separating ftests from datacard writing
# this part will create a ws.root file containing the histograms, pdfs (all orders), and fit results


from ROOT import RooRealVar, RooDataHist, RooArgList, RooGenericPdf, RooBernstein, RooExtendPdf, RooCmdArg, RooWorkspace, RooFit, RooDataSet, RooArgSet, RooCategory, RooFitResult, RooCurve, RooParametricShapeBinPdf
import os, sys
from array import array
import copy, math, pickle
import collections
import numpy as np
from numpy import ndarray
from bruteForce import VarInfo, PdfInfo, varToInfo, bruteForce, silence
silence()

mopt = ROOT.Math.MinimizerOptions()
mopt.SetMaxFunctionCalls(100000)
mopt.SetMaxIterations(100000)

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


#*******************************************************#
#                                                       #
#                      Datacard                         #
#                                                       #
#*******************************************************#

def getCard(sig, ch, ifilename, npool = 1, initvals = [1.0], bias = False, verbose = False):
       try:
              ifile = ROOT.TFile.Open(ifilename)
       except IOError:
              print "Cannot open ", ifilename
       else:
              print "Opening file ",  ifilename
              ifile.cd()


       print "BIAS?", bias

       hist_filename = os.getcwd()+"/"+ifilename
       hist = getHist(ch, sig, ifile)

       #*******************************************************#
       #                                                       #
       #                   Generate workspace                  #
       #                                                       #
       #*******************************************************#
       mode = "template"
       doModelling = True
       if(mode == "template"):

              ###ATT: CHECK BKG AND DATA NORMALIZATION AND DISTRIBUTION
              histBkgData = getHist(ch, "Bkg", ifile)
              histData = getHist(ch, "data_obs", ifile)
              print "channel ", ch             
              print "signal ", sig

              xvarmin = 1850.
              xvarmax = 8050.
              mT = RooRealVar(  "mH"+ch,    "m_{T}", xvarmin, xvarmax, "GeV")
              binMin = histData.FindBin(xvarmin)
              binMax = histData.FindBin(xvarmax)
              obsData = RooDataHist("data_obs", "Data",  RooArgList(mT), histData, 1.)
              print "Bkg Integral: ", histData.Integral() 
              nBkgEvts = histBkgData.Integral(binMin, binMax)
              nDataEvts = histData.Integral(binMin, binMax)

              print "channel: ", ch
              normBkg = RooRealVar("Bkg_"+ch+"_norm", "Number of background events", nBkgEvts, 0., 2.e4)
              normData = RooRealVar("Data_"+ch+"_norm", "Number of background events", nDataEvts, 0., 2.e4)
              ch_red = ch[:-5]
              modelName = "Bkg_"+ch
              modelAltName =  "Bkg_Alt_"+ch
              
              if(doModelling):
                     print "channel: ", ch_red
                     lowerLimit = -50
                     upperLimit = 150
                     if "lowCut" in ch_red:
                            lowerLimit = -60
                            upperLimit = 150
                     if "highCut" in ch_red:
                            lowerLimit = -50
                            upperLimit = 150
                     if "lowSVJ2" in ch_red:
                            lowerLimit = -80
                            upperLimit = 190
                     if "highSVJ2" in ch_red:
                            lowerLimit = -20
                            upperLimit = 100


                     p1_1 = RooRealVar(ch_red + "_p1_1", "p1", 1., lowerLimit, upperLimit)
                     p1_2 = RooRealVar(ch_red + "_p1_2", "p1", 1., lowerLimit, upperLimit)
                     p1_3 = RooRealVar(ch_red + "_p1_3", "p1", 1., lowerLimit, upperLimit)
                     p1_4 = RooRealVar(ch_red + "_p1_4", "p1", 1., lowerLimit, upperLimit)

                     p2_1 = RooRealVar(ch_red + "_p2_1", "p2", 1., lowerLimit, upperLimit)
                     p2_2 = RooRealVar(ch_red + "_p2_2", "p2", 1., lowerLimit, upperLimit)
                     p2_3 = RooRealVar(ch_red + "_p2_3", "p2", 1., lowerLimit, upperLimit)
                     p2_4 = RooRealVar(ch_red + "_p2_4", "p2", 1., lowerLimit, upperLimit)

                     p3_2 = RooRealVar(ch_red + "_p3_2", "p3", 1., lowerLimit, upperLimit)
                     p3_3 = RooRealVar(ch_red + "_p3_3", "p3", 1., lowerLimit, upperLimit)
                     p3_4 = RooRealVar(ch_red + "_p3_4", "p3", 1., lowerLimit, upperLimit)

                     p4_3 = RooRealVar(ch_red + "_p4_3", "p4", 1., lowerLimit, upperLimit)
                     p4_4 = RooRealVar(ch_red + "_p4_4", "p4", 1., lowerLimit, upperLimit)

                     p5_4 = RooRealVar(ch_red + "_p5_4", "p5", 1., lowerLimit, upperLimit)


                     #Function from Theorists, combo testing, sequence E, 1, 11, 12, 22
                     # model NM has N params on 1-x and M params on x. exponents are (p_i + p_{i+1} * log(x))
                     # these are the RooGenericPdf verisons, convert to RooParametricShapeBinPdf below
                     modelBkg1_rgp = RooGenericPdf(modelName+"1_rgp", "Thry. fit (11)", "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2))",                                                                  RooArgList(mT, p1_1, p2_1))
                     modelBkg2_rgp = RooGenericPdf(modelName+"2_rgp", "Thry. fit (12)", "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)))",                                                 RooArgList(mT, p1_2, p2_2, p3_2))
                     #modelBkg3_rgp = RooGenericPdf(modelName+"3_rgp", "Thry. fit (13)", "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2)))",                         RooArgList(mT, p1_3, p2_3, p3_3, p4_3))
                     #modelBkg4_rgp = RooGenericPdf(modelName+"4_rgp", "Thry. fit (14)", "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2)+@5*pow(log(@0/13000),3)))", RooArgList(mT, p1_4, p2_4, p3_4, p4_4, p5_4))
                     modelBkg3_rgp = RooGenericPdf(modelName+"3_rgp", "Thry. fit (22)", "pow(1 - @0/13000, @1+@2*log(@0/13000)) * pow(@0/13000, -(@3+@4*log(@0/13000)))",                         RooArgList(mT, p1_3, p2_3, p3_3, p4_3))
                     modelBkg4_rgp = RooGenericPdf(modelName+"4_rgp", "Thry. fit (32)", "pow(1 - @0/13000, @1+@2*log(@0/13000)+@3*pow(log(@0/13000),2)) * pow(@0/13000, -(@4+@5*log(@0/13000)))", RooArgList(mT, p1_4, p2_4, p3_4, p4_4, p5_4))
                     #modelBkg4_rgp = RooGenericPdf(modelName+"4_rgp", "Thry. fit (41)", "pow(1 - @0/13000, @1+@2*log(@0/13000)+@3*pow(log(@0/13000),2)+@4*pow(log(@0/13000),3)) * pow(@0/13000, -@5)", RooArgList(mT, p1_4, p2_4, p3_4, p4_4, p5_4))
                     modelBkg = [
                        RooParametricShapeBinPdf(modelName+"1", "Thry. Fit (11)", modelBkg1_rgp, mT, RooArgList(p1_1, p2_1), histBkgData),
                        RooParametricShapeBinPdf(modelName+"2", "Thry. Fit (12)", modelBkg2_rgp, mT, RooArgList(p1_2, p2_2, p3_2), histBkgData),
                        RooParametricShapeBinPdf(modelName+"3", "Thry. Fit (22)", modelBkg3_rgp, mT, RooArgList(p1_3, p2_3, p3_3, p4_3), histBkgData),
                        RooParametricShapeBinPdf(modelName+"4", "Thry. Fit (32)", modelBkg4_rgp, mT, RooArgList(p1_4, p2_4, p3_4, p4_4, p5_4), histBkgData),
                     ]

                     fitrange = "Full"
                     fitRes = [modelBkg[i].fitTo(obsData, RooFit.Extended(False), RooFit.Save(1), RooFit.SumW2Error(True), RooFit.Strategy(2), RooFit.Minimizer("Minuit2"), RooFit.PrintLevel(-1 if not verbose else 2), RooFit.Range(fitrange)) for i in range(len(modelBkg))]

                     #**********************************************************
                     #                    ALTERNATIVE MODEL                    *
                     #**********************************************************
                     if bias:

                            normAlt = RooRealVar("Bkg_"+ch+"alt_norm", "Number of background events", nBkgEvts, 0., 2.e4)
                            normData = RooRealVar("Data_"+ch+"alt_norm", "Number of background events", nDataEvts, 0., 2.e4) 

                            lowAlt = -100
                            highAlt = 100
                            pdfsAlt = [
                                PdfInfo(modelAltName+"1", "Alt. Fit 1par", "exp(@1*(@0/13000))", hist=histBkgData,
                                    x = varToInfo(mT, True),
                                    pars = [
                                        VarInfo(ch_red + "_p1_1_alt", "p1", 1., lowAlt, highAlt, 0, "", False),
                                    ],
                                ),
                                PdfInfo(modelAltName+"2", "Alt. Fit 2par", "exp(@1*(@0/13000)) * pow(@0/13000,@2)", hist=histBkgData,
                                    x = varToInfo(mT, True),
                                    pars = [
                                        VarInfo(ch_red + "_p1_2_alt", "p1", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p2_2_alt", "p2", 1., lowAlt, highAlt, 0, "", False),
                                    ],
                                ),
                                PdfInfo(modelAltName+"3", "Alt. Fit 3par", "exp(@1*(@0/13000)) * pow(@0/13000,@2*(1+@3*log(@0/13000)))", hist=histBkgData,
                                    x = varToInfo(mT, True),
                                    pars = [
                                        VarInfo(ch_red + "_p1_3_alt", "p1", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p2_3_alt", "p2", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p3_3_alt", "p3", 1., lowAlt, highAlt, 0, "", False),
                                    ],
                                ),
                                PdfInfo(modelAltName+"4", "Alt. Fit 4par", "exp(@1*(@0/13000)) * pow(@0/13000,@2*(1+@3*log(@0/13000)*(1+@4*log(@0/13000))))", hist=histBkgData,
                                    x = varToInfo(mT, True),
                                    pars = [
                                        VarInfo(ch_red + "_p1_4_alt", "p1", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p2_4_alt", "p2", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p3_4_alt", "p3", 1., lowAlt, highAlt, 0, "", False),
                                        VarInfo(ch_red + "_p4_4_alt", "p4", 1., lowAlt, highAlt, 0, "", False),
                                    ],
                                ),
                            ]
                            modelAlt = []
                            objsAlt = []
                            fitResAlt = []
                            for pdfAlt in pdfsAlt:
                                print "fit ",pdfAlt.name
                                # truncate brute force at 3 param (4 takes too long, has some segfaults)
                                #mtmp, otmp, ftmp = bruteForce(pdfAlt, obsData, initvals, npool, [3,1.], 2)
                                mtmp, otmp, ftmp = bruteForce(pdfAlt, obsData, initvals, npool, None, 2)
                                modelAlt.append(mtmp)
                                objsAlt.append(otmp)
                                fitResAlt.append(ftmp)


                     #*******************************************************#
                     #                                                       #
                     #                  Saving RooWorkspace                  #
                     #                                                       #
                     #*******************************************************#

                     rfilename = "fitResults_{}.root".format(ch_red)
                     rfile_ = ROOT.TFile.Open(rfilename, "RECREATE")
                     w_ = RooWorkspace("FitWS", "workspace")
                     getattr(w_, "import")(obsData)
                     for i in range(len(modelBkg)):
                        getattr(w_, "import")(modelBkg[i])
                        fitRes[i].Write()
                     if bias:
                        for i in range(len(modelAlt)):
                            getattr(w_, "import")(modelAlt[i])
                            fitResAlt[i].Write()
                     rfile_.Close()

                     wsfilename = "ws_allFits_{}.root".format(ch_red)
                     wstatus = w_.writeToFile(wsfilename, True)
