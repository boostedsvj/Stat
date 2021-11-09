import ROOT
import json

#changed to no longer need the settings.py file
# -syst list is created here
# -rateParams list is created here
# -processes list is created here
# files now saved to base directory and exported to EOS after jobs completion


from ROOT import RooRealVar, RooDataHist, RooArgList, RooGenericPdf, RooBernstein, RooExtendPdf, RooCmdArg, RooWorkspace, RooFit, RooDataSet, RooArgSet, RooCategory, RooFitResult, RooCurve, RooParametricShapeBinPdf
from Stat.Limits.bruteForce import remakePdf
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
# uncertainties apply to signal only
# data doesn't get any (cause its data)
# and MC bkg is indirectly used (via the fits)
# and the uncertainties on the fit parameters take care of that
syst["lumi"] = ("lnN", "sig", 1.026)
syst["trig"] = ("lnN", "sig", 1.020)
syst["scale"] = ("lnN", "sig", 1.15)

'''syst["mcstat"] = ("shape", ["sig"])

syst["MC2016JEC"] = ("shape",["sig"])
syst["MC2016JER"] = ("shape",["sig"])
syst["MC2016puunc"] = ("shape",["sig"])
syst["MC2016trigfnunc"] = ("shape",["sig"])

syst["MC2017JEC"] = ("shape",["sig"])
syst["MC2017JER"] = ("shape",["sig"])
syst["MC2017puunc"] = ("shape",["sig"])
syst["MC2017trigfnunc"] = ("shape",["sig"])

syst["MC2018JEC"] = ("shape",["sig"])
syst["MC2018JER"] = ("shape",["sig"])
syst["MC2018puunc"] = ("shape",["sig"])
syst["MC2018trigfnunc"] = ("shape",["sig"])

syst["MCRun2JES"] = ("shape",["sig"])
syst["MCRun2pdfallunc"] = ("shape",["sig"])
syst["MCRun2psfsrunc"] = ("shape",["sig"])
syst["MCRun2psisrunc"] = ("shape",["sig"])'''


processes = ["QCD"]

#*******************************************************#
#                                                       #
#                     Utility Functions                 #
#                                                       #
#*******************************************************#

isData = True

def getRate(ch, process, ifile):
       hName = ch + "/"+ process
       h = ifile.Get(hName)
       #return h.Integral(1,h.GetXaxis().GetNbins()-1)
       return h.Integral(1,h.GetXaxis().GetNbins())

def getHist(ch, process, ifile):
       hName = ch + "/"+ process
       print"ch is: %s, process is: %s"%(ch, process)
       print "Histo Name ", hName
       h = ifile.Get(hName)
       h.SetDirectory(0)
       return h


#*******************************************************#
#                                                       #
#                      Datacard                         #
#                                                       #
#*******************************************************#

#def getCard(sig, ch, ifilename, wsdir, doModelling, mode = "histo", bias = False, verbose = False, doSys = True):
def getCard(sig, ch, ifilename, wsdir, doModelling, mode = "template", bias = True, verbose = False, doSys = False):
       try:
              ifile = ROOT.TFile.Open(ifilename)
       except IOError:
              print "Cannot open ", ifilename
       else:
              print "Opening file ",  ifilename
              ifile.cd()


       print "BIAS?", bias
       carddir = ""

       ch = 'bsvj'
       hist_filename = os.getcwd()+"/"+ifilename
       hist = getHist(ch, sig, ifile)

       


       #*******************************************************#
       #                                                       #
       #                   Generate workspace                  #
       #                                                       #
       #*******************************************************#
    
       if(mode == "template"):

              #ch_red = ch[:-5]
	      ch_red = ch
              modelName = "Bkg_"+ch
              modelAltName =  "Bkg_Alt_"+ch
              histData = getHist(ch, "data_obs", ifile)
              histSig = getHist(ch, sig, ifile)
              xvarmin = 190.
	      xvarmax = 502.
              binMin = histData.FindBin(xvarmin)
              binMax = histData.FindBin(xvarmax)-1
              nDataEvts = histData.Integral(binMin, binMax)
              print("binMin: %d, binMax: %d")%(binMin, binMax)
              print("number of data events from the integral of the histogram: ", nDataEvts)
              nSigEvts = histSig.Integral(binMin, binMax)
              mT = RooRealVar(  "mH"+ch,    "m_{T}", xvarmin, xvarmax, "GeV")
              obsData = RooDataHist("data_obs", "Data",  RooArgList(mT), histData, 1.)
              sigData = RooDataHist("sigdata", "MC Sig",  RooArgList(mT), histSig, 1.)
              #*******************************************************#
              #                                                       #
              #                    Datacard Writing                   #
              #                                                       #
              #*******************************************************#
              wfile = ROOT.TFile.Open(wsdir+"/ws_{}.root".format(ch_red))
              wBkg =  wfile.Get("BackgroundWS")
              print "workspace ", wBkg.Print()
              print "Model name: ", modelName
              modelBkg = wBkg.pdf(modelName)
              modelBkg, objsBkg = remakePdf(modelBkg)
              parSet = modelBkg.getParameters(obsData)
              parSet.Print()
              argList = ROOT.RooArgList(parSet)
              parNames = [ argList[i].GetName() for i in xrange(0, len(parSet))]                                 


              for i in xrange(0, len(parSet)):
                    print  argList[i].GetName(), "   ", argList[i].getVal()

              if bias:
                     modelAlt = wBkg.pdf(modelAltName)
                     modelAlt, objsAlt = remakePdf(modelAlt)
                     ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
                     from ROOT import RooMultiPdf
                     pdf_index_string = "pdf_index_%s" % (ch)
                     cat = RooCategory(pdf_index_string, "Index of Pdf which is active")
                     pdfs = RooArgList()
                     pdfs.add(modelBkg)
                     pdfs.add(modelAlt)
                     roomultipdf = RooMultiPdf("roomultipdf", "All Pdfs", cat, pdfs)
                     #normulti = RooRealVar("roomultipdf_norm", "Number of background events", nDataEvts, 0., 1.e6)
                     normulti = RooRealVar("roomultipdf_norm", "Number of background events", 1.0, 0., 1.e6)

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

              '''for i in xrange(hist.GetNbinsX()):
                     for year in ["2016","2017","2018"]:
                            mcstatSysName = "mcstat_%s_%s_MC%sbin%d"  % (ch, sig, year, i+1)
                            mcstatSigUp = getHist(ch, sig + "_" + mcstatSysName + "Up", ifile)

                            mcstatSigDown = getHist(ch, sig + "_" + mcstatSysName + "Down", ifile)
                            mcstatSigHistUp = RooDataHist(sig + "_" + mcstatSysName + "Up", "Data (MC sig)",  RooArgList(mT), mcstatSigUp, 1.)
                            mcstatSigHistDown = RooDataHist(sig + "_" + mcstatSysName + "Down", "Data (MC sig)",  RooArgList(mT), mcstatSigDown, 1.)
                            getattr(w, "import")(mcstatSigHistUp, RooFit.Rename(sig + "_" + mcstatSysName + "Up") )
                            getattr(w, "import")(mcstatSigHistDown, RooFit.Rename(sig + "_" + mcstatSysName + "Down") )

              for sysName,sysValue  in syst.iteritems():
                     if not ((doSys == True) or (sysName == "lumi")):
                            continue # always include lumi, and only do others if doSys is True
                     if(sysValue[0]=="shape" and "mcstat" not in sysName):              
                            sysUp =  getHist(ch, sig + "_" + sysName + "Up", ifile)
                            sysDown =  getHist(ch, sig + "_" + sysName + "Down", ifile)
                            sysSigHistUp = RooDataHist(sig + "_" + sysName + "Up", sysName + " uncertainty",  RooArgList(mT), sysUp, 1.)
                            sysSigHistDown = RooDataHist(sig + "_" + sysName + "Down", sysName + " uncertainty",  RooArgList(mT), sysDown, 1.)
                            getattr(w, "import")(sysSigHistUp, RooFit.Rename(sig + "_" + sysName + "Up") )
                            getattr(w, "import")(sysSigHistDown, RooFit.Rename(sig + "_" + sysName + "Down") )'''
                            


              getattr(w, "import")(modelBkg, (modelBkg.GetName()))
              if bias:
                     getattr(w, "import")(modelAlt, (modelAltName))
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

              processes = ["Bkg"]
              rates["Bkg"] = nDataEvts #nBkgEvts, Kevin asked to change to data yield, limits converge better
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
       rates[sig] = getRate(ch, sig, ifile)
       if mode == "template":
              print "TEST sig: ", rates[sig], nSigEvts, binMax
              print "TEST bkg/obs: ", rates["data_obs"], nDataEvts, binMax




       card  = "imax 1 number of channels \n"
       card += "jmax * number of backgrounds \n"
       card += "kmax * number of nuisance parameters\n"
       card += "-----------------------------------------------------------------------------------\n"

       if(mode == "template"):
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
              print(sysName, sysValue, "testSystematicValues")
              if not ((doSys == True) or (sysName == "lumi")):
                  continue # always include lumi, and only do others if doSys is True
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
                            if("sig" in sysValue[1]):
                                   try:
                                          if(getRate(ch, sig, ifile) != 0.): sigSys = abs((getRate(ch, sig+hsysNameUp, ifile) - getRate(ch, sig+hsysNameDown, ifile))/ (2* getRate(ch, sig, ifile)))
                                   except AttributeError:
                                          sigSys = sysValue[2]
                                   if(sigSys < 1. and sigSys > 0.): sigSys = sigSys + 1
                                   card += "%-20s" % (sigSys)
                            else:  card += "%-20s" % ("-")
                            for p in processes:
                                   print(p, processes)
                                   if (p in sysValue[1]):
                                          if (getRate(ch, p, ifile) != 0.): bkgSys = abs((getRate(ch, p+hsysNameUp, ifile) - getRate(ch, p+hsysNameDown, ifile))/ (2* getRate(ch, p, ifile)) )
                                          else: bkgSys = 1
                                          if(bkgSys<1. and bkgSys >0.): bkgSys = bkgSys + 1
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
                                          for year in ["2016","2017","2018"]:
                                                 sysName = "mcstat_%s_%s_MC%sbin%d      "  % (ch, sampName, year, i+1)
                                                 card += "%-20s   shape   " % (sysName)
                                                 card += line   
                                                 card += "\n"    
       
              card += "\n"
              
       if mode == "template":
              for par in parNames: card += "%-20s%-20s\n" % (par, "flatParam")

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

       return card
