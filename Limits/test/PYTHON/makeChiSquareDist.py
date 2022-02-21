import ROOT as rt
import os, sys
import subprocess
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)
rt.gStyle.SetPalette(rt.kRainBow)
if len(sys.argv) < 2:
	print("must include EOSDIR option. Try again. (/store/user/cfallon/<INPUT>/)")
	exit()
eosDir = sys.argv[1]
eosArea = "root://cmseos.fnal.gov//store/user/cfallon/"+eosDir+"/"


listOfParams1 = [
['1500', '20', '03', 'peak'],
['1700', '20', '03', 'peak'],
['1900', '20', '03', 'peak'],
['2100', '20', '03', 'peak'],
['2300', '20', '03', 'peak'],
['2500', '20', '03', 'peak'],
['2700', '20', '03', 'peak'],
['2900', '20', '03', 'peak'],
['3100', '20', '03', 'peak'],
['3300', '20', '03', 'peak'],
['3500', '20', '03', 'peak'],
['3700', '20', '03', 'peak'],
['3900', '20', '03', 'peak'],
['4100', '20', '03', 'peak'],
['4300', '20', '03', 'peak'],
['4500', '20', '03', 'peak'],
['4700', '20', '03', 'peak'],
['4900', '20', '03', 'peak'],
['5100', '20', '03', 'peak']]

baseline = [['3100', '20', '03', 'peak']]

#regions = ["highCut"]
#regions = ["lowSVJ2","highCut","highSVJ2"]
regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
#expSig = ""# "" for excpSig = 0, "_expSig1" for expSig = 1, add "_extra" for SVJ options, or nothing for Dijet options
#n = 0 if expSig == "" else int(expSig.split("_")[1][-1:])


pae = ""
pae2 = ""

i = 0
# make histograms for each Signal Injection (0,1) and each genFunc (Main, Alt) combo
# that are distributed in mZ
"""
makeNew = True
if makeNew:
	outputROOTFile = rt.TFile.TFile("massScan.root","recreate")
	hist_0M_lc = rt.TH1F("hist_0M_lc","Gaussian Mean for lowCut, r_{inj} = 0, genMain", 19, 1500,5200)
	hist_1M_lc = rt.TH1F("hist_1M_lc","Gaussian Mean for lowCut, r_{inj} = 1, genMain", 19, 1500,5200)
	hist_0A_lc = rt.TH1F("hist_0A_lc","Gaussian Mean for lowCut, r_{inj} = 0, genAlt", 19, 1500,5200)
	hist_1A_lc = rt.TH1F("hist_1A_lc","Gaussian Mean for lowCut, r_{inj} = 1, genAlt", 19, 1500,5200)
	hist_0M_l2 = rt.TH1F("hist_0M_l2","Gaussian Mean for lowSVJ2, r_{inj} = 0, genMain", 19, 1500,5200)
	hist_1M_l2 = rt.TH1F("hist_1M_l2","Gaussian Mean for lowSVJ2, r_{inj} = 1, genMain", 19, 1500,5200)
	hist_0A_l2 = rt.TH1F("hist_0A_l2","Gaussian Mean for lowSVJ2, r_{inj} = 0, genAlt", 19, 1500,5200)
	hist_1A_l2 = rt.TH1F("hist_1A_l2","Gaussian Mean for lowSVJ2, r_{inj} = 1, genAlt", 19, 1500,5200)
	hist_0M_hc = rt.TH1F("hist_0M_hc","Gaussian Mean for highCut, r_{inj} = 0, genMain", 19, 1500,5200)
	hist_1M_hc = rt.TH1F("hist_1M_hc","Gaussian Mean for highCut, r_{inj} = 1, genMain", 19, 1500,5200)
	hist_0A_hc = rt.TH1F("hist_0A_hc","Gaussian Mean for highCut, r_{inj} = 0, genAlt", 19, 1500,5200)
	hist_1A_hc = rt.TH1F("hist_1A_hc","Gaussian Mean for highCut, r_{inj} = 1, genAlt", 19, 1500,5200)
	hist_0M_h2 = rt.TH1F("hist_0M_h2","Gaussian Mean for highSVJ2, r_{inj} = 0, genMain", 19, 1500,5200)
	hist_1M_h2 = rt.TH1F("hist_1M_h2","Gaussian Mean for highSVJ2, r_{inj} = 1, genMain", 19, 1500,5200)
	hist_0A_h2 = rt.TH1F("hist_0A_h2","Gaussian Mean for highSVJ2, r_{inj} = 0, genAlt", 19, 1500,5200)
	hist_1A_h2 = rt.TH1F("hist_1A_h2","Gaussian Mean for highSVJ2, r_{inj} = 1, genAlt", 19, 1500,5200)
	listOfParams = listOfParams1
else:
	outputROOTFile = rt.TFile.TFile("massScan.root","update")
	listOfParams = listOfParams2
"""
choices = subprocess.check_output(["eos","root://cmseos.fnal.gov","ls",eosArea.split("//")[2]])

for sigPars in listOfParams1:
	#continue
	SVJNAME = "SVJ_mZprime{}_mDark{}_rinv{}_alpha{}".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3])
	if not (SVJNAME in choices): continue
	print("************************",SVJNAME,"************************")
	print("files to be opened from: " + eosArea + SVJNAME)
	for region in regions:
		print("/ws_"+SVJNAME+"_{}_2018_template.root".format(region))
		wsFile = rt.TFile.Open(eosArea+SVJNAME+"/ws_"+SVJNAME+"_{}_2018_template.root".format(region),"read")
		# signal is a RooDataHist within SVJ RooWorkspace
		svjWS = wsFile.Get("SVJ")
		sigDataSet = rt.RooDataHist()
		sigDataSet = svjWS.data(SVJNAME)
		sigHist = sigDataSet.createHistogram("mH"+region+"_2018")
		#for expSig in ["Sig0GenAltFitMain","Sig1GenAltFitMain","Sig0GenMainFitMain","Sig1GenMainFitMain"]:#FitMain
		for expSig in ["Sig0GenAltFitAlt","Sig1GenAltFitAlt","Sig0GenMainFitAlt","Sig1GenMainFitAlt"]:#FitAlt
			for combineOpts in [""]:#["OptS","OptD"]:
				bigName = region+combineOpts+expSig
				rGenCode = expSig[3]+expSig[7]
				n = int(expSig[3])
				print("/fitDiagnostics"+bigName+".mH120.123456.root")
				print("/higgsCombine"+bigName+".FitDiagnostics.mH120.123456.root")
				print("/higgsCombine"+bigName.split("Fit")[0]+".MultiDimFit.mH120.root")
				fitDiagFile = rt.TFile.Open(eosArea + SVJNAME + "/fitDiagnostics"+bigName+".mH120.123456.root","read")
				fitOnlyFile = rt.TFile.Open(eosArea + SVJNAME + "/higgsCombine"+bigName+".FitDiagnostics.mH120.123456.root","read")
				mdfFile = rt.TFile.Open(eosArea + SVJNAME + "/higgsCombine"+bigName.split("Fit")[0]+".MultiDimFit.mH120.root","read")
				#dataObsFile = rt.TFile.Open(eosArea + SVJNAME + "/ws_"+SVJNAME+"_"+region+"_2018_template.root","read")
				if ( type(fitDiagFile) == type(rt.TFile())) or ( type(fitOnlyFile) == type(rt.TFile())) or (type(mdfFile) == type(rt.TFile())):
				#if (type(fitDiagFile) == type(rt.TFile())) or (type(fitOnlyFile) == type(rt.TFile())):
					continue
				wsMDF = mdfFile.Get("w")
				snapshotArgList = wsMDF.getSnapshot("MultiDimFit")
				#svjWs = dataObsFile.Get("SVJ")
				data = rt.RooDataHist()
				data = svjWS.data("data_obs")
				#genPdf = svjWs.pdf("Bkg_"+region+"_2018")
				tree = fitDiagFile.Get("tree_fit_sb")
				limit = fitOnlyFile.Get("limit")
				limit.Print()
				mT = rt.RooRealVar("mH"+region+"_2018","m_{T}", 1500., 8000., "GeV")
				dataHist = data.createHistogram("histBkg",mT,rt.RooFit.Binning(65, 1500, 8000))
				numNonZeroBinsInBkg = 0
				for iBin in range(dataHist.GetNbinsX()):
					if dataHist.GetBinContent(iBin) >= 0.01: numNonZeroBinsInBkg +=1
				print("Number of NonZero Bins in Bkg Histogram: {}".format(numNonZeroBinsInBkg))
				
				chi2Div = 1.55


				# updated nPar and nParAlt 2Feb21
				#
				#      | lC | l2 | hC | h2 
				# main |  2 |  2 |  3 |  2
				# alt  |  2 |  2 |  3 |  1
				if region == "highCut":
					regCode = "hc"
					nPar = 3
					nParAlt = 3
					eitLimitUp = 12000 # events in toy
					eitLimitDown = 10000
					rDownLim = -4
					rUpLim = 4
					downLimitrErr = 0
					upLimitrErr = 2
					downLimitRMU = -4
					upLimitRMU = 4
					parLimDown1 = -100
					parLimUp1 = 100
					parLimDown2 = -100
					parLimUp2 = 100
					parLimDown3 = -100
					parLimUp3 = 100
					parLimDown4 = -100
					parLimUp4 = 100
					parLimDown5 = -100
					parLimUp5 = 100
					chi2LimitDown = 0
					chi2LimitUp = 200
				elif region == "highSVJ2":
					regCode = "h2"
					nPar = 2
					nParAlt = 2
					eitLimitUp = 1000 # events in toy
					eitLimitDown = 0
					rDownLim = -2
					rUpLim = 2
					downLimitrErr = 0
					upLimitrErr = 1
					downLimitRMU = -4
					upLimitRMU = 4
					parLimDown1 = -1000
					parLimUp1 = 1000
					parLimDown2 = -1000
					parLimUp2 = 1000
					parLimDown3 = -1000
					parLimUp3 = 1000
					parLimDown4 = -1000
					parLimUp4 = 1000
					parLimDown5 = -1000
					parLimUp5 = 1000
					chi2LimitDown = 0
					chi2LimitUp = 200
				elif region == "lowCut":
					regCode = "lc"
					nPar = 2
					nParAlt = 2
					eitLimitUp = 80000 # events in toy
					eitLimitDown = 70000
					rDownLim = -5
					rUpLim = 5
					downLimitrErr = 0
					upLimitrErr = 2
					downLimitRMU = -5
					upLimitRMU = 5
					parLimDown1 = -100
					parLimUp1 = 100
					parLimDown2 = -100
					parLimUp2 = 100
					parLimDown3 = -100
					parLimUp3 = 100
					parLimDown4 = -100
					parLimUp4 = 100
					parLimDown5 = -100
					parLimUp5 = 100
					chi2LimitDown = 0
					chi2LimitUp = 200
				elif region == "lowSVJ2":
					regCode = "l2"
					nPar = 2
					nParAlt = 2
					eitLimitUp = 2000 # events in toy
					eitLimitDown = 0
					rDownLim = -2
					rUpLim = 2
					downLimitrErr = 0
					upLimitrErr = 1
					downLimitRMU = -6
					upLimitRMU = 6
					parLimDown1 = -100
					parLimUp1 = 100
					parLimDown2 = -100
					parLimUp2 = 100
					parLimDown3 = -100
					parLimUp3 = 100
					parLimDown4 = -100
					parLimUp4 = 100
					parLimDown5 = -100
					parLimUp5 = 100
					chi2LimitDown = 0
					chi2LimitUp = 200

				
				nBins = 50
				rHist = rt.TH1F("rHist","r, all Events",nBins,rDownLim,rUpLim)
				rHistChi2Up = rt.TH1F("rHistChi2Up","r, chi2 > {}".format(chi2Div),nBins,rDownLim,rUpLim)
				rHistChi2Down = rt.TH1F("rHistChi2Down","r, chi2 < {}".format(chi2Div),nBins,rDownLim,rUpLim)

				rErrHist = rt.TH1F("rErrHist","rErr, all Events",nBins,downLimitrErr,upLimitrErr)
				rErrHistChi2Up = rt.TH1F("rErrHistChi2Up","rErr, chi2 > {}".format(chi2Div),nBins,downLimitrErr,upLimitrErr)
				rErrHistChi2Down = rt.TH1F("rErrHistChi2Down","rErr, chi2 < {}".format(chi2Div),nBins,downLimitrErr,upLimitrErr)

				rmuHist = rt.TH1F("rmuHist","(r-{})/rErr, all Events".format(n),nBins,downLimitRMU,upLimitRMU)
				rmuHistChi2Up = rt.TH1F("rmuHistChi2Up","(r-{})/rErr, chi2 > {}".format(n,chi2Div),nBins,downLimitRMU,upLimitRMU)
				rmuHistChi2Down = rt.TH1F("rmuHistChi2Down","(r-{})/rErr, chi2 < {}".format(n,chi2Div),nBins,downLimitRMU,upLimitRMU)

				p1Hist = rt.TH1F("p1Hist","p1, all events", nBins,parLimDown1,parLimUp1)
				p1HistChi2Up = rt.TH1F("p1HistChi2Up","p1, chi2 > {}".format(chi2Div), nBins,parLimDown1,parLimUp1)
				p1HistChi2Down = rt.TH1F("p1HistChi2Down","p1, chi2 < {}".format(chi2Div), nBins,parLimDown1,parLimUp1)

				p2Hist = rt.TH1F("p2Hist","p2, all events", nBins,parLimDown2,parLimUp2)
				p2HistChi2Up = rt.TH1F("p2HistChi2Up","p2, chi2 > {}".format(chi2Div), nBins,parLimDown2,parLimUp2)
				p2HistChi2Down = rt.TH1F("p2HistChi2Down","p2, chi2 < {}".format(chi2Div), nBins,parLimDown2,parLimUp2)

				p3Hist = rt.TH1F("p3Hist","p3, all events", nBins,parLimDown3,parLimUp3)
				p3HistChi2Up = rt.TH1F("p3HistChi2Up","p3, chi2 > {}".format(chi2Div), nBins,parLimDown3,parLimUp3)
				p3HistChi2Down = rt.TH1F("p3HistChi2Down","p3, chi2 < {}".format(chi2Div), nBins,parLimDown3,parLimUp3)

				p4Hist = rt.TH1F("p4Hist","p4, all events", nBins,parLimDown4,parLimUp4)
				p4HistChi2Up = rt.TH1F("p4HistChi2Up","p4, chi2 > {}".format(chi2Div), nBins,parLimDown4,parLimUp4)
				p4HistChi2Down = rt.TH1F("p4HistChi2Down","p4, chi2 < {}".format(chi2Div), nBins,parLimDown4,parLimUp4)

				p5Hist = rt.TH1F("p5Hist","p5, all events", nBins,parLimDown5,parLimUp5)
				p5HistChi2Up = rt.TH1F("p5HistChi2Up","p5, chi2 > {}".format(chi2Div), nBins,parLimDown5,parLimUp5)
				p5HistChi2Down = rt.TH1F("p5HistChi2Down","p5, chi2 < {}".format(chi2Div), nBins,parLimDown5,parLimUp5)

				chi2Hist = rt.TH1F("chi2Hist","chi2, all events", nBins, chi2LimitDown, chi2LimitUp)
				chi2HistChi2Up = rt.TH1F("chi2HistChi2Up","chi2, chi2 > {}".format(chi2Div), nBins, chi2LimitDown, chi2LimitUp)
				chi2HistChi2Down = rt.TH1F("chi2HistChi2Down","chi2, chi2 < {}".format(chi2Div), nBins, chi2LimitDown, chi2LimitUp)

				chiRHist = rt.TH1F("chiRHist","reduced chi2, all events", nBins, 0, 10)
				chiRHistChi2Up = rt.TH1F("chiRHistChi2Up","reduced chi2, chi2 > {}".format(chi2Div), nBins, 0, 10)
				chiRHistChi2Down = rt.TH1F("chiRHistChi2Down","reduced chi2, chi2 < {}".format(chi2Div), nBins, 0, 10)

				eitHist = rt.TH1F("eitHist","Events in Toy, all events", nBins, eitLimitDown, eitLimitUp)
				eitHistChi2Up = rt.TH1F("eitHistChi2Up","eit, chi2 > {}".format(chi2Div), nBins, eitLimitDown, eitLimitUp)
				eitHistChi2Down = rt.TH1F("eitHistChi2Down","eit, chi2 < {}".format(chi2Div), nBins, eitLimitDown, eitLimitUp)

				#p1vRHist = rt.TH2F("p1vRHist","p1, all events",nBins,rDownLim,rUpLim, nBins,parLimDown1,parLimUp1)
				#p2vRHist = rt.TH2F("p2vRHist","p2, all events",nBins,rDownLim,rUpLim, nBins,parLimDown2,parLimUp2)
				#p3vRHist = rt.TH2F("p3vRHist","p3, all events",nBins,rDownLim,rUpLim, nBins,parLimDown3,parLimUp3)
				#p4vRHist = rt.TH2F("p4vRHist","p4, all events",nBins,rDownLim,rUpLim, nBins,parLimDown4,parLimUp4)
				#p5vRHist = rt.TH2F("p5vRHist","p5, all events",nBins,rDownLim,rUpLim, nBins,parLimDown5,parLimUp5)

				rHist.SetLineColor(rt.kBlack)
				rHistChi2Up.SetLineColor(rt.kBlue)
				rHistChi2Down.SetLineColor(rt.kRed)
				rErrHist.SetLineColor(rt.kBlack)
				rErrHistChi2Up.SetLineColor(rt.kBlue)
				rErrHistChi2Down.SetLineColor(rt.kRed)
				rmuHist.SetLineColor(rt.kBlack)
				rmuHistChi2Up.SetLineColor(rt.kBlue)
				rmuHistChi2Down.SetLineColor(rt.kRed)
				p1Hist.SetLineColor(rt.kBlack)
				p1HistChi2Up.SetLineColor(rt.kBlue)
				p1HistChi2Down.SetLineColor(rt.kRed)
				p2Hist.SetLineColor(rt.kBlack)
				p2HistChi2Up.SetLineColor(rt.kBlue)
				p2HistChi2Down.SetLineColor(rt.kRed)
				p3Hist.SetLineColor(rt.kBlack)
				p3HistChi2Up.SetLineColor(rt.kBlue)
				p3HistChi2Down.SetLineColor(rt.kRed)
				p4Hist.SetLineColor(rt.kBlack)
				p4HistChi2Up.SetLineColor(rt.kBlue)
				p4HistChi2Down.SetLineColor(rt.kRed)
				p5Hist.SetLineColor(rt.kBlack)
				p5HistChi2Up.SetLineColor(rt.kBlue)
				p5HistChi2Down.SetLineColor(rt.kRed)
				chi2Hist.SetLineColor(rt.kBlack)
				chi2HistChi2Up.SetLineColor(rt.kBlue)
				chi2HistChi2Down.SetLineColor(rt.kRed)
				chiRHist.SetLineColor(rt.kBlack)
				chiRHistChi2Up.SetLineColor(rt.kBlue)
				chiRHistChi2Down.SetLineColor(rt.kRed)
				eitHist.SetLineColor(rt.kBlack)
				eitHistChi2Up.SetLineColor(rt.kBlue)
				eitHistChi2Down.SetLineColor(rt.kRed)

				rHist.SetLineWidth(2)
				rHistChi2Up.SetLineWidth(2)
				rHistChi2Down.SetLineWidth(2)
				rErrHist.SetLineWidth(2)
				rErrHistChi2Up.SetLineWidth(2)
				rErrHistChi2Down.SetLineWidth(2)
				rmuHist.SetLineWidth(2)
				rmuHistChi2Up.SetLineWidth(2)
				rmuHistChi2Down.SetLineWidth(2)
				p1Hist.SetLineWidth(2)
				p1HistChi2Up.SetLineWidth(2)
				p1HistChi2Down.SetLineWidth(2)
				p2Hist.SetLineWidth(2)
				p2HistChi2Up.SetLineWidth(2)
				p2HistChi2Down.SetLineWidth(2)
				p3Hist.SetLineWidth(2)
				p3HistChi2Up.SetLineWidth(2)
				p3HistChi2Down.SetLineWidth(2)
				p4Hist.SetLineWidth(2)
				p4HistChi2Up.SetLineWidth(2)
				p4HistChi2Down.SetLineWidth(2)
				chi2Hist.SetLineWidth(2)
				chi2HistChi2Up.SetLineWidth(2)
				chi2HistChi2Down.SetLineWidth(2)
				chiRHist.SetLineWidth(2)
				chiRHistChi2Up.SetLineWidth(2)
				chiRHistChi2Down.SetLineWidth(2)
				eitHist.SetLineWidth(2)
				eitHistChi2Up.SetLineWidth(2)
				eitHistChi2Down.SetLineWidth(2)
				


				#snapshotArgList
				# var names are <region>_p<i>_<nPar>
				mdfParList = [0,0,0,0]
				mdfParListAlt = [0,0,0,0]
				for i in [1,2,3,4]:
					if i <= nPar:
						mdfParList[i-1] = snapshotArgList.getRealValue(region+"_p"+str(i)+"_"+str(nPar))
					if i <= nParAlt:
						mdfParListAlt[i-1] = snapshotArgList.getRealValue(region+"_p"+str(i)+"_"+str(nParAlt)+"_alt")

				pae+="{} {} {} {} {} {} {} {} {} {} {}\n".format(sigPars, region, expSig,mdfParList[0],mdfParList[1],mdfParList[2],mdfParList[3],mdfParListAlt[0],mdfParListAlt[1],mdfParListAlt[2],mdfParListAlt[3] )
					

				listOfFuncs = []
				nPass = 0
				nTotal = 0
				#canv_Plots = rt.TCanvas("canv_Plots","canv_Plots",1200,1200)
				#canv_Plots.SetLogy()
				#data.Draw()
				sigHist.SetMaximum(10000)
				#sigHist.Draw()
				toyToDraw = -1
				ndfTot = 0.
				nToysTot = 0.
				for iEvt in range(limit.GetEntries()):
					#print(iEvt)
					limit.GetEvent(iEvt)
					if limit.quantileExpected != -1:
						continue
					tree.GetEvent(limit.iToy)
					toy = fitOnlyFile.Get("toys/toy_{}".format(limit.iToy))
					#toy.Print()
					#toyExp = getattr(tree,"n_exp_final_bin"+region+"_2018_proc_roomultipdf")
					toyExp = toy.sumEntries()
					#toy.Draw("same")

					rHistVal = tree.r
					rErrHistVal = tree.rErr
					if rErrHistVal == 0:
						continue
					rmuHistVal = (tree.r-int(n))/tree.rErr
					if toyToDraw == -1 and rmuHistVal > 2:
						toyToDraw = limit.iToy
					if bigName[-3:] != "Alt": #for fitting at Main func
						if nPar == 2:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_1")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_1")
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [1]) * pow(x/13000, -[2])",1500,8000))
						elif nPar == 3:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_2")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_2")
							p3 = getattr(limit,"trackedParam_"+region+"_p3_2")
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [1]) * pow(x/13000, -([2]+[3]*log(x/13000)))",1500,8000))
						elif nPar == 5:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_4")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_4")
							p3 = getattr(limit,"trackedParam_"+region+"_p3_4")
							p4 = getattr(limit,"trackedParam_"+region+"_p4_4")
							p5 = getattr(limit,"trackedParam_"+region+"_p5_4")
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [1]) * pow(x/13000, -[2]-[3]*log(x/13000)-[4]*pow(log(x/13000),2)-[5]*pow(log(x/13000),3))",1500,8000))#func14
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [1]+[2]*log(x/13000)+[3]*pow(log(x/13000),2)) * pow(x/13000, -[4]-[5]*log(x/13000))",1500,8000))#func32
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [1]+[2]*log(x/13000)+[3]*pow(log(x/13000),2)+[4]*pow(log(x/13000),3)) * pow(x/13000, [5])",1500,8000))#func41
						else:
							print("nPar is not 2 nor 5!")
							exit(0)
					elif bigName[-3:] == "Alt": #for fitting with Alt Func
						if nParAlt == 1:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_1_alt")
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*x/13000)",1500,8000))
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*x/13000)",1500,8000))
						elif nParAlt == 2:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_2_alt")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_2_alt")
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*x/13000 + [2]*log(x/13000))",1500,8000))
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*(x/13000)) * pow(x/13000,[2])",1500,8000))
						elif nParAlt == 3:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_3_alt")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_3_alt")
							p3 = getattr(limit,"trackedParam_"+region+"_p3_3_alt")
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*x/13000 + [2]*log(x/13000) + [3]*pow(log(x/13000),2))",1500,8000))
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*(x/13000)) * pow(x/13000,[2]*(1+[3]*log(x/13000)))",1500,8000)) #reparam, decor
						elif nParAlt == 4:
							p1 = getattr(limit,"trackedParam_"+region+"_p1_4_alt")
							p2 = getattr(limit,"trackedParam_"+region+"_p2_4_alt")
							p3 = getattr(limit,"trackedParam_"+region+"_p3_4_alt")
							p4 = getattr(limit,"trackedParam_"+region+"_p4_4_alt")
							#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*x/13000 + [2]*log(x/13000) + [3]*pow(log(x/13000),2) + [4]*pow(log(x/13000),3))",1500,8000))
							listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * exp([1]*(x/13000)) * pow(x/13000,[2]*(1+[3]*log(x/13000)*(1+[4]*log(x/13000)))",1500,8000))# reparam, decor
						else:
							print("nParAlt is not 1, 2, 3, nor 4!")
							exit(0)

					if bigName[-3:] == "Alt": nPar = nParAlt
					listOfFuncs[-1].SetParameter(1,p1)
					if nPar >= 2:
						listOfFuncs[-1].SetParameter(2,p2)
					if nPar >= 3:
						listOfFuncs[-1].SetParameter(3,p3)
					if nPar >= 4:
						listOfFuncs[-1].SetParameter(4,p4)
					if nPar >= 5:
						listOfFuncs[-1].SetParameter(5,p5)
					# pdf's in RooFit do not have a normlization factor (i.e., Sum(allspace) of pdf = 1)
					# so we need to scale it to our dataset.
					# first, set norm to 1, and we know the number of events in the toy)
					# then scale the function by numEvents divided by intergral of function without a normialztion
					# i.e., first normalize to 1 (divide by integral), then scale to (multi by) numEvents
					# factor of 100 is becaue our bins are 100 GeV wide
					norm = 1
					listOfFuncs[-1].SetParameter(0,norm)
					denom = listOfFuncs[-1].Integral(1500,8000) 
					if denom > 0:
						norm = (toyExp-rHistVal*sigHist.GetEffectiveEntries())*100/(denom)
					else:
						continue
					listOfFuncs[-1].SetParameter(0,norm)
					listOfFuncs[-1].Draw("same")
					chi2 = 0
					ndf = -nPar-2
					print(ndf)
					#ndf = 0
					for iBin in range(toy.numEntries()):
						toy.get(iBin)
						x = toy.get(iBin).getRealValue("mH"+region+"_2018")
						Ei = listOfFuncs[-1].Eval(x) + rHistVal*sigHist.GetBinContent(iBin)
						#Ei = listOfFuncs[-1].Eval(x)
						Oi = toy.weight()
						#print(x, Ei, Oi)
						if (Ei > 0.) and (Oi > 0.):#Ei <= 0 means divide by 0, Oi <= 0 means bin is empty and should be excluded
							ndf+=1
							print(ndf, Ei, Oi, x)
							chi2 += ((Oi-Ei)**2)/Ei
					if ndf < 0:
						continue
					ndfTot += ndf
					nToysTot += 1
					chiR = chi2/float(ndf)
					#ndf -= (nPar - 2)
					#print(ndf, toy.numEntries())

					chi2Hist.Fill(chi2)
					chiRHist.Fill(chiR)
					#print(chi2, ndf, chi2/ndf)

					rHist.Fill(rHistVal)
					rErrHist.Fill(rErrHistVal)
					rmuHist.Fill(rmuHistVal)
					p1Hist.Fill(p1)
					#p1vRHist.Fill(rHistVal,p1)
					if nPar >= 2:
						p2Hist.Fill(p2)
						#p2vRHist.Fill(rHistVal,p2)
					if nPar >= 3:
						p3Hist.Fill(p3)
						#p3vRHist.Fill(rHistVal,p3)
					if nPar >= 4:
						p4Hist.Fill(p4)
						#p4vRHist.Fill(rHistVal,p4)
					if nPar >= 5:
						p5Hist.Fill(p5)
						#p5vRHist.Fill(rHistVal,p5)
					eitHist.Fill(toyExp)

					#varCheck = (rt.Math.chisquared_cdf_c(chi2, ndf) > 0.05) # pass means chi2 is lower than crit value of alpha = 0.05
					varCheck = (rt.Math.chisquared_cdf_c(chi2, 65-2-nPar) > 0.05) # pass means chi2 is lower than crit value of alpha = 0.05
					legText = "#chi^2"
					#if region == "highSVJ2":
					#	varCheck = (rmuHistVal > 2)
					#	legText = "FOM"
					nTotal += 1
					nPass += int(varCheck)
					#varCheck = (chi2/ndf < chi2Div)

					rHistChi2Up.Fill(rHistVal, varCheck)
					rErrHistChi2Up.Fill(rErrHistVal, varCheck)
					rmuHistChi2Up.Fill(rmuHistVal, varCheck)
					p1HistChi2Up.Fill(p1, varCheck)
					if nPar >= 2:
						p2HistChi2Up.Fill(p2, varCheck)
					if nPar >= 3:
						p3HistChi2Up.Fill(p3, varCheck)
					if nPar >= 4:
						p4HistChi2Up.Fill(p4, varCheck)
					if nPar >= 5:
						p5HistChi2Up.Fill(p5, varCheck)
					chi2HistChi2Up.Fill(chi2, varCheck)
					chiRHistChi2Up.Fill(chiR, varCheck)
					eitHistChi2Up.Fill(toyExp, varCheck)

					rHistChi2Down.Fill(rHistVal, not varCheck)
					rErrHistChi2Down.Fill(rErrHistVal, not varCheck)
					rmuHistChi2Down.Fill(rmuHistVal, not varCheck)
					p1HistChi2Down.Fill(p1, not varCheck)
					if nPar >= 2:
						p2HistChi2Down.Fill(p2, not varCheck)
					if nPar >= 3:
						p3HistChi2Down.Fill(p3, not varCheck)
					if nPar >= 4:
						p4HistChi2Down.Fill(p4, not varCheck)
					if nPar >= 5:
						p5HistChi2Down.Fill(p5, not varCheck)
					chi2HistChi2Down.Fill(chi2, not varCheck)
					chiRHistChi2Down.Fill(chiR, not varCheck)
					eitHistChi2Down.Fill(toyExp, not varCheck)
				#canv_Plots.SaveAs("plots/manyFits_"+region+"_"+expSig+".png")
				#try: ndf = ndfTot/nToysTot
				#ndf = ndfTot/nToysTot
				ndf = 65-2-nPar
				#except ZeroDivisionError: continue 
				#ndf = 65 - nPar - 2
				#pae += "\n {} {} {} {} {}".format(sigPars, region, expSig, nTotal, nPass)
				c1 = rt.TCanvas("c{}".format(i),"c{}".format(i),1500,1500)
				i+=1
				c1.Divide(3,3)

				iPlot1 = 5
				iPlot2 = 6
				iPlot3 = 7
				iPlot4 = 8
				iPlot5 = 9
				if nPar < 4:
					doLeg = True
				else:
					doLeg = False
				c1.cd(1)
				rHist.Draw("hist")
				rHistChi2Up.Draw("hist same")
				rHistChi2Down.Draw("hist same")
				c1.cd(2)
				rErrHist.Draw("hist")
				rErrHistChi2Up.Draw("hist same")
				rErrHistChi2Down.Draw("hist same")
				c1.cd(3)
				rmuHist.Draw("hist")
				rmuHistChi2Up.Draw("hist same")
				rmuHistChi2Down.Draw("hist same")
				#c1.cd(4)
				#eitHist.Draw("hist")
				#eitHistChi2Up.Draw("hist same")
				#eitHistChi2Down.Draw("hist same")
				c1.cd(iPlot1)
				p1Hist.Draw("hist")
				p1HistChi2Up.Draw("hist same")
				p1HistChi2Down.Draw("hist same")
				if nPar >= 2:
					c1.cd(iPlot2)
					p2Hist.Draw("hist")
					p2HistChi2Up.Draw("hist same")
					p2HistChi2Down.Draw("hist same")
				if nPar >= 3:
					c1.cd(iPlot3)
					p3Hist.Draw("hist")
					p3HistChi2Up.Draw("hist same")
					p3HistChi2Down.Draw("hist same")
				if nPar >= 4:
					c1.cd(iPlot4)
					p4Hist.Draw("hist")
					p4HistChi2Up.Draw("hist same")
					p4HistChi2Down.Draw("hist same")
				if nPar >= 5:
					c1.cd(iPlot5)
					p5Hist.Draw("hist")
					p5HistChi2Up.Draw("hist same")
					p5HistChi2Down.Draw("hist same")
				if nPar < 3:
					c1.cd(7)
					chiRHist.Draw("hist")
					chiRHistChi2Up.Draw("hist same")
					chiRHistChi2Down.Draw("hist same")
					c1.cd(8)
					rt.gPad.SetLogy()
					dataHist.Draw()
					if len(listOfFuncs) > 0:
						try: listOfFuncs[toyToDraw].Draw("same")
						except IndexError: pass
				c1.cd(4)
				chi2Hist.Draw("hist")
				chi2HistChi2Up.Draw("hist same")
				chi2HistChi2Down.Draw("hist same")
				chi2Func_test = rt.TF1("chi2Func","ROOT::Math::chisquared_pdf(x,{},0)".format(ndf),chi2LimitDown,chi2LimitUp)
				chi2Norm = chi2Hist.GetBinContent(chi2Hist.GetMaximumBin())/chi2Func_test.GetMaximum()
				#chi2Func = rt.TF1("chi2Func","{}*({}/{}) * ROOT::Math::chisquared_pdf(x,{},0)".format(50,10,nBins,ndf,ndf),chi2LimitDown,chi2LimitUp)
				chi2Func = rt.TF1("chi2Func","{}*ROOT::Math::chisquared_pdf(x,{},0)".format(chi2Norm,ndf),chi2LimitDown,chi2LimitUp)
				chi2Func.SetLineColor(rt.kGreen)
				chi2Func.Draw("same")
				if doLeg:
					c1.cd(9)
					leg = rt.TLegend(0.0,0.0,1.0,1.0)
					chi2Func.Draw("axis")
					leg.SetHeader("Legend","C")
					leg.AddEntry(rHist,"All Events","l")
					leg.AddEntry(rHistChi2Up,"Low "+legText+" Events","l")
					leg.AddEntry(rHistChi2Down,"High "+legText+" Events","l")
					leg.AddEntry(chi2Func, "Ideal #chi^2 Distribtuion","l")
					leg.Draw()
				#c1.SaveAs("../condorTests/chi2Dist/"+SVJNAME+"_"+region+expSig+".png")
				c1.SaveAs("plots/"+SVJNAME+"_"+bigName+"_chi2Diagnostic.png")
				c1.Delete()
				#c3 = rt.TCanvas("c3","c3",1500,1500)
				#c3.Divide(2,3)
				#c3.cd(1)
				#p1vRHist.Draw("colz")
				#c3.cd(2)
				#p2vRHist.Draw("colz")
				#c3.cd(3)
				#p3vRHist.Draw("colz")
				#c3.cd(4)
				#p4vRHist.Draw("colz")
				#c3.cd(5)
				#p5vRHist.Draw("colz")
				#c3.SaveAs("plots/"+SVJNAME+"_"+bigName+"_2D_paramVr.png")
				#c2 = rt.TCanvas("c{}".format(i),"c{}".format(i),1500,1500)
				#i += 1
				#rmuHist.Draw("hist")
				#if rmuHist.GetEntries() > 0:
				#	rmuHist.Fit("gaus",'','PLC',-5,5)
				#	rmuHist.GetFunction("gaus").SetLineColor(rt.kBlack)
				#	rmuHist.GetFunction("gaus").SetLineStyle(9)
				#	rmuHist.GetFunction("gaus").SetLineWidth(8)
				#	outputROOTFile.Get("hist_"+rGenCode+"_"+regCode).Fill(sigPars[0],rmuHist.GetFunction("gaus").GetParameter(1))
				#	pae2 += "{} {} {} {} {} {} {} {} {}\n".format(region, sigPars[0],sigPars[1],sigPars[2],sigPars[3],expSig, rmuHist.GetFunction("gaus").GetParameter(0),rmuHist.GetFunction("gaus").GetParameter(1),rmuHist.GetFunction("gaus").GetParameter(2) )
				#rmuHist.Draw("func same")
				#rmuHistChi2Up.Draw("hist same")
				#rmuHistChi2Down.Draw("hist same")
				#c2.SaveAs("../Oct1/plots/"+SVJNAME+"_"+bigName+"_pullOnly.png")
				#c2.Delete()
				fitDiagFile.Close()
				fitOnlyFile.Close()
				mdfFile.Close()
				#dataObsFile.Close()
		wsFile.Close()
"""
if makeNew:
	outputROOTFile.cd()
	hist_0M_lc.Write()
	hist_1M_lc.Write()
	hist_0A_lc.Write()
	hist_1A_lc.Write()
	hist_0M_l2.Write()
	hist_1M_l2.Write()
	hist_0A_l2.Write()
	hist_1A_l2.Write()
	hist_0M_hc.Write()
	hist_1M_hc.Write()
	hist_0A_hc.Write()
	hist_1A_hc.Write()
	hist_0M_h2.Write()
	hist_1M_h2.Write()
	hist_0A_h2.Write()
	hist_1A_h2.Write()

outputROOTFile.Write()
outputROOTFile.Close()
"""


#print(pae) 
#print(pae2)


















