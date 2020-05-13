import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function

#regions = ["lowCut","lowSVJ1","lowSVJ2","highCut","highSVJ1","highSVJ2"]
regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
#regions = ["highCut"]
# Mar30_modified
# Region Thry/Main Dijet/Alt
# highCut 2 3
# highSVJ2 1 1
# lowCut 2 2
# lowSVJ2 2 2

for region in regions:
	for toyFunc in ["Main","Alt"]:
		for sigStr in ["Sig0"]:
			for fitFunc in ["Main"]:
				halfName = toyFunc+sigStr
				fullName = toyFunc+sigStr+fitFunc
				n = sigStr[3]
				fitDiagFile = rt.TFile.Open("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/fitDiagnostics"+fullName+".root","read")
				fitOnlyFile = rt.TFile.Open("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+fullName+".FitDiagnostics.mH120.123456.root","read")
				genOnlyFile = rt.TFile.Open("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+halfName+".GenerateOnly.mH125.123456.root","read")
				dataObsFile = rt.TFile.Open("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","read")
				svjWs = dataObsFile.Get("SVJ")
				data = rt.RooDataHist()
				data = svjWs.data("data_obs")
				if toyFunc == "Main":
					genPdf = svjWs.pdf("Bkg_"+region+"_2018")
				else:
					genPdf = svjWs.pdf("Bkg_Alt_"+region+"_2018")
				tree = fitDiagFile.Get("tree_fit_sb")
				limit = fitOnlyFile.Get("limit")
				
				chi2Div = 1.55
				#if region == "highSVJ1" and toyFunc == "Alt":
				#	chi2Div = -10
				#if region == "highSVJ2" and toyFunc == "Main":
				#	chi2Div = -2
				#if region == "lowSVJ1" and toyFunc == "Main":
				#	chi2Div = -2

				if ((region == "highCut") and (toyFunc == "Alt")):
					rDownLim = 0.9
					rUpLim = 1.1
					downLimitrErr = 2.82
					upLimitrErr = 2.83
					downLimitRMU = 0.3
					upLimitRMU = 0.4
					parLimDown1 = 5.8
					parLimUp1 = 7.3
					parLimDown2 = 7.37
					parLimUp2 = 11.438
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "highCut") and (toyFunc == "Main")):
					rDownLim = 0.9
					rUpLim = 1.1
					downLimitrErr = 2.82
					upLimitrErr = 2.83
					downLimitRMU = 0.3 
					upLimitRMU = 0.4
					parLimDown1 = 6.3
					parLimUp1 = 6.31
					parLimDown2 = 11.437
					parLimUp2 = 11.438
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "highSVJ2") and (toyFunc == "Alt")):
					rDownLim = 0
					rUpLim = 20
					downLimitrErr = 0
					upLimitrErr = 35
					downLimitRMU = 0 
					upLimitRMU = 9
					parLimDown1 = 7
					parLimUp1 = 7.4
					parLimDown2 = -1000
					parLimUp2 = 1000
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "highSVJ2") and (toyFunc == "Main")):
					rDownLim = -7
					rUpLim = 4
					downLimitrErr = 1
					upLimitrErr = 3
					downLimitRMU = -4 
					upLimitRMU = 4
					parLimDown1 = 6.8
					parLimUp1 = 7.2
					parLimDown2 = -1000
					parLimUp2 = 1000
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "lowCut") and (toyFunc == "Alt")):
					rDownLim = -1
					rUpLim = 3
					downLimitrErr = 2.75
					upLimitrErr = 2.85
					downLimitRMU = -0.5 
					upLimitRMU = 1.2
					parLimDown1 = 5.8
					parLimUp1 = 5.81
					parLimDown2 = 7.36
					parLimUp2 = 7.4
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "lowCut") and (toyFunc == "Main")):
					rDownLim = 1
					rUpLim = 4
					downLimitrErr = 16.4
					upLimitrErr = 16.8
					downLimitRMU = 0 
					upLimitRMU = 0.3
					parLimDown1 = 5.6
					parLimUp1 = 6
					parLimDown2 = 7.36
					parLimUp2 = 7.4
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "lowSVJ2") and (toyFunc == "Alt")):
					rDownLim = -1
					rUpLim = 2
					downLimitrErr = 2.81
					upLimitrErr = 2.83
					downLimitRMU = -0.5
					upLimitRMU = 1
					parLimDown1 = 5.8
					parLimUp1 = 6.3
					parLimDown2 = 10
					parLimUp2 = 12
					chi2LimitDown = 0
					chi2LimitUp = 30
				elif ((region == "lowSVJ2") and (toyFunc == "Main")):
					rDownLim = -1
					rUpLim = 5
					downLimitrErr = 6
					upLimitrErr = 10
					downLimitRMU = -0.5 
					upLimitRMU = 1
					parLimDown1 = 5.8
					parLimUp1 = 6.4
					parLimDown2 = 10
					parLimUp2 = 12
					chi2LimitDown = 0
					chi2LimitUp = 30
				else:
					rDownLim = -20
					rUpLim = -10
					upLimitrErr = 40
					downLimitrErr = 0
					upLimitRMU = 40 
					downLimitRMU = -40 
					parLimDown1 = 0.5
					parLimUp1 = 4
					parLimDown2 = -50
					parLimUp2 = -30
					chi2LimitDown = 0
					chi2LimitUp = 30

				
				nBins = 50
				rHist = rt.TH1F("rHist","r, all Events",nBins,rDownLim,rUpLim)
				rHistChi2Up = rt.TH1F("rHistChi2Up","r, chi2 > {}".format(chi2Div),nBins,rDownLim,rUpLim)
				rHistChi2Down = rt.TH1F("rHistChi2Down","r, chi2 < {}".format(chi2Div),nBins,rDownLim,rUpLim)
				rErrHist = rt.TH1F("rErrHist","rErr, all Events",nBins,downLimitrErr,upLimitrErr)
				rErrHistChi2Up = rt.TH1F("rErrHistChi2Up","rErr, chi2 > {}".format(chi2Div),nBins,downLimitrErr,upLimitrErr)
				rErrHistChi2Down = rt.TH1F("rErrHistChi2Down","rErr, chi2 < {}".format(chi2Div),nBins,downLimitrErr,upLimitrErr)
				rmuHist = rt.TH1F("rmuHist","(r-mu)/rErr, all Events",nBins,downLimitRMU,upLimitRMU)
				rmuHistChi2Up = rt.TH1F("rmuHistChi2Up","(r-mu)/rErr, chi2 > {}".format(chi2Div),nBins,downLimitRMU,upLimitRMU)
				rmuHistChi2Down = rt.TH1F("rmuHistChi2Down","(r-mu)/rErr, chi2 < {}".format(chi2Div),nBins,downLimitRMU,upLimitRMU)
				p1Hist = rt.TH1F("p1Hist","p1, all events", nBins,parLimDown1,parLimUp1)
				p1HistChi2Up = rt.TH1F("p1HistChi2Up","p1, chi2 > {}".format(chi2Div), nBins,parLimDown1,parLimUp1)
				p1HistChi2Down = rt.TH1F("p1HistChi2Down","p1, chi2 < {}".format(chi2Div), nBins,parLimDown1,parLimUp1)
				p2Hist = rt.TH1F("p2Hist","p2, all events", nBins,parLimDown2,parLimUp2)
				p2HistChi2Up = rt.TH1F("p2HistChi2Up","p2, chi2 > {}".format(chi2Div), nBins,parLimDown2,parLimUp2)
				p2HistChi2Down = rt.TH1F("p2HistChi2Down","p2, chi2 < {}".format(chi2Div), nBins,parLimDown2,parLimUp2)
				chi2Hist = rt.TH1F("chi2Hist","chi2, all events", nBins, chi2LimitDown, chi2LimitUp)
				chi2HistChi2Up = rt.TH1F("chi2HistChi2Up","chi2, chi2 > {}".format(chi2Div), nBins, chi2LimitDown, chi2LimitUp)
				chi2HistChi2Down = rt.TH1F("chi2HistChi2Down","chi2, chi2 < {}".format(chi2Div), nBins, chi2LimitDown, chi2LimitUp)

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
				chi2Hist.SetLineColor(rt.kBlack)
				chi2HistChi2Up.SetLineColor(rt.kBlue)
				chi2HistChi2Down.SetLineColor(rt.kRed)

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
				chi2Hist.SetLineWidth(2)
				chi2HistChi2Up.SetLineWidth(2)
				chi2HistChi2Down.SetLineWidth(2)

				hist2D_chi2_rmu = rt.TH2F("hist2D_chi2_rmu","Chi2VsRMuErr;chi2;rmuerr",nBins, 0, 10, nBins, downLimitRMU,10)
				
				if "low" in region:
					nPar = 2
				else:
					nPar = 1
				
				listOfFuncs = []
				for iEvt in range(limit.GetEntries()):
					limit.GetEvent(iEvt)
					if limit.quantileExpected != -1:
						continue
					tree.GetEvent(limit.iToy)
					toy = genOnlyFile.Get("toys/toy_{}".format(limit.iToy))
					toyExp = getattr(tree,"n_exp_final_bin"+region+"_2018_proc_roomultipdf")

					rHistVal = tree.r
					rErrHistVal = tree.rErr
					rmuHistVal = (tree.r-int(n))/tree.rErr
					#if region != "highSVJ2":
					if nPar == 1:
						p1 = getattr(limit,"trackedParam_"+region+"_p1_1")
						listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(x, -[1])",1500,6000))
					elif nPar == 2:
						p1 = getattr(limit,"trackedParam_"+region+"_p1_2")
						p2 = getattr(limit,"trackedParam_"+region+"_p2_2")
						listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x, [2]) * pow(x, -[1])",1500,6000))
					else:
						print("nPar is not 1 or 2!")
						exit(0)
 
					listOfFuncs[-1].SetParameter(1,p1)
					if nPar >= 2:
						listOfFuncs[-1].SetParameter(2,p2)
					# pdf's in RooFit do not have a normlization facotr (i.e., Sum(allspace) of pdf = 1)
					# so we need to scale it to our dataset.
					# first, set norm to 1, and we know the number of events in the toy)
					# then scale the function by numEvents divided by intergral of function without a normialztion
					# i.e., first normalize to 1 (divide by integral), then scale to numEvents
					# factor of 100 is becaue our bins are 100 GeV wide
					norm = 1
					listOfFuncs[-1].SetParameter(0,norm)
					listOfFuncs[-1].SetParameter(0,toyExp/listOfFuncs[-1].Integral(1500,6000)*50)
					chi2 = 0
					for iBin in range(toy.numEntries()):
						toy.get(iBin)
						x = toy.get(iBin).getRealValue("mH")
						Oi = toy.weight()
						Ei = listOfFuncs[-1].Eval(x)
						print(x, Oi, Ei)
						try:
							chi2 += ((Oi-Ei)**2)/Ei
						except ZeroDivisionError:
							chi2 += 0
					ndf = 90 - nPar
					chi2Hist.Fill(chi2/ndf)
					#else:
					#	p1 = getattr(limit,"trackedParam_"+region+"_p1_1")
					#	p2 = -9
					#	listOfFuncs.append(rt.TF1("iToy_"+str(iEvt-1),"[1]*pow(1 - x/13000, [0])",1500,3800))
					#	listOfFuncs[-1].SetParameter(0,p1)
					#	listOfFuncs[-1].SetParameter(1,1)
					#	# pdf's in RooFit do not have a normlization facotr (i.e., Sum(allspace) of pdf = 1)
					#	# so we need to scale it to our dataset.
					#	# first (above line), set norm to 1, and we know the number of events in the toy)
					#	# then scale the function by numEvents divided by intergral of function without a normialztion
					#	# i.e., first normalize to 1 (divide by integral), then scale to numEvents
					#	# factor of 100 is becaue our bins are 100 GeV wide
					#	listOfFuncs[-1].SetParameter(1,toyExp/listOfFuncs[-1].Integral(1500,3800)*100)
					#	chi2 = 0
					#	for iBin in range(toy.numEntries()):
					#		toy.get(iBin)
					#		x = toy.get(iBin).getRealValue("mH")
					#		Oi = toy.weight()
					#		Ei = listOfFuncs[-1].Eval(x)
					#		chi2 += ((Oi-Ei)**2)/Ei
					#	ndf = 22
					#	chi2Hist.Fill(chi2/ndf)
					chi2Func = rt.TF1("chi2Func","{}*({}/{}) * ROOT::Math::chisquared_pdf(x,{},0)/{}".format(50,10,nBins,ndf,ndf),0,10)
					chi2Func.SetLineColor(rt.kGreen)

					rHist.Fill(rHistVal)
					rErrHist.Fill(rErrHistVal)
					rmuHist.Fill(rmuHistVal)
					p1Hist.Fill(p1)
					p2Hist.Fill(p2)

					hist2D_chi2_rmu.Fill(chi2/ndf,rmuHistVal)

					varCheck = (chi2/ndf < chi2Div)

					rHistChi2Up.Fill(rHistVal, varCheck)
					rErrHistChi2Up.Fill(rErrHistVal, varCheck)
					rmuHistChi2Up.Fill(rmuHistVal, varCheck)
					p1HistChi2Up.Fill(p1, varCheck)
					p2HistChi2Up.Fill(p2, varCheck)
					chi2HistChi2Up.Fill(chi2/ndf, varCheck)

					rHistChi2Down.Fill(rHistVal, not varCheck)
					rErrHistChi2Down.Fill(rErrHistVal, not varCheck)
					rmuHistChi2Down.Fill(rmuHistVal, not varCheck)
					p1HistChi2Down.Fill(p1, not varCheck)
					p2HistChi2Down.Fill(p2, not varCheck)
					chi2HistChi2Down.Fill(chi2/ndf, not varCheck)

				c1 = rt.TCanvas("c1","c1",2500,1000)
				c1.Divide(3,2)
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
				c1.cd(4)
				p1Hist.Draw("hist")
				p1HistChi2Up.Draw("hist same")
				p1HistChi2Down.Draw("hist same")
				if not (region == "highSVJ2"):
					c1.cd(5)
					p2Hist.Draw("hist")
					p2HistChi2Up.Draw("hist same")
					p2HistChi2Down.Draw("hist same")
				c1.cd(6)
				chi2Hist.Draw("hist")
				chi2HistChi2Up.Draw("hist same")
				chi2HistChi2Down.Draw("hist same")
				#chi2Func.Draw("same")
				c1.SaveAs("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/AA_"+region+"_"+fullName+"_10plots.png")
	
				c2 = rt.TCanvas("c2","c2",1000,1000)
				hist2D_chi2_rmu.Draw("colz")
				c2.SaveAs("cards_Mar30_modified/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/BB_"+region+"_"+fullName+"_2dPlot.png")

				fitDiagFile.Close()
				genOnlyFile.Close()
				dataObsFile.Close()
























