import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function

regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
#regions = ["lowSVJ1","lowSVJ2","highSVJ1","highSVJ2"]
#regions = ["lowSVJ2"]

# Mar3_ThryDij
# Region Thry/Main Dijet/Alt
# highCut 2 3
# highSVJ2 1 1
# lowCut 2 2
# lowSVJ2 2 2

# Generate 4 plots:
#	distribution of param1, for all and for qunatExp = 0.5
#	distribution of param2, for all and for quantExp = 0.5

copyToExcel = "Region & Toy & Sig & p1 & p2\n"

for region in regions:
	for toyFunc in ["Main","Alt"]:
		for sigStr in ["Sig0"]:
			for fitFunc in ["Main"]:
				halfName = toyFunc+sigStr
				fullName = toyFunc+sigStr+fitFunc
				n = sigStr[3]
				fitDiagFile = rt.TFile.Open("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/fitDiagnostics"+fullName+".root","read")
				fitOnlyFile = rt.TFile.Open("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+fullName+".FitDiagnostics.mH120.123456.root","read")
				genOnlyFile = rt.TFile.Open("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+halfName+".GenerateOnly.mH125.123456.root","read")
				dataObsFile = rt.TFile.Open("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","read")
				svjWs = dataObsFile.Get("SVJ")
				toysDirect = genOnlyFile.Get("toys")
				#print(type(svjWs),type(toysDirect))
				#toy = rt.RooDataHist()
				#data = rt.RooDataHist()
				data = svjWs.data("data_obs")
				#norm = svjWs.var("roomultipdf_norm").getValV()
				dataHist = data.createHistogram("mH")
				fitPdf = svjWs.pdf("Bkg_"+region+"_2018")
				if float(fitDiagFile.GetSize()) < 20000:
					print("="*25)
					print("{} {} could not be drawn".format(region, fullName))
					print("="*25)
					continue
				tree = fitDiagFile.Get("tree_fit_sb")
				limit = fitOnlyFile.Get("limit")
				c1 = rt.TCanvas("c1","c1",1500,1000)
				c1.Divide(3,2)
				quantCut = "quantileExpected==-1"
				if region != "highSVJ2":
					c1.cd(1)
					#limit.Draw("trackedParam_"+region+"_p1_2>>h1(100,-200,100)","")
					#rt.gDirectory.Get("h1").SetLineColor(rt.kBlack)
					#rt.gDirectory.Get("h1").SetLineWidth(2)
					limit.Draw("trackedParam_"+region+"_p1_2>>h2(100,-200,100)",quantCut)
					rt.gDirectory.Get("h2").SetLineColor(rt.kRed)
					rt.gDirectory.Get("h2").SetLineWidth(2)
					rt.gDirectory.Get("h2").SetAxisRange(0,rt.gDirectory.Get("h2").GetMaximum()*1.1,"Y")
					rt.gPad.Update()
					c1.cd(2)
					#limit.Draw("trackedParam_"+region+"_p2_2>>h3(100,-500,500)","")
					#rt.gDirectory.Get("h3").SetLineColor(rt.kBlack)
					#rt.gDirectory.Get("h3").SetLineWidth(2)
					limit.Draw("trackedParam_"+region+"_p2_2>>h4(100,-500,500)",quantCut)
					rt.gDirectory.Get("h4").SetLineColor(rt.kRed)
					rt.gDirectory.Get("h4").SetLineWidth(2)
					rt.gDirectory.Get("h4").SetAxisRange(0,rt.gDirectory.Get("h4").GetMaximum()*1.1,"Y")
					rt.gPad.Update()
				else:
					c1.cd(1)
					#limit.Draw("trackedParam_"+region+"_p1_1>>h1","")
					#rt.gDirectory.Get("h1").SetLineColor(rt.kBlack)
					#rt.gDirectory.Get("h1").SetLineWidth(2)
					limit.Draw("trackedParam_"+region+"_p1_1>>h2",quantCut)
					rt.gDirectory.Get("h2").SetLineColor(rt.kRed)
					rt.gDirectory.Get("h2").SetLineWidth(2)
					rt.gDirectory.Get("h2").SetAxisRange(0,rt.gDirectory.Get("h1").GetMaximum()*1.1,"Y")
					rt.gPad.Update()
				#c1.cd(3)
				#limit.Draw("trackedParam_"+region+"_p3_3>>h5(100,-100,100)","")
				#rt.gDirectory.Get("h5").SetLineColor(rt.kBlack)
				#rt.gDirectory.Get("h5").SetLineWidth(2)
				#limit.Draw("trackedParam_"+region+"_p3_3>>h6(100,-100,100)",quantCut,"same")
				#rt.gDirectory.Get("h6").SetLineColor(rt.kRed)
				#rt.gDirectory.Get("h6").SetLineWidth(2)
				#rt.gDirectory.Get("h5").SetAxisRange(0,rt.gDirectory.Get("h5").GetMaximum()*1.1,"Y")
				#rt.gPad.Update()
				#c1.cd(4)
				#limit.Draw("trackedParam_"+region+"_p4_4>>hi(100,-100,100)","")
				#rt.gDirectory.Get("hi").SetLineColor(rt.kBlack)
				#rt.gDirectory.Get("hi").SetLineWidth(2)
				#limit.Draw("trackedParam_"+region+"_p4_4>>hj(100,-100,100)",quantCut,"same")
				#rt.gDirectory.Get("hj").SetLineColor(rt.kRed)
				#rt.gDirectory.Get("hj").SetLineWidth(2)
				#rt.gDirectory.Get("hi").SetAxisRange(0,rt.gDirectory.Get("hi").GetMaximum()*1.1,"Y")
				#rt.gPad.Update()
				c1.cd(5)
				#limit.Draw("limitErr>>h8","")
				#rt.gDirectory.Get("h8").SetLineColor(rt.kBlack)
				#rt.gDirectory.Get("h8").SetLineWidth(2)
				limit.Draw("limitErr>>h9",quantCut)
				rt.gDirectory.Get("h9").SetLineColor(rt.kRed)
				rt.gDirectory.Get("h9").SetLineWidth(2)
				rt.gDirectory.Get("h9").SetAxisRange(0,rt.gDirectory.Get("h9").GetMaximum()*1.1,"Y")
				rt.gPad.Update()
				c1.cd(6)
				#limit.Draw("(limit-"+n+")/limitErr>>h10","")
				#rt.gDirectory.Get("h10").SetLineColor(rt.kBlack)
				#rt.gDirectory.Get("h10").SetLineWidth(2)
				limit.Draw("(limit-"+n+")/limitErr>>h11",quantCut)
				rt.gDirectory.Get("h11").SetLineColor(rt.kRed)
				rt.gDirectory.Get("h11").SetLineWidth(2)
				rt.gDirectory.Get("h11").SetAxisRange(0,rt.gDirectory.Get("h11").GetMaximum()*1.1,"Y")
				rt.gPad.Update()
				c1.cd(0)
				c1.SaveAs("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/"+region+"_fitParams_"+fullName+".png")
				c2 = rt.TCanvas("c2","c2",1000,500)
				c2.Divide(2,1)
				c2.cd(1)
				rt.gPad.SetLogy()
				listOfFuncs = []
				chi2hist = rt.TH1F("hist2Hist","Chi2",80,0,8)
				MaxChi2Dict = {} # index is key, chi2 is value
				MinChi2Dict = {} # index is key, chi2 is value
				minFuncVal = 1000
				maxFuncVal = 0.001
				for iEvt in range(limit.GetEntries()):
					limit.GetEvent(iEvt)
					#print(limit.iToy)
					if limit.quantileExpected != -1:
						continue
					tree.GetEvent(limit.iToy)
					rVal = tree.r
					toy = genOnlyFile.Get("toys/toy_{}".format(limit.iToy))
					toyExp = getattr(tree,"n_exp_final_bin"+region+"_2018_proc_roomultipdf")
					norm = 1
					#if region != "highSVJ2":
					p1 = getattr(limit,"trackedParam_"+region+"_p1_3")
					p2 = getattr(limit,"trackedParam_"+region+"_p2_3")
					p3 = getattr(limit,"trackedParam_"+region+"_p3_3")
					#p4 = getattr(limit,"trackedParam_"+region+"_p4_4")
					#copyToExcel += "{} {} {} {} {}\n".format(region, toyFunc, sigStr, p1, p2)
					#listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[2]*pow(1 - x/13000, [0]+[1]*log(x/13000)) * pow(x/13000,-[3]-[4]*log(x/13000))",1500,3800))
					listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[2]*pow(1 - x/13000, [0]+[1]*log(x/13000)) * pow(x/13000,-[3])",1500,3800))
					if rVal < 0:
						listOfFuncs[-1].SetLineColor(rt.kRed)
					else:
						listOfFuncs[-1].SetLineColor(rt.kBlue)
					listOfFuncs[-1].SetParameter(0,p1)
					listOfFuncs[-1].SetParameter(1,p2)
					listOfFuncs[-1].SetParameter(2,norm)
					listOfFuncs[-1].SetParameter(3,p3)
					#listOfFuncs[-1].SetParameter(4,p4)
					# pdf's in RooFit do not have a normlization facotr (i.e., Sum(allspace) of pdf = 1)
					# so we need to scale it to our dataset.
					# first (above line), set norm to 1, and we know the number of events in the toy)
					# then scale the function by numEvents divided by intergral of function without a normialztion
					# i.e., first normalize to 1 (divide by integral), then scale to numEvents
					# factor of 100 is becaue our bins are 100 GeV wide
					listOfFuncs[-1].SetParameter(2,toyExp/listOfFuncs[-1].Integral(1500,3800)*100)
					funcMax = listOfFuncs[-1].GetMaximum(1500,3800)
					funcMin = listOfFuncs[-1].GetMinimum(1500,3800)
					if funcMax > maxFuncVal:
						maxFuncVal = funcMax
					if funcMin < minFuncVal:
						minFuncVal = funcMin
					chi2 = 0
					for iBin in range(toy.numEntries()):
						toy.get(iBin)
						x = toy.get(iBin).getRealValue("mH")
						Oi = toy.weight()
						Ei = listOfFuncs[-1].Eval(x)
						try:
							chi2 += ((Oi-Ei)**2)/Ei
						except ZeroDivisionError:
							chi2 += 0
					chi2hist.Fill(chi2/20) # ndf = nBins - nParams = 23 - 2 = 21
					#else:
					#	p1 = getattr(limit,"trackedParam_"+region+"_p1_1")
					#	listOfFuncs.append(rt.TF1("iToy_"+str(iEvt-1),"[1]*pow(1 - x/13000, [0])",1500,3800))
					#	listOfFuncs[-1].SetParameter(0,p1)
					#	listOfFuncs[-1].SetParameter(1,norm)
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
					#	chi2hist.Fill(chi2/22)# ndf = nBins - nParams = 23 - 1 = 22
					MaxChi2Dict[limit.iToy] = chi2
					MinChi2Dict[limit.iToy] = chi2
				while len(MaxChi2Dict) > 20:
					minVal = min(list(MaxChi2Dict.values()))
					for index, Chi2val in MaxChi2Dict.items():
						if Chi2val == minVal:
							del MaxChi2Dict[index]
				while len(MinChi2Dict) > 20:
					maxVal = max(list(MinChi2Dict.values()))
					for index, Chi2val in MinChi2Dict.items():
						if Chi2val == maxVal:
							del MinChi2Dict[index]
				#print("++++++++++++++++++++++++++++++++++++++++++++++")
				#print(MaxChi2Dict.keys())
				#print(MinChi2Dict.keys())
				listOfFuncs[0].GetYaxis().SetRangeUser(0.1*minFuncVal,10*maxFuncVal)
				listOfFuncs[0].Draw()
				for thing in listOfFuncs[1:]:
					thing.Draw("same")
				dataHist.Draw("same")
				c2.Update()
				c2.cd(2)
				chi2hist.Draw()
				c2.SaveAs("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/"+region+"_manyFits_"+fullName+".png")
				c3 = rt.TCanvas("c3","c3",1000,1000)
				c3.SetLogy()
				opts = ""
				chi2MinText = rt.TPaveText(0.5,0.6,0.75,0.9)
				chi2MinText.AddText("Min Index and Vals")
				chi2MaxText = rt.TPaveText(0.75,0.6,1.0,0.9)
				chi2MaxText.AddText("Max Index and Vals")
				for iToy, chi2Val in MinChi2Dict.items() + MaxChi2Dict.items():
					if iToy in MinChi2Dict.keys():
						minMax = "Min"
						#print(len(listOfFuncs), iToy-1)
						listOfFuncs[iToy-1].SetLineColor(rt.kBlue)
						chi2MinText.AddText("{} {}".format(iToy, chi2Val))
					else:
						minMax = "Max"
						#print(len(listOfFuncs), iToy-1)
						listOfFuncs[iToy-1].SetLineColor(rt.kRed)
						chi2MaxText.AddText("{} {}".format(iToy, chi2Val))
					listOfFuncs[iToy-1].GetYaxis().SetRangeUser(0.01,100000)
					listOfFuncs[iToy-1].Draw(opts)
					opts = "same"
				chi2MaxText.Draw("same")
				chi2MinText.Draw("same")
				c3.SaveAs("cards_Mar3_ThryDij/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/minMax"+region+"/minMax_"+fullName+".png")
				fitDiagFile.Close()
				genOnlyFile.Close()
				dataObsFile.Close()
#print(copyToExcel)
























