import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function

regions = ["lowCut","lowSVJ1","lowSVJ2","highCut","highSVJ1","highSVJ2"]
#regions = ["lowSVJ1","lowSVJ2","highSVJ1","highSVJ2"]
#regions = ["lowSVJ2"]

# Generate 4 plots:
#	distribution of param1, for all and for qunatExp = 0.5
#	distribution of param2, for all and for quantExp = 0.5

copyToExcel = "Region & Toy & Sig & Fit & HistMean & HistStdDev & FitMean & FitStdDev\\\\\n"

for region in regions:
	for toyFunc in ["Main","Alt"]:
		for sigStr in ["Sig0"]:
			for fitFunc in ["Main"]:
				halfName = toyFunc+sigStr
				fullName = toyFunc+sigStr+fitFunc
				n = sigStr[3]
				fitDiagFile = rt.TFile.Open("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/fitDiagnostics"+fullName+".root","read")
				fitOnlyFile = rt.TFile.Open("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+fullName+".FitDiagnostics.mH120.123456.root","read")
				genOnlyFile = rt.TFile.Open("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+halfName+".GenerateOnly.mH125.123456.root","read")
				dataObsFile = rt.TFile.Open("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","read")
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
				listOfFuncs = []
				chi2hist = rt.TH1F("hist2Hist","Chi2",80,0,8)
				resHist = rt.TH2F("resHist","Residual;Bin;Toy",23,1,24,500,1,501)
				resProfile = rt.TProfile("resHist","Profile of Residual vs Bin;binNum (1500 < MT < 3800);Average(500 toys)",23,1,24)
				resProfileNorm = rt.TProfile("resHistNorm","Profile of Residual/funcValue vs Bin;binNum (1500 < MT < 3800);Average(500 toys)",23,1,24)
				resProfileNorm.SetMinimum(-0.4)
				resProfileNorm.SetMaximum(0.8)
				listOfChi2 = [0]
				for iEvt in range(limit.GetEntries()):
					limit.GetEvent(iEvt)
					if limit.quantileExpected != -1:
						continue
					tree.GetEvent(limit.iToy)
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
					#print(listOfFuncs[-1].Integral(1500,3800))
					chi2 = 0
					for iBin in range(toy.numEntries()):
						toy.get(iBin)
						x = toy.get(iBin).getRealValue("mH")
						Oi = toy.weight()
						Ei = listOfFuncs[-1].Eval(x)
						chi2 += ((Oi-Ei)**2)/Ei
						resHist.Fill(iBin,limit.iToy,Oi-Ei)
						resProfile.Fill(iBin, Oi-Ei)
						resProfileNorm.Fill(iBin, (Oi-Ei)/Ei)
					chi2hist.Fill(chi2/20) # ndf = nBins - nParams = 23 - 2 = 21
					listOfChi2.append(chi2/20)
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
					#		resHist.Fill(iBin,limit.iToy,Oi-Ei)
					#		resProfile.Fill(iBin, Oi-Ei)
					#		resProfileNorm.Fill(iBin, (Oi-Ei)/Ei)
					#		print(Oi-Ei)
					#	chi2hist.Fill(chi2/22)# ndf = nBins - nParams = 23 - 1 = 22
					#	listOfChi2.append(chi2/22)
				#for iToy in range(1,20):
				#	c3 = rt.TCanvas("c3","c3",1000,1000)
				#	frame = rt.RooPlot(svjWs.var("mH"),1500,3800,23)
				#	#frame.GetYaxis().SetRangeUser(0.0000001,10000000000)
				#	c3.SetLogy()
				#	toy = genOnlyFile.Get("toys/toy_{}".format(iToy))
				#	toy.plotOn(frame,rt.RooFit.LineColor(rt.kRed))
				#	if toyFunc == "Main":
				#		genPdf = svjWs.pdf("Bkg_"+region+"_2018")
				#	elif toyFunc == "Alt":
				#		genPdf = svjWs.pdf("Bkg_Alt_"+region+"_2018")
				#	genPdf.plotOn(frame,rt.RooFit.LineColor(rt.kBlue))
				#	frame.Draw()
				#	listOfFuncs[iToy].Draw("same")
				#	text = rt.TPaveText(0.5,0.6,0.9,0.9,"NDC ARC NB")
				#	genPdf.getParameters(toy).Print()
				#	genPdf.Print()
				#	if toyFunc == "Alt":
				#		ext = "_alt"
				#	else:
				#		ext = ""
				#	text.AddText("p1, p2 of GenFunc: {} {}".format(genPdf.getParameters(toy).getRealValue(region+"_p1_2"+ext),genPdf.getParameters(toy).getRealValue(region+"_p2_2"+ext)))
				#	text.AddText("p1, p2 of FitFunc: {} {}".format(listOfFuncs[iToy].GetParameter(0),listOfFuncs[iToy].GetParameter(1)))
				#	text.AddText("chi2 of fitFunc to Toy: {}".format(listOfChi2[iToy]))
				#	text.Draw("same")
				#	c3.SaveAs("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/"+region+"_toy{}_".format(iToy)+fullName+".png")
				c4 = rt.TCanvas("c4","c4",2000,1000)
				c4.Divide(2,1)
				rt.gStyle.SetOptStat(0)
				#c4.SetLogz()
				c4.cd(1)
				resProfile.SetLineWidth(3)
				resProfile.Draw("e1")
				c4.cd(2)
				resProfileNorm.SetLineWidth(3)
				resProfileNorm.Draw("e1")
				c4.SaveAs("cards_Feb12_3params/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/"+region+"_resHist_"+fullName+".png")
				fitDiagFile.Close()
				genOnlyFile.Close()
				dataObsFile.Close()
print(copyToExcel)
























