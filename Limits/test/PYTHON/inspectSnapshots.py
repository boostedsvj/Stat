import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function

regions = ["lowSVJ1","lowSVJ2","highSVJ1","highSVJ2"]
#regions = ["highSVJ2"]

# Generate 4 plots:
#	distribution of param1, for all and for qunatExp = 0.5
#	distribution of param2, for all and for quantExp = 0.5

copyToExcel = "Region & Toy & Sig & Fit & HistMean & HistStdDev & FitMean & FitStdDev\\\\\n"

for region in regions:
	for toyFunc in ["Main","Alt"]:
		for sigStr in ["Sig1", "Sig0"]:
			for fitFunc in ["Main"]:
				halfName = toyFunc+sigStr
				fullName = toyFunc+sigStr+fitFunc
				n = sigStr[3]
				fitDiagFile = rt.TFile.Open("cards_21Jan/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/fitDiagnostics"+fullName+".root","read")
				fitOnlyFile = rt.TFile.Open("cards_21Jan/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+fullName+".FitDiagnostics.mH120.123456.root","read")
				genOnlyFile = rt.TFile.Open("cards_21Jan/SVJ_mZprime3000_mDark20_rinv03_alphapeak/"+region+"/higgsCombine"+halfName+".GenerateOnly.mH125.123456.root","read")
				dataObsFile = rt.TFile.Open("cards_21Jan/SVJ_mZprime3000_mDark20_rinv03_alphapeak/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template.root","read")
				svjWs = dataObsFile.Get("SVJ")
				toysDirect = genOnlyFile.Get("toys")
				#print(type(svjWs),type(toysDirect))
				#toy = rt.RooDataHist()
				#data = rt.RooDataHist()
				data = svjWs.data("data_obs")
				print("data",type(data))
				dataHist = data.createHistogram("mH")
				fitPdf = svjWs.pdf("Bkg_"+region+"_2018")
				if float(fitDiagFile.GetSize()) < 20000:
					print("="*25)
					print("{} {} could not be drawn".format(region, fullName))
					print("="*25)
					continue
				tree = fitDiagFile.Get("tree_fit_sb")
				limit = fitOnlyFile.Get("limit")
				for iToy in range(1,20):
					#c3 = rt.TCanvas("c3","c3",1000,1000)
					#frame = rt.RooPlot(svjWs.var("mH"),1500,3800,23)
					#c3.SetLogy()
					toy = genOnlyFile.Get("toys/toy_{}".format(iToy))
					toySH = genOnlyFile.Get("toys/toy_{}_snapshot".format(iToy))
					toy.Print()
					toySH.Print()
					print(toySH.getRealValue("lumi_In"))
					#toy.plotOn(frame,rt.RooFit.LineColor(rt.kRed))
					#if toyFunc == "Main":
					#	genPdf = svjWs.pdf("Bkg_"+region+"_2018")
					#elif toyFunc == "Alt":
					#	genPdf = svjWs.pdf("Bkg_Alt_"+region+"_2018")
					#genPdf.plotOn(frame,rt.RooFit.LineColor(rt.kBlue))
					#frame.Draw()
					#listOfFuncs[iToy].Draw("same")
					#c3.SaveAs("cards_21Jan/SVJ_mZprime3000_mDark20_rinv03_alphapeak/plots/"+region+"_toy{}_".format(iToy)+fullName+".png")
				fitDiagFile.Close()
				genOnlyFile.Close()
				dataObsFile.Close()
print(copyToExcel)
























