import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function
eosArea = "root://cmsxrootd.fnal.gov//store/user/cfallon/biasStudies/"

listOfParams = [
['3000', '20', '03', 'peak'],
['1500', '20', '03', 'peak'],
['1600', '20', '03', 'peak'],
['1700', '20', '03', 'peak'],
['1800', '20', '03', 'peak'],
['1900', '20', '03', 'peak'],
['2000', '20', '03', 'peak'],
['2100', '20', '03', 'peak'],
['2200', '20', '03', 'peak'],
['2300', '20', '03', 'peak'],
['2400', '20', '03', 'peak'],
['2500', '20', '03', 'peak'],
['2600', '20', '03', 'peak'],
['2700', '20', '03', 'peak'],
['2800', '20', '03', 'peak'],
['2900', '20', '03', 'peak'],
['3100', '20', '03', 'peak'],
['3200', '20', '03', 'peak'],
['3300', '20', '03', 'peak'],
['3400', '20', '03', 'peak'],
['3500', '20', '03', 'peak'],
['3600', '20', '03', 'peak'],
['3700', '20', '03', 'peak'],
['3800', '20', '03', 'peak'],
['3900', '20', '03', 'peak'],
['4000', '20', '03', 'peak'],
['4100', '20', '03', 'peak'],
['4200', '20', '03', 'peak'],
['4300', '20', '03', 'peak'],
['4400', '20', '03', 'peak'],
['4500', '20', '03', 'peak'],
['3000', '1', '03', 'peak'],
['3000', '5', '03', 'peak'],
['3000', '10', '03', 'peak'],
['3000', '30', '03', 'peak'],
['3000', '40', '03', 'peak'],
['3000', '50', '03', 'peak'],
['3000', '60', '03', 'peak'],
['3000', '70', '03', 'peak'],
['3000', '80', '03', 'peak'],
['3000', '90', '03', 'peak'],
['3000', '100', '03', 'peak'],
['3000', '20', '0', 'peak'],
['3000', '20', '01', 'peak'],
['3000', '20', '02', 'peak'],
['3000', '20', '04', 'peak'],
['3000', '20', '05', 'peak'],
['3000', '20', '06', 'peak'],
['3000', '20', '07', 'peak'],
['3000', '20', '08', 'peak'],
['3000', '20', '09', 'peak'],
['3000', '20', '1', 'peak'],
['3000', '20', '03', 'low'],
['3000', '20', '03', 'high']
]
#regions = ["lowCut"]#,"lowSVJ1","lowSVJ2","highCut","highSVJ1","highSVJ2"]
regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
#regions = ["lowSVJ2"]
expSig = "_expSig1_extra"# "" for excpSig = 0, "_expSig1" for expSig = 1
# Generate plots of fit_status, r, rErr, and (r-STR)/rErr for all, fit_status <300, and fit_status==0
pae = "mZ\tmD\trI\taD\tregion\tpullMean\tpullStdDev\n"

outFile = open("../condorTests/gausFits"+expSig+".txt","w+")
outFile.write(pae)


#copyToExcel = "Region & Toy & Sig & Fit & HistMean & HistStdDev & FitMean & FitStdDev\\\\\n"
for sigPars in listOfParams:
	SVJNAME = "SVJ_mZprime{}_mDark{}_rinv{}_alpha{}".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3])
	print("************************",SVJNAME,"************************")
	for region in regions:
		print("************************",region, "************************")
		#if ()
		fitDiagFile = rt.TFile.Open(eosArea + SVJNAME + "/fitDiagnostics"+region+expSig+".root","read")
		if type(fitDiagFile) == type(rt.TFile()):
			outFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3],region, "DNE","DNE"))
			continue
		dataObsFile = rt.TFile.Open(eosArea + SVJNAME + "/ws_"+SVJNAME+"_"+region+"_2018_template.root","read")
		svjWs = dataObsFile.Get("SVJ")
		data = rt.RooDataHist()
		data = svjWs.data("data_obs")
		fitPdf = svjWs.pdf("Bkg_"+region+"_2018")
		print(type(fitDiagFile))
		if float(fitDiagFile.GetSize()) < 20000:
			print("="*25)
			print("{} could not be drawn".format(region))
			print("="*25)
			continue
		tree = fitDiagFile.Get("tree_fit_sb")
		c1 = rt.TCanvas("c1","c1",1000,1000)
		#c1.Divide(2,2)
		#c1.cd(1)
		selection = "" #"fit_status==0||fit_status==1"
		#tree.Draw("r>>h1(60,-20,20)",selection)
		#rt.gDirectory.Get("h1").SetLineColor(rt.kBlack)
		#rt.gDirectory.Get("h1").SetLineWidth(2)
		#tree.Draw("r>>h1p","r>0","same")
		#rt.gDirectory.Get("h1p").SetLineColor(rt.kRed)
		#rt.gDirectory.Get("h1").SetAxisRange(0,rt.gDirectory.Get("h1").GetMaximum()*1.1,"Y")
		#rt.gPad.Update()
		#c1.cd(2)
		#tree.Draw("rErr>>h2(100,0,10)",selection)
		#rt.gDirectory.Get("h2").SetLineColor(rt.kBlack)
		#rt.gDirectory.Get("h2").SetLineWidth(2)
		#tree.Draw("rErr>>h2p","r>0","same")
		#rt.gDirectory.Get("h2p").SetLineColor(rt.kRed)
		#rt.gDirectory.Get("h2").SetAxisRange(0,rt.gDirectory.Get("h2").GetMaximum()*1.1,"Y")
		#rt.gPad.Update()
		#c1.cd(3)
		#frame = rt.RooPlot(svjWs.var("mH"),1500.,3800.,23)
		#rame = rt.RooPlot(svjWs.var("mH"),1500.,6000.,45)
		#frame.SetTitle("data, inital fit")
		#rt.gPad.SetLogy()
		#data.plotOn(frame, rt.RooFit.MarkerColor(rt.kBlack))
		#fitPdf.plotOn(frame,rt.RooFit.LineColor(rt.kMagenta))
		#frame.SetAxisRange(0.001,50000,"Y")
		#frame.Draw()
		#c1.cd(4)
		injSig = 0 if expSig == "" else 1
		tree.Draw("(r-{})/rErr>>h3(100,-20,20)".format(injSig),selection)
		rt.gDirectory.Get("h3").SetLineColor(rt.kBlack)
		rt.gDirectory.Get("h3").SetLineWidth(2)
		#tree.Draw("(r-"+n+")/rErr>>h3p","r>0","same")
		#rt.gDirectory.Get("h3p").SetLineColor(rt.kRed)
		rt.gDirectory.Get("h3").SetAxisRange(0,rt.gDirectory.Get("h3").GetMaximum()*1.1,"Y")
		rt.gPad.Update()
		gaus = rt.TF1("gaus","gaus(0)", rt.gDirectory.Get("h3").GetMean()-2*rt.gDirectory.Get("h3").GetStdDev(), rt.gDirectory.Get("h3").GetMean()+2*rt.gDirectory.Get("h3").GetStdDev())
		rt.gDirectory.Get("h3").Fit("gaus","R")
		newLine = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3],region, gaus.GetParameter(1), gaus.GetParameter(2)) 
		outFile.write(newLine)
		pae += newLine
		gaus.SetLineColor(rt.kBlue)
		gaus.Draw("same")
		#copyToExcel += "{} & {} & {} & {} & {} & {} & {} & {} \\\\\n".format(region, toyFunc,sigStr,fitFunc,rt.gDirectory.Get("h3").GetMean(),rt.gDirectory.Get("h3").GetStdDev(),gaus.GetParameter(1),gaus.GetParameter(2))
		#c1.cd(0)
		#title = rt.TPaveText(0.3,0.5,0.7,0.53)
		title = rt.TPaveText(0.1,0.9,0.4,0.93)
		title.AddText(SVJNAME + "{}".format(region))
		title.Draw("same")
		c1.SaveAs("~/nobackup/SVJ/biasStudies2/CMSSW_10_2_13/src/Stat/Limits/test/condorTests/"+SVJNAME+"_"+region+"_fit_status_plots"+expSig+".png")

		fitDiagFile.Close()
		dataObsFile.Close()
print(pae)
outFile.close()
























