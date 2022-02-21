import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# runover each SR to fit all 100 toys to the main function

#regions = ["lowSVJ1"]
regions = ["highCut"]

# Generate plots of fit_status, r, rErr, and (r-STR)/rErr for all, fit_status <300, and fit_status==0
iToy = 1
fileError 			= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_Error/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")
fileErrorDoFit		= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_ErrorDoFit/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")
fileNoError 		= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_NoError/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")
fileNoErrorDoFit	= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_NoErrorDoFit/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")
fileBigError 		= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_BigError/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")
fileBigErrorDoFit	= rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/highCut_BigErrorDoFit/higgsCombineMainSig0.GenerateOnly.mH125.123456.root","read")

dataObsFile = rt.TFile.Open("cards_Feb21_Thry_23bins/SVJ_mZprime3000_mDark20_rinv03_alphapeak/ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_highCut_2018_template.root","read")
svjWs = dataObsFile.Get("SVJ")
c1 = rt.TCanvas("c1","c1",1000,1000)
frameA = rt.RooPlot(svjWs.var("mH"),1500.,3800.,23)
frameB = rt.RooPlot(svjWs.var("mH"),1500.,3800.,23)

data_obs = svjWs.data("data_obs")

toyError 		= fileError.Get("toys/toy_{}".format(iToy))
toyErrorDoFit 		= fileErrorDoFit.Get("toys/toy_{}".format(iToy))
toyNoError 		= fileNoError.Get("toys/toy_{}".format(iToy))
toyNoErrorDoFit 		= fileNoErrorDoFit.Get("toys/toy_{}".format(iToy))
toyBigError 		= fileBigError.Get("toys/toy_{}".format(iToy))
toyBigErrorDoFit 		= fileBigErrorDoFit.Get("toys/toy_{}".format(iToy))

toyBigError.plotOn(frameA,rt.RooFit.MarkerColor(rt.kBlue), rt.RooFit.MarkerSize(3))
toyBigErrorDoFit.plotOn(frameB,rt.RooFit.MarkerColor(rt.kBlue+2), rt.RooFit.MarkerSize(3))
toyError.plotOn(frameA,rt.RooFit.MarkerColor(rt.kGreen), rt.RooFit.MarkerSize(2))
toyErrorDoFit.plotOn(frameB,rt.RooFit.MarkerColor(rt.kGreen+2), rt.RooFit.MarkerSize(2))
toyNoError.plotOn(frameA,rt.RooFit.MarkerColor(rt.kRed), rt.RooFit.MarkerSize(1))
toyNoErrorDoFit.plotOn(frameB,rt.RooFit.MarkerColor(rt.kRed+2), rt.RooFit.MarkerSize(1))

print(type(toyBigError), toyBigError.sumEntries("mH<=1600"))
print(type(data_obs), data_obs.sumEntries("mH<=1600"))

SF = toyBigError.sumEntries("mH<=1600")/data_obs.sumEntries("mH<=1600")
print(SF)
data_obs.plotOn(frameA,rt.RooFit.Rescale(SF))
data_obs.plotOn(frameB)


rt.gPad.SetLogy()
frameA.Draw()
c1.SaveAs("ErrorVsNoError_noFit_highCut.png")
frameB.Draw()
c1.SaveAs("ErrorVsNoError_doFit_highCut.png")

































