import ROOT as rt
import os
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)

# plot:
# snapshot with Toys
# snapshot with Fits
# toys with Fits
eosArea = "root://cmsxrootd.fnal.gov//store/user/cfallon/biasStudies_ThreeSteps/"


listOfParams1 = [
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
['2900', '20', '03', 'peak']]

listOfParams2 = [
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
['3000', '1', '03', 'peak']]

listOfParams3 = [
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
['3000', '20', '06', 'peak']]

listOfParams4 = [
['3000', '20', '07', 'peak'],
['3000', '20', '08', 'peak'],
['3000', '20', '09', 'peak'],
['3000', '20', '1', 'peak'],
['3000', '20', '03', 'low'],
['3000', '20', '03', 'high']]

listOfParams5 = [
['2200', '20', '03', 'peak'],
['2400', '20', '03', 'peak']]

baseline = [['3000', '20', '03', 'peak']]

#regions = ["lowSVJ2","highSVJ2"]
regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
#expSig = ""# "" for excpSig = 0, "_expSig1" for expSig = 1, add "_extra" for SVJ options, or nothing for Dijet options
#n = 0 if expSig == "" else int(expSig.split("_")[1][-1:])
for sigPars in baseline:
	SVJNAME = "SVJ_mZprime{}_mDark{}_rinv{}_alpha{}".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3])
	print("************************",SVJNAME,"************************")
	for region in regions:
		for expSig in ["Sig0GenMainFitMain","Sig1GenMainFitMain","Sig0GenAltFitMain","Sig1GenAltFitMain"]:
			for combineOpts in ["OptS","OptD"]:
				nameNoFit = region+combineOpts+expSig[:-7]
				bigName = region+combineOpts+expSig
				n = int(expSig[3])
				# files: 
				# bigName : fitDiag, hC.FD
				# nameNoFit : hC.GO, hC.MDF
				fitDiagFile = rt.TFile.Open(eosArea + SVJNAME + "/fitDiagnostics"+bigName+".root","read")
				fitOnlyFile = rt.TFile.Open(eosArea + SVJNAME + "/higgsCombine"+bigName+".FitDiagnostics.mH120.123456.root","read")
				genOnlyFile = rt.TFile.Open(eosArea + SVJNAME + "/higgsCombine"+nameNoFit+".GenerateOnly.mH120.123456.root","read")
				multiDFFile = rt.TFile.Open(eosArea + SVJNAME + "/higgsCombine"+nameNoFit+".MultiDimFit.mH120.root","read")
				#dataObsFile = rt.TFile.Open(eosArea + SVJNAME + "/ws_"+SVJNAME+"_"+region+"_2018_template.root","read")
				#if (type(fitDiagFile) == type(rt.TFile())) or (type(fitOnlyFile) == type(rt.TFile())) or(type(dataObsFile) == type(rt.TFile())):
				if (type(fitDiagFile) == type(rt.TFile())) or (type(fitOnlyFile) == type(rt.TFile()) or (type(genOnlyFile) == type(rt.TFile()) or (type(multiDFFile) == type(rt.TFile())):
					continue
				#svjWs = dataObsFile.Get("SVJ")
				mdfWS = multiDFFile.Get("w")
				mdfWS.loadSnapshot("MultiDimFit") # load the snapshot parameters from the fit to bkg MC only, non-Gen parameters are 0
				
				#data = rt.RooDataHist()
				#data = svjWs.data("data_obs")
				#genPdf = svjWs.pdf("Bkg_"+region+"_2018")
				tree = fitDiagFile.Get("tree_fit_sb")
				limit = fitOnlyFile.Get("limit")
				
				if "low" in region:
					nPar = 2
					nParGen = 2
				elif "highCut" == region:
					nPar = 3
					nParGen = (2 if ("GenAlt" in expSig) else 3 )
				elif "highSVJ2" == region:
					nPar = 1
					nParGen = 1
				if "GenAlt" in expSig:
					genExtra = "_alt"
				else:
					genExtra = ""

				if nGenPar == 1:
					p1Gen = getattr(limit,"trackedParam_"+region+"_p1_1"+genExtra)
					if genExtra = "_alt"
						genFunc = rt.TF1("GenFuncAlt","[0] * exp(@1*(@0/13000))",1500,6000))
					else:
						genFunc = rt.TF1("GenFuncMain","[0] * pow(x/13000, -[1])",1500,6000))
				elif nGenPar == 2:
					p1Gen = getattr(limit,"trackedParam_"+region+"_p1_2"+genExtra)
					p2Gen = getattr(limit,"trackedParam_"+region+"_p2_2"+genExtra)
					if genExtra = "_alt"
						genFunc = rt.TF1("GenFuncAlt","[0] * exp(@1*(@0/13000)+@2*log(@0/13000))",1500,6000))
					else:
						genFunc = rt.TF1("GenFuncMain","[0] * pow(1 - x/13000, [2]) * pow(x/13000, -[1])",1500,6000))
				elif nGenPar == 3:
					p1Gen = getattr(limit,"trackedParam_"+region+"_p1_3"+genExtra)
					p2Gen = getattr(limit,"trackedParam_"+region+"_p2_3"+genExtra)
					p3Gen = getattr(limit,"trackedParam_"+region+"_p3_3"+genExtra)
					if genExtra = "_alt"
						genFunc = rt.TF1("GenFuncAlt","[0] * exp(@1*(@0/13000)+@2*log(@0/13000)+@3*pow(log(@0/13000),2))",1500,6000))
					else:
						genFunc = rt.TF1("GenFuncMain","[0] * pow(1 - x/13000, [2]) * pow(x/13000, -[1]-[3]*log(x/13000))",1500,6000))
				else:
					print("nPar is not 1, 2, nor 3!")
					exit(0)
				genFunc.SetParameter(0, 
				genFunc.SetParameter(1,p1Gen)
				if nParGen >= 2:
					genFunc[-1].SetParameter(2,p2Gen)
				if nParGen >= 3:
					genFunc[-1].SetParameter(3,p3Gen)
				
				listOfFuncs = []
				for iEvt in range(limit.GetEntries()):
					limit.GetEvent(iEvt)
					if limit.quantileExpected != -1:
						continue
					tree.GetEvent(limit.iToy)
					toy = fitOnlyFile.Get("toys/toy_{}".format(limit.iToy))
					#toyExp = getattr(tree,"n_exp_final_bin"+region+"_2018_proc_roomultipdf")
					toyExp = toy.sumEntries()


					rHistVal = tree.r
					rErrHistVal = tree.rErr
					if rErrHistVal == 0:
						continue
					rmuHistVal = (tree.r-int(n))/tree.rErr
					if nPar == 1:
						p1 = getattr(limit,"trackedParam_"+region+"_p1_1")
						listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(x/13000, -[1])",1500,6000))
					elif nPar == 2:
						p1 = getattr(limit,"trackedParam_"+region+"_p1_2")
						p2 = getattr(limit,"trackedParam_"+region+"_p2_2")
						listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [2]) * pow(x/13000, -[1])",1500,6000))
					elif nPar == 3:
						p1 = getattr(limit,"trackedParam_"+region+"_p1_3")
						p2 = getattr(limit,"trackedParam_"+region+"_p2_3")
						p3 = getattr(limit,"trackedParam_"+region+"_p3_3")
						listOfFuncs.append(rt.TF1("iToy_"+str(iEvt),"[0] * pow(1 - x/13000, [2]) * pow(x/13000, -[1]-[3]*log(x/13000))",1500,6000))
					else:
						print("nPar is not 1, 2, nor 3!")
						exit(0)

					listOfFuncs[-1].SetParameter(1,p1)
					if nPar >= 2:
						listOfFuncs[-1].SetParameter(2,p2)
					if nPar >= 3:
						listOfFuncs[-1].SetParameter(3,p3)
					# pdf's in RooFit do not have a normlization factor (i.e., Sum(allspace) of pdf = 1)
					# so we need to scale it to our dataset.
					# first, set norm to 1, and we know the number of events in the toy)
					# then scale the function by numEvents divided by intergral of function without a normialztion
					# i.e., first normalize to 1 (divide by integral), then scale to numEvents
					# factor of 50 is becaue our bins are 50 GeV wide
					norm = 1
					listOfFuncs[-1].SetParameter(0,norm)
					denom = listOfFuncs[-1].Integral(1500,6000) 
					if denom > 0:
						norm = toyExp/listOfFuncs[-1].Integral(1500,6000)*50
					else:
						continue
					listOfFuncs[-1].SetParameter(0,norm)


				c1 = rt.TCanvas("c1","c1",1500,500)
				c1.Divide(3,1)
				c1.cd(1) # snapshot func with toys
				c1.cd(2) # snapshot Func with fitFuncs
				c1.cd(3) # toys with fit funcs
				
				c1.SaveAs("../condor_ThreeSteps/"+SVJNAME+"_"+bigName+".png")

				fitDiagFile.Close()
				fitOnlyFile.Close()
				#dataObsFile.Close()
























