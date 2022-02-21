from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os,sys
import subprocess

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--eosDir", dest="eosDir", type=str, required=True, help="eos directory for input files (/store/...)")
parser.add_argument("-l", "--doLimit", dest="doLimit", default=False, action="store_true", help="limit mode (use combined regions)")
parser.add_argument("-m", "--maxVal", dest="maxVal", type=float, default=0, required=False, help="maximum allowed output (default = 0, no max value)")
parser.add_argument("-v", "--varyParam", dest="varyParam", type=str, default="mZ", required=False, help="Which param to vary? (mZ, mD, rI, aD)")
parser.add_argument("-Z", "--baselineZ", dest="baselineZ", type=str, default="3100", required=False, help="If mZ isn't varied, what value to use? (default 3100)")
parser.add_argument("-D", "--baselineD", dest="baselineD", type=str, default="20", required=False, help="If mD isn't varied, what value to use? (default 20)")
parser.add_argument("-R", "--baselineR", dest="baselineR", type=str, default="03", required=False, help="If rI isn't varied, what value to use? (default 03)")
parser.add_argument("-A", "--baselineA", dest="baselineA", type=str, default="peak", required=False, help="If aD isn't varied, what value to use? (default peak)")

args = parser.parse_args()

# runover each SR to fit all 100 toys to the main function
eosArea = "root://cmseos.fnal.gov/"+args.eosDir+"/"

import ROOT as rt
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(1)
rt.gStyle.SetOptFit(1011)
rt.gROOT.SetBatch(True)

#generate list of parameters to run over.
# only available for mZ and rI at the moment

if args.varyParam == "mZ":
	listOfParams = [
	["1500", args.baselineD, args.baselineR, args.baselineA],
	["1700", args.baselineD, args.baselineR, args.baselineA],
	["1900", args.baselineD, args.baselineR, args.baselineA],
	["2100", args.baselineD, args.baselineR, args.baselineA],
	["2300", args.baselineD, args.baselineR, args.baselineA],
	["2500", args.baselineD, args.baselineR, args.baselineA],
	["2700", args.baselineD, args.baselineR, args.baselineA],
	["2900", args.baselineD, args.baselineR, args.baselineA],
	["3100", args.baselineD, args.baselineR, args.baselineA],
	["3300", args.baselineD, args.baselineR, args.baselineA],
	["3500", args.baselineD, args.baselineR, args.baselineA],
	["3700", args.baselineD, args.baselineR, args.baselineA],
	["3900", args.baselineD, args.baselineR, args.baselineA],
	["4100", args.baselineD, args.baselineR, args.baselineA],
	["4300", args.baselineD, args.baselineR, args.baselineA],
	["4500", args.baselineD, args.baselineR, args.baselineA],
	["4700", args.baselineD, args.baselineR, args.baselineA],
	["4900", args.baselineD, args.baselineR, args.baselineA],
	["5100", args.baselineD, args.baselineR, args.baselineA]]
if args.varyParam == "rI":
	listOfParams = [
	[args.baselineZ, args.baselineD, "01", args.baselineA],
	[args.baselineZ, args.baselineD, "03", args.baselineA],
	[args.baselineZ, args.baselineD, "05", args.baselineA],
	[args.baselineZ, args.baselineD, "07", args.baselineA],
	[args.baselineZ, args.baselineD, "09", args.baselineA]]

choices = subprocess.check_output(["eos","root://cmseos.fnal.gov/","ls",args.eosDir])

if args.doLimit:
	regions = ["cut","bdt"]
	sigList = ["Sig3M"]
else:
	regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
	sigList = ["Sig0","Sig1"]
for expSig in sigList:#["Sig1","Sig1.5","Sig2","Sig3","Sig5","Sig7","Sig10"]:#,"SigM"]:#["SigM"]:#"Sig0","Sig1",
	#for funcs in ["GenMainFitMain","GenAltFitMain"]:
	for funcs in ["GenAltFitAlt", "GenMainFitAlt"]:
		#setup four TGraphErrors, one for each varying-variable
		if args.doLimit:
			vZ_m = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyZ_mean"}
			vZ_s = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyZ_stdev"}
			vD_m = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyD_mean"}
			vD_s = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyD_stdev"}
			vR_m = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyR_mean"}
			vR_s = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyR_stdev"}
			vA_m = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyA_mean"}
			vA_s = {"cut":[[],[],[]],"bdt":[[],[],[]], "name":"varyA_stdev"}
		else:
			vZ_m = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyZ_mean"}
			vZ_s = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyZ_stdev"}
			vZ_c = {"lowCut":[[],[]],"lowSVJ2":[[],[]],"highCut":[[],[]],"highSVJ2":[[],[]], "name":"varyZ_chi2"}
			vD_m = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyD_mean"}
			vD_s = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyD_stdev"}
			vR_m = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyR_mean"}
			vR_s = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyR_stdev"}
			vA_m = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyA_mean"}
			vA_s = {"lowCut":[[],[],[]],"lowSVJ2":[[],[],[]],"highCut":[[],[],[]],"highSVJ2":[[],[],[]], "name":"varyA_stdev"}


		for sigPars in listOfParams:
			SVJNAME = "SVJ_mZprime{}_mDark{}_rinv{}_alpha{}".format(sigPars[0],sigPars[1],sigPars[2],sigPars[3])
			if not (SVJNAME in choices): continue
			print("************************",SVJNAME,"************************")
			#(re)set varyBools to all be false
			varyZ = args.varyParam=="mZ"
			varyD = args.varyParam=="mD"
			varyR = args.varyParam=="rI"
			varyA = args.varyParam=="aD"
			zMass = int(sigPars[0])
			dMass = int(sigPars[1])
			rInvs = 1 if sigPars[2] == "1" else int(sigPars[2])*0.1
			if sigPars[3] == "peak":
				aDark = 1
			elif sigPars[3] == "low":
				aDark = 0.5
			elif sigPars[3] == "high":
				aDark = 1.5
			for region in regions:
				print("************************",region, "************************")
				fileName = eosArea + SVJNAME + "/fitDiagnostics"+region+expSig+funcs+".mH120.123456.root"
				# check if file exists, skipping if it doesnt
				fileDir = subprocess.check_output(["eos","root://cmseos.fnal.gov","ls",eosArea.split("//")[2]+SVJNAME])
				if not ("fitDiagnostics"+region+expSig+funcs+".mH120.123456.root" in fileDir): continue
				fitDiagFile = rt.TFile.Open(fileName,"read")
				if float(fitDiagFile.GetSize()) < 20000:
					fitDiagFile.Close()
					print("="*25)
					print("{} could not be drawn".format(fileName))
					print("="*25)
					continue
				tree = fitDiagFile.Get("tree_fit_sb")
				#print(tree.Print())
				c1 = rt.TCanvas("c1","c1",1500,500)
				c1.Divide(3,1)
				selection = "fit_status==0"
				if "M" in expSig:
					_file = rt.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Limits/datacards_07tsb_sys/limit_{}AltManualBFInitSyst.root".format(region),"read")
					limitTree = _file.Get("limit")
					for event in limitTree:
						if limitTree.trackedParam_mZprime == int(sigPars[0]) and limitTree.quantileExpected == 0.5:
							injSig = limitTree.limit
					_file.Close()
					if expSig=="Sig3M": injSig *= 3
					if args.maxVal>0 and injSig>args.maxVal: injSig = args.maxVal
				else:
					try:
						injSig = int(expSig[3:])
					except ValueError:
						injSig = float(expSig[3:])
				print("INJSIG CHECK ************ ", str(injSig))
				c1.cd(3)
				minPull,maxPull = -5,5
				nToys = tree.GetEntries()
				nBins = int((nToys/5.0))
				tree.Draw("(r-{})/rErr>>h3({},{},{})".format(injSig,nBins,minPull,maxPull),selection)
				#tree.Draw("(r-{})/rErr>>h4({},{},{})".format(injSig,nBins,minPull,maxPull),selection+"&&rErr<100","same")
				#tree.Draw("(r-{})/rErr>>h4(50,-5,5)".format(injSig),"fit_status==0","same")
				#rt.gDirectory.Get("h4").SetLineColor(rt.kRed)
				rt.gDirectory.Get("h3").SetLineColor(rt.kBlack)
				rt.gDirectory.Get("h3").SetLineWidth(2)
				rt.gDirectory.Get("h3").SetAxisRange(0,rt.gDirectory.Get("h3").GetMaximum()*1.1,"Y")
				#rt.gDirectory.Get("h4").SetLineColor(rt.kRed)
				#rt.gDirectory.Get("h4").SetLineWidth(2)
				#rt.gDirectory.Get("h4").SetAxisRange(0,rt.gDirectory.Get("h4").GetMaximum()*1.1,"Y")
				rt.gPad.Update()
				gaus = rt.TF1("gaus","gaus(0)",minPull,maxPull)
				fitResult = rt.gDirectory.Get("h3").Fit("gaus","RS")
				gaus.SetLineColor(rt.kBlue)
				gaus.Draw("same")
				c1.cd(1)
				tree.Draw("r>>h1({0},-{1},{1})".format(nBins, max(injSig*3,5000)),selection)
				#tree.Draw("r>>h5(160,-40,40)".format(injSig),"fit_status==0","same")
				#rt.gDirectory.Get("h5").SetLineColor(rt.kRed)
				c1.cd(2)
				tree.Draw("rErr>>h2({},0,{})".format(nBins, max(5000,injSig*5)),selection)
				c1.SaveAs("./plots/"+SVJNAME+"_"+region+expSig+funcs+"_3PlotsR.pdf")
				if gaus.GetNDF() == 0:
					continue
				if varyZ:
					vZ_m[region][0].append(zMass)
					vZ_m[region][1].append(gaus.GetParameter(1))
					vZ_m[region][2].append(gaus.GetParError(1))
					if not args.doLimit:
						vZ_c[region][0].append(zMass)
						vZ_c[region][1].append(fitResult.Chi2()/fitResult.Ndf())
					vZ_s[region][0].append(zMass)
					vZ_s[region][1].append(gaus.GetParameter(2))
					vZ_s[region][2].append(gaus.GetParError(2))
				if varyD:
					vD_m[region][0].append(dMass)
					vD_m[region][1].append(gaus.GetParameter(1))
					vD_m[region][2].append(gaus.GetParError(1))
					vD_s[region][0].append(dMass)
					vD_s[region][1].append(gaus.GetParameter(2))
					vD_s[region][2].append(gaus.GetParError(2))
				if varyR:
					vR_m[region][0].append(rInvs)
					vR_m[region][1].append(gaus.GetParameter(1))
					vR_m[region][2].append(gaus.GetParError(1))
					vR_s[region][0].append(rInvs)
					vR_s[region][1].append(gaus.GetParameter(2))
					vR_s[region][2].append(gaus.GetParError(2))
				if varyA:
					vR_m[region][0].append(aDark)
					vR_m[region][1].append(gaus.GetParameter(1))
					vR_m[region][2].append(gaus.GetParError(1))
					vR_s[region][0].append(aDark)
					vR_s[region][1].append(gaus.GetParameter(2))
					vR_s[region][2].append(gaus.GetParError(2))
					

				fitDiagFile.Close()
		for region in regions:
			if args.doLimit:
				BASELINE = ""
				vecs = [vZ_m, vZ_s]
			if args.varyParam == "mZ":
				BASELINE = args.baselineR+"_"
				vecs = [vZ_m, vZ_s]
			elif args.varyParam == "rI":
				BASELINE = args.baselineZ+"_"
				vecs = [vR_m, vR_s]

			for vec in vecs:
				out = open("./plots/"+vec["name"].replace("_","/")+"_"+BASELINE+region+expSig+funcs+".txt","w")
				out.write("#varVal mean/stdev error\n")
				for i in range(len(vec[region][0])):
					if len(vec[region]) == 3:
						out.write("{} {} {}\n".format(vec[region][0][i], vec[region][1][i], vec[region][2][i]))
					if len(vec[region]) == 2:
						out.write("{} {}\n".format(vec[region][0][i], vec[region][1][i]))
				out.close()










