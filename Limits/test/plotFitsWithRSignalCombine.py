import sys
region = sys.argv[-2]
mass = sys.argv[-1]
if len(sys.argv)>3:
    verbose = len(sys.argv[-3])>0
else:
    verbose = False

import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)
rt.gStyle.SetStatH(0.30)
rt.gStyle.SetStatW(0.20)
import math
import numpy as np
from array import array
from scipy.special import gammaln

# plot dataObs + r*Sig for atleast r = 0 and r = min(deltaNLL)

#regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
regions = [region]
#masses = ["2000","3000","4000"]
masses = [mass]

# disable root garbage collection
def ngc(obj):
    rt.SetOwnership(obj,False)

def NLL(r,s,b,d):
    # L_s+b/L_b
#    return r*s - d*math.log(1+r*s/b)
    # L_s+b
    if r*s+b<=0: return 0
#    return (r*s+b)-d*math.log(r*s+b)-gammaln(d+1)
    return (r*s+b)-d*math.log(r*s+b)

def poisErr(val):
    return rt.Math.gamma_quantile_c((1-0.6827)/2.,val+1,1)-val

def makeRHist(r,region,mass,entry,hist_data,hist_sig,ybounds,draw=False,hist0=None):
        rStr = "{:.3f}".format(r).replace(".","p").replace("-","n")
        hist_dataE = hist_data.Clone("data from entry"); ngc(hist_dataE)
        hist_R = hist_data.Clone("data + "+(str(r) if hist0 is not None else "")+"*signal"); ngc(hist_R)
        hist_NLL = hist_R.Clone("NLL "+hist_R.GetName()); ngc(hist_NLL)

        for iBin in range(1,hist_data.GetNbinsX()+1):
            hist_dataE.SetBinContent(iBin,entry["HIST"][iBin-1][3] if iBin-1 < len(entry["HIST"]) else 0)
            hist_dataE.SetBinError(iBin,poisErr(hist_dataE.GetBinContent(iBin)))
            hist_R.SetBinContent(iBin,entry["HIST"][iBin-1][2] if iBin-1 < len(entry["HIST"]) else 0)
            hist_NLL.SetBinContent(iBin,-1*entry["HIST"][iBin-1][4] if iBin-1 < len(entry["HIST"]) else 0)
            hist_NLL.SetBinError(iBin,0)

        # set the norm of hist_R
        coeff = sum(entry["COEF"])
        hist_R.Scale(coeff/hist_R.Integral())
        for iBin in range(1,hist_data.GetNbinsX()+1):
            hist_R.SetBinError(iBin,poisErr(hist_R.GetBinContent(iBin)))

        hist_sigr = hist_sig.Clone()
        hist_sigr.GetYaxis().SetTitle("residual")
        hist_sigr.Scale(r)
        hist_sigr.SetLineColor(rt.kCyan + 2)
        hist_sigr.SetMarkerColor(rt.kCyan + 2)

        hist_R.SetMarkerStyle(20)
        hist_R.SetMarkerColor(rt.kBlack)
        hist_R.SetLineColor(rt.kBlack)

#        # remove negative bins (and bins w/ very small # events from signal)
#        for iBin in range(hist_R.GetNbinsX()+1):
#            if hist_R.GetBinContent(iBin)<0.9:
#                hist_R.SetBinContent(iBin,0)
#                hist_R.SetBinError(iBin,0)

        hist_R.SetTitle(region+" "+mass+" "+rStr)
        hist_R.GetYaxis().SetRangeUser(ybounds[0],ybounds[1])
        
        if "high" in region:
            fitFunc = rt.TF1("fitFunc","[0]*pow(x/13000, -[1])",1500,6000); ngc(fitFunc)
            fitFunc.SetParameter(0,1)
            fitFunc.FixParameter(1,entry["PAR"][0])
        elif "low" in region:
            fitFunc = rt.TF1("fitFunc","[0]*pow(1 - x/13000, [2]) * pow(x/13000, -[1])",1500,6000); ngc(fitFunc)
            fitFunc.SetParameter(0,1)
            fitFunc.FixParameter(1,entry["PAR"][0])
            fitFunc.FixParameter(2,entry["PAR"][1])
        else:
            exit(0)
        # set the fit func normalization
        fitFunc.FixParameter(0,coeff/fitFunc.Integral(1500,6000)*hist_R.GetBinWidth(1))
        # get chisquare
        hist_R.Fit(fitFunc,"Q0")
#        fitFunc.SetChisquare(hist_R.Chisquare(fitFunc))

        hist_res = hist_R.Clone("res")
        hist_res.SetLineColor(rt.kRed)
        hist_res.SetMarkerColor(rt.kRed)
        for iBin in range(hist_res.GetNbinsX()+1):
            x = hist_res.GetXaxis().GetBinCenter(iBin)
            hist_res.SetBinContent(iBin,hist_R.GetBinContent(iBin)-fitFunc.Eval(x))

        if hist0 is not None:
            hist0["fit"].SetLineColor(rt.kMagenta+2)
            hist0["res"].SetLineColor(rt.kBlue)
            hist0["res"].SetMarkerColor(rt.kBlue)

        #hist_NLL.Reset()
        hist_NLL.SetName("hist_NLL")
        hist_NLL.SetTitle(";m_{T};deltaNLL")
        hist_NLL.SetStats(0)
        hist_NLLP = hist_NLL.Clone(); ngc(hist_NLLP)
        hist_NLLP.SetStats(0)
        hist_NLLP.SetLineColor(rt.kBlack)
        hist_NLLP.SetMarkerColor(rt.kBlack)
        hist_NLLP.SetMarkerStyle(20)
        for iBin in range(1,hist_NLLP.GetNbinsX()+1):
            x = hist_NLLP.GetXaxis().GetBinCenter(iBin)
            d = hist_data.GetBinContent(iBin)
            s = hist_sig.GetBinContent(iBin)
            b = fitFunc.Eval(x)
            hist_NLLP.SetBinContent(iBin, NLL(r,s,b,d))
            if verbose: print "DEBUGNLL: {} {} {} {} {} {} {}".format(iBin,r,s,b,d,hist_NLLP.GetBinContent(iBin),hist_NLLP.GetBinContent(iBin)-hist0["nllp"].GetBinContent(iBin) if hist0 is not None else "")
            hist_NLLP.SetBinError(iBin, 0)
        if hist0 is not None:
            hist_NLL.Add(hist0["nll"],-1)
            hist_NLLP.Add(hist0["nllp"],-1)
        nll = hist_NLL.Integral(1,hist_NLL.GetNbinsX())
        nllp = hist_NLLP.Integral(1,hist_NLLP.GetNbinsX())
        hist_NLL.GetYaxis().SetRangeUser(-5,10)
#        hist_NLL.Print("all")
        
        if draw:
            print("data: {}, sig: {}, d+s: {}, r: {}, chi2/ndf: {}, nll: {}".format(hist_data.Integral(),hist_sig.Integral(),hist_R.Integral(),r,fitFunc.GetChisquare()/fitFunc.GetNDF(),nll))

            c1 = rt.TCanvas("c1","c1",1000,1300)
            p1 = rt.TPad("p1",region+" "+mass+" r_{inj}="+str(r),0.0,0.5,1.0,1.0)
            p2 = rt.TPad("p2","res",0.0,0.25,1.0,0.5)
            p3 = rt.TPad("p3","NLL",0.0,0.0,1.0,0.25)
            p1.SetBottomMargin(0.0001)
            p1.SetBorderMode(0)
            p1.SetLogy(1)
            p2.SetTopMargin(0.0001)
            p2.SetBorderMode(0)
            p2.SetBottomMargin(0.0001)
            p3.SetTopMargin(0.0001)
            p3.SetBorderMode(0)

            p1.cd()
            # get y bounds on first try
            if hist0 is None:
                hist_data.Draw()
                hist_R.Draw("same")
            else:
                hist_R.Draw()
                hist_data.Draw("same")
                hist_R.Draw("same")
                fitFunc.Draw("same")
                hist0["fit"].Draw("same")

            p2.cd()
            if "SVJ" in region: hist_sigr.GetYaxis().SetRangeUser(-10,10)
            else: hist_sigr.GetYaxis().SetRangeUser(-50,50)
            hist_sigr.Draw("p")
            hist_res.Draw("p same")
            if hist0 is not None: hist0["res"].Draw("p same")

            p3.cd()
            hist_NLL.Draw("p")
            hist_NLLP.Draw("p same")
            c1.cd()
            p1.Draw()
            if hist0 is None:
                ybounds[0] = math.pow(10,p1.GetUymin())
                ybounds[1] = math.pow(10,p1.GetUymax())
            p2.Draw("same")
            p3.Draw("same")
            #c1.SaveAs("plots_BestRFit_res_combine/dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")
            c1.SaveAs("plots_BestRFit_res_combine/sigMod_dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")

        return nll, {"nll": hist_NLL, "nllp": hist_NLLP, "fit": fitFunc, "res": hist_res}

for region in regions:
    for mass in masses:
        print(region, mass)
        scan = []
        with open("log_mdf_debug_"+str(mass)+"_all.log",'r') as logfile:
            found_region = False
            entry = {}
            for line in logfile:
                linesplit = line.split()
                if len(linesplit)==2 and linesplit[0]=="REGION:":
                    if linesplit[1]==region: found_region = True
                    else: found_region = False
                if not found_region: continue
                if line.startswith("DEBUG"):
                    ltype = linesplit[0].replace("DEBUGNLL","").replace(":","")
                    # marks the beginning of a new entry
                    if ltype=="COEF":
                        if len(entry)>0:
                            entry['r'] = entry["RET"][0]
                            scan.append(entry)
                        entry = {}
                    tmp = [float(x) for x in linesplit[1:] if not x[0].isalpha()]
                    if ltype=="HIST":
                        if ltype not in entry: entry[ltype] = []
                        entry[ltype].append(tmp)
                    else:
                        entry[ltype] = tmp

        #rootFile = rt.TFile.Open("datacard2.root","r") # make sure to also change what the plots are saved as!
        rootFile = rt.TFile.Open("datacard3.root","r")

        tDirName = region+"_2018"
        
        dataName = "data_obs" # TFirectoryFile in rootFile
        sigName = "SVJ_mZprime"+mass+"_mDark20_rinv03_alphapeak" # TH1F in dataName

        hist_data = rootFile.Get(tDirName+"/"+dataName); ngc(hist_data)
        hist_sig = rootFile.Get(tDirName+"/"+sigName); ngc(hist_sig)

        print("nEvents Data: {} {}".format(hist_data.Integral(),hist_data.GetEntries()))
        print("nEvents Sig: {} {}".format(hist_sig.Integral(),hist_sig.GetEntries()))

        hist_data.SetMarkerStyle(24)
        hist_data.SetMarkerColor(rt.kBlue)
        hist_data.SetLineColor(rt.kBlue)

        ybounds = [-1,-1]
        entry0 = None
        for entry in scan:
            if abs(entry['r']+0.015)<1e-10:
                entry0 = entry
                break
        nll0, hist0 = makeRHist(0,region,mass,entry0,hist_data,hist_sig,ybounds,draw=True)
        nlls = []
        rvals = []
        for entry in scan:
            rvals.append(entry['r'])
            nll, hist = makeRHist(entry['r'],region,mass,entry,hist_data,hist_sig,ybounds,draw=True,hist0=hist0)
            nlls.append(nll)

        # compute deltaNLL
        rvals = np.array(rvals)
        nlls = np.array(nlls)
        print "min(deltaNLL) = "+str(np.min(nlls))+" at r = "+str(rvals[np.argmin(nlls)])

        # make graph
        g = rt.TGraph(len(rvals),rvals,nlls)
        g.SetTitle("")
        g.GetXaxis().SetTitle("r")
        g.GetYaxis().SetTitle("deltaNLL")
        c1 = rt.TCanvas("c1","c1",1000,700)
        p1 = rt.TPad("p1","deltaNLL "+region+" "+mass,0.0,0.0,1.0,1.0)
        p1.SetBorderMode(0)
        p1.cd()
        g.Draw("ap")
        c1.cd()
        p1.Draw()
        c1.SaveAs("plots_BestRFit_res_combine/deltaNLL_sigMod_dataAndSigWithFits_"+region+"_"+mass+".png")












