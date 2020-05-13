import sys
region = sys.argv[2]

import ROOT as rt
import random
rt.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
rt.gStyle.SetOptStat(111111)
rt.gStyle.SetOptFit(1011)
import math
import numpy as np
from array import array
from scipy.special import gammaln

# plot dataObs + r*Sig for atleast r = 0 and r = min(deltaNLL)

#regions = ["lowCut","lowSVJ2","highCut","highSVJ2"]
regions = [region]

endOfScriptOutput = "Bins with err/val > 0.5\nRegion\tMass\tiBin\tValue\tError\n"

def NLL(r,s,b,d):
    # L_s+b/L_b
#    return r*s - d*math.log(1+r*s/b)
    # L_s+b
    if r*s+b<=0: return 0
#    return (r*s+b)-d*math.log(r*s+b)-gammaln(d+1)
    return (r*s+b)-d*math.log(r*s+b)

for region in regions:
#    for mass in ["2000","3000","4000"]:
    for mass in ["3000"]:
        print(region, mass)
        #rootFile = rt.TFile.Open("datacard2.root","r") # make sure to also change what the plots are saved as!
        rootFile = rt.TFile.Open("datacard3.root","r")

        tDirName = region+"_2018"
        
        dataName = "data_obs" # TFirectoryFile in rootFile
        sigName = "SVJ_mZprime"+mass+"_mDark20_rinv03_alphapeak" # TH1F in dataName

        hist_data = rootFile.Get(tDirName+"/"+dataName)
        hist_sig = rootFile.Get(tDirName+"/"+sigName)
        #print("Signal {} {}".format(region, mass))
        #print("iBin\tval\terr\tval/err")
        for iBin in range(hist_sig.GetNbinsX()):
            if hist_sig.GetBinContent(iBin) != 0:
                if (hist_sig.GetBinError(iBin)/hist_sig.GetBinContent(iBin) > 0.5):
                    endOfScriptOutput += ("{}\t{}\t{}\t{}\t{}\t{}\n".format(region, mass, iBin, hist_sig.GetBinContent(iBin),hist_sig.GetBinError(iBin), hist_sig.GetBinError(iBin)/hist_sig.GetBinContent(iBin)))

        print("nEvents Data: {} {}".format(hist_data.Integral(),hist_data.GetEntries()))
        print("nEvents Sig: {} {}".format(hist_sig.Integral(),hist_sig.GetEntries()))

        rvals = np.linspace(-3,3,99)
        nlls = []
        pads = []
        for r in rvals:
#        for r in [0, -0.46, -0.11, -1.91, -0.07]:
            rStr = "{:.3f}".format(r).replace(".","p").replace("-","n")
            hist_R = hist_data.Clone()
            hist_R.Add(hist_data, hist_sig, 1, r)

            # remove negative bins (and bins w/ very small # events from signal)
            for iBin in range(hist_R.GetNbinsX()+1):
                if hist_R.GetBinContent(iBin)<0.5:
                    hist_R.SetBinContent(iBin,0)
                    hist_R.SetBinError(iBin,0)

            #if hist_R.GetMinimum() < 0:
            #    print("NEGATIVE BIN: {} {} {} {}".format(region, mass, r, hist_R.GetMinimum()))
            #    continue
            hist_R.SetTitle(region+" "+mass+" "+rStr)
            hist_NLL = hist_R.Clone()
            
            if "high" in region:
                fitFunc = rt.TF1("fitFunc","[0]*pow(x/13000, -[1])",1500,6000)
                if "SVJ" in region:
                    fitFunc.SetParameter(0,4.3e-6)
                    fitFunc.SetParameter(1,7.3)
                else:
                    fitFunc.SetParameter(0,0.0018)
                    fitFunc.SetParameter(1,6.4)
                fitFunc.SetParLimits(0, 0.5e-6, 0.5)
                fitFunc.SetParLimits(1, 1, 30)
            elif "low" in region:
                fitFunc = rt.TF1("fitFunc","[0]*pow(1 - x/13000, [2]) * pow(x/13000, -[1])",1500,6000)
                if "SVJ" in region:
                    fitFunc.SetParameter(0,4.4e-5)
                    fitFunc.SetParameter(1,7.5)
                    fitFunc.SetParameter(2,4.2)
                else:
                    fitFunc.SetParameter(0,0.08)
                    fitFunc.SetParameter(1,6.0)
                    fitFunc.SetParameter(2,6.7)
                fitFunc.SetParLimits(0, 0.5e-6, 0.5)
                fitFunc.SetParLimits(1, 1, 30)
                fitFunc.SetParLimits(2, 1, 16)
            else:
                exit(0)
            hist_R.Fit(fitFunc,"Q")

            #hist_NLL.Reset()
            hist_NLL.SetName("hist_NLL")
            hist_NLL.SetTitle(";m_{T};NLL")
            hist_NLL.SetStats(0)
            for iBin in range(hist_NLL.GetNbinsX()+1):
                x = hist_NLL.GetXaxis().GetBinCenter(iBin)
                d = hist_data.GetBinContent(iBin)
                s = hist_sig.GetBinContent(iBin)
                b = fitFunc.Eval(x)
                hist_NLL.SetBinContent(iBin, NLL(r,s,b,d))
                hist_NLL.SetBinError(iBin, 0)
            nlls.append(hist_NLL.Integral(-1,-1))
            
            print("data: {}, sig: {}, d+s: {}, r: {}, chi2/ndf: {}, nll: {}".format(hist_data.Integral(),hist_sig.Integral(),hist_R.Integral(),r,fitFunc.GetChisquare()/fitFunc.GetNDF(),nlls[-1]))

            c1 = rt.TCanvas("c1","c1",1000,1000)
            p1 = rt.TPad("p1",region+" "+mass+" r_{inj}="+str(r),0.0,0.35,1.0,1.0)
            p2 = rt.TPad("p2","NLL",0.0,0.0,1.0,0.35)
            p1.SetBottomMargin(0.0001)
            p1.SetBorderMode(0)
            p1.SetLogy(1)
            p2.SetTopMargin(0.0001)
            p2.SetBorderMode(0)
            p1.cd()
            hist_R.Draw()
            p2.cd()
            hist_NLL.Draw("p")
            c1.cd()
            p1.Draw()
            p2.Draw("same")
            #c1.SaveAs("plots_BestRFit_res/dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")
            c1.SaveAs("plots_BestRFit_res/sigMod_dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")
            pads.extend([c1,p1,p2])

#            p1.cd()
#            p1.SetLogy(0)
#            hist_R.Draw()
#            p2.cd()
#            hist_NLL.Draw("p")
#            c1.cd()
#            p1.Draw()
#            p2.Draw("same")
#            #c1.SaveAs("plots_BestRFit_res/Linear_dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")
#            c1.SaveAs("plots_BestRFit_res/LinearsigMod_dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")

        ind0 = np.where(rvals==0.)[0][0]
        # compute deltaNLL
        nlls = np.array(nlls)
        nlls = nlls - nlls[ind0]
        print "min(deltaNLL) = "+str(np.min(nlls))+" at r = "+str(rvals[np.argmin(nlls)])

        # make graph
        g = rt.TGraph(len(rvals),rvals,nlls)
        g.SetTitle("")
        g.GetXaxis().SetTitle("r")
        g.GetYaxis().SetTitle("deltaNLL")
        c1 = rt.TCanvas("c1","c1",1000,700)
        p1 = rt.TPad("p1","deltaNLL "+region+" "+mass+" r_{inj}="+str(r),0.0,0.0,1.0,1.0)
        p1.SetBorderMode(0)
        p1.cd()
        g.Draw("ap")
        c1.cd()
        p1.Draw()
        c1.SaveAs("plots_BestRFit_res/deltaNLL_sigMod_dataAndSigWithFits_"+region+"_"+mass+"_"+rStr+"_chi2.png")



print(endOfScriptOutput)









