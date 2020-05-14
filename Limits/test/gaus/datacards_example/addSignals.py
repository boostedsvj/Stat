import os
import ROOT as r
from array import array
from string import Template

def fill_template(inname, outname, **kwargs):
    if inname==outname:
        raise ValueError("Attempt to overwrite template file: "+inname)
    with open(inname,'r') as temp:
        old_lines = Template(temp.read())
        new_lines = old_lines.substitute(**kwargs)
    with open(outname,'w') as temp:
        temp.write(new_lines)

r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

fname = "dijet_combine_qq_1900_lumi-137.500_PFDijet2017"

infile = r.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/gaus/datacards_example/"+fname+".root","READ")

w = infile.Get("wPFDijet2017")

th1x = w.var("th1x")
bins = th1x.getBinning().array()
bins.SetSize(th1x.getBinning().numBins())
bins = array('d',bins)

# gaussian fit to original signal: (not very good fit)
# Constant     2.94156e+02
# Mean         4.21729e+00
# Sigma        1.82464e+00
# vary sigma, keep others
# also check SVJ-like mean: (2141-1500)/(6000-1500)*42 = 5.983
# use normalized gaussian to maintain integral as signal rate

gaus = r.TF1("gaus","gausn",min(bins),max(bins))
gaus.SetParameter(0, 1140.858)
means = [4.21729, 14.987]
sigmas = [1.0,2.0,3.0,4.0,5.0]
for iMean,mean in enumerate(means):
    for iSigma,sigma in enumerate(sigmas):
        print iMean,iSigma
        gaus.SetParameter(1,mean)
        gaus.SetParameter(2,sigma)
        hname = "gaus_m{}_s{}".format(iMean,iSigma)
        hist = r.TH1F(hname,hname,len(bins)-1,bins)
        for b in range(1,len(bins)):
            binc = gaus.Eval((bins[b-1]+bins[b])/2.)
            hist.SetBinContent(b,binc)
            hist.SetBinError(b,r.Math.gamma_quantile_c((1-0.6827)/2.,binc+1,1)-binc)
        hyield = hist.Integral()
        # make datacard
        fill_template(
            "dijet_combine_sig_template.txt",
            "dijet_combine_"+hname+".txt",
            signame = hname,
            sigrate = hyield
        )
        # make pdf
        pdf = r.RooDataHist(hname,hname, r.RooArgList(th1x), hist, 1.)
        # import in workspace
        getattr(w,"import")(pdf)

# write out modified workspace
mfname = fname+"_mod.root"
os.remove(mfname)
w.writeToFile(mfname)
