import os, sys

region = sys.argv[1]

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

fname = "ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_"+region+"_2018_template"

infile = r.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/gaus/svj/"+fname+".root","READ")

w = infile.Get("SVJ")

mH = w.var("mH")
bins = mH.getBinning().array()
bins.SetSize(mH.getBinning().numBins())
bins = array('d',bins)

# gaussian fit to original signal: (not very good fit)
# Constant     2.94156e+02
# Mean         4.21729e+00
# Sigma        1.82464e+00
# vary sigma, keep others
# use normalized gaussian to maintain integral as signal rate

gaus = r.TF1("gaus","gausn",min(bins),max(bins))
norms = {
"highCut" : 156.003928,
"highSVJ2": 120.191287,
"lowCut": 188.853831,
"lowSVJ2": 145.396926,
}
gaus.SetParameter(0, norms[region])
means = [2141.1]
sigfac = (6000.-1500.)/42.
sigmas = [1.0*sigfac,2.0*sigfac,3.0*sigfac,4.0*sigfac,5.0*sigfac]
for iMean,mean in enumerate(means):
    for iSigma,sigma in enumerate(sigmas):
        print iMean,iSigma
        gaus.SetParameter(1,mean)
        gaus.SetParameter(2,sigma)
        hname = "gaus_m{}_s{}".format(iMean,iSigma)
        hist = r.TH1F(hname,hname,len(bins)-1,bins)
        for b in range(1,len(bins)):
            center = (bins[b-1]+bins[b])/2.
            # include bin width
            binc = gaus.Eval(center)*(bins[b]-bins[b-1])
            hist.SetBinContent(b,binc)
            hist.SetBinError(b,r.Math.gamma_quantile_c((1-0.6827)/2.,binc+1,1)-binc)
        hyield = hist.Integral()
        # make datacard
        fill_template(
            "svj_"+region+"_2018_template_bias.txt",
            "svj_"+region+"_2018_"+hname+".txt",
            signame = hname,
            sigrate = hyield
        )
        # make pdf
        pdf = r.RooDataHist(hname,hname, r.RooArgList(mH), hist, 1.)
        # import in workspace
        getattr(w,"import")(pdf)

# write out modified workspace
mfname = fname+"_mod.root"
os.remove(mfname)
w.writeToFile(mfname)
