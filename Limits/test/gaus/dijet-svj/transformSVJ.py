import os,sys

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

infile = r.TFile.Open("dijet_combine_qq_1900_lumi-137.500_PFDijet2017.root","READ")

w = infile.Get("wPFDijet2017")

th1x = w.var("th1x")
bins = th1x.getBinning().array()
bins.SetSize(th1x.getBinning().numBins())
bins = array('d',bins)

yields = {
"dijet": 6557850,
"highCut" : 10941,
"highSVJ2": 166,
"lowCut": 74920,
"lowSVJ2": 1225,
}

sigyields = {
"dijet": 1140.858,
"highCut" : 156.003928,
"highSVJ2": 120.191287,
"lowCut": 188.853831,
"lowSVJ2": 145.396926,
}

sigratios = {key : sigyields[key]/yields[key] for key in yields}

gaus = r.TF1("gaus","gausn",min(bins),max(bins))
# use SVJ-like mean
mean = 5.9827
sigmas = [2.0,4.0]

data_obs = w.data("data_obs")
if region=="dijet":
    data_obs_gen = data_obs
else:
    # generate toy data
    pdf_data = r.RooHistPdf("pdf_data","pdf_data",r.RooArgSet(th1x),data_obs)
    data_obs_gen = pdf_data.generateBinned(r.RooArgSet(th1x),yields[region])

ofname = "dijet_combine_{}_PFDijet2017.root".format(region)

# make new workspace, import contents from prev (but take new data)
w2 = r.RooWorkspace("wPFDijet2018")
getattr(w2,"import")(w.allVars())
getattr(w2,"import")(w.allPdfs())
getattr(w2,"import")(data_obs_gen,r.RooFit.Rename("data_obs"))

for sregion,ratio in sigratios.iteritems():
    # use normalized gaussian to maintain integral as signal rate
    gaus.SetParameter(0, ratio*yields[region])
    for sigma in sigmas:
        gaus.SetParameter(1,mean)
        gaus.SetParameter(2,sigma)
        hname = "{}_gaus_r{}_s{}".format(region,sregion,sigma)
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
            sigrate = hyield,
            datrate = yields[region],
            wsfile = ofname
        )
        # make pdf
        pdf = r.RooDataHist(hname,hname, r.RooArgList(th1x), hist, 1.)
        # import in workspace
        getattr(w2,"import")(pdf)

# write out modified workspace
if os.path.exists(ofname): os.remove(ofname)
w2.writeToFile(ofname)
