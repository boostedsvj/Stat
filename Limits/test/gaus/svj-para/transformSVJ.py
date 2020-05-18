import sys

region = sys.argv[1]

import ROOT as r
import array
r.gSystem.Load('libHiggsAnalysisCombinedLimit')

npars = {
    "lowCut": 2,
    "highCut": 3,
    "lowSVJ2": 2,
    "highSVJ2": 1,
}

f_ws = r.TFile.Open('ws_SVJ_mZprime3000_mDark20_rinv03_alphapeak_'+region+'_2018_template.root')
workspace = f_ws.Get('SVJ')

pdf = workspace.pdf('Bkg_'+region+'_2018')
pars = r.RooArgList()
for i in range(npars[region]):
    pname = region+'_p'+str(i+1)+'_'+str(npars[region])
    pars.add(workspace.var(pname))

pars.Print('v')

mT = workspace.var('mH')
binning = mT.getBinning('')
x = array.array('d',[])
for i in range(0,binning.numBins()):
    x.append(binning.binLow(i))
x.append(binning.binHigh(binning.numBins()-1))
print "binning:", x

shape = r.TH1D('empty','empty',len(x)-1,x)

new_pdf = r.RooParametricShapeBinPdf('Bkg2_'+region+'_2018','Bkg2_'+region+'_2018',pdf,mT,pars,shape)
getattr(workspace,'import')(new_pdf,r.RooFit.RecycleConflictNodes())

workspace.Print('v')
workspace.writeToFile('ws2_SVJ_mZprime3000_mDark20_rinv03_alphapeak_'+region+'_2018_template.root',True)
