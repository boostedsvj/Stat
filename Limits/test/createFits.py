import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

#channels = ["lowSVJ2", "highSVJ2", "highCut", "lowCut"]
channels = ["highCut", "lowCut"]
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', dest='ifile', type=str, default= "root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig4/sigfull/",help='Where can I find input histos? trig4/sigfull = new (24 July 2020) files created by Kevin')
parser.add_argument("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_argument("-n","--npool",dest="npool",type=int,default=0,help="number of parallel processes for brute force method (0 = parallelism disabled)")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False)
parser.add_argument("-I", "--initvals", dest="initvals", type=float, default=[-10.0,-1.0,-0.1,0.1,1.0,10.0], nargs='+', help="list of allowed initial values for brute force method")
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.allFits import *

ifilename = opt.ifile + "datacard_final_SVJ_2900_20_0.3_peak.root"

signals = []

sigpoints.append(["2900","20","03","peak"])

for p in sigpoints:

    mZprime=p[0]
    mDark=p[1]
    rinv=p[2]
    alpha=p[3]

    print "Creating datacards for mZprime = ", mZprime, " GeV, mDark = ", mDark, " GeV, rinv = ", rinv, " , alpha = ", alpha
    signal  = "SVJ_mZprime%s_mDark%s_rinv%s_alpha%s" % (mZprime, mDark, rinv, alpha) 
    signals.append(signal)

try:
    ifile = ROOT.TFile.Open(ifilename)
except IOError:
    print "Cannot open ", ifilename
else:
    print "Opening file ",  ifilename
    ifile.cd()
    
    r = ROOT.gDirectory.GetListOfKeys()[0]
    
    r_years = [r.ReadObj().GetName()[-4:] for r in ROOT.gDirectory.GetListOfKeys() ]
    
    years =  list(set(r_years))
    

ch_year = []

print "====> CHANNELS: ", channels
for y in years:
    channels_years = [ch + '_' + y for ch in channels ]
    ch_year= ch_year + channels_years
    

print "====> CHANNELS + YEAR: ", ch_year

for s in signals:
    for ch in ch_year:
        getCard(s, ch, ifilename, opt.npool, opt.initvals, opt.bias, opt.verbose)


