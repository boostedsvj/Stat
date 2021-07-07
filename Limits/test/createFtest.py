import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

channels = ["highCut", "lowCut"]
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', dest='idir', type=str, default= "root://cmseos.fnal.gov//store/user/cfallon/datacards_2Feb/",help='Location of fit output files ws_allFits_{}.root, fitResults_{}.root')
parser.add_argument("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_argument("-x", "--useChi2", dest="useChi2", action="store_true", default=False, help="use chi2 in F-test")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False)
parser.add_argument("-p", "--noplots", dest="doplots", action="store_false", default=True)
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.ftest import *

years = ['2018']
ch_year = []

print "====> CHANNELS: ", channels
for y in years:
    channels_years = [ch + '_' + y for ch in channels ]
    ch_year= ch_year + channels_years
    

print "====> CHANNELS + YEAR: ", ch_year

for ch in ch_year:
    getCard(ch, opt.idir, opt.bias, opt.useChi2, opt.verbose, opt.doplots)


