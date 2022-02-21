import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

channels = ['bsvj_0p0', 'bsvj_0p1', 'bsvj_0p2', 'bsvj_0p3']
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input',   dest='idir',   type=str,             default= "/home/saranabili/LimitStudy/September16_latestSetUp/CMSSW_10_2_13/src/Stat/Limits/test",help='the file for bsvj is in hepcms')
parser.add_argument("-t", "--test",    dest="bias",    action="store_true",  default=False)
parser.add_argument("-x", "--useChi2", dest="useChi2", action="store_true",  default=True, help="use chi2 in F-test")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",  default=True)
parser.add_argument("-p", "--noplots", dest="doplots", action="store_false", default=True)
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.ftest import *


print "====> CHANNELS: ", channels
ch_year = ['bsvj_0p0', 'bsvj_0p1', 'bsvj_0p2', 'bsvj_0p3']
print "====> CHANNELS + YEAR: ", ch_year

for ch in ch_year:
    print"these are the parameters:"
    print"	ch: ",ch
    print"      opt.idir: ",opt.idir
    print"      opt.bias: ",opt.bias
    print"      opt.useChi2: ",opt.useChi2
    print"      opt.verbose: ",opt.verbose
    print"      opt.doplots: ",opt.doplots
    getCard(ch, opt.idir, opt.bias, opt.useChi2, opt.verbose, opt.doplots)
