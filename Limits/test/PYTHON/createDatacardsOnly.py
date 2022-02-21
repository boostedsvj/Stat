import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

channels = ["lowSVJ2", "highSVJ2", "highCut", "lowCut"]
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', dest='ifile', type=str, default= "root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig7/sigfull/",help='Where can I find input histos? trig7/sigfull = new (13 Sep 2021) files created by Kevin')
parser.add_argument('-w', '--workspaceDir', dest='workspaceDir', type=str, default="root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Limits/datacards_Sep13/",help='Location of F-test output files ws_{}.root')
parser.add_argument("-m","--mode",dest="mode",type=str,default="hist",help="Kind of shape analysis: parametric fit or fit to histos?")
parser.add_argument("-Z", "--zMass", dest="mZ", type=str,help="Mass [GeV] of the Z' in MC signal. range: [500, 4400] in steps of 100, inclusive", default='2900')
parser.add_argument("-D", "--dMass", dest="mD", type=str, help="Mass [GeV] of dark quarks in MC signal", default = '20')
parser.add_argument("-R", "--rInv", dest="rI", type=str, help="Fraction of invisible particles in MC signal", default = '03')
parser.add_argument("-A", "--aDark", dest="aD", type=str, help="alphaDark value in MC signal. Options: 'low', 'peak', 'high'", default = "peak")
parser.add_argument("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_argument("-s", "--noSys",dest="doSys",action='store_false', default=True)
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.datacardsOnly import *

ifilename = opt.ifile + "datacard_final_SVJ_"+opt.mZ+"_"+opt.mD+"_"+(opt.rI if len(opt.rI)==1 else opt.rI[0]+"."+opt.rI[1:])+"_"+opt.aD+".root"

signals = []

sigpoints.append([opt.mZ, opt.mD, opt.rI, opt.aD])

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

cmd = "rm ws.root"
os.system(cmd)

for s in signals:
    doModelling = True # need to evaluate Fisher test for every batch
    for ch in ch_year:
        getCard(s, ch, ifilename, opt.workspaceDir, doModelling, opt.mode, opt.bias, True, opt.doSys)
