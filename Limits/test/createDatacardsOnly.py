import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

#channels = ['bsvj_0p0']
#channels = ['bsvj_0p1']
#channels = ['bsvj_0p2']
#channels = ['bsvj_0p3']
channels = ['bsvj_0p0', 'bsvj_0p1', 'bsvj_0p2']
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', dest='ifile', type=str, default= "/home/saranabili/LimitStudy/September16_latestSetUp/CMSSW_10_2_13/src/Stat/Limits/test/",help='Location of F-test output files ws_{}.root')
parser.add_argument('-w', '--workspaceDir', dest='workspaceDir', type=str, default="/home/saranabili/LimitStudy/September16_latestSetUp/CMSSW_10_2_13/src/Stat/Limits/test/",help='Location of F-test output files ws_{}.root')
parser.add_argument("-m","--mode",dest="mode",type=str,default="template",help="Kind of shape analysis: parametric fit or fit to histos?")
parser.add_argument("-Z", "--zMass", dest="mZ", type=str,help="Mass [GeV] of the Z' in MC signal. range: [250, 450] in steps of 50, inclusive", default='350')
#parser.add_argument("-Z", "--zMass", dest="mZ", type=str, help="Mass [GeV] of the Z' in MC signal. range: [250, 450] in steps of 50, inclusive", default = '[250 300]')
parser.add_argument("-D", "--dMass", dest="mD", type=str, help="Mass [GeV] of dark quarks in MC signal", default = '10')
parser.add_argument("-R", "--rInv", dest="rI", type=str, help="Fraction of invisible particles in MC signal", default = '03')
parser.add_argument("-A", "--aDark", dest="aD", type=str, help="alphaDark value in MC signal. Options: 'low', 'peak', 'high'", default = "peak")
parser.add_argument("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_argument("-s", "--noSys",dest="doSys",action='store_false', default=False)
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.datacardsOnly import *

#ifilename = opt.ifile + "mt_Feb01_bdtcuts0p4.root"
#ifilename = opt.ifile + "mt_Feb03_bdtcuts0p04.root"
ifilename = opt.ifile + "datacard_final_SVJ_"+opt.mZ+"_"+opt.mD+"_"+(opt.rI if len(opt.rI)==1 else opt.rI[0]+"."+opt.rI[1:])+"_"+opt.aD+".root"
#ifilename = opt.ifile + "datacard_final_SVJ_250_10_0.3_peak.root"
#ifilename = opt.ifile + "mt_Jan23_bdtcuts0p40p45.root"
#ifilename = opt.ifile + "mt_Jan23_bdtcuts0p50p60.root"
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
 
#ch_year = ['bsvj_0p0'] 
#ch_year = ['bsvj_0p1']
#ch_year = ['bsvj_0p2']
#ch_year = ['bsvj_0p3']
ch_year = ['bsvj_0p0', 'bsvj_0p1', 'bsvj_0p2']
print "====> CHANNELS + YEAR: ", ch_year

cmd = "rm ws.root"
os.system(cmd)

for s in signals:
    doModelling = True # need to evaluate Fisher test for every batch
    for ch in ch_year:
        getCard(s, ch, ifilename, opt.workspaceDir, doModelling, opt.mode, opt.bias, True, opt.doSys)
