import ROOT
import os, sys
import optparse
from Stat.Limits.datacardsUnified import *
#changed to no longer need the settings.py file
# -signal parameters are now command-line input
# -Fisher testing is only done on the baseline (3000, 20, 03, peak) signal
# -channels list is now created here

channels = ["lowSVJ2", "highSVJ2", "highCut", "lowCut"]
sigpoints = []
print "====> CHANNELS: ", channels

path = "/uscms_data/d3/cfallon/SVJ/bs7/CMSSW_10_2_13/src/Stat/Limits/test"

usage = 'usage: %prog -p histosPath -o outputFile'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--input', dest='ifile', type='string', default= "root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig4/sigfull/",help='Where can I find input histos? Default is new (24 July 2020) files created by Kevin')
parser.add_option("-d","--outdir",dest="outdir",type="string",default="outdir",help="Name of the output directory where to store datacards. Default is outdir")
parser.add_option("-m","--mode",dest="mode",type="string",default="hist",help="Kind of shape analysis: parametric fit or fit to histos?. Default is hist")
parser.add_option("-Z", "--zMass", dest="mZ", type='str',help="str: Mass [GeV] of the Z' in MC signal. Default = '3100', range: [500, 4400] in steps of 100, inclusive", default='3100')
parser.add_option("-D", "--dMass", dest="mD", type='str', help="str: Mass [GeV] of dark quarks in MC singal. Default = '20'", default = '20')
parser.add_option("-R", "--rInv", dest="rI", type='str', help="str: Fraction of invisible particles in MC signal. Default = '03'", default = '03')
parser.add_option("-A", "--aDark", dest="aD", type='str', help="str: alphaDark value in MC signal. Default = 'peak'. Options: 'low', 'peak', 'high'", default = "peak")
parser.add_option("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_option("-u","--unblind",dest="unblind",action='store_true', default=False)
parser.add_option("-s", "--doSys",dest="doSys",action='store_false', default=True)

(opt, args) = parser.parse_args()
sys.argv.append('-b')

if len(opt.rI) > 1:
	ifilename = opt.ifile + "datacard_final_SVJ_"+opt.mZ+"_"+opt.mD+"_"+opt.rI[0]+"."+opt.rI[1]+"_"+opt.aD+".root"
else:
	ifilename = opt.ifile + "datacard_final_SVJ_"+opt.mZ+"_"+opt.mD+"_"+opt.rI+"_"+opt.aD+".root"
outdir = opt.outdir
mode = opt.mode
unblind = opt.unblind

bias = opt.bias

signals = []

sigpoints.append([opt.mZ, opt.mD, opt.rI, opt.aD])

print "====> CHANNELS: ", channels

for p in sigpoints:

    mZprime=p[0]
    mDark=p[1]
    rinv=p[2]
    alpha=p[3]

    print "Creating datacards for mZprime = ", mZprime, " GeV, mDark = ", mDark, " GeV, rinv = ", rinv, " , alpha = ", alpha
    signal  = "SVJ_mZprime%s_mDark%s_rinv%s_alpha%s" % (mZprime, mDark, rinv, alpha) 
    #if(bias):     signal  = "SVJ_mZprime%s_mDark%s_rinv%s_alpha%s_bias" % (mZprime, mDark, rinv, alpha) 
    signals.append(signal)



#signals = ["SVJ_mZprime3000_mDark20_rinv03_alphapeak"]



    

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

cmd = "rm Efficiencies.txt"
os.system(cmd)

effs = {}


ch_eff = ["lowSVJ0", "lowSVJ1", "lowSVJ2", "highSVJ0", "highSVJ1", "highSVJ2"]
ch_eff_year = []
for y in years:
    ch_eff_years = [ch + '_' + y for ch in ch_eff ]
    ch_eff_year= ch_eff_year + ch_eff_years


for s in signals:

    effs[s] = getEfficiency(s, ch_eff_year, ifilename)

y = json.dumps(effs)
outname =  "Efficiencies.txt"

efile = open(outname, 'w')
efile.write(y)
efile.close()

for s in signals:

    #if s=="SVJ_mZprime3000_mDark20_rinv03_alphapeak":
    #    doModelling = True
    #else:
    #    doModelling = False
    doModelling = True # need to evaluate Fisher test for every batch
    # as an alternatative, we could set things up to F-test once, then save
    # those functions and import them to every submission, TODO?
    for ch in ch_year:
        getCard(s, ch, ifilename, outdir, doModelling, mode, bias, True, opt.doSys)


