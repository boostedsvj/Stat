import os, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

channels = ['bsvj']
sigpoints = []

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', dest='ifile', type=str, default= "/home/saranabili/LimitStudy/CMSSW_10_2_13/src/Stat/Limits/test/",help='the file for bsvj is in hepcms')
parser.add_argument("-t", "--test", dest="bias", action="store_true", default=False)
parser.add_argument("-n","--npool",dest="npool",type=int,default=0,help="number of parallel processes for brute force method (0 = parallelism disabled)")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=True)
#parser.add_argument("-I", "--initvals", dest="initvals", type=float, default=[-8.0,-0.8,-0.1,0.1,0.8,8.0], nargs='+', help="list of allowed initial values for brute force method")
parser.add_argument("-I", "--initvals", dest="initvals", type=float, default=[-10.0,-1.0,-0.1,0.1,1.0,10.0], nargs='+', help="list of allowed initial values for brute force method")
opt = parser.parse_args()
sys.argv.append('-b')

import ROOT
from Stat.Limits.allFits import *

#ifilename = opt.ifile + "tmva_MThisto_June25.root"
#ifilename = opt.ifile + "tmva_MThisto_August6_systematicup.root"
#ifilename = opt.ifile + "qcd_mz250_tight.root"
#ifilename = opt.ifile + "histograms_3masspoints_qcdttjets_Sep21.root"
#ifilename = opt.ifile + "Sept27_3masspoints_qcd.root"
#ifilename = opt.ifile + "Sept27_3masspoints_qcd_bdt50.root"
#ifilename = opt.ifile + "Sept28_3masspoints_bdt50.root"
#ifilename = opt.ifile + "tmva_MThisto_August13_hacksys.root"
#ifilename = opt.ifile + "Sept29_qcdttjets_3mz.root"
#ifilename =  opt.ifile + "Oct04_3mass_allbkg.root"
#ifilename = opt.ifile + "histograms_3masspoints_allbkgs_highermassrange_Oct07.root"
#ifilename = opt.ifile + "histograms_Oct13.root"
ifilename = opt.ifile + "test.root" #the last file with the good results: Sept29_qcdttjets_3mz.root

signals = []

sigpoints.append(["250","10","03","peak"])
sigpoints.append(["300","10","03","peak"])
sigpoints.append(["350","10","03","peak"])
sigpoints.append(["400","10","03","peak"])
sigpoints.append(["450","10","03","peak"])

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
    print"r_years are: ", r_years

    print "years are:  ", years
    


print "====> CHANNELS: ", channels
ch_year = ['bsvj']

print "====> CHANNELS + YEAR: ", ch_year

for s in signals:
    for ch in ch_year:
        getCard(s, ch, ifilename, opt.npool, opt.initvals, opt.bias, opt.verbose)
        print"=====> these are the characteristics:"
	print"		s is: ",s
	print"		ch is: ",ch
	print"		ifilename is:  ",ifilename
	print"		npool is: ",opt.npool
	print"		bias is : ",opt.bias
	print"		verbose is: ",opt.verbose
	print"		initial values are: ", opt.initvals
	print"**************************************************************************************"


