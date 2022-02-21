#! /usr/bin/env python
import sys
import os
import commands
import string
import optparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import ROOT as rt
from Stat.Limits.settings import *

# open root file from Kevin's EOS area
# print limits
# exit
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--file", dest="fName", type=str, required=True, help="fileName")
parser.add_argument("-z", "--mZprime", dest="mZprime", type=int, required=True, help="mZprime")
parser.add_argument("-x", "--expr", dest="expr", type=str, default="M", help="mathematical expression to compute (M = median expected limit)")
parser.add_argument("-m", "--max", dest="max", type=float, default=10., help="max value to return")
args = parser.parse_args()

_file = rt.TFile.Open("root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Limits/datacards_07tsb_sys/{}".format(args.fName),"read")

limitTree = _file.Get("limit")

expr = args.expr.replace("M","{}")

for iEvt in range(limitTree.GetEntries()):
	limitTree.GetEvent(iEvt)
	if limitTree.quantileExpected == 0.5 and limitTree.trackedParam_mZprime == args.mZprime:
		result = eval(expr.format(limitTree.limit))
		print min(result,args.max) if args.max>0 else result
		break
