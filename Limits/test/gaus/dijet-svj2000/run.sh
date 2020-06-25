#!/bin/bash

mkdir -p reuse
mkdir -p onTheFly
mkdir -p bypass

combine dijet_combine_highSVJ2_SVJ_mZprime2000_mDark20_rinv03_alphapeak.txt -M MultiDimFit --saveWorkspace  --rMin -80 --rMax 80

combine higgsCombineTest.MultiDimFit.mH120.root --snapshotName MultiDimFit -M FitDiagnostics --robustFit=1 --setRobustFitTolerance=1. -n testFit -t 1 --toysFrequentist --saveToys --expectSignal 1 --rMin -80 --rMax 80 --saveShapes --plots -s 1 --saveWorkspace >& log_onTheFly.log && mv *.png *Diag*.root onTheFly/


combine higgsCombineTest.MultiDimFit.mH120.root --snapshotName MultiDimFit -M FitDiagnostics --robustFit=1 --setRobustFitTolerance=1. -n testFit -t 1 --toysFrequentist --saveToys --expectSignal 1 --rMin -80 --rMax 80 --saveShapes --plots -s 1 --saveWorkspace --bypassFrequentistFit >& log_bypass.log && mv *.png *Diag*.root bypass/



combine dijet_combine_highSVJ2_SVJ_mZprime2000_mDark20_rinv03_alphapeak.txt -M MultiDimFit --saveWorkspace  --rMin -80 --rMax 80

combine higgsCombineTest.MultiDimFit.mH120.root --snapshotName MultiDimFit -M GenerateOnly -t 1 --saveToys --toysFrequentist --expectSignal 1 -s 1 

combine higgsCombineTest.MultiDimFit.mH120.root --snapshotName MultiDimFit -M FitDiagnostics --robustFit=1 --setRobustFitTolerance=1. -n testFit --toysFile higgsCombineTest.GenerateOnly.mH120.1.root -t 1 --toysFrequentist --saveToys --expectSignal 1 --rMin -80 --rMax 80 --saveShapes --plots -s 1 --saveWorkspace >& log_reuse.log && mv *.png *Diag*.root reuse/

