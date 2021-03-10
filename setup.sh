#!/bin/bash

scram project CMSSW_10_2_13
cd CMSSW_10_2_13/src
eval `scramv1 runtime -sh`
git clone github:kpedro88/HiggsAnalysis-CombinedLimit -b faster_para_plus HiggsAnalysis/CombinedLimit
git clone github:kpedro88/CombineHarvester -b fixes
git clone github:kpedro88/Analysis -b SVJ2018
git clone github:CTFallon/Stat
scram b -j 4
cd Analysis
./setupXML.sh
./setupScripts.sh
./recompile.sh

