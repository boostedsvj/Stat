#!/bin/bash

#give names to parameters
mZ=${1}
mD=${2}
rI=${3}
aD=${4}
REGION=${5}
nTOYS=${6}
expSig=${7}
genFunc=${8} # 0 for bkgFunc, 1 for altFunc
fitFunc=${9} # 0 for bkgFunc, 1 for altFunc
cores=${10}
eosArea=${11}
rVal=${12}

SVJ_NAME="SVJ_mZprime${mZ}_mDark${mD}_rinv${rI}_alpha${aD}"
DC_NAME="${SVJ_NAME}_${REGION}_2018_template_bias.txt"
WS_NAME="ws_${SVJ_NAME}_${REGION}_2018_template.root"

if [ -n "$eosArea" ]; then
	export EOSDIR=${eosArea}/${SVJ_NAME}
	xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${DC_NAME} .
	xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${WS_NAME} .
fi

if [ -z "$rVal" ]; then
#	rVal=10
	rVal=50000
fi

rMin=-${rVal}
rMax=${rVal}

#parOptsMain is for using the Main function, we want to ignore the Alt parameters
#parOptsAlt is the opposite
#parOptsTrack<X> is used in the FitDiag Command to track the parameters we use
#  i.e., if we use parOptsMain, we should also use parOptsTrackMain, but only in FitDiag

eval $($CMSSW_BASE/src/Stat/Limits/test/getBiasArgs.py -r ${REGION} -n ${genFunc} -s Gen -f ${WS_NAME})
eval $($CMSSW_BASE/src/Stat/Limits/test/getBiasArgs.py -r ${REGION} -n ${fitFunc} -s Fit -f ${WS_NAME})

#print obtained variables
echo "SetArgGen: $SetArgGen"
echo "FrzArgGen: $FrzArgGen"
echo "TrkArgGen: $TrkArgGen"
echo "SetArgFit: $SetArgFit"
echo "FrzArgFit: $FrzArgFit"
echo "TrkArgFit: $TrkArgFit"

if [ ${genFunc} -eq 0 ]
then
  genName="${REGION}${optName}Sig${expSig}GenMain"
  genPdfName="Bkg_${REGION}_2018"
elif [ ${genFunc} -eq 1 ]
then
  genName="${REGION}${optName}Sig${expSig}GenAlt"
  genPdfName="Bkg_Alt_${REGION}_2018"
fi

if [ ${fitFunc} -eq 0 ]
then
  fitName="${genName}FitMain"
  pdfName="Bkg_${REGION}_2018"
elif [ ${fitFunc} -eq 1 ]
then
  fitName="${genName}FitAlt"
  pdfName="Bkg_Alt_${REGION}_2018"
fi

declare -A InitVals
InitVals["highCut"]="highCut_p1_3_alt=26.6638050328,highCut_p2_3_alt=-23.6201008859,highCut_p3_3_alt=0.146132892342"
InitVals["lowCut"]="lowCut_p1_2_alt=-9.91152443735,lowCut_p2_2_alt=-5.75832620307"
InitVals["highSVJ2"]="highSVJ2_p1_2_alt=-3.37387807418,highSVJ2_p2_2_alt=-6.59422151231"
InitVals["lowSVJ2"]="lowSVJ2_p1_2_alt=-11.1315816159,lowSVJ2_p2_2_alt=-6.43882169901"
if [ "$genFunc" -ne "$fitFunc" ]; then
	echo "Using precomputed brute force initial values"
	SetArgFit="${SetArgFit},${InitVals[$REGION]}"
	echo "SetArgFit: $SetArgFit"
fi
if false; then
	echo "Running brute force to get initial parameter values"
	cmdBF="python $CMSSW_BASE/src/Stat/Limits/python/bruteForce.py -f ${WS_NAME} -p ${pdfName} -n ${cores} -g ${genPdfName}"
	echo ${cmdBF}
	if [ -z "$DRYRUN" ]; then
		${cmdBF} >& bf.log
		cat bf.log
		InitVals=$(grep "setParameters" bf.log | tail -n 1 | cut -d' ' -f2)
		SetArgFit="${SetArgFit},${InitVals}"
		echo "SetArgFit: $SetArgFit"
	fi
fi

#cmdMDFa="combine ${DC_NAME} -M MultiDimFit -n ${genName} --saveWorkspace --rMin -80 --rMax 80 --setParameters $SetArgGen --freezeParameters $FrzArgGen"
cmdGen="combine ${DC_NAME} -M GenerateOnly -n ${genName} -t ${nTOYS} --toysFrequentist --saveToys --expectSignal ${expSig} --bypassFrequentistFit --saveWorkspace -s 123456 -v -1 --setParameters $SetArgGen --freezeParameters $FrzArgGen"

#cmdMDFb="combine ${DC_NAME} -M MultiDimFit -n ${fitName} --saveWorkspace --rMin -80 --rMax 80 --X-rtd MINIMIZER_MaxCalls=100000 --setParameters $SetArgFit --freezeParameters $FrzArgFit"
cmdFit="combine ${DC_NAME} -M FitDiagnostics -v -1 -n ${fitName} --toysFile higgsCombine${genName}.GenerateOnly.mH120.123456.root -t ${nTOYS} --toysFrequentist --saveToys --expectSignal ${expSig} --rMin ${rMin} --rMax ${rMax} --savePredictionsPerToy --bypassFrequentistFit --X-rtd MINIMIZER_MaxCalls=100000 --setParameters $SetArgFit --freezeParameters $FrzArgFit --trackParameters $TrkArgFit"

echo "combine commands:"
#echo ${cmdMDFa} # | tee /dev/stderr
#if [ -z "$DRYRUN" ]; then ${cmdMDFa}; fi
echo ${cmdGen} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdGen}; fi
#echo ${cmdMDFb} # | tee /dev/stderr
#if [ -z "$DRYRUN" ]; then ${cmdMDFb}; fi
echo ${cmdFit} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdFit}; fi
