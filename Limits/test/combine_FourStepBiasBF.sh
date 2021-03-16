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

SVJ_NAME="SVJ_mZprime${mZ}_mDark${mD}_rinv${rI}_alpha${aD}"
DC_NAME="${SVJ_NAME}_${REGION}_2018_template_bias.txt"
WS_NAME="ws_${SVJ_NAME}_${REGION}_2018_template.root"

if [ -n "$eosArea" ]; then
	export EOSDIR=${eosArea}/${SVJ_NAME}
	xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${DC_NAME} .
	xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${WS_NAME} .
fi

rMin=-10
rMax=10

#parOptsMain is for using the Main function, we want to ignore the Alt parameters
#parOptsAlt is the opposite
#parOptsTrack<X> is used in the FitDiag Command to track the parameters we use
#  i.e., if we use parOptsMain, we should also use parOptsTrackMain, but only in FitDiag

eval $(./getBiasArgs.py -r ${REGION} -n ${genFunc} -s Gen -f ${WS_NAME})
eval $(./getBiasArgs.py -r ${REGION} -n ${fitFunc} -s Fit -f ${WS_NAME})

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

if [ "$genFunc" -ne "$fitFunc" ]; then
	echo "Running brute force to get initial parameter values"
	cmdBF="python $CMSSW_BASE/src/Stat/Limits/python/bruteForce.py -f ${WS_NAME} -p ${pdfName} -n ${cores} -g ${genPdfName}"
	echo ${cmdBF}
	if [ -z "$DRYRUN" ]; then
		${cmdBF} >& bf.log
		cat bf.log
		InitVals=$(grep "setParameters" bf.log | tail -n 1 | cut -d' ' -f2)
		SetArgFit="${SetArgFit},${InitVals}"
	fi
fi

cmdMDFa="combine ${DC_NAME} -M MultiDimFit -n ${genName} --saveWorkspace --rMin -80 --rMax 80 --setParameters $SetArgGen --freezeParameters $FrzArgGen"
cmdGen="combine higgsCombine${genName}.MultiDimFit.mH120.root --snapshotName MultiDimFit -M GenerateOnly -n ${genName} -t ${nTOYS} --toysFrequentist --saveToys --expectSignal ${expSig} --bypassFrequentistFit --saveWorkspace -s 123456 -v -1 --setParameters $SetArgGen --freezeParameters $FrzArgGen"

cmdMDFb="combine ${DC_NAME} -M MultiDimFit -n ${fitName} --saveWorkspace --rMin -80 --rMax 80 --X-rtd MINIMIZER_MaxCalls=100000 --setParameters $SetArgFit --freezeParameters $FrzArgFit"
cmdFit="combine higgsCombine${fitName}.MultiDimFit.mH120.root --snapshotName MultiDimFit -M FitDiagnostics -v -1 -n ${fitName} --toysFile higgsCombine${genName}.GenerateOnly.mH120.123456.root -t ${nTOYS} --toysFrequentist --saveToys --expectSignal ${expSig} --rMin ${rMin} --rMax ${rMax} --savePredictionsPerToy --bypassFrequentistFit --setParameters $SetArgFit --freezeParameters $FrzArgFit --trackParameters $TrkArgFit"

echo "combine commands:"
echo ${cmdMDFa} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdMDFa}; fi
echo ${cmdGen} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdGen}; fi
echo ${cmdMDFb} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdMDFb}; fi
echo ${cmdFit} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdFit}; fi
