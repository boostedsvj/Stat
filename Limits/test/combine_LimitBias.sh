#!/bin/bash -x

#give names to parameters
mZ=${1}
mD=${2}
rI=${3}
aD=${4}
COMBO=${5}
nTOYS=${6}
expSig=${7} # 0, 1, or M (unsure if decimals are OK for naming scheme)
genFunc=${8} # 0 for bkgFunc, 1 for altFunc
fitFunc=${9} # 0 for bkgFunc, 1 for altFunc
cores=${10}
eosArea=${11}

SVJ_NAME="SVJ_mZprime${mZ}_mDark${mD}_rinv${rI}_alpha${aD}"
if [ -n "$eosArea" ]; then
	export EOSDIR=${eosArea}/${SVJ_NAME}
fi

declare -A regions
regions[cut]="highCut lowCut"
regions[bdt]="highSVJ2 lowSVJ2"

DC_NAMES=""
WS_NAMES=""
SetArgGenAll=""
FrzArgGenAll=""
TrkArgGenAll=""
SetArgFitAll=""
FrzArgFitAll=""
TrkArgFitAll=""
for REGION in ${regions[$COMBO]}; do
	DC_NAME=${SVJ_NAME}_${REGION}_2018_template_bias.txt
	DC_NAMES="${DC_NAMES}${DC_NAMES:+ }${DC_NAME}"
	WS_NAME=ws_${SVJ_NAME}_${REGION}_2018_template.root
	WS_NAMES="${WS_NAMES}${WS_NAMES:+ }${WS_NAME}"
	if [ -n "$eosArea" ]; then
		xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${DC_NAME} .
		xrdcp -s root://cmseos.fnal.gov/${EOSDIR}/${WS_NAME} .
	fi
	eval $(./getBiasArgs.py -r ${REGION} -n ${genFunc} -s Gen -f ${WS_NAME})
	SetArgGenAll="${SetArgGenAll}${SetArgGenAll:+,}${SetArgGen}"
	FrzArgGenAll="${FrzArgGenAll}${FrzArgGenAll:+,}${FrzArgGen}"
	TrkArgGenAll="${TrkArgGenAll}${TrkArgGenAll:+,}${TrkArgGen}"
	eval $(./getBiasArgs.py -r ${REGION} -n ${fitFunc} -s Fit -f ${WS_NAME})
	SetArgFitAll="${SetArgFitAll}${SetArgFitAll:+,}${SetArgFit}"
	FrzArgFitAll="${FrzArgFitAll}${FrzArgFitAll:+,}${FrzArgFit}"
	TrkArgFitAll="${TrkArgFitAll}${TrkArgFitAll:+,}${TrkArgFit}"
done
DC_NAME_ALL=datacard_${mZ}_${COMBO}.txt
combineCards.py $DC_NAMES > $DC_NAME_ALL


#print obtained variables
echo "SetArgGenAll: $SetArgGenAll"
echo "FrzArgGenAll: $FrzArgGenAll"
echo "TrkArgGenAll: $TrkArgGenAll"
echo "SetArgFitAll: $SetArgFitAll"
echo "FrzArgFitAll: $FrzArgFitAll"
echo "TrkArgFitAll: $TrkArgFitAll"

if [ ${genFunc} -eq 0 ]
then
  genName="${COMBO}${optName}Sig${expSig}GenMain"
  genPdfName="Bkg_${REGION}_2018"
elif [ ${genFunc} -eq 1 ]
then
  genName="${COMBO}${optName}Sig${expSig}GenAlt"
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

# this is not fully updated for combos yet
if false && [ "$genFunc" -ne "$fitFunc" ]; then
	echo "Running brute force to get initial parameter values"
	cmdBF="python $CMSSW_BASE/src/Stat/Limits/python/bruteForce.py -f ${WS_NAME} -p ${pdfName} -n ${cores} -g ${genPdfName}"
	echo ${cmdBF}
	if [ -z "$DRYRUN" ]; then
		${cmdBF} >& bf.log
		cat bf.log
		InitVals=$(grep "setParameters" bf.log | tail -n 1 | cut -d' ' -f2)
		SetArgFitAll="${SetArgFitAll},${InitVals}"
		echo "SetArgFitAll: $SetArgFitAll"
	fi
fi

#setting expSig to median value if expSig equals M
maxValrinj=10
if [ ${expSig} = "M" ]; then
    cmdRE="python readREXT.py -f limit_${COMBO}AltManualBFInitSyst.root -z ${mZ} -m ${maxValrinj}"
    ${cmdRE} >& rext.log
    expSig=$(cat rext.log)
    echo $cmdRE
    echo "Median Extracted R is ${expSig} max ${maxValrinj}"
    #if [ ${expSig%.*} -gt 9 ]; then
    #    echo "Median Extracted R is ${expSig}. Setting to r_inj to 10."
    #    expSig=10
    #fi
fi

rMax=10

smallMax=$(bc <<< "${rMax}*1.0 < 3.0*${expSig}")
if [[ ${smallMax} > 0 ]]; then
    rMax=$(bc <<< "3.0*${expSig}")
fi

rMin=-${rMax}

cmdGen="combine ${DC_NAME_ALL} -M GenerateOnly -n ${genName} -t ${nTOYS} --toysFrequentist --saveToys --expectSignal ${expSig} --bypassFrequentistFit --saveWorkspace -s 123456 -v -1 --setParameters $SetArgGenAll --freezeParameters $FrzArgGenAll"

cmdFit="combine ${DC_NAME_ALL} -M FitDiagnostics -n ${fitName} --toysFile higgsCombine${genName}.GenerateOnly.mH120.123456.root -t ${nTOYS} -v -1 --toysFrequentist --saveToys --expectSignal ${expSig} --rMin ${rMin} --rMax ${rMax} --savePredictionsPerToy --bypassFrequentistFit --X-rtd MINIMIZER_MaxCalls=100000 --setParameters $SetArgFitAll --freezeParameters $FrzArgFitAll --trackParameters $TrkArgFitAll"

echo "combine commands:"
echo ${cmdGen} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdGen}; fi
echo ${cmdFit} # | tee /dev/stderr
if [ -z "$DRYRUN" ]; then ${cmdFit}; fi
