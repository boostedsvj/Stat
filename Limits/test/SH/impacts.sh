#!/bin/bash

#DRYRUN=echo
DRYRUN=""

MASS=$1
if [ -z "$MASS" ]; then
	exit 1
fi
eosArea=$2
TOYARG=$3
TOYNAME=$4
COMBOIN=$5
if [ -z "$COMBOIN" ]; then
	COMBOS=(cut bdt)
else
	COMBOS=($COMBOIN)
fi

SVJ_NAME=SVJ_mZprime${MASS}_mDark20_rinv03_alphapeak
if [ -n "$eosArea" ]; then
	export EOSDIR=${eosArea}/${SVJ_NAME}
fi

declare -A regions
regions[cut]="highCut lowCut"
regions[bdt]="highSVJ2 lowSVJ2"

for COMBO in ${COMBOS[@]}; do
	DC_NAMES=""
	WS_NAMES=""
	SetArgAll=""
	FrzArgAll=""
	TrkArgAll=""
	for REGION in ${regions[$COMBO]}; do
		DC_NAME=${SVJ_NAME}_${REGION}_2018_template_bias.txt
		DC_NAMES="${DC_NAMES}${DC_NAMES:+ }${DC_NAME}"
		WS_NAME=ws_${SVJ_NAME}_${REGION}_2018_template.root
		WS_NAMES="${WS_NAMES}${WS_NAMES:+ }${WS_NAME}"
		if [ -n "$eosArea" ]; then
			xrdcp -f -s root://cmseos.fnal.gov/${EOSDIR}/${DC_NAME} .
			xrdcp -f -s root://cmseos.fnal.gov/${EOSDIR}/${WS_NAME} .
		fi
		eval $(./getBiasArgs.py -r ${REGION} -n 1 -s "" -f ${WS_NAME})
		SetArgAll="${SetArgAll}${SetArgAll:+,}${SetArg}"
		FrzArgAll="${FrzArgAll}${FrzArgAll:+,}${FrzArg}"
		TrkArgAll="${TrkArgAll}${TrkArgAll:+,}${TrkArg}"
	done
	DC_NAME_ALL=datacard_${MASS}_${COMBO}
	$DRYRUN combineCards.py $DC_NAMES > ${DC_NAME_ALL}.txt
	TESTNAME=Test${TOYNAME}${COMBO}
	ARGS="--rMin -10 -n ${TESTNAME} --setParameters $SetArgAll --freezeParameters $FrzArgAll --trackParameters $TrkArgAll"
	if [ -n "$TOYARG" ]; then
		ARGS="$ARGS $(echo "$TOYARG" | sed 's~COMBO~'${COMBO}'~')"
	fi

	OUTNAME=impacts_${COMBO}${TOYNAME:+_}${TOYNAME}
	(set -x
	$DRYRUN text2workspace.py ${DC_NAME_ALL}.txt
	$DRYRUN combineTool.py -M Impacts -d ${DC_NAME_ALL}.root --doInitialFit -m 125 $ARGS
	$DRYRUN combineTool.py -M Impacts -d ${DC_NAME_ALL}.root -m 125 --doFits --parallel 16 $ARGS
	$DRYRUN combineTool.py -M Impacts -d ${DC_NAME_ALL}.root -o ${OUTNAME}.json -m 125 -n ${TESTNAME}
	$DRYRUN plotImpacts.py -i ${OUTNAME}.json -o ${OUTNAME}
	# subset
	$DRYRUN python excludeImpacts.py -i ${OUTNAME}.json -m mcstat -x $(echo $FrzArgAll | tr ',' ' ')
	$DRYRUN plotImpacts.py -i ${OUTNAME}_include.json -o ${OUTNAME}_include
	$DRYRUN plotImpacts.py -i ${OUTNAME}_include.json -o ${OUTNAME}_include_multi --per-page 4 --height 300 --label-size -1 --width 400
	)
done
