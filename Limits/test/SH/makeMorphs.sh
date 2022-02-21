#!/bin/bash

MASS=$1
SVJNAME=SVJ_mZprime${MASS}_mDark20_rinv005_alphapeak
SVJNAME2=SVJ_${MASS}_20_0.05_peak
for REGION in highCut lowCut highSVJ2 lowSVJ2; do
	STATS=""
	for YEAR in 2016 2017 2018; do for BIN in {1..65}; do for VAR in Up Down; do STATS="${STATS} mcstat_${REGION}_2018_${SVJNAME}_MC${YEAR}bin${BIN}${VAR}"; done; done; done
	for UNC in "" MC2016JECUp MC2016JECDown MC2016JERUp MC2016JERDown MC2016puuncUp MC2016puuncDown MC2016trigfnuncUp MC2016trigfnuncDown MC2017JECUp MC2017JECDown MC2017JERUp MC2017JERDown MC2017puuncUp MC2017puuncDown MC2017trigfnuncUp MC2017trigfnuncDown MC2018JECUp MC2018JECDown MC2018JERUp MC2018JERDown MC2018puuncUp MC2018puuncDown MC2018trigfnuncUp MC2018trigfnuncDown MCRun2JESUp MCRun2JESDown MCRun2pdfalluncUp MCRun2pdfalluncDown MCRun2psfsruncUp MCRun2psfsruncDown MCRun2psisruncUp MCRun2psisruncDown $STATS; do
		python testMorph.py -s $MASS 20 005 peak -p rinv 0 01 -m integral -r $REGION -u "$UNC"
	done
	mkdir -p $SVJNAME
	hadd -f ${SVJNAME}/datacard_${SVJNAME2}_${REGION}.root ${SVJNAME}_*.root
	rm ${SVJNAME}_*.root
done
hadd -f ${SVJNAME}/datacard_final_${SVJNAME2}.root ${SVJNAME}/datacard_${SVJNAME2}_*.root root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig5/sigfull/datacard_bkg_data.root
rm ${SVJNAME}/datacard_${SVJNAME2}_*.root
python createDatacardsOnly.py -t -m template -Z ${MASS} -D 20 -R 005 -A peak -i ${SVJNAME}/
mv ws_${SVJNAME}_*.root ${SVJNAME}_*.txt $SVJNAME/
