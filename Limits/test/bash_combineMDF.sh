#!/bin/bash
# to be run : ./bash_combineMDF.sh <mZprime>
# script to run MultiDimFit on four regions for Likelihood scan 


# sample command for highCut 3000
#combine -M MultiDimFit SVJ_mZprime3000_mDark20_rinv03_alphapeak_highCut_2018_template_bias.txt --cminDefaultMinimizerType="Minuit" --freezeParameters highCut_p1_3_alt,highCut_p2_3_alt,highCut_p3_3_alt,pdf_index_highCut_2018 --setParameters pdf_index_highCut_2018=0 --algo grid --points 200 --setParameterRanges r=-10,10 --autoRange 2 --squareDistPoiStep

if [ -z "$1" ]; then
	echo "must provide mass"
	exit 1
fi

testDIR="/uscms_data/d3/cfallon/SVJ/biasStudies/CMSSW_10_2_13/src/Stat/Limits/test/"
testDIR=`pwd`/
wsDIR="cards_Mar30_modified/SVJ_mZprime${1}_mDark20_rinv03_alphapeak/"
wsPATH="${testDIR}${wsDIR}"
echo $wsPATH
for region in lowCut lowSVJ2
do
	cd ${wsPATH}${region}
	echo "REGION: $region"
	combine -M MultiDimFit ${wsPATH}SVJ_mZprime${1}_mDark20_rinv03_alphapeak_${region}_2018_template_bias.txt -n ${region} --cminDefaultMinimizerType="Minuit" --freezeParameters pdf_index_${region}_2018,${region}_p1_2_alt,${region}_p2_2_alt --setParameters r=0,pdf_index_${region}_2018=0,${region}_p1_2_alt=0,${region}_p2_2_alt=0 --algo grid --points 200 --setParameterRanges r=-3,3:${region}_p1_2=1,30:${region}_p2_2=1,16 --trackParameters ${region}_p1_2,${region}_p2_2
# --verbose 4
# --cminDefaultMinimizerAlgo Simplex --cminDefaultMinimizerStrategy 0 --robustFit 1 
# --skipInitialFit
done

for region in highSVJ2 highCut
do
	cd ${wsPATH}${region}
	echo "REGION: $region"
	combine -M MultiDimFit ${wsPATH}SVJ_mZprime${1}_mDark20_rinv03_alphapeak_${region}_2018_template_bias.txt -n ${region} --cminDefaultMinimizerType="Minuit" --freezeParameters pdf_index_${region}_2018,${region}_p1_1_alt --setParameters r=0,pdf_index_${region}_2018=0,${region}_p1_1_alt=0 --algo grid --points 200 --setParameterRanges r=-3,3:${region}_p1_1=1,30 --trackParameters ${region}_p1_1
# --verbose 4
# --cminDefaultMinimizerAlgo Simplex --cminDefaultMinimizerStrategy 0 --robustFit 1 
# --skipInitialFit
done
