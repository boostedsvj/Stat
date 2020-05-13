#!/bin/bash

extra=$1
year=2018

EXTRA="--robustFit=1 --setRobustFitTolerance=1."
ENAME=""
if [ -n "$extra" ]; then
  EXTRA="--robustFit=1 --minos none --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex"
  ENAME="svjOpt_"
fi

for region in highCut lowCut highSVJ2 lowSVJ2; do
if [[ "$region" == "highCut" ]]; then
  PARS="--setParameters pdf_index_${region}_${year}=0,${region}_p1_2_alt=0,${region}_p2_2_alt=0 --freezeParameters pdf_index_${region}_${year},${region}_p1_2_alt,${region}_p2_2_alt --trackParameters ${region}_p1_3,${region}_p2_3,${region}_p3_3"
elif [[ "$region" == "lowCut" ]]; then
  PARS="--setParameters pdf_index_${region}_${year}=0,${region}_p1_2_alt=0,${region}_p2_2_alt=0 --freezeParameters pdf_index_${region}_${year},${region}_p1_2_alt,${region}_p2_2_alt --trackParameters ${region}_p1_2,${region}_p2_2,"
elif [[ "$region" == "highSVJ2" ]]; then
  PARS="--setParameters pdf_index_${region}_${year}=0,${region}_p1_1_alt=0 --freezeParameters pdf_index_${region}_${year},${region}_p1_1_alt --trackParameters ${region}_p1_1"
elif [[ "$region" == "lowSVJ2" ]]; then
  PARS="--setParameters pdf_index_${region}_${year}=0,${region}_p1_2_alt=0,${region}_p2_2_alt=0 --freezeParameters pdf_index_${region}_${year},${region}_p1_2_alt,${region}_p2_2_alt --trackParameters ${region}_p1_2,${region}_p2_2,"
else
  exit 1
fi

	echo "region: $region"
	combine -M FitDiagnostics SVJ_mZprime3000_mDark20_rinv03_alphapeak_${region}_2018_bias.txt $EXTRA -t 300 --toysFrequentist --saveToys --expectSignal 0.0 --rMin -80 --rMax 80 -n bias_${ENAME}${region} $PARS -v 3
done
