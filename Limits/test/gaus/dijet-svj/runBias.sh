#!/bin/bash

region=$1
extra=$2

EXTRA="--robustFit=1 --setRobustFitTolerance=1."
ENAME=""
if [ -n "$extra" ]; then
  EXTRA="--robustFit=1 --minos none --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex"
  ENAME="svjOpt_"
fi

for SREGION in dijet highCut lowCut highSVJ2 lowSVJ2; do
	for SIGMA in 2.0 4.0; do
		NAME=${region}_gaus_r${SREGION}_s${SIGMA}
		echo "Signal: $NAME"
		combine -M FitDiagnostics dijet_combine_${NAME}.txt $EXTRA -t 300 --toysFrequentist --saveToys --expectSignal 0.0 --rMin -80 --rMax 80 -n bias_$ENAME$NAME
	done
done
