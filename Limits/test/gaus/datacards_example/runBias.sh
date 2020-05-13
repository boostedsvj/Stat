#!/bin/bash

extra=$1

EXTRA="--robustFit=1 --setRobustFitTolerance=1."
ENAME=""
if [ -n "$extra" ]; then
  EXTRA="--robustFit=1 --minos none --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex"
  ENAME="svjOpt_"
fi

for MU in 0 1; do
	for SIGMA in 0 1 2 3 4; do
		NAME=gaus_m${MU}_s${SIGMA}
		echo "Signal: $NAME"
		combine -M FitDiagnostics dijet_combine_${NAME}.txt $EXTRA -t 300 --toysFrequentist --saveToys --expectSignal 0.0 --rMin -80 --rMax 80 -n bias_$ENAME$NAME
	done
done
