#!/bin/bash

for MU in 0 1; do
	for SIGMA in 0 1 2 3 4; do
		NAME=gaus_m${MU}_s${SIGMA}
		echo "Signal: $NAME"
		combine -M MultiDimFit dijet_combine_${NAME}.txt -n ${NAME} --algo grid --points 200 --setParameterRanges r=-3,3
	done
done
