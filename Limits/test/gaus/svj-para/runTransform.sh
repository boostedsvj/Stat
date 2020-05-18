#!/bin/bash

for region in lowCut highCut lowSVJ2 highSVJ2; do
	sed 's/roomultipdf/Bkg2_'$region'_2018/g; s/ws_/ws2_/g; s/pdf_index/#pdf_index/g;' SVJ_mZprime3000_mDark20_rinv03_alphapeak_${region}_2018_template_bias.txt > SVJ_mZprime3000_mDark20_rinv03_alphapeak_${region}_2018_template_bias2.txt
	python transformSVJ.py $region
done

