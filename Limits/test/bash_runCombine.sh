#!/bin/bash
# this should be run as ./bash_runCombine.py x y from the test dir.
# where x == 1 will force GEN to run
# and y == 1 will force FIT to run
runGen=$1
runFit=$2

zMass=$3

if [ $runGen == 1 ]
then
	echo "Doing GEN"
fi
if [ $runFit == 1 ]
then
	echo "Doing FIT"
fi

# Mar3_ThryDij
# Region Thry/Main Dijet/Alt
# highCut 2 3
# highSVJ2 1 1
# lowCut 2 2
# lowSVJ2 2 2

testDIR="/uscms_data/d3/cfallon/SVJ/biasStudies/CMSSW_10_2_13/src/Stat/Limits/test/"
wsDIR="cards_Mar30_modified/SVJ_mZprime${3}_mDark20_rinv03_alphapeak/"
wsPATH="${testDIR}${wsDIR}"
#echo $wsPATH
for region in lowCut lowSVJ2
do
	regionPATH="${testDIR}${wsDIR}${region}"
	for year in 2018
	do
		echo "Performing bias test for region $region"
		echo $wsPATH
		echo $regionPATH
		if [ ! -d $regionPATH ]
		then
			echo "Destination directory for $region doesn't exits. Making..."
			mkdir $regionPATH
		fi
		cd $regionPATH
		if [ $runGen == 1 ]
		then
			echo "Running GEN for MainSig0"
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M GenerateOnly --setParameters pdf_index_${region}_${year}=0 --toysFrequentist -t 1000 --expectSignal 0 -n MainSig0 --saveToys -m 125 --freezeParameters pdf_index_${region}_${year},${region}_p2_2_alt,${region}_p1_2_alt,${region}_p2_2,${region}_p1_2 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --bypassFrequentistFit > /dev/null
		fi
		if [ $runGen == 1 ]
		then
			echo "Running GEN for AltSig0"
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M GenerateOnly --setParameters pdf_index_${region}_${year}=1 --toysFrequentist -t 1000 --expectSignal 0 -n AltSig0 --saveToys -m 125 --freezeParameters pdf_index_${region}_${year},${region}_p2_2_alt,${region}_p1_2_alt,${region}_p2_2,${region}_p1_2 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --bypassFrequentistFit > /dev/null
		fi
		if [ $runFit == 1 ]
		then
			echo "M0M"		
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M FitDiagnostics --setParameters pdf_index_${region}_${year}=0,${region}_${year}_p2_2_alt=0,${region}_${year}_p1_2_alt=0 --toysFile $regionPATH/higgsCombineMainSig0.GenerateOnly.mH125.123456.root -n MainSig0Main -t 1000  --rMin -20 --rMax 20 --minos none --freezeParameters pdf_index_${region}_${year},${region}_p1_2_alt,${region}_p2_2_alt --cminDefaultMinimizerStrategy 0 --robustFit 1 --verbose 4 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --savePredictionsPerToy --trackParameters ${region}_p1_2,${region}_p2_2 --saveNormalizations --bypassFrequentistFit > /dev/null
		fi
		if [ $runFit == 1 ]
		then
			echo "A0M"
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M FitDiagnostics --setParameters pdf_index_${region}_${year}=0 --toysFile $regionPATH/higgsCombineAltSig0.GenerateOnly.mH125.123456.root -n AltSig0Main -t 1000 --rMin -20 --rMax 20 --minos none --freezeParameters pdf_index_${region}_${year} --cminDefaultMinimizerStrategy 0 --robustFit 1 --verbose 4 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --savePredictionsPerToy --trackParameters ${region}_p1_2,${region}_p2_2 --saveNormalizations --bypassFrequentistFit > /dev/null
		fi
		cd ..
	done
done


for region in highSVJ2 highCut
do
	regionPATH="${testDIR}${wsDIR}${region}"
	for year in 2018
	do
		echo "Performing bias test for region $region"
		echo $wsPATH
		echo $regionPATH
		if [ ! -d $regionPATH ]
		then
			echo "Destination directory for $region doesn't exits. Making..."
			mkdir $regionPATH
		fi
		cd $regionPATH
		if [ $runGen == 1 ]
		then
			echo "Running GEN for MainSig0"
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M GenerateOnly --setParameters pdf_index_${region}_${year}=0 --toysFrequentist -t 1000 --expectSignal 0 -n MainSig0 --saveToys -m 125 --freezeParameters pdf_index_${region}_${year},${region}_p1_1_alt,${region}_p1_1 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --bypassFrequentistFit > /dev/null
		fi
		if [ $runGen == 1 ]
		then
			echo "Running GEN for AltSig0"
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M GenerateOnly --setParameters pdf_index_${region}_${year}=1 --toysFrequentist -t 1000 --expectSignal 0 -n AltSig0 --saveToys -m 125 --freezeParameters pdf_index_${region}_${year},${region}_p1_1_alt,${region}_p1_1 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --bypassFrequentistFit > /dev/null
		fi
		if [ $runFit == 1 ]
		then
			echo "M0M"		
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M FitDiagnostics --setParameters pdf_index_${region}_${year}=0,${region}_${year}_p1_1_alt=0 --toysFile $regionPATH/higgsCombineMainSig0.GenerateOnly.mH125.123456.root -n MainSig0Main -t 1000  --rMin -20 --rMax 20 --minos none --freezeParameters pdf_index_${region}_${year},${region}_p1_1_alt --cminDefaultMinimizerStrategy 0 --robustFit 1 --verbose 4 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --savePredictionsPerToy --trackParameters ${region}_p1_1 --saveNormalizations --bypassFrequentistFit > /dev/null
		fi
		if [ $runFit == 1 ]
		then
			echo "A0M"		
			combine $wsPATH/SVJ_mZprime${3}_mDark20_rinv03_alphapeak_${region}_${year}_template_bias.txt -M FitDiagnostics --setParameters pdf_index_${region}_${year}=0 --toysFile $regionPATH/higgsCombineAltSig0.GenerateOnly.mH125.123456.root -n AltSig0Main -t 1000  --rMin -20 --rMax 20 --minos none --freezeParameters pdf_index_${region}_${year},${region}_p1_1_alt --cminDefaultMinimizerStrategy 0 --robustFit 1 --verbose 4 --cminDefaultMinimizerType Minuit --cminDefaultMinimizerAlgo Simplex --savePredictionsPerToy --trackParameters ${region}_p1_1 --saveNormalizations --bypassFrequentistFit > /dev/null
		fi
		cd ..
	done
done

