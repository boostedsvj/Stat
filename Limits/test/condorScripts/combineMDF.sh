#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh

xrdcp -s root://cmseos.fnal.gov//store/user/cfallon/CMSSW_10_2_13.tgz .
tar -xf CMSSW_10_2_13.tgz
rm CMSSW_10_2_13.tgz
export SCRAM_ARCH=slc6_amd64_gcc700
cd CMSSW_10_2_13/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

echo "CMSSW: "$CMSSW_BASE

cd Stat/Limits/test

#give names to parameters
mZ=${1}
mD=${2}
rI=${3}
aD=${4}
REGION=${5}

SVJ_NAME="SVJ_mZprime${mZ}_mDark${mD}_rinv${rI}_alpha${aD}"

xrdcp -s root://cmseos.fnal.gov//store/user/cfallon/biasStudies/${SVJ_NAME}/${SVJ_NAME}_${REGION}_2018_template_bias.txt ${SVJ_NAME}_${REGION}_2018_template_bias.txt
xrdcp -s root://cmseos.fnal.gov//store/user/cfallon/biasStudies/${SVJ_NAME}/ws_${SVJ_NAME}_${REGION}_2018_template.root ws_${SVJ_NAME}_${REGION}_2018_template.root

if [ ${REGION} = "highCut" ]
then
  frzPar="pdf_index_${REGION}_2018,${REGION}_p1_2_alt,${REGION}_p2_2_alt"
  setPar="r=0,pdf_index_${region}_2018=0,${REGION}_p1_2_alt=0,${REGION}_p2_2_alt=0"
  rngPar="r=-3,3:${REGION}_p1_3=1,30:${REGION}_p2_3=-30,-1:${REGION}_p3_3=1,10"
  trkPar="${REGION}_p1_3,${REGION}_p2_3,${REGION}_p3_3"
elif [ ${REGION} = "highSVJ2" ]
then
  frzPar="pdf_index_${REGION}_2018,${REGION}_p1_1_alt"
  setPar="r=0,pdf_index_${region}_2018=0,${REGION}_p1_1_alt=0"
  rngPar="r=-3,3:${REGION}_p1_1=1,15"
  trkPar="${REGION}_p1_1"
elif [ ${REGION} = "lowCut" ]
then
  frzPar="pdf_index_${REGION}_2018,${REGION}_p1_2_alt,${REGION}_p2_2_alt"
  setPar="r=0,pdf_index_${region}_2018=0,${REGION}_p1_2_alt=0,${REGION}_p2_2_alt=0"
  rngPar="r=-3,3:${REGION}_p1_2=1,15:${REGION}_p2_2=1,15"
  trkPar="${REGION}_p1_2,${REGION}_p2_2"
elif [ ${REGION} = "lowSVJ2" ]
then
  frzPar="pdf_index_${REGION}_2018,${REGION}_p1_2_alt,${REGION}_p2_2_alt"
  setPar="r=0,pdf_index_${region}_2018=0,${REGION}_p1_2_alt=0,${REGION}_p2_2_alt=0"
  rngPar="r=-3,3:${REGION}_p1_2=1,15:${REGION}_p2_2=1,15"
  trkPar="${REGION}_p1_2,${REGION}_p2_2"
fi



cmd="combine -M MultiDimFit ${SVJ_NAME}_${REGION}_2018_template_bias.txt -n ${REGION} --cminDefaultMinimizerType=Minuit --freezeParameters ${frzPar} --setParameters ${setPar} --algo grid --points 200 --setParameterRanges ${rngPar} --trackParameters ${trkPar}"
echo "combine command:"
echo "${cmd}"
$cmd

# export items to EOS
echo "List all root files = "
ls *.root
echo "List all files"
ls 
echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/cfallon/biasStudies/${SVJ_NAME}
echo "xrdcp output for condor"
for FILE in *.root
do
  if [ ! "ws_${SVJ_NAME}_${REGION}_2018_template.root" = ${FILE} ]
  then
    echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
    xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  fi
  rm ${FILE}
done
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_13
