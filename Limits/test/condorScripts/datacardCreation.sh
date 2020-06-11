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

ls -la Stat/Limits/python
cd Stat/Limits/test
echo "Arguments passed to this script are:"
echo "Name of output directory : ${1}"
echo "Mode: ${2}"
echo "Signal Parameters: ${3} ${4} ${5} ${6}"
python createDatacardsCondor.py -d ${1} -m ${2} -t -Z ${3} -D ${4} -R ${5} -A ${6}

# export items to EOS
echo "List all root files = "
ls *.root
echo "List all files"
ls 
echo "*******************************************"
OUTDIR=root://cmseos.fnal.gov//store/user/cfallon/biasStudies/${1}
echo "xrdcp output for condor"
for FILE in *.root *.pdf *.txt
do
  echo "xrdcp -f ${FILE} ${OUTDIR}/${FILE}"
  xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
  rm ${FILE}
done
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_13
