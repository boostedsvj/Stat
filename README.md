# Stat

## Updated method

Run [setup.sh](./setup.sh) to install the correct branches of all dependencies.

Fisher testing and datacard writing are split into three jobs.

First, do all the fits:
```
python createFits.py
```

Second, to do the F-tests and generate the workspaces:
```
python createFtest.py
```

Third, create the datacards:
```
python createDatacardsOnly.py
```

Fourth, run a bias test:
```
./combine_FourStepBiasBF.sh 2100 20 03 peak highCut 300 0 0 1 8
```

for Condor submission:
make sure to change directories in the following files:
* Stat/Limits/test/condorScripts/scramTarEos.sh : lines 19, 27
* Stat/Limits/test/condorScripts/ftest.sh : lines 7, 37
* Stat/Limits/test/condorScripts/ftest.jdl : line 10 
* Stat/Limits/test/condorScripts/allFits.sh : lines 7, 37
* Stat/Limits/test/condorScripts/allFits.jdl : line 11 
* Stat/Limits/test/condorScripts/datacardsOnly.sh : lines 7, 48 (the jdl doesn't require any change)
* Stat/Limits/test/condorScripts/condor_FourStepBiasBF.jdl : line 21
* Stat/Limits/test/condorScripts/condor_LimitBias.jdl : line 22
* Stat/Limits/test/condorScripts/condor_LimitBias.sh : line 7

step 1, do all fits:
```
condor_submit allFits.jdl
```

step 2, do F-tests:
```
condor_submit ftest.jdl
```

step 3, create datacards:
```
condor_submit datacardsOnly.jdl
```

step 4, run the combine commands to do the bias testing:
```
condor_submit condor_FourStepBiasBF.jdl
```

step 5, run the limit bias tests:
```
condor_submit condor_LimitBias.jdl
```
