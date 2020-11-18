# Miscellaneous useful Combine scripts

## `printws.C(file,name)`
* ROOT macro to print the contents of RooFit workspace `name` from `file`
* useful for viewing large workspaces (redirect output to text file)

## `getParamsTracked.py`
* gets SVJ background fit parameters, signal & background normalizations, r values from Combine output tree `limit`
* can specify which quantile to use

## `runLimitsPool.py`
* make combined datacards for (highCut, lowCut) or (highSVJ2, lowSVJ2)
* make Combine commands (`AsymptoticLimits` or `manualCLs`, below) w/ all parameters set, tracked, frozen as needed
* run Combine commands locally using multiprocessing
* hadd output files from Combine commands

## `manualCLs.py`
* "manual" implementation of asymptotic CLs procedure using likelihood scans:
1. Estimate r range using AsymptoticLimits w/ systematics disabled
2. Run likelihood scans using MultiDimFit
3. Compute CLs values (can also make diagnostic plots)
4. Optional: do MultiDimFit for each quantile r value (for prefit/postfit plots)
5. Make new `limit` tree w/ results from step 3 or 4
* Combine output (steps 1,2,4) can be reused when rerunning command

Some features in these scripts require one of the following Combine branches:
1. Minor fixes: https://github.com/kpedro88/HiggsAnalysis-CombinedLimit/tree/faster_para_plus
2. above + ad-hoc procedures to try to improve false minima: https://github.com/kpedro88/HiggsAnalysis-CombinedLimit/tree/faster_para_plus_debug_improve
