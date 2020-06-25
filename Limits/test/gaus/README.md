# Bias studies

Folder : bias study description:
* [datacards_example](./datacards_example): inject Gaussian signal w/ varying width (dijet search)
* [dijet](./dijet): multiply dijet signal rate by 100
* [svj](./svj): inject Gaussian signal w/ varying width (svj search)
* [dijet-svj](./dijet-svj): Rescale dijet data to SVJ-like yield, vary signal width and signal/data ratio
* [svj3000](./svj3000): svj search with MT < 3000
* [svj-dijet](./svj-dijet): svj search w/ dijet bins
* [svj-dijet2](./svj-dijet2): svj search w/ dijet bins and RooParametricShapeBinPdf
* [svj-para](./svj-para): svj search with RooParametricShapeBinPdf
* [dijet-svj2000](./dijet-svj2000): dijet search, data rescaled to SVJ-like yield, 2 TeV SVJ signal

Input workspaces in ROOT files are available in the corresponding folders at:
```
root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/gaus/
```

Instructions:
1. Run either `addSignals.py` or `transformSVJ.py` (sometimes the region might need to be provided: dijet, highCut, highSVJ2, lowCut, or lowSVJ2)
2. Run `runBias.sh` (possible arguments: region, "extra" = SVJ options rather than dijet options)
3. Run `plotBias.py` (again providing region, "extra"), `plotSignals.py`
