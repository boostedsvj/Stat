import os,sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from Stat.Limits.bruteForce import makeVarInfoList

# make status messages useful
def fprint(msg):
    import sys
    print(msg)
    sys.stdout.flush()

def alphaVal(val):
    result = 0
    if val=="peak": result = -2
    elif val=="high": result = -1
    elif val=="low": result = -3
    else: result = float(val)
    return result

def getCombos():
    combos = {
        "cut": ["highCut","lowCut"],
        "bdt": ["highSVJ2","lowSVJ2"],
    }
    return combos

def getBranches(tree, matches=None, exact=False):
    import ROOT as r
    if not exact and not isinstance(matches,list): matches = [matches]
    elif exact and isinstance(matches,list): matches = matches[0]
    branches = []
    for b in tree.GetListOfBranches():
        bname = b.GetName()
        leaf = b.GetLeaf(bname)
        if matches is None or (exact and matches==b.GetName()) or (not exact and all(m in bname for m in matches)):
            branches.append(bname)
    return branches

def getParamNames():
    params = ["mZprime", "mDark", "rinv", "alpha"]
    return params

def makeSigDict(values,names=None):
    if names is None: names=getParamNames()
    sigdict = {names[j]:values[j] for j in range(len(values))}
    return sigdict

def getSigname(sig):
    params = getParamNames()
    return "SVJ_"+"_".join("{}{}".format(key,sig[key]) for key in params)

def getWname(sig,region):
    signame = sig if isinstance(sig,str) else getSigname(sig)
    wname = "ws_{}_{}_2018_template.root".format(signame,region)
    if not os.path.basename(os.getcwd())==signame: wname = signame+"/"+wname
    return wname

def getPname(region, alt=False):
    pname = "Bkg{}_{}_2018".format("_Alt" if alt else "", region)
    return pname

def getFname(name, method, combo, sig=None, prefix="higgsCombine", seed=None):
    if sig is not None:
        signame = getSigname(sig)
    fname = "{}{}.{}.mH120.ana{}{}.root".format(prefix,name,method,combo,"."+str(seed) if seed is not None else "")
    if sig is not None and not os.path.basename(os.getcwd())==signame: fname = signame+"/"+fname
    return fname

def getParamsText(params):
    return ["{}={}".format(p.replace('trackedParam_','').replace('trackedError_',''),v) for p,v in sorted(params.iteritems()) if any(x in p for x in ['high','low','shapeBkg'])]

def getParamsTracked(fname, quantile, includeParam=True, includeErr=False, extraCondition=""):
    import ROOT as r

    condition = "abs(quantileExpected-{})<0.001".format(quantile)
    if len(extraCondition)>0: condition += "&&{}".format(extraCondition)
    results = {}
    if not os.path.exists(fname): return results
    file = r.TFile.Open(fname)
    tree = file.Get("limit")

    matches = []
    if includeParam and not includeErr: matches = ["trackedParam"]
    elif not includeParam and includeErr: matches = ["trackedError"]
    params = []
    # background normalization factors
    params.extend(getBranches(tree, matches+["shapeBkg"]))
    # background & signal final normalizations
    params.extend(getBranches(tree, matches+["n_exp_final"]))
    # signal strength (from MDF)
    params.extend(getBranches(tree, "r", exact=True))

    for region in ["high","low"]:
        # background fit parameters
        params.extend(sorted(getBranches(tree, matches+[region])))

    # deliver dict of params:values
    n = tree.Draw(':'.join(params),condition,"goff")
    if n<=0: return results

    results = {p:tree.GetVal(i)[0] for i,p in enumerate(params)}
    return results

def getInitFromBF(fname, wsname, pdfname, region=None):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    file = r.TFile.Open(fname)
    ws = file.Get(wsname)
    pdf = ws.pdf(pdfname)
    pars = makeVarInfoList(pdf.getPars())
    npars = len(pars)

    if region is not None:
        fname2 = "fitResults_{}.root".format(region)
        if not os.path.isfile(fname2): fname2 = "../"+fname2
        file2 = r.TFile.Open(fname2)
        result = file2.Get("fitresult_{}{}_data_obs".format(pdfname,npars))

        setargs = {p.name:result.floatParsFinal().find(p.name).getValV() for p in pars}
        errargs = {p.name:result.floatParsFinal().find(p.name).getError() for p in pars}
        return setargs, errargs
    else:
        setargs = ["{}={}".format(p.name,p.val) for p in pars]
        return setargs

