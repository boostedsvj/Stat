import os,sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pprint import pprint

def getBranches(tree, match="", exact=False):
    import ROOT as r
    branches = []
    for b in tree.GetListOfBranches():
        leaf = b.GetLeaf(b.GetName())
        if len(match)==0 or (exact and match==b.GetName()) or (not exact and match in b.GetName()):
            branches.append(b.GetName())
    return branches

def getParamsAuto(tree, region):
    return sorted(getBranches(tree, region))

# retained for comparison w/ automatic approach above
def getParams(fn, region, prefix="trackedParams_"):
    params = dict(
        main = dict(
            highCut = (4,5),
            lowCut = (1,2),
            highSVJ2 = (1,2),
            lowSVJ2 = (1,2),
        ),
        alt = dict(
            highCut = (3,3),
            lowCut = (3,3),
            highSVJ2 = (2,2),
            lowSVJ2 = (2,2),
        )
    )

    order = params[fn][region][0]
    n = params[fn][region][1]
    pnames = ["{}{}_p{}_{}".format(prefix,region,i+1,order)+("_alt" if "Alt" in fn else "") for i in range(n)]
    return pnames

def getFname(mass, name, method, combo, quiet=True):
    signame = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(mass)
    if not quiet: print signame
    fname = "higgsCombine{}.{}.mH120.ana{}.root".format(name,method,combo)
    if mass!=0 and not os.path.basename(os.getcwd())==signame: fname = signame+"/"+fname
    return fname

def getParamsText(params):
    return ["{}={}".format(p.replace('trackedParam_',''),v) for p,v in sorted(params.iteritems()) if any(x in p for x in ['high','low','shapeBkg'])]

def getAll(mass, name, method, quantile, do_set, do_param, quiet):
    combos = {
        "cut": ["highCut","lowCut"],
        "bdt": ["highSVJ2","lowSVJ2"],
    }

    results = {}
    for combo,regions in combos.iteritems():
        results[combo] = getParamsTracked(getFname(mass, name, method, combo, quiet), quantile)
        if not quiet:
            pprint(results[combo])
        if do_set:
            output = getParamsText(results[combo])
            if not quiet: print ','.join(output)
        # todo: restore do_param case (consider min and max from all expected limit values, print in gaussian constrained param format)
    return results

def getParamsTracked(fname, quantile):
    import ROOT as r

    condition = "abs(quantileExpected-{})<0.001".format(quantile)
    results = {}
    if not os.path.exists(fname): return results
    file = r.TFile.Open(fname)
    tree = file.Get("limit")

    params = []
    # background normalization factors
    params.extend(getBranches(tree, "shapeBkg"))
    # background & signal final normalizations
    params.extend(getBranches(tree, "n_exp_final"))
    # signal strength (from MDF)
    params.extend(getBranches(tree, "r", exact=True))

    for region in ["high","low"]:
        # background fit parameters
        params.extend(getParamsAuto(tree, region))

    # deliver dict of params:values
    n = tree.Draw(':'.join(params),condition,"goff")
    if n<=0: return results

    results = {p:tree.GetVal(i)[0] for i,p in enumerate(params)}
    return results

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--mass", dest="mass", type=int, required=True, help="Zprime mass")
    parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[name])")
    parser.add_argument("-M", "--method", dest="method", type=str, default="AsymptoticLimits", help="method name (higgsCombineTest.[method])")
    parser.add_argument("-q", "--quantile", dest="quantile", type=float, default=-1, help="quantile to obtain tracked params")
    parser.add_argument("--set", dest="set", default=False, action="store_true", help="print in setParameters format")
    parser.add_argument("--param", dest="param", default=False, action="store_true", help="print in param format")
    parser.add_argument("--quiet", dest="quiet", default=False, action="store_true", help="suppress printouts")
    args = parser.parse_args()

    getAll(args.mass, args.name, args.method, args.quantile, args.set, args.param, args.quiet)

