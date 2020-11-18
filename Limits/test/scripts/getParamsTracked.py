import os,sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-m", "--mass", dest="mass", type=int, required=True, help="Zprime mass")
parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[name])")
parser.add_argument("-M", "--method", dest="method", type=str, default="AsymptoticLimits", help="method name (higgsCombineTest.[method])")
parser.add_argument("-q", "--quantile", dest="quantile", type=float, default=-1, help="quantile to obtain tracked params")
args = parser.parse_args()

import ROOT as r

combos = {
    "cut": ["highCut","lowCut"],
    "bdt": ["highSVJ2","lowSVJ2"],
}
main_params = {
    "highCut": (4,5),
    "lowCut": (1,2),
    "highSVJ2": (1,2),
    "lowSVJ2": (1,2),
}
if "Alt" in args.name:
    main_params = {
        "highCut": (3,3),
        "lowCut": (3,3),
        "highSVJ2": (2,2),
        "lowSVJ2": (2,2),
    }

# todo: list all tree branches to find # params automatically
def get_params(region):
    order = main_params[region][0]
    n = main_params[region][1]
    params = ["{}_p{}_{}".format(region.replace("_2018",""),i+1,order)+("_alt" if "Alt" in args.name else "") for i in range(n)]
    return params

condition = "(quantileExpected-{})<0.001".format(args.quantile)
signame = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(args.mass)
print signame
for combo,regions in combos.iteritems():
    fname = "higgsCombine{}.{}.mH120.ana{}.root".format(args.name,args.method,combo)
    if not os.path.basename(os.getcwd())==signame: fname = signame+"/"+fname
    if not os.path.exists(fname): continue
    file = r.TFile.Open(fname)
    tree = file.Get("limit")
    for region in regions:
        p = get_params(region)
        p = ['trackedParam_'+pp for pp in p]
        n = tree.Draw(':'.join(p),condition,"goff")
        if n<=0: continue
        v = [str(tree.GetVal(i)[0]) for i in range(len(p))]

        ch = "ch1" if "high" in region else "ch2"
        p2 = ['trackedParam_shapeBkg_roomultipdf_'+ch+'__norm','trackedParam_n_exp_final_bin'+ch+'_proc_roomultipdf','trackedParam_n_exp_final_bin'+ch+'_proc_'+signame]
        n2 = tree.Draw(':'.join(p2),condition,"goff")
        v2 = [str(tree.GetVal(i)[0]) for i in range(len(p2))]

		# todo: print in useful format for postfit plots
        print '{}: {}\n  {}'.format(region, ','.join(v), ','.join(v2))

    n = tree.Draw("limit",condition,"goff")
    if n<=0: continue
    print 'r: {}'.format(tree.GetVal(0)[0])
