import os,sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def getParamsTracked(mass, name, method, quantile, do_set, do_param, quiet):
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
    if "Alt" in name:
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
        params = ["{}_p{}_{}".format(region.replace("_2018",""),i+1,order)+("_alt" if "Alt" in name else "") for i in range(n)]
        return params

    condition = "abs(quantileExpected-{})<0.001".format(quantile)
    signame = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(mass)
    results = {}
    if not quiet: print signame
    for combo,regions in combos.iteritems():
        fname = "higgsCombine{}.{}.mH120.ana{}.root".format(name,method,combo)
        if mass!=0 and not os.path.basename(os.getcwd())==signame: fname = signame+"/"+fname
        if not os.path.exists(fname): continue
        results[combo] = []
        file = r.TFile.Open(fname)
        tree = file.Get("limit")
        for region in regions:
            # todo: refactor this part into separate function that just delivers dict of params:values
            p = get_params(region)
            p = ['trackedParam_'+pp for pp in p]
            n = tree.Draw(':'.join(p),condition,"goff")
            if n<=0: continue
            v = [str(tree.GetVal(i)[0]) for i in range(len(p))]

            ch = "ch1" if "high" in region else "ch2"
            p2 = ['trackedParam_shapeBkg_roomultipdf_'+ch+'__norm','trackedParam_n_exp_final_bin'+ch+'_proc_roomultipdf']
            if mass!=0: p2.append('trackedParam_n_exp_final_bin'+ch+'_proc_'+signame)
            n2 = tree.Draw(':'.join(p2),condition,"goff")
            v2 = [str(tree.GetVal(i)[0]) for i in range(len(p2))]

            # todo: print in useful format for postfit plots
            if not quiet: print '{}: {}\n  {}'.format(region, ','.join(v), ','.join(v2))

            if do_set:
                output = ["{}={}".format(p[i].replace('trackedParam_',''),v[i]) for i in range(len(p))]+["{}={}".format(p2[i].replace('trackedParam_',''),v2[i]) for i in range(len(p2)) if "n_exp" not in p2[i]]
                results[combo].extend(output)
                if not quiet: print ','.join(output)
            if do_param:
                for pp in [p,p2]:
                    vv = v if pp==p else v2
                    for i in range(len(pp)):
                        if "n_exp" in pp[i]: continue
                        # consider min and max from all expected limit values
                        n = tree.Draw(pp[i],"quantileExpected>0","goff")
                        v3 = tree.GetVal(0)
                        v3.SetSize(tree.GetSelectedRows())
                        stdev = max(abs(float(vv[i])-max(v3)),abs(float(vv[i])-min(v3)))
                        if not quiet: print '{} param {} {}'.format(pp[i].replace('trackedParam_',''),vv[i],stdev)

        n = tree.Draw("limit",condition,"goff")
        if n<=0: continue
        if not quiet: print 'r: {}'.format(tree.GetVal(0)[0])

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

    getParamsTracked(args.mass, args.name, args.method, args.quantile, args.set, args.param, args.quiet)
