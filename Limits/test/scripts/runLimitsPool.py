import os,sys,subprocess,shlex
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool

# make status messages useful
def fprint(msg):
    import sys
    print(msg)
    sys.stdout.flush()

def alpha_val(val):
    result = 0
    if val=="peak": result = -2
    elif val=="high": result = -1
    elif val=="low": result = -3
    else: result = float(val)
    return result

with open('dict_xsec_Zprime.txt','r') as xfile:
    xsecs = {int(xline.split('\t')[0]): float(xline.split('\t')[1]) for xline in xfile}

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-r", "--regions", dest="regions", type=str, default=["cut","bdt"], nargs="+", help="list of combined regions")
parser.add_argument("-n", "--npool", dest="npool", type=int, default=6, help="number of processes")
parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
parser.add_argument("-f", "--freezeNorm", dest="freezeNorm", default=False, action="store_true", help="freeze bkg normalization to data")
parser.add_argument("-m", "--mod", dest="mod", type=str, default=[], choices=["F12","Alt","Minfix","Robust","S0","Nostat"], nargs="*", help="modification(s)")
parser.add_argument("-j", "--just-hadd", dest="just_hadd", default=False, action="store_true", help="don't run any combine commands, just hadd")
parser.add_argument("-M", "--manualCLs", dest="manualCLs", default=False, action='store_true', help="use manual CLs algorithm")
parser.add_argument("--extra", dest="extra", type=str, default="", help="extra args for manual CLs")
args = parser.parse_args()

pwd = os.getcwd()

main_params = {
    "highCut": (4,5),
    "lowCut": (1,2),
    "highSVJ2": (1,2),
    "lowSVJ2": (1,2),
}
alt_params = {
    "highCut": (3,3),
    "lowCut": (3,3),
    "highSVJ2": (2,2),
    "lowSVJ2": (2,2),
}
def get_params(region,pdict,suff=""):
    order = pdict[region][0]
    n = pdict[region][1]
    params = ["{}_p{}_{}{}".format(region.replace("_2018",""),i+1,order,'_'+suff if len(suff)>0 else "") for i in range(n)]
    return params

def runCmd(args):
    output = ""
    try:
        output += subprocess.check_output(shlex.split(args),stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output += e.output
    return output

def doLimit(mass):
    signame = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(mass)
    params = {
        "xsec": xsecs[mass],
        "mZprime": mass,
        "mDark": 20,
        "rinv": 0.3,
        "alpha": alpha_val("peak"),
    }
    setargs = []
    trkargs = []
    extargs = ""
    for p,v in params.iteritems():
        setargs.append(p+"="+str(v))
        trkargs.append(p)
        extargs = extargs+p+" extArg "+str(v)+"\n"
    frzargs = trkargs[:]

    indarg = "pdf_index_{}_2018"
    indparams = {
        "cut": [indarg.format("highCut"),indarg.format("lowCut")],
        "bdt": [indarg.format("highSVJ2"),indarg.format("lowSVJ2")],
    }
    if "Alt" in args.mod:
        setargs.extend([p+"=1" for p in indparams[combo]])
    else:
        setargs.extend([p+"=0" for p in indparams[combo]])
    frzargs.extend(indparams[combo])

    fitparams = {
        "cut": get_params("highCut",main_params)+get_params("lowCut",main_params),
        "bdt": get_params("highSVJ2",main_params)+get_params("lowSVJ2",main_params),
    }
    altparams = {
        "cut": get_params("highCut",alt_params,"alt")+get_params("lowCut",alt_params,"alt"),
        "bdt": get_params("highSVJ2",alt_params,"alt")+get_params("lowSVJ2",alt_params,"alt"),
    }
    if "Alt" in args.mod:
        trkargs.extend(altparams[combo])
        frzargs.extend(fitparams[combo])
    else:
        trkargs.extend(fitparams[combo])
        frzargs.extend(altparams[combo])

    if "F12" in args.mod and combo=="cut":
        setargs.extend("highCut_p4_4=0,highCut_p5_4=0,lowCut_p3_4=0,lowCut_p4_4=0,lowCut_p5_4=0".split(','))
        frzargs.extend("highCut_p4_4,highCut_p5_4,lowCut_p3_4,lowCut_p4_4,lowCut_p5_4".split(','))

    for ch in ["ch1","ch2"]:
        trkargs.extend(["n_exp_bin{}_proc_roomultipdf".format(ch),"shapeBkg_roomultipdf_{}__norm".format(ch),"n_exp_final_bin{}_proc_roomultipdf".format(ch),"n_exp_final_bin{}_proc_SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(ch,mass)])
        if args.freezeNorm: frzargs.append("shapeBkg_roomultipdf_{}__norm".format(ch))

    if "S0" in args.mod:
        frzargs.append("allConstrainedNuisances")
    if "Nostat" in args.mod:
        frzargs.append("rgx{mcstat_.*}")

    os.chdir(os.path.join(pwd,signame))
    cargs = "--setParameters "+','.join(setargs)+" --freezeParameters "+','.join(frzargs)+" --trackParameters "+','.join(trkargs)+" --keyword-value ana="+combo+" -n "+cname
    if "Minfix" in args.mod:
        cargs += " --X-rtd improveFalseMinima"
    if "Robust" in args.mod:
        cargs += " --X-rtd allowRobustBisection1"
    datacards = []
    for region in regions:
        datacards.append(signame+"_{}_2018_template_bias_toy.txt".format(region))
    dcfname = "datacard_{}_{}.txt".format(mass,combo)

    outputs = []
    fprint("Calculating limit for {}...".format(mass))
    # combine cards
    command = "combineCards.py "+" ".join(datacards)+" > "+dcfname
    outputs.append(command)
    if not args.dry_run:
        os.system(command)
        # add ext args
        with open(dcfname,'a') as dcfile:
            dcfile.write(extargs)
    # run combine
    cargs += " -d "+dcfname
    if args.manualCLs:
        command = 'python ../manualCLs.py {} -a "{}" -n {}'.format(args.extra,cargs,combo+"_"+cname)
    else:
        command = "combine -M AsymptoticLimits "+cargs
    outputs.append(command)
    if not args.dry_run: 
        outputs.append(runCmd(command))
    os.chdir(pwd)

    return outputs

masses = [
1500,
1700,
1900,
2100,
2300,
2500,
2700,
2900,
3100,
3300,
3500,
3700,
3900,
4100,
4300,
4500,
4700,
4900,
5100,
]

cname = "Test"
if len(args.mod)>0: cname += ''.join(args.mod)
if args.freezeNorm: cname += "Frz"
if args.manualCLs: cname += "Manual"

combos = {
"cut": ["highCut","lowCut"],
"bdt": ["highSVJ2","lowSVJ2"],
}
for combo,regions in combos.iteritems():
    if not combo in args.regions: continue
    if not args.just_hadd:
        p = Pool(args.npool if not args.dry_run else 1)
        for outputs in p.imap_unordered(doLimit, masses):
            fprint('\n'.join(outputs))
        p.close()
        p.join()

    outfiles = []
    for mass in masses:
        signame = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(mass)
        mname = "ManualCLs" if args.manualCLs else "AsymptoticLimits"
        fname = signame+"/higgsCombine"+cname+"."+mname+".mH120.ana"+combo+".root"
        append = False
        if args.dry_run:
            append = True
        else:
            # check if limit converged
            from ROOT import TFile, TTree
            f = TFile.Open(fname)
            if f==None: continue
            t = f.Get("limit")
            if t==None: continue
            # 5 expected + 1 observed (+ prefit sometimes)
            append = t.GetEntries() >= 6
        if append: outfiles.append(fname)
        else: fprint("Warning: {} limit for mZprime = {} did not converge".format(combo, mass))

    # combine outfiles
    os.chdir(pwd)
    outname = "limit_"+combo+cname[4:]+".root"
    command = "hadd -f2 "+outname+''.join(" "+ofn for ofn in outfiles)
    fprint(command)
    if not args.dry_run: os.system(command)
