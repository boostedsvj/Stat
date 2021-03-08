import os,sys,subprocess,shlex
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
import getBiasArgs
from Stat.Limits.bruteForce import makeVarInfoList
from collections import OrderedDict
from paramUtils import fprint, alphaVal, getSigname, getFname, getWname, getPname, getCombos, getInitFromBF

def runCmd(args):
    output = ""
    try:
        output += subprocess.check_output(shlex.split(args),stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output += e.output
    return output

def doLimit(info):
    args = info["args"]
    sig = info["sig"]
    signame = getSigname(sig)
    os.chdir(os.path.join(args.pwd,signame))

    params = {key:float(val) if key!="alpha" else alphaVal(val) for key,val in sig.iteritems()}
    setargs = []
    trkargs = []
    treargs = []
    extargs = ""
    for p,v in params.iteritems():
        setargs.append(p+"="+str(v))
        trkargs.append(p)
        extargs = extargs+p+" extArg "+str(v)+"\n"
    frzargs = trkargs[:]

    fn = 1 if "Alt" in args.mod else 0
    for region in args.combo_regions:
        fname = getWname(signame,region)
        biasargs = getBiasArgs.main(fname, args.default_ws, region, fn, verbose=False)
        setargs.extend(biasargs['SetArg'])
        frzargs.extend(biasargs['FrzArg'])
        trkargs.extend(biasargs['TrkArg'])
        treargs.extend(biasargs['TrkArg'])
        if args.init:
            if isinstance(args.init,dict): setargs.extend(args.init[region].split(','))
            else: setargs.extend(getInitFromBF(fname, args.default_ws, getPname(region, "Alt" in args.mod)))

    for ch in ["ch1","ch2"]:
        normargs = ["n_exp_bin{}_proc_roomultipdf".format(ch),"shapeBkg_roomultipdf_{}__norm".format(ch),"n_exp_final_bin{}_proc_roomultipdf".format(ch),"n_exp_final_bin{}_proc_{}".format(ch,signame)]
        trkargs.extend(normargs)
        treargs.extend(normargs)
        if args.freezeNorm: frzargs.append("shapeBkg_roomultipdf_{}__norm".format(ch))

    if "S0" in args.mod:
        frzargs.append("allConstrainedNuisances")
    if "Nostat" in args.mod:
        frzargs.append("rgx{mcstat_.*}")

    cargs = args.args
    cargs += " --setParameters {} --freezeParameters {} --trackParameters {} --trackErrors {} --keyword-value ana={} -n {}".format(
        ','.join(setargs), ','.join(frzargs), ','.join(trkargs), ','.join(treargs), args.combo, args.cname
    )
    if "Calls" in args.mod:
        cargs += " --X-rtd MINIMIZER_MaxCalls=100000"
    if len(args.toyfile)>0:
        cargs += " --toysFile {} -t {} --toysFrequentist".format(args.toyfile.format(args.combo),-1 if args.asimov else 1)
    datacards = []
    for region in args.combo_regions:
        datacards.append(signame+"_{}_2018_template_bias{}.txt".format(region,args.suff))
    dcfname = "datacard_{}_{}.txt".format('_'.join([val for key,val in sig.iteritems() if key!="xsec"]),args.combo)

    outputs = []
    fprint("Calculating limit for {}...".format(signame))
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
        extra = args.extra[:]
        # to use initCLs and manualCLs together, specify params to extract
        if args.initCLs:
            extra += " -i shapeBkg high low"
        command = 'python ../manualCLs.py {} -a "{}" -n {}'.format(extra,cargs,args.combo+"_"+args.cname)
    else:
        command = "combine -M AsymptoticLimits "+cargs
    outputs.append(command)
    if not args.dry_run: 
        outputs.append(runCmd(command))
    os.chdir(args.pwd)

    return outputs

def main(args):
    if args.newbf:
        sys.path.append(os.getcwd())
        args.init = getattr(__import__(args.newbf,fromlist=["bf"]),"bf")

    cname = args.name[:]
    if len(args.mod)>0: cname += ''.join(args.mod)
    if args.freezeNorm: cname += "Frz"
    if args.manualCLs and not "-A" in args.extra: cname += "Manual"
    if args.init: cname += "BFInit"
    if args.initCLs: cname += "Init"
    if "-b" in args.extra: cname += "Bonly"
    if "-s" in args.extra: cname += "Syst"
    args.cname = cname

    combos = getCombos()
    for combo,regions in combos.iteritems():
        if not combo in args.regions: continue
        args.combo = combo
        args.combo_regions = regions

        if not args.just_hadd:
            if args.npool==0:
                for outputs in [doLimit({"args":args,"sig":sig}) for sig in args.signals]:
                    fprint('\n'.join(outputs))
            else:
                p = Pool(args.npool if not args.dry_run else 1)
                for outputs in p.imap_unordered(doLimit, [{"args":args,"sig":sig} for sig in args.signals]):
                    fprint('\n'.join(outputs))
                p.close()
                p.join()

        outfiles = []
        for sig in args.signals:
            mname = "ManualCLs" if args.manualCLs and not "-A" in args.extra else "AsymptoticLimits"
            sname = "StepA" if args.manualCLs and "-A" in args.extra else ""
            seedname = None
            if len(args.toyfile)>0 and not args.asimov:
                seedname = args.toyfile.split('.')[-2]
            fname = getFname(sname+cname,mname,combo,sig=sig,seed=seedname)
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
            else: fprint("Warning: {} limit for {} did not converge".format(combo, signame))

        # combine outfiles
        if not args.no_hadd:
            os.chdir(args.pwd)
            outname = "limit_"+combo+cname[4:]+".root"
            command = "hadd -f2 "+outname+''.join(" "+ofn for ofn in outfiles)
            fprint(command)
            if not args.dry_run: os.system(command)

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--regions", dest="regions", type=str, default=["cut","bdt"], nargs="+", help="list of combined regions")
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=6, help="number of processes")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-f", "--freezeNorm", dest="freezeNorm", default=False, action="store_true", help="freeze bkg normalization to data")
    parser.add_argument("-m", "--mod", dest="mod", type=str, default=[], choices=["Alt","S0","Nostat","Calls"], nargs="*", help="modification(s)")
    parser.add_argument("-j", "--just-hadd", dest="just_hadd", default=False, action="store_true", help="don't run any combine commands, just hadd")
    parser.add_argument("--no-hadd", dest="no_hadd", default=False, action="store_true", help="don't hadd")
    parser.add_argument("-M", "--manualCLs", dest="manualCLs", default=False, action='store_true', help="use manual CLs algorithm")
    init_group = parser.add_mutually_exclusive_group()
    init_group.add_argument("-i", "--initCLs", dest="initCLs", default=False, action='store_true', help="use initialized CLs algorithm")
    init_group.add_argument("-I", "--init", dest="init", default=False, action='store_true', help="use existing initial values of parameters")
    parser.add_argument("--newbf", dest="newbf", type=str, default="", help="file containing bf dict to import for initial values")
    parser.add_argument("--extra", dest="extra", type=str, default="", help="extra args for manual CLs")
    sig_group = parser.add_mutually_exclusive_group()
    sig_group.add_argument("--signal", dest="signals", metavar=("mZprime","mDark","rinv","alpha"), type=str, default=[], nargs=4, help="signal parameters")
    sig_group.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    parser.set_defaults(signals="default_signals.txt")
    parser.add_argument("-N", "--name", dest="name", type=str, default="Test", help="name for combine files")
    parser.add_argument("-s", "--suff", dest="suff", type=str, default="", help="suffix to pick different version of datacards")
    parser.add_argument("-t", "--toyfile", dest="toyfile", type=str, default="", help="toy file ({} in filename will be substituted with combined region)")
    parser.add_argument("--asimov", dest="asimov", default=False, action="store_true", help="toy file contains asimov dataset")
    parser.add_argument("-a", "--args", dest="args", type=str, default="", help="extra args for combine")
    args = parser.parse_args()

    # parse signal info
    with open('dict_xsec_Zprime.txt','r') as xfile:
        xsecs = {xline.split('\t')[0]: float(xline.split('\t')[1]) for xline in xfile}
    param_names = ["mZprime", "mDark", "rinv", "alpha", "xsec"]
    param_values = []
    if isinstance(args.signals,list):
        param_values.append(args.signals)
        param_values[-1].append(xsecs[param_values[-1][0]])
    else:
        with open(args.signals,'r') as sfile:
            for line in sfile:
                line = line.rstrip()
                if len(line)==0: continue
                param_values.append(line.split())
                param_values[-1].append(xsecs[param_values[-1][0]])
    args.signals = [OrderedDict([(param_names[j],param_values[i][j]) for j in range(len(param_values[i]))]) for i in range(len(param_values))]

    # pass some defaults

    args.default_ws = "SVJ"

    args.pwd = os.getcwd()

    main(args)

