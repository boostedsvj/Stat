import os,sys,shlex,traceback
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from multiprocessing import Pool
from StringIO import StringIO
import getBiasArgs
from Stat.Limits.bruteForce import makeVarInfoList
from paramUtils import fprint, runCmd, alphaVal, makeSigDict, getParamNames, getSigname, getSignameCheck, getFname, getWname, getDname, getDCname, getPname, getCombos, getInitFromBF, paramVal
import plotParamsScan
from makePostfitPlot import makePostfitPlot

def getXsecs(fname):
    with open(fname,'r') as xfile:
        xsecs = {xline.split('\t')[0]: float(xline.split('\t')[1]) for xline in xfile}
    return xsecs

def getFromToyfile(toyfile,match,replace=True,delim='.'):
    def match_text(item,match):
        return match in item
    def match_int(item,match):
        try:
            int(item)
            return True
        except:
            return False

    # special case
    if match=="#": matcher = match_int
    else: matcher = match_text
    result = next(item for item in toyfile.split(delim) if matcher(item,match))
    if replace and match!="#": result = result.replace(match,"")
    return result

def append_or_print(outputs,line,append=True):
    if append:
        outputs.append(line)
    else:
        fprint(line)
    return outputs

def try_plot_command(func, repr, pargs, outputs, append):
    # in order to restore stdout after capturing printouts from calls in plotting script/function
    backup_stdout = sys.stdout
    outputs = append_or_print(outputs, repr.format(','.join('"{}"'.format(p) if isinstance(p,str) else str(p) for p in pargs)), append)
    sys.stdout = StringIO()
    try:
	print('func: ', func)
        func(*pargs)
        this_output = sys.stdout.getvalue()
        sys.stdout = backup_stdout
        outputs = append_or_print(outputs, this_output, args.npool!=0)
    except Exception as e:
        sys.stdout = backup_stdout
        outputs = append_or_print(outputs, traceback.format_exc(), args.npool!=0)
    return outputs

def doLimit(info):
    args = info["args"]
    sig = info["sig"]
    signame = getSigname(sig)
    os.chdir(os.path.join(args.pwd,signame))

    params = {key:paramVal(key,val) for key,val in sig.iteritems()}
    setargs = []
    trkargs = []
    treargs = []
    extargs = ""
    for p,v in params.iteritems():
	print(p+"="+str(v))
        setargs.append(p+"="+str(v))
        trkargs.append(p)
        extargs = extargs+p+" extArg "+str(v)+"\n"
    frzargs = trkargs[:]

    fn = 1 if "Alt" in args.mod else 0
    for region in args.combo_regions:
        fname = getWname(signame,region)
        biasargs = getBiasArgs.main(fname, args.default_ws, region, fn, verbose=False)
        print(region, args.default_ws, fname)
        setargs.extend(biasargs['SetArg'])
        frzargs.extend(biasargs['FrzArg'])
        trkargs.extend(biasargs['TrkArg'])
        treargs.extend(biasargs['TrkArg'])
        if args.init:
            if isinstance(args.init,dict): setargs.extend(args.init[region].split(','))
            else: setargs.extend(getInitFromBF(fname, args.default_ws, getPname(region, "Alt" in args.mod)))

    #for ch in ["ch1","ch2"]:
    '''for ch in ["ch1"]:
        #normargs = ["n_exp_bin{}_proc_roomultipdf".format(ch),"shapeBkg_roomultipdf_{}__norm".format(ch),"n_exp_final_bin{}_proc_roomultipdf".format(ch),"n_exp_final_bin{}_proc_{}".format(ch,signame)]
        normargs = ["roomultipdf_norm"]
        trkargs.extend(normargs)
        treargs.extend(normargs)
        #if args.freezeNorm: frzargs.append("shapeBkg_roomultipdf_{}__norm".format(ch))
        if args.freezeNorm: frzargs.append("roomultipdf_norm")'''

    if "S0" in args.mod:
        frzargs.append("allConstrainedNuisances")
    if "Nostat" in args.mod:
        frzargs.append("rgx{mcstat_.*}")

    cargs = args.args
    if len(cargs)>0: cargs += " "
    cargs += "--setParameters {} --freezeParameters {} --trackParameters {} --trackErrors {} --keyword-value ana={} -n {}".format(
        ','.join(setargs), ','.join(frzargs), ','.join(trkargs), ','.join(treargs), args.combo, args.cname
    )
    if "Calls" in args.mod:
        cargs += " --X-rtd MINIMIZER_MaxCalls=100000"
    if len(args.toyfile)>0:
        cargs += " --toysFile {} -t {} --toysFrequentist".format(args.toyfile.format(args.combo),-1 if args.asimov else 1)
    datacards = []
    for region in args.combo_regions:
        datacards.append(getDname(signame,region).replace(".txt","{}.txt".format(args.suff)))
    dcfname = getDCname(sig,args.combo)

    outputs = []
    fprint("Calculating limit for {}...".format(signame))
    # combine cards
    command = "combineCards.py "+" ".join(datacards)+" > "+dcfname
    outputs = append_or_print(outputs, command, args.npool!=0)
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
            #extra += " -i shapeBkg high low"
	    extra += " -i shapeBkg bsvj"
        #command = 'python ../manualCLs.py {} -a="{}" -n {}'.format(extra,cargs,args.combo+"_"+args.cname)
        command = 'python ../manualCLs.py {} -p -f -s -a="{}" -n {}'.format(extra,cargs,args.combo+"_"+args.cname)
    else:
        command = "combine -M AsymptoticLimits "+cargs
    outputs = append_or_print(outputs, command, args.npool!=0)
    if not args.dry_run:
        # in this case, runCmd will print directly if not in pool mode, so just append to output (which will go unused)
        outputs = append_or_print(outputs, runCmd(command, args.npool==0))
        if args.plots:
            for step in ["Asimov","Observed"]:
                pargs = [sig,args.cname,step,args.combo,args.seedname,args.init]
                try_plot_command(plotParamsScan.main, "plotParamsScan.main({})", pargs, outputs, args.npool!=0)
            # currently, postfit files only created for ManualCLs
            if args.manualCLs:
                obs = len(args.toyfile)==0
                dfile = "hists.{}.{}.root".format(getFromToyfile(args.toyfile,"higgsCombine").replace("/",""), getFromToyfile(args.toyfile,"#")) if not obs else ""
                injected = int(getFromToyfile(args.toyfile,"mZprime",delim='_')) if "sigtoy" in args.toyfile else 0
                for q in [-3, -2, -1, 0.5]:
                    for region in args.combo_regions:
                        pargs = [sig,args.cname,"ManualCLsFit",q,dfile,args.datacards,obs,injected,args.combo,region,None,None,True]
                        try_plot_command(makePostfitPlot, "makePostfitPlot({})", pargs, outputs, args.npool!=0)
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
                for sig in args.signals:
                    # in non-pool mode, outputs are printed as they occur
                    doLimit({"args":args,"sig":sig})
            else:
                p = Pool(args.npool if not args.dry_run else 1)
                for outputs in p.imap_unordered(doLimit, [{"args":args,"sig":sig} for sig in args.signals]):
                    fprint('\n'.join(outputs))
                p.close()
                p.join()

        import ROOT as r
        if len(args.updateXsec)>0:
            xsecs = getXsecs(args.updateXsec[0])
            try:
                r.xsec_t
            except:
                r.gROOT.ProcessLine("struct xsec_t { Float_t mZprime; Float_t xsec; Double_t limit; };")
        outtrees = []
        outtreesroot = r.TList()
        for sig in args.signals:
            mname = "ManualCLs" if args.manualCLs and not "-A" in args.extra else "AsymptoticLimits"
            sname = "StepA" if args.manualCLs and "-A" in args.extra else ""
            fname = getFname(sname+cname,mname,combo,sig=sig,seed=args.seedname)
            if len(args.hadd_dir)>0: fname = args.hadd_dir+"/"+fname
            append = False
            if args.dry_run:
                append = True
                outtrees.append(fname)
            else:
                # check if limit converged
                f = r.TFile.Open(fname)
		print("filename: ", fname)
                if f!=None:
                    t = f.Get("limit")
                    if t!=None:
                        # 5 expected + 1 observed (+ prefit sometimes)
                        if t.GetEntries() >= 6:
                            append = True
                            t.SetDirectory(0)
                            # disable signal normalization branches: won't hadd properly
                            #t.SetBranchStatus("*{}".format(getSignameCheck(sig)),0)
			    #t.SetBranchStatus("{}".format(getSignameCheck(sig)),0)
                            t.SetBranchStatus("ana",0)
                            # can't update hadded tree because of ROOT bug, so update first
                            if len(args.updateXsec)>0:
                                xobj = r.xsec_t()
                                t.SetBranchAddress("trackedParam_mZprime",r.AddressOf(xobj,'mZprime'))
                                t.SetBranchAddress("trackedParam_xsec",r.AddressOf(xobj,'xsec'))
                                t.SetBranchAddress("limit",r.AddressOf(xobj,'limit'))
                                nt = t.CloneTree(0)
                                nt.SetDirectory(0)
                                for i in range(t.GetEntries()):
                                    t.GetEntry(i)
                                    xsec_orig = xobj.xsec
                                    xobj.xsec = xsecs[str(int(xobj.mZprime))]
                                    xobj.limit = xobj.limit*xsec_orig/xobj.xsec
                                    nt.Fill()
                                outtrees.append(nt)
                                outtreesroot.Add(nt)
                            else:
                                outtrees.append(t)
                                outtreesroot.Add(t)
                if not append: fprint("Warning: {} limit for {} did not converge".format(combo, getSigname(sig)))

        # combine trees
        if not args.no_hadd:
            if args.dry_run: fprint(outtrees)
            else:
                os.chdir(args.pwd)
                outname = "limit_"+combo+cname[4:]+(args.updateXsec[1] if len(args.updateXsec)>0 else "")+".root"
                outfile = r.TFile.Open(outname,"RECREATE")
                outtree = r.TTree.MergeTrees(outtreesroot)
                outtree.Write()
                outfile.Close()

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--regions", dest="regions", type=str, default=["cut","bdt"], nargs="+", help="list of combined regions")
    parser.add_argument("-n", "--npool", dest="npool", type=int, default=6, help="number of processes")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-f", "--freezeNorm", dest="freezeNorm", default=False, action="store_true", help="freeze bkg normalization to data")
    parser.add_argument("-m", "--mod", dest="mod", type=str, default=[], choices=["Alt","S0","Nostat","Calls"], nargs="*", help="modification(s)")
    parser.add_argument("-j", "--just-hadd", dest="just_hadd", default=False, action="store_true", help="don't run any combine commands, just hadd")
    parser.add_argument("--no-hadd", dest="no_hadd", default=False, action="store_true", help="don't hadd")
    parser.add_argument("--hadd-dir", dest="hadd_dir", type=str, default="", help="directory for files to be hadded if not local")
    #parser.add_argument("-M", "--manualCLs", dest="manualCLs", default=False, action='store_true', help="use manual CLs algorithm")
    parser.add_argument("-M", "--manualCLs", dest="manualCLs", default=True, action='store_true', help="use manual CLs algorithm") 
    init_group = parser.add_mutually_exclusive_group()
    init_group.add_argument("-i", "--initCLs", dest="initCLs", default=False, action='store_true', help="use initialized CLs algorithm")
    init_group.add_argument("-I", "--init", dest="init", default=False, action='store_true', help="use existing initial values of parameters")
    parser.add_argument("--newbf", dest="newbf", type=str, default="", help="file containing bf dict to import for initial values")
    parser.add_argument("--extra", dest="extra", type=str, default="", help="extra args for manual CLs")
    sig_group = parser.add_mutually_exclusive_group()
    sig_group.add_argument("--signal", dest="signals", metavar=tuple(getParamNames()), type=str, default=[], nargs=4, help="signal parameters")
    sig_group.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    parser.set_defaults(signals="default_signals.txt")
    parser.add_argument("-N", "--name", dest="name", type=str, default="Test", help="name for combine files")
    parser.add_argument("-s", "--suff", dest="suff", type=str, default="", help="suffix to pick different version of datacards")
    parser.add_argument("-t", "--toyfile", dest="toyfile", type=str, default="", help="toy file ({} in filename will be substituted with combined region)")
    parser.add_argument("--asimov", dest="asimov", default=False, action="store_true", help="toy file contains asimov dataset")
    parser.add_argument("-a", "--args", dest="args", type=str, default="", help="extra args for combine")
    parser.add_argument("-p", "--plots", dest="plots", default=False, action="store_true", help="make plots")
    #parser.add_argument("-p", "--plots", dest="plots", default=True, action="store_true", help="make plots")
    #parser.add_argument("--datacards", dest="datacards", type=str, default="root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig7/sigfull/", help="datacard histogram location (for postfit plots)")
    parser.add_argument("--datacards", dest="datacards", type=str, default="/home/saranabili/LimitStudy/September16_latestSetUp/CMSSW_10_2_13/src/Stat/Limits/test", help="datacard histogram location (for postfit plots)")
    parser.add_argument("-u", "--update-xsec", dest="updateXsec", type=str, metavar=('filename','suffix'), default=[], nargs=2, help="info to update cross sections when hadding")
    args = parser.parse_args()

    args.seedname = None
    if len(args.toyfile)>0:
        args.seedname = args.toyfile.split('.')[-2]
        args.args = " ".join([args.args,"-s {}".format(args.seedname),"--seedInName"])

    # parse signal info
    xsecs = getXsecs('dict_xsec_Zprime.txt')
    param_names = getParamNames()+["xsec"]
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
    args.signals = [makeSigDict(param_values[i],param_names) for i in range(len(param_values))]

    # pass some defaults

    args.default_ws = "SVJ"

    args.pwd = os.getcwd()

    main(args)

