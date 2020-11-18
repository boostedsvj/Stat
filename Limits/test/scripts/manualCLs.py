import os,sys,subprocess,shlex,uuid
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# make status messages useful
def fprint(msg):
    import sys
    print(msg)
    sys.stdout.flush()

# flags should be different representations of the same flag, e.g. -a/--args
def updateArg(args, flags, val, sep=""):
    # make a copy
    args = args[:]
    if not isinstance(flags,list): flags = [flags]
    for flag in flags:
        if flag+" " in args:
            args = args.replace(flag+" ", flag+" "+val+sep, 1)
            return args
    args += " "+flag+" "+val
    return args

def runCmd(cmd, log, opt='w'):
    with open(log,opt) as logfile:
        subprocess.call(shlex.split(cmd),stdout=logfile,stderr=logfile)

def getOutfile(log):
    ofname = ""
    indicator = "COMBINE_OUTPUT_FILE: "
    failure = "WARNING: MultiDimFit failed"
    success = True
    with open(log,'r') as logfile:
        for line in logfile:
            # stop if found both checks
            if len(ofname)>0 and not success: break
        
            if indicator in line:
                ofname = line.rstrip().replace(indicator,"")
            elif failure in line:
                success = False

    if len(ofname)==0: raise RuntimeError("Could not find output file name from log: {}".format(log))
    return ofname, success

def getRange(dry_run, ofname1, count_lower, count_upper):
    # get r range
    # (default vals provided for dryrun printouts)
    rmin = 0.
    rmax = 10.
    factor = 10.
    factor_min = factor*pow(2,count_lower)
    factor_max = factor*pow(2,count_upper)
    npts = 100
    # increase npts by 1 to include both endpoints
    npts += 1
    if not dry_run:
        import ROOT as r
        file1 = r.TFile.Open(ofname1)
        limit1 = file1.Get("limit")
        n = limit1.Draw("limit","abs(quantileExpected-0.975)<0.001||abs(quantileExpected-0.025)<0.001","goff")
        if n!=2: raise RuntimeError("Malformed limit tree in "+ofname1)
        vals = [limit1.GetVal(0)[0],limit1.GetVal(0)[1]]
        rmin = min(vals)/factor_min
        rmax = max(vals)*factor_max

    return rmin,rmax,npts

def step1(args):
    # run AsymptoticLimits w/ nuisances disabled
    args1 = updateArg(args.args, ["--freezeParameters"], "allConstrainedNuisances", ',')
    args1 = updateArg(args1, ["-n","--name"], "Step1")
    cmd1 = "combine -M AsymptoticLimits "+args1
    fprint(cmd1)
    logfname1 = "log_step1_{}.log".format(args.name)
    ofname1 = ""
    if not args.dry_run:
        if "step1" not in args.reuse: runCmd(cmd1,logfname1)
        ofname1, _ = getOutfile(logfname1)
    return ofname1

def step2impl(args, name, lname, rmin, rmax, npts, extra=""):
    # run MDF likelihood scan
    args2 = updateArg(args.args, ["-n","--name"], name)
    cmd2 = "combine -M MultiDimFit --redefineSignalPOIs r --setParameterRanges r={},{} --algo grid --points {} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --alignEdges 1 --saveNLL {} {} {}".format(rmin,rmax,npts,extra,args2,args.fitopts)
    fprint(cmd2)
    logfname2 = "log_{}_{}.log".format(lname,args.name)
    ofname2 = ""
    if not args.dry_run:
        if "step2" not in args.reuse: runCmd(cmd2,logfname2)
        ofname2, _ = getOutfile(logfname2)
    return ofname2

def step2(args, ofname1, count_lower, count_upper):
    # get rmin, rmax from step1
    rmin, rmax, npts = getRange(args.dry_run,ofname1,count_lower,count_upper)

    # observed
    ofname2d = step2impl(args,"Observed","step2d",rmin,rmax,npts)

    # expected (asimov)
    ofname2a = step2impl(args,"Asimov","step2a",rmin,rmax,npts,"-t -1 --toysFreq")
    
    return ofname2d, ofname2a

def step3(args, scans):
    # based on https://gitlab.cern.ch/cms-hcg/cadi/hig-19-003/-/blob/master/HJMINLO/plot_cls.py
    interpolate = True
    inter_npoints = 1000
    quantiles = [0.025, 0.16, 0.50, 0.84, 0.975]

    if args.dry_run:
        # return dummy dictionary to be used in step4 commands
        return {q:0.0 for q in quantiles+[-1]}

    import ROOT as r
    # put root in batch mode
    r.gROOT.SetBatch()
    file_data = r.TFile.Open(scans[0])
    file_asimov = r.TFile.Open(scans[1])
    tree_data = file_data.Get('limit')
    tree_asimov = file_asimov.Get('limit')

    # in case of failures, only keep r values that succeeded for both data and asimov
    r_dict = {}
    r_data_bestfit = -1
    r_asimov_bestfit = -1

    for e in tree_data:
        if e.quantileExpected == -1:
            r_data_bestfit = e.r
        elif e.quantileExpected > -1:
            r_dict[e.r] = [e.deltaNLL,None]

    for e in tree_asimov:
        if e.quantileExpected == -1:
            r_asimov_bestfit = e.r
        elif e.quantileExpected > -1:
            if e.r in r_dict: r_dict[e.r][1] = e.deltaNLL

    r_dict = {k:v for k,v in r_dict.iteritems() if v[1] is not None}
    r_data = list(sorted(r_dict))
    rmin = min(r_data)
    rmax = max(r_data)
    dnll_data = [r_dict[k][0] for k in r_data]
    dnll_asimov = [r_dict[k][1] for k in r_data]

    clsb = []
    clb = []
    cls = []
    cls_exp = {}
    for q in quantiles:
        cls_exp[q] = []

    # some plots are needed for interp, so just setup all
    cls_graph = r.TGraph(len(r_data))
    clsb_graph = r.TGraph(len(r_data))
    clb_graph = r.TGraph(len(r_data))
    cls_exp_graph = {}
    for q in quantiles:
        cls_exp_graph[q] = r.TGraph(len(r_data))
    dn2ll_data_graph = r.TGraph(len(r_data))
    dn2ll_asimov_graph = r.TGraph(len(r_data))

    for i in range(len(r_data)):
        rval = r_data[i]
        if rval < r_data_bestfit:
            dnll_data_min = min([dnll_data[j] for j in range(0,i+1)])
            dnll_data_constrained = dnll_data[i] - dnll_data_min
        else: 
            dnll_data_constrained = dnll_data[i]

        qmu = 2*max(dnll_data_constrained,0.)
        qA = 2*max(dnll_asimov[i],0.)
        s_exp = {}
        for q in quantiles:
            n = r.Math.normal_quantile(q, 1.0)
            s_exp[q] = r.Math.normal_cdf_c( r.TMath.Sqrt(qA) - n, 1.)/q

        if qmu >= 0 and qmu <= qA:
            sb = r.Math.normal_cdf_c( r.TMath.Sqrt(qmu) )
            b = r.Math.normal_cdf( r.TMath.Sqrt(qA) - r.TMath.Sqrt(qmu) )
        elif qmu > qA:
            sb = r.Math.normal_cdf_c( (qmu + qA) / (2*r.TMath.Sqrt(qmu)) )
            b = r.Math.normal_cdf_c( (qmu - qA) / (2*r.TMath.Sqrt(qmu)) )
        else:
            raise ValueError("Unexpected q values: q_mu = {}, q_A = {}".format(q_mu,q_A))

        s = sb/b
        clsb.append(sb)
        clb.append(b)
        cls.append(s)
        for q in quantiles:
            cls_exp[q].append(s_exp[q])

        # fill plots
        cls_graph.SetPoint(i, rval, s)
        clsb_graph.SetPoint(i, rval, sb)
        for q in quantiles:
            cls_exp_graph[q].SetPoint(i, rval, s_exp[q])
        clb_graph.SetPoint(i, rval, b)
        dn2ll_data_graph.SetPoint(i, rval, 2*dnll_data_constrained)
        dn2ll_asimov_graph.SetPoint(i, rval, 2*dnll_asimov[i])

    import numpy as np
    def find_crossing(array, value):
        array = np.asarray(array)
        array = array - value
        across = (array < 0)
        idx = across.argmax() if across.any() else -1
        # check if previous point is closer to zero
        if idx>0 and abs(array[idx-1]) < abs(array[idx]): idx = idx-1
        return idx

    if interpolate:
        cls = []
        r_data = []
        for q in quantiles:
            cls_exp[q] = []
        for rval in np.linspace(rmin, rmax, inter_npoints+1):
            r_data.append(rval)
            cls.append(cls_graph.Eval(rval,0,'S'))
            for q in quantiles: 
                cls_exp[q].append(cls_exp_graph[q].Eval(rval,0,'S'))

    limits = {}
    quantiles_at_lower_boundary = []
    quantiles_at_upper_boundary = []
    alpha = 0.05
    idx = find_crossing(cls,alpha)
    limits[-1] = r_data[idx]
    fprint("Observed Limit: r < {:.2f}".format(limits[-1]))
    if idx==0: quantiles_at_lower_boundary.append(-1)
    if idx==len(r_data)-1 or idx==-1: quantiles_at_upper_boundary.append(-1)
    for q in quantiles:
        qidx = find_crossing(cls_exp[q],alpha)
        limits[q] = r_data[qidx]
        fprint("Expected {:3.1f}%: r < {:.2f}".format(q*100., limits[q]))
        if qidx==0: quantiles_at_lower_boundary.append(q)
        if qidx==len(r_data)-1 or qidx==-1: quantiles_at_upper_boundary.append(q)

    # repeat step 2 w/ wider range if this happens
    at_lower_boundary = len(quantiles_at_lower_boundary)>0
    at_upper_boundary = len(quantiles_at_upper_boundary)>0
    if at_lower_boundary or at_upper_boundary:
        fprint("WARNING: found limits for quantiles {} at boundary".format(','.join([str(q) for q in quantiles_at_lower_boundary+quantiles_at_upper_boundary])))

    # draw plots
    if args.plots:
        r.gStyle.SetOptStat(0)
        r.gStyle.SetOptTitle(0)

        hist = r.TH1D('hist','hist',100,rmin,rmax)
        hist.GetXaxis().SetTitle('#mu')
        c = r.TCanvas('c','c',700,550)
        hist.Draw()
        hist.SetMaximum(1)
        hist.SetMinimum(0)
        cls_graph.Draw('l')
        cls_graph.SetLineColor(r.kBlack)
        cls_graph.SetLineStyle(7)
        clsb_graph.Draw('l')
        clsb_graph.SetLineColor(r.kRed)
        clb_graph.Draw('l')
        clb_graph.SetLineColor(r.kBlue)
        clb_graph.SetLineStyle(7)
        for q in quantiles:                                       
            cls_exp_graph[q].Draw('l')
            if q==0.5: cls_exp_graph[q].SetLineColor(r.kGray+2)
            if q==0.16 or q==0.84: cls_exp_graph[q].SetLineColor(r.kGray+1)
            if q==0.025 or q==0.975: cls_exp_graph[q].SetLineColor(r.kGray)
        line = r.TLine(rmin, 0.05, rmax, 0.05)
        line.SetLineColor(r.kGreen)
        line.Draw()

        leg = r.TLegend(0.68,0.17,0.89,0.36)

        leg.SetTextFont(42)
        leg.SetFillColor(r.kWhite)
        leg.SetLineColor(r.kWhite)
        leg.SetFillStyle(0)
        leg.SetLineWidth(0)
        leg.AddEntry(clb_graph, "CL_{b}","l")
        leg.AddEntry(clsb_graph, "CL_{s+b}","l")
        leg.AddEntry(cls_graph, "CL_{s}","l")
        leg.AddEntry(cls_exp_graph[0.5], "CL_{s} exp.","l")
        leg.Draw()
        pname1 = "cls_{}.{}"
        for pformat in ["png","pdf"]:
            c.Print(pname1.format(args.name,pformat))

        hist.Draw()
        hist.SetMaximum(10)
        dn2ll_data_graph.Draw('l')
        dn2ll_data_graph.SetLineColor(r.kBlack)
        dn2ll_asimov_graph.Draw('l')
        dn2ll_asimov_graph.SetLineColor(r.kBlue)

        leg = r.TLegend(0.68,0.17,0.89,0.28)
        leg.SetTextFont(42)
        leg.SetFillColor(r.kWhite)
        leg.SetLineColor(r.kWhite)
        leg.SetFillStyle(0)
        leg.SetLineWidth(0)
        leg.AddEntry(dn2ll_asimov_graph, "#tilde{q}_{#mu,A}","l")
        leg.AddEntry(dn2ll_data_graph, "#tilde{q}_{#mu}","l")
        leg.Draw()

        pname2 = "dn2ll_{}.{}"
        for pformat in ["png","pdf"]:
            c.Print(pname2.format(args.name,pformat))

    return limits, at_lower_boundary, at_upper_boundary

def step4(args, limits):
    # run MDF for each r value to get output tree w/ proper fit params, normalizations, etc.
    # include prefit (bkg-only) as quantile=-2 w/ r=0
    limits[-2] = 0.
    
    ofnames = {}
    index = 0
    no_reuse = "step4" not in args.reuse
    for q, rval in sorted(limits.iteritems()):
        args4 = updateArg(args.args, ["--setParameters"], "r={}".format(rval), ',')
        args4 = updateArg(args4, ["--freezeParameters"], "r", ',')
        args4 = updateArg(args4, ["-n","--name"], "Postfit{}".format(index))
        cmd4 = "combine -M MultiDimFit {} {}".format(args4,args.fitopts)
        fprint(cmd4)
        logfname4 = "log_step4q{}_{}.log".format(q,args.name)
        ofname4 = ""
        if not args.dry_run:
            if no_reuse: runCmd(cmd4,logfname4)
            ofname4, success = getOutfile(logfname4)
            retries = 0
            max_retries = 10
            rval_old = rval
            rval_new = rval
            while no_reuse and not success and retries <= max_retries:
                retries += 1
                rval_old = rval_new
                rval_new = rval-0.00000000001*retries
                fprint("MultiDimFit failed for quantile {}, trying a perturbation: r={}".format(q,rval_new))
                cmd4 = cmd4.replace("r={}".format(rval_old),"r={}".format(rval_new))
                runCmd(cmd4,logfname4)
                ofname4, success = getOutfile(logfname4)
            if not success:
                fprint("WARNING: MultiDimFit failed {} times, giving up".format(retries))
            ofnames[q] = ofname4
        index += 1
        
    return ofnames
    
def step5(args, limits, fits, ofname1):
    # should this be a separate step?
    if not args.dry_run:
        import ROOT as r
        r.gROOT.ProcessLine("struct quantile_t { Float_t quantile; Double_t limit; };")
        qobj = r.quantile_t()
        
        # combine trees, setting quantile values
        if len(fits)>0:
            trees = r.TList()
            for q,ofname in sorted(fits.iteritems()):
                file_q = r.TFile.Open(ofname)
                tree_q = file_q.Get("limit")
                tree_q.SetDirectory(0)
                tree_q.SetBranchAddress("quantileExpected",r.AddressOf(qobj,'quantile'))
                tree_q.SetBranchAddress("limit",r.AddressOf(qobj,'limit'))
                tree_q_new = tree_q.CloneTree(0)
                tree_q_new.SetDirectory(0)
                tree_q.GetEntry(0)
                qobj.quantile = q
                qobj.limit = limits[q]
                tree_q_new.Fill()
                trees.Add(tree_q_new)
            newtree = r.TTree.MergeTrees(trees)
        # reuse step1 tree
        else:
            file1 = r.TFile.Open(ofname1)
            tree1 = file1.Get("limit")
            tree1.SetDirectory(0)
            tree1.SetBranchAddress("quantileExpected",r.AddressOf(qobj,'quantile'))
            tree1.SetBranchAddress("limit",r.AddressOf(qobj,'limit'))
            newtree = tree1.CloneTree(0)
            for i in range(tree1.GetEntries()):
                tree1.GetEntry(i)
                for q,rval in limits.iteritems():
                    if abs(q-qobj.quantile)<0.01:
                        qobj.limit = rval
                        break
                newtree.Fill()
        
        # output
        ofile = r.TFile.Open(ofname1.replace("AsymptoticLimits","ManualCLs").replace("Step1",""),"RECREATE")
        ofile.cd()
        newtree.Write()
        ofile.Close()

def manualCLs(args):
    # 1. estimate r range
    ofname1 = step1(args)

    # repeat steps 2 and 3 if boundaries are hit
    # todo: add option to "refine" (smaller range)
    at_lower = True
    at_upper = True
    count_lower = 0
    count_upper = 0
    while at_lower or at_upper:
        # 2. run likelihood scans
        scans = step2(args, ofname1, count_lower, count_upper)

        # 3. compute CLs from likelihood scans
        limits, at_lower, at_upper = step3(args, scans)
        if at_lower: count_lower += 1
        if at_upper: count_upper += 1

    # 4. run MDF for each r value to get fit params
    fits = []
    if args.fit:
        fits = step4(args, limits)
    
    # 5. make new limit tree from step 4 MDF runs
    step5(args, limits, fits, ofname1)

if __name__=="__main__":
    reusable_steps = ["step1","step2","step4"]

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-a", "--args", dest="args", type=str, required=True, help="input arguments for combine")
    parser.add_argument("-n", "--name", dest="name", type=str, default="", help="name for output files")
    parser.add_argument("-p", "--plots", dest="plots", default=False, action='store_true', help="make likelihood plots")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-r", "--reuse", dest="reuse", type=str, default=[], nargs='*', choices=reusable_steps + ["all"], help="reuse Combine results from specified steps")
    parser.add_argument("-f", "--fit", dest="fit", default=False, action='store_true', help="run MDF for prefit and postfit")
    parser.add_argument("-x", "--extra", dest="extra", default=False, action='store_true', help="enable extra fit options for MDF")
    args = parser.parse_args()

    if "all" in args.reuse: args.reuse = reusable_steps[:]
    args.reuse = set(args.reuse)

    if len(args.name)==0: args.name = uuid.uuid4()

    args.fitopts = '--setRobustFitStrategy 0 --setRobustFitTolerance 0.1 --robustFit 1 --cminPreScan --cminPreFit 1 --cminFallbackAlgo "Minuit2,0:0.2"  --cminFallbackAlgo "Minuit2,0:1.0" --cminFallbackAlgo "Minuit2,0:10.0" --cminOldRobustMinimize 0 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_MaxCalls=9999999' if args.extra else ''

    manualCLs(args)

