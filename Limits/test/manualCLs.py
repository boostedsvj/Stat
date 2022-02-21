import os,sys,subprocess,shlex,uuid
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict

#########################
# helper functions
#########################

# make status messages useful
def fprint(msg):
    import sys
    print(msg)
    sys.stdout.flush()

# common structure for interpolations
class GenericInterpolator(object):
    def __init__(self, x, y):
        from array import array
        self.x = array('d', x)
        self.y = array('d', y)
    def eval(self, val):
        pass

class LinearInterpolator(GenericInterpolator):
    def __init__(self, x, y):
        super(LinearInterpolator,self).__init__(x,y)
        # just use TGraph
        import ROOT as r
        self.graph = r.TGraph(len(self.x),self.x,self.y)
    def eval(self, val):
        return self.graph.Eval(val,0,"")

class CubicInterpolator(GenericInterpolator):
    def __init__(self, x, y):
        # use TSpline3 (requires sorted data)
        self.x, self.y = zip(*sorted(zip(x,y),key=lambda p: p[0]))
        super(CubicInterpolator,self).__init__(x,y)
        import ROOT as r
        self.spline = r.TSpline3("", self.x, self.y, len(self.x))
    def eval(self, val):
        return self.spline.Eval(val)

class SteffenInterpolator(GenericInterpolator):
    def __init__(self, x, y):
        super(SteffenInterpolator,self).__init__(x,y)
        # use gsl via ROOT pythonization
        import ROOT as r
        def pythonize_gsl():
            try:
                r.SteffenWrapper
            except:
                r.gROOT.ProcessLine("""#include <gsl/gsl_spline.h>
class SteffenWrapper {
public:
SteffenWrapper(int n, double* x, double* y){
    spline_ = gsl_spline_alloc(gsl_interp_steffen,n);
    gsl_spline_init(spline_,x,y,n);
}
~SteffenWrapper() { gsl_spline_free(spline_); }
double Eval(double x) { return gsl_spline_eval(spline_, x, 0); }
gsl_spline* spline_;
};""")
        pythonize_gsl()
        self.spline = r.SteffenWrapper(len(self.x), self.x, self.y)
    def eval(self, val):
        return self.spline.Eval(val)

class PchipInterpolator(GenericInterpolator):
    def __init__(self, x, y):
        super(PchipInterpolator,self).__init__(x,y)
        # use scipy
        from scipy import interpolate
        self.spline = interpolate.PchipInterpolator(x,y)
    def eval(self, val):
        return self.spline(val)

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

def replaceArg(args, flags, val):
    # make a copy
    argsplit = args.split()
    if not isinstance(flags,list): flags = [flags]
    for flag in flags:
        for iarg,arg in enumerate(argsplit):
            if flag==arg:
                argsplit[iarg+1] = val
                break
    return ' '.join(argsplit)

def removeArg(args, matches, before=0, after=0, exact=True):
    # make a copy
    argsplit = args.split()
    if not isinstance(matches,list): matches = [matches]
    toremove = ""
    for match in matches:
        for iarg,arg in enumerate(argsplit):
            if (exact and match==arg) or (not exact and match in arg):
                toremove = ' '.join(argsplit[iarg-before:iarg+after+1])
                break
    args = args.replace(toremove,"",1)
    return args

def removeToyArgs(args):
    args = removeArg(args, "--toysFile", after=1)
    args = removeArg(args, "-t", after=1)
    args = removeArg(args, "--toysFrequentist")
    return args

def runCmd(cmd, log, opt='w'):
    if not isinstance(log,int): logfile = open(log,opt)
    else: logfile = log
    subprocess.call(shlex.split(cmd),stdout=logfile,stderr=logfile)
    if not isinstance(log,int): logfile.close()

def getOutfile(log):
    ofname = ""
    indicator = "COMBINE_OUTPUT_FILE: "
    # MDF and FD have different failure messages
    failure = ["WARNING: MultiDimFit failed","Fit failed."]
    success = True
    with open(log,'r') as logfile:
        for line in logfile:
            # stop if found both checks
            if len(ofname)>0 and not success: break
        
            if indicator in line:
                ofname = line.rstrip().replace(indicator,"")
            elif any(f in line for f in failure):
                success = False

    if len(ofname)==0: raise RuntimeError("Could not find output file name from log: {}".format(log))
    return ofname, success

def getRange(dry_run, ofname1, nuisances, permissive):
    # get r range
    # (default vals provided for dryrun printouts)
    rmin = 0.
    rmax = 10.
    factor = 10.
    if nuisances: factor = 2.
    npts = 100
    # increase npts by 1 to include both endpoints
    npts += 1
    if not dry_run:
        import ROOT as r
        file1 = r.TFile.Open(ofname1)
        limit1 = file1.Get("limit")
	#print("filename is: "file1.GetName())
        entries = limit1.GetEntries()
        n = limit1.Draw("limit","abs(quantileExpected-0.975)<0.001||abs(quantileExpected-0.025)<0.001","goff")
        vals = []
        if n!=2:
            if not permissive or entries==0: raise RuntimeError("Malformed limit tree in "+ofname1)
            elif permissive:
                n = limit1.Draw("limit","","goff")
                vals = [limit1.GetVal(0)[0]]
        if len(vals)==0: vals = [limit1.GetVal(0)[0],limit1.GetVal(0)[1]]
        rmax = max(vals)*factor

    return rmin,rmax,npts

def getBranches(tree, matches=None, exact=False):
    import ROOT as r
    if not exact and not isinstance(matches,list): matches = [matches]
    elif exact and isinstance(matches,list): matches = matches[0]
    branches = []
    print('tree in getBranches: ', tree)
    for b in tree.GetListOfBranches():
        bname = b.GetName()
        leaf = b.GetLeaf(bname)
        if matches is None or (exact and matches==b.GetName()) or (not exact and all(m in bname for m in matches)):
            branches.append(bname)
    return branches

def replaceMulti(string,replacements):
    for old,new in replacements.iteritems():
        string = string.replace(old,new)
    return string

def getParams(dry_run, ofname1, pnames, quantile, matches=[]):
    if dry_run: return

    import ROOT as r

    condition = "abs(quantileExpected-{})<0.001".format(quantile)

    file1 = r.TFile.Open(ofname1)
    limit1 = file1.Get("limit")

    params = []
    for p in pnames:
        params.extend(getBranches(limit1, matches+[p]))

    # deliver string of param=value
    limit1.Draw(':'.join(params),condition,"goff")
    results = {p:limit1.GetVal(i)[0] for i,p in enumerate(params)}
    return results

def getParamsTxt(params):
    return ','.join(["{}={}".format(replaceMulti(p,{'trackedParam_':'','trackedError_':''}),v) for p,v in sorted(params.iteritems())])

def getParamConstraints(params,errors):
    constraints = {}
    psuff = 'trackedParam_'
    esuff = 'trackedError_'
    for ename in errors:
        name = ename.replace(esuff,'')
        pname = psuff+name
        constraints[name] = (params[pname]-errors[ename],params[pname]+errors[ename])
    return constraints

def getConstraintsTxt(constraints):
    return ':'.join(["{}={},{}".format(p,v[0],v[1]) for p,v in sorted(constraints.iteritems())])

def getInitArgs(args, ofname, quantile):
    # optional: initialize parameter values from step1 (initCLs) or step0 (bonly)
    init_args = ""
    constraint_args = ""
    params = getParams(args.dry_run,ofname,args.initCLs,quantile,matches=['trackedParam'])
    init_args = getParamsTxt(params)
    if len(args.initConstraints)>0:
        errors = getParams(args.dry_run,ofname,args.initConstraints,quantile,matches=['trackedError'])
        constraints = getParamConstraints(params,errors)
        constraint_args = getConstraintsTxt(constraints)
    init_args = (init_args,constraint_args)
    return init_args

def handleInitArgs(args, init_args):
    if init_args is not None:
        params = init_args[0]
        constraints = init_args[1]
        if len(params)>0: args = updateArg(args, ['--setParameters'], params, ',')
        if len(constraints)>0: args = updateArg(args, ['--setParameterRanges'], constraints, ':')
    return args

#########################
# step functions
#########################

def step0(args, products):
    # run MultiDimFit w/ fixed r=0
    args0 = updateArg(args.args, ["--freezeParameters"], "r", ',')
    args0 = updateArg(args0, ["--setParameters"], "r=0", ',')
    args0 = updateArg(args0, ["-n","--name"], "Step0")
    cmd0 = "combine -M MultiDimFit --saveShapes --saveNLL {} {}".format(args0,args.fitopts)
    fprint(cmd0)
    logfname0 = "log_step0_{}.log".format(args.name)
    products["ofname0"] = ""
    products["init_args"] = None
    if not args.dry_run:
        if "step0" not in args.reuse: runCmd(cmd0,logfname0)
        products["ofname0"], _ = getOutfile(logfname0)
        products["init_args"] = getInitArgs(args, products["ofname0"], -1)
    return products

def step1(args, products, nuisances=False, stepname="Step1"):
#def step1(args, products, nuisances=False, stepname="Step"):
    # run AsymptoticLimits w/ nuisances disabled (default)
    args1 = updateArg(args.args, ["-n","--name"], stepname)
    print("args1: ", args1)
    print("args: ", args)
    print("products: ", products)
    if not nuisances: args1 = updateArg(args1, ["--freezeParameters"], "allConstrainedNuisances", ',')
    args1 = handleInitArgs(args1, products["init_args"])
    cmd1 = "combine -M AsymptoticLimits "+args1
    fprint(cmd1)
    logfname1 = "log_{}_{}.log".format(stepname[0].lower()+stepname[1:],args.name)
    print('logfname1: ', logfname1)
    products["ofname1"] = ""
    if not args.dry_run:
        if "step1" not in args.reuse: runCmd(cmd1,logfname1)
        products["ofname1"], _ = getOutfile(logfname1)
        if len(args.initCLs)>0 and products["init_args"] is None:
            products["init_args"] = getInitArgs(args, products["ofname1"], 0.5)
    return products

def stepA(args, products):
    return step1(args, products, True, "StepA")

def step2impl(args, products, name, lname, ofname, extra=""):
    # rename file from previous r range
    ofname_old = None
    if not args.dry_run and products[ofname] is not None and "step2" not in args.reuse:
        ofname_old = products[ofname].replace(".","0.",1)
        os.rename(products[ofname], ofname_old)

    # run MDF likelihood scan
    args2 = updateArg(args.args, ["-n","--name"], name)
    if "toysFile" in extra:
        args2 = removeToyArgs(args2)
    args2 = removeArg(args2, "--saveToys")
    args2 = handleInitArgs(args2, products["init_args"])
    args2 = updateArg(args2, ['--setParameterRanges'], "r={},{}".format(products["rmin"],products["rmax"]), ':')
    # addl args for second r ranges
    npts = products["npts"]
    if ofname_old is not None:
        npts -= 1
        extra += " --snapshotName MultiDimFit --skipInitialFit"
        args2 = replaceArg(args2, ['-d', '--datacard'], ofname_old)
    cmd2 = "combine -M MultiDimFit --redefineSignalPOIs r --algo grid --points {} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --alignEdges 1 --saveWorkspace --saveNLL {} {} {}".format(products["npts"],extra,args2,args.fitopts)
    fprint(cmd2)

    logfname2 = "log_{}_{}.log".format(lname,args.name)
    products[ofname] = ""
    if not args.dry_run:
        if "step2" not in args.reuse: runCmd(cmd2,logfname2)
        products[ofname], _ = getOutfile(logfname2)

        # combine w/ previous file
        if ofname_old is not None:
            ofname_new = products[ofname].replace(".","1.",1)
            os.rename(products[ofname], ofname_new)
            cmdh = "hadd {} {} {}".format(products[ofname], ofname_old, ofname_new)
            fprint(cmdh)
            if "step2" not in args.reuse:
                runCmd(cmdh, subprocess.PIPE)
                for ofname_ in [ofname_old,ofname_new]: os.remove(ofname_)
    print("product: ", products)
    return products
    #print("product: ", products)

def step2(args, products):
    asimov_args = "-t -1 --toysFreq"

    if products["count_upper"]==0:
        # get rmin, rmax from step1
        products["rmin"], products["rmax"], products["npts"] = getRange(args.dry_run, products["ofname1"], args.syst, args.permissive)

        # get asimov dataset separately (for some reason, hadding MultiDimFit output files crashes if both --saveWorkspace and --saveToys are used)
        argsG = updateArg(args.args, ["-n","--name"], "Asimov")
        argsG = removeToyArgs(argsG)
        argsG = removeArg(argsG, "--saveToys")
        argsG = handleInitArgs(argsG, products["init_args"])
        cmdG = "combine -M GenerateOnly {} --saveToys {}".format(asimov_args, argsG)
        fprint(cmdG)
        logfnameG = "log_{}_{}.log".format("step2g",args.name)
        products["ofname2g"] = ""
        if not args.dry_run:
            if "step2" not in args.reuse: runCmd(cmdG, logfnameG)
            products["ofname2g"], _ = getOutfile(logfnameG)
    else:
        # go to next r range
        range_factor = 2
        products["rmin"] = products["rmax"]
        products["rmax"] = products["rmax"]*range_factor
        # avoid rerunning previous rmax
        products["rmin"] += (products["rmax"]-products["rmin"])/float(products["npts"])

    # observed
    products = step2impl(args, products, "Observed", "step2d", "ofname2d")

    # expected (asimov)
    products = step2impl(args, products, "Asimov", "step2a", "ofname2a", extra="{} --toysFile {}".format(asimov_args, products["ofname2g"]))

    return products

def step3(args, products):
    # based on https://gitlab.cern.ch/cms-hcg/cadi/hig-19-003/-/blob/master/HJMINLO/plot_cls.py
    interpolate = args.spline is not None
    inter_npoints = args.npoints
    quantiles = [0.025, 0.16, 0.50, 0.84, 0.975]

    if args.dry_run:
        # return dummy dictionary to be used in step4 commands
        products["limits"] = {q:0.0 for q in quantiles+[-1]}
        products["at_upper"] = False
        return products

    import ROOT as r
    # put root in batch mode
    r.gROOT.SetBatch()
    file_data = r.TFile.Open(products["ofname2d"])
    file_asimov = r.TFile.Open(products["ofname2a"])
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

    fprint("INFO: r_data_bestfit = {}, r_asimov_bestfit = {}".format(r_data_bestfit,r_asimov_bestfit))
    if r_data_bestfit<0: fprint("WARNING: data bestfit failed!")
    if r_asimov_bestfit<0: fprint("WARNING: asimov bestfit failed!")
    r_dict = {k:v for k,v in r_dict.iteritems() if v[1] is not None}
    r_data = list(sorted(r_dict))
    rmin = min(r_data)
    rmax = max(r_data)
    dnll_data = [r_dict[k][0] for k in r_data]
    dnll_asimov = [r_dict[k][1] for k in r_data]
    len_r_data = len(r_data)
    if products["count_upper"]>0 and len_r_data==products["len_r_data"]:
        raise RuntimeError("ERROR: no fits converged when extending r range to avoid limits at boundary")

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

    dn2ll_min = 1e10
    dn2ll_max = -1e10
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
        dn2ll_min = min(dn2ll_min,2*dnll_data_constrained,2*dnll_asimov[i])
        dn2ll_max = max(dn2ll_max,2*dnll_data_constrained,2*dnll_asimov[i])

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
        r_data_interp = []
        cls_helper = args.spline(r_data,cls)
        cls = []
        cls_exp_helpers = {}
        for q in quantiles:
            cls_exp_helpers[q] = args.spline(r_data,cls_exp[q])
            cls_exp[q] = []
        # fill in more points between each pair of r values
        for rind in range(len(r_data)-1):
            r1 = r_data[rind]
            r2 = r_data[rind+1]
            for rval in np.linspace(r1, r2, inter_npoints+1):
                r_data_interp.append(rval)
                cls.append(cls_helper.eval(rval))
                for q in quantiles:
                    cls_exp[q].append(cls_exp_helpers[q].eval(rval))
        r_data = r_data_interp

    products["limits"] = {}
    quantiles_at_upper_boundary = []
    alpha = 0.05
    idx = find_crossing(cls,alpha)
    products["limits"][-1] = r_data[idx]
    fprint("Observed Limit: r < {:.2f}".format(products["limits"][-1]))
    if idx==len(r_data)-1 or idx==-1: quantiles_at_upper_boundary.append(-1)
    for q in quantiles:
        qidx = find_crossing(cls_exp[q],alpha)
        products["limits"][q] = r_data[qidx]
        fprint("Expected {:3.1f}%: r < {:.2f}".format(q*100., products["limits"][q]))
        if qidx==len(r_data)-1 or qidx==-1: quantiles_at_upper_boundary.append(q)

    # store best fit values
    products["limits"][-3] = r_data_bestfit
    products["limits"][-4] = r_asimov_bestfit

    # repeat step 2 w/ wider range if this happens
    products["at_upper"] = len(quantiles_at_upper_boundary)>0
    if products["at_upper"]:
        fprint("WARNING: found limits for quantiles {} at boundary".format(','.join([str(q) for q in quantiles_at_upper_boundary])))
        products["len_r_data"] = len_r_data

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
        cls_graph.SetLineWidth(2)
        cls_graph.SetLineStyle(7)
        clsb_graph.Draw('l')
        clsb_graph.SetLineColor(r.kRed)
        clsb_graph.SetLineWidth(2)
        clb_graph.Draw('l')
        clb_graph.SetLineColor(r.kBlue)
        clb_graph.SetLineWidth(2)
        clb_graph.SetLineStyle(7)
        for q in quantiles:                                       
            cls_exp_graph[q].Draw('l')
            if q==0.5: 
              cls_exp_graph[q].SetLineColor(r.kGray+2)
              cls_exp_graph[q].SetLineWidth(2)
            if q==0.16 or q==0.84: 
	      cls_exp_graph[q].SetLineColor(r.kGray+1)
              cls_exp_graph[q].SetLineWidth(2)
            if q==0.025 or q==0.975: 
	      cls_exp_graph[q].SetLineColor(r.kGray)
  	      cls_exp_graph[q].SetLineWidth(2)
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
        hist.GetYaxis().SetRangeUser(dn2ll_min,dn2ll_max)
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

    return products

def step4(args, products):
    # run MDF for each r value to get output tree w/ proper fit params, normalizations, shapes, etc.
    # include: prefit (bkg-only) as quantile=-2 w/ r=0
    #          bestfit (obs) as quantile=-3
    #          bestfit (asimov) as quantile=-4
    products["limits"][-2] = 0.

    def checkFitRes(ofname,success):
        if not success: return success
        import ROOT as r
        ofname = ofname.replace("MultiDimFit.","").replace("higgsCombine","multidimfit")
        f = r.TFile.Open(ofname)
        fitres = f.Get("fit_mdf")
        return fitres!=None
    
    products["ofitnames"] = {}
    no_reuse = "step4" not in args.reuse
    for q, rval in sorted(products["limits"].iteritems()):
        args4 = args.args[:]
        extra = ""
        if q==-3 or q==-4:
            extra = "--X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL"
        if q==-4 or q>0:
            extra += " -t -1 --toysFreq"
            args4 = removeToyArgs(args4)
        else:
            args4 = updateArg(args4, ["--setParameters"], "r={}".format(rval), ',')
            args4 = updateArg(args4, ["--freezeParameters"], "r", ',')
        args4 = updateArg(args4, ["-n","--name"], "Postfit{:.3f}".format(q))
        # not needed when params already known
        args4 = removeArg(args4,"MINIMIZER_MaxCalls",before=1,exact=False)
        # avoid default cutoff
        args4 += " --rMax {}".format(rval*2)
        args4 = handleInitArgs(args4, products["init_args"])
        cmd4 = "combine -M MultiDimFit --saveShapes --saveWorkspace --saveFitResult {} {} {}".format(extra,args4,args.fitopts)
        fprint(cmd4)
        logfname4 = "log_step4q{}_{}.log".format(q,args.name)
        ofname4 = ""
        if not args.dry_run:
            if no_reuse: runCmd(cmd4,logfname4)
            ofname4, success = getOutfile(logfname4)
            success = checkFitRes(ofname4,success)
            retries = 0
            max_retries = 10
            rval_old = rval
            rval_new = rval
            while no_reuse and not success and retries <= max_retries and q!=-2:
                retries += 1
                rval_old = rval_new
                rval_new = rval-0.00000000001*retries
                fprint("MultiDimFit failed for quantile {}, trying a perturbation: r={}".format(q,rval_new))
                cmd4 = cmd4.replace("r={}".format(rval_old),"r={}".format(rval_new))
                runCmd(cmd4,logfname4)
                ofname4, success = getOutfile(logfname4)
                success = checkFitRes(ofname4,success)
            if not success:
                fprint("WARNING: MultiDimFit failed {} times, giving up".format(retries))
            products["ofitnames"][q] = ofname4
        
    return products
    
def step5(args, products, title="ManualCLs"):
    if not args.dry_run:
        import ROOT as r
        try:
            r.quantile_t
        except:
            r.gROOT.ProcessLine("struct quantile_t { Float_t quantile; Double_t limit; };")
        qobj = r.quantile_t()
        
        # combine trees, setting quantile values
        if title=="ManualCLsFit":
            if products["ofitnames"] is None:
                raise RuntimeError("Can't run step5 w/ title {} w/o step4 fits".format(title))
            trees = r.TList()
            for q,ofname in sorted(products["ofitnames"].iteritems()):
                file_q = r.TFile.Open(ofname)
                tree_q = file_q.Get("limit")
                tree_q.SetDirectory(0)
                tree_q.SetBranchAddress("quantileExpected",r.AddressOf(qobj,'quantile'))
                tree_q.SetBranchAddress("limit",r.AddressOf(qobj,'limit'))
                tree_q_new = tree_q.CloneTree(0)
                tree_q_new.SetDirectory(0)
                tree_q.GetEntry(0)
                qobj.quantile = q
                qobj.limit = products["limits"][q]
                tree_q_new.Fill()
                trees.Add(tree_q_new)
            newtree = r.TTree.MergeTrees(trees)
        # reuse step1 tree
        elif title=="ManualCLs":
            file1 = r.TFile.Open(products["ofname1"])
            tree1 = file1.Get("limit")
            tree1.SetDirectory(0)
            tree1.SetBranchAddress("quantileExpected",r.AddressOf(qobj,'quantile'))
            tree1.SetBranchAddress("limit",r.AddressOf(qobj,'limit'))
            newtree = tree1.CloneTree(0)
            # todo: handle permissive case here
            for i in range(tree1.GetEntries()):
                tree1.GetEntry(i)
                for q,rval in products["limits"].iteritems():
                    if abs(q-qobj.quantile)<0.01:
                        qobj.limit = rval
                        break
                newtree.Fill()
        else:
            raise ValueError("Unknown step5 title: {}".format(title))
        
        # output
        products["ofname5"] = products["ofname1"].replace("AsymptoticLimits",title).replace("Step1","")
	#products["ofname5"] = products["ofname1"].replace("AsymptoticLimits",title).replace("Step","")
        ofile = r.TFile.Open(products["ofname5"], "RECREATE")
        ofile.cd()
        newtree.Write()
        ofile.Close()
        return products

#########################
# main functions
#########################

def manualCLs(args):
    products = defaultdict(lambda: None)
    #print('product is: ', products)

    # 0. b-only fit if requested for initCLs
    if args.bonly:
        products = step0(args, products)

    # 1. estimate r range (& params for initCLs)
    if not args.asymptotic:
        products = step1(args, products, args.syst)
    # alternate path: just do asymptotic and exit
    else:
        products = stepA(args, products)
        return

    # repeat steps 2 and 3 if boundaries are hit
    # todo: add option to "refine" (smaller range)
    products["at_upper"] = True
    products["count_upper"] = 0
    while products["at_upper"]:
        # 2. run likelihood scans
        products = step2(args, products)

        # 3. compute CLs from likelihood scans
        products = step3(args, products)
        if products["at_upper"]: products["count_upper"] += 1

        # stop loop if reusing step 2
        if "step2" in args.reuse: products["at_upper"] = False

    # 4. run MDF for each r value to get fit params etc.
    if args.fit:
        products = step4(args, products)
    
    # 5. make new limit tree from step 1 + step 3 (and/or step 4 MDF runs, if available)
    products = step5(args, products)
    if args.fit:
        products = step5(args, products, "ManualCLsFit")

# usage notes:
# patch Combine with https://github.com/kpedro88/HiggsAnalysis-CombinedLimit/commit/dc34ebf8bd4db11814e417d54c74514528bb6b47
# also currently relies on https://github.com/kpedro88/HiggsAnalysis-CombinedLimit/commit/c56f7cac98b5ee90b293eb094bb797b34571b495,
# but will be updated to https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/pull/642
if __name__=="__main__":
    reusable_steps = ["step0","step1","step2","step4"]
    allowed_splines = {
        "none" : None,
        "linear" : LinearInterpolator,
        "cubic" : CubicInterpolator,
        "steffen" : SteffenInterpolator,
        "pchip" : PchipInterpolator,
    }

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-a", "--args", dest="args", type=str, required=True, help="input arguments for combine")
    parser.add_argument("-n", "--name", dest="name", type=str, default="", help="name for output files")
    parser.add_argument("-p", "--plots", dest="plots", default=False, action='store_true', help="make likelihood plots")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-r", "--reuse", dest="reuse", type=str, default=[], nargs='*', choices=reusable_steps + ["all"], help="reuse Combine results from specified steps")
    parser.add_argument("-f", "--fit", dest="fit", default=False, action='store_true', help="run MDF for prefit and postfit")
    parser.add_argument("-x", "--extra", dest="extra", default=False, action='store_true', help="enable extra fit options for MDF")
    parser.add_argument("-i", "--initCLs", dest="initCLs", type=str, default=[], nargs='*', help="use initCLs for specified parameters")
    parser.add_argument("-b", "--bonly", dest="bonly", default=False, action="store_true", help="use b-only fit for initCLs")
    parser.add_argument("-s", "--syst", dest="syst", default=False, action="store_true", help="enable systematics for step1")
    parser.add_argument("-I", "--initConstraints", dest="initConstraints", type=str, default=[], nargs='*', help="use errors from initCLs to constrain ranges for specified parameters")
    parser.add_argument("--robustHesse", dest="robustHesse", default=False, action='store_true', help="enable robustHesse for MDF")
    parser.add_argument("-A", "--asymptotic", dest="asymptotic", default=False, action="store_true", help="just run AsymptoticLimits after init step")
    parser.add_argument("-S", "--spline", dest="spline", type=str, default="linear", choices=list(allowed_splines), help="spline to use for interpolation of CLs curves (none: disable interpolation)")
    parser.add_argument("--npoints", dest="npoints", type=int, default=10, help="number of points to use in interpolation")
    parser.add_argument("--permissive", dest="permissive", default=False, action='store_true', help="be more permissive about finding range in step1")
    args = parser.parse_args()

    if "all" in args.reuse: args.reuse = reusable_steps[:]
    args.reuse = set(args.reuse)

    if len(args.name)==0: args.name = uuid.uuid4()

    args.fitopts = ''
    if args.extra:
        args.fitopts = '--setRobustFitStrategy 0 --setRobustFitTolerance 0.1 --robustFit 1 --cminPreScan --cminPreFit 1 --cminFallbackAlgo "Minuit2,0:0.2"  --cminFallbackAlgo "Minuit2,0:1.0" --cminFallbackAlgo "Minuit2,0:10.0" --cminOldRobustMinimize 0 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_MaxCalls=9999999'
    elif args.robustHesse:
        args.fitopts = '--robustHesse 1'

    # get the class type directly
    args.spline = allowed_splines[args.spline]

    manualCLs(args)


