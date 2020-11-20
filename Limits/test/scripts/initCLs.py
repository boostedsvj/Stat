import os,sys,subprocess,shlex,uuid
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from getParamsTracked import getParamsTracked

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

def step2(args, ofname1):
    # obtain parameter values from step1 (median exp)
    params = getParamsTracked(0, ofname1.split('.')[0].replace("higgsCombine",""), ofname1.split('.')[1], 0.5, True, False, True)
    combo = ofname1.split('.')[3].replace("ana","")
    
    # run AsymptoticLimits w/ initial values from above
    args2 = updateArg(args.args, ["--setParameters"], ','.join(params[combo]))
    args2 = updateArg(args2, ["-n","--name"], "Step2")
    cmd2 = "combine -M AsymptoticLimits "+args2
    fprint(cmd2)
    logfname2 = "log_step2_{}.log".format(args.name)
    ofname2 = ""
    if not args.dry_run:
        if "step2" not in args.reuse: runCmd(cmd2,logfname2)
        ofname2, _ = getOutfile(logfname2)
    if args.doprint: fprint(','.join(params[combo]))
    return ofname2

def initCLs(args):
    # 1. get initial values
    ofname1 = step1(args)

    # 2. use initial values
    ofname2 = step2(args, ofname1)

if __name__=="__main__":
    reusable_steps = ["step1","step2","step4"]

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-a", "--args", dest="args", type=str, required=True, help="input arguments for combine")
    parser.add_argument("-n", "--name", dest="name", type=str, default="", help="name for output files")
    parser.add_argument("-D", "--dry-run", dest="dry_run", default=False, action='store_true', help="dry run (print commands but don't execute)")
    parser.add_argument("-p", "--print", dest="doprint", default=False, action='store_true', help="print parameter values")
    parser.add_argument("-r", "--reuse", dest="reuse", type=str, default=[], nargs='*', choices=reusable_steps + ["all"], help="reuse Combine results from specified steps")
    args = parser.parse_args()

    if "all" in args.reuse: args.reuse = reusable_steps[:]
    args.reuse = set(args.reuse)

    if len(args.name)==0: args.name = uuid.uuid4()

    initCLs(args)

