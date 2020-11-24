import os,sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from getParamsTracked import getParamsTracked, getFname

input_template = """INPUT
input/input_svj_stack_dijetmtdetahad_2017.txt
input/input_svj_mt_postfit_options.txt
input/input_svj_mt_hist_full.txt
"""

options_template = """OPTION
string+:printsuffix[{psuff}]
string:exthisto_dir[{hdir}]
vstring:extra_text[{etxt}]
vstring:fits[{fitlist}]
vstring+:chosensets[{signames}]
vstring+:numers[{signames}]
string:rootfile[{ofile}]
"""

fit_template = "{fitname}\ts:fn[{fnname}]\tvd:pars[1,{pvals}]\td:yield[{yieldval}]\ts:legname[{legname}]\tin:input[input/input_svj_mt_fit_opt.txt]\tb:legpars[0]\tc:linecolor[{fitcol}]"

set_template = """hist\tmc\t{signamefull}\ts:legname[{legname}]\tc:color[{sigcol}]\ti:linestyle[7]\ti:panel[0]\tvs:fits[]\td:yieldnormval[{signorm}]
\tbase\text\t{signamefull}\ts:extfilename[{sigfile}]\tvs:exthisto_in[{signamesafe}]\tvs:exthisto_out[MTAK8]"""

# todo: handle signal-injection toys in legname (& psuff)
data_template = """hist\tdata\tdata\ti:panel[1]\ts:legname[toy data (no signal)]\tvs:extra_text[{fitleg}]\tb:yieldnorm[0]
\tbase\text\tdata\ts:extfilename[{dfile}]\tvs:exthisto_in[Bkg_toy]\tvs:exthisto_out[MTAK8]"""

quantile_info = {
    -3: {"legname": "bestfit", "name": "bestfit", "color": "kOrange + 2", "sigcolor": "kMagenta"},
    -2: {"legname": "b-only", "name": "bonly", "color": "kRed"},
    -1: {"legname": "postfit (obs)", "name": "postfitobs", "color": "kBlue", "sigcolor": "kCyan + 2"},
}

function_info = {
    ("alt",2):  {"formula": "[0]*exp([1]*x/13000+[2]*log(x/13000))", "legname": "g_{2}(x)", "name": "g"},
    ("alt",3):  {"formula": "[0]*exp([1]*x/13000+[2]*log(x/13000)+[3]*log(x/13000)^2)", "legname": "g_{3}(x)", "name": "g"},
    ("main",2): {"formula": "([0]*(1-x/13000)^[1])*((x/13000)^(-[2]))", "legname": "f_{1,1}(x)", "name": "f"},
    ("main",5): {"formula": "([0]*(1-x/13000)^([1]+[2]*log(x/13000)+[3]*log(x/13000)^2))*((x/13000)^(-([4]+[5]*log(x/13000))))", "legname": "f_{3,2}(x)", "name": "f"},
}

region_info = {
    "highCut": {"alt": 3, "main": 5, "legname": "high-^{}R_{T}"},
    "lowCut": {"alt": 3, "main": 2, "legname": "low-^{}R_{T}"},
    "highSVJ2": {"alt": 2, "main": 2, "legname": "high-SVJ2"},
    "lowSVJ2": {"alt": 2, "main": 2, "legname": "low-SVJ2"},
}

def makePostfitPlot(mass, name, method, quantiles, data_file, datacard_dir, combo, region):
    ch = "ch1" if "high" in region else "ch2"

    iname = "input_svj_mt_fit_toy_{region}_{name}_mZprime{mass}.txt".format(region=region,mass=mass,name=name)
    signame = "SVJ_{}_20_0.3_peak".format(mass)
    signamesafe = "SVJ_mZprime{}_mDark20_rinv03_alphapeak".format(mass)
    rinfo = region_info[region]
    ftype = ""
    finfo = None

    fits = OrderedDict()
    sigs = OrderedDict()
    for quantile in quantiles:
        qinfo = quantile_info[quantile]

        params = getParamsTracked(getFname(mass, name, method, combo), quantile)
        if len(params)==0: return ""

        pfit = {p:v for p,v in params.iteritems() if region in p}
        pvals = [str(v) for p,v in sorted(pfit.iteritems())]
        # only pick finfo once, because it should be consistent for the region
        if finfo is None:
            ftype = "alt" if any("_alt" in p for p in pfit) else "main"
            finfo = function_info[(ftype,rinfo[ftype])]
        fname = "{}_{}".format(finfo["name"],qinfo["name"])

        fits[fname] = fit_template.format(
            fitname = fname,
            fnname = finfo["formula"],
            pvals = ','.join(pvals),
            # yieldval must be multiplied by bin width
            yieldval = str(params["trackedParam_n_exp_final_bin{}_proc_roomultipdf".format(ch)]*100),
            legname = qinfo["legname"],
            fitcol = qinfo["color"],
        )

        # no need to show signal for b-only
        if qinfo["name"]=="bonly": continue

        signamefull = "{}_{}".format(signame,qinfo["name"])
        legname = "{} (r = {:.2g})".format(signame,params["r"])

        sigs[signamefull] = set_template.format(
            signamefull = signamefull,
            legname = legname,
            sigcol = qinfo["sigcolor"],
            signorm = str(params["trackedParam_n_exp_final_bin{}_proc_{}".format(ch,signamesafe)]),
            sigfile = "{}/datacard_{}.root".format(datacard_dir,signame),
            signame = signame,
            signamesafe = signamesafe,
        )

    data = data_template.format(
        fitleg = finfo["legname"],
        dfile = data_file,
    )

    options = options_template.format(
        # should quantiles be included here?
        psuff = "_toy_fit_{}_mZprime{}".format(ftype,mass),
        hdir = "{}_2018".format(region),
        etxt = rinfo["legname"],
        fitlist = ','.join(fits),
        signames = ','.join(sigs),
        ofile = "test/fit_toy_mZprime{}_{}".format(mass,region),
    )

    with open(iname,'w') as ifile:
        lines = [
            input_template,
            options,
            "FIT",
            '\n'.join(fits.values()),
            '',
            "SET",
            '\n'.join(sigs.values()),
            data,
        ]
        ifile.write('\n'.join(lines))

    return iname

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--mass", dest="mass", type=int, required=True, help="Zprime mass")
    parser.add_argument("-n", "--name", dest="name", type=str, default="Test", help="test name (higgsCombine[name])")
    parser.add_argument("-M", "--method", dest="method", type=str, default="AsymptoticLimits", help="method name (higgsCombineTest.[method])")
    parser.add_argument("-d", "--data", dest="data", type=str, default="test/trig4_sigfull_datacard.root", help="data file name")
    parser.add_argument("-D", "--datacards", dest="datacards", type=str, default="root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/Datacards/trig4/sigfull/", help="datacard location")
    parser.add_argument("-q", "--quantiles", dest="quantiles", type=float, default=[], nargs='+', help="quantile(s) to plot fits")
    parser.add_argument("-c", "--combos", dest="combos", type=str, default=[], nargs='+', choices=["cut","bdt"], help="combo(s) to plot")
    args = parser.parse_args()

    combos = {
        "cut": ["highCut","lowCut"],
        "bdt": ["highSVJ2","lowSVJ2"],
    }
    args.combos = {c:combos[c] for c in args.combos}

    input_files = []
    for combo,regions in args.combos.iteritems():
        for region in regions:
            tmp = makePostfitPlot(args.mass,args.name,args.method,args.quantiles,args.data,args.datacards,combo,region)
            input_files.append(tmp)
    print ' '.join(input_files)

