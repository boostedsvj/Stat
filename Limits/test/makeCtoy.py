import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import getBiasArgs
from paramUtils import fprint, makeSigDict, getParamNames, getSigname, getSignameShort, getFname, getWname, getDname, getDCname, getCombos, runCmds, getChannel, PoissonErrorUp

def main(args):
    import ROOT as r
    r.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

    current_dir = os.getcwd()
    signame = getSigname(args.signal) if len(args.signal)>0 else ""
    combos = getCombos()
    toyname = "{}{}toy{}".format("a" if args.asimov else "b", "sig" if len(args.inject)!=0 else "", "_"+signame if len(args.inject)!=0 else "")
    outname = "hists.{}.{}.root".format(toyname,args.seed)
    outfile = r.TFile.Open(outname,"RECREATE")
    outputs = []
    for combo,regions in combos.iteritems():
        os.chdir(signame)
        setargs = []
        frzargs = []
        datacards = []
        for region in regions:
            wname = getWname(signame,region)
            biasargs = getBiasArgs.main(wname, args.default_ws, region, args.function, verbose=False)
            setargs.extend(biasargs['SetArg'])
            frzargs.extend(biasargs['FrzArg'])
            datacards.append(getDname(signame, region))
        dcfname = getDCname(args.signal, combo)
        toyfname = getFname(toyname,"GenerateOnly",combo,seed=args.seed)

        # make the toy, put in base directory
        if len(args.inject)==0: inject = 0
        elif len(args.inject)==1: inject = args.inject[0]
        else: inject = args.inject[0 if combo=="cut" else 1]
        command = "combineCards.py "+" ".join(datacards)+" > "+dcfname
        fprint(command)
        os.system(command)
        commands = [
            "combine -M GenerateOnly -n {} -t {} -s {} --setParameters {} --freezeParameters {} --keyword-value ana={} -d {} --expectSignal {} -v -1 --toysFrequentist --saveToys --bypassFrequentistFit --saveWorkspace".format(
                toyname, -1 if args.asimov else 1, args.seed, ','.join(setargs), ','.join(frzargs), combo, dcfname, inject
            ),
            "mv {} ../".format(toyfname),
        ]
        runCmds(commands)

        # escape from roofit
        os.chdir(current_dir)
        toyfile = r.TFile.Open(toyfname)
        for region in regions:
            channel = getChannel(region)
            #fullregion = "{}_2018".format(region)
            fullregion = "{}".format(region)
            hdata = toyfile.Get("toys/toy_{}".format("asimov" if args.asimov else "1")).reduce("CMS_channel==CMS_channel::{}".format(channel))
            hist = r.RooAbsData.createHistogram(hdata, "mH{}".format(fullregion))
            # fix error bars
            for b in range(hist.GetNbinsX()):
                bin = b+1
                content = hist.GetBinContent(bin)
                if content>0:
                    hist.SetBinError(bin,PoissonErrorUp(content))
            tmpdir = outfile.mkdir(fullregion)
            tmpdir.cd()
            hist.SetName("data_toy")
            hist.Write()

        outputs.append(toyfname)
    outputs.append(outname)
    fprint('\n'.join(outputs))

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--signal", dest="signal", metavar=tuple(getParamNames()), type=str, default=["3100","20","03","peak"], nargs=4, help="signal parameters")
    parser.add_argument("-i", "--inject", dest="inject", type=float, default=[], nargs='+', help="injected signal strength (one arg, or two args: cut bdt)")
    parser.add_argument("-a", "--asimov", dest="asimov", default=False, action="store_true", help="use asimov dataset")
    parser.add_argument("-n", "--function", dest="function", type=int, required=True, choices=[0,1], help="function (0: main, 1: alt)")
    parser.add_argument("-S", "--seed", dest="seed", type=int, default=123456, help="random seed")
    args = parser.parse_args()

    args.signal = makeSigDict(args.signal)
    # pass some defaults
    args.default_ws = "SVJ"

    main(args)

