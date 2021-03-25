#!/usr/bin/env python

from __future__ import print_function
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from paramUtils import OpenFile

def get_param_names(ws, region, fn):
    allVars = ws.allVars()
    iter = allVars.createIterator()
    v = iter.Next()
    names = []
    while v:
        vname = v.GetName()
        if vname.startswith(region) and ((fn==0 and "alt" not in vname) or (fn==1 and "alt" in vname)):
            names.append(vname)
        v = iter.Next()
    return names

def main(fname, wsname, region, fn, suff="", verbose=True):
    from ROOT import TFile, RooWorkspace
    file = OpenFile(fname)
    ws = file.Get(wsname)

    other_fn = int(not fn)
    indarg = "pdf_index_{}_2018".format(region)

    setargs = ["{}={}".format(indarg,fn)]
    frzargs = [indarg]+get_param_names(ws,region,other_fn)
    trkargs = get_param_names(ws,region,fn)

    if verbose:
        args = [
            'SetArg{}="{}"'.format(suff,','.join(setargs)),
            'FrzArg{}="{}"'.format(suff,','.join(frzargs)),
            'TrkArg{}="{}"'.format(suff,','.join(trkargs)),
        ]
        print('\n'.join(args))
    else:
        return {'SetArg': setargs, 'FrzArg': frzargs, 'TrkArg': trkargs}

if __name__=="__main__":
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--region", dest="region", type=str, required=True, help="region name")
    parser.add_argument("-n", "--function", dest="function", type=int, required=True, choices=[0,1], help="function (0: main, 1: alt)")
    parser.add_argument("-s", "--suffix", dest="suffix", type=str, required=True, help="suffix for variable names")
    parser.add_argument("-f", "--file", dest="file", type=str, required=True, help="file containing workspace")
    parser.add_argument("-w", "--workspace", dest="workspace", type=str, default="SVJ", help="workspace name")
    args = parser.parse_args()

    main(args.file, args.workspace, args.region, args.function, args.suffix)
