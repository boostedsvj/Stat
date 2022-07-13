"""
Scripts using building blocks in boosted_fits.py to create datacards
"""

import argparse, inspect, os, os.path as osp
from boosted_fits import *
from time import strftime
from copy import copy

import ROOT
ROOT.RooMsgService.instance().setSilentMode(True)

_scripts = {}
def is_script(fn):
    _scripts[fn.__name__] = fn
    return fn

# _______________________________________________________________________



def fit_scipy(pdf_type, npars, histogram):
    return fit_expr_to_histogram_robust(pdf_expression(pdf_type, npars), histogram)

def fit_roofit(pdf_type, npars, histogram, init_vals=None):
    mt = get_mt_from_th1(histogram)
    pdf = make_pdf(pdf_type, npars, histogram, mt, name=uid())
    data_obs = th1_to_datahist(histogram, mt)
    return fit_pdf_to_datahist(pdf, data_obs, init_vals=init_vals)


@is_script
def plot_scipy_fits():
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-o', '--plotdir', type=str, default='plots_bkgfits_%b%d')
    parser.add_argument('-b', '--bdtcut', type=float, default=None)
    parser.add_argument('-n', '--npars', type=int, nargs='*')
    parser.add_argument('-p', '--pdftype', type=str, default=None, choices=['main', 'alt'])

    args = parser.parse_args()
    plotdir = strftime(args.plotdir)
    if not osp.isdir(plotdir): os.makedirs(plotdir)

    import matplotlib.pyplot as plt
    mpl_fontsizes()

    with open_root(args.rootfile) as tf:

        def do_plot(tdir_name):
            tdir = tf.Get(tdir_name)
            bkg_hist = tdir.Get('Bkg')

            for pdf_type in ['main', 'alt']:
                if args.pdftype and pdf_type != args.pdftype: continue
                logger.info('Fitting pdf_type=%s, tdir_name=%s', pdf_type, tdir_name)
                fig = plt.figure(figsize=(8,8))
                ax = fig.gca()
                binning, counts = th1_binning_and_values(tdir.Get('Bkg'))
                bin_centers = np.array([.5*(l+r) for l, r in zip(binning[:-1], binning[1:])])
                # Bkg histogram
                ax.step(binning[:-1], counts, where='post', label='bkg {}'.format(tdir_name))
                # Fits
                if args.npars is not None and len(args.npars):
                    npars_iter = list(args.npars)
                else:
                    npars_iter = list(range(1,5) if pdf_type == 'alt' else range(2,6))
                for npars in npars_iter:
                    logger.info('Fitting pdf_type=%s, tdir_name=%s, npars=%s', pdf_type, tdir_name, npars)
                    res = fit_scipy(pdf_type, npars, bkg_hist)
                    y_pdf = eval_expression(pdf_expression(pdf_type, npars), [bin_centers] + list(res.x))
                    y_pdf = y_pdf/y_pdf.sum() * counts.sum()
                    chi2 = ((y_pdf-counts)**2 / y_pdf).sum() / (len(bin_centers) - npars)
                    label = '{}_npars{}, chi2={:.3f}, {}'.format(
                        pdf_type, npars, chi2,
                        ', '.join(['p{}={:.3f}'.format(i, v) for i, v in enumerate(res.x)])
                        )
                    ax.plot(bin_centers, y_pdf, label=label)
                ax.legend()
                ax.set_xlabel(r'$m_{T}$ (GeV)')
                ax.set_ylabel(r'$N_{events}$')
                ax.set_yscale('log')
                plt.savefig(osp.join(plotdir, tdir_name + '_' + pdf_type + '.png'), bbox_inches='tight')        

        bdtcut = None
        if args.bdtcut is not None:
            tdir_name = 'bsvj_{:.1f}'.format(args.bdtcut).replace('.', 'p')
            do_plot(tdir_name)
        else:
            for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
                do_plot(tdir_name)


@is_script
def plot_roofit_fits():
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-o', '--plotdir', type=str, default='plots_bkgfits_%b%d')
    parser.add_argument('-b', '--bdtcut', type=float, default=None)
    parser.add_argument('-n', '--npars', type=int, nargs='*')
    parser.add_argument('-p', '--pdftype', type=str, default=None, choices=['main', 'alt'])
    args = parser.parse_args()
    plotdir = strftime(args.plotdir)
    if not osp.isdir(plotdir): os.makedirs(plotdir)

    with open_root(args.rootfile) as tf:

        def do_plot(tdir_name):
            tdir = tf.Get(tdir_name)
            bkg_hist = tdir.Get('Bkg')
            mt = get_mt_from_th1(bkg_hist)
            for pdf_type in ['main', 'alt']:
                if args.pdftype and pdf_type != args.pdftype: continue
                logger.info('Fitting pdf_type=%s, tdir_name=%s', pdf_type, tdir_name)
                if args.npars is not None and len(args.npars):
                    npars_iter = list(args.npars)
                else:
                    npars_iter = list(range(1,5) if pdf_type == 'alt' else range(2,6))
                for npars in npars_iter:
                    logger.info('Fitting pdf_type=%s, tdir_name=%s, npars=%s', pdf_type, tdir_name, npars)
                    res_scipy = fit_scipy(pdf_type, npars, bkg_hist)
                    if len(res_scipy.x) != npars:
                        raise Exception(
                            'Wrong number of fitted parameters.'
                            ' Found {} parameters in scipy fit result, but npars is {}.'
                            ' Scipy fit result:\n{}'
                            .format(len(res_scipy.x), npars, res_scipy)
                            )
                    res_roofit_only = fit_roofit(pdf_type, npars, bkg_hist)
                    res_roofit_wscipy = fit_roofit(pdf_type, npars, bkg_hist, init_vals=res_scipy.x)
                    plot_pdf_for_various_fitresults(
                        make_pdf(pdf_type, npars, bkg_hist, mt=mt, name=uid()),
                        [res_scipy, res_roofit_only, res_roofit_wscipy],
                        th1_to_datahist(bkg_hist, mt=mt),
                        osp.join(plotdir, '{0}_{1}_npar{2}.png'.format(tdir.GetName(), pdf_type, npars)),
                        labels=['Scipy only', 'RooFit only', 'RooFit init w/ scipy']
                        )
                    print('-'*60)
                    print('Summary of varous fit strategies')
                    print('\nScipy only:')
                    print(res_scipy)
                    print('\nRooFit with initial parameters from Scipy:')
                    res_roofit_wscipy.Print()
                    print('\nRooFit with initial parameters set to 1.:')
                    res_roofit_only.Print()

        if args.bdtcut is not None:
            tdir_name = 'bsvj_{:.1f}'.format(args.bdtcut).replace('.', 'p')
            do_plot(tdir_name)
        else:
            for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
                do_plot(tdir_name)



def this_fn_name():
    """
    Returns the name of whatever function this function was called from.
    (inspect.stack()[0][3] would be "this_fn_name"; [3] is just the index for names)
    """
    return inspect.stack()[1][3]


@is_script
def gen_datacard_ul():
    """
    Generate datacards for UL
    """
    parser = argparse.ArgumentParser(this_fn_name())
    parser.add_argument('jsonfile', type=str)
    parser.add_argument('-b', '--bdtcut', type=float, default=.1)
    parser.add_argument('-i', '--injectsignal', action='store_true')
    args = parser.parse_args()
    bdt_str = '{:.1f}'.format(args.bdtcut).replace('.', 'p')

    import json
    with open(args.jsonfile, 'r') as f:
        d = json.load(f)

    mt_binning = d['mt']

    # Cut off left part of mt distr: Find first bin >300 GeV
    for i in range(len(mt_binning)-1):
        if mt_binning[i] >= 300.:
            i_bin_min = i
            break
    else:
        raise Exception()

    logger.info('Starting mt from {}'.format(mt_binning[i_bin_min]))
    mt_binning = mt_binning[i_bin_min:]
    n_bins = len(mt_binning)-1
    mt = get_mt(mt_binning[0], mt_binning[-1], n_bins, name='mt')

    def construct_th1_from_hist(name, hist):
        th1 = ROOT.TH1F(name, name, n_bins, array('f', mt_binning))
        ROOT.SetOwnership(th1, False)
        vals = hist['vals'][i_bin_min:]
        errs = hist['errs'][i_bin_min:]
        assert len(vals) == n_bins
        assert len(errs) == n_bins
        for i in range(n_bins):
            th1.SetBinContent(i+1, vals[i])
            th1.SetBinError(i+1, errs[i])
        return th1

    # Construct the bkg TH1
    bkg_hist = d['histograms']['{:.1f}/bkg'.format(args.bdtcut)]
    bkg_th1 = construct_th1_from_hist('bkg', bkg_hist)

    # data RooDataHist: just the bkg now
    data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_th1, 1.)

    def fit(pdf):
        res_scipy = fit_scipy(pdf.pdftype, pdf.npars, bkg_th1)
        res_roofit_wscipy = fit_roofit(pdf.pdftype, pdf.npars, bkg_th1, init_vals=res_scipy.x)
        return res_roofit_wscipy

    winner_pdfs = []
    for pdf_type in ['main', 'alt']:
        pdfs = get_pdfs(pdf_type, bkg_th1, mt)
        ress = [ fit(pdf) for pdf in pdfs ]
        i_winner = do_fisher_test(mt, data_datahist, pdfs)
        winner_pdfs.append(pdfs[i_winner])
        plot_fits(pdfs, ress, data_datahist, pdf_type + '.pdf')

    systs = [['lumi', 'lnN', 1.026, '-']]

    for mz in [450]:
        sig_name = 'mz{}'.format(mz)
        sig_hist = d['histograms']['{:.1f}/{}'.format(args.bdtcut, sig_name)]
        sig_th1 = construct_th1_from_hist(sig_name, sig_hist)
        sig_datahist = ROOT.RooDataHist(sig_name, sig_name, ROOT.RooArgList(mt), sig_th1, 1.)

        if args.injectsignal:
            logger.info('Injecting signal in data_obs')
            data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_th1+sig_th1, 1.)

        compile_datacard_macro(
            winner_pdfs, data_datahist, sig_datahist,
            strftime('dc_%b%d/dc_mz{}_bdt{}.txt'.format(mz, bdt_str)),
            systs=systs
            )


@is_script
def simple_test_fit():
    """
    Runs a simple AsymptoticLimits fit on a datacard, without many options
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('datacard', type=str)
    parser.add_argument('-c', '--chdir', type=str, default=None)
    args = parser.parse_args()

    cmd = CombineCommand(args.datacard, 'AsymptoticLimits')
    cmd.track_parameters.extend(['r'])
    cmd.args.add('--saveWorkspace')
    cmd.set_parameter('pdf_index', 1)
    cmd.freeze_parameters.extend([
        'pdf_index',
        'bsvj_bkgfitmain_npars4_p1', 'bsvj_bkgfitmain_npars4_p2', 'bsvj_bkgfitmain_npars4_p3',
        'bsvj_bkgfitmain_npars4_p4',
        # 'bsvj_bkgfitalt_npars3_p1', 'bsvj_bkgfitalt_npars3_p2', 'bsvj_bkgfitalt_npars3_p3'
        ])
    run_combine_command(cmd, args.chdir)


@is_script
def multidimfit():
    """
    Runs a single MultiDimFit on a datacard
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('datacard', type=str)
    parser.add_argument('-c', '--chdir', type=str, default=None)
    parser.add_argument('-a', '--asimov', action='store_true')
    args = parser.parse_args()

    cmd = CombineCommand(args.datacard, 'MultiDimFit')
    cmd.track_parameters.extend(['r'])
    cmd.args.add('--saveWorkspace')
    cmd.args.add('--saveNLL')
    if args.asimov:
        cmd.kwargs['-t'] = '-1'
        cmd.args.add('--toysFreq')
    cmd.set_parameter('pdf_index', 1)
    cmd.freeze_parameters.extend(['pdf_index', r'rgx{bsvj_bkgfitmain_npars.*}'])

    cmd.redefine_signal_pois.append('r')
    cmd.kwargs['--X-rtd'] = 'REMOVE_CONSTANT_ZERO_POINT=1'
    cmd.track_parameters.extend(['r'])

    run_combine_command(cmd, args.chdir)


@is_script
def likelihood_scan():
    """
    Runs a likelihood scan on a datacard
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('datacard', type=str)
    parser.add_argument('-c', '--chdir', type=str, default=None)
    parser.add_argument('-a', '--asimov', action='store_true')
    parser.add_argument('--injectsignal', action='store_true')
    parser.add_argument('-n', '--npoints', type=int, default=51)
    parser.add_argument('-r', '--range', type=float, default=[-.7, .7], nargs=2)
    parser.add_argument('-v', '--verbosity', type=int, default=0)
    args = parser.parse_args()

    cmd = CombineCommand(args.datacard, 'MultiDimFit')

    cmd.redefine_signal_pois.append('r')
    if args.injectsignal: cmd.kwargs['--expectSignal'] = 1
    cmd.add_range('r', args.range[0], args.range[1])
    cmd.track_parameters.extend(['r'])

    cmd.args.add('--saveWorkspace')
    cmd.args.add('--saveNLL')
    cmd.kwargs['--algo'] = 'grid'
    cmd.kwargs['--points'] = args.npoints
    cmd.kwargs['--X-rtd'] = 'REMOVE_CONSTANT_ZERO_POINT=1'
    cmd.kwargs['--alignEdges'] = 1
    cmd.kwargs['-v'] = args.verbosity

    if args.asimov:
        cmd.kwargs['-t'] = '-1'
        cmd.args.add('--toysFreq')
        cmd.kwargs['-n'] = 'Asimov'
    else:
        cmd.kwargs['-n'] = 'Observed'

    if args.injectsignal: cmd.kwargs['-n'] += 'InjectedSig'

    cmd.set_parameter('pdf_index', 1)
    cmd.freeze_parameters.extend([
        'pdf_index',
        # r'rgx{bsvj_bkgfitmain_npars.*}'
        'bsvj_bkgfitmain_npars4_p1', 'bsvj_bkgfitmain_npars4_p2', 'bsvj_bkgfitmain_npars4_p3',
        'bsvj_bkgfitmain_npars4_p4',
        # 'bsvj_bkgfitalt_npars3_p1', 'bsvj_bkgfitalt_npars3_p2', 'bsvj_bkgfitalt_npars3_p3'
        ])

    run_combine_command(cmd, args.chdir)


@is_script
def printws():
    """
    Prints a workspace contents
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-w', '--workspace', type=str)
    args = parser.parse_args()
    with open_root(args.rootfile) as f:
        ws = get_ws(f, args.workspace)
    ws.Print()
    return ws



def test_chi2():
    import time

    with open_root('ws_allFits_bsvj.root') as ws_file:
        ws = ws_file.Get('FitWS')
    pdfs = [ ws.pdf('Bkg_bsvj'+str(i+1)) for i in range(4) ]
    data_obs = ws.data('data_obs')

    with open_root('fitResults_bsvj.root') as fit_file:
        fit_results = [ fit_file.Get('fitresult_{}_data_obs'.format(pdf.GetName())) for pdf in pdfs ]

    mt_min = 160.
    mt_max = 504.
    mt = ROOT.RooRealVar("mHbsvj", "m_{T}", mt_min, mt_max, "GeV")
    mt.setBins(43) # Exactly 8 GeV

    pdf = pdfs[0]
    fit_result = fit_results[0]
    n_fit_parameters = fit_result.floatParsFinal().getSize()
    
    t0 = time.time()
    chi2_viaframe = get_chi2_viaframe(mt, pdf, data_obs, n_fit_parameters, verbose=True)
    t1 = time.time()
    chi2_manual = get_chi2(mt, pdf, data_obs, n_fit_parameters, verbose=True)
    t2 = time.time()
    things = [ 'chi2/ndf', 'chi2', 'prob', 'ndf' ]

    for thing, manual, viaframe in zip(things, chi2_manual, chi2_viaframe):
        print('{thing:8s}: manual={manual:9.3f}, viaframe={viaframe:9.3f}'.format(**locals()))

    print('viaframe took {}s, manual {}s'.format((t1-t0), (t2-t1)))


def test_rss():
    with open_root('ws_allFits_bsvj.root') as ws_file:
        ws = ws_file.Get('FitWS')
    pdfs = [ ws.pdf('Bkg_bsvj'+str(i+1)) for i in range(4) ]
    data_obs = ws.data('data_obs')

    mt_min = 160.
    mt_max = 504.
    mt = ROOT.RooRealVar("mHbsvj", "m_{T}", mt_min, mt_max, "GeV")
    mt.setBins(43) # Exactly 8 GeV

    pdf = pdfs[0]
    rss_manual = get_rss(mt, pdf, data_obs, verbose=True)
    rss_viaframe = get_rss_viaframe(mt, pdf, data_obs, verbose=True)
    print('rss_manual: {rss_manual:.2f}; rss_viaframe: {rss_viaframe:.2f}'.format(**locals()))



if __name__ == '__main__':
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument('script', type=str, choices=list(_scripts.keys()))
    parser.add_argument('-v', '--verbose', action='store_true')
    args, other_args = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + other_args
    if args.verbose: debug()
    r = _scripts[args.script]()