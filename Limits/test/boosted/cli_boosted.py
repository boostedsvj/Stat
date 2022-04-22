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


@is_script
def run_bkg_fit():
    bdt_score = .1
    bdt_str = '{:.1f}'.format(bdt_score).replace('.', 'p')
    mt = get_mt(220, 800, 36)
    with open_root('newmtbinning_dataobs.root') as tf:
        tdir = tf.Get('bsvj_'+bdt_str)
        bkg_hist = tdir.Get('Bkg')
        pdfs, fit_results, data_obs = fit_pdfs_to_histograms(mt, bkg_hist, bkg_hist)
    plot_fits(pdfs[:2], fit_results[:2], data_obs, 'test.pdf')

@is_script
def test_eval_expression():
    assert eval_expression('pow(@0, 2)', [2.]) == 4
    import numpy as np
    np.testing.assert_array_equal(
        eval_expression('pow(@0, 2)', [np.array([2., 4.])]),
        np.array([4., 16.])
        )
    assert add_normalization('@0*@1') == '@2*(@0*@1)'

@is_script
def python_fit():
    from pprint import pprint
    import numpy as np
    bdt_score = .1
    bdt_str = '{:.1f}'.format(bdt_score).replace('.', 'p')
    mt = get_mt(220, 800, 36)
    with open_root('newmtbinning_dataobs.root') as tf:
        tdir = tf.Get('bsvj_'+bdt_str)
        bkg_hist = tdir.Get('Bkg')
        pdf = get_main_pdfs(bkg_hist, mt)[3]
    # res = fit_pdf_expression_to_histogram_python(pdf.expression, bkg_hist)
    # print(res)

    from pprint import pprint

    pprint(dir(pdf.getPdf()))

    print(pdf.getPdf().dependents())

    # argset = pdf.getComponents()

    # fwdit = argset.fwdIterator()
    # for i in range(argset.getSize()):
    #     a = fwdit.next()
    #     print a

    # print a.GetName()
    # pprint(dir(a))

    # print a.getVariables()
    
    # fwdit = a.getVariables().fwdIterator()
    # for i in range(argset.getSize()):
    #     av = fwdit.next()
    #     print av



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
                    print '-'*60
                    print 'Summary of varous fit strategies'
                    print '\nScipy only:'
                    print res_scipy
                    print '\nRooFit with initial parameters from Scipy:'
                    res_roofit_wscipy.Print()
                    print '\nRooFit with initial parameters set to 1.:'
                    res_roofit_only.Print()

        if args.bdtcut is not None:
            tdir_name = 'bsvj_{:.1f}'.format(args.bdtcut).replace('.', 'p')
            do_plot(tdir_name)
        else:
            for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
                do_plot(tdir_name)




    #         if mt_pointer[0] is None:
    #             left = bkg_hist.GetBinLowEdge(1)
    #             right = bkg_hist.GetBinLowEdge(bkg_hist.GetNbinsX()+1)
    #             n_bins = bkg_hist.GetNbinsX()
    #             logger.info('Making mt variable: %s to %s with %s bins', left, right, n_bins)
    #             mt_pointer[0] = get_mt(left, right, n_bins)
    #         mt = mt_pointer[0]
    #         data_obs = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_hist, 1.)

    #         for pdf_set_name, pdfs in [
    #             ('alt', get_alt_pdfs(bkg_hist, mt)),
    #             ('main', get_main_pdfs(bkg_hist, mt)),
    #             ]:
    #             for pdf in pdfs:
    #                 data_th1 = data_obs.createHistogram(str(uuid.uuid4()), mt)
    #                 res_purepy = fit_pdf_expression_to_histogram_python(pdf.expression, data_th1)
    #                 res_purepy = fit_pdf_expression_to_histogram_python(pdf.expression, data_th1, res_purepy.x, tol=1e-6)
    #                 res_pureroot = fit_pdf_to_datahist(rebuild_rpsbp(pdf), data_obs, init_with_python_fit=False)
    #                 res_rootinitpy = fit_pdf_to_datahist(rebuild_rpsbp(pdf), data_obs, init_with_python_fit=True)
    #                 plot_pdf_for_various_fitresults(
    #                     pdf,
    #                     [res_purepy, res_pureroot, res_rootinitpy],
    #                     data_obs,
    #                     osp.join(
    #                         plotdir,
    #                         '{0}_{1}_npar{2}.png'.format(
    #                             tdir.GetName(), pdf_set_name, count_parameters(pdf.expression)
    #                             )
    #                         ),
    #                     labels=['Scipy only', 'RooFit only', 'RooFit init w/ scipy']
    #                     )

    # with open_root(args.rootfile) as tf:
    #     for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
    #         make_fits(tf.Get(tdir_name))



@is_script
def gen_datacard():
    """
    Intending to make full datacard
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-b', '--bdtcut', type=float, default=.1)
    args = parser.parse_args()

    bdt_str = '{:.1f}'.format(args.bdtcut).replace('.', 'p')
    signal_name = lambda mz: 'SVJ_mZprime{}_mDark10_rinv03_alphapeak'.format(mz)
    with open_root(args.rootfile) as tf:
        tdir = tf.Get('bsvj_'+bdt_str)
        bkg_hist = tdir.Get('Bkg')
        mt = get_mt_from_th1(bkg_hist, name='mt')
        data_hist = tdir.Get('data_obs')
        data_obs = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), tdir.Get('data_obs'), 1.)
        signals = [
            ROOT.RooDataHist(
                signal_name(mz), signal_name(mz), ROOT.RooArgList(mt),
                tdir.Get(signal_name(mz)), 1.
                )
            for mz in [250, 300, 350]
            ]

    def fit(pdf):
        res_scipy = fit_scipy(pdf.pdftype, pdf.npars, bkg_hist)
        res_roofit_wscipy = fit_roofit(pdf.pdftype, pdf.npars, bkg_hist, init_vals=res_scipy.x)
        return res_roofit_wscipy

    winner_pdfs = []
    for pdf_type in ['main', 'alt']:
        pdfs = get_pdfs(pdf_type, bkg_hist, mt)
        ress = [ fit(pdf) for pdf in pdfs ]
        i_winner = do_fisher_test(mt, data_obs, pdfs)
        winner_pdfs.append(pdfs[i_winner])
        plot_fits(pdfs, ress, data_obs, pdf_type + '.pdf')

    # print(signals)
    # return

    systs = [['lumi', 'lnN', 1.026, '-']]
    for mz, sig in zip([250, 300, 350], signals):
        compile_datacard_macro(winner_pdfs, data_obs, sig, strftime('dc_%b%d/dc_mz{}_bdt{}.txt'.format(mz, bdt_str)), systs=systs)



@is_script
def test_Feb18():
    bdt_score = .1
    bdt_str = '{:.1f}'.format(bdt_score).replace('.', 'p')
    mt = get_mt(220, 800, 36)
    signal_name = lambda mz: 'SVJ_mZprime{}_mDark10_rinv03_alphapeak'.format(mz)
    with open_root('newmtbinning_dataobs.root') as tf:
        tdir = tf.Get('bsvj_'+bdt_str)
        bkg_hist = tdir.Get('Bkg')
        data_hist = tdir.Get('data_obs')
        data_obs = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), tdir.Get('data_obs'), 1.)
        signals = [
            ROOT.RooDataHist(
                "sig", signal_name(mz), ROOT.RooArgList(mt),
                tdir.Get(signal_name(mz)), 1.
                )
            for mz in [250, 300, 350]
            ]

    pdf = get_alt_pdfs(bkg_hist, mt)[2]

    data_th1 = data_obs.createHistogram(str(uuid.uuid4()), mt)
    res_purepy = fit_pdf_expression_to_histogram_python(pdf.expression, data_th1)
    res_purepy = fit_pdf_expression_to_histogram_python(pdf.expression, data_th1, res_purepy.x, tol=1e-6)
    res_pureroot = fit_pdf_to_datahist(rebuild_rpsbp(pdf), data_obs, init_with_python_fit=False)
    res_rootinitpy = fit_pdf_to_datahist(rebuild_rpsbp(pdf), data_obs, init_with_python_fit=True)

    plot_pdf_for_various_fitresults(
        pdf,
        # [res_purepy],
        [res_purepy, res_pureroot, res_rootinitpy],
        data_obs, 'multires.pdf'
        )
    
    print(res_purepy)
    res_pureroot.Print()
    res_rootinitpy.Print()

    return

    # pdf2 = get_alt_pdfs(bkg_hist, mt)[0]
    # res_without_python = fit_pdf_to_datahist(pdf2, data_obs, init_with_python_fit=False)
    # plot_fits([pdf2], [res_without_python], data_obs, strftime('main_pdf_fit_wopy_%b%d.pdf'))
    # pdf1 = get_alt_pdfs(bkg_hist, mt)[0]
    # res_with_python = fit_pdf_to_datahist(pdf1, data_obs, init_with_python_fit=True, unc_multiplier=100.)
    # plot_fits([pdf1], [res_with_python], data_obs, strftime('main_pdf_fit_wipy_%b%d.pdf'))
    # return

    pdfs = get_alt_pdfs(bkg_hist, mt)
    fit_results_nopy = [
        fit_pdf_to_datahist(pdf, data_obs, init_with_python_fit=False) for pdf in pdfs
        ]
    plot_fits(pdfs, fit_results_nopy, data_obs, strftime('main_pdf_fit_nopy_%b%d.pdf'))
    # return

    pdfs = get_alt_pdfs(bkg_hist, mt)
    fit_results_wipy = [
        fit_pdf_to_datahist(pdf, data_obs, init_with_python_fit=True) for pdf in pdfs
        ]
    plot_fits(pdfs, fit_results_wipy, data_obs, strftime('main_pdf_fit_wipy_%b%d.pdf'))

    print('No py results:')
    for res in fit_results_nopy: res.Print()
    print('With py results:')
    for res in fit_results_wipy: res.Print()

    return


    pdfs = get_main_pdfs(bkg_hist, mt)
    unc_multipliers = [ 10. for i in range(len(pdfs)) ]
    unc_multipliers[2] = 20. # PDF 3 needs a little wider ranges to converge
    fit_results = [
        fit_pdf_to_datahist(pdf, data_obs, unc_multiplier=unc_multiplier) \
        for pdf, unc_multiplier in zip(pdfs, unc_multipliers)
        ]
    for res in fit_results: res.Print()
    
    plot_fits(pdfs, fit_results, data_obs, 'test.pdf')
    n_fit_parameters_per_pdf = [ pdf.getVariables().getSize()-1 for pdf in pdfs ]
    winner = do_fisher_test(mt, data_obs, pdfs, n_fit_parameters_per_pdf)

    systs = [
        ['lumi', 'lnN', 1.026, '-']
        ]

    for mz, sig in zip([250, 300, 350], signals):
        compile_datacard_macro(pdfs[winner], data_obs, sig, strftime('dc_%b%d/dc_mz{}_bdt{}.txt'.format(mz, bdt_str)), systs=systs)



# ______________________________________________________________
# Development on Sara's histograms

def reimplemented_getCard():
    verbose = False
    sig = 'SVJ_mZprime{}_mDark{}_rinv{}_alpha{}'.format('250', '20', '03', 'peak')
    ch = 'bsvj'

    with open_root('tmva_MThisto_August6_systematicup.root') as ifile:
        histBkgData = ifile.Get('bsvj/Bkg')
        histData = ifile.Get('bsvj/data_obs')
        histSig = ifile.Get('bsvj/SVJ_mZprime250_mDark20_rinv03_alphapeak')
        histBkgData.SetDirectory(0)
        histData.SetDirectory(0)
        histSig.SetDirectory(0)

        mt_min = 160.
        mt_max = 500.
        mT = ROOT.RooRealVar('mH'+ch, 'm_{T}', mt_min, mt_max, 'GeV')
        mT.setBins(30)
        binMin = histData.FindBin(mt_min)
        binMax = histData.FindBin(mt_max)
        data_obs = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mT), histData, 1.)
        sig = ROOT.RooDataHist("sig", "Signal", ROOT.RooArgList(mT), histSig, 1.)
        print "Total data integral: ", histData.Integral()
        nBkgEvts = histBkgData.Integral(binMin, binMax)
        nDataEvts = histData.Integral(binMin, binMax)
        print 'Ranged data integral: ', nDataEvts

        normBkg = ROOT.RooRealVar("Bkg_"+ch+"_norm", "Number of background events", nBkgEvts, 0., 2.e4)
        normData = ROOT.RooRealVar("Data_"+ch+"_norm", "Number of background events", nDataEvts, 0., 2.e4)

        pdfs, keep = get_pdfs(histBkgData, mT, prefix=ch)

        fit_results = []
        for i, rpsbp in enumerate(pdfs):
            print '\n' + '-'*70 + '\nFitting pdf {}'.format(i+1)
            print(data_obs)
            print(rpsbp)
            print(mT)
            try:
                res = rpsbp.fitTo(
                    data_obs,
                    ROOT.RooFit.Extended(False),
                    ROOT.RooFit.Save(1),
                    ROOT.RooFit.SumW2Error(True),
                    ROOT.RooFit.Strategy(2),
                    ROOT.RooFit.Minimizer("Minuit2"),
                    ROOT.RooFit.PrintLevel(-1 if not verbose else 2),
                    ROOT.RooFit.Range('Full')
                    )
            except:
                print 'Problem fitting pdf {}'.format(i+1)
                raise
            fit_results.append(res)

        # Print out fit results
        for i, res in enumerate(fit_results):
            print '\nResult of fit {}'.format(i+1)
            res.Print()

        # Create output workspace and return
        ws_name = 'FitWS'
        ws = ROOT.RooWorkspace(ws_name, 'workspace')
        commit = getattr(ws, 'import')
        commit(data_obs)
        commit(sig)
        for pdf in pdfs: commit(pdf)

        dump_fits_to_file('fitResults_bsvj.root', fit_results)
        dump_ws_to_file('ws_allFits_bsvj.root', ws)

        return ws

def test_plot_fits():
    with open_root('ws_allFits_bsvj.root') as ws_file:
        ws = ws_file.Get('FitWS')
    pdfs = [ ws.pdf('Bkg_bsvj'+str(i+1)) for i in range(4) ]
    data_obs = ws.data('data_obs')
    with open_root('fitResults_bsvj.root') as fit_file:
        fit_results = [ fit_file.Get('fitresult_{}_data_obs'.format(pdf.GetName())) for pdf in pdfs ]
    plot_fits(pdfs, fit_results, data_obs)



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
        print '{thing:8s}: manual={manual:9.3f}, viaframe={viaframe:9.3f}'.format(**locals())

    print 'viaframe took {}s, manual {}s'.format((t1-t0), (t2-t1))


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
    print 'rss_manual: {rss_manual:.2f}; rss_viaframe: {rss_viaframe:.2f}'.format(**locals())



def test_fisher():
    with open_root('ws_allFits_bsvj.root') as ws_file:
        ws = ws_file.Get('FitWS')
    pdfs = [ ws.pdf('Bkg_bsvj'+str(i+1)) for i in range(4) ]
    data_obs = ws.data('data_obs')
    sig = ws.data('sig')

    n_fit_parameters_per_pdf = [ pdf.getVariables().getSize()-1 for pdf in pdfs ]

    mt_min = 160.
    mt_max = 504.
    mt = ROOT.RooRealVar("mHbsvj", "m_{T}", mt_min, mt_max, "GeV")
    mt.setBins(43) # Exactly 8 GeV

    winner = do_fisher_test(mt, data_obs, pdfs, n_fit_parameters_per_pdf)

    compile_datacard_macro(pdfs[winner], data_obs, sig)


if __name__ == '__main__':
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument('script', type=str, choices=list(_scripts.keys()))
    parser.add_argument('-v', '--verbose', action='store_true')
    args, other_args = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + other_args
    if args.verbose: debug()
    _scripts[args.script]()