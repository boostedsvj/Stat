"""
Scripts using building blocks in boosted_fits.py to create datacards
"""


from boosted_fits import *
from time import strftime


_scripts = {}
def is_script(fn):
    _scripts[fn.__name__] = fn
    return fn

@is_script
def test_Feb18():
    mt = get_mt(220, 800, 36)
    with open_root('newmtbinning_dataobs.root') as tf:
        tdir = tf.Get('bsvj_0p0')
        bkg_hist = tdir.Get('Bkg')
        pdfs, fit_results, data_obs, _ = fit_pdfs_to_histograms(mt, bkg_hist, bkg_hist)
        signal_name = lambda mz: 'SVJ_mZprime{}_mDark10_rinv03_alphapeak'.format(mz)
        signals = [
            ROOT.RooDataHist(
                "sig", signal_name(mz), ROOT.RooArgList(mt),
                tdir.Get(signal_name(mz)), 1.
                )
            for mz in [250, 300, 350]
            ]
    plot_fits(pdfs[:2], fit_results[:2], data_obs, 'test.pdf')
    n_fit_parameters_per_pdf = [ pdf.getVariables().getSize()-1 for pdf in pdfs ]
    winner = do_fisher_test(mt, data_obs, pdfs, n_fit_parameters_per_pdf)
    for mz, sig in zip([250, 300, 350], signals):
        compile_datacard_macro(pdfs[winner], data_obs, sig, strftime('dc_%b%d/dc_mz{}.txt'.format(mz)))



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
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('script', type=str, choices=list(_scripts.keys()))
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    if args.verbose: debug()
    _scripts[args.script]()