import ROOT
import boosted_fits as bsvj
import numpy as np
import inspect

TEST_ROOTFILE = 'example_fit_result.root'

def get_ws(rootfile=None):
    if rootfile is None: rootfile = TEST_ROOTFILE
    with bsvj.open_root(rootfile) as f:
        return bsvj.get_ws(f)


def test_roofit_get_y_values():
    # Values from text datacard used to create the example ws
    n_expected_data = 825105
    n_expected_sig = 24666

    ws = get_ws()
    ws.loadSnapshot('MultiDimFit')

    mt = ws.var('mt')
    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[1:]+mt_binning[:-1])

    # Test getting from RooDataSet
    data = ws.data('data_obs')
    x_data, y_data = bsvj.roodataset_values(data)

    assert abs(y_data.sum() - n_expected_data) <= 1.
    np.testing.assert_array_equal(mt_bin_centers, x_data)
    assert len(mt_binning)-1 == len(y_data)

    # Compare with values from createHistogram
    data_th1 = ROOT.RooDataHist(bsvj.uid(), '', ROOT.RooArgSet(mt), data).createHistogram(bsvj.uid(), mt)
    data_th1_binning, y_data_th1, errs_data_th1 = bsvj.th1_binning_and_values(data_th1, True)

    np.testing.assert_array_equal(mt_binning, data_th1_binning)
    np.testing.assert_array_equal(y_data, y_data_th1)
    np.testing.assert_array_equal(y_data, errs_data_th1)


    # Test getting from RooMultiPdf/RooAbsPdf
    bkg = ws.pdf('shapeBkg_roomultipdf_bsvj')
    y_bkg = bsvj.pdf_values(bkg, x_data)
    bkg_norm = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()

    # Compare with values from createHistogram
    bkg_th1 = bkg.createHistogram(bsvj.uid(), mt)
    bkg_th1_binning, y_bkg_th1 = bsvj.th1_binning_and_values(bkg_th1)

    assert len(y_bkg) == len(y_data)
    np.testing.assert_almost_equal(bkg_th1_binning, mt_binning)
    np.testing.assert_almost_equal(y_bkg, y_bkg_th1)

    # Test getting from RooDataHist
    sig = ws.embeddedData('shapeSig_sig_bsvj')
    x_sig, y_sig = bsvj.roodataset_values(sig)

    np.testing.assert_array_equal(mt_bin_centers, x_sig)
    assert abs(y_sig.sum() - n_expected_sig) <= 1.

    bsvj.logger.info(inspect.stack()[0][3] + ' finished succesfully')


def test_eval_expression():
    assert bsvj.eval_expression('pow(@0, 2)', [2.]) == 4
    import numpy as np
    np.testing.assert_array_equal(
        bsvj.eval_expression('pow(@0, 2)', [np.array([2., 4.])]),
        np.array([4., 16.])
        )
    assert bsvj.add_normalization('@0*@1') == '@2*(@0*@1)'


def test_chi2():
    import time
    from scipy import stats

    ws = get_ws()

    mt = ws.var('mt')
    pdf = ws.pdf('bsvj_bkgfitalt_npars1_rpsbp')
    data = ws.data('data_obs')

    res = pdf.fitTo(data, ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Save(True))
    bsvj.logger.info('Fit result: %s', res)

    n_fit_parameters = res.floatParsFinal().getSize()
    bsvj.logger.info('n_fit_parameters: %s', n_fit_parameters)

    t0 = time.time()
    chi2_viaframe = bsvj.get_chi2_viaframe(mt, pdf, data, n_fit_parameters)
    t1 = time.time()
    rss_viaframe = bsvj.get_rss_viaframe(mt, pdf, data)
    t2 = time.time()

    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[:-1] + mt_binning[1:])
    _, y_data = bsvj.roodataset_values(data)
    y_pdf = bsvj.pdf_values(pdf, mt_bin_centers) * y_data.sum()
    raw_chi2 = ((y_pdf-y_data)**2 / y_pdf).sum()
    ndf = len(mt_bin_centers) - n_fit_parameters
    chi2 = raw_chi2 / ndf
    prob = stats.distributions.chi2.sf(raw_chi2, ndf) # or cdf?
    rss = np.sqrt(((y_pdf-y_data)**2).sum())
    t3 = time.time()

    bsvj.logger.info('chi2_viaframe: %s, took %s sec', chi2_viaframe, t1-t0)
    bsvj.logger.info(
        'chi2_manual: %s, took %s sec',
        (chi2, raw_chi2, prob, ndf),
        t3-t2
        )

    # ________________________________________________________
    # Extra test: Fit normalization manually to data again

    def fom(norm):
        y = y_pdf * norm
        return ((y-y_data)**2 / y).sum() / ndf
    from scipy.optimize import minimize
    res = minimize(fom, 1.)
    
    y_pdf = res.x * y_pdf
    raw_chi2 = ((y_pdf-y_data)**2 / y_pdf).sum()
    ndf = len(mt_bin_centers) - n_fit_parameters
    chi2 = raw_chi2 / ndf
    prob = stats.distributions.chi2.sf(raw_chi2, ndf) # or cdf?
    bsvj.logger.info('chi2_manual after fitting norm: %s',(chi2, raw_chi2, prob, ndf))

    # ________________________________________________________
    # RSS

    bsvj.logger.info('rss_viaframe: %s, took %s sec', rss_viaframe, t2-t1)
    bsvj.logger.info('rss_manual: %s', rss)



def test_combine_command():
    cmd = bsvj.CombineCommand('bla.txt')
    cmd.kwargs['--test'] = .1
    cmd.kwargs['--test2'] = 'NO'
    cmd.args.add('--saveNLL')
    cmd.freeze_parameters.extend(['x', 'y'])
    assert cmd.parameters == {}
    assert cmd.freeze_parameters == ['x', 'y']
    assert cmd.dc == 'bla.txt'
    assert cmd.str == 'combine -M MultiDimFit bla.txt --saveNLL --test 0.1 --test2 NO --freezeParameters x,y'
    cmd.add_range('x', .4, 1.6)
    cmd.add_range('y', 100, 101)
    assert cmd.str == (
        'combine -M MultiDimFit bla.txt --saveNLL --test 0.1 --test2 NO --freezeParameters x,y'
        ' --setParameterRanges x=0.4,1.6:y=100,101'
        )


if __name__ == '__main__':
    # test_roofit_get_y_values()
    # test_eval_expression()
    # test_chi2()
    test_combine_command()