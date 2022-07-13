import ROOT
import boosted_fits as bsvj
import numpy as np
import inspect

TEST_ROOTFILE = 'example_fit_result.root'

def get_ws(rootfile=None):
    if rootfile is None: rootfile = TEST_ROOTFILE
    with bsvj.open_root(rootfile) as f:
        return bsvj.get_ws(f)

def test_get_y_values():
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



if __name__ == '__main__':
    test_get_y_values()