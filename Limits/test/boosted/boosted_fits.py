"""
Building blocks to create the boosted SVJ analysis datacards
"""

import ROOT
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetPadColor(0)

from contextlib import contextmanager
from array import array
from math import sqrt
import itertools, re, logging, os, os.path as osp
from collections import OrderedDict


DEFAULT_LOGGING_LEVEL = logging.INFO

def setup_logger(name='boosted'):
    if name in logging.Logger.manager.loggerDict:
        logger = logging.getLogger(name)
        logger.info('Logger %s is already defined', name)
    else:
        fmt = logging.Formatter(
            fmt = '\033[33m%(levelname)7s:%(asctime)s:%(module)s:%(lineno)s\033[0m %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
            )
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logger = logging.getLogger(name)
        logger.setLevel(DEFAULT_LOGGING_LEVEL)
        logger.addHandler(handler)
    return logger
logger = setup_logger()

def debug(flag=True):
    """Sets the logger level to debug (for True) or warning (for False)"""
    logger.setLevel(logging.DEBUG if flag else DEFAULT_LOGGING_LEVEL)


@contextmanager
def open_root(path, mode='READ'):
    '''Context manager that takes care of closing the ROOT file properly'''
    try:
        tfile = ROOT.TFile.Open(path, mode)
        yield tfile
    finally:
        tfile.Close()

def dump_fits_to_file(filename, results):
    logger.info('Dumping fit results to ' + filename)
    with open_root(filename, 'RECREATE') as tf:
        for result in results: result.Write()

def dump_ws_to_file(filename, ws):
    logger.info('Dumping pdfs ws to ' + filename)
    wstatus = ws.writeToFile(filename, True)
    return wstatus


def get_mt(mt_min=160., mt_max=500., n_bins=43):
    """
    Sensible defaults for the mt axis
    """
    mt = ROOT.RooRealVar('mHbsvj', 'm_{T}', mt_min, mt_max, 'GeV')
    mt.setBins(n_bins)
    # Manually add the boundaries to it as python attributes for easy access
    mt.mt_min = mt_min
    mt.mt_max = mt_max
    return mt


def get_pdfs(bkg_hist, mt, prefix='bsvj'):
    '''
    Creates the bkg pdfs to be fitted. Returns a list of RooParametricShapeBinPdf.
    '''
    # Function from Theorists, combo testing, sequence E, 1, 11, 12, 22
    # model NM has N params on 1-x and M params on x. exponents are (p_i + p_{i+1} * log(x))
    # these are the ROOT.RooGenericPdf verisons, convert to ROOT.RooParametricShapeBinPdf below
    model_name = 'Bkg_' + prefix
    logger.info('Doing get_pdfs for %s', model_name)

    logger.info('Making pdf for fn 1')
    p1_1 = ROOT.RooRealVar(prefix + "_p1_1", "p1", 1., -450., 450.)
    p2_1 = ROOT.RooRealVar(prefix + "_p2_1", "p2", 1., -100., 100.)
    modelBkg1_rgp = ROOT.RooGenericPdf(
        model_name+"1_rgp", "Thry. fit (11)",
        "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2))",
        ROOT.RooArgList(mt, p1_1, p2_1)
        )
    modelBkg1_rpsbp = ROOT.RooParametricShapeBinPdf(
        model_name+"1", "Thry. Fit (11)",
        modelBkg1_rgp, mt, ROOT.RooArgList(p1_1, p2_1), bkg_hist
        )

    logger.info('Making pdf for fn 2')
    p1_2 = ROOT.RooRealVar(prefix + "_p1_2", "p1", 1., -450., 450.)
    p2_2 = ROOT.RooRealVar(prefix + "_p2_2", "p2", 1., -100., 100.)
    p3_2 = ROOT.RooRealVar(prefix + "_p3_2", "p3", 1., -15, 15)
    modelBkg2_rgp = ROOT.RooGenericPdf(
        model_name+"2_rgp", "Thry. fit (12)",
        "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)))",
        ROOT.RooArgList(mt, p1_2, p2_2, p3_2)
        )
    modelBkg2_rpsbp = ROOT.RooParametricShapeBinPdf(
        model_name+"2", "Thry. Fit (12)",
        modelBkg2_rgp, mt, ROOT.RooArgList(p1_2, p2_2, p3_2), bkg_hist
        )

    logger.info('Making pdf for fn 3')
    p1_3 = ROOT.RooRealVar(prefix + "_p1_3", "p1", 1., -450., 450.)
    p2_3 = ROOT.RooRealVar(prefix + "_p2_3", "p2", 1., -100., 50.)
    p3_3 = ROOT.RooRealVar(prefix + "_p3_3", "p3", 1., -50., 20.)
    p4_3 = ROOT.RooRealVar(prefix + "_p4_3", "p4", 1., -20., 20.)
    modelBkg3_rgp = ROOT.RooGenericPdf(
        model_name+"3_rgp","Thry. fit (22)",
        "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2)))",
        # "pow(1 - @0/13000, @1+@2*log(@0/13000)) * pow(@0/13000, -(@3+@4*log(@0/13000)))",
        ROOT.RooArgList(mt, p1_3, p2_3, p3_3, p4_3)
        )
    # modelBkg3_rgp = ROOT.RooGenericPdf(
    #     model_name+"3_rgp", "Thry. fit (13)",
    #     "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2)))",
    #     ROOT.RooArgList(mt, p1_3, p2_3, p3_3, p4_3)
    #     )
    modelBkg3_rpsbp = ROOT.RooParametricShapeBinPdf(
        model_name+"3", "Thry. Fit (22)",
        modelBkg3_rgp, mt, ROOT.RooArgList(p1_3, p2_3, p3_3, p4_3), bkg_hist
        )

    logger.info('Making pdf for fn 4')
    p1_4 = ROOT.RooRealVar(prefix + "_p1_4", "p1", 1., -450., 450.)
    p2_4 = ROOT.RooRealVar(prefix + "_p2_4", "p2", 1., -15., 15.)
    p3_4 = ROOT.RooRealVar(prefix + "_p3_4", "p3", 1., -450., 450.)
    p4_4 = ROOT.RooRealVar(prefix + "_p4_4", "p4", 1., -15.0, 100.)
    p5_4 = ROOT.RooRealVar(prefix + "_p5_4", "p5", 1., -15., 15.)
    modelBkg4_rgp = ROOT.RooGenericPdf(
        model_name+"4_rgp", "Thry. fit (32)",
        "pow(1 - @0/13000, @1+@2*log(@0/13000)+@3*pow(log(@0/13000),2)) * pow(@0/13000, -(@4+@5*log(@0/13000)))",
        ROOT.RooArgList(mt, p1_4, p2_4, p3_4, p4_4, p5_4)
        )
    # modelBkg4_rgp = ROOT.RooGenericPdf(
    #     model_name+"4_rgp", "Thry. fit (14)",
    #     "pow(1 - @0/13000, @1) * pow(@0/13000, -(@2+@3*log(@0/13000)+@4*pow(log(@0/13000),2)+@5*pow(log(@0/13000),3)))",
    #     ROOT.RooArgList(mt, p1_4, p2_4, p3_4, p4_4, p5_4)
    #     )
    # modelBkg4_rgp = ROOT.RooGenericPdf(
    #     model_name+"4_rgp", "Thry. fit (41)",
    #     "pow(1 - @0/13000, @1+@2*log(@0/13000)+@3*pow(log(@0/13000),2)+@4*pow(log(@0/13000),3)) * pow(@0/13000, -@5)",
    #     ROOT.RooArgList(mt, p1_4, p2_4, p3_4, p4_4, p5_4)
    #     )
    modelBkg4_rpsbp = ROOT.RooParametricShapeBinPdf(
        model_name+"4", "Thry. Fit (32)",
        modelBkg4_rgp, mt, ROOT.RooArgList(p1_4, p2_4, p3_4, p4_4, p5_4), bkg_hist
        )

    pdfs = [modelBkg1_rpsbp, modelBkg2_rpsbp, modelBkg3_rpsbp, modelBkg4_rpsbp]
    # Return a list of needed input parameters as well, to prevent ROOT from garbage
    # collecting them
    keep = [
        p1_1, p2_1, p1_2, p2_2, p3_2, p1_3, p2_3, p3_3, p4_3, p1_4, p2_4, p3_4, p4_4, p5_4,
        modelBkg1_rgp, modelBkg2_rgp, modelBkg3_rgp, modelBkg4_rgp
        ]
    logger.info('Created the following (unfitted) pdfs: %s', pdfs)
    return pdfs, keep



def fit_pdfs_to_histograms(mt, bkg_hist, data_hist):
    """
    mt: RooRealVar representing the mT axis
    bkg_hist: TH1 of the MC bkg
    data_hist: TH1 of the data (replaced by the MC bkg during development)
    """
    logger.info('Fitting pdfs to data_hist (%s) in variable mT (%s)', repr(data_hist), repr(mt))
    logger.info("Total data integral: %s", data_hist.Integral())

    # Get (at this point unfitted) pdfs
    pdfs, keep = get_pdfs(bkg_hist, mt)

    data_obs = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), data_hist, 1.)

    i_bin_min = data_hist.FindBin(mt.mt_min)
    i_bin_max = data_hist.FindBin(mt.mt_max)
    n_events_bkg = bkg_hist.Integral(i_bin_min, i_bin_max)
    n_events_data = data_hist.Integral(i_bin_min, i_bin_max)

    logger.info('data_hist {}: total events: {}; events from {} to {}: {}'.format(
        data_hist.GetName(), data_hist.Integral(), mt.mt_min, mt.mt_max, n_events_data
        ))

    # norm_bkg = ROOT.RooRealVar("Bkg_bsvj_norm", "Number of background events", n_events_bkg, 0., 2.e4)
    # norm_data = ROOT.RooRealVar("Data_bsvj_norm", "Number of background events", n_events_data, 0., 2.e4)

    # Fit the pdfs to data
    fit_results = []
    for i, rpsbp in enumerate(pdfs):
        logger.info(
            'Fitting pdf %s %s to data_obs (%s) in variable mT (%s)',
            i, repr(rpsbp), repr(data_obs), repr(mt)
            )
        try:
            res = rpsbp.fitTo(
                data_obs,
                ROOT.RooFit.Extended(False),
                ROOT.RooFit.Save(1),
                ROOT.RooFit.SumW2Error(True),
                ROOT.RooFit.Strategy(2),
                ROOT.RooFit.Minimizer("Minuit2"),
                ROOT.RooFit.PrintLevel(2 if logger.level <= logging.DEBUG else -1),
                ROOT.RooFit.Range('Full')
                )
        except:
            logger.error('Problem fitting pdf {}'.format(i+1))
            raise
        fit_results.append(res)

    # Print out fit results
    for i, res in enumerate(fit_results):
        logger.info('Result of fit {}'.format(i+1))
        res.Print()

    return pdfs, fit_results, data_obs, keep


def get_variables(rooabsarg):
    """
    Returns a list of all variables a RooAbsArg depends on
    """
    argset = ROOT.RooArgList(rooabsarg.getVariables())
    return [argset.at(i) for i in range(argset.getSize())]


def plot_fits(pdfs, fit_results, data_obs, outfile='test.pdf'):
    """
    Plots the fitted bkg pdfs on top of the data histogram.
    """
    # First find the mT Roo variable in one of the pdfs
    for variable in get_variables(pdfs[0]):
        if variable.GetName().lower().startswith('mh'):
            break
    else:
        raise Exception('Could not determine mt variable')
    mT = variable
    mT_min = variable.getMin()
    mT_max = variable.getMax()

    # Open the frame
    xframe = mT.frame(ROOT.RooFit.Title("extended ML fit example"))
    c1 = ROOT.TCanvas()
    c1.cd()

    # Plot the data histogram
    data_obs.plotOn(xframe, ROOT.RooFit.Name("data_obs"))

    # Plot the pdfs (its parameters already at fit result)
    colors = [ROOT.kPink+6, ROOT.kBlue-4, ROOT.kRed-4, ROOT.kGreen+1]
    for pdf, color in zip(pdfs, colors):
        pdf.plotOn(
            xframe,
            ROOT.RooFit.Name(pdf.GetName()),
            ROOT.RooFit.LineColor(color),
            ROOT.RooFit.Range("Full")
            )

    # Add the fit result text labels
    for i, fit_result in enumerate(fit_results):
        n_fit_pars = len(fit_result.floatParsFinal())
        chi2 = xframe.chiSquare(pdfs[i].GetName(), "data_obs", n_fit_pars)
        txt = ROOT.TText(
            .12, 0.12+i*.05,
            "model {}, nP {}, chi2: {}".format(i, n_fit_pars, chi2)
            )
        txt.SetNDC()
        txt.SetTextSize(0.04)
        txt.SetTextColor(colors[i])
        xframe.addObject(txt) 
        txt.Draw()

    xframe.SetMinimum(0.0002)
    xframe.Draw()
    c1.SetLogy()
    c1.SaveAs(outfile)


def pdf_ploton_frame(frame, pdf, norm):
    pdf.plotOn(
        frame,
        ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),
        ROOT.RooFit.LineColor(ROOT.kBlue),
        ROOT.RooFit.FillColor(ROOT.kOrange),
        ROOT.RooFit.FillStyle(1001),
        ROOT.RooFit.DrawOption("L"),
        ROOT.RooFit.Name(pdf.GetName()),
        ROOT.RooFit.Range("Full")
        )
    pdf.paramOn(
        frame,
        ROOT.RooFit.Label(pdf.GetTitle()),
        ROOT.RooFit.Layout(0.45, 0.95, 0.94),
        ROOT.RooFit.Format("NEAU")
        )

def data_ploton_frame(frame, data, is_data=True):
    data_graph = data.plotOn(
        frame,
        ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson if is_data else ROOT.RooAbsData.SumW2),
        ROOT.RooFit.DrawOption("PE0"),
        ROOT.RooFit.Name(data.GetName())
        )
    return data_graph


def get_chi2_viaframe(mt, pdf, data, n_fit_parameters):
    """
    Get the chi2 value of the fitted pdf on data by plotting it on
    a temporary frame.

    For some reason this is much faster than the non-plotting method.
    """
    logger.debug('Using plotOn residuals')
    frame = mt.frame(ROOT.RooFit.Title(""))
    pdf_ploton_frame(frame, pdf, norm=data.sumEntries())
    data_graph = data_ploton_frame(frame, data)
    roochi2 = frame.chiSquare(pdf.GetName(), data.GetName(), n_fit_parameters)
    # Number of degrees of freedom: data will contain zeros of mt binning
    # is finer than original data binning; don't count those zeros
    dhist = frame.findObject(data.GetName(),ROOT.RooHist.Class())
    n_bins = 0
    for i in range(dhist.GetN()):
        x = ROOT.Double(0.)
        y = ROOT.Double(0.)
        dhist.GetPoint(i,x,y)
        if y!=0: n_bins += 1
    ndf = n_bins - n_fit_parameters
    chi2 = roochi2 * ndf
    roopro = ROOT.TMath.Prob(chi2, ndf)
    return roochi2, chi2, roopro, ndf

def get_chi2(mt, pdf, data, n_fit_parameters, norm=None):
    """
    Calculate chi2 of pdf on data by reading bin-by-bin the fit
    value and the data value.

    This is much slower than the `get_chi2_viaframe` function.
    Do not use this.
    """
    h_pdf = pdf.createHistogram('dummy', mt)
    h_data = data.createHistogram('dummy2', mt)
    raw_chi2 = 0.
    n_bins = 0
    for i in range(h_pdf.GetNbinsX()):
        val_pdf = h_pdf.GetBinContent(i) * (data.sumEntries() if norm is None else norm)
        val_data = h_data.GetBinContent(i)
        if val_data == 0.: continue
        raw_chi2 += (val_pdf-val_data)**2 / val_pdf
        n_bins += 1
    ndf = n_bins - n_fit_parameters
    chi2 = raw_chi2 / ndf
    from scipy import stats
    prob = stats.distributions.chi2.sf(raw_chi2, ndf) # or cdf?
    return chi2, raw_chi2, prob, ndf


def get_rss(mt, pdf, data, norm=None):
    """
    Manual RSS calculation avoiding plotting.
    For some reason this is much slower than the plotting-way.
    Do not use this.
    """
    h_pdf = pdf.createHistogram('dummy', mt)
    h_data = data.createHistogram('dummy2', mt)
    rss = 0.
    for i in range(h_pdf.GetNbinsX()):
        val_pdf = h_pdf.GetBinContent(i) * (data.sumEntries() if norm is None else norm)
        val_data = h_data.GetBinContent(i)
        if val_data == 0.: continue
        if logger.level <= logging.DEBUG:
            left_pdf = h_pdf.GetBinLowEdge(i)
            right_pdf = left_pdf + h_pdf.GetBinWidth(i)
            left_data = h_data.GetBinLowEdge(i)
            right_data = left_data + h_data.GetBinWidth(i)
            logger.debug(
                '{i}:'
                '\n  pdf  ({left_pdf:8.1f}:{right_pdf:8.1f}): {val_pdf:8.3f}'
                '\n  data ({left_data:8.1f}:{right_data:8.1f}): {val_data:8.3f}'
                .format(**locals())
                )
        rss += (val_pdf-val_data)**2
    rss = sqrt(rss)
    logger.debug('rss_manual:', rss)
    return rss


def get_rss_viaframe(mt, pdf, data, norm=None, return_n_bins=False):
    """
    Get the Residual Sum of Squares (RSS) of the fitted pdf on data by plotting it on
    a temporary frame.

    if `return_n_bins` is True, also the number of bins that were used to calculate the
    RSS.

    For some reason this is much faster than the non-plotting method.
    """
    logger.info('Calculating RSS using plotOn residuals')
    rss = 0.
    frame = mt.frame(ROOT.RooFit.Title(""))
    pdf_ploton_frame(frame, pdf, norm=(data.sumEntries() if norm is None else norm))
    data_graph = data_ploton_frame(frame, data)

    hist = data_graph.getHist()
    residuals = frame.residHist(data.GetName(), pdf.GetName(), False, True) # this is y_i - f(x_i)
    xmin, xmax = array('d', [0.]), array('d', [0.])
    data.getRange(mt, xmin, xmax)

    n_bins = 0
    for i in xrange(0, hist.GetN()):
        x, y = hist.GetX()[i], hist.GetY()[i]
        res_y = residuals.GetY()[i]
        left  = x - hist.GetErrorXlow(i)
        right = x + hist.GetErrorXhigh(i)
        if left > xmax[0] and right > xmax[0]: continue
        elif y <= 0.: continue
        if logger.level <= logging.DEBUG:
            y_pdf = y - res_y
            logger.debug(
                '{i} ({left:.2f} to {right:.2f}):'
                # '\n  pdf  : {val_pdf:8.3f}'
                '\n  data : {y:8.3f}'
                '\n  residual : {res_y:8.3f}'
                '\n  pdf : {y_pdf:8.3f}'
                .format(**locals())
                )
        rss += res_y**2
        n_bins += 1
    rss = sqrt(rss)
    logger.info('rss_viaframe: {}'.format(rss))
    return (rss, n_bins) if return_n_bins else rss


def do_fisher_test(mt, data, pdfs, n_fit_parameters_per_pdf, a_crit=.07):
    """

    """
    rsss = [ get_rss_viaframe(mt, pdf, data, return_n_bins=True) for pdf in pdfs ]
    # Compute test values of all combinations beforehand
    cl_vals = {}
    for i, j in itertools.combinations(range(len(pdfs)), 2):
        n1 = n_fit_parameters_per_pdf[i]
        n2 = n_fit_parameters_per_pdf[j]
        rss1, _      = rsss[i]
        rss2, n_bins = rsss[j]
        f = ((rss1-rss2)/(n2-n1)) / (rss2/(n_bins-n2))
        cl = 1.-ROOT.TMath.FDistI(f, n2-n1, n_bins-n2)
        cl_vals[(i,j)] = cl
    # Get the winner index
    get_winner = lambda i, j: i if cl_vals[(i,j)] > a_crit else j
    winner = get_winner(0,1)
    for i in range(2,len(pdfs)):
        winner = get_winner(winner, i)
    if logger.level <= logging.INFO:
        # Print the table
        logger.info(
            'Winner is pdf {} with {} parameters'
            .format(winner, n_fit_parameters_per_pdf[winner])
            )

        table = [[''] + range(1,len(pdfs))]
        for i in range(len(pdfs)-1):
            table.append(
                [i] + ['{:6.4f}'.format(cl_vals[(i,j)]) for j in range(i+1,len(pdfs))]
                )
        logger.info('alpha values of pdf i vs j:\n' + tabelize(table))

    return winner
    

class Datacard:
    def __init__(self):
        self.shapes = [] # Expects len-4 or len-5 iterables as elements
        self.channels = [] # Expects len-2 iterables as elements, (name, rate)
        self.rates = OrderedDict() # Expects str as key, OrderedDict as value
        self.systs = []


def parse_dc(dc):
    '''
    Very basic datacard formatter
    '''
    txt = ''
    line = '\n' + '-'*83

    txt += (
        'imax {} number of channels'
        '\njmax * number of backgrounds'
        '\nkmax * number of nuisance parameters'
        ).format(len(dc.channels))
    txt += line
    txt += '\n' + tabelize([['shapes']+s for s in dc.shapes])
    txt += line
    txt += '\n' + tabelize(transpose([('bin', 'observation')] + list(dc.channels)))
    txt += line

    # Format the bin/process table
    # Format per column, and only transpose at str format time
    proc_nr_dict = {}
    proc_nr_counter = [0]
    def proc_nr(proc):
        if not proc in proc_nr_dict:
            proc_nr_dict[proc] = proc_nr_counter[0]
            proc_nr_counter[0] += 1
        return proc_nr_dict[proc]
    table = [['bin', 'process', 'process', 'rate']]
    for bin in dc.rates:
        for proc in dc.rates[bin]:
            table.append([bin, proc, proc_nr(proc), int(dc.rates[bin][proc])])
    txt += '\n' + tabelize(transpose(table))

    txt += line
    txt += '\n' + tabelize(dc.systs)
    return txt


def transpose(l):
    '''Transposes a list of lists'''
    return map(list, zip(*l))


def tabelize(data):
    '''
    Formats a list of lists to a single string (no seps).
    Rows need not be of same length.
    '''
    # Ensure data is strings
    data = [ [ str(i) for i in row ] for row in data ]
    # Determine the row with the most columns
    n_columns = max(map(len, data))
    # Determine how wide each column should be (max)
    col_widths = [0 for i in range(n_columns)]
    for row in data:
        for i_col, item in enumerate(row):
            if len(item) > col_widths[i_col]:
                col_widths[i_col] = len(item)
    # Format
    return '\n'.join(
        ' '.join(
            format(item, str(w)) for item, w in zip(row, col_widths)
            )
        for row in data
        )

def compile_datacard_macro(bkg_pdf, data_obs, sig, outfile='dc_bsvj.txt'):
    w = ROOT.RooWorkspace("SVJ", "workspace")
    commit = getattr(w, 'import')
    commit(data_obs)
    commit(bkg_pdf, ROOT.RooFit.RenameVariable(bkg_pdf.GetName(), 'bkg'))
    commit(sig)

    wsfile = outfile.replace('.txt', '.root')
    dump_ws_to_file(wsfile, w)

    # Write the dc
    dc = Datacard()
    dc.shapes.append(['bkg', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.shapes.append([sig.GetName(), 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.shapes.append(['data_obs', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.channels.append(('bsvj', int(data_obs.sumEntries())))
    dc.rates['bsvj'] = OrderedDict()
    dc.rates['bsvj']['bkg'] = data_obs.sumEntries()
    dc.rates['bsvj'][sig.GetName()] = sig.sumEntries()
    # Freely floating bkg parameters
    for par in get_variables(bkg_pdf):
        if re.search(r'p\d_\d', par.GetName()):
            dc.systs.append([par.GetName(), 'flatParam'])
    txt = parse_dc(dc)

    logger.info('txt datacard:\n%s', txt)
    logger.info('Dumping txt to ' + outfile)
    if not osp.isdir(osp.dirname(outfile)): os.makedirs(osp.dirname(outfile))
    with open(outfile, 'w') as f:
        f.write(txt)


