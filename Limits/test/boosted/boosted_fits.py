"""
Building blocks to create the boosted SVJ analysis datacards
"""

import uuid
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
import numpy as np
import itertools, re, logging, os, os.path as osp, copy
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

def mpl_fontsizes(small=14, medium=18, large=24):
    import matplotlib.pyplot as plt
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=large)   # fontsize of the figure title


def uid():
    return str(uuid.uuid4())

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

def eval_expression(expression, pars):
    """
    Evaluates a ROOT TFormula expression in python.
    Only a limited amount of keywords are implemented (pow, log, sqrt, exp).
    """
    # Load keywords in local scope
    from numpy import log, sqrt, exp
    def pow(base, exponent):
        return base ** exponent
    # Python variables can't start with '@'; replace with some keyword
    expression = expression.replace('@', 'PARAMETER')
    # Plug parameters in local scope
    par_dict = {'PARAMETER'+str(i) : p for i, p in enumerate(pars)}
    locals().update(par_dict)
    try:
        return eval(expression)
    except NameError:
        logger.error(
            'Missing variables for expression:\n{0}\nAvailable parameters: {1}'
            .format(expression, list(par_dict.keys()))
            )
        raise

def eval_pdf_python(pdf, parameters, mt_array=None):
    if mt_array is None:
        mt = pdf.parameters[0]
        binning = mt.getBinning()
        mt_array = np.array([ binning.binCenter(i) for i in range(binning.numBins()) ])    
    parameters = list(copy.copy(parameters))
    parameters.insert(0, mt_array)
    return eval_expression(pdf.expression, parameters)

def count_parameters(expr):
    """Returns the number of parameters in an expression (i.e. highest @\d"""
    return max(map(int, re.findall(r'@(\d+)', expr))) + 1

def add_normalization(expr):
    """
    Takes an expression string, and basically adds "@NORM*(...)" around it.
    """
    return '@{0}*('.format(count_parameters(expr)) + expr + ')'

def build_rss(expr, h):
    """
    Builds a residual-sum-of-squares function between a pdf (expression)
    and a histogram.
    """
    binning, counts = th1_binning_and_values(h)
    # counts /= counts.sum() # normalize to 1
    bin_centers = [.5*(l+r) for l, r in zip(binning[:-1], binning[1:])]
    mtarray = np.array(bin_centers)
    def rss(parameters):
        # Insert mT array as first parameter
        parameters = list(copy.copy(parameters))
        parameters.insert(0, mtarray)
        y_pdf = eval_expression(expr, parameters)
        # Normalize pdf to counts too before evaluating, so as the compare only shape
        y_pdf = (y_pdf/y_pdf.sum()) * counts.sum()
        return np.sqrt(np.sum((counts-y_pdf)**2))
    return rss

def build_chi2(expr, h):
    """
    Builds a chi2 function between a pdf (expression) and a histogram.
    """
    binning, counts = th1_binning_and_values(h)
    # Use the bin centers as the mT array
    mt_array = np.array(.5*(binning[:-1]+binning[1:]))
    def chi2(parameters):
        # Construct the parameters of the expression:
        # [ @0 (mt), @1, ... @N (pdf parameters) ]
        parameters = list(copy.copy(parameters))
        parameters.insert(0, mt_array)
        y_pdf = eval_expression(expr, parameters)
        # Normalize pdf to counts too before evaluating, so as the compare only shape
        y_pdf = (y_pdf/y_pdf.sum()) * counts.sum()
        return np.sum((counts-y_pdf)**2 / y_pdf)
    return chi2


FIT_CACHE_FILE = 'fit_cache.pickle'


def _fit_hash(expression, histogram, init_vals=None, **minimize_kwargs):
    import hashlib
    m = hashlib.sha256()
    def add_floats_to_hash(floats):
        for number in floats:
            s = '{:.5f}'.format(number)
            m.update(s)
    m.update(expression)
    binning, values = th1_binning_and_values(histogram)
    add_floats_to_hash(binning)
    add_floats_to_hash(values)
    if init_vals is not None: add_floats_to_hash(init_vals)
    if 'tol' in minimize_kwargs: m.update('{:.3f}'.format(minimize_kwargs['tol']))
    if 'method' in minimize_kwargs: m.update(minimize_kwargs['method'])
    return m.hexdigest()


def _fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=None, hash=None, **minimize_kwargs):
    """
    The actual entry point to the scipy fit
    """
    n_fit_pars = count_parameters(expression) - 1 # -1 because par 0 is mT
    logger.info('Fitting {0} with {1} parameters'.format(expression, n_fit_pars))
    from scipy.optimize import minimize
    chi2 = build_chi2(expression, histogram)
    if init_vals is None: init_vals = np.ones(n_fit_pars)
    res = minimize(chi2, init_vals, **minimize_kwargs)
    # Save some extra information in the result
    res.x_init = np.array(init_vals)
    res.expression = expression
    res.hash = hash
    # Set approximate uncertainties; see https://stackoverflow.com/a/53489234
    # Assume ftol ~ function value
    # try:
    #     res.dx = np.sqrt(res.fun * np.diagonal(res.hess_inv))
    # except:
    #     logger.error('Failed to set uncertainties; using found function values as proxies')
    #     res.dx = res.x.copy()
    return res

def _read_fit_cache():
    import pickle
    if osp.isfile(FIT_CACHE_FILE):
        with open(FIT_CACHE_FILE, 'rb') as f:
            return pickle.load(f)
    else:
        return {}

def _write_fit_cache(cache_dict):
    import pickle
    with open(FIT_CACHE_FILE, 'wb') as f:
        pickle.dump(cache_dict, f)

def _get_from_fit_cache(fit_hash):
    return _read_fit_cache().get(fit_hash, None)

def _add_one_to_fit_cache(key, result):
    d = _read_fit_cache()
    d[key] = result
    _write_fit_cache(d)

def brute_force_init_vals(npars, values):
    import itertools
    return np.array(list(itertools.product(*[values for i in range(npars)])))


def fit_expr_to_histogram_robust(expression, histogram):
    """
    Heuristic around `fit_pdf_expression_to_histogram_python`.
    First attempts a single fit, and only goes for the bruteforce if the single fit
    did not converge properly.
    """
    fit_hash = _fit_hash(expression, histogram)
    res = _get_from_fit_cache(fit_hash)
    if res is not None:
        logger.warning('Returning cached robust fit')
        return res
    res = fit_pdf_expression_to_histogram_python(
        expression, histogram,
        tol=1e-3, method='BFGS'
        )
    # Refit with output from first fit
    res = fit_pdf_expression_to_histogram_python(
        expression, histogram,
        init_vals=res.x,
        tol=1e-6, method='Nelder-Mead'
        )
    if not res.success:
        logger.info('Fit did not converge with single try; brute forcing it')
        results = []
        for method in ['BFGS', 'Nelder-Mead']:
            results = fit_pdf_expression_to_histogram_python(
                expression, histogram,
                init_vals=brute_force_init_vals(count_parameters(expression)-1, [-1., 1.]),
                tol=1e-3, method=method
                )
        results = [ r for r in results if not(np.isnan(r.fun) or np.isposinf(r.fun) or np.isneginf(r.fun)) ]
        if len(results) == 0: raise Exception('Not a single fit of the brute force converged!')
        i_min = np.argmin([r.fun for r in results])        
        res = results[i_min]
    logger.info('Final scipy fit result:\n%s', res)
    _add_one_to_fit_cache(fit_hash, res)
    return res


def fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=None, cache=True, cache_dict=None, **minimize_kwargs):
    """
    Fits a a background pdf expression to a TH1 histogram with scipy.optimize.minimize.
    Assumes @0 in the expression is mT, and the histogram is binned in mT.

    If `cache` is True, it tries to save the fit result to a file. If init_vals is
    a 2-dimensional array, the fit is repeated for each initial value.

    If `cache_dict` is given, no read/write to the cache file is performed.
    """
    _write_cache = False
    if cache and cache_dict is None:
        # If cache is enabled, and no cache dict was specified, treat this as one fit result
        # and do the read/write IO of the cache file
        cache_dict = _read_fit_cache()
        _write_cache = True
    
    if init_vals is not None:
        init_vals = np.array(init_vals)
        if len(init_vals.shape) > 1:
            logger.info('Will run fit for %s different initial values', init_vals.shape[0])
            results = [
                fit_pdf_expression_to_histogram_python(
                    expression, histogram, init_vals=x_init,
                    cache=cache, cache_dict=cache_dict, **minimize_kwargs
                    ) for x_init in init_vals
                ]
            _write_fit_cache(cache_dict)
            return results

    fit_hash = _fit_hash(expression, histogram, init_vals, **minimize_kwargs)

    if cache and fit_hash in cache_dict:
        # Nothing new to save, so return immediately
        logger.warning('Returning cached fit')
        return cache_dict[fit_hash]
    else:
        res = _fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=init_vals, hash=fit_hash, **minimize_kwargs)
        cache_dict[fit_hash] = res
        if _write_cache: _write_fit_cache(cache_dict)
        return res


def get_mt(mt_min=160., mt_max=500., n_bins=43, name='mHbsvj'):
    """
    Sensible defaults for the mt axis
    """
    mt = ROOT.RooRealVar(name, 'm_{T}', mt_min, mt_max, 'GeV')
    mt.setBins(n_bins)
    # Manually add the boundaries to it as python attributes for easy access
    mt.mt_min = mt_min
    mt.mt_max = mt_max
    return mt

def get_mt_from_th1(histogram, name=None):
    """
    Returns mT from the x axis of a TH1 histogram.
    Min and max are simply the left/right boundary of the first/last bin,
    and bin width is copied.
    """
    mt = get_mt(
        histogram.GetBinLowEdge(1),
        histogram.GetBinLowEdge(histogram.GetNbinsX()+1),
        histogram.GetNbinsX(),
        name = uid() if name is None else name
        )
    object_keeper.add(mt)
    return mt


class ROOTObjectKeeper:
    """
    Keeps ROOT objects in it so ROOT doesn't garbage clean them.
    """
    def __init__(self):
        self.objects = {}
    
    def add(self, thing):
        try:
            key = thing.GetName()
        except AttributeError:
            key = str(uuid.uuid4())
        if key in self.objects:
            logger.warning('Variable %s (%s) already in object keeper', thing.GetName(), thing.GetTitle())
        self.objects[key] = thing

    def add_multiple(self, things):
        for t in things: self.add(t)


object_keeper = ROOTObjectKeeper()


def build_rpsbp(name, expression, mt, parameters, histogram):
    '''Builds a RooParametricShapeBinPdf'''
    logger.info(
        'Building pdf {}, expression "{}", mt.GetName()={}, parameter_names={}'
        .format(name, expression, mt.GetName(), ', '.join([p.GetName() for p in parameters]))
        )
    generic_pdf = ROOT.RooGenericPdf(
        name+'_rgp', name+'_rgp',
        expression, ROOT.RooArgList(mt, *parameters)
        )
    object_keeper.add(generic_pdf)
    parametric_shape_bin_pdf = ROOT.RooParametricShapeBinPdf(
        name+'_rpsbp', name+'_rpsbp',
        generic_pdf, mt, ROOT.RooArgList(*parameters), histogram
        )
    parametric_shape_bin_pdf.expression = expression # Tag it onto the instance
    parametric_shape_bin_pdf.parameters = [mt] + parameters
    parametric_shape_bin_pdf.npars = len(parameters)
    parametric_shape_bin_pdf.histogram = histogram
    parametric_shape_bin_pdf.pdftype = 'main' if 'main' in name else 'alt'
    object_keeper.add(parametric_shape_bin_pdf)
    logger.info(
        'Created RooParametricShapeBinPdf {} with {} parameters on histogram {}'
        .format(
            parametric_shape_bin_pdf.GetName(),
            len(parameters), histogram.GetName()
            )
        )
    return parametric_shape_bin_pdf

def rebuild_rpsbp(pdf):
    name = uid()
    def remake_parameter(parameter):
        variable = ROOT.RooRealVar(
            name + '_' + parameter.GetName(), parameter.GetTitle(),
            1., parameter.getMin(), parameter.getMax()
            )
        object_keeper.add(variable)
        return variable
    return build_rpsbp(
        name, pdf.expression, pdf.parameters[0],
        [remake_parameter(p) for p in pdf.parameters[1:]], pdf.histogram
        )

def pdf_expression(pdf_type, npars, mt_scale='1000'):
    if pdf_type == 'main':
        if npars == 2:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2))'
        elif npars == 3:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})))'
        elif npars == 4:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})+@4*pow(log(@0/{0}),2)))'
            # Alternatives to 22:
            # 13: pow(1 - @0/{0}, @1+@2*log(@0/{0})) * pow(@0/{0}, -(@3+@4*log(@0/{0})))
        elif npars == 5:
            expression = 'pow(1 - @0/{0}, @1+@2*log(@0/{0})+@3*pow(log(@0/{0}),2)) * pow(@0/{0}, -(@4+@5*log(@0/{0})))'
            # Alternatives to 32:
            # 14: pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})+@4*pow(log(@0/{0}),2)+@5*pow(log(@0/{0}),3)))
            # 41: pow(1 - @0/{0}, @1+@2*log(@0/{0})+@3*pow(log(@0/{0}),2)+@4*pow(log(@0/{0}),3)) * pow(@0/{0}, -@5)
        else:
            raise Exception('Unavailable npars for main: {0}'.format(npars))
    elif pdf_type == 'alt':
        if npars == 1:
            expression = 'exp(@1*(@0/{0}))'
        elif npars == 2:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2)'
        elif npars == 3:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2*(1+@3*log(@0/{0})))'
        elif npars == 4:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2*(1+@3*log(@0/{0})*(1+@4*log(@0/{0}))))'
        else:
            raise Exception('Unavailable npars for alt: {0}'.format(npars))
    else:
        raise Exception('Unknown pdf type {0}'.format(pdf_type))
    return expression.format(mt_scale)


def make_pdf(pdf_type, npars, bkg_hist, mt=None, name=None, mt_scale='1000', **kwargs):
    if mt is None: mt = get_mt_from_th1(bkg_hist)
    if pdf_type == 'alt':
        return make_alt_pdf(npars, bkg_hist, mt, name=name, mt_scale=mt_scale, **kwargs)
    elif pdf_type == 'main':
        return make_main_pdf(npars, bkg_hist, mt, name=name, mt_scale=mt_scale, **kwargs)
    else:
        raise Exception('Unknown pdf type {0}'.format(pdf_type))

def make_main_pdf(npars, bkg_hist, mt, name=None, mt_scale='1000'):
    # Function from Theorists, combo testing, sequence E, 1, 11, 12, 22
    # model NM has N params on 1-x and M params on x. exponents are (p_i + p_{i+1} * log(x))
    pdf_name = 'bsvj_bkgfitmain_npars'+str(npars) if name is None else name
    expression = pdf_expression('main', npars, mt_scale)
    if npars == 2:
        parameters = [
            ROOT.RooRealVar(pdf_name + "_p1", "p1", 1., -45., 45.),
            ROOT.RooRealVar(pdf_name + "_p2", "p2", 1., -10., 10.)
            ]
    elif npars == 3:
        parameters = [
            ROOT.RooRealVar(pdf_name + "_p1", "p1", 1., -45., 45.),
            ROOT.RooRealVar(pdf_name + "_p2", "p2", 1., -10., 10.),
            ROOT.RooRealVar(pdf_name + "_p3", "p3", 1., -15, 15),
            ]
    elif npars == 4:
        parameters = [
            ROOT.RooRealVar(pdf_name + "_p1", "p1", 1., -95., 95.),
            ROOT.RooRealVar(pdf_name + "_p2", "p2", 1., -25., 20.),
            ROOT.RooRealVar(pdf_name + "_p3", "p3", 1., -2., 2.),
            ROOT.RooRealVar(pdf_name + "_p4", "p4", 1., -2., 2.),
            ]
    elif npars == 5:
        parameters = [
            ROOT.RooRealVar(pdf_name + "_p1", "p1", 1., -15., 15.),
            ROOT.RooRealVar(pdf_name + "_p2", "p2", 1., -95., 95.),
            ROOT.RooRealVar(pdf_name + "_p3", "p3", 1., -25., 25.),
            ROOT.RooRealVar(pdf_name + "_p4", "p4", 1., -5., 5.),
            ROOT.RooRealVar(pdf_name + "_p5", "p5", 1., -1.5, 1.5),
            ]
    object_keeper.add_multiple(parameters)
    return build_rpsbp(pdf_name, expression.format(mt_scale), mt, parameters, bkg_hist)

def make_alt_pdf(npars, bkg_hist, mt, name=None, mt_scale='1000', par_lo=-50., par_up=50.):
    pdf_name = 'bsvj_bkgfitalt_npars'+str(npars) if name is None else name
    parameters = [
        ROOT.RooRealVar(pdf_name + '_p{0}'.format(i+1), '', 1., par_lo, par_up) \
        for i in range(npars)
        ]
    expression = pdf_expression('alt', npars, mt_scale)
    object_keeper.add_multiple(parameters)
    return build_rpsbp(pdf_name, expression.format(mt_scale), mt, parameters, bkg_hist)


def get_main_pdfs(bkg_hist, mt, mt_scale='1000'):
    return [make_main_pdf(npars, bkg_hist, mt, mt_scale=mt_scale) for npars in range(2,6)]

def get_alt_pdfs(bkg_hist, mt, mt_scale='1000'):
    return [make_alt_pdf(npars, bkg_hist, mt, mt_scale=mt_scale) for npars in range(1,5)]

def get_pdfs(pdftype, *args, **kwargs):
    if pdftype == 'main':
        return get_main_pdfs(*args, **kwargs)
    elif pdftype == 'alt':
        return get_alt_pdfs(*args, **kwargs)
    else:
        raise Exception('No such pdftype: {}'.format(pdftype))


def fit_pdf_to_datahist(pdf, data_hist, init_vals=None, init_ranges=None):
    logger.info('Fitting pdf {0} to data_hist {1} with RooFit'.format(pdf, data_hist))

    if init_vals is not None:
        if len(init_vals) != len(pdf.parameters)-1:
            raise Exception('Expected {} values; got {}'.format(len(pdf.parameters)-1, len(init_vals)))
        for par, value in zip(pdf.parameters[1:], init_vals):
            old_range = max(abs(par.getMin()), abs(par.getMax()))
            if abs(value) > old_range:
                new_range = 1.1*abs(value)
                logger.info(
                    'Increasing range for {} ({}) from ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})'
                    .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax(), -new_range, new_range)
                    )
                par.setRange(-new_range, new_range)
            elif abs(value) / old_range < .1:
                new_range = 1.1*abs(value)
                logger.info(
                    'Decreasing range for {} ({}) from ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})'
                    .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax(), -new_range, new_range)
                    )
                par.setRange(-new_range, new_range)
            # par.setRange(value-0.01*abs(value), value+0.01*abs(value))
            par.setVal(value)
            logger.info(
                'Setting {0} ({1}) value to {2}, range is {3} to {4}'
                .format(par.GetName(), par.GetTitle(), value, par.getMin(), par.getMax())
                )
    if init_ranges is not None:
        if len(init_ranges) != len(pdf.parameters)-1:
            raise Exception('Expected {} values; got {}'.format(len(pdf.parameters)-1, len(init_ranges)))
        for par, (left, right) in zip(pdf.parameters[1:], init_ranges):
            par.setRange(left, right)
            logger.info(
                'Setting {0} ({1}) range to {2} to {3}'
                .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax())
                )

    print '><' * 40
    print pdf.parameters[0].GetName(), data_hist.mt.GetName()

    try:
        res = pdf.fitTo(
            data_hist,
            ROOT.RooFit.Extended(False),
            ROOT.RooFit.Save(1),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Minimizer("Minuit2"),
            ROOT.RooFit.PrintLevel(2 if logger.level <= logging.DEBUG else -1),
            ROOT.RooFit.Range('Full'),
            ROOT.RooFit.PrintEvalErrors(-1)
            )
    except:
        logger.error('Problem fitting pdf {}'.format(pdf.GetName()))
        raise
    return res


def to_list(rooarglist):
    return [rooarglist.at(i) for i in range(rooarglist.getSize())]

def get_variables(rooabsarg):
    """
    Returns a list of all variables a RooAbsArg depends on
    """
    argset = ROOT.RooArgList(rooabsarg.getVariables())
    return [argset.at(i) for i in range(argset.getSize())]


def set_pdf_to_fitresult(pdf, res):
    """
    Sets the parameters of a pdf to the fit result. 
    """
    def set_par(par, value):
        par.setRange(value-10., value+10.)
        par.setVal(value)
    import scipy
    if isinstance(res, ROOT.RooFitResult):
        vals = []
        for p_fit, p_pdf in zip(to_list(res.floatParsFinal()), pdf.parameters[1:]):
            set_par(p_pdf, p_fit.getVal())
            vals.append(p_fit.getVal())
        return vals
    elif isinstance(res, scipy.optimize.optimize.OptimizeResult):
        for val, p_pdf in zip(res.x, pdf.parameters[1:]):
            set_par(p_pdf, val)
        return res.x


def plot_pdf_for_various_fitresults(pdf, fit_results, data_obs, outfile='test.pdf', labels=None, title=''):
    """
    Plots the fitted bkg pdfs on top of the data histogram.
    """
    # First find the mT Roo variable in one of the pdfs
    mt = pdf.parameters[0]
    mt_min = mt.getMin()
    mt_max = mt.getMax()

    # Open the frame
    xframe = mt.frame(ROOT.RooFit.Title(title))
    c1 = ROOT.TCanvas(str(uuid.uuid4()), '', 1000, 800)
    c1.cd()

    # Plot the data histogram
    data_obs.plotOn(xframe, ROOT.RooFit.Name("data_obs"))
    norm = data_obs.sumEntries()

    # Plot the pdfs (its parameters already at fit result)
    colors = [ROOT.kPink+6, ROOT.kBlue-4, ROOT.kRed-4, ROOT.kGreen+1]

    py_chi2 = build_chi2(pdf.expression, data_obs.createHistogram(mt.GetName()))

    base_pdf = pdf
    for i, res in enumerate(fit_results):
        pdf = rebuild_rpsbp(base_pdf)

        print '>'*80
        print mt.GetName(), pdf.parameters[0].GetName(), data_obs.mt.GetName()

        vals = set_pdf_to_fitresult(pdf, res)
        logger.info(
            'i=%s; Manual chi2=%.5f, chi2_via_frame=%.5f',
            i, py_chi2(vals), get_chi2_viaframe(mt, pdf, data_obs, len(vals))[1]
            )
        pdf.plotOn(
            xframe,
            ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),
            ROOT.RooFit.LineColor(colors[i]),
            # ROOT.RooFit.FillColor(ROOT.kOrange),
            ROOT.RooFit.FillStyle(1001),
            ROOT.RooFit.DrawOption("L"),
            ROOT.RooFit.Name(pdf.GetName()),
            ROOT.RooFit.Range("Full")
            )
        chi2 = xframe.chiSquare(pdf.GetName(), "data_obs", len(pdf.parameters)-1)
        par_value_str = ', '.join(['p{}={:.3f}'.format(iv, v) for iv, v in enumerate(vals)])
        label = labels[i] if labels else 'fit'+str(i)
        txt = ROOT.TText(
            .13, 0.13+i*.045,
            "{}, chi2={:.4f}, {}".format(label, chi2, par_value_str)
            )
        txt.SetNDC()
        txt.SetTextSize(0.03)
        txt.SetTextColor(colors[i])
        xframe.addObject(txt) 
        txt.Draw()

    xframe.SetMinimum(0.002)
    xframe.Draw()
    c1.SetLogy()
    c1.SaveAs(outfile)
    if outfile.endswith('.pdf'): c1.SaveAs(outfile.replace('.pdf', '.png'))
    del xframe, c1


def plot_fits(pdfs, fit_results, data_obs, outfile='test.pdf'):
    """
    Plots the fitted bkg pdfs on top of the data histogram.
    """
    # First find the mT Roo variable in one of the pdfs
    mT = pdfs[0].parameters[0]
    mT_min = mT.getMin()
    mT_max = mT.getMax()

    # Open the frame
    xframe = mT.frame(ROOT.RooFit.Title("extended ML fit example"))
    c1 = ROOT.TCanvas()
    c1.cd()

    # Plot the data histogram
    data_obs.plotOn(xframe, ROOT.RooFit.Name("data_obs"))

    # Set to fitresult
    for pdf, res in zip(pdfs, fit_results): set_pdf_to_fitresult(pdf, res)

    # Plot the pdfs
    colors = [ROOT.kPink+6, ROOT.kBlue-4, ROOT.kRed-4, ROOT.kGreen+1]
    colors.extend(colors)
    colors.extend(colors)
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

        par_values = [ 'p{}={:.3f}'.format(i, v.getVal()) for i, v in enumerate(pdfs[i].parameters[1:])]
        par_value_str = ', '.join(par_values)
        
        txt = ROOT.TText(
            .12, 0.12+i*.05,
            "model {}, nP {}, chi2: {:.4f}, {}".format(i, n_fit_pars, chi2, par_value_str)
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
    c1.SaveAs(outfile.replace('.pdf', '.png'))
    del xframe, c1


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


def do_fisher_test(mt, data, pdfs, a_crit=.07):
    """

    """
    rsss = [ get_rss_viaframe(mt, pdf, data, return_n_bins=True) for pdf in pdfs ]
    # Compute test values of all combinations beforehand
    cl_vals = {}
    for i, j in itertools.combinations(range(len(pdfs)), 2):
        n1 = pdfs[i].npars
        n2 = pdfs[j].npars
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
            .format(winner, pdfs[winner].npars)
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

    # Careful with workspace path: Should be relative path to DC
    shapes = copy.copy(dc.shapes)
    for i in range(len(shapes)):
        shapes[i][2] = osp.basename(shapes[i][2])

    txt += '\n' + tabelize([['shapes']+s for s in shapes])
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
    txt += '\n'
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

def make_multipdf(pdfs, name='roomultipdf'):
    # ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    cat = ROOT.RooCategory('pdf_index', "Index of Pdf which is active")
    pdf_arglist = ROOT.RooArgList()
    for pdf in pdfs:
        pdf_arglist.add(pdf)
    multipdf = ROOT.RooMultiPdf(name, "All Pdfs", cat, pdf_arglist)
    multipdf.cat = cat
    multipdf.pdfs = pdfs
    norm = ROOT.RooRealVar(name+'_norm', "Number of background events", 1.0, 0., 1.e6)
    object_keeper.add(multipdf)
    object_keeper.add(norm)
    return multipdf, norm

def compile_datacard_macro(bkg_pdf, data_obs, sig, outfile='dc_bsvj.txt', systs=None):
    do_syst = systs is not None
    w = ROOT.RooWorkspace("SVJ", "workspace")

    def commit(thing, *args, **kwargs):
        print(thing)
        name = thing.GetName() if hasattr(thing, 'GetName') else '?'
        print('Importing {} ({})'.format(name, thing))
        logger.info('Importing {} ({})'.format(name, thing))
        getattr(w, 'import')(thing, *args, **kwargs)

    # Bkg pdf: May be multiple
    is_multipdf = hasattr(bkg_pdf, '__len__')
    if is_multipdf:
        # Multipdf
        multipdf, norm = make_multipdf(bkg_pdf)
        for pdf in multipdf.pdfs:
            commit(pdf)
        commit(multipdf.cat)
        commit(norm)
        commit(multipdf)
    else:
        commit(bkg_pdf, ROOT.RooFit.RenameVariable(bkg_pdf.GetName(), 'bkg'))

    commit(data_obs)
    commit(sig, ROOT.RooFit.RenameVariable(sig.GetName(), 'sig'))

    wsfile = outfile.replace('.txt', '.root')
    dump_ws_to_file(wsfile, w)

    # Write the dc
    dc = Datacard()
    dc.shapes.append(['roomultipdf' if is_multipdf else 'bkg', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.shapes.append(['sig', 'bsvj', wsfile, 'SVJ:$PROCESS'] + (['SVJ:$PROCESS_$SYSTEMATIC'] if do_syst else []))
    dc.shapes.append(['data_obs', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.channels.append(('bsvj', int(data_obs.sumEntries())))
    dc.rates['bsvj'] = OrderedDict()
    dc.rates['bsvj']['sig'] = sig.sumEntries()
    dc.rates['bsvj']['roomultipdf' if is_multipdf else 'bkg'] = data_obs.sumEntries()
    # Freely floating bkg parameters
    def systs_for_pdf(pdf):
        for par in pdf.parameters[1:]:
            dc.systs.append([par.GetName(), 'flatParam'])
    [systs_for_pdf(p) for p in multipdf.pdfs] if is_multipdf else systs_for_pdf(bkg_pdf)
    # Rest of the systematics
    if is_multipdf: dc.systs.append([multipdf.cat.GetName(), 'discrete'])
    if do_syst: dc.systs.extend(systs)
    txt = parse_dc(dc)

    logger.info('txt datacard:\n%s', txt)
    logger.info('Dumping txt to ' + outfile)
    if not osp.isdir(osp.dirname(outfile)): os.makedirs(osp.dirname(outfile))
    with open(outfile, 'w') as f:
        f.write(txt)


def th1_binning_and_values(h):
    """
    Returns the binning and values of the histogram.
    Does not include the overflows.
    """
    n_bins = h.GetNbinsX()
    # GetBinLowEdge of the right overflow bin is the high edge of the actual last bin
    binning = np.array([h.GetBinLowEdge(i) for i in range(1,n_bins+2)])
    values = np.array([h.GetBinContent(i) for i in range(1,n_bins+1)])
    return binning, values

def th1_to_datahist(histogram, mt=None):
    if mt is None: mt = get_mt_from_th1(histogram)
    datahist = ROOT.RooDataHist(uid(), '', ROOT.RooArgList(mt), histogram, 1.)
    datahist.mt = mt
    return datahist
