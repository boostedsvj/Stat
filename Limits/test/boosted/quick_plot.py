"""
Some scripts to quickly plot basic outcomes from combine scans
"""
from __future__ import print_function
import ROOT
from time import strftime
import argparse

# Add the directory of this file to the path so the boosted tools can be imported
import sys, os, os.path as osp, pprint, re, traceback
sys.path.append(osp.dirname(osp.abspath(__file__)))
import boosted_fits as bsvj
logger = bsvj.setup_logger('quickplot')

import numpy as np
import matplotlib.pyplot as plt

def set_mpl_fontsize(small=16, medium=22, large=26):
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=medium)    # legend fontsize
    plt.rc('figure', titlesize=large)  # fontsize of the figure title
set_mpl_fontsize()

_scripts = {}
def is_script(fn):
    _scripts[fn.__name__] = fn
    return fn

def cmd_exists(executable):
    """
    Checks if a command can be found on the system path.
    Not a very smart implementation but does the job usually.
    See https://stackoverflow.com/a/28909933/9209944 .
    """
    return any(os.access(os.path.join(path, executable), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

def name_from_combine_rootfile(rootfile):
    return osp.basename(rootfile).split('.',1)[0].replace('higgsCombine','')

def quickplot_parser(outfile='test.png'):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile', type=str, default=outfile)
    parser.add_argument('-b', '--batch', action='store_true')
    return parser

def namespace_to_attrdict(args):
    return bsvj.AttrDict(**vars(args))


def extract_scan(rootfile, correct_minimum=False):
    with bsvj.open_root(rootfile) as tf:
        limit = tf.Get('limit')
        logger.info(limit)
        mus = []
        deltanlls = []
        for _ in limit:
            mus.append(limit.r)
            deltanlls.append(limit.deltaNLL)

    mus = np.array(mus)
    deltanlls = np.array(deltanlls)
    order = np.argsort(mus)
    mus = mus[order]
    deltanlls = deltanlls[order]

    if correct_minimum:
        minimum = np.min(deltanlls)
        logger.warning('Shifting curve by {0:.4f}'.format(minimum))
        deltanlls = deltanlls - minimum

    return mus, deltanlls


@is_script
def muscan(args):
    parser = quickplot_parser('muscan.png')
    parser.add_argument('rootfiles', type=str, nargs='+')
    parser.add_argument('--xmin', type=float)
    parser.add_argument('--xmax', type=float)
    parser.add_argument('--ymin', type=float)
    parser.add_argument('--ymax', type=float)
    parser.add_argument('--correctminimum', action='store_true')
    args = parser.parse_args(args)

    fig = plt.figure(figsize=(8,8))
    ax = fig.gca()

    min_mu = 1e6
    max_mu = -1e6
    for rootfile in args.rootfiles:
        mus, deltanlls = extract_scan(rootfile, args.correctminimum)
        min_mu = min(min_mu, np.min(mus))
        max_mu = max(max_mu, np.max(mus))
        ax.plot(mus, deltanlls, '-o', label=name_from_combine_rootfile(rootfile))

    ax.plot([min_mu, max_mu], [.0, .0], color='lightgray')
    ax.plot([min_mu, max_mu], [.5, .5], label='$1\sigma$')
    ax.plot([min_mu, max_mu], [1., 1.], label='$2\sigma$')

    ax.set_xlabel('$\mu$')
    ax.set_ylabel('$\Delta NLL$')
    if args.xmax: ax.set_xlim(right=args.xmax)
    if args.xmin: ax.set_xlim(left=args.xmin)
    if args.ymax: ax.set_ylim(top=args.ymax)
    if args.ymin: ax.set_ylim(bottom=args.ymin)
    ax.legend()

    plt.savefig(args.outfile, bbox_inches='tight')
    if not(args.batch) and cmd_exists('imgcat'): os.system('imgcat ' + args.outfile)



@is_script
def mtdist(args):
    from scipy.interpolate import make_interp_spline
    if not isinstance(args, bsvj.AttrDict):
        parser = quickplot_parser('mtdist.png')
        parser.add_argument('rootfile', type=str)
        args = parser.parse_args(args)

    with bsvj.open_root(args.rootfile) as f:
        ws = bsvj.get_ws(f)
    
    ws.loadSnapshot('MultiDimFit')

    mu = ws.var('r').getVal()

    mt = ws.var('mt')
    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[1:]+mt_binning[:-1])
    mt_bin_widths = mt_binning[1:] - mt_binning[:-1]

    data = ws.data('data_obs')
    _, y_data = bsvj.roodataset_values(data)
    errs_data = np.sqrt(y_data)
    # errs_data = y_data

    bkg = ws.pdf('shapeBkg_roomultipdf_bsvj')
    y_bkg = bsvj.pdf_values(bkg, mt_bin_centers)
    bkg_norm = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()

    sig = ws.embeddedData('shapeSig_sig_bsvj')
    _, y_sig = bsvj.roodataset_values(sig)


    bsvj.logger.info(
        'Found {} entries in data, {} entries in signal (should match with datacard!)'
        .format(y_data.sum(), y_sig.sum())
        )

    fig = plt.figure(figsize=(10,10))
    ax = fig.gca()

    ax.plot([], [], ' ', label=name_from_combine_rootfile(args.rootfile))

    ax.errorbar(
        mt_bin_centers, y_data,
        xerr=.5*mt_bin_widths, yerr=errs_data,
        fmt='o', c='black', label='Data'
        )

    mt_fine = np.linspace(mt_binning[0], mt_binning[-1], 100) # For fine plotting
    spl = make_interp_spline(mt_bin_centers, y_bkg, k=3)  # type: BSpline
    y_bkg_fine = spl(mt_fine)
    ax.plot(mt_fine, y_bkg_fine*bkg_norm, label='Bkg (fit)', c='b')

    y_sb = bkg_norm*y_bkg + mu*y_sig
    ax.step(
        mt_binning[:-1], y_sb, where='post', c='r',
        label=r'B+$\mu_{{fit}}$S ($\mu_{{fit}}$={0:.1f})'.format(mu)
        )

    ax.step(mt_binning[:-1], y_sig, where='post', label=r'S ($\mu$=1)', c='g')

    ax.legend()
    ax.set_ylabel('$N_{events}$')
    ax.set_xlabel(r'$m_{T}$ (GeV)')
    ax.set_yscale('log')

    plt.savefig(args.outfile, bbox_inches='tight')
    if not(args.batch) and cmd_exists('imgcat'): os.system('imgcat ' + args.outfile)



def get_cls(obs_arrays, asimov_arrays):
    from scipy.stats import norm
    quantiles = np.array([0.025, 0.16, 0.50, 0.84, 0.975])

    def get_mu_dnll_quantiles(dct):
        mu = dct['r']
        dnll = dct['deltaNLL']
        dnll -= np.min(dnll) # Correct for bad first fit
        quantile = dct['quantileExpected']
        # Take out the bestfit
        is_bestfit = quantile==-1.
        i_bestfit = is_bestfit.argmax()
        mu_best = mu[i_bestfit]
        dnll_best = dnll[i_bestfit]
        mu = mu[~is_bestfit]
        dnll = dnll[~is_bestfit]
        quantile = quantile[~is_bestfit]
        # Sort by ascending mu
        order = np.argsort(mu)
        mu = mu[order]
        dnll = dnll[order]
        quantile = quantile[order]
        return bsvj.AttrDict(mu_best=mu_best, dnll_best=dnll_best, mu=mu, dnll=dnll, quantile=quantile, n=mu.shape[0])

    def align_mu_values(obs, asimov):
        """
        Make sure all arrays in obs and asimov concern the same mu values
        """
        i_obs_in_asimov = np.isin(obs.mu, asimov.mu)
        i_asimov_in_obs = np.isin(asimov.mu, obs.mu)
        for key in ['mu', 'dnll', 'quantile']:
            obs[key] = obs[key][i_obs_in_asimov]
            asimov[key] = asimov[key][i_asimov_in_obs]
        obs.n = obs.mu.shape[0]
        asimov.n = asimov.mu.shape[0]

    obs = get_mu_dnll_quantiles(obs_arrays)
    asimov = get_mu_dnll_quantiles(asimov_arrays)
    align_mu_values(obs, asimov)
    np.testing.assert_array_equal(obs.mu, asimov.mu)

    # I do not understand why not simply q_obs = 2.*dnll_obs?
    # Or is this just to offset a potentially bad best fit?
    # If so, why not just shift the whole dnll array so its minimum is at 0...
    q_obs = []
    for i, mu in enumerate(obs.mu):
        if mu < obs.mu_best:
            dnll_obs_min = np.min(obs.dnll[:i+1])  # Why?
            dnll_obs_constrained = obs.dnll[i] - dnll_obs_min
        else:
            dnll_obs_constrained = obs.dnll[i]
        q_obs.append(2.*max(dnll_obs_constrained, 0.))
    q_obs = np.array(q_obs)
    assert q_obs.shape == (obs.n,)

    q_A = 2. * asimov.dnll
    q_A[q_A < 0.] = 0.  # Set negative values to 0

    # Also this formula I don't fully understand
    s_exp = { q : (1.-norm.cdf(np.sqrt(q_A) - norm.ppf(q))) / q for q in quantiles}

    assert np.all(  ((q_obs >= 0.) & (q_obs <= q_A)) | (q_obs > q_A)  )

    # This is just black magic
    sb = np.where(
        q_obs <= q_A,
        1. - norm.cdf( np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs+q_A) , np.sqrt(q_obs)) )
        )
    b = np.where(
        q_obs <= q_A,
        norm.cdf( np.sqrt(q_A)-np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs-q_A) , np.sqrt(q_obs)) )
        )
    s = sb / b
    return bsvj.AttrDict(s=s, b=b, sb=sb, q_obs=q_obs, q_A=q_A, obs=obs, asimov=asimov, s_exp=s_exp)


def interpolate_95cl_limit(cls):
    mu = cls.obs.mu
    def interpolate(cl):
        # print('Interpolating')
        select = ((cl < .15) & (cl > .01) & (mu>0))
        # print('  {} values left'.format(select.sum()))
        order = np.argsort(cl[select])
        # print('  {:14s}  {:14s}'.format('cl', 'mu'))
        # for c, m in zip(cl[select][order], mu[select][order]):
        #     print('  {:+14.7f}  {:+14.7f}'.format(c, m))
        res = np.interp(.05, cl[select][order], mu[select][order])
        # print('Interpolation result: cl=0.05, mu={}'.format(res))
        return res

    d = bsvj.AttrDict()
    d['twosigma_down'] = interpolate(cls.s_exp[0.975])
    d['onesigma_down'] = interpolate(cls.s_exp[0.84])
    d['expected'] = interpolate(cls.s_exp[0.5])
    d['onesigma_up'] = interpolate(cls.s_exp[0.16])
    d['twosigma_up'] = interpolate(cls.s_exp[0.025])
    d['observed'] = interpolate(cls.s)
    return d


def safe_divide(a, b):
    return np.divide(a, b, out=np.zeros_like(a), where=b!=0)

@is_script
def cls(args):
    if not isinstance(args, bsvj.AttrDict):
        parser = quickplot_parser('cls.png')
        parser.add_argument('observed', type=str)
        parser.add_argument('asimov', type=str)
        args = parser.parse_args(args)

    obs = bsvj.get_arrays(args.observed)
    asimov = bsvj.get_arrays(args.asimov)

    cls = get_cls(obs, asimov)
    limit = interpolate_95cl_limit(cls)

    fig = plt.figure(figsize=(14,14))
    ax = fig.gca()

    mu = cls.obs.mu
    ax.plot([mu[0], mu[-1]], [.05, .05], label='95%', c='purple')
    ax.plot(mu, cls.s, label='s', c='black')
    ax.plot(mu, cls.b, label='b', c='blue')
    ax.plot(mu, cls.sb, label='sb', c='red')
    ax.plot([cls.obs.mu_best, cls.obs.mu_best], [0., 1.05], c='xkcd:cobalt', label='bestfit', alpha=.3)

    # Expected
    ax.fill_between(mu, cls.s_exp[0.975], cls.s_exp[0.84], color="yellow", alpha=0.25)
    ax.fill_between(mu, cls.s_exp[0.84], cls.s_exp[0.16], color="green", alpha=0.25)
    ax.fill_between(mu, cls.s_exp[0.16], cls.s_exp[0.025], color="yellow", alpha=0.25)
    ax.plot(mu, cls.s_exp[0.5], c='black', linestyle='--', label='s_exp')
    
    # Limit points
    s = 45
    ax.scatter([limit.twosigma_down, limit.twosigma_up], [.05, .05], c='xkcd:dark yellow', s=s)
    ax.scatter([limit.onesigma_down, limit.onesigma_up], [.05, .05], c='green', s=s)
    ax.scatter([limit.expected, limit.observed], [.05, .05], c='black', s=s)

    ax.legend()
    ax.set_xlim(0.)
    ax.set_ylim(0., 1.05)
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel('CL')

    plt.savefig(args.outfile, bbox_inches='tight')
    if not(args.batch) and cmd_exists('imgcat'): os.system('imgcat ' + args.outfile)


@is_script
def brazil(args):
    if not isinstance(args, bsvj.AttrDict):
        parser = quickplot_parser('brazil.png')
        parser.add_argument('rootfiles', type=str, nargs='+')
        parser.add_argument('--xmin', type=float)
        parser.add_argument('--xmax', type=float)
        parser.add_argument('--ymin', type=float)
        parser.add_argument('--ymax', type=float)
        args = parser.parse_args(args)
    rootfiles = args.rootfiles

    get_mz = lambda rootfile: int(re.search(r'mz(\d+)', osp.basename(rootfile)).group(1))
    rootfiles.sort(key=get_mz)
    obs_rootfiles = [f for f in args.rootfiles if 'Observed' in osp.basename(f)]
    asimov_rootfiles = [f for f in args.rootfiles if 'Asimov' in osp.basename(f)]

    assert [get_mz(f) for f in obs_rootfiles] == [get_mz(f) for f in asimov_rootfiles]

    points = []
    for obs_rootfile, asimov_rootfile in zip(obs_rootfiles, asimov_rootfiles):
        mz = get_mz(obs_rootfile)
        assert mz == get_mz(asimov_rootfile)
        obs = bsvj.get_arrays(obs_rootfile)
        asimov = bsvj.get_arrays(asimov_rootfile)
        try:
            cls = get_cls(obs, asimov)
            limit = interpolate_95cl_limit(cls)
            points.append(bsvj.AttrDict(
                mz = mz,
                limit = limit,
                cls = cls
                ))
        except Exception:
            bsvj.logger.error('Error for mz {}:'.format(mz))
            traceback.print_exc()

    print(
        '{:<5s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} | {:>8s}'
        .format(
            'mz',
            '2s down',
            '1s down',
            'exp',
            '1s up',
            '2s up',
            'obs'
            )
        )
    for p in points:
        print(
            '{:<5.0f} {:+8.3f} {:+8.3f} {:+8.3f} {:+8.3f} {:+8.3f} | {:+8.3f}'
            .format(
                p.mz,
                p.limit.twosigma_down,
                p.limit.onesigma_down,
                p.limit.expected,
                p.limit.onesigma_up,
                p.limit.twosigma_up,
                p.limit.observed
                )
            )

    fig = plt.figure(figsize=(12,10))
    ax = fig.gca()
    
    mzs = [p.mz for p in points]

    ax.fill_between(
        mzs,
        [p.limit.twosigma_down for p in points],
        [p.limit.twosigma_up for p in points],
        color="yellow"
        )
    ax.fill_between(
        mzs,
        [p.limit.onesigma_down for p in points],
        [p.limit.onesigma_up for p in points],
        color="green"
        )    
    ax.plot(mzs, [p.limit.expected for p in points], c='black', linestyle='--', label='Exp')
    ax.plot(mzs, [p.limit.observed for p in points], c='black', linestyle='-', label='Obs')

    ax.set_xlabel(r'$m_{Z\prime}$ (GeV)')
    ax.set_ylabel(r'$\mu$')
    if hasattr(args, 'xmax') and args.xmax: ax.set_xlim(right=args.xmax)
    if hasattr(args, 'xmin') and args.xmin: ax.set_xlim(left=args.xmin)
    if hasattr(args, 'ymax') and args.ymax: ax.set_ylim(top=args.ymax)
    if hasattr(args, 'ymin') and args.ymin: ax.set_ylim(bottom=args.ymin)
    ax.legend()


    plt.savefig(args.outfile, bbox_inches='tight')
    if not(args.batch) and cmd_exists('imgcat'): os.system('imgcat ' + args.outfile)



@is_script
def allplots(args):
    if not isinstance(args, bsvj.AttrDict):
        parser = quickplot_parser()
        parser.add_argument('rootfiles', type=str, nargs='+')
        args = parser.parse_args(args)

    outdir = strftime('plots_%b%d')
    if not osp.isdir(outdir): os.makedirs(outdir)

    d = namespace_to_attrdict(args)
    d.batch = True

    for rootfile in args.rootfiles:
        mtdist(bsvj.AttrDict(
            d,
            rootfile=rootfile,
            outfile=osp.join(outdir, 'mtdist_{}.png'.format(name_from_combine_rootfile(rootfile)))
            ))

    brazil(bsvj.AttrDict(d, outfile=osp.join(outdir, 'brazil.png')))






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('script', type=str, choices=list(_scripts.keys()))
    parser.add_argument('-v', '--verbose', action='store_true')
    args, remaining_args = parser.parse_known_args()
    if args.verbose: bsvj.debug()
    r = _scripts[args.script](remaining_args)