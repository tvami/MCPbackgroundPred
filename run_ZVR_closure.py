# 2DAlphabet closure in the Z->mumu validation region (ZVR): Z probes carry no signal,
# so the pass region above SIGSTART is fit-masked, predicted from the sideband, and compared.
# usage: python3 run_ZVR_closure.py <workingArea>
import os, sys
from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction

workingArea = sys.argv[1]
configJSON = os.environ.get("ZVR_CONFIG", "config_BinningZv1_InputZVRv1_ZVR_closure.json")
SIGNAL = "Signal_S0p5_M1000_Q30_ZVR"
BKG = 'CMS_MCP_Background'

# transfer functions in the normalized x in [0, 1]; R_P/F rises steeply with probQ
_rpf_options = {
    '0x0':   {'form': '@0', 'constraints': {0: {'MIN': 0, 'MAX': 100, 'NOMINAL': 5}}},
    '1x0':   {'form': '@0+@1*x', 'constraints': {0: {'MIN': 0, 'MAX': 100, 'NOMINAL': 5}, 1: {'MIN': 0, 'MAX': 3000, 'NOMINAL': 500}}},
    'expo':  {'form': 'exp(@0+@1*x)', 'constraints': {0: {'MIN': 0, 'MAX': 6, 'NOMINAL': 2}, 1: {'MIN': 0, 'MAX': 30, 'NOMINAL': 5}}},
    'expo2': {'form': 'exp(@0+@1*x+@2*x*x)', 'constraints': {0: {'MIN': 0, 'MAX': 6, 'NOMINAL': 2}, 1: {'MIN': -30, 'MAX': 30, 'NOMINAL': 5},
                                                            2: {'MIN': -30, 'MAX': 30, 'NOMINAL': 0}}},
    'expo3': {'form': 'exp(@0+@1*x+@2*x*x+@3*x*x*x)', 'constraints': {0: {'MIN': 0, 'MAX': 6, 'NOMINAL': 2}, 1: {'MIN': -30, 'MAX': 30, 'NOMINAL': 5},
                                                                     2: {'MIN': -30, 'MAX': 30, 'NOMINAL': 0}, 3: {'MIN': -30, 'MAX': 30, 'NOMINAL': 0}}},
}

def _select_signal(row, args):
    signame, tf = args
    if row.process_type == 'SIGNAL':
        return signame in row.process
    if 'Background' in row.process:
        return row.process in (BKG + '_' + tf, BKG)
    return True

def make_workspace():
    twoD = TwoDAlphabet(workingArea, configJSON, loadPrevious=False)
    bkg_hists = twoD.InitQCDHists()
    for f in [r for r in twoD.ledger.GetRegions() if 'fail' in r]:
        p = f.replace('fail', 'pass')
        binning_f, _ = twoD.GetBinningFor(f)
        fail_name = BKG + '_' + f
        bkg_f = BinnedDistribution(fail_name, bkg_hists[f], binning_f, constant=False)
        twoD.AddAlphaObj(BKG, f, bkg_f)
        for opt_name, opt in _rpf_options.items():
            rpf = ParametricFunction(fail_name.replace('fail', 'rpf') + '_' + opt_name, binning_f, opt['form'], opt['constraints'])
            bkg_p = bkg_f.Multiply(fail_name.replace('fail', 'pass') + '_' + opt_name, rpf)
            twoD.AddAlphaObj(BKG + '_' + opt_name, p, bkg_p, title='Background')
    twoD.Save()

def fit_and_plot(tf, rMax=10):
    twoD = TwoDAlphabet(workingArea, '{}/runConfig.json'.format(workingArea), loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, SIGNAL, tf)
    area = '{}-{}_area'.format(SIGNAL, tf)
    twoD.MakeCard(subset, area)
    twoD.MLfit(area, rMin=0, rMax=rMax, verbosity=1, defMinStrat=int(os.environ.get('ZVR_STRAT', '0')), extra='--robustHesse 1')
    try:  # the s+b fit is irrelevant for the closure; its plots may fail
        import json
        lj = os.path.join(os.path.dirname(os.path.abspath(configJSON)), json.load(open(configJSON))['GLOBAL']['path'], 'lumi_shown.json')
        lumi = r'%.1f $fb^{-1}$ (13.6 TeV)' % json.load(open(lj))['lumi_shown_fb'] if os.path.exists(lj) else '(13.6 TeV)'
        twoD.StdPlots(area, subset, lumiText=lumi,
                      pf_slice_str={'fail': 'sizeX residual #leq 1.25', 'pass': 'sizeX residual > 1.25'}, units='')
        plot.plot_transfer_funcs(workingArea, area)
    except Exception as e:
        print('plotting failed for %s: %s' % (tf, e))
    closure(area)

def closure(area):
    """predicted (b-only postfit) vs observed per window bin of the pass region"""
    import ROOT
    f = ROOT.TFile('%s/%s/postfitshapes_b.root' % (workingArea, area))
    tot_p = tot_o = 0.
    for reg in ['pass_SIG', 'pass_HIGH']:
        d = f.Get(reg + '_postfit')
        if not d: continue
        b, o = d.Get('TotalBkg').ProjectionX(), d.Get('data_obs').ProjectionX()
        for i in range(1, b.GetNbinsX() + 1):
            pb, ob = b.GetBinContent(i), o.GetBinContent(i)
            tot_p += pb; tot_o += ob
            print('CLOSURE %s %-9s x [%.2f,%.2f) predicted %9.1f +- %7.1f  observed %6.0f  ratio %.2f' % (
                area, reg, b.GetBinLowEdge(i), b.GetBinLowEdge(i + 1), pb, b.GetBinError(i), ob, pb / max(ob, 1)))
    print('CLOSURE %s window total predicted %.1f observed %.0f ratio %.2f' % (area, tot_p, tot_o, tot_p / max(tot_o, 1)))

def _n_fit_bins(twoD):
    b = twoD.binnings['default']
    nx = len(b.xbinList) - 1
    nwin = len(b.xbinByCat['SIG']) - 1 + len(b.xbinByCat['HIGH']) - 1
    return 2 * nx - nwin  # fail everywhere + pass sideband

def ftest(tf1, tf2):
    """F-test on the saturated GoF of the data, fitted bins only"""
    import ROOT
    from TwoDAlphabet.helpers import cd, execute_cmd
    twoD = TwoDAlphabet(workingArea, '{}/runConfig.json'.format(workingArea), loadPrevious=True)
    chi, npar = {}, {}
    for tf in (tf1, tf2):
        area = '{}-{}_area'.format(SIGNAL, tf)
        pars = twoD.ledger.select(_select_signal, SIGNAL, tf).alphaParams
        npar[tf] = len(pars[pars['name'].str.contains('rpf')].index)
        with cd('{}/{}'.format(workingArea, area)):
            if not os.path.exists('higgsCombine_gof_data.GoodnessOfFit.mH120.root'):
                execute_cmd('text2workspace.py -b card.txt -o gofws.root --channel-masks --X-no-jmax')
                execute_cmd('combine -M GoodnessOfFit -d gofws.root --algo=saturated -n _gof_data '
                            '--setParameters mask_pass_SIG=1,mask_pass_HIGH=1')
            tf_ = ROOT.TFile('higgsCombine_gof_data.GoodnessOfFit.mH120.root'); t = tf_.Get('limit'); t.GetEntry(0)
            chi[tf] = t.limit
    n = _n_fit_bins(twoD); p1, p2 = npar[tf1], npar[tf2]
    F = (chi[tf1] - chi[tf2]) / max(p2 - p1, 1) / (chi[tf2] / (n - p2))
    pval = 1. - ROOT.Math.fdistribution_cdf(F, p2 - p1, n - p2) if F > 0 else 1.
    print('FTEST %s (%d par, chi2 %.1f) vs %s (%d par, chi2 %.1f), n=%d: F=%.2f p=%.3g' % (tf1, p1, chi[tf1], tf2, p2, chi[tf2], n, F, pval))

def gof(tf, ntoys=300):
    """saturated GoF with local toys (blinding handled by 2DAlphabet)"""
    twoD = TwoDAlphabet(workingArea, '{}/runConfig.json'.format(workingArea), loadPrevious=True)
    area = '{}-{}_area'.format(SIGNAL, tf)
    twoD.GoodnessOfFit(area, ntoys=ntoys, freezeSignal=0, condor=False)
    plot.plot_gof(workingArea, area, condor=False)
    res = open('{}/{}/gof_results.txt'.format(workingArea, area)).read().replace('\n', ' ')
    print('GOF %s %s' % (area, res))

if __name__ == '__main__':
    make_workspace()
    for tf in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['expo', 'expo2', '1x0']):
        fit_and_plot(tf)
