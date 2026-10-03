# 2DAlphabet inputs for the Z->mumu closure (ZVR): pass/fail = sizeX residual > 1.25,
# x = -log10(probQ_pixel) (saturated tracks at 12.5), y = |eta| (one bin in the config).
# usage: python3 make_hists_ZVR.py <outdir>
import glob, json, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
OUT = sys.argv[1] if len(sys.argv) > 1 else 'histograms_for_2DAlphabet_ZVRv1'
SD = os.path.expanduser('~/HSCP-MCP/study_data')
SX, XSAT = 1.25, 12.5
NX, XMAX, NY = 260, 13.0, 20  # fine binning: 0.05 in x, 0.05 in |eta|
PROBE = ('tpOS && tpMass>81 && tpMass<101 && muMatched && highPurity && hasDeDx && pixSizeXresidual>-90'
         ' && !(passMET||passJet) && abs(eta)<1')
SIG = '(passMET||passJet||passTau) && highPurity && hasDeDx && pixSizeXresidual>-90 && abs(eta)<1 && genMatched'

def done(f):
    lg = os.path.join(os.path.dirname(f), 'cmsrun_' + os.path.basename(f)).replace('.root', '.log')
    return not os.path.exists(lg) or 'exit 0' in open(lg).read()[-200:]

def fill(files, sel, satcol, name, scale_to=None):
    ch = ROOT.TChain('MCPAnalyzer/tracks')
    for f in files: ch.Add(f)
    d = {k: np.asarray(v, float) for k, v in ROOT.RDataFrame(ch).Filter(sel).AsNumpy(
        ['probQ_pixel', satcol, 'pixSizeXresidual', 'eta']).items()}
    pq = d['probQ_pixel']
    x = np.where(pq > 0, -np.log10(np.clip(pq, 1e-30, 1)), np.nan)
    x = np.where((pq < 0) & (d[satcol] > 0), XSAT, x)
    ok = ~np.isnan(x)
    x, sx, ae = np.clip(x[ok], 0, XMAX - 1e-6), d['pixSizeXresidual'][ok], np.abs(d['eta'][ok])
    hp = ROOT.TH2D('hpass', '', NX, 0, XMAX, NY, 0, 1); hf = ROOT.TH2D('hfail', '', NX, 0, XMAX, NY, 0, 1)
    for h, m in [(hp, sx > SX), (hf, sx <= SX)]:
        h.Sumw2()
        for a, b in zip(x[m], ae[m]): h.Fill(a, b)
    if scale_to:
        s = scale_to / (hp.Integral() + hf.Integral()); hp.Scale(s); hf.Scale(s)
        for h in (hp, hf):  # empty-bin sentinel for templates: combine rejects zero-norm channels
            for i in range(1, NX + 1):
                for j in range(1, NY + 1):
                    if h.GetBinContent(i, j) <= 0: h.SetBinContent(i, j, 1e-8)
    fo = ROOT.TFile('%s/MCP_%s.root' % (OUT, name), 'RECREATE'); hp.Write(); hf.Write(); fo.Close()
    print('%-28s pass %10.1f  fail %10.1f  (pass x>=3.4: %.1f, fail x>=3.4: %.1f)' % (
        name, hp.Integral(), hf.Integral(), hp.Integral(hp.GetXaxis().FindBin(3.401), NX, 1, NY),
        hf.Integral(hf.GetXaxis().FindBin(3.401), NX, 1, NY)))

os.makedirs(OUT, exist_ok=True)

def shown_lumi(files):  # 109 fb^-1 x events read / all events of the 2024 ZMu skim (both halves)
    L = json.load(open(SD + '/zmu/lumi_ledger.json'))
    n = sum(e for f in files for _, e in L['jobs'].get(os.path.basename(f)[:-5], []))
    return L['lumi_total_fb'] * n / L['total_events']
dfiles = [f for f in glob.glob(SD + '/zmu/output/z*_*.root') if done(f) and not os.path.basename(f).startswith('dyaod')]
fill(dfiles, PROBE, 'nPixQFloor', 'Data_ZVR')
json.dump({'lumi_shown_fb': shown_lumi(dfiles), 'n_files': len(dfiles)}, open('%s/lumi_shown.json' % OUT, 'w'))
print('data shown: %.2f fb^-1 (%d files)' % (shown_lumi(dfiles), len(dfiles)))
fill([f for f in glob.glob(SD + '/zmu/output/dyaod_*.root') if done(f)], PROBE, 'nPixQFloor', 'DYMC_ZVR')
fill([SD + '/acc/output/acc_S0p5_M1000_Q30.root'], SIG, 'nPixHitsUsed', 'Signal_S0p5_M1000_Q30_ZVR', scale_to=50.)
