# x >= 5 category: signal vs background (JVR pass tracks, x >= 5) in candidate discriminants, from the high-x skim.
# Prints the background efficiency at fixed signal efficiency per variable and benchmark; draws normalized overlays.
# usage: python3 hix_study.py <skim dir> <plot dir>
import glob, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
SKIM, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
SD = os.path.expanduser('~/HSCP-MCP/study_data')
sys.path.insert(0, SD); import hscp_style as S
SIG = '(passMET||passJet||passTau) && highPurity && hasDeDx && pixSizeXresidual>1.25 && abs(eta)<1 && genMatched'
BENCH = [('S0p5_M1000_Q24', 'M = 1 TeV, |Q| = 8e', ROOT.kRed + 1), ('S0p5_M1000_Q48', 'M = 1 TeV, |Q| = 16e', ROOT.kOrange + 7),
         ('S0p5_M2000_Q90', 'M = 2 TeV, |Q| = 30e', ROOT.kViolet + 1)]
# name, expression (numpy on the column dict), bins, lo, hi, title, expected signal side (cuts try both)
VARS = [('clch', lambda d: d['pixClCharge'] / 1e3, 50, 0, 500, 'Mean pixel cluster charge [ke]', '>'),
        ('ih', lambda d: d['ih_pixel'], 50, 0, 50, 'I_{h} (pixel) [MeV/cm]', '>'),  # diagnostic only, never a selection cut
        ('sizex', lambda d: d['pixSizeXresidual'], 50, 1.25, 11.25, 'Pixel sizeX residual', '>'),
        ('sizey', lambda d: d['pixClSizeY'], 40, 0, 20, 'Mean pixel cluster sizeY', '>'),
        ('clmax', lambda d: d['pixClSizeMax'], 40, 0, 40, 'Largest pixel cluster size', '>'),
        ('fsat', lambda d: (d['x'] > 12) * 1., 2, 0, 2, 'Saturated (all pixel hits at the floor)', '>'),
        ('pt', lambda d: d['pt'], 50, 0, 1000, 'Track p_{T} [GeV]', '<'),
        ('relpterr', lambda d: d['ptError'] / d['pt'], 50, 0, 1, '#sigma_{p_{T}} / p_{T}', '>'),
        ('chi2', lambda d: d['normChi2'], 50, 0, 10, 'Track #chi^{2}/ndof', '>'),
        ('calo', lambda d: d['caloEmEnergy'] + d['caloHadEnergy'], 50, 0, 500, 'Calo energy near the track [GeV]', '<'),
        ('npix', lambda d: d['nValidPixelHits'] * 1., 8, 0, 8, 'Valid pixel hits', '>')]
COLS = ['pixClCharge', 'ih_pixel', 'pixSizeXresidual', 'pixClSizeY', 'pixClSizeMax', 'pt', 'ptError', 'normChi2',
        'caloEmEnergy', 'caloHadEnergy', 'nValidPixelHits', 'probQ_pixel']

def load(chain_name, files, sel, satcol):
    ch = ROOT.TChain(chain_name)
    for f in files: ch.Add(f)
    cols = COLS + ([satcol] if satcol else ['x'])
    d = {k: np.asarray(v, float) for k, v in ROOT.RDataFrame(ch).Filter(sel).AsNumpy(cols).items()}
    if satcol:  # signal: build x as in the templates
        pq = d['probQ_pixel']
        x = np.where(pq > 0, -np.log10(np.clip(pq, 1e-30, 1)), np.nan)
        d['x'] = np.where((pq < 0) & (d[satcol] > 0), 12.5, x)
        keep = d['x'] >= 5
        d = {k: v[keep] for k, v in d.items()}
    return d

bkg = load('t', glob.glob(SKIM + '/hix_*.root'), 'reg == 0 && pas && x >= 5', None)
srf = load('t', glob.glob(SKIM + '/hix_*.root'), 'reg > 0 && !pas && x >= 5', None)
sig = {k: load('MCPAnalyzer/tracks', [SD + '/acc/output/acc_%s.root' % k], SIG, 'nPixHitsUsed') for k, _, _ in BENCH}
print('background (JVR pass, x>=5) %d tracks; SR fail x>=5 %d; signal %s' % (len(bkg['x']), len(srf['x']), {k: len(v['x']) for k, v in sig.items()}))

S.lumi(None, run='2024'); S.CMS.SetEnergy(13.6)
print('%-10s' % 'var' + ''.join('%24s' % b[0] for b in BENCH) + '   (bkg eff at sig eff 0.9 / 0.7)')
for name, fn, nb, lo, hi, xt, side in VARS:
    vb = fn(bkg); row = '%-10s' % name
    for k, _, _ in BENCH:
        vs = fn(sig[k]); cells = []
        for es in (0.9, 0.7):  # best of the two cut directions
            eb = min((vb >= np.quantile(vs, 1 - es)).mean(), (vb <= np.quantile(vs, es)).mean())
            cells.append('%.3f' % eb)
        row += '%24s' % ' / '.join(cells)
    print(row)
    hs = []
    for lab, arr, col, sty in [('JetMET VR, pass, x #geq 5', vb, ROOT.kBlack, 'pe')] + \
                              [(t, fn(sig[k]), c, 'hist') for k, t, c in BENCH]:
        h = ROOT.TH1D('h_%s_%d' % (name, len(hs)), '', nb, lo, hi); h.Sumw2()
        for v in np.clip(arr, lo + 1e-9, hi - 1e-9): h.Fill(v)
        h.Scale(1. / max(h.Integral(), 1)); hs.append((h, lab, col, sty)); S.keep.append(h)
    c = S.canvas('c_' + name, lo, hi, 1e-4, 30, xt, 'Fraction of tracks', logy=True)
    next(o for o in c.GetListOfPrimitives() if o.InheritsFrom('TH1')).GetXaxis().SetNdivisions(505)
    leg = S.legend(0.52, 0.70, 0.93, 0.88, 0.028)
    for h, lab, col, sty in hs:
        if sty == 'pe': S.CMS.cmsDraw(h, 'pe', **S.DATA); leg.AddEntry(h, lab, 'pe')
        else: S.CMS.cmsDraw(h, 'hist', lcolor=col, lwidth=2, fstyle=0); leg.AddEntry(h, lab, 'l')
    S.note('|#eta| < 1, sizeX residual > 1.25', 0.255, 0.64, 0.03)
    S.save(c, '%s/hix_%s' % (OUT, name))
# do the SR fail tracks look like the JVR fail ones? (background composition check)
jf = load('t', glob.glob(SKIM + '/hix_*.root'), 'reg == 0 && !pas && x >= 5', None)
for name, fn, *_ in VARS:
    print('fail x>=5 median %-9s JVR %.3g  SR %.3g' % (name, np.median(fn(jf)), np.median(fn(srf))))
