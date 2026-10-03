# 2DAlphabet signal templates over the |Q| >= 8e grid (central 2024 MiniAODv6, DY production):
# pass/fail = sizeX residual > 1.25, x = -log10(probQ_pixel) (saturated at 12.5), y = |eta|,
# normalized to sigma(DY, Run 3) x 109 fb^-1 (yield in tracks). Also a summary CSV and efficiency maps.
# usage: python3 make_hists_SR_signal.py [outdir]
import csv, glob, os, re, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
OUT = sys.argv[1] if len(sys.argv) > 1 else 'histograms_for_2DAlphabet_SRv1'
SD = os.path.expanduser('~/HSCP-MCP/study_data')
XS = os.path.expanduser('~/HSCP-MCP/AN/AN-26-065/Figures/DataAndSignal')
LUMI_PB = 109.0e3
ACC = os.environ.get('ACCDIR', SD + '/acc/output')  # signal ntuples
SAT = os.environ.get('SATCOL', 'nPixHitsUsed')  # saturated: probQ < 0 and SAT > 0 (nPixQFloor as in data)
IHMIN = float(os.environ.get('IHMIN', '-1e9'))  # high-side requirement on pixel Ih [MeV/cm], default none
SX, XSAT, XWIN, XEND = 1.25, 12.5, 3.4, 5.0
NX, XMAX, NY = 260, 13.0, 20
SIG = '(passMET||passJet||passTau) && highPurity && hasDeDx && pixSizeXresidual>-90 && abs(eta)<1 && genMatched && ih_pixel>%g' % IHMIN

def xsec_table(fn):  # (mass, charge in e/3) -> pb
    return {(int(float(r['mass'])), int(float(r['charge']))): float(r['xsec_pb']) for r in csv.DictReader(open(fn))}
XSEC = {'S0p5': xsec_table(XS + '/heco_cross_sections_run3.csv'), 'S0': xsec_table(XS + '/heco_cross_sections_spin0.csv')}

def templates(f, tag, norm):
    ch = ROOT.TChain('MCPAnalyzer/tracks'); ch.Add(f)
    d = {k: np.asarray(v, float) for k, v in ROOT.RDataFrame(ch).Filter(SIG).AsNumpy(
        ['probQ_pixel', SAT, 'pixSizeXresidual', 'eta']).items()}
    pq = d['probQ_pixel']
    x = np.where(pq > 0, -np.log10(np.clip(pq, 1e-30, 1)), np.nan)
    x = np.where((pq < 0) & (d[SAT] > 0), XSAT, x)
    ok = ~np.isnan(x)
    x, sx, ae = np.clip(x[ok], 0, XMAX - 1e-6), d['pixSizeXresidual'][ok], np.abs(d['eta'][ok])
    hp = ROOT.TH2D('hpass', '', NX, 0, XMAX, NY, 0, 1); hf = ROOT.TH2D('hfail', '', NX, 0, XMAX, NY, 0, 1)
    for h, m in [(hp, sx > SX), (hf, sx <= SX)]:
        h.Sumw2()
        for a, b in zip(x[m], ae[m]): h.Fill(a, b)
        h.Scale(norm)
    out = {}
    for k, h in [('pass', hp), ('fail', hf)]:  # yields in the fit window and the x > 5 category
        bw, be = h.GetXaxis().FindBin(XWIN + 1e-6), h.GetXaxis().FindBin(XEND - 1e-6)
        out[k + '_win'] = h.Integral(bw, be, 1, NY); out[k + '_hi'] = h.Integral(be + 1, NX, 1, NY)
        out[k + '_all'] = h.Integral()
    for h in (hp, hf):  # empty-bin sentinel: combine rejects zero-norm channels
        for i in range(1, NX + 1):
            for j in range(1, NY + 1):
                if h.GetBinContent(i, j) <= 0: h.SetBinContent(i, j, 1e-8)
    fo = ROOT.TFile('%s/MCP_Signal_%s.root' % (OUT, tag), 'RECREATE'); hp.Write(); hf.Write(); fo.Close()
    return out

os.makedirs(OUT, exist_ok=True)
rows = []
for f in sorted(glob.glob(ACC + '/acc*_*.root')):
    m = re.search(r'acc2?_(S0p5|S0)_M(\d+)_Q(\d+)\.root', f)
    sp, M, Q3 = m.group(1), int(m.group(2)), int(m.group(3))
    if Q3 < 24: continue
    tf = ROOT.TFile(f); t = tf.Get('MCPAnalyzer/events'); nev = t.GetEntries() if t else 0; tf.Close()
    if nev < 1000: continue  # low-stat points left out, as in the acceptance maps
    xs = XSEC[sp].get((M, Q3))
    if xs is None:  # DY xsec scales as Q^2: take the nearest charge at this mass
        qs = [q for m_, q in XSEC[sp] if m_ == M]
        if not qs: print('no xsec for', sp, M, Q3); continue
        qn = min(qs, key=lambda q: abs(q - Q3)); xs = XSEC[sp][(M, qn)] * (Q3 / qn) ** 2
    tag = '%s_M%d_Q%d' % (sp, M, Q3)
    y = templates(f, tag, xs * LUMI_PB / nev)
    rows.append(dict(spin='1/2' if sp == 'S0p5' else '0', M=M, Q=Q3 / 3., nev=nev, xsec_pb=xs, ntot=xs * LUMI_PB, **y))
    print('%-18s sigma %.3g pb  pass win %.3g  pass x>5 %.3g  fail win %.3g' % (tag, xs, y['pass_win'], y['pass_hi'], y['fail_win']), flush=True)
with open('%s/signal_yields.csv' % OUT, 'w') as fo:
    w = csv.DictWriter(fo, list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# efficiency maps (per produced event, tracks): pass in the window, pass x > 5, and the sum
sys.path.insert(0, SD); import hscp_style as S
MS = sorted({r['M'] for r in rows}); QS = sorted({r['Q'] for r in rows})
S.CMS.SetLumi(None, run='Simulation 2024'); S.CMS.SetEnergy(0, unit=''); S.CMS.SetExtraText('Work in Progress')
for key, title in [('pass_win', 'tracks in pass, 3.4 #leq x < 5'), ('pass_hi', 'tracks in pass, x #geq 5'),
                   ('pass_sr', 'tracks in pass, x #geq 3.4')]:
    for sp, stag in [('1/2', 'spin12'), ('0', 'spin0')]:
        h = ROOT.TH2D('h_%s_%s' % (key, stag), '', len(MS), 0, len(MS), len(QS), 0, len(QS))
        for r in rows:
            if r['spin'] != sp: continue
            v = (r['pass_win'] + r['pass_hi'] if key == 'pass_sr' else r[key]) / r['ntot']
            h.SetBinContent(MS.index(r['M']) + 1, QS.index(r['Q']) + 1, max(v, 1e-4))
        for i, m in enumerate(MS): h.GetXaxis().SetBinLabel(i + 1, str(m))
        for j, q in enumerate(QS): h.GetYaxis().SetBinLabel(j + 1, '%g' % q)
        h.SetMinimum(0); h.SetMaximum(max(1., h.GetMaximum())); S.keep.append(h)
        c = S.canvas('c_%s_%s' % (key, stag), 0, len(MS), 0, len(QS), 'M [GeV]', '|Q| [e]', z=True, iPos=0, yoff=1.1)
        c.SetRightMargin(0.2); c.SetBottomMargin(0.14)
        ROOT.gStyle.SetPaintTextFormat('.2f')
        h.GetXaxis().SetTitle('M [GeV]'); h.GetYaxis().SetTitle('|Q| [e]')
        h.GetZaxis().SetTitle('Spin-%s: %s per event' % (sp, title)); h.GetZaxis().SetTitleOffset(1.45)
        for ax in (h.GetXaxis(), h.GetYaxis()): ax.SetLabelSize(0.032); ax.SetTitleSize(0.045)
        h.GetXaxis().LabelsOption('v'); h.GetXaxis().SetTitleOffset(1.35); h.GetYaxis().SetTitleOffset(1.0)
        h.GetZaxis().SetTitleSize(0.035); h.GetZaxis().SetLabelSize(0.03); h.SetMarkerSize(0.85)
        h.Draw('COLZ TEXT')
        S.CMS.CMS_lumi(c, 0)
        S.save(c, '%s/sig_%s_%s' % (OUT, key, stag))
print('%d points' % len(rows))
