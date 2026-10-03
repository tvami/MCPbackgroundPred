# x >= 5 counting category (includes the saturated bin): background = fail(x>=5, SR) x R, R = pass/fail at x>=5
# measured in JVR (JetMET group, pfMET < 100). Checks R vs pfMET, era and pT in JVR, and compares the sideband
# pass/fail ratio (2 <= x < 3.4) between the JetMET and Tau-only streams. SR pass at x >= 3.4 is never read.
# usage: python3 xgt5_category.py [out.json]
import glob, json, sys
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.EnableImplicitMT(8)
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'xgt5_category.json'
TRK = 'highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1'
ERAS = [('C', 379415, 380255), ('D', 380255, 380948), ('E', 380948, 381944), ('F', 381944, 383780),
        ('G', 383780, 385814), ('H', 385814, 386409), ('I', 386409, 387200)]

def frame(pd, sel):
    ch = ROOT.TChain('MCPAnalyzer/tracks')
    for f in sorted(glob.glob(NT + '/%s/*/*/*/*.root' % pd)): ch.Add(f)
    df = ROOT.RDataFrame(ch).Filter(sel + ' && ' + TRK)
    df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5 : -1.)').Filter('x >= 0')
    return ch, df.Define('pas', 'pixSizeXresidual > 1.25')

def book(df, groups, blind):
    """per group: counts in the sideband [2, 3.4), window [3.4, 5), high x >= 5; pass only if not blind"""
    out = {}
    for g, cut in groups.items():
        d = df.Filter(cut)
        out[g] = {r + '_' + k: d.Filter('%s && %s' % (rc, kc)).Count()
                  for r, rc in [('sb', 'x >= 2 && x < 3.4'), ('win', 'x >= 3.4 && x < 5'), ('hi', 'x >= 5')]
                  for k, kc in [('fail', '!pas'), ('pass', 'pas')] if not (blind and r != 'sb' and k == 'pass')}
    return out

jv = {'JVR': 'pfMET < 100'}
jv.update({'JVR_met%d' % lo: 'pfMET >= %d && pfMET < %d' % (lo, lo + 25) for lo in (0, 25, 50, 75)})
jv.update({'JVR_era' + e: 'pfMET < 100 && run >= %d && run < %d' % (a, b) for e, a, b in ERAS})
jv.update({'JVR_pt%d' % lo: 'pfMET < 100 && pt >= %d && pt < %d' % (lo, hi) for lo, hi in ((0, 50), (50, 100), (100, 200), (200, 100000))})
c1, dj = frame('JetMET*', '(passMET || passJet)')
c2, dt = frame('Tau', '!(passMET || passJet) && passTau')
res = book(dj, jv, False)
res.update(book(dj, {'SR_JetMET': 'pfMET >= 100'}, True))
res.update(book(dt, {'SR_Tau': 'pfMET >= 0', 'SR_Tau_met100': 'pfMET < 100'}, True))
ROOT.RDF.RunGraphs([v for g in res.values() for v in g.values()])
res = {g: {k: int(v.GetValue()) for k, v in d.items()} for g, d in res.items()}
json.dump(res, open(OUT, 'w'), indent=1)

def ratio(p, f):
    return (p / f, p / f * (1. / max(p, 1) + 1. / max(f, 1)) ** 0.5) if f > 0 else (float('nan'), float('nan'))
print('%-14s %9s %9s %10s %7s %7s %9s %9s' % ('group', 'R_sb', '+-', 'fail_win', 'fail_hi', 'pass_hi', 'R_hi', '+-'))
for g, d in res.items():
    rs, es = ratio(d['sb_pass'], d['sb_fail'])
    rh, eh = ratio(d.get('hi_pass', 0), d['hi_fail']) if 'hi_pass' in d else (float('nan'), float('nan'))
    print('%-14s %9.1f %9.1f %10d %7d %7s %9.0f %9.0f' % (g, rs, es, d['win_fail'], d['hi_fail'], d.get('hi_pass', 'blind'), rh, eh))
R, eR = ratio(res['JVR']['hi_pass'], res['JVR']['hi_fail'])
for g in ('SR_JetMET', 'SR_Tau'):
    f = res[g]['hi_fail']
    print('prediction %s: fail(x>=5) = %d -> pass(x>=5) = %.0f +- %.0f (stat: R %.0f%%, fail %.0f%%)' % (
        g, f, f * R, f * R * ((eR / R) ** 2 + 1. / max(f, 1)) ** 0.5, 100 * eR / R, 100 / max(f, 1) ** 0.5))
