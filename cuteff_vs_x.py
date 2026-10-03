# Efficiency of the x >= 5 category cuts (calo < 200 GeV, pixel Ih > 8) vs x = -log10 probQ, per region and pass/fail,
# from ntuples v1 (|eta| < 1, x >= 2). SR pass at x >= 3.4 is never counted. usage: python3 cuteff_vs_x.py <chunk> <n> <outdir>
import glob, json, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
k, n, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
XB = [2, 2.5, 3, 3.4, 4, 5, 13]
CUTS = {'c200i8': '(caloEmEnergy + caloHadEnergy) < 200 && ih_pixel > 8', 'c150i10': '(caloEmEnergy + caloHadEnergy) < 150 && ih_pixel > 10'}
ch = ROOT.TChain('MCPAnalyzer/tracks')
for f in sorted(glob.glob(NT + '/*/*/*/*/*.root'))[k::n]: ch.Add(f)
df = ROOT.RDataFrame(ch).Filter('highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1')
df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5f : -1.f)').Filter('x >= 2')
df = df.Define('pas', 'pixSizeXresidual > 1.25').DefinePerSample('fname', 'rdfsampleinfo_.AsString()')
jm = '(passMET || passJet) && TString(fname).Contains("/JetMET")'
df = df.Define('reg', '(%s) ? (pfMET < 100 ? 0 : 1) : ((!(passMET || passJet) && passTau && TString(fname).Contains("/Tau/")) ? 2 : -1)' % jm)
df = df.Filter('reg >= 0 && !(reg > 0 && pas && x >= 3.4)')
res = {}
for reg in (0, 1, 2):
    for p in (0, 1):
        for lo, hi in zip(XB[:-1], XB[1:]):
            if reg > 0 and p == 1 and lo >= 3.4: continue
            d = df.Filter('reg == %d && pas == %d && x >= %g && x < %g' % (reg, p, lo, hi))
            key = '%d_%d_%g' % (reg, p, lo)
            res[key] = {'all': d.Count()}
            res[key].update({c: d.Filter(e).Count() for c, e in CUTS.items()})
ROOT.RDF.RunGraphs([v for r in res.values() for v in r.values()])
os.makedirs(OUT, exist_ok=True)
json.dump({kk: {c: int(v.GetValue()) for c, v in r.items()} for kk, r in res.items()}, open('%s/cuteff_%03d.json' % (OUT, k), 'w'))
print('chunk %d done' % k)
