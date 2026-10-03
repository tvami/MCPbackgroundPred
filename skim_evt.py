# Track skim for the event-level analysis: every selected track (HP, DeDxHitInfo, |eta| < 1, sizeX defined) with its
# event id, region and quality variables, so the quality cuts and the one-candidate-per-event choice are applied offline.
# SR pass tracks are kept only to choose the candidate; they must never be histogrammed at x >= 3.4 (blinding).
# usage: python3 skim_evt.py <chunk> <nchunks> <outdir>
import glob, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
k, n, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
COLS = ['run', 'lumi', 'event', 'pt', 'eta', 'pfMET', 'probQ_pixel', 'nPixQFloor', 'pixSizeXresidual', 'pixClCharge',
        'ih_pixel', 'nValidPixelHits', 'nTrackerLayers', 'validFraction', 'normChi2', 'ptError', 'nStrip',
        'caloEmEnergy', 'caloHadEnergy', 'x', 'pas', 'reg']
ch = ROOT.TChain('MCPAnalyzer/tracks')
for f in sorted(glob.glob(NT + '/*/*/*/*/*.root'))[k::n]: ch.Add(f)
df = ROOT.RDataFrame(ch).Filter('highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1')
df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5f : -1.f)').Filter('x >= 0')
df = df.Define('pas', 'pixSizeXresidual > 1.25').DefinePerSample('fname', 'rdfsampleinfo_.AsString()')
jm = '(passMET || passJet) && TString(fname).Contains("/JetMET")'
df = df.Define('reg', '(%s) ? (pfMET < 100 ? 0 : 1) : ((!(passMET || passJet) && passTau && TString(fname).Contains("/Tau/")) ? 2 : -1)' % jm)
df = df.Filter('reg >= 0')
os.makedirs(OUT, exist_ok=True)
df.Snapshot('t', '%s/evt_%03d.root' % (OUT, k), COLS)
print('chunk %d: done' % k)
