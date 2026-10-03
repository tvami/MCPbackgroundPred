# Skim of the high-x tracks (x = -log10 probQ >= 3.4, |eta| < 1) from ntuples v1 for the x >= 5 category studies.
# Kept: JVR (JetMET group, pfMET < 100) pass and fail; SR (JetMET pfMET >= 100, Tau-only) FAIL only, SR pass stays blind.
# usage: python3 skim_hix.py <chunk> <nchunks> <outdir>
import glob, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
k, n, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
COLS = ['run', 'lumi', 'event', 'pt', 'eta', 'phi', 'pfMET', 'passMET', 'passJet', 'passTau', 'probQ_pixel', 'probXY_pixel',
        'probQ_pixelNoL1', 'nPixQFloor', 'nPixHitsUsed', 'pixSizeXresidual', 'pixSizeXpred', 'pixClSizeX', 'pixClSizeY',
        'pixClCharge', 'pixClSizeMax', 'ih_pixel', 'dedxPixel_builtin', 'nValidPixelHits', 'nPixXYvalid', 'nPixClusters',
        'nonL1PixHits', 'normChi2', 'ptError', 'validFraction', 'nTrackerLayers', 'caloEmEnergy', 'caloHadEnergy',
        'pcCaloFrac', 'charge', 'x', 'pas', 'reg']
files = sorted(glob.glob(NT + '/*/*/*/*/*.root'))[k::n]
ch = ROOT.TChain('MCPAnalyzer/tracks')
for f in files: ch.Add(f)
jm = '(passMET || passJet) && TString(fname).Contains("/JetMET")'
df = ROOT.RDataFrame(ch).Filter('highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1')
df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5f : -1.f)').Filter('x >= 3.4')
df = df.Define('pas', 'pixSizeXresidual > 1.25').DefinePerSample('fname', 'rdfsampleinfo_.AsString()')
# reg: 0 = JVR, 1 = SR JetMET, 2 = SR Tau-only, -1 = other (Tau PD events of the JetMET group: duplicates)
df = df.Define('reg', '(%s) ? (pfMET < 100 ? 0 : 1) : ((!(passMET || passJet) && passTau && TString(fname).Contains("/Tau/")) ? 2 : -1)' % jm)
df = df.Filter('reg == 0 || (reg > 0 && !pas)')
os.makedirs(OUT, exist_ok=True)
df.Snapshot('t', '%s/hix_%03d.root' % (OUT, k), COLS)
print('chunk %d: %d files done' % (k, len(files)))
