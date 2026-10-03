# 2DAlphabet inputs per region from ntuples v1 (one chunk of files): JVR (JetMET group, pfMET < 100), SR_JetMET (pfMET >= 100),
# SR_Tau (Tau-only), each without and with the high-side requirement ih_pixel > 4 MeV/cm. SR pass is filled only below x = 3.4.
# usage: python3 make_hists_regions.py <chunk> <nchunks> <outdir>; then hadd per region
import glob, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
k, n, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
ch = ROOT.TChain('MCPAnalyzer/tracks')
for f in sorted(glob.glob(NT + '/*/*/*/*/*.root'))[k::n]: ch.Add(f)
df = ROOT.RDataFrame(ch).Filter('highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1')
df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5f : -1.f)').Filter('x >= 0')
df = df.Define('xc', 'x < 13.f ? x : 12.999f').Define('ae', 'abs(eta)').DefinePerSample('fname', 'rdfsampleinfo_.AsString()')
jm = '(passMET || passJet) && TString(fname).Contains("/JetMET")'
df = df.Define('reg', '(%s) ? (pfMET < 100 ? 0 : 1) : ((!(passMET || passJet) && passTau && TString(fname).Contains("/Tau/")) ? 2 : -1)' % jm)
df = df.Filter('reg >= 0')
hs = {}
for reg, rname in [(0, 'JVR'), (1, 'SRJetMET'), (2, 'SRTau')]:
    for ih, iname in [('1', ''), ('ih_pixel > 4', '_ih4')]:
        d = df.Filter('reg == %d && %s' % (reg, ih))
        for p, cut in [('pass', 'pixSizeXresidual > 1.25'), ('fail', 'pixSizeXresidual <= 1.25')]:
            if reg > 0 and p == 'pass': cut += ' && x < 3.4'  # SR stays blind
            hs[(rname + iname, p)] = d.Filter(cut).Histo2D(('h%s' % p, '', 260, 0, 13, 20, 0, 1), 'xc', 'ae')
ROOT.RDF.RunGraphs(list(hs.values()))
for (r, p), h in hs.items():
    os.makedirs('%s/%s' % (OUT, r), exist_ok=True)
    fo = ROOT.TFile('%s/%s/part_%s_%03d.root' % (OUT, r, p, k), 'RECREATE'); h.GetValue().Write('h' + p); fo.Close()
print('chunk %d done' % k)
