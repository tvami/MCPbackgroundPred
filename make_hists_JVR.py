# 2DAlphabet inputs for the JetMET validation region (JVR): JetMET PD events triggered by the JetMET group with
# pfMET < 100 GeV (signal <= 3% of triggered signal). Tracks: HP, DeDxHitInfo, |eta| < 1, sizeX defined.
# pass/fail = sizeX residual > 1.25, x = -log10(probQ_pixel) (saturated at 12.5).
import glob, json, os, sys
import ROOT
ROOT.gROOT.SetBatch(True); ROOT.EnableImplicitMT(8)
OUT = sys.argv[1] if len(sys.argv) > 1 else 'histograms_for_2DAlphabet_JVRv1'
NT = '/ceph/cms/store/user/tvami/HSCP_MCP/Ntuples_v1'
os.makedirs(OUT, exist_ok=True)
files = sorted(glob.glob(NT + '/JetMET*/*/*/*/*.root'))
ch = ROOT.TChain('MCPAnalyzer/tracks')
for f in files: ch.Add(f)
df = ROOT.RDataFrame(ch).Filter('(passMET || passJet) && pfMET < 100 && highPurity && hasDeDx && pixSizeXresidual > -90 && abs(eta) < 1')
df = df.Define('x', 'probQ_pixel > 0 ? -log10(probQ_pixel) : ((nPixQFloor > 0) ? 12.5 : -1.)').Filter('x >= 0') \
       .Define('xc', 'x < 13. ? x : 12.999').Define('ae', 'abs(eta)')
m = ('h', '', 260, 0, 13, 20, 0, 1)
hp = df.Filter('pixSizeXresidual > 1.25').Histo2D(('hpass',) + m[1:], 'xc', 'ae')
hf = df.Filter('pixSizeXresidual <= 1.25').Histo2D(('hfail',) + m[1:], 'xc', 'ae')
fo = ROOT.TFile('%s/MCP_Data_JVR.root' % OUT, 'RECREATE'); hp.GetValue().Write('hpass'); hf.GetValue().Write('hfail'); fo.Close()
# signal template: copy the ZVR one (same selection on the signal track)
os.system('cp histograms_for_2DAlphabet_ZVRv1/MCP_Signal_S0p5_M1000_Q30_ZVR.root %s/MCP_Signal_S0p5_M1000_Q30_ZVR.root' % OUT)
json.dump({'lumi_shown_fb': 109.0, 'n_files': len(files)}, open('%s/lumi_shown.json' % OUT, 'w'))
h1, h2 = hp.GetValue(), hf.GetValue()
print('files %d; pass %.0f fail %.0f; window (3.4-5): pass %.0f fail %.0f; x>5: pass %.0f fail %.0f' % (
    len(files), h1.Integral(), h2.Integral(), h1.Integral(69, 100, 1, 20), h2.Integral(69, 100, 1, 20),
    h1.Integral(101, 260, 1, 20), h2.Integral(101, 260, 1, 20)))
