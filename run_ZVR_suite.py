# ZVR closure suite: window bins from the fail statistics, data and DY MC as data_obs, x up to 5 and full range,
# all transfer functions, F-tests and local GoF toys. usage: python3 run_ZVR_suite.py <tag> [min_fail]
import json, os, subprocess, sys
import ROOT

TAG = sys.argv[1]
MINFAIL = float(sys.argv[2]) if len(sys.argv) > 2 else 10.
HDIR = os.environ.get('ZVR_HDIR', 'histograms_for_2DAlphabet_ZVRv1')
TFS = os.environ.get('ZVR_TFS', '0x0,1x0,expo,expo2').split(',')
XLO = float(os.environ.get('ZVR_XLO', '0'))
SIDE = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.4] if XLO == 0 else [round(XLO + 0.25 * i, 2) for i in range(int(round((3.4 - XLO) / 0.25)))] + [3.4]
FINE = 0.05

def window_bins(sample, xend):
    """window edges from 3.4: close a bin once it holds MINFAIL fail entries (last bin takes the rest)"""
    f = ROOT.TFile('%s/MCP_%s.root' % (HDIR, sample)); px = f.Get('hfail').ProjectionX()
    edges, acc, x = [3.4], 0., 3.4
    while x < xend - 1e-9:
        b = px.FindBin(x + 1e-6); acc += px.GetBinContent(b); x = round(x + FINE, 4)
        if acc >= MINFAIL and xend - x > 0.3: edges.append(x); acc = 0.
    if edges[-1] != xend: edges.append(xend)
    if len(edges) < 3:  # 2DA needs a non-empty HIGH: split the window
        edges = [3.4, round((3.4 + xend) / 2, 2), xend]
    return edges

def config(sample, xend):
    c = json.load(open('config_BinningZv1_InputZVRv1_ZVR_closure.json'))
    w = window_bins(sample, xend)
    c['BINNING']['default']['X']['BINS'] = SIDE + w[1:]
    c['BINNING']['default']['X']['SIGSTART'] = 3.4
    c['BINNING']['default']['X']['SIGEND'] = w[-2]
    c['PROCESSES']['data_obs']['ALIAS'] = sample
    c['GLOBAL']['path'] = './' + HDIR
    name = 'config_%s_%s_x%s.json' % (TAG, sample, str(xend).replace('.', 'p'))
    json.dump(c, open(name, 'w'), indent=1)
    return name, w

def run(cfg, area):
    env = dict(os.environ, ZVR_CONFIG=cfg)
    log = open('%s.log' % area, 'w')
    subprocess.call(['python3', '-c', """
import sys; sys.argv=['x','%s']
import run_ZVR_closure as m
m.make_workspace()
for tf in %r:
    try: m.fit_and_plot(tf)
    except Exception as e: print('FIT_FAILED', tf, e)
for a, b in [(t1, t2) for t1, t2 in zip(%r[:-1], %r[1:])]:
    try: m.ftest(a, b)
    except Exception as e: print('FTEST_FAILED', a, b, e)
for tf in %r:
    try: m.gof(tf)
    except Exception as e: print('GOF_FAILED', tf, e)
""" % (area, TFS, TFS, TFS, TFS)], env=env, stdout=log, stderr=subprocess.STDOUT)

summary = {}
for sample in os.environ.get('ZVR_SAMPLES', 'Data_ZVR,DYMC_ZVR').split(','):
    for xend in [float(v) for v in os.environ.get('ZVR_XEND', '5.0,13.0').split(',')]:
        cfg, w = config(sample, xend)
        area = 'rpf_%s_%s_x%s' % (TAG, sample, str(xend).replace('.', 'p'))
        os.system('rm -rf ' + area)
        print('running', area, 'window', w, flush=True)
        run(cfg, area)
        summary[area] = dict(window=w, closure=[l.strip() for l in open(area + '.log') if l.startswith(('CLOSURE', 'FTEST', 'GOF', 'FIT_FAILED'))])
json.dump(summary, open('zvr_suite_%s.json' % TAG, 'w'), indent=1)
for a, v in summary.items():
    print('==', a, v['window']); [print('  ', l) for l in v['closure']]
