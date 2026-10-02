# MCPbackgroundPred: 2DAlphabet background prediction for the HSCP-MCP search

Data-driven background estimate and statistical analysis for the search for multiply charged
particles (MCPs, HSCP-MCP) with Run-3 data. The background is predicted with the **Alphabet**
method as implemented in [2DAlphabet](https://github.com/ammitra/2DAlphabet) (vendored here), and
the fits and limits are done with `combine`.

The ntuples come from [MCPAnalyzer](https://github.com/tvami/MCPAnalyzer); this repository starts
from the histograms made out of them.

## The method in one paragraph

Tracks (events) are split by a **tagger** into a **pass** and a **fail** region. The background shape in the
fail region is taken from data, bin by bin, as free parameters. The pass region is predicted as

```
N_pass(x, y) = R_P/F(x, y) * N_fail(x, y)
```

where the transfer function `R_P/F` is a low-order polynomial in the two analysis variables, fitted to
data **outside the blinded signal window** together with the signal template. This is ABCD
generalized to a smooth, fitted ratio: it does not require the tagger to be independent of `x` and
`y`, it only requires the ratio to vary slowly, which the polynomial order (chosen by an F-test)
captures. Signal enters through its pass and fail templates with a common signal strength `r`.

## Planned MCP setup

The MCP search targets **|Q| >= 8e** (Q label `24` and up). In that range the MCP fires neither the
IsoMu24 nor the Mu50 trigger (efficiency below 0.003 over the grid), so the Muon PD is free to serve as
a control sample.

| Ingredient | MCP choice |
|---|---|
| Pass / fail | pixel **sizeX residual** (measured minus angle-predicted cluster size): pass above the working point (currently 1.25), fail below |
| Fit variable `x` | **`-log10(probQ_pixel)`**; MIPs fall steeply, the MCP sits at large values. Saturated tracks (probQ undefined, all pixel hits at the template floor, `nPixQFloor`) go in the last `x` bin |
| `y` (optional) | eta region (barrel / endcap), see below |
| Signal region | JetMET (+ Tau) triggered data, pass region, `x` inside `[SIGSTART, SIGEND)` (the probQ tail, e.g. above 3.4) blinded |
| Validation regions | the same fit in a JetMET sideband with negligible signal (low MET, or prescaled jet triggers), and in the `x` sideband below `SIGSTART` |
| MIP reference | Z->mumu tag-and-probe from the Run-3 `ZMu` RAW-RECO skim (`/Muon{0,1}/Run2024*-ZMu-*/RAW-RECO`): `R_P/F(x)` for MIPs, data/MC, stability |
| Signal | DY (and photon fusion) MCP grid, both spins, `M = 100-2400 GeV` |

So the fit predicts the probQ spectrum of sizeX-tagged tracks from the probQ spectrum of untagged
tracks: `N_pass(x) = R_P/F(x) * N_fail(x)`.

Why Alphabet and not a plain ABCD in the (probQ, sizeX) plane: the two variables are correlated for
MIPs (rho = 0.2-0.3; the joint tail rate over the product of the single rates is 1.2 in the barrel and
3.5 in the endcaps), so ABCD does not close. Here that correlation is exactly the `x` dependence of
`R_P/F`, which the polynomial fits, and the Z probes measure it for MIPs directly.

Things to settle before the first fit:
- the sizeX residual has a barrel/endcap offset for MIPs (median about 1.56 in the barrel, 0.9-1.1 in
  the endcaps), so with one cut `R_P/F` mostly measures the eta mix. Either use an eta-dependent
  working point (or a residual corrected per layer), or split the fit in eta with `y`;
- the binning of `x` at large probQ, where the MIP background runs out, and the treatment of the
  saturated bin.

## Set up

```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git clone git@github.com:tvami/MCPbackgroundPred.git .
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v10.0.1
cd ../..
scram b -j
```

### Create the python environment (only once)
```
python3 -m virtualenv twoD-env
source twoD-env/bin/activate
cd 2DAlphabet
python3 setup.py develop
cd ..
```

### Use the environment
```
cd /path/to/CMSSW_14_1_0_pre4/src
cmsenv
source twoD-env/bin/activate
```

**OS note (uaf):** this release is built for el8. On an el9 host run inside the el8 container, and
always pass the command after `--` (without it the first token is read as an image name):
```
cmssw-el8 -- bash -c 'cd /path/to/CMSSW_14_1_0_pre4/src && eval `scram runtime -sh` && source twoD-env/bin/activate && <command>'
```
All `combine` based commands below assume the el8 environment.

## Repository layout

| Path | What it is |
|---|---|
| `2DAlphabet/` | vendored 2DAlphabet, installed with `setup.py develop`, so edits under `2DAlphabet/TwoDAlphabet/` take effect immediately |
| `CombineHarvester/` | vendored CombineHarvester (used for impacts) |
| `HiggsAnalysis/CombinedLimit/` | `combine`, cloned at setup, not tracked |
| `twoD-env/` | the python virtualenv, created at setup, not tracked |

`twoDalphabetMod.py`, `alphawrapMod.py` and `binningMod.py` are experimental variants next to the
original modules. Run scripts should import the non-`Mod` modules
(`from TwoDAlphabet.twoDalphabet import ...`); check the imports before assuming which one is in use.

## Inputs

### Histograms
One ROOT file per process and region, holding 2D histograms of `(x, y)` for the pass and fail
regions, plus the systematic variations of the signal:

```
histograms_for_2DAlphabet_v<N>/
    MCP_Data_<REGION>.root          # data_obs, pass and fail
    MCP_Signal_<NAME>_<REGION>.root # one per signal point, nominal + up/down templates
    MCP_BkgMC_<REGION>.root         # optional: simulated background, for closure tests
```

with `<NAME>` such as `S0p5_M1000_Q30` (spin, mass in GeV, `Q` label = 3|Q|/e). Input versions
(`v<N>`) and binning versions are bumped independently and each one gets a line in this README
(see "Versioning").

### Config JSON
Each fit is driven by a JSON config (see `2DAlphabet/example_config.json` for the full schema):

- `GLOBAL.path`, `GLOBAL.FILE`, `GLOBAL.HIST`: where the histograms are and how they are named
  (`$process` and `$region` are substituted).
- `GLOBAL.SIGNAME`: the list of signal points.
- `PROCESSES`: `data_obs` (`TYPE: DATA`), the signal (`TYPE: SIGNAL`) and its `SYSTEMATICS`.
- `REGIONS`: the pass and fail regions and which processes enter them.
- `BINNING.<name>.X/Y`: `MIN`, `MAX`, `NBINS` (or explicit `BINS`), and on `X` the blinding window
  `SIGSTART`/`SIGEND`. The last bin edge must lie inside the input histogram range.
- `SYSTEMATICS`: `CODE 0` = lnN (`VAL`), shape systematics point to up/down templates.
- `OPTIONS`: blinding and plotting switches (see "Blinding").

Name configs as `config_Binningv<M>_Inputv<N>_<REGION>_<SIGNAL|BkgMC>.json`.

## Running

A run script takes the working-area name as its first argument and runs the chain for one region:

1. build the workspace and the cards from the config (`make_workspace`);
2. fit each transfer-function order (`R_P/F` polynomial, e.g. `0x0`, `1x0`, `1x1`, `2x0`);
3. make the post-fit, transfer-function and pull plots;
4. run the goodness of fit (toys on condor), harvest and plot it;
5. run the F-test between neighboring orders to pick the transfer function.

```
python3 run_<REGION>.py <workingArea>
```

Keep the config JSON **hardcoded in the run script** (`configJSON = ...`), so a working area always
corresponds to exactly one config; the command line only names the output directory.

Before any condor round, unit test the workspace build (no fit, no proxy). It catches out-of-range
bins and histogram-name mismatches:
```
python3 -c "import sys; sys.argv=['x','UNITTEST']; import run_SR as m; m.make_workspace()"
rm -rf UNITTEST
```

### Re-plotting without re-running fits
Run from `CMSSW_14_1_0_pre4/src` with the environment active; these read existing output only:
```
python3 -c "from TwoDAlphabet import plot; plot.plot_gof('<workingArea>','<signal>-<tf>_area', condor=True)"
python3 -c "from TwoDAlphabet import plot; plot.plot_transfer_funcs('<workingArea>','<signal>-<tf>_area')"
```
Use `condor=True` for condor toy tarballs (`*_gof_toys_output_*.tgz`), `condor=False` for one local
toy file.

### Condor
The GoF toys, the per-signal limits and the impacts run on condor (T2_US_UCSD, rhel8 image). A valid
grid proxy is needed (`voms-proxy-init -voms cms`). For many signal points use a templated config
(`..._InputTemplate_SR_Blind.json`) plus a single-signal run script, with one job per point and
transfer-function order.

## Blinding

- `OPTIONS.blindedFit: ["pass"]` removes the blinded part of the pass region from the **likelihood**
  (it is still built and can be plotted); `blindedPlots: ["pass"]` hides it in the **plots** only.
- `blindedFitSubregions` selects which `x` sub-regions are masked: default `["SIG","HIGH"]` (everything
  at or above `SIGSTART`); `["HIGH"]` masks only the bins at or above `SIGEND`, e.g. to drop a pile-up
  overflow bin while keeping it in the plots.
- A key with an unknown name is silently ignored. Writing `blindedFitNA` instead of `blindedFit` is
  therefore a way to switch blinding **off**: double check the key before unblinding anything.

The pass region of the signal region stays blinded until the validation regions and the MIP
reference are understood.

## Systematics

Only the **signal** carries named nuisances: the background is data driven, so its parameters are the
transfer-function coefficients and the free fail-region bin yields. Nuisances are named
`CMS_<CADI>_<name>`.

Expected signal systematics: luminosity (lnN), trigger efficiency (JetMET leg from the PFMETNoMu turn-on
in Muon data, ditau leg), the sizeX working point (pass/fail migration) and the probQ shape (both from
data/MC of the Z probes), pixel cluster
saturation modeling at high charge, track pT scale (pT = pT_true / Q), pileup, and signal theory
(scale, PDF; the photon-fusion PDF weights are not meaningful).

**Implement efficiency systematics of the sizeX tagger as a shift of the working point**, refilling pass
and fail, not as a flat scale of the pass histogram. Scaling pass alone breaks `pass + fail = total`,
and the transfer function absorbs the difference. Check that nominal and both variations give the same
`pass + fail`.

**When adding a systematic**, every list in the histogram-making script that loops over templates
(normalization, empty-bin sentinel, ...) must be extended, not only the booking. A template left
un-normalized builds a valid card and only shows up as a wildly pulled nuisance.

## Post-fit validation

- **Goodness of fit**: saturated test statistic with toys on condor; `gof_plot.{png,pdf}` and
  `gof_results.txt` (p-value) in the `*_area` directory. A p-value below about 0.05 means the background
  model does not describe the data.
- **Pulls and correlations**: `nuisance_pulls.{pdf,root}` and `plots_fit_{b,s}/correlation_matrix.*`
  from `fitDiagnosticsTest.root`. The bin-by-bin fail yields are left out of the matrix.
- **Impacts**: run per signal point in three modes: Asimov with signal injected (`-t -1 --expectSignal 1`),
  Asimov background-only, and (after unblinding) data. Seed the transfer-function parameters from the
  b-only fit (`rpf_params_*_fitb.txt`) in **both** `--doInitialFit` and `--doFits`, since each `-t -1`
  call regenerates its Asimov. Rank by impact on `r`: the fail-bin "pulls" are yields in event units,
  not sigma pulls.
- **Closure**: run the full chain with simulated background as `data_obs` (`..._BkgMC` configs) and on
  the validation regions before looking at the signal region.

## Limits

`combine -M AsymptoticLimits` per signal point gives the limit on `r`; multiplied by the signal cross
section of that point (DY plus photon fusion where available) it is the cross-section limit. Results are
shown as `sigma` vs `M` for each charge, and as excluded regions in the `(M, |Q|)` plane, separately for
spin-0 and spin-1/2.

## Versioning

Bump the input version (`histograms_for_2DAlphabet_v<N>`) when the histograms change and the binning
version (`Binningv<M>`) when the analysis binning changes. For each new version add a line below with
the ntuple production it comes from and what changed, and copy the previous configs and run scripts,
editing only `GLOBAL.path`, the `BINS`, and the hardcoded `configJSON`.

### Input versions
- (none yet)

### Binning versions
- (none yet)
