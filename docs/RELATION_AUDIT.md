# relation(): pre-flight audit (interface step 3)

Status: for agreement. Nothing is implemented. Builds on `docs/INTERFACE_DESIGN.md` (decisions 1 to 7) and `docs/TREESAMPLE_AUDIT.md`.

## 1. The scientific question

For one set of parameters, what relation does the model predict between two galaxy properties at z = 0 (for example M_BH against M\*), for the population of central galaxies whose abundance follows the halo mass function, in a form that can be compared with an observed sample selected on M\*?

That means: the distribution of y *given* x, in bins of x, with its scatter and its occupation fraction (the fraction of galaxies that have y > 0), and an estimate of how much of the result is noise from the finite number of trees.

## 2. Assumptions the relation step makes

Physical
- P1. One tree is one central galaxy, with its z = 0 properties read from the last time step (z = 0 exactly).
- P2. The population is the weighted tree sample (`TreeSample.weight`), so no environment, no satellites, no merging companions.
- P3. y given x is what an observation selected on x measures. An observation selected on y (a BH sample selected by detectability) measures something different; that is handled by restricting the range, not by modelling it (design decision 2).

Statistical
- S1. Bins are fixed intervals in log10 x. Within each bin: weighted median, 16th and 84th percentiles of log10 y, from `weighted_quantile`.
- S2. Scatter is half the 16th-to-84th range of log10 y, in dex.
- S3. Occupation fraction is the weighted fraction of trees with y > 0 in the bin. The median, percentiles and scatter use the occupied trees only; empty ones are never dropped silently (design decision 3).
- S4. Uncertainty comes from a bootstrap over trees (below), giving a standard error on the median, scatter and occupation fraction in each bin.
- S5. No scatter is added beyond what the model produces (assembly history and stochastic inflow if present). Adding intrinsic scatter belongs to the fit.

Numerical
- N1. A bin is reported only if its Kish effective number of trees is at least a minimum; otherwise its values are NaN and the bin is marked.
- N2. Each bin carries flags rather than being removed: *low statistics*, *range-limited*, *BH-unreliable* (below).

## 3. What you have specified

Design decisions 1 (select on the model's M\* afterwards, weighted), 2 (restricted range, caveat declared), 3 (zero-BH separately), 4 (median and scatter), 7 (z = 0). The restriction of BH relations to halo masses of about 1e11 Msun and above (the resolution result, `TREESAMPLE_AUDIT.md` section 10).

## 4. What I have inferred (please check)

- I1. Bin width 0.25 dex in log10 x.
- I2. Minimum effective number of trees per bin: 20.
- I3. Range-limited flag: a bin is flagged if more than 2 per cent of its weight comes from trees within 0.3 dex of either end of the sampled halo-mass range (then the selected population plausibly extends beyond what was sampled, as it does when scatter moves haloes across the edge).
- I4. BH-unreliable flag: a bin is flagged if more than 10 per cent of its weight comes from trees with final halo mass below 1e11 Msun (`TreeSample.bh_reliable`). Requesting a BH quantity prints a warning that lists the flagged bins.
- I5. Bootstrap resampling in blocks of adjacent trees (sorted by M0, blocks of about 10), with each tree keeping its own weight. The sample is stratified in log M0, so resampling single trees would overstate the noise from the halo-mass dimension.
- I6. 200 bootstrap resamples by default, reproducible from a seed. They reuse the one model run: no extra model evaluations.
- I7. x and y are named outputs of the run (`stars_mass`, `bh_mass`, `gas_mass`, `sfr`, ...) or a function of the run's results; both are taken at the last step.
- I8. Flagged bins are returned and marked. A fit uses only unflagged bins by default.

## 5. Where the existing code is inconsistent with this

- E1. **Stellar mass counts all mass ever formed.** Ashvini's `stars_mass` has no stellar mass loss: gas returned by dying stars is not removed from the stellar mass. Observed stellar masses are present-day masses after loss. With instantaneous recycling at R = 0.44 (the value used by GRUMPY) the model's M\* would be about 0.25 dex too high against the same galaxy observed. Separately, observed stellar masses assume a Chabrier or Kroupa IMF, while the Config's IMF is Salpeter, which would add up to another 0.2 to 0.25 dex. Both shift M\* in M_BH–M\* directly, and by the slope of the relation shift the apparent BH normalisation. The SF efficiency in a fit can partly absorb it, but then it is wrong for the right reason only by accident. I did not correct for either in the earlier SHMR comparisons either (noted in `WIND_SCALING_AUDIT.md`); that comparison is affected the same way.
- E2. The earlier scripts computed relations conditional on halo mass (percentiles in bins of M_halo), which does not need weights. The new step is conditional on M\*, which does. Nothing earlier is invalid, but none of it is a prediction for an M\*-selected sample.
- E3. `space_common` reduces runs to medians per tree-mass stratum. It cannot be reused for this.

## 6. Ambiguities that change the result materially

- A1. **Which stellar mass.** (a) Use `stars_mass` as the model produces it and state that it is formed mass. (b) Multiply by a stellar-mass retention factor (default 1.0, user-settable) and by an IMF offset (default 1.0), both recorded in the relation's metadata. (c) Apply R = 0.44 and a Chabrier conversion by default. I recommend (b) with defaults of 1.0 and a visible note: (c) bakes an unexamined recycling model into the comparison, and the right default depends on whether the SF efficiency and the observations are on the same footing, which is your call.
- A2. **Scatter definition.** Half the 16-to-84 range as above, or the standard deviation. The first is robust to the long tails seen at low mass; I recommend it.
- A3. **Flagged bins in a fit.** Exclude by default (recommended) or include with their larger errors.
- A4. **x, y from satellites' perspective.** Not applicable in v1 (centrals only), recorded here so it is not forgotten.

## 7. Minimum tests that separate a correct implementation from one that merely runs

All on a fake galaxy model whose answer is analytic. A function `relation_from_arrays(x, y, weight, m0, ...)` does the work, so no trees are needed.

1. **Known power law.** y = g x^b with log-normal scatter sigma_y independent of everything else: the recovered median slope, normalisation and scatter match within the bootstrap error.
2. **Eddington bias through the relation.** M\* = f M_h with scatter and y depending on M_h: the recovered median of y given M\* equals the analytic value including the shift 0.9 ln(10) sigma^2 (as in the TreeSample test), and equal weights miss it.
3. **Occupation.** A known probability of having y > 0 as a function of x is recovered; the median and scatter of the occupied are unchanged by the empty ones.
4. **Low statistics.** A bin with fewer than the minimum effective trees is NaN and flagged, not reported.
5. **Range flag.** Fires when the halo range is cut inside the selection (the truncated-range case), and does not fire for an adequate range.
6. **BH flag.** Fires when more than 10 per cent of the weight comes from trees below 1e11 Msun.
7. **Bootstrap.** Standard errors fall as 1 over the square root of the number of trees (four times the trees, half the error, within noise); reproducible from the seed; unchanged by multiplying all weights by a constant.
8. **Scaling.** Multiplying x by a factor shifts the relation by exactly that factor in log x (the check for the stellar-mass factor of A1).
9. **Integration.** On a small real tree set with the fiducial `Model`, `sample.run` then `relation` returns bins of the expected shape, its occupation fraction agrees with a direct count, and the flagged bins are the ones expected from the halo mass range.

## 8. Interface (proposed)

```python
rel = model.relation(sample, theta, y="bh_mass", x="stars_mass", bins=np.arange(8.0, 11.5, 0.25),
                     x_scale=1.0, n_boot=200, seed=0)
rel.table()        # per bin: x centre, effective n, occupation fraction, median, 16th, 84th, scatter, errors, flags
rel.plot()         # median and 16-84 band, occupation fraction below, flagged bins marked
rel.metadata       # the settings, x_scale, bin edges, thresholds: so a result says what it was
```

`x_scale` is the retention and IMF factor of A1 (a single multiplier on x), recorded in the metadata.

## 9. What I need from you

- A1: the default for the stellar-mass factor: 1.0 and a visible note (my recommendation), or R = 0.44 and a Chabrier offset applied?
- I1 to I5: the bin width, the thresholds (20 trees, 2 per cent at the range edge, 10 per cent below 1e11 Msun) and the block bootstrap: accept, or give different values.
- A3: flagged bins excluded from fits by default?

## 10. Implemented (3 Nov 2026): `ashvini/relation.py`, `Model.relation`

Agreed: `x_scale` default 1.0 with the note in the output and metadata (A1). The other inferred values (I1 to I8) were used as proposed and are arguments (`bins`, `min_ess`, `edge_dex`, `edge_fraction`, `bh_fraction`, `block`, `n_boot`); flagged bins are returned and marked, and `Relation.usable` gives the unflagged ones for a fit (A3).

Tests (`tests/test_relation.py`, 10): known power law with its scatter; the Eddington shift through the relation (equal weights miss it by more than 0.3 dex); occupation fraction recovered, and empty trees leave the occupied statistics unchanged; low-statistics bins are NaN and flagged; range flag fires on a truncated range and not on an adequate one; BH flag and warning; bootstrap error scales as 1/sqrt(n), is reproducible, and is unchanged by rescaling the weights; `x_scale` shifts the bins exactly; the stars note; an integration test on a real small tree set.

### First real case: fiducial model, 1000 trees, M_BH against M\*

`TreeSample.for_stellar_masses((8.0, 10.5), n_trees=1000, seed=1)` (halo range 10^9.57 to 10^12.95 Msun, effective number of trees 290), `Model("fiducial").relation(sample, y="bh_mass", x="stars_mass")`. Run and relation take 2 s once the trees exist (the trees took 46 s to build). Medians in the unflagged bins (`log10 x`, `log10 M_BH`, scatter in dex):

| log M\* | log M_BH | scatter |
|---|---|---|
| 8.87 | 7.66 | 0.09 |
| 9.12 | 7.93 | 0.08 |
| 9.36 | 8.13 | 0.09 |
| 9.63 | 8.42 | 0.09 |
| 9.87 | 8.65 | 0.09 |
| 10.11 | 8.92 | 0.09 |

(the bin at log M\* = 10.36 is range-limited; those below 8.7 are BH-unreliable). Every galaxy in these bins has a black hole. Observations from this run:
- The model has M_BH of about 5 to 6 per cent of M\* across the range, with a slope of about 1.0. The fiducial BH parameters are untuned, so this is the offset a fit would have to remove; observed ratios at these masses are of order 0.1 per cent or less. I have not compared against a specific data set.
- The scatter in M_BH at fixed M\* is only 0.08 to 0.12 dex, from assembly history alone.
- The flags behaved as predicted: bins below M\* of about 4e8 (log 8.6) are BH-unreliable, matching the resolution result, and the top bin is range-limited, because the fiducial model forms fewer stars than the Moster relation the halo range was built from. `for_stellar_masses` with a larger `margin_dex` (or the range from a build with an upper halo mass higher than 10^12.95) fixes it.
- The tree builder's timestep check warns for this sample (expected splits per step 1.4 at its most massive trees, 9e12 Msun, against 0.47 for the stored 1e12 sets and the 0.1 guideline). The most massive bins are therefore less well resolved in assembly time than the earlier sets; I expect the bins shown to be fed mostly by haloes below about 2e12 Msun, where the check is closer to the stored sets, but I have not verified that.
