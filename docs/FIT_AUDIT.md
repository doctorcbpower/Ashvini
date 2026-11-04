# Fitting and mock recovery: pre-flight audit (interface steps 4 and 5)

Status: for agreement. Nothing is implemented. Builds on `INTERFACE_DESIGN.md`, `TREESAMPLE_AUDIT.md` and `RELATION_AUDIT.md`.

## 1. The scientific question

Which values of the free knobs make the model's relation (for example M_BH against M\*) agree with an observed relation, and how well does the data constrain each knob? A posterior over the knobs, not only a best fit.

Before any real data: if the model itself generates a relation at known parameters, does the fit recover them? (A fit that cannot do this says nothing about real data.)

## 2. Assumptions the fit makes

Physical
- P1. The model's relation (median, scatter, occupation of y given x) from `Model.relation`, on a fixed tree sample, is the prediction. All the limits of the earlier steps carry over: centrals only; z = 0; BH results trustworthy only for halo masses above about 1e11 Msun; the stellar mass is formed mass unless `x_scale` is changed.
- P2. The observed sample is selected on x (M\*), or its selection is declared and the range restricted (design decision 2). A BH sample selected on BH detectability is not modelled.
- P3. The data are detections of galaxies that have a black hole, with upper limits where given. The occupation fraction is not fitted in v1 (it is reported, and `f_occ` below 1 in a bin where data exist is a warning).

Statistical
- S1. Each data point is Gaussian in log10 y about the model median at its x, with variance
  `sigma_i^2 + s_model(x)^2 + s_int^2 + sigma_model(x)^2`:
  the measurement error, the scatter the model produces, an extra intrinsic scatter `s_int` that is a fitted parameter, and the bootstrap error of the model median (from the finite tree sample). Errors in x are carried into y by the local slope of the model relation.
- S2. The model median at x is interpolated linearly in log x between the centres of the unflagged bins. A data point outside the range of unflagged bins is excluded and counted, never extrapolated.
- S3. An upper limit contributes the Gaussian cumulative probability that the model lies below it.
- S4. Priors are uniform in the knob's sampling variable (log for knobs that span decades) within the registry bounds; `s_int` is uniform in log between 0.01 and 1 dex.
- S5. The tree sample is fixed during a fit, so the likelihood is a deterministic function of the knobs. The tree-sample noise is not removed, only measured (it enters through sigma_model).

Numerical
- N1. emcee ensemble sampler, 4 x (number of parameters) walkers at least, parallelised across processes.
- N2. Convergence is judged from the integrated autocorrelation time (chain length at least 50 times it), the acceptance fraction, and agreement between two independent runs started from different places.
- N3. The bootstrap is switched off inside the likelihood (the reference errors come from one run at a starting point and are held fixed), because it is the larger part of the cost.

## 3. What you have specified

Decisions 1 to 7 of the design note. The workflow in your first message: generate trees covering the stellar-mass range, run with a set of knobs, and use emcee against empirical data.

## 4. What I have inferred (please check)

- I1. The Gaussian-in-log-y likelihood above, including `s_int` as a free parameter.
- I2. Mock data of 40 points, uniform in log M\* over the unflagged range, with measurement errors of 0.3 dex and `s_int` of 0.2 dex. These are typical of the M_BH–M\* literature but are my choice.
- I3. Inclusion of the bootstrap model error in the variance.
- I4. Linear interpolation of the model median between bin centres.

## 5. Where the existing code is inconsistent with this, or cannot yet do it

- E1. **No sampler installed.** `emcee` (and `corner`, for the plots) are not in the environment. Installing them downloads packages; I need your permission, and they would be listed as optional dependencies (`requirements-dev.txt`, not `requirements.txt`).
- E2. **Cost.** A model run plus relation took about 2 s for 1000 trees on a laptop (with 100 bootstrap resamples, which the fit does not need). A fit of 4 knobs with 32 walkers and 3000 steps is about 100 000 evaluations, roughly 1 CPU-day at 1 s each, and about an hour on a 32-core node. It belongs on OzSTAR. The emulator of design decision 6 comes after a timed direct run.
- E3. **Bounds may not reach the data.** The fiducial model gives M_BH about 5 per cent of M\*, against 0.1 per cent or less observed. The accretion-efficiency bounds in the registry (3e-4 to 1e-2, with default 1e-3) were taken from the earlier exploration and may not span the factor of 50 needed, so a fit could pile up against a bound. A one-at-a-time scan of the BH knobs (below) answers this before any sampling.
- E4. **The M\* relation is not what a BH-only data set constrains.** M_BH against M\* data constrain how black holes grow for a given stellar mass, not how stars form. The star formation and wind knobs change which haloes host a given M\* and so move M_BH(M\*) only indirectly. I expect strong degeneracies, or posteriors equal to the prior, if they are freed together with the BH knobs.
- E5. **Seeding.** Below the BH-reliable mass the black hole mass depends on the tree resolution, so the usable M\* range starts at about 4e8 Msun, and the seed-mass knobs may be nearly unconstrained by data above it: the black hole has lost memory of its seed by z = 0 at 1e11 and 1e12 (`TREESAMPLE_AUDIT.md`, section 10).

## 6. Ambiguities that change the result materially

- A1. **Which knobs are free in the first fit.** (a) BH knobs only (accretion efficiency, feedback coupling, and the seed masses if they turn out to matter), with star formation and wind fixed at a set that already reproduces the stellar-to-halo relation. (b) Everything, constrained by data that include the stellar-to-halo relation (Moster/Behroozi, Nadler) as a second term. I recommend (a) first: it is cheap, identifiable, and cleanly separates what the data can tell us. (b) is the right end point but needs a joint likelihood and an audit of its own.
- A2. **Which stellar-mass fixed set for (a).** One of the accepted parameter sets from the exploration (`PARAMETER_SPACE_AUDIT.md`), or the preset values. The accepted sets were chosen to match the stellar-to-halo relation, the preset ones were not.
- A3. **Which data first.** The real data set (Reines & Volonteri 2015 for AGN in low-mass galaxies; Greene, Strader & Ho 2020; Kormendy & Ho 2013 for bulges) decides the mass range and the selection caveat. For the mock stage no data are needed.
- A4. **Same trees or different trees for the mock.** (a) Generate the mock from the sample the fit then uses: tests only the inference machinery ("inverse crime"). (b) Generate from an independent tree realisation: also tests that tree noise does not bias the answer. I recommend both, (a) first.
- A5. **Intrinsic scatter.** Free parameter (I1), or fixed. I recommend free, with its prior stated.

## 7. Minimum tests that separate a correct implementation from one that merely runs

1. **Likelihood by hand.** For a fake model relation with a closed form, the log-likelihood equals a direct evaluation of the Gaussian terms, including an upper limit, the excluded points count, and the interpolation between bins.
2. **Prior.** Samples respect the bounds; log and linear knobs map correctly.
3. **Coverage on a cheap synthetic model.** With a model whose relation is a closed-form function of its parameters, 50 or more mocks are fitted; the fraction of truths within the 68 per cent interval is within the binomial error of 0.68, and the rank statistics are uniform. This validates the fitting code without the cost of the tree runs.
4. **Recovery on the real model.** Three mocks (different seeds, true values chosen from the preset's neighbourhood, one near a bound): the truth lies inside the 95 per cent region for each identifiable parameter, and a parameter the data cannot constrain returns its prior, not a number.
5. **Convergence diagnostics are produced**, and the fit reports non-convergence when chains are cut short.
6. **Same input, same chain** (seeded).
7. **Bound pile-up.** A mock whose true value is outside the bounds produces a posterior piled against the bound and a visible warning.

## 8. Order of work

1. One-at-a-time scan of the BH knobs for their effect on M_BH at fixed M\*, with the tree sample and relation we have (cheap, local): decides E3, E4 and which knobs enter the first fit.
2. `Data`, the likelihood, and the synthetic-model coverage test (tests 1 to 3). No sampler needed for the likelihood tests.
3. `fit()` with emcee, diagnostics, and the real-model mock recovery on OzSTAR (tests 4 to 7).

## 9. What I need from you

- Permission to install `emcee` and `corner` into the project venv (E1).
- A1 and A2: BH knobs only with a fixed star-formation set, and which set.
- A3: the first real data set, if you have one in mind.
- I1 to I4 and A4, A5: accept or change.

## 10. Result of the one-at-a-time scan (3 Nov 2026, `scripts/bh_knob_scan.py`)

Fiducial preset, 1000 trees (the sample of `RELATION_AUDIT.md` section 10), the six unflagged bins at log M\* = 8.9 to 10.1, each knob at the ends of its registry range, others at the preset. Change in the median log10 M_BH at fixed M\* (the same in every bin to within 0.03 dex, so the slope is unchanged unless stated):

| knob | range | change in log M_BH at fixed M\* |
|---|---|---|
| bh_accretion_eff | 3e-4 to 1e-2 (preset 1e-3) | -0.51 to +0.98: M_BH proportional to the efficiency |
| sf_efficiency | 0.003 to 0.1 (preset 0.015) | +0.69 to -0.83: M_BH inversely proportional to it |
| bh_feedback_coupling (eta_agn) | 0.05 to 2 | 0.00 to 0.01 |
| bh_seed_mass_pop3, _direct_collapse, _threshold | full ranges | at most 0.03 |
| wind_loading | 1 to 25 | at most 0.05 |
| z_reion | 6 to 9 | at most 0.01 |

No knob changes the scatter in M_BH at fixed M\* (0.09 dex throughout, from assembly history alone) or the occupation fraction (1.0 in these bins).

Findings
1. **Only one combination matters here: bh_accretion_eff divided by sf_efficiency.** The two one-at-a-time exponents are +1 and -1, so I expect only the ratio to be constrained (not tested jointly). Black holes accrete in proportion to the gas, and at fixed M\* a lower star formation efficiency leaves more gas. A fit to M_BH against M\* alone fixes the normalisation of that ratio and nothing else; the accretion efficiency is meaningful only relative to the star formation efficiency it was fitted with.
2. **The slope is fixed by the model** at about 1.0 over this range, whatever these knobs do. If the data's slope differs significantly, no setting of these knobs can reproduce it. I have not compared with published slopes.
3. **AGN feedback (eta_agn), the seed masses, the wind and the reionisation redshift are unconstrained** by this relation above M\* of about 7e8. A fit would return their priors. This agrees with the earlier finding that the BH parameters had zero Sobol index for M\*.
4. **The registry range of bh_accretion_eff does not reach a ratio of 0.1 per cent.** M_BH/M\* at log M\* = 9.5 is 10^-0.23 at efficiency 1e-2, 10^-1.21 at 1e-3, 10^-2.19 at 1e-4, 10^-2.71 at 3e-5 and 10^-3.18 at 1e-5. At the bottom of the registry range (3e-4) it is 10^-1.7 (2 per cent). Reaching 0.1 per cent needs an efficiency near 2e-5, so the lower bound would have to move down by about 1.5 decades. (This assumes a target ratio of that order, which comes from my recollection of the literature and not from a data set we have chosen.)
5. **Other options (fixed values, not registry knobs):**
   - Eddington multiplier 0.3: slope 0.80, offset 0.24 dex lower, scatter 0.90 dex. Growth then hits the Eddington cap and the scatter explodes. 3.0: no change from the fiducial.
   - M-sigma feedback switched on alone: no change. Switched on with its growth cap: offset 1.2 dex lower (M_BH about 0.4 per cent of M\*), slope 1.07, scatter 0.13. The only option here that reaches a realistic ratio without moving the accretion efficiency far below its range.
   - Feedback delay of 0.3 Gyr: no change.
   - The slim-disc growth model at its default values: M_BH of order 1e11 Msun in 1e12 Msun haloes and M\* suppressed by about 3 dex (M\* of 2e7 against 1.5e10), so no bin of the fiducial M\* range is populated. It needs its own retuning together with star formation before it can be compared with anything. The slim-disc knobs in the registry should not be fitted until then.

Consequences for the plan
- The first BH-only fit is effectively one parameter (the log of the accretion-to-star-formation efficiency ratio) plus the intrinsic scatter that the data require. That is a one-dimensional problem and does not need MCMC; the grid of the likelihood is enough and better.
- A fit with more than one live BH parameter needs knobs that act on the slope or the scatter: the Eddington multiplier and the M-sigma switch with its growth cap. Neither is in the registry now.
- The MCMC machinery is still the right tool once the star formation, wind and BH parameters are fitted together against the stellar-to-halo relation as well (A1 option b), where the degeneracy above becomes a property of the joint posterior. That is a larger project.

## 11. Registry changes after the scan (3 Nov 2026)

Agreed: add the Eddington multiplier and the M-sigma cap as knobs, extend the accretion-efficiency range. Done in `ashvini/parameters.py` and `ashvini/model.py`:
- `bh_accretion_eff`: range now 1e-6 to 1e-2 (provisional; the exploration used 3e-4 to 1e-2).
- `bh_eddington_multiplier` (new): 0.1 to 10, log; pznk11_freefall growth model only.
- `bh_msigma_baryon_fraction` (new): 0.05 to 0.3, linear; needs the M-sigma growth cap, which the new preset `m_sigma_cap` switches on (`bh.sigma_feedback.enabled` and `growth_cap_enabled`). Its physical value is 0.16, so freeing it makes a theory prediction into a calibration; I kept it because it is the only continuous parameter of the cap, and noted that in its description.

What the new knobs do (fiducial trees, `scripts/bh_knob_scan.py`; log10 M_BH/M\* at log M\* = 9.5, slope, scatter in dex):

| setting | log M_BH/M\* | slope | scatter |
|---|---|---|---|
| fiducial, accretion efficiency 2e-5 | -2.89 | 0.98 | 0.08 |
| fiducial, Eddington multiplier 0.1 | -3.21 | 0.00 | 1.28 |
| 0.3 | -1.46 | 0.80 | 0.90 |
| 0.5 | -1.29 | 0.98 | 0.14 |
| 1 or 10 | -1.22, -1.19 | 1.01, 0.99 | 0.09 |
| m_sigma_cap, f_g = 0.05 / 0.16 / 0.30 | -2.93 / -2.42 / -2.15 | 1.07 | 0.13 |
| m_sigma_cap, accretion efficiency 1e-3 / 1e-5 | -2.42 / -3.20 | 1.07 / 0.99 | 0.13 / 0.09 |

Findings
- **The Eddington multiplier is a threshold knob.** At or above 1 it does nothing. Below about 0.5 the Eddington cap binds, the scatter grows, and at 0.1 black holes barely grow from their seeds (slope 0, scatter 1.3 dex), which is not a regime the data would accept. A flat posterior above about 0.5 and a sharp wall below is what to expect.
- **With the M-sigma cap, M_BH at the cap is proportional to f_g** (0.5 dex per factor 3.2) and independent of the accretion efficiency, as long as the cap is reached. At an efficiency of 1e-5 the cap is not reached and the result is accretion-limited, the same as without the cap. So the cap and the accretion efficiency are two regimes of one relation: M_BH = the smaller of the two. A fit that frees both will have a knee-shaped, non-Gaussian posterior along their trade-off. That is the first place in this problem where an MCMC (rather than a grid) is justified.
- The cap changes the slope slightly (1.07) and the scatter (0.13).

## 12. Implemented: Data, likelihood, prior, grid posterior (3 Nov 2026): `ashvini/fit.py`, `tests/test_fit.py`

Agreed: I1 to I4, A4 and A5 as proposed (Gaussian in log y with `s_int` free, bootstrap model error in the variance, linear interpolation between unflagged bin centres, intrinsic scatter free and uniform in log between 0.01 and 1 dex). The star-formation fixed set and the first data set are still open; neither is needed yet.

- `Data(log_x, log_y, sigma_y, upper_limit, sigma_x)`, `Curve` (median, scatter and median error at nodes; `Curve.from_relation` keeps unflagged bins and can take fixed reference errors), `log_likelihood(data, curve, s_int)`, `n_excluded`.
- `Problem`: free parameters with uniform priors in the sampling variable (log10 for log parameters, `s_int` last), `log_prob`, `sample_prior`. `predict(theta)` returns a Curve, so the closed-form model of the tests and the real `Model` use the same code.
- `grid_posterior`, exact up to the grid, for one or two parameters (section 10: the BH-only problem is one-dimensional). The sampler comes next and reuses `Problem.log_prob`.

Tests (9): the likelihood by hand (interpolation, upper limit, x errors carried by the slope); out-of-range points excluded and counted, never extrapolated; prior mapping and bounds; coverage on the closed-form model; centring on the truth for a large sample; all-outside data is an error; determinism. Coverage test: 150 mocks of 30 points with truths drawn from the prior, half of them upper limits from a detection threshold that does not depend on the truth; the truth lies in the central 68 per cent interval in the right fraction of mocks and the ranks are uniform (Kolmogorov-Smirnov). A 300-mock run gave 0.70 for both the offset and `s_int` (binomial error 0.027), rank-test p = 0.47 and 0.98.

A first version of the mock had a detection threshold built from the true offset, so the limits carried information about it. It gave coverage of 0.76 and p = 0.02 for the offset. That was a flaw in the test, not the likelihood; the threshold is now independent of the truth.

Not yet done: `fit()` with emcee and its convergence diagnostics; the binding to the real `Model` (reference errors held fixed, bootstrap off); recovery on the real model (tests 4 to 7, on OzSTAR).

## 13. Implemented: `fit()` with emcee and the Model binding (4 Nov 2026)

`ashvini/fit.py`:
- `fit(problem, n_steps, n_walkers, n_runs=2, seed, pool)`: `n_runs` independent emcee ensembles (default 2, at least 16 walkers, 4 per dimension), each started from different prior draws with a finite posterior; `pool` takes a `multiprocessing.Pool` or an MPI pool. Returns a `FitResult` (chain, tau, acceptance, R-hat, `summary()`, `physical()`, `corner()`).
- Convergence (N2): kept steps at least 50 times the integrated autocorrelation time; acceptance fraction between 0.1 and 0.9; split R-hat between the ensembles below 1.1. Any failure is in `FitResult.warnings` and issued as a warning.
- Bound pile-up (test 7): more than 10 per cent of the posterior within 2 per cent of the prior range at one end (a flat posterior puts 2 per cent there) gives a warning naming the parameter and the side.
- `problem_for_model(model, sample, data, y, x, bins, x_scale, reference_theta, n_boot)`: parameters are the free knobs; the bootstrap errors of the model medians come from one run at the reference point and are held fixed (interpolated at each bin), so each likelihood evaluation is one model run and no bootstrap (N3).

Tests (tests 5 to 7 and the sampler version of test 3, 15 in the file): recovery of a known offset, with the sampler's 16th to 84th percentiles agreeing with the exact grid posterior; a short chain is reported, not trusted; same seed gives the same chain; a truth outside the bounds piles up with a warning; data the model cannot reach is an error; and a stand-in model confirming one bootstrap in total and recovery of a known knob. A 1500-step chain was flagged as too short (tau 29 to 38 steps); 3000 steps were needed.

Still to do: recovery on the real `Model` (test 4: three mocks, one near a bound; needs OzSTAR and a choice of fixed star-formation set, for which the preset values are the default) and the first real data set.

## 14. Recovery on the real model: design (4 Nov 2026) `scripts/fit_recovery.py`, `hpc/submit_fit_recovery.sbatch`

Question: if the model itself generates the data, does `fit` return the parameters that made them, and does a parameter the relation cannot constrain come back as its prior?

Setup (my choices; change any):
- Model: `fiducial` preset, so star formation, winds and reionisation sit at the preset values (A2 default). Free: `bh_accretion_eff` (1e-6 to 1e-2, log) and `bh_feedback_coupling` (eta_agn, 0.05 to 2, log), plus `s_int`. eta_agn is there because section 10 found it changes nothing: it should return its prior.
- Tree sample: `TreeSample.for_stellar_masses((8.0, 10.5), n_trees=1000, seed=1)`.
- Mock: 40 points uniform in log M\* over the unflagged bins, measurement error 0.3 dex, `s_int` 0.2 dex (I2), drawn from the model relation at the truth (median, plus model scatter, `s_int` and the measurement error in quadrature).
- Truths (accretion efficiency, eta_agn): case 0: (2e-4, 0.5); case 1: (1e-3, 0.1), the preset; case 2: (3e-6, 1.0), a factor of three above the lower bound of the accretion efficiency. Cases 3 to 5 repeat these with the mock drawn from an independent tree realisation (seed 2) while the fit keeps seed 1 (A4b).
- Sampler: two ensembles of 16 walkers, 4000 steps (1000 discarded), 32 forked workers. Cost: the local smoke run (600 trees, 4 workers on a laptop, 60 steps) took 2.0 s per evaluation per worker, so 128 000 evaluations is about 2.2 hours per case on 32 cores at 600 trees, perhaps 3.5 hours at 1000 (not measured on the cluster; the job asks for 8 hours).

What would count as success: for each case, the accretion efficiency truth in the 95 per cent interval (and in the 68 per cent interval in roughly two thirds of the six cases); `s_int` the same; eta_agn flat across its prior range (Kolmogorov-Smirnov test against uniform, and about half the posterior in the central half of the range); no convergence warnings. Case 2 may show a mild pile-up if the lower bound is felt; a case where the truth sits outside the 95 per cent interval for the accretion efficiency in a same-trees mock would be a defect in the machinery, and in an independent-trees mock a sign that tree noise biases the answer. With six cases, one or two misses at 68 per cent are expected and mean nothing.

Not covered: the fit uses the same bins and the same reference errors as the mock, so it does not test sensitivity to the flagged-bin thresholds.
