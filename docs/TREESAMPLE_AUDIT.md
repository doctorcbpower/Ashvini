# TreeSample: pre-flight audit (interface step 2)

Status: for agreement. Nothing is implemented. Facts about the stored trees were checked on 3 Nov 2026; items marked "unverified" were not.

## 1. The scientific question

Given a model, what is the relation between two galaxy properties (for example M_BH and M*) at redshift z = 0, for a *population* of galaxies that can be compared with an observed sample selected on stellar mass?

A merger tree is a halo history rooted at one halo mass. Observed samples are selected on M*. So the sample of trees has to be turned into a population in M*, with each tree carrying the right abundance. Getting this wrong biases the relation even if the physics is right.

## 2. Assumptions the TreeSample would make

Physical
- P1. Each tree is an isolated central halo's main-progenitor history. No environment, no subhaloes, no satellites (stripping, quenching).
- P2. Abundance of haloes follows a halo mass function n(M_h) for the same cosmology and mass definition as the trees.
- P3. The stellar-mass selection is applied to the model's own M*, after the run.

Statistical
- S1. The set of trees is a stratified sample of the halo population: strata in M_h(z=0), weights = n(M_h) x stratum width / number of trees in the stratum.
- S2. Weighted medians, percentiles and scatter are computed with those weights.
- S3. The tree set is fixed during a fit (deterministic likelihood); the selection and the weights are recomputed at each parameter value, because M* changes.

Numerical
- N1. Merger trees from foraois, Zhang-Hui algorithm, z_max = 30, dz = 0.005, kmax override 3000, M_res = 1e-4 x M0, as in the existing stored sets.
- N2. Trees cached on disk under a hash of all of these settings and the seed.

## 3. What you have specified

Decisions 1, 3, 4, 7 of the design note: select on the model's M* afterwards and weight by the halo mass function; zero-BH haloes reported separately; median and scatter; z = 0 only.

## 4. What I have inferred (please check)

- I1. Strata evenly spaced in log M_h, with the same number of trees each.
- I2. The halo-mass range is set from an assumed M*(M_h) plus a generous margin, only to decide how far to extend; the answer should not depend on it (tested, see 7).
- I3. The halo mass function is Sheth-Tormen for the Planck 2018 cosmology in the foraois config. Unverified: I have not yet checked whether foraois or pymctrees already provides a mass function, or which mass definition its trees use. The weights depend on this.
- I4. Centrals only (P1) is acceptable for v1.

## 5. Where the existing code is inconsistent with this

- E1. **Stored trees are not a population.** The 9 stored sets have all 100 trees at exactly the same M0 (checked: min = max for 1e8, 1e10 and 1e12), spaced 0.5 dex apart, with equal weight. Any relation computed *conditional on M\** from them over-weights rare haloes. The conditional-on-halo-mass numbers in earlier audits are not affected.
- E2. **Eddington bias.** With a steep mass function and scatter in M* at fixed M_h, the haloes at fixed M* are mostly lower-mass haloes scattered up. Equal-weight strata hide this. Weighted selection reproduces it only if the strata are fine enough that the scatter (0.3 to 0.5 dex at low mass from assembly alone) is sampled smoothly, so 0.5 dex spacing is marginal.
- E3. **Resolution scales with M0** (M_res = 1e-4 x M0). A 1e12 tree resolves haloes only above 1e8 Msun; a 1e8 tree resolves down to 1e4. The Pop III (minihalo, 1e6 Msun) and direct-collapse (1e7 Msun) seeding thresholds therefore cannot trigger early in massive trees, because the progenitor is not resolved when it would seed. The stored 1e12 trees do show BH masses of 100 Msun for some haloes and 1e5 for others, so some seeding does occur, but I have not measured how the seeding channel and seed time depend on M_res. This affects M_BH directly and is the largest risk to a BH relation. Unverified how far it goes.
- E4. From my notes, the foraois compiled kernel clamps the variance at small scales (audit in ~/MyCodes/foraois_sigma_audit, fix-before-release awaiting your decision). It affects first-resolved quantities, which is what seeding uses. I have not rechecked whether the version that built the stored trees has the fix.
- E5. Tree generation takes the final mass as `np.full(n, M0)` in the adapter, although the underlying foraois call takes an array. I have not checked that a per-tree array works end to end.

## 6. Ambiguities that change the result materially

- A1. **Strata or continuous.** (a) Keep strata (grid of M0, N each, weights), with a finer grid (e.g. 0.1 to 0.15 dex, about 30 to 40 strata over 4 dex). (b) Draw every tree's M0 from a continuous distribution (the array form). (b) removes the lumps and the choice of grid. It needs E5 checked and gives weights that are all equal if M0 is drawn from the mass function itself, which is the simplest statistics. Cost is the same per tree; (b) needs more trees in the sparse high-mass tail, so draw from a flattened distribution and weight back.
- A2. **Resolution policy.** (a) Keep M_res = 1e-4 x M0. (b) Fix the absolute M_res (say 1e5 Msun) for every tree, so the seeding thresholds are resolved everywhere. (b) is more physical for BH seeding and costs much more for massive trees. A measurement of how M_BH(1e12) changes between the two should decide it.
- A3. **Centrals only or include satellites.** Observed BH samples (for example AGN in dwarfs) include satellites, and satellites are stripped. Centrals only will be wrong at some level in the dwarf regime. For v1 I would state it, not model it.
- A4. **Halo mass function source and definition** (I3).

## 7. Minimum tests that separate a correct implementation from one that merely runs

All with a fake galaxy model where M* is an exact known function of M_h with a known log-normal scatter and the BH mass is an exact function of M*, so the right answer is analytic:

1. **Abundance.** The weighted number of trees in any M_h range equals the integral of the mass function over it.
2. **Recovery of Eddington bias.** The weighted mean log M_h at fixed M* matches the analytic conditional mean for the chosen mass function and scatter; equal weights fail this test.
3. **Grid independence.** Halving the stratum spacing (or doubling the number of continuous draws) changes the weighted median and scatter of the relation by less than a stated tolerance.
4. **Range independence.** Widening the halo-mass range by 0.5 dex does not change the relation inside the selected M* window. A completeness check must warn when the selected window's edge weights are not negligible.
5. **Determinism and cache.** The same inputs give identical samples; changing any tree setting changes the hash and rebuilds.
6. **Zero-BH accounting.** Weighted occupation fraction plus weighted occupied fraction equals one.
7. **Real-model check (not synthetic).** For the fiducial Model, the weighted relation does not depend on the stratum grid, and the dependence of M_BH on M_res is measured (A2) and reported before anything is fitted.

## 8. What I need from you

- A1: strata with a finer grid, or continuous M0?
- A2: fixed fractional resolution (as now), or fixed absolute resolution for BH work? If unsure, I will measure first (test 7) and report.
- I3/A4: whether you have a preferred halo mass function and mass definition, or I use the one consistent with foraois after I check it.
- E4: whether the sigma-clamp fix is in the version that should build the new trees.

## 9. Update, 3 Nov 2026

- **A1 decided: continuous halo masses.**
- **E5 checked.** foraois `build_forest_numpy` accepts a per-tree array of final masses. A 30-tree test (final masses log-uniform 1e9 to 1e11 Msun, M_res = 1e6, dz = 0.05, z_max = 10) ran, the histories decrease monotonically back in time, and all trees are resolved at z = 10. The Ashvini adapter (`build_forest_live`) still takes a single mass, so it needs a variant that takes an array.
- **New coupling.** One call has one `M_res`. With continuous masses a fractional resolution (M_res = 1e-4 x M0) would need the trees to be built in batches by mass range. A fixed absolute M_res needs no batching. This pushes toward the fixed-resolution option of A2.
- **E4.** The sigma-clamp fix is in foraois (CHANGELOG v0.1.2, "P1": low-mass sigma(M)/alpha(M) table clamp). New trees built with the installed foraois get it. I have not established whether the stored trees were built before or after it.
- **I3.** foraois has a `MassFunctions` class (`mass_function_utils.py`). I have not yet read which mass functions it provides or its mass definition.
- **A2: fixed absolute resolution chosen, then found not feasible; reopened.** The tree builders register at most one split per step, so the expected number of splits per step must stay below about 0.1 (`_warn_if_timestep_noncompliant`). Checked with the foraois diagnostic (worst case over redshift, evaluated at M0, so conservative), z_max = 30:

  | M0 (Msun) | M_res (Msun) | M_res/M0 | expected splits per step, dz = 0.005 | dz = 0.002 |
  |---|---|---|---|---|
  | 1e8 | 1e4 | 1e-4 | 0.47 | 0.19 |
  | 1e11 | 1e5 | 1e-6 | 40 | 16 |
  | 1e12 | 1e6 | 1e-6 | 50 | 20 |
  | 1e12 | 1e5 | 1e-7 | 360 | 140 |
  | 1e12 | 1e4 | 1e-8 | 2500 | 990 |

  Even the stored sets (ratio 1e-4) are at 0.47, above the 0.1 guideline, so they were already the finest affordable. A fixed absolute resolution of 1e4 to 1e6 Msun for the massive trees would need dz of order 1e-4 or smaller, which is millions of steps per tree. **The resolution ratio M_res/M0 cannot go much below 1e-4, so the resolution has to be fractional.** The earlier statement that fixed resolution "costs more" understated this: it is not available with this builder.

## 10. Measured: black hole results against tree resolution (3 Nov 2026)

`scripts/resolution_bh_scan.py`: 200 Zhang-Hui trees per case, fiducial Model (the run_params.yaml defaults, not a tuned set), M_res/M0 = 1e-3 and 1e-4 (the finest allowed at dz = 0.005), and a second seed at 1e-4 as the noise estimate. Medians; "occupied" means haloes with a black hole.

| M0 (Msun) | M_res/M0 | with BH | channel that fired | z seeded | log M_BH(z=0), occupied | log M* |
|---|---|---|---|---|---|---|
| 1e9 | 1e-3 | 0.09 | Pop III | 16.0 | 2.79 | 4.16 |
| 1e9 | 1e-4 | 0.01 / 0.00 (two seeds) | Pop III | 15.5 | 2.84 | 4.34 / 4.35 |
| 1e10 | 1e-3 | 0.68 | direct collapse (1e5 Msun) | 11.7 | **5.08** | 5.27 |
| 1e10 | 1e-4 | 0.57 / 0.59 | Pop III (100 Msun) | 17.0 / 16.7 | **3.36 / 3.50** | 5.35 / 5.48 |
| 1e11 | 1e-3 | 1.00 | threshold 0.54, direct collapse 0.47 | 4.0 | 7.23 | 8.48 |
| 1e11 | 1e-4 | 1.00 | direct collapse 0.70, Pop III 0.28 | 13.9 / 14.1 | 7.29 / 7.30 | 8.51 / 8.53 |
| 1e12 | 1e-3 | 1.00 | threshold 0.93 | 5.2 | 8.59 | 9.84 |
| 1e12 | 1e-4 | 1.00 | direct collapse 0.85 | 12.1 / 12.1 | 8.67 / 8.67 | 9.88 / 9.88 |

Findings
1. **M\* is converged** (differences below 0.1 dex, comparable to the seed noise).
2. **At 1e10 Msun, M_BH changes by about 1.7 dex with resolution.** At 1e-3 the tree first resolves the halo at z of about 11, below the Pop III redshift limit (15), so the direct-collapse channel seeds 1e5 Msun black holes. At 1e-4 the halo is resolved at z of about 15.5, the Pop III channel (checked first) fires, and a 100 Msun seed that has grown only to a few thousand Msun by z = 0 blocks the direct-collapse channel. Which channel fires depends on when the tree starts, not only on the physics.
3. **At 1e9 Msun the occupation fraction is not converged** (0.09 to 0.01 to 0.00, noise-limited at the low end).
4. **At 1e11 and 1e12 Msun, M_BH at z = 0 is robust**: 0.06 and 0.08 dex between the two resolutions, and 0.01 and 0.00 dex between the two seeds. The seeding channel and redshift change a great deal (z of 4 against 14 at 1e11), but the black hole has lost memory of its seed by z = 0.
5. M_BH/M* at 1e11 Msun is about 0.06 in the fiducial model, far above observed relations. The fiducial BH parameters are untuned, so this is what a fit would have to correct. I have not examined it further.

Limits of this measurement: 200 trees per case, one decade of resolution (1e-4 is the finest the builder allows), the default seeding thresholds and fiducial parameters. Whether M_BH at 1e11 and 1e12 stays converged for other BH parameters is not shown.

Consequence for TreeSample
- Trees are built with fractional resolution M_res/M0 = 1e-4, in batches by mass range.
- M_BH relations are trustworthy for M_h(z=0) of about 1e11 Msun and above (M* of about 3e8 and above in the fiducial model). Below that, the black hole mass and occupation depend on the tree resolution through the priority order of the seeding channels. `TreeSample` should warn when a BH relation is requested below that mass.
- Making the dwarf regime resolution-independent needs a change to the seeding prescription itself, which is a modelling decision and needs its own audit.

## 11. Implemented (3 Nov 2026): `ashvini/tree_sample.py`

Decisions applied: continuous halo masses (stratified, even in log10 M), fixed resolution ratio M_res/M0 = 1e-4 built in batches of 0.1 dex (never finer than 1e-4), Sheth-Tormen weights from foraois, z = 0 only, BH relations flagged below 1e11 Msun (`TreeSample.bh_reliable`; the warning itself belongs to the relation step, which knows which haloes were selected).

Tests (`tests/test_tree_sample.py`): weights sum to the analytic number density; with a known steep mass function and log-normal scatter the weighted selection at fixed M* recovers the analytic Eddington shift (0.33 dex for sigma = 0.4) and equal weights do not; the selected relation is unchanged by widening the halo range and is changed by cutting it inside the selection; cache hit, cache miss on a changed seed, trees end at their own masses, `run` works with a `Model`.

Open
- **Mass function checked (3 Nov 2026).** The foraois Sheth-Tormen function agrees with an independent calculation (CAMB P(k) with the same parameters, my own top-hat sigma(M), the standard Sheth-Tormen form with delta_c = 1.686) to within 0.2 per cent from 2e8 to 3e14 Msun, in shape and normalisation (`test_foraois_mass_function_matches_an_independent_sheth_tormen`). sigma8 comes out 0.849. The mass fraction in haloes between 1e8 and 10^14.9 Msun/h is 0.57 in both; my earlier expectation of 0.7 to 0.8 was wrong. What this does not test is whether Sheth-Tormen is the right mass function: against N-body mass functions it is known to be off by order 10 to 30 per cent at some masses, which would tilt the weights by a similar amount across 1e9 to 1e12 Msun. That is an accuracy limit of the weights, not an implementation error.
- **Mass definition.** foraois masses are linear-theory top-hat masses, not M200c or Mvir; they differ from the masses of published abundance-matching relations by roughly 0.1 to 0.2 dex.
- **Halo range.** `for_stellar_masses` sets the range from the Moster+13 relation with a 1 dex margin. The fiducial model forms fewer stars than Moster at 1e12 Msun, so the upper edge may be too low for M* of about 3e10; the relation step must check that edge weights are negligible and say so if not.
- The lowest halo mass is 10^8.2 Msun (the mass function's table starts at 1e8 Msun/h).
