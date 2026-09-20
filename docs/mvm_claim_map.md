# MVM claim map

> **The central result is a conditional statement about the initial BH mass required to reach a specified
> M_BH/M_star,tot target within this baryon-cycle model, not a prediction of how black-hole seeds are formed.**

A map of what the frozen model supports, what it only qualifies, and what must not be claimed. Every row points to an entry of
`docs/mvm_numerical_crib_sheet.md` (C-numbers) and, where relevant, a figure. Numbers are quoted from the crib sheet; none is new here
except where marked. Nothing in this file is manuscript text.

## 0. Provenance chain

| Step | Repository | Commit | Content |
|---|---|---|---|
| 1 | foraois | `ed28326` | collapse barrier normalised to the P(k) redshift, with tests: the dependency used to generate the production trees |
| 2 | foraois | `9d5fdae` | PCH08 diagnostic note, reproduction script, ROADMAP entry (no source change) |
| 3 | Ashvini | `695b114` | frozen MVM, pre-MVM reference, production ensemble, figures, original diagnostic scripts (its logs were omitted by mistake; see 4) |
| 4 | Ashvini | `8ffb61e` | C17 checks, PCH08 diagnosis, all 27 diagnostic logs, crib-sheet update; corrects 3 (which did not contain the logs) |
| 5 | Ashvini | (this commit) | claim map, accessibility-reach diagnostic and log, crib-sheet C9b |

* Between foraois `1ba7073` (the last commit before this work) and `9d5fdae`, the only source change is `src/foraois/cosmo_utils.py` (the barrier fix), so the
  foraois code at `ed28326` is the state that generated the production trees.
* Frozen model: `ashvini/reservoir_stock.py`, sha256 `31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63`, recorded in
  `scripts/paper_figures/output/mvm_production_results.json` (2026-09-20T15:44:29). Every C17 script asserts this hash before running.
* The numba tree sampler is not seed-reproducible: the stored JSON and the logs are the record of the ensembles.

## 1. Headline results

| ID | Claim | Evidence | Numbers to quote | Qualification |
|---|---|---|---|---|
| H1 | The critical seed boundary rises monotonically with M_halo(z=5); it is not a power law and has a kink at M_hot | C1, C3; Fig 1 | 9.45e7 [7.61e7, 1.06e8] (3e10), 1.39e9 [1.03e9, 1.61e9] (3e11), 4.85e9 [4.54e9, 5.17e9] (3e13) Msun; local slope 1.17 -> 0.09 | conditional on M_hot (R3); 3e13 dz and M_res untested (Q1) |
| H2 | The boundary is set predominantly by M_star,tot: M_seed,crit = f_BH M_star,tot / G_BH with G_BH close to 1 | C4; Fig 1 lower panel | 1/G_BH = 0.973 [0.971, 0.975] -> 0.929 [0.922, 0.935]; M_star,tot/M_halo = 6.5e-3, 9.7e-3, 3.5e-4 | the decomposition is an identity; its physical reading is model-conditional; G_BH - 1 shrinks with dt (Q2), so quote "close to unity", not "grows by 3 to 8 per cent" |
| H3 | At the fiducial the critical seed is supply-limited: the Eddington cap does not bind | C5 | 0% of growth steps and 0% of mass gained Eddington-limited; 0.028, 0.044, 0.074 of about 20.8 e-folds used; peak Mdot_acc/Mdot_Edd 0.04 (3e10), 0.07 (3e11) | fiducial only; a statement about the critical boundary, not "Eddington is irrelevant" (light seeds: S1) |
| H4 | The inverse question: at the fiducial a fixed seed ends at M_BH/M_star,tot of about M_seed/M_star,tot | C6; Fig 2 | at 3e10: 4.0e-7 (1e2), 4.1e-6 (1e3), 8.2e-4 (2e5), 4.2e-2 (1e7); target 0.5 | 1e7 is a benchmark seed mass, not a channel; ratios uncertain at the few per cent of C13 |
| H5 | The critical seed is heavier than its host's baryons when the host is first resolved, and stays so until z of about 10 | C10; Fig 3 | seed/(f_b M_halo) at first resolution 5.1e3 [3.3e3, 1.7e4], 4.9e4 [4.8e3, 1.3e5], 3.8e3 [8.5e2, 2.6e4]; median M_BH/(f_b M_halo) crosses 1 at z = 11.0, 9.7, 12.2 and is 0.021, 0.031, 0.0011 at z = 5 | a diagnostic; the seed is placed at z_seed regardless of its host; interpretation is model-conditional |

## 2. Secondary results (nuclear accessibility and the seed channels)

| ID | Claim | Evidence | Numbers to quote | Qualification |
|---|---|---|---|---|
| S1 | Light seeds can enter an Eddington-limited transition when nuclear accessibility is increased | C7, C8, C9; Fig 2 dashed | 1e3 seed: peak Mdot_acc/Mdot_Edd 0.09 at the fiducial (3e10), crossing 1 between R_nuc 125 and 150 pc (sigma_j = 0.5) or sigma_j 0.7 and 0.8 (R_nuc = 100); at R_nuc = 250 pc median G = 13 (3e10), 272 (3e11) | threshold location approximate: 20 trees, cells of 25 to 50 pc and 0.1 in sigma_j, mass-dependent; the R_nuc = 250 pc case is a labelled experiment, not the model |
| S2 | In the R_nuc = 250 pc test light seeds remain far below the target | C7; Fig 2 dashed | highest 84th-percentile ratio 5.9e-4 (1e3 seed), 2.93 dex below f_BH = 0.5 | does not generalise to stronger accessibility (S3) |
| S3 | Accessible growth approaches a characteristic BH-to-stellar mass ratio of order 0.1 to 0.2 in the explored models, with weak dependence on seed mass; the critical seed falls by up to a factor of about 9 but remains supply-limited | C9, C9b | sigma_j = 1.5, R_nuc = 100: 1e3 seed 0.072 [0.052, 0.094] (3e10), 0.092 [0.069, 0.12] (3e11); 1e7 seed 0.42 (3e10), 0.17 (3e11); critical seed ratio to fiducial 0.14 (3e10), 0.11 (3e11); peak Mdot_acc/Mdot_Edd at the critical seed at most 0.48 | **Not established as a universal saturation scale.** It may reflect the momentum-driven feedback prescription, the accessible angular-momentum distribution, the nuclear radius, the star-formation law, the range of sigma_j explored, or their interaction; accessibility range and feedback/nuclear prescriptions remain model-dependent. Do not write "BH growth saturates at 0.1 to 0.2 M_star". 60 trees (C9b) and 20 trees (C9); no cell has more than 5% of trees at f_BH = 0.5; Eddington-limited critical seeds only with eta_acc of 0.05 or more (C9), outside the physically motivated range at 0.5 |
| S4 | f_BH rescales the boundary and does not change its shape | C11 | f_BH = 0.1: x0.24, x0.23, x0.25; f_BH = 0.9: x1.58, x1.65, x1.59 (3e10, 3e11, 3e13) | near-linear because G_BH is close to 1 |

## 3. Robustness, model dependence and sensitivity

| ID | Claim | Evidence | Numbers to quote | Qualification |
|---|---|---|---|---|
| R1 | The starting epoch z_seed = 15 to 35 has negligible effect | C17 | paired ratio to z_seed = 25 within 0.4% in the median, 3.4% in the worst tree | fixed step; nothing happens before the host exists; the MVM starts with an empty galaxy at z_seed |
| R2 | A stochastic ensemble against a smooth mean assembly history changes the boundary by 6 to 10% | C17 | x0.897 (3e10), x0.908 (3e11), x0.939 (3e13) | the Fakhouri et al. (2010) fit is extrapolated beyond its calibration (z below about 2); only one alternative description |
| R3 | M_hot controls the high-mass end of the boundary and the kink | C17 | none at 3e10; at 3e12 to 3e13, x0.44 to 0.47 (M_hot halved), x2.0 to 2.2 (doubled) | the step is imposed; its location was varied, not its sharpness |
| R4 | Star-formation parameters shift the boundary by tens of per cent, through M_star,tot | C12, C17 | eps_sf x0.5 / x2: 0.77 / 1.20 (3e10), 0.69 / 1.37 (3e11); n_rd 5 / 20: x1.06 to 1.23 / x0.71 to 0.80; eta_SN and wind off: 0.67 to 1.74 | n_rd enters through the galaxy-scale clock, not through accessibility |
| R5 | The boundary depends on nuclear accessibility through sigma_j much more than through R_nuc | C12, C9b | R_nuc 50 / 200 pc: x1.10 to 1.11 / x0.87 to 0.89; sigma_j 0.75: x0.67 to 0.69; sigma_j 1.0: x0.39 to 0.40; sigma_j 1.5: x0.11 to 0.14 | different tree sets in C12 and C9b; R_nuc beyond 250 pc only at 20 trees |
| R6 | Accretion efficiency and radiative efficiency | C12 | eta_acc x0.1 / x10: 1.17 to 1.27 / 0.52 to 0.75; radiative efficiency 0.057 / 0.32: 1.04 to 1.05 / 0.87 to 0.91 | radiative efficiency is weak because the Eddington limit does not bind at the critical seed |
| R7 | The simplification from the exploratory model to the MVM changes the boundary by 0.6 to 3.4% | C15 | per-tree ratio 1.034 (3e10), 1.006 (3e11) | 20 trees, 401 steps, identical trees |

## 4. Numerical qualifications

| ID | Statement | Evidence | Numbers | Qualification |
|---|---|---|---|---|
| Q1 | Numerical uncertainty on M_seed,crit is a few per cent at 3e10 and 3e11 | C13 | sampling 0.001 to 0.008 dex (240 trees); 801 steps about 0.5%, 1% below the refined limit; dz halving +1.3%, +3.3%; M_res at most 2.4% | at 3e13 the time-step offset is about 4% and dz and M_res were not tested |
| Q2 | Nuclear share of stars, star formation before z = 10, G_BH - 1 and the feedback-limited fraction are resolution-dependent | C14 | e.g. 3e10 nuclear share 0.139, 0.069, 0.041, 0.028 at 201, 401, 801, 1601 steps | quote only as qualitative; not convergence claims |
| Q3 | The model is feedback-limited before z of about 10 and supply-rich after z of about 7 | C14, C16 | 38.4%, 43.7%, 26.3% of feedback steps at 801 steps | qualitative statement is stable across the resolutions tested; the percentages are not |
| Q4 | Tree-to-tree scatter | C2 | 16-84 width 0.142, 0.195, 0.057 dex at 3e10, 3e11, 3e13 (range 0.057 to 0.195) | define the measure; not comparable with the old paper's 0.12 to 0.31 dex |

## 5. Methodological limitation

| ID | Statement | Evidence | Qualification |
|---|---|---|---|
| M1 | A PCH08 comparison is not valid at the required resolution | C17; foraois `docs/PCH08_HIGH_Z_DIAGNOSTIC.md` | whether this is an implementation issue or a regime limitation is not established; do not claim algorithm agreement; the old "PCH08 within 0.02 to 0.14 dex" cannot be reproduced |

## 6. Claims we explicitly do not make

1. A prediction of how black hole seeds form, or of any seed-formation channel or mass.
2. A universal critical seed mass independent of the baryon-cycle model.
3. An Eddington-timescale argument: the Eddington cap does not bind at the fiducial critical seed.
4. That light seeds cannot grow: at the fiducial they do not, but a modest increase in accessibility moves them to Eddington-limited growth (S1).
5. That the boundary is unchanged by accessibility, or that accessibility cannot approach the target: at sigma_j of 1 to 1.5 the critical seed falls to 0.1 to 0.4 of its fiducial value and light seeds reach 0.07 to 0.12 of M_star,tot (S3).
5b. That 0.1 to 0.2 of M_star,tot is a physical saturation scale or ceiling for black hole growth (S3): it is a result of the explored models.
6. Physical bimodality of black hole growth: the bimodality in Figure 2 belongs to the accessibility test at high halo mass; the median drop there is not a decline in individual growth.
7. A converged prediction of nuclear or early star formation (Q2).
8. Tree-algorithm convergence: only a smooth mean history was compared, and PCH08 is not usable (R2, M1).
9. That the kink or the low stellar mass at 3e13 is an emergent black hole phenomenon: it follows from the imposed cold/hot step at M_hot (R3).
10. That the model has been comprehensively parameter-tested: z_reion and UV suppression, lambda, f_b, the NFW concentration and the sharpness of the cold/hot step are untested in the MVM.
11. That the R_nuc = 250 pc experiment is a second physical model or a candidate fiducial.
12. That 1e7 Msun is a seed-formation channel (it is a benchmark seed mass).

## 7. Figures and claims

| Figure | File | Claims it carries |
|---|---|---|
| 1 | `output/mvm_fig1_boundary.png` | H1, H2 (lower panel), R3 (kink location), Q4 (band) |
| 2 | `output/mvm_fig2_inverse.png` | H4, S1, S2 (dashed curves and band); caption must say the accessibility test is an experiment and that the median is not the individual growth |
| 3 | `output/mvm_fig3_hostbaryons.png` | H5 |

## 8. Corrections to the working sentences (found by checking against the data) and the agreed wording

The two working sentences of 2026-09-20 were overstated; both were tested by `c17_accessibility_reach.log` (C9b). The wording below was agreed on 2026-09-20.

* Replaced: "the nuclear-accessibility prescription ... does not materially change the critical boundary over the tested accessibility range."
  It holds for R_nuc of 50 to 250 pc (x0.82 to 1.11) but not for sigma_j: the boundary falls to x0.39 to 0.40 at sigma_j = 1.0 and x0.11 to 0.14 at sigma_j = 1.5 (S3, R5).
  **Agreed:** *Increasing nuclear accessibility can substantially reduce the critical seed mass, particularly through the width of the angular-momentum distribution, but the critical solution remains supply-limited over the tested range.*
* Replaced: "increased accessibility alone does not reach f_BH = 0.5."
  True for the R_nuc = 250 pc test (light seeds at least 2.9 dex short) but not for stronger accessibility: at sigma_j = 1.5 light seeds reach 0.07 to 0.12 of M_star,tot, a factor of about 4 to 7 short of 0.5, and a 1e7 Msun seed reaches 0.42 at 3e10 (5% of trees at or above 0.5).
  **Agreed:** *Increasing nuclear accessibility can drive orders-of-magnitude growth of light seeds, but in the explored models the resulting BH mass saturates at a substantial fraction of the stellar mass rather than generically reaching the adopted M_BH/M_star = 0.5 target.*
  (Read "saturates" as a description of the explored models, not a claim of a universal ceiling: see S3.)
* What this says: accessibility matters, but it matters by changing how much black hole growth can be achieved from a given seed, not by turning the critical solution into a conventional Eddington-growth problem.

## 9. Open narrative decisions

1. Which figures go into the manuscript, and whether the inverse figure (Fig 2) is presented as central or interpretive.
2. S3 is retained as a model result with the qualification beside it. Open: whether the 0.1 to 0.2 scale survives one carefully chosen physical stress test (for example the feedback and star-formation prescriptions at sigma_j = 1.5), which would decide whether it belongs in a paper at all.
3. How to describe M_hot: as a modelling dependence of the high-mass boundary, with the location varied but not the sharpness.
4. The old-model numbers (crib sheet, section 3) must not appear.
