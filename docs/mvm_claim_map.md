# MVM claim map

> **The central result is a conditional statement about the initial BH mass required to reach a specified
> M_BH/M_star,tot target within this baryon-cycle model, not a prediction of how black-hole seeds are formed.**

A map of what the frozen model supports, what it only qualifies, and what must not be claimed. Every row points to an entry of
`docs/mvm_numerical_crib_sheet.md` (C-numbers) and, where relevant, a figure. Numbers are quoted from the crib sheet; none is new here
except where marked. Nothing in this file is manuscript text.

## 0. Provenance chain

The foraois hashes in this section are those of the foraois history after its message-only rewrite; the earlier hashes are listed in `docs/PRODUCTION_PROVENANCE.md`, section 2.

| Step | Repository | Commit | Content |
|---|---|---|---|
| 1 | foraois | `ec1a66f` | collapse barrier normalised to the P(k) redshift, with tests. The production trees were generated from `4370a4a` plus this change while it was still uncommitted (inferred from file times); `ec1a66f` committed it afterwards and is the recommended reproducibility pin, not the run-time commit (`docs/PRODUCTION_PROVENANCE.md`) |
| 2 | foraois | `767398d` | PCH08 diagnostic note, reproduction script, ROADMAP entry (no source change) |
| 3 | Ashvini | `695b114` | frozen MVM, pre-MVM reference, production ensemble, figures, original diagnostic scripts (its logs were omitted by mistake; see 4) |
| 4 | Ashvini | `8ffb61e` | C17 checks, PCH08 diagnosis, all 27 diagnostic logs, crib-sheet update; corrects 3 (which did not contain the logs) |
| 5 | Ashvini | `fdb2138`, `59418eb` | claim map, accessibility-reach diagnostic and log, crib-sheet C9b; agreed wording |
| 6 | Ashvini | (the C18 commit) | C18 saturation/mechanism tests, phase-space diagnostic, efficiency grid, their logs, data and figures; efficiency literature check; claim-map and crib-sheet updates |

* Between foraois `4370a4a` (the last commit before this work) and `767398d`, the only source change is `src/foraois/cosmo_utils.py` (the barrier fix). The production run used `4370a4a` plus that change while it was still uncommitted (inferred from file times); it was committed afterwards as `ec1a66f`, which is the recommended reproducibility pin, not the run-time commit. The exact run-time commit is not known (`docs/PRODUCTION_PROVENANCE.md`).
* Frozen model: `ashvini/reservoir_stock.py`, sha256 `31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63`, recorded in
  `scripts/paper_figures/output/mvm_production_results.json` (2026-09-20T15:44:29). `c17_diagnostics.py`, `c17_accessibility_reach.py`, `c17_smooth_assembly.py` and the three `c18_*` scripts assert this hash before running; the `c17_pch08_check*.py` scripts and the `mvm_*.py` diagnostics do not (`docs/PRODUCTION_PROVENANCE.md`, section 3).
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
| S3 | In the accessible, supply-limited regime the model rapidly reduces seed-mass dependence, because black hole growth and nuclear star formation compete for the same gas reservoir. In the no-feedback limit the resulting M_BH/M_star,tot scales approximately with the ratio of the adopted BH accretion and star-formation efficiencies; AGN feedback further suppresses it. The critical seed falls by up to a factor of about 9 but remains supply-limited | C9, C9b, C18 | 3e10 and 3e11 only; five decades of seed mass compress to a factor 5.3 (3e10) and 2.05 (3e11) by z = 5 (seeds 1e2 to 1e7), or 1.2 and 1.4 for the three lightest; no-feedback R/(eta_acc/eps_sf) = 0.6 to 0.99 (3e10), 0.3 to 0.94 (3e11), log-log slope 0.79 to 0.98; AGN on: slope 0.56 to 0.69 and R/(eta_acc/eps_sf) of 0.08 to 0.44 (3e10, seed 1e3); critical seed ratio to fiducial 0.14 (3e10), 0.11 (3e11) at sigma_j = 1.5 | **The resulting M_BH/M_star of about 0.1 to 0.2 at the fiducial efficiencies is a model-dependent consequence of the adopted nuclear sink efficiencies and feedback prescription, not an independent prediction of a characteristic BH-to-stellar mass scale.** Finite-time relaxation to z = 5, not an asymptote; not at 3e13 (max/min across seeds 51); two masses, 100 trees, one accessibility cell plus a fiducial control |
| S4 | f_BH rescales the boundary and does not change its shape | C11 | f_BH = 0.1: x0.24, x0.23, x0.25; f_BH = 0.9: x1.58, x1.65, x1.59 (3e10, 3e11, 3e13) | near-linear because G_BH is close to 1 |
| S5 | The restoring tendency of R: d ln R/dt is positive below and negative above a moving zero crossing in the accessible regime, and has no crossing at the fiducial accessibility | C18 (phase space); `sat_fig7_phase.png` | median crossings (z <= 7): 0.108 (3e10), 0.154 (3e11) with AGN on; 0.264, 0.267 with AGN off; 0.13, 0.164 with the stellar wind off; drift from 0.048 to 0.108 (3e10) with time; no crossing at fiducial accessibility | describe as restoring behaviour / finite-time relaxation, not an attractor: it follows from d ln R/dt = Mdot_BH/M_BH - Mdot_star/M_star with both drawing on the same reservoir; the z > 10 band is in the resolution-dependent feedback-limited regime and is not read |

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
| R8 | The accessible-regime ratio depends on the adopted efficiencies: the no-feedback floor scales as about eta_acc/eps_sf, but not exactly as the ratio alone | C18 (efficiency grid); `sat_fig8_eta_eps.png`; `docs/mvm_efficiency_literature_check.md` | over a factor of about 80 in the ratio, log-log slope 0.79 to 0.98 (no feedback); the two grid points with equal ratio 0.667, (0.005, 0.0075) and (0.02, 0.03), give R = 0.546 and 0.401 (3e10, seed 1e3), 0.373 and 0.241 (3e11) | absolute efficiencies matter because they change the gas trajectory and the Eddington limitation (for seed 1e3, 20 to 62% of accreting steps are Eddington-capped in every cell of the grid, carrying 5 to 37% of the BH mass gained without feedback); AGN feedback adds a further dependence; heavy seeds at low ratio retain their seed memory. eps_sf = 0.015 is inherited (applied per local free-fall time in the MVM, per Hubble-scale clock in the source) and eta_acc = 0.005 has no independent calibration in the sources checked: the ratio is not independently motivated |

## 4. Numerical qualifications

| ID | Statement | Evidence | Numbers | Qualification |
|---|---|---|---|---|
| Q1 | Observed ensemble-to-ensemble variation of the fiducial median critical seed is about 5% at 3e10 (independent ensembles of 60 to 240 trees give medians from 9.0e7 to 9.9e7); it is not established at other masses | C13; `docs/PRODUCTION_PROVENANCE.md` | time step: 801 nodes about 0.5% (3e10) and 1% (3e11) below the refined limit; dz halving +1.3%, +3.3%; M_res at most 2.4% (earlier repeated draws gave 0.001 to 0.008 dex, smaller than the observed variation) | the Zhang and Hui trees at dz = 0.05 are outside the practical single-split timestep-compliance regime, so the dz and M_res tests are not convergence tests of the tree construction and the trees are not timestep-converged; the effect on the black-hole calculation is not quantified; at 3e13 the time-step offset is about 4% and dz and M_res were not tested |
| Q2 | Nuclear share of stars, star formation before z = 10, G_BH - 1 and the feedback-limited fraction are resolution-dependent | C14 | e.g. 3e10 nuclear share 0.139, 0.069, 0.041, 0.028 at 201, 401, 801, 1601 steps | quote only as qualitative; not convergence claims |
| Q3 | The model is feedback-limited before z of about 10 and supply-rich after z of about 7 | C14, C16 | 38.4%, 43.7%, 26.3% of feedback steps at 801 steps | qualitative statement is stable across the resolutions tested; the percentages are not |
| Q4 | Tree-to-tree scatter | C2 | 16-84 width 0.142, 0.195, 0.057 dex at 3e10, 3e11, 3e13 (range 0.057 to 0.195) | define the measure; not comparable with the old paper's 0.12 to 0.31 dex |

## 5. Methodological limitation

| ID | Statement | Evidence | Qualification |
|---|---|---|---|
| M1 | A PCH08 comparison is not valid at the required resolution | C17; foraois `docs/PCH08_HIGH_Z_DIAGNOSTIC.md` | the foraois analysis attributes the near-deterministic PCH08 histories to timestep non-compliance (expected splits per step far above about 0.1 at fixed dz and small M_res/M0); this identifies the cause, not the size of any effect on the black-hole calculation; the same single-split criterion applies to the Zhang and Hui trees at the production settings; do not claim algorithm agreement; the old "PCH08 within 0.02 to 0.14 dex" cannot be reproduced |

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
13. That M_BH/M_star of 0.1 to 0.2 is a universal physical ceiling, an attractor, or an emergent characteristic scale. It is the model's nuclear supply partition (about eta_acc/eps_sf) reduced by feedback.
14. That eta_acc and eps_sf are independently calibrated, or that the black hole-to-stellar partition is a prediction. It is conditional on the adopted efficiencies.
15. That the critical seed is close to the target unconditionally: it requires poor accessibility (the fiducial) or eta_acc/eps_sf times the feedback factor below f_BH.

## 7. Figures and claims

| Figure | File | Claims it carries |
|---|---|---|
| 1 | `output/mvm_fig1_boundary.png` | H1, H2 (lower panel), R3 (kink location), Q4 (band) |
| 2 | `output/mvm_fig2_inverse.png` | H4, S1, S2 (dashed curves and band); caption must say the accessibility test is an experiment and that the median is not the individual growth |
| 3 | `output/mvm_fig3_hostbaryons.png` | H5 |
| C18 | `output/sat_fig1_trajectories.png`, `sat_fig2_final_vs_seed.png`, `sat_fig3_fmom.png`, `sat_fig4_fmom_trajectories.png`, `sat_fig5_all_accessible.png`, `sat_fig6_test4.png`, `sat_fig7_phase.png`, `sat_fig8_eta_eps.png` | S3, S5, R8 (diagnostic figures; not proposed for the manuscript without a decision) |

## 8. Corrections to the working sentences (found by checking against the data) and the agreed wording

The two working sentences of 2026-09-20 were overstated; both were tested by `c17_accessibility_reach.log` (C9b). The wording below was agreed on 2026-09-20.

* Replaced: "the nuclear-accessibility prescription ... does not materially change the critical boundary over the tested accessibility range."
  It holds for R_nuc of 50 to 250 pc (x0.82 to 1.11) but not for sigma_j: the boundary falls to x0.39 to 0.40 at sigma_j = 1.0 and x0.11 to 0.14 at sigma_j = 1.5 (S3, R5).
  **Agreed:** *Increasing nuclear accessibility can substantially reduce the critical seed mass, particularly through the width of the angular-momentum distribution, but the critical solution remains supply-limited over the tested range.*
* Replaced: "increased accessibility alone does not reach f_BH = 0.5."
  True for the R_nuc = 250 pc test (light seeds at least 2.9 dex short) but not for stronger accessibility: at sigma_j = 1.5 light seeds reach 0.07 to 0.12 of M_star,tot, a factor of about 4 to 7 short of 0.5, and a 1e7 Msun seed reaches 0.42 at 3e10 (5% of trees at or above 0.5).
  **Agreed:** *Increasing nuclear accessibility can drive orders-of-magnitude growth of light seeds, but in the explored models the resulting BH mass saturates at a substantial fraction of the stellar mass rather than generically reaching the adopted M_BH/M_star = 0.5 target.*
  (Read "saturates" as a description of the explored models, not a claim of a universal ceiling: see S3.)
* After C18 the second sentence has a mechanism: the light-seed ratio approaches a level set by the competition between black hole growth and nuclear star formation (S3, R8), which for the fiducial efficiencies lies below f_BH. "Saturates" describes the explored models and is conditional on eta_acc/eps_sf.
* What this says: accessibility matters, but it matters by changing how much black hole growth can be achieved from a given seed, not by turning the critical solution into a conventional Eddington-growth problem.

## 9. Open narrative decisions

1. Which figures go into the manuscript, and whether the inverse figure (Fig 2) is presented as central or interpretive.
2. S3 is retained as a mechanism, not a physical scale (see the wording in S3 and the not-claimed items 13 to 15). The stress test that could have made it a result of the paper has been run (C18) and reduced it to the model's efficiencies; no further sweeps are planned.
3. How to describe M_hot: as a modelling dependence of the high-mass boundary, with the location varied but not the sharpness.
4. The old-model numbers (crib sheet, section 3) must not appear.
5. How prominently the efficiency dependence (R8) and its provenance appear, given that eta_acc has no independent calibration in the sources checked (`docs/mvm_efficiency_literature_check.md`).

## 10. Proposed paper pitch (2026-09-20) against the evidence

Status: a proposal from the author, assessed here; not adopted wording. The manuscript is untouched. The pitch reframes the paper from a critical-seed calculation to a mechanism paper: growth is gated first by nuclear accessibility, then becomes a competition with star formation for a shared reservoir, and the critical seed follows from that. Related working notes: `docs/mvm_feeding_efficiency_discussion.md` (efficiency provenance; torque-equation comparison) and `docs/mvm_discussion_alignment.md` (old-model sentences that conflict with the MVM).

**Pitch sentence, tightened to the evidence:** *In a baryon-cycle model of early black-hole growth, whether a seed grows is set first by whether gas can reach the nucleus. Once it can, growth becomes a competition with star formation for a shared reservoir, which erases much of the seed memory within the available time at the halo masses tested.*

The agreed central sentence still governs the framing: the result is a conditional statement about the initial BH mass required to reach a specified M_BH/M_star,tot target within this baryon-cycle model, not a prediction of how black-hole seeds are formed.

### 10.1 Pitch statements

| Pitch statement | Status | Evidence | Constraint on the wording |
|---|---|---|---|
| Growth begins with access, not with the seed | Supported as a statement about the model | H2, H3, R5, S1; light-seed transition near R_nuc 125 to 150 pc (sigma_j 0.5) or sigma_j 0.7 to 0.8 (R_nuc 100 pc) | The accessibility parameterisation (lognormal j, memoryless galaxy cut, sigma_j, R_nuc) has no independent calibration; the paper cannot say which regime real z > 7 nuclei occupy |
| The critical seed grows by only a few per cent | Supported at the fiducial only | H2, H3; G_BH 1.03 to 1.08 | Quote "close to unity" (Q2: G_BH - 1 shrinks with dt); the fiducial is Regime I |
| The critical seed is set mostly by the stellar mass the host builds | Supported, and near-definitional at the fiducial | H2: M_seed,crit = f_BH M_star,tot / G_BH | If the BH barely grows, the boundary approaches the target BH mass by construction; say so explicitly |
| Regime I (inaccessible): the BH remembers its seed | Supported | H4: ratio about M_seed / M_star,tot; 4.0e-7 (1e2), 4.1e-6 (1e3), 4.2e-2 (1e7) at 3e10 | Frozen fiducial, chosen as such; the light-seed threshold is a result of the model, not a physical threshold |
| Regime II (accessible): small seeds grow by orders of magnitude | Supported | C9, C18 test 1: G median 5.1e4 (1e2 seed, 3e10), 1.0e6 (3e11) at sigma_j = 1.5, R_nuc = 300 | Transition map rests on the 20 and 60-tree grids (C9, C9b), not the production ensemble |
| Regime III: shared-reservoir competition, finite-time relaxation | Supported in mechanism (S5, R8) | C18: d ln R/dt = Mdot_BH/M_BH - Mdot_star/M_star; restoring crossing at 0.108 (3e10) and 0.154 (3e11) for z <= 7 with AGN on | One accessibility cell (sigma_j = 1.5, R_nuc = 300 pc), two masses, 100 trees; z = 5 is not an asymptote; not "attractor" |
| Seeds five decades apart end within a factor of a few | Partly | S3, C18 test 1 | Three lightest seeds (1e2 to 1e5): x1.2 (3e10), x1.4 (3e11). All four seeds (to 1e7): x5.3, x2.05. Not at 3e13 (x51), where hot-mode supply limits growth; the 1e7 seed is still falling at z = 5 |
| The seed sets how long the system takes to forget it | Consistent with, not measured | dispersion across seeds, median ratio max/min: 1.45e3 (z = 10), 48.6 (z = 7), 5.3 (z = 5) at 3e10; 1.06e3, 12.3, 2.2 at 3e11 | No relaxation timescale was measured; do not quote one |
| M_BH/M_star of 0.1 to 0.2 is not a universal scale | Supported | S3, R8 | The ratio is approximately eta_acc/eps_sf times a feedback factor, and not a function of the ratio alone (equal-ratio pairs differ) |
| The fiducial M_BH/M_star is a high-side estimate | **Not yet** | `docs/mvm_feeding_efficiency_discussion.md` section 2: eta_acc = 0.005 is about x2.4 above the eps_T = 5 equivalent for one nucleus | Rests on a torque equation taken from the author's message, not checked by me, and one nucleus; high-side only in eta_acc at fixed eps_sf, and eps_sf is uncertain in both directions; keep out of the abstract until eq. 2 is checked |
| Early BH growth is controlled by two bottlenecks (access, then competition) | Supported as a statement within the model | as above | Keep the conditional form ("within this baryon-cycle model"); omitted: super-Eddington growth, BH mergers, nuclear dark matter, recycling |
| An observed high M_BH/M_star is not uniquely diagnostic of a heavy seed | Conditional on the regime | S3, R8, H4 | If nuclei are inaccessible, only heavy seeds are visible; if accessible, the ratio reflects eta_acc/eps_sf and feedback more than the seed. Observations do not yet distinguish these. Light seeds do not reach ratios of 0.5 in the explored models, so the inference cuts both ways. The model ends at z = 5 and uses M_star,tot, not an aperture stellar mass |

### 10.2 What the paper must say plainly if this pitch is adopted

1. The fiducial is in the inaccessible regime by construction. The critical seed there is about 4e3 to 5e4 times the host's baryons at first resolution (H5: medians 5.1e3, 4.9e4, 3.8e3), so no physically formed seed reaches the target in that configuration. This is a property of the frozen accessibility, not a result about nature.
2. Regime placement of real nuclei is unconstrained. sigma_j, R_nuc and eta_acc are uncalibrated (R8; provenance note).
3. The accessible-regime evidence is one cell, two masses, finite time. It does not hold at 3e13.
4. The critical seed is insensitive to eta_acc (R6: x1.17 to 1.27 for eta_acc x0.1; tested only to x0.1) while the accessible-regime ratio is not (R8).

### 10.3 Evidence gaps for the pitch, in order of value

* Seed convergence at an intermediate accessibility (for example sigma_j = 1.0), which would show whether Regime III is a single cell or a range. Optional; not run; the author's instruction is not to expand the sweeps.
* Verification of the torque equation (Angles-Alcazar et al. 2017, eq. 2 and the eps_T range), on which the "high-side" statement and the feeding-efficiency discussion rest.
* Evaluation of the torque expression along stored MVM nuclear states against the model's Mdot_acc. This needs no new simulation. It depends on the previous item.

The transition map (Regime II) would need to be drawn from the 60-tree reach grid (C9b) or rerun at production tree number if it is to be a headline figure.

### 10.4 Open decisions from the pitch

* Whether the paper is reframed as a mechanism paper or remains a critical-seed paper with a mechanism section (affects which of Figs 1 to 3 and `sat_fig1, 2, 7, 8` are headline).
* Title: the candidate "Black-hole growth in the early Universe: nuclear accessibility and competition for a shared gas reservoir" is acceptable; a subtitle naming the model ("in a baryon-cycle model") avoids promising a result about nature.
* Whether the JWST implication is stated at all before the regime placement can be discussed.
