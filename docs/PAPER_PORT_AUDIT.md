# Porting the discreteness paper into Ashvini: pre-flight audit (2026-11-04; for agreement, no code written)

## 1. Objective
Ashvini becomes the single codebase for the paper "discrete versus continuous supernova feedback in the smallest star-forming haloes". Every number and figure is regenerated with Ashvini, the fork is tagged, and Stochastic_Occupation (private) is retired. Data are generated with Ashvini and available on request; a Zenodo deposit of summary arrays is optional.

## 2. What Ashvini has and what the paper needs
Has (M0 to M6, pushed to the fork): event integrator with exact substeps, stochastic and K-quadrature smooth samplers, energy and momentum wind models, stellar-lifetime delays, chemistry, and the options cooling threshold, inflow lag, recycling, energy bank, access fraction, coupling switch, wind travel time. Agreement with Stochastic_Occupation's smooth control about 1e-3.

Missing:
| Id | Item | Paper section it feeds |
|---|---|---|
| M7 | Supply read from a cumulative table of f_b f_ret dM_h per trajectory, identical in both treatments (`events.supply = "table"`; the old start-of-substep rule kept as `"start"`) | Appendix A; every number |
| M8 | Packet sampler (k supernovae per event at 1/k the rate) and switch-off samplers (fixed or Poisson counts; random or stratified delays) | event-size sweep, fixed counts and delays |
| M9 | Diagnostics: young-star gas removal (rate lambda, 4 Myr window), lifetime scale (zero lifetimes), gas limit off, and cumulative records of ejected mass, potential removal and realised removal on the output grid | radiation test, delay test, gas-limit test, surplus mechanism |
| M10 | Ensemble runner: independent random draws for every run in every tree (shared draws as an option), 5 Myr and 1 Myr output grids, resume-safe shards, OzSTAR Slurm scripts | all ensembles |
| M11 | Paper scripts and a README mapping each figure and number to one command | reproducibility |
| M12 | Regenerate everything with Ashvini; tag a release on the fork; archive Stochastic_Occupation | submission |

## 3. Known differences that limit parity
- Ashvini's tabulated Hubble time and constants, and the substep cap (min of 0.1 t_ff and 1 Myr against 0.1 t_dyn and 1 Myr): about 1e-3.
- Random draw order differs between the two codes, so per-tree stochastic output will not match; parity is statistical.
- The paper's control uses K = 200 (Ashvini default 100); the runner sets K = 200 and a test pins it.
- The stored M4/M5 reference numbers used the start-of-substep supply. They are regenerated against Stochastic_Occupation run with `STOCHOCC_ACC=1`.
- The supply convention at mergers (suppression at the interval midpoint mass) is a choice common to both treatments; it is recorded in the documentation.

## 4. Acceptance targets (registered before any porting)
- **A1** Supply conservation: total baryons (gas + stars + ejected) agree between discrete and smooth runs to 1e-7 of the total, for 20 trees per mass (now 8e-4 with the old rule).
- **A2** Smooth control against Stochastic_Occupation with `STOCHOCC_ACC=1`: stars and gas at z = 5 within 1e-3 on 5 trees at each of 1e7, 1e8, 1e9.
- **A3** Discrete ensemble, 100 trees x 50 runs, independent draws, chemistry on, same trees as `results_ens_ts_acc/fid`: median offsets at z = 5 within two bootstrap standard errors of 148.3, 135.4, 59.9, 14.3 and 0.1 per cent (5e6, 1e7, 3e7, 1e8, 1e9); offsets at 20 and 50 Myr within 3 points of the stored values.
- **A4** Mechanism metrics (nochem, 100 trees x 25 runs): late loading ratio within 0.05 of 0.29, 0.65, 0.90 (1e7, 3e7, 1e8) and 1.00 at 1e9; ejection deficit equals the ejected-mass deficit to 1e-3.
- **A5** Sampler tests: packet sampler conserves expected energy per unit stellar mass for k = 0.01 to 10; switch-off modes reproduce the mean count and the quantile delays.
- **A6** Event-size sweep, fixed counts and delays, threshold, momentum, radiation and zero-lifetime runs agree with the corrected OzSTAR repeats of Stochastic_Occupation (results_acc) within sampling error, as each set arrives.
- **A7** Same seed gives identical output; different seeds give independent output.

## 5. Process
- Work on a branch `paper-port` in the fork, one commit per milestone, merged to `main` when A1 to A7 hold. Nothing is pushed to `origin` (Anand's repository).
- The uncommitted files already in the Ashvini working tree are not touched or committed.
- The OzSTAR repeats of Stochastic_Occupation continue; they are the parity reference.
- Paper numbers are replaced by the Ashvini-generated ones only after A3 to A6 pass. Where they differ from the current text by more than sampling error, that is reported.

## 6. Data
Generated with Ashvini; available on request. A Zenodo deposit of summary arrays (offsets per tree and time series, not the raw per-run files) is optional.

## 7. Decisions assumed (to confirm)
1. Default supply becomes `"table"` with the goldens regenerated; `"start"` stays available.
2. Statistical, not per-tree, parity for stochastic runs.
3. The fork is tagged at M12 (name to be agreed).
4. Stochastic_Occupation is archived as a tarball and no longer developed.
