# MVM v2: revised reservoir model, validation and production

Written 2026-11-14. Frozen model: `ashvini/reservoir_stock.py` (SHA-256 31c701dd..., unchanged). Revised model:
`ashvini/reservoir_stock_v2.py`. Rationale and the evidence for each change: `model-verification-audit-2026-11-14.md` in the
paper folder (High-Redshift-Black-Hole-Growth).

## What changed (each is a switch; legacy values reproduce the frozen model exactly, asserted in the tests)

| Switch | Frozen | Revised | Why |
|---|---|---|---|
| `dm_nuclear` | False | True | NFW dark matter inside R_nuc added to M_enc for the nuclear capture threshold and free-fall time. The omitted term exceeds the nuclear baryons of a light seed by orders of magnitude and changes the light-seed regime. |
| `inflow` | `rate` (GRUMPY spline rate times dt) | `conservative` (mass the halo gains on the grid times the suppression factors times f_b) | The rate estimator over-delivers baryons on merger trees: cumulative inflow over f_b M_halo(z=5) was 1.18 (median, 3e10), up to 1.6 for individual trees; exactly 1.00 now at 3e10 and 0.93 at 3e11, where the cold/hot factor legitimately suppresses. |
| `cosmology` | None (Planck18 tables for time and sigma_h, 67.8/0.308 for halo structure) | `PaperCosmology` (one flat LCDM, closed-form time) | The two differ by 4.4 Myr in the z = 25 to 5 duration. |
| `n_sub` | 1 | 8 (variant F only) | Consumption and feedback were first-order accurate in the step; with 25 to 45 per cent of live steps gas-capped by feedback, the nuclear stellar fraction and light-seed growth factor drifted with the step. |

Unchanged: the order of operations in a step, all parameters, the galaxy cut, the nuclear transfer with history, the exact
Eddington step, proportional sharing, the form of both outflows, the critical-seed definition and bisection.

## Tests (`tests/test_reservoir_stock_v2.py`, 10 tests, about 8 s)

Legacy switches reproduce the frozen model bit for bit (smooth and merger-like histories, fiducial and high-delivery
parameters); the closed-form cosmology matches astropy to 1e-7 and inverts; the conservative scheme delivers exactly f_b
times the halo gain; baryon conservation to 1e-12 for every switch combination; nuclear dark matter shortens the nuclear
free-fall time and raises delivery; one sign change of F on a smooth history; n_sub=1 is the frozen step, n_sub>1 conserves
baryons and approaches the refined-grid answer.

## Local validation (small ensembles, indicative only; the production run replaces these numbers)

Paired median of the per-tree ratio of the critical seed to the frozen model, same trees:

| M_halo(z=5) | trees | dark matter only | conservative inflow only | + one cosmology | + dark matter (E) |
|---|---|---|---|---|---|
| 3e10 | 30 | 0.975 | 0.868 | 0.869 | 0.841 |
| 3e11 | 20 | 0.997 | 0.868 | 0.870 | 0.867 |
| 3e13 | 20 | 0.997 | 0.922 | 0.923 | 0.919 |

* The cosmology change moves the critical seed by under 0.1 per cent at these masses (critical seed is in the G_BH of about 1
  regime, set by the stellar mass, not by the exponential growth time).
* Dark matter barely moves the critical seed but moves light seeds: at 3e11 a 1e3 Msun seed has median G_BH of 298 in the
  revised model against 1.0 in the frozen model (the paper's separate test gave 322); the black hole to stellar ratio of that
  seed rises from 3e-7 to 1e-4, still far below the target.
* Step convergence (`diagnostics/v2_timestep_convergence.py`, 8 trees, 3e10, nodes:n_sub): the critical seed agrees to 0.5
  per cent between 401 and 1601 nodes with or without sub-stepping. Nuclear stellar fraction: 0.074 / 0.041 / 0.029 for
  401 / 801 / 1601 nodes (n_sub 1), but 0.024 at 801:4, 0.022 at 801:16 and 0.022 at 1601:16. Light-seed G_BH: 167 / 145 / 133
  (n_sub 1), 132 at 801:4, 132 at 801:16, 129 at 1601:16. Fraction of live steps limited by gas: 0.42 / 0.38 / 0.35 (n_sub 1)
  against 0.30 (4 sub-steps), 0.12 (16) and 0.04 (1601:16): the gas-limited fraction is a step artefact.
* Single crossing of F(M_seed): 0 of 30 trees at 3e10 (17-point grid); the production run checks every tree on 33 points.

## Production on OzSTAR/NT

```
git pull                                   # on the cluster, after pushing from the Mac
export FORAOIS_ROOT=/path/to/foraois       # v0.1.2 or later
mkdir -p /fred/oz009/cpower/Ashvini/logs /fred/oz009/cpower/Ashvini/mvm_v2
sbatch --array=0-12 hpc/submit_mvm_v2.sbatch
python scripts/paper_figures/merge_mvm_v2_production.py /fred/oz009/cpower/Ashvini/mvm_v2
```
One task per halo mass (13); each runs the tests first, then samples 240 Zhang and Hui trees (stored in `trees_M<idx>.npz`;
the numba sampler is not seed-reproducible, so these files are the record), and runs six variants on them: A frozen
(reproduction), B dark matter only, C conservative inflow only, D conservative inflow plus one cosmology, E all three physics
changes, F = E with 8 sub-steps (the revised model). F also gets the fixed-seed fates (seeds 1e2, 1e3, 2e5, 1e7 in the
fiducial, R_nuc = 250 pc and high-delivery configurations) and the single-crossing check for every tree. The JSON records the
hashes of the frozen module, v2 module and script and the git state of both repositories. Estimated 20 to 40 minutes per task
from the local timings (30 trees: about 14 s per variant at 801 nodes, 68 s for F including fixed seeds); the 8 hour request is generous.

## Not changed, and open

* Tree sampler: the Zhang and Hui binary-per-step sampler is not step-compliant at these masses (about 6 to 7 per cent high
  at 3e10 and 3e11). `foraois` has a Brownian-bridge constrained branch sampler, but it conditions on reaching a specified
  mass at a specified redshift, so it is not a drop-in replacement for unconstrained main-progenitor histories.
* Assumptions kept deliberately and to be stated in the paper: no atomic-cooling threshold for gas entering the galaxy at
  z above 10, stellar mass is mass formed (no return fraction), instantaneous stellar outflow, fixed NFW concentration, the
  memoryless galaxy cut.

## Tree algorithm (added 2026-11-14; tests run with FORAOIS_ROOT set, otherwise the PCH08 test is skipped)

The production script now defaults to `--tree-algorithm pch08`: the main progenitor from compliant, adaptive Parkinson, Cole
& Helly (2008) trees, with the step chosen per branch so that the expected number of splits per step is at most
`--target-nupper` (0.1), as in the Stochastic_Occupation paper. Implementation: `ashvini/pch08_main_branch.py` (secondary
branches are discarded, which does not affect the main branch because the step is chosen per branch). Output is on the same
uniform redshift grid as before (`--dz` is now only the output spacing), so the variants A to F run unchanged. Trees are seeded
as 1000*(mass index + 1) + tree index and are reproducible; `trees_M<idx>.npz` also stores `tree_stats` (steps, resolved
splits, and steps at the step floor; the last must be zero).

`--tree-algorithm zh` keeps the Zhang & Hui sampler at fixed `--dz` as a bracketing comparison (not compliant at these masses).

Local checks (M0 = 3e10, 4 trees): reproducible, no steps at the floor, 0.8 to 1.5 million steps per tree; median M(z)/M0 is
0.13 at z = 10 and 0.007 at z = 15 (Zhang & Hui: about 0.044 at z = 10). Cost per tree at z = 5 to 25 is seconds (3e10) to about
half a minute (3e13) of compiled time, plus tens of seconds of table set-up per task. The sbatch script now uses 4 CPUs.
