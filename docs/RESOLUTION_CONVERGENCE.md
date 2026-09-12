# Resolution/timestep convergence of the z=0 stellar-to-halo mass relation

## Summary

While computing a stellar-to-halo mass relation (SHMR) with Ashvini run on
`foraois`-generated PCH08 merger trees, the low-mass end (`M_halo(z=0) <~
1e11 Msun`) turned out to be sensitive to two tree-generation parameters
that `run_params.yaml`/`build_forest_live` leave at fixed defaults:

- `dz`, the tree-growth timestep in redshift (default `0.05`)
- `m_res`, the mass resolution below which a progenitor is not tracked
  (default `1e-3 * mass_bin`)

A one-parameter-at-a-time test suggested both mattered independently, in
*opposite* directions (finer `dz` alone lowered the low-mass `M_star`;
finer `m_res` alone raised it). A joint grid test resolved this: **`dz`
has essentially no independent effect once `m_res` is resolved finely
enough.** The three `m_res_fraction=1e-5` curves at `dz=0.05/0.01/0.005`
agreed to within ~2% at every mass bin tested (`1e8`-`1e12 Msun`). The
earlier apparent `dz` sensitivity was a confound of testing `dz` only at
the coarse default `m_res_fraction=1e-3`, which is itself far from
converged.

**`m_res` is the real convergence axis, and it is not yet fully converged
even at the finest value tested.** Refining `m_res_fraction` by one decade
(`1e-3 -> 1e-4 -> 1e-5`) changed the lowest-mass bin's median `M_star` by
roughly 7-8x *each time*, with no sign of the change shrinking by the
finest value tested (`1e-5`). This means:

- `run_params.yaml`'s default `m_res_fraction=1e-3` is **not usable** for
  a low-mass (`<~1e11 Msun`) SHMR result -- it is not a "slightly coarse
  but roughly right" choice, it is off by close to two orders of
  magnitude from wherever the true converged answer lies.
- `m_res_fraction=1e-5` is a better floor but has **not been verified
  converged** -- a follow-up test extending to `1e-6`/`1e-7` is needed
  before trusting an absolute low-mass `M_star` value, only relative
  comparisons at matched (and adequately fine) resolution.
- **Follow-up (resolved)**: the initial `z_max=15` used for this test
  turned out to have a real, separate artifact -- at the finest `m_res`
  tested (`1e-5*M0`), median halo *formation redshift* (the first step a
  halo's tracked mass exceeds `M_res`) was pinned essentially exactly at
  `z_max=15` for several mass bins, meaning the true M_res-crossing point
  for those halos lies beyond `z_max` and was invisible to the test --
  those halos were handed the maximum possible integration time the box
  allowed rather than their real (even earlier) formation time. Re-running
  the same grid at `z_max=30` confirmed formation redshift really does
  push much higher with finer `m_res` (e.g. median z_form 11.0 -> 20.5 for
  the `M_halo=1e8` bin going from `m_res_fraction=1e-3` to `1e-5`) --
  exactly the CDM small-scale-structure behaviour one would expect (no
  cutoff in the linear power spectrum, so ever-finer resolution keeps
  revealing earlier-forming, lower-mass structure). **However, `M_star`
  itself barely changed between the `z_max=15` and `z_max=30` runs**
  (e.g. `2.51e4 -> 2.71e4 Msun` for that same bin/resolution) despite
  formation redshift moving substantially. This rules out "the model
  keeps crediting more integration time at ever-earlier cosmic epochs" as
  the mechanism for `M_star`'s non-convergence -- extending the visible
  time window barely mattered once it was extended. The `M_star` growth
  with finer `m_res` is therefore better explained as **resolving more
  sub-structure/mergers within the already-accessible time window
  (`z <~ 15`)** rather than extending the assembly history backward in
  time -- i.e. genuinely closer to the user's original "more small-scale
  structure gets resolved as you refine, the same way an N-body
  convergence study would show" intuition, just confirmed by ruling out
  the competing "artificial extra formation time" explanation rather than
  assumed. The `z_max=15` pinning was real but turned out to be a
  secondary artifact, not the actual cause -- a useful reminder that a
  confirmed artifact isn't automatically *the* explanation for the
  symptom it was found while investigating.

## What this does NOT affect

The reionization on/off and gamma/omega sensitivity tests (see the
adjacent `stellar_to_halo_mass_relation.py` and
`reionization_gamma_omega_sensitivity.py` scripts) were run at the
default `m_res_fraction=1e-3`, so their *qualitative* conclusions (UVB
suppression drives the low-mass plateau/kink; BH growth has negligible
effect on `M_star`; `omega` matters more than `gamma` for the transition
shape) should still hold, since all of those comparisons held `m_res`
fixed across the compared runs -- but their *absolute* low-mass `M_star`
values should not be treated as final for the same reason as above.

## Reproducing

See `scripts/resolution_timestep_sensitivity.py` (one-parameter-at-a-time,
superseded by the joint test below but kept for the record of how the
confound was found) and `scripts/joint_resolution_convergence.py` (the
joint `dz x m_res` grid that produced the conclusion above).

## Code changes made as a result

`ashvini/pymctrees_adapter.py`'s `build_forest_live()` now calls
`_warn_if_mres_unconverged()`, which emits a `UserWarning` whenever
`m_res/mass_bin >= 1e-3` (the empirically-known-bad default) is used. This
is a "known insufficient" floor, not a "verified converged" ceiling --
passing a smaller `m_res` silences the warning but is not by itself proof
of convergence for your specific mass range; rerun a convergence check
like `joint_resolution_convergence.py` for any new regime before trusting
absolute low-mass results.

## Follow-up work

1. Extend the `m_res_fraction` grid to `1e-6`/`1e-7` (at `z_max=30`, now
   that the ceiling artifact is understood and avoided) to find where, if
   anywhere, `M_star` actually stops changing.
2. ~~Investigate *why* ever-finer `m_res` keeps adding stellar mass~~ --
   partially answered above: it is resolving more sub-structure/mergers
   within the already-accessible time window, not extending the
   integration window backward. Still open: whether that resolved
   sub-structure's contribution to `M_star` itself converges (point 1
   above) or whether a physical minimum halo mass/virial-temperature
   floor for star formation is needed independent of the numerical `m_res`
   choice.
3. Consider whether `run_params.yaml`'s default `m_res_fraction=1e-3`
   should be changed, given it is now known to be badly unconverged for
   any low-mass science case, not just this SHMR calculation.
4. `z_max=15` (the value used elsewhere in Ashvini's own notebooks/demos
   for similar low-mass work) is itself a trap at fine `m_res` -- any
   future low-mass, fine-resolution study should default to a
   substantially higher `z_max` (`>=30`) and verify formation redshift
   isn't pinned at the ceiling, the same check that caught this here.
