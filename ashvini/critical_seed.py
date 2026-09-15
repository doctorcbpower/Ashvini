"""
Critical black-hole seed mass M_seed,crit(M_halo): the seed mass at which
strict Eddington-limited growth exactly reaches a specified degree of
overmassiveness by the anchor redshift (2026 "Differential Growth" paper,
Section "The critical seed mass"):

    M_BH(z_anchor; M_seed,crit) = f_BH * M_star(z_anchor),   f_BH = 0.5 fiducial

M_seed < M_seed,crit means Eddington-limited growth alone cannot reach the
target (a super-Eddington episode, earlier formation, or a lower effective
radiative efficiency would be needed); M_seed >= M_seed,crit means it can.

This was previously only a standalone analysis-script function
(seed_halo_boundary.py's critical_seed, in the paper's now-lost ephemeral
session sandbox -- see docs/2026_paper_session_code_catalogue.md). Rebuilt
here, on top of run_forest()'s "hobbs_slimdisk" growth model and its
a_star/r_crit/epsilon_f/eta_acc overrides, as a general, reusable utility
(not paper-specific glue) -- any caller with a halo assembly history can
ask "what seed mass reaches this target".

Seeding mechanism: bisection needs to seed every halo at its own,
independently-scanned trial mass, at the very first redshift step each
halo actually exists (its tree's own formation time) -- not at a fixed
halo mass threshold tied to any one halo's own resolution limit. This is
done by (temporarily) reconfiguring the existing
black_holes.seeding.halo_mass_threshold channel with M_halo_min set below
any halo's own mass resolution (so every halo seeds at its first nonzero
step) and M_seed set to a per-halo NumPy array of trial seed masses --
seeding_mask_and_mass's own np.where(mask, hm.M_seed, seed_mass)
broadcasts an (N,)-shaped M_seed against the (N,)-shaped mask exactly as
if it were the scalar the dataclass field was designed for; no core
package change was needed for this. pop3/direct_collapse are disabled for
the duration so the scanned seed mass is the only channel that can fire.
The original seeding config is restored afterward unconditionally.
"""

import numpy as np

from . import main
from .run_params import PARAMS

STRICT_EDDINGTON_R_CRIT = 1e12
# r_crit=1.0 does NOT give strict Eddington-limited growth, despite an
# earlier version of this module's own docstring claiming it does. Found
# and fixed 2026-09-15 in the sibling paper_reservoir.py module (see
# STRICT_EDDINGTON_CHI_CRIT there for the full derivation and a direct
# empirical check against bh_growth_step_slimdisk); the same bug was
# independently present here, unfixed until now, since critical_seed()
# reimplements the identical bisection pattern against main.run_forest()
# instead of paper_reservoir.run_reservoir_paper(). Summary: the
# closed-form growth update's three regimes have boundaries
# M1=A_bh/(kappa_edd*r_crit) and M2=A_bh/kappa_edd=M1*r_crit; at r_crit=1,
# M1=M2 exactly, so the genuine Eddington-limited (exponential,
# kappa_edd-dependent) regime between them has zero width, and growth
# equals the raw, uncapped supply rate A_bh at every black hole mass,
# completely independent of kappa_edd/epsilon. r_crit -> infinity is the
# correct limit for "no super-Eddington cap ever engages" (verified
# directly against bh_growth_step_slimdisk in the sibling module and its
# own test suite, which uses r_crit=1e12 for exactly this limit) -- not
# r_crit=1, which instead removes the Eddington cap entirely.


def _seed_all_and_run(halo_mass, halo_mass_rate, redshift, seed_mass, m_halo_min, **run_forest_kwargs):
    seeding = PARAMS.bh.seeding
    saved = (
        seeding.pop3.enabled,
        seeding.direct_collapse.enabled,
        seeding.halo_mass_threshold.enabled,
        seeding.halo_mass_threshold.M_halo_min,
        seeding.halo_mass_threshold.M_seed,
    )
    try:
        seeding.pop3.enabled = False
        seeding.direct_collapse.enabled = False
        seeding.halo_mass_threshold.enabled = True
        seeding.halo_mass_threshold.M_halo_min = m_halo_min
        seeding.halo_mass_threshold.M_seed = np.asarray(seed_mass, dtype=float)
        return main.run_forest(halo_mass, halo_mass_rate, redshift, **run_forest_kwargs)
    finally:
        (
            seeding.pop3.enabled,
            seeding.direct_collapse.enabled,
            seeding.halo_mass_threshold.enabled,
            seeding.halo_mass_threshold.M_halo_min,
            seeding.halo_mass_threshold.M_seed,
        ) = saved


def critical_seed(
    halo_mass,
    halo_mass_rate,
    redshift,
    f_bh=0.5,
    seed_mass_lo=1.0,
    seed_mass_hi=1.0e8,
    n_iter=50,
    m_halo_min=1.0,
    r_crit=STRICT_EDDINGTON_R_CRIT,
    **overrides,
):
    """
    Vectorised bisection for M_seed,crit, one value PER HALO (row) of
    halo_mass/halo_mass_rate (shape (N, n)) -- e.g. one per tree in a
    merger-tree ensemble. `redshift` (shape (n,)) is shared across all
    haloes and must be chronological ascending-in-cosmic-time (Ashvini's
    own convention, see pymctrees_adapter.build_forest_for_bin) -- the
    anchor redshift is redshift[-1].

    r_crit defaults to STRICT_EDDINGTON_R_CRIT (strict Eddington-limited
    growth: no super-Eddington cap ever engages). See that constant's own
    comment for why r_crit=1.0 -- this parameter's former default --
    does NOT give strict Eddington-limited growth despite an earlier
    version of this docstring claiming it does: it instead removes the
    Eddington cap entirely, making growth equal the raw, uncapped nuclear
    supply rate at every black hole mass. Passing a smaller, finite
    r_crit computes a *different* quantity (a seed-mass boundary under a
    permitted super-Eddington episode) -- valid, but not "M_seed,crit" in
    the paper's sense.

    Growth is monotonically non-decreasing in seed mass at
    r_crit=STRICT_EDDINGTON_R_CRIT (a genuine two-regime --
    Eddington-exponential then supply-limited -- closed form in y0, see
    bh_growth_step_slimdisk's docstring), so a straight bisection (not a
    general root-find) is well-posed here, provided M_star(z_anchor) does
    not vary enough with M_seed via AGN feedback to overwhelm that
    monotonicity -- true for any realistic epsilon_f (the coupling is a
    small perturbation on the gas budget, not a dominant term -- see
    MODELS.md's AGN feedback section). This module has not been
    re-exercised against real data since the r_crit fix (unlike the
    sibling paper_reservoir.critical_seed_paper, which was stress-tested
    extensively) -- callers relying on it should re-verify their own
    results were computed after this fix, not before.

    Any additional **overrides (a_star, epsilon_f, eta_acc) are forwarded
    to run_forest() unchanged, alongside growth_model="hobbs_slimdisk"
    (always set here -- this boundary is only defined for that model).

    Returns
    -------
    M_seed_crit : (N,) array, Msun. NaN where seed_mass_hi does not reach
        the target (never_reaches_target) or where seed_mass_lo already
        exceeds it (already_above_at_lo) -- see those two returned masks.
    never_reaches_target : (N,) bool array -- even the largest seed mass
        tried (seed_mass_hi) does not reach f_bh*M_star(z_anchor) under
        strict Eddington growth. This is the light-seed-problem boundary
        itself: widen seed_mass_hi if this fires for seed masses you
        actually care about resolving, rather than trusting the NaN.
    already_above_at_lo : (N,) bool array -- even the smallest seed mass
        tried (seed_mass_lo) already exceeds the target (M_seed,crit is
        below your bracket). Widen seed_mass_lo (downward) if this fires.
    """
    halo_mass = np.atleast_2d(halo_mass)
    halo_mass_rate = np.atleast_2d(halo_mass_rate)
    N = halo_mass.shape[0]

    run_forest_kwargs = dict(growth_model="hobbs_slimdisk", r_crit=r_crit, **overrides)

    def excess(seed_mass):
        result = _seed_all_and_run(
            halo_mass, halo_mass_rate, redshift, seed_mass, m_halo_min, **run_forest_kwargs
        )
        return result["bh_mass"][:, -1] - f_bh * result["stars_mass"][:, -1]

    lo = np.full(N, float(seed_mass_lo))
    hi = np.full(N, float(seed_mass_hi))

    excess_lo = excess(lo)
    excess_hi = excess(hi)
    already_above_at_lo = excess_lo > 0.0
    never_reaches_target = excess_hi < 0.0

    for _ in range(n_iter):
        mid = np.sqrt(lo * hi)  # bisect in log-space
        excess_mid = excess(mid)
        go_higher = excess_mid < 0.0
        lo = np.where(go_higher, mid, lo)
        hi = np.where(go_higher, hi, mid)

    M_seed_crit = np.sqrt(lo * hi)
    unresolved = never_reaches_target | already_above_at_lo
    M_seed_crit = np.where(unresolved, np.nan, M_seed_crit)
    return M_seed_crit, never_reaches_target, already_above_at_lo
