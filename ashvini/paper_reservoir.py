"""
Dedicated reservoir integrator for the 2026 "Differential Growth" paper's
own bespoke halo-baryon-accretion + gas/star ODEs (the paper's Eqs.
halo_accretion, inflow, gas, stars) -- NOT the same computation as
main.run_forest(), which implements Ashvini's own separate,
general-purpose gas_inflow_rate/star_formation_rate, calibrated for a
different set of production runs (the Menon & Power 2024 reproduction,
SFR/dust/metals modelling, etc.) that this paper's model was never
validated against and was not intended to share.

Only the pieces that WERE actually ported into Ashvini and validated
against this paper are reused here: black_holes_growth_slimdisk's
closed-form BH growth + graded super-Eddington cap, spin's Novikov-Thorne
epsilon(a*), and the King (2003) AGN wind term. The halo-accretion/
gas-inflow/star-formation equations below are new, implementing the
paper's own equations directly -- see docs/2026_paper_session_code_catalogue.md
for why the original standalone prototype that first implemented them
(ashvini_halo_model.py) no longer exists anywhere to import instead.

Section 2 of the paper (2026-09 revision) settles the two pieces that were
previously undefined:

    zeta(M_halo, z) = zeta_UV(M_halo, z) * zeta_ch(M_halo): the product of
        a UV-suppression window (uv_suppression_mass_only below, the
        literature-standard Okamoto et al. 2008 mass/z-only step function
        -- NOT reionization.uv_suppression's own full form, whose extra
        ratio_term correction needs a genuinely differentiable dM_halo/dt
        that a real merger tree's raw mass history does not have at a
        resolved merger; ablation-tested directly to change M_seed,crit
        by at most +0.04 dex across this paper's full halo-mass grid, so
        dropped in favour of the simpler, non-fragile literature form)
        and a cold/hot-mode transition (cold_hot_mode_suppression below,
        same step-function form reused with different arguments, per
        Kravtsov & Manwadkar 2022's implementation choice) -- see Eq.
        coldhot/stepfunc. Applied to the *full* accreted mass (smooth +
        merger, undifferentiated).

    t_star = t_ff/epsilon_sf, where t_ff is the *same* nuclear free-fall
        time (black_holes_growth_slimdisk.time_freefall_enclosed,
        M_enc = M_BH+M_gas+M_star, R_nuc=100 pc) that sets the black
        hole's own nuclear supply -- not a separate galactic-scale
        dynamical time. An earlier attempt using Ashvini's own
        time_freefall(z)=0.141*Hubble_time(z) gave M_BH ~50x M_star at the
        fiducial halo mass; a second attempt using the paper's own
        orbital-timescale approximation (t_orb ~ 8.8 Myr *
        (r/r_200/1e-3)*(1+z)^-1.5) got within an order of magnitude but
        not exactly reproducing the target -- both superseded by this
        shared-clock picture, which keeps stars and the black hole
        competing for literally the same reservoir on the same timescale,
        differing only by efficiency (epsilon_sf vs eta_acc).

Any number produced by this module that is quoted in the paper should
still be cross-checked against Ashvini's own record of the original run,
if one turns up, rather than trusted blind -- see
docs/2026_paper_session_code_catalogue.md.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import erf

from . import black_holes_growth_slimdisk as bh_slimdisk
from . import reionization as reion
from .constants import G, Msun_g, pc_cm
from .reionization import s as _okamoto_step
from .spin import epsilon_from_spin
from .utils import time_at_z, z_at_time
from .paper_reservoir_params import PAPER_PARAMS

# Sourced from paper_reservoir_params.yaml (via PAPER_PARAMS) rather than
# hardcoded here, so that every fiducial value used across this module's
# many function-signature defaults stays in exactly one place -- a stale,
# independently-hardcoded eta_acc default drifted out of sync with the
# paper's own stated fiducial for a period before being caught (see
# docs/2026_paper_session_code_catalogue.md). See paper_reservoir_params.py
# for the full parameter set and its documentation.
F_B_FIDUCIAL = PAPER_PARAMS.f_b  # cosmic baryon fraction
M_HOT_FIDUCIAL = PAPER_PARAMS.M_hot  # Msun; cold/hot-mode transition mass (Eq. coldhot)
PHI_COLDHOT_FIDUCIAL = PAPER_PARAMS.phi_coldhot  # dimensionless step sharpness (Eq. coldhot)
LAMBDA_MEDIAN_FIDUCIAL = PAPER_PARAMS.lam_median  # Bullock et al. (2001) median halo spin parameter
DELTA_VIR_FIDUCIAL = PAPER_PARAMS.delta_vir  # overdensity (relative to rho_crit) defining R_vir


def virial_radius_velocity(M_halo_msun, z, omega_m, omega_lambda, H0_km_s_mpc, delta_vir=DELTA_VIR_FIDUCIAL):
    """
    R_vir (cm), V_vir (cm/s) for a halo of mass M_halo_msun at redshift z,
    via M_halo = (4/3)*pi*delta_vir*rho_crit(z)*R_vir^3 (delta_vir=200,
    relative to the critical density) and V_vir = sqrt(G*M_halo/R_vir).
    """
    H0_cgs = H0_km_s_mpc * 1.0e5 / (1.0e6 * pc_cm)  # km/s/Mpc -> 1/s
    H_z = H0_cgs * np.sqrt(omega_m * (1.0 + z) ** 3 + omega_lambda)
    rho_crit = 3.0 * H_z**2 / (8.0 * np.pi * G)
    M_halo_g = np.asarray(M_halo_msun, dtype=float) * Msun_g
    R_vir_cm = (3.0 * M_halo_g / (4.0 * np.pi * delta_vir * rho_crit)) ** (1.0 / 3.0)
    V_vir_cms = np.sqrt(G * M_halo_g / R_vir_cm)
    return R_vir_cm, V_vir_cms


def low_spin_available_fraction(
    M_enc_msun, M_halo_msun, z, R_nuc_pc=100.0, lam_median=LAMBDA_MEDIAN_FIDUCIAL, sigma_lnj=0.5,
    omega_m=None, omega_lambda=None, H0_km_s_mpc=67.8,
):
    """
    f_avail(M_enc, M_halo, z): the fraction of a galaxy's gas whose
    specific angular momentum j is low enough to circularize within
    R_nuc, given the enclosed mass M_enc there -- the mechanism that
    replaces an ad hoc galaxy-to-nucleus transport timescale with an
    angular-momentum selection.

    A gas parcel with specific angular momentum j circularizes at
    R_circ = j^2/(G*M_enc) (Keplerian, v_circ^2=GM_enc/R_circ,
    j=v_circ*R_circ); only parcels with R_circ <= R_nuc can reach the
    nuclear region, i.e. j <= j_thresh = sqrt(G*M_enc*R_nuc). Assuming
    specific angular momentum across gas parcels is log-normal with
    median j_gas = lambda*sqrt(2)*V_vir*R_vir (the halo's own bulk
    specific angular momentum, Bullock et al. 2001 normalization,
    lambda=0.035 the fiducial median) and log-normal width sigma_lnj (a
    free parameter -- unlike lambda's halo-to-halo scatter, which is a
    different, better-constrained quantity, the parcel-to-parcel spread
    within one halo's gas is not; treated as a free parameter here, to be
    explored rather than assumed), the fraction below threshold is the
    log-normal CDF:

        f_avail = Phi( ln(j_thresh/j_gas) / sigma_lnj )

    Critically, j_thresh grows only as sqrt(M_enc), so f_avail shrinks as
    the nuclear reservoir grows -- this is what breaks the runaway
    feedback loop that resulted from feeding the *entire* galaxy-scale
    M_gas into the R_nuc=100 pc free-fall time directly (t_ff crashing as
    M_enc grew unboundedly): only ever a self-limiting fraction of the
    galaxy's gas is treated as nuclear-available at all.
    """
    omega_m = OMEGA_M_FIDUCIAL if omega_m is None else omega_m
    omega_lambda = OMEGA_LAMBDA_FIDUCIAL if omega_lambda is None else omega_lambda

    M_halo_msun = np.asarray(M_halo_msun, dtype=float)
    has_halo = M_halo_msun > 0
    M_halo_safe = np.where(has_halo, M_halo_msun, 1.0)  # dummy value, never used below (see has_halo guard)

    R_vir_cm, V_vir_cms = virial_radius_velocity(M_halo_safe, z, omega_m, omega_lambda, H0_km_s_mpc)
    j_gas = lam_median * np.sqrt(2.0) * V_vir_cms * R_vir_cm

    M_enc_g = np.maximum(np.asarray(M_enc_msun, dtype=float), 0.0) * Msun_g
    R_nuc_cm = R_nuc_pc * pc_cm
    j_thresh = np.sqrt(G * M_enc_g * R_nuc_cm)

    j_gas_safe = np.where(j_gas > 0, j_gas, 1.0)
    j_thresh_safe = np.where(j_thresh > 0, j_thresh, 1e-300)
    x = np.log(j_thresh_safe / j_gas_safe) / sigma_lnj
    f_avail = 0.5 * (1.0 + erf(x / np.sqrt(2.0)))
    # No host halo yet (e.g. a tree's pre-formation frozen prefix, where
    # pymctrees_adapter zeroes M_halo): there is no established potential
    # to set j_gas at all, so nothing should be treated as nuclear-available
    # yet, rather than falling back on whatever j_gas=nan/0 happens to
    # produce through the safety guards above.
    return np.where(has_halo, f_avail, 0.0)


def cold_hot_mode_suppression(M_halo_msun, M_hot=M_HOT_FIDUCIAL, phi=PHI_COLDHOT_FIDUCIAL):
    """
    zeta_ch(M_halo) = 1 - s(M_halo/M_hot, phi) -- Eq. coldhot: the soft
    step-function suppression of baryonic inflow once the halo grows
    massive enough to sustain a virial shock and transition from cold- to
    hot-mode accretion (Birnboim & Dekel 2003; Keres et al. 2005, 2009;
    Dekel & Birnboim 2006; Dekel et al. 2009), reusing the same step
    function s(x,y) already implemented for UV suppression
    (reionization.s) rather than a second, independent formula (cf.
    Kravtsov & Manwadkar 2022). M_hot=4e11 Msun lies within this paper's
    own 3e10-3e13 Msun halo-mass grid, so this term is expected to matter
    mainly for the more massive halos in the sample.
    """
    M_halo_msun = np.asarray(M_halo_msun, dtype=float)
    mu_hot = np.where(M_halo_msun > 0, M_halo_msun / M_hot, 1e-30)
    return 1.0 - _okamoto_step(mu_hot, phi)


def uv_suppression_mass_only(z_val, m_halo):
    """
    zeta_UV(M_halo, z) = s(mu_c(z), omega) -- the literature-standard
    Okamoto et al. (2008) UV-suppression step function, a function of
    M_halo and z only. Used here instead of reionization.uv_suppression's
    full form, which adds a correction (`ratio_term`) requiring a
    genuinely differentiable dM_halo/dt -- something a real merger tree's
    raw mass history does not have at a resolved merger (see
    grumpy_halo_growth_rate's docstring for why estimating one is
    unavoidably fragile there).

    Ablation-tested directly against the full ratio_term-corrected form
    for this paper's own halo-mass grid (3e10-3e13 Msun, N=200 trees per
    point): dropping ratio_term changes M_seed,crit by at most +0.04 dex
    (fiducial mass) and +0.00 dex (high mass) -- indistinguishable from
    noise. ratio_term was also the *only* place in this entire
    calculation where a rate estimate is divided into, rather than
    multiplied into (as Mdot_b = f_b*halo_growth_rate is): a rate
    multiplied by dt only ever gives a bounded mass increment, so
    grumpy_halo_growth_rate's own individual-tree rate spikes (checked
    directly: up to ~80x the characteristic rate for some trees, even at
    the fiducial dz=0.05) pose no numerical risk there. Dropping
    ratio_term therefore removes the one demonstrably fragile piece of
    this calculation in exchange for a physically inconsequential
    correction, rather than continuing to manage its numerical behaviour.
    """
    z_val = np.asarray(z_val, dtype=float)
    m_halo = np.asarray(m_halo, dtype=float)
    suppression = np.ones_like(z_val)
    mask = z_val <= 10
    mh_masked = m_halo[mask]
    resolved = mh_masked > 0
    mh_safe = np.where(resolved, mh_masked, 1.0)
    vals = _okamoto_step(reion.mu_c(z_val[mask], mh_safe), reion.omega)
    suppression[mask] = np.where(resolved, vals, 0.0)
    return suppression

# Menon & Power (2024) cosmology, matching the tree ensembles' own config
# (foraois/config/menon_power_2024.yml) -- see smooth_halo_growth_rate.
OMEGA_M_FIDUCIAL = PAPER_PARAMS.omega_m
OMEGA_LAMBDA_FIDUCIAL = PAPER_PARAMS.omega_lambda


def smooth_halo_growth_rate(M_halo_msun, z, omega_m=OMEGA_M_FIDUCIAL, omega_lambda=OMEGA_LAMBDA_FIDUCIAL):
    """
    dM_halo/dt (Msun/Gyr): the Fakhouri, Ma & Boylan-Kolchin (2010, MNRAS
    406, 2267) mean halo mass accretion rate fit (their Eq. 2),

        dM/dt [Msun/yr] = 46.1 * (M_halo/1e12)^1.1 * (1+1.11z)
                           * sqrt(Omega_m*(1+z)^3 + Omega_Lambda)

    NOT the paper's own Eq. halo_accretion (a much weaker, near-z-
    independent specific rate that, taken at face value, barely grows a
    halo at all between z=25 and z=5 -- checked directly: integrating it
    gives M_halo(z=25)/M_halo(z=5) ~ 0.94 for a fiducial 5e10 Msun halo,
    two orders of magnitude too flat for real hierarchical assembly, e.g.
    against this FMBK10 rate's own implied growth or basic H(z) scaling).
    ashvini.pymctrees_adapter.build_forest_for_bin's own module docstring,
    and docs/2026_paper_session_code_catalogue.md's description of the
    lost standalone prototype's ashvini_cosmology.py, both independently
    name FMBK10 (via Dave, Finlator & Oppenheimer 2012) as what the
    original calculation actually used -- Eq. halo_accretion in the
    current paper draft is very likely a transcription error from an
    earlier drafting pass, not what produced the paper's existing
    reported numbers; worth fixing in the paper text itself, separately
    from this module.
    """
    M_halo_msun = np.asarray(M_halo_msun, dtype=float)
    rate_per_yr = 46.1 * (M_halo_msun / 1.0e12) ** 1.1 * (1.0 + 1.11 * z) * np.sqrt(
        omega_m * (1.0 + z) ** 3 + omega_lambda
    )
    return rate_per_yr * 1.0e9  # Msun/yr -> Msun/Gyr


def halo_baryon_accretion_rate(M_halo_msun, z, f_b=F_B_FIDUCIAL, **fmbk10_kwargs):
    """Mdot_b (Msun/Gyr): the baryonic share (universal fraction f_b) of smooth_halo_growth_rate."""
    return f_b * smooth_halo_growth_rate(M_halo_msun, z, **fmbk10_kwargs)


def integrate_smooth_halo_trajectory(M_halo_seed_msun, z_seed, z_anchor, n_steps=400, **fmbk10_kwargs):
    """
    Forward-integrates the smooth halo trajectory M_halo(z) from
    (M_halo_seed_msun, z_seed) down to z_anchor via smooth_halo_growth_rate
    (RK4, uniform in cosmic time).

    Returns (M_halo_msun, redshift, cosmic_time), each shape (n_steps+1,),
    chronological (ascending cosmic time / descending redshift -- index 0
    = z_seed, index -1 = z_anchor), matching main.run_forest()'s own
    halo_mass/redshift convention.
    """
    t_seed = float(time_at_z(np.array([z_seed]))[0])
    t_anchor = float(time_at_z(np.array([z_anchor]))[0])
    t_grid = np.linspace(t_seed, t_anchor, n_steps + 1)

    def rate(t, m):
        z = float(z_at_time(np.array([t]))[0])
        return smooth_halo_growth_rate(max(m, 0.0), z, **fmbk10_kwargs)

    M = np.empty(n_steps + 1)
    M[0] = M_halo_seed_msun
    for j in range(n_steps):
        dt = t_grid[j + 1] - t_grid[j]
        t0, t1 = t_grid[j], t_grid[j] + dt
        k1 = rate(t0, M[j])
        k2 = rate(t0 + 0.5 * dt, M[j] + 0.5 * dt * k1)
        k3 = rate(t0 + 0.5 * dt, M[j] + 0.5 * dt * k2)
        k4 = rate(t1, M[j] + dt * k3)
        M[j + 1] = M[j] + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    redshift = z_at_time(t_grid)
    return M, redshift, t_grid


def shoot_smooth_halo_trajectory(
    M_halo_target_msun, z_seed=25.0, z_anchor=5.0, n_steps=400,
    seed_lo=1.0, seed_hi=None, tol=1e-3, max_iter=60, **fmbk10_kwargs,
):
    """
    Bisects on the z_seed initial halo mass so forward integration of
    smooth_halo_growth_rate lands within `tol` (fractional) of
    M_halo_target_msun at z_anchor -- the paper's own "shooting" method
    for constructing a smooth assembly history anchored to a specified
    z=5 halo mass (Section "Smooth halo assembly histories").

    Growth is monotonically increasing in the initial mass (a larger
    seed mass can only grow into an equal-or-larger mass under this
    positive-definite growth rate), so M_halo_target_msun itself is a
    safe upper bracket by default.

    Returns the same (M_halo_msun, redshift, cosmic_time) as
    integrate_smooth_halo_trajectory, for the converged initial mass --
    plus a bool `converged` (False if max_iter was exhausted without
    reaching `tol`; the last trial is still returned).
    """
    if seed_hi is None:
        seed_hi = M_halo_target_msun

    lo, hi = seed_lo, seed_hi
    converged = False
    M, redshift, t_grid = None, None, None
    for _ in range(max_iter):
        mid = np.sqrt(lo * hi)
        M, redshift, t_grid = integrate_smooth_halo_trajectory(mid, z_seed, z_anchor, n_steps, **fmbk10_kwargs)
        final = M[-1]
        if abs(final / M_halo_target_msun - 1.0) < tol:
            converged = True
            break
        if final < M_halo_target_msun:
            lo = mid
        else:
            hi = mid

    return M, redshift, t_grid, converged


def interpolate_tree_onto_grid(tree_mass_msun, tree_redshift, target_redshift):
    """
    Interpolates a tree's main-progenitor mass history (`tree_mass_msun`,
    shape (N, n_tree), chronological -- ascending cosmic time -- against
    `tree_redshift`, shape (n_tree,)) onto `target_redshift` (e.g. the
    uniform grid a smooth trajectory's own reservoir integration uses),
    in cosmic time (not redshift -- the natural variable for this
    integration, see the paper's own "the exponential black-hole growth
    law and the galaxy reservoir evolution are evaluated on the same
    physical clock"). Linear in cosmic time; masses below the tree's own
    resolution limit are already 0 (see pymctrees_adapter's own
    pre-formation masking) and interpolate smoothly into the first
    resolved value, not discontinuously -- no separate masking needed here.

    Returns an (N, len(target_redshift)) array. `target_redshift`'s own
    endpoints must lie within [tree_redshift.min(), tree_redshift.max()]
    (raises via np.interp's own extrapolation-as-clipping otherwise --
    callers should choose target_redshift bracketed by the tree's own
    z_seed/z_anchor).
    """
    tree_mass_msun = np.atleast_2d(tree_mass_msun)
    t_tree = time_at_z(tree_redshift)
    t_target = time_at_z(target_redshift)
    return np.array([np.interp(t_target, t_tree, row) for row in tree_mass_msun])


def grumpy_halo_growth_rate(tree_mass_msun, tree_redshift, target_redshift):
    """
    Per-halo smoothed dM_halo/dt (Msun/Gyr), following the mass-accretion-
    history smoothing recipe of Kravtsov & Manwadkar (2022, GRUMPY): fit
    each tree's own raw main-progenitor mass history (mergers included,
    not split out) with a cubic spline of log M_halo vs. log t; take that
    spline's analytic derivative to get dlnM/dlnt at the tree's own native
    time nodes; clip negative values to zero (their monotonicity
    condition -- Mdot should read zero during a real mass *decrease*,
    e.g. post-merger relaxation or stripping, never negative); re-spline
    this clipped dlnM/dlnt as a function of cosmic time t; then evaluate
    Mdot(t) = (M_halo(t)/t) * dlnM/dlnt(t) wherever needed.

    This replaces two earlier, unsuccessful attempts at getting a
    well-behaved rate out of a real (bursty) merger-tree trajectory: a
    single-step finite difference (noise-amplifying at a merger jump) and
    a fixed-physical-time trailing window (still let M_seed,crit swing
    76x between two reasonable window choices). GRUMPY's own reasoning
    for reionization/cold-hot suppression (reionization.uv_suppression's
    ratio_term, cold_hot_mode_suppression) requires a genuinely
    differentiable dM_halo/dt -- not "the accretion rate that feeds the
    reservoir" in a loose sense, but literally the time-derivative that
    falls out of differentiating the equilibrium relation
    M_gas = f_bar(M_halo, z) * M_halo. A tree's raw main-progenitor mass
    has real jump discontinuities at every resolved merger, where
    dM_halo/dt does not exist in the classical sense at all -- no
    smoothing scheme applied *after* the fact can manufacture a
    derivative a genuine discontinuity doesn't have. GRUMPY's own use
    case sidesteps this by fitting the (fixed, simulation-snapshot-
    cadence) sampled mass history with a smooth interpolant first and
    differentiating *that*, which is exactly what this function does.

    Because this is still exact interpolation through the tree's own raw
    (jumpy) mass nodes, it inherits the same fundamental limit as any
    such recipe: as the tree's own step size dz shrinks, a genuine merger
    event's mass jump spans an ever-smaller time interval, and the local
    interpolated slope there grows without bound in exactly that limit --
    checked directly (see docs/2026_paper_session_code_catalogue.md for
    the dz=0.0125 case where one tree's merger jump right at z_anchor
    produced a >100x outlier in M_seed,crit). This is not a numerical bug
    to patch away; it reflects that "smooth accretion rate" is not a
    resolution-independent quantity for a stochastic assembly history --
    exactly GRUMPY's own point in fixing their node spacing to their
    simulation's snapshot cadence rather than refining it indefinitely.
    `dz` here is therefore a genuine, disclosed resolution choice (like
    M_res), not a numerical-error control parameter to converge to zero;
    see the paper's own Section 3 for the adopted fiducial value.

    Parameters
    ----------
    tree_mass_msun : (N, n_tree) array, Msun -- a tree's raw main-
        progenitor mass history (pymctrees_adapter.build_forest_for_bin's
        `halo_masses`, NOT split into smooth/merger channels), 0 during
        any pre-formation prefix.
    tree_redshift : (n_tree,) array -- the tree's own native z_steps grid.
    target_redshift : (n,) array -- where to evaluate Mdot (e.g. the
        reservoir integration's own uniform-in-cosmic-time grid).

    Returns
    -------
    rate : (N, n) array, Msun/Gyr -- 0 wherever target_redshift falls
        outside a given tree's resolved (post-formation) range.
    """
    tree_mass_msun = np.atleast_2d(tree_mass_msun)
    N = tree_mass_msun.shape[0]
    t_native = time_at_z(tree_redshift)
    t_target = time_at_z(target_redshift)
    n_out = len(t_target)
    rate = np.zeros((N, n_out))

    for i in range(N):
        m = tree_mass_msun[i]
        alive = m > 0
        if alive.sum() < 4:
            continue  # too few resolved points to spline meaningfully
        t_alive = t_native[alive]
        m_alive = m[alive]
        order = np.argsort(t_alive)
        t_alive = t_alive[order]
        m_alive = m_alive[order]

        cs_logm = CubicSpline(np.log(t_alive), np.log(m_alive))
        dlnM_dlnt_nodes = np.clip(cs_logm(np.log(t_alive), 1), 0.0, None)
        cs_rate = CubicSpline(t_alive, dlnM_dlnt_nodes)

        in_range = (t_target >= t_alive.min()) & (t_target <= t_alive.max())
        tt = t_target[in_range]
        m_t = np.exp(cs_logm(np.log(tt)))
        dlnM_dlnt_t = np.clip(cs_rate(tt), 0.0, None)
        rate[i, in_range] = (m_t / tt) * dlnM_dlnt_t

    return rate


def run_reservoir_paper(
    halo_mass, redshift, seed_mass, f_b=F_B_FIDUCIAL,
    epsilon_sf=PAPER_PARAMS.epsilon_sf, a_star=None, epsilon=None,
    chi_crit=PAPER_PARAMS.chi_crit, epsilon_f=PAPER_PARAMS.epsilon_f, eta_acc=PAPER_PARAMS.eta_acc,
    R_nuc_pc=PAPER_PARAMS.R_nuc_pc, compaction_boost=PAPER_PARAMS.compaction_boost,
    M_hot=M_HOT_FIDUCIAL, phi_coldhot=PHI_COLDHOT_FIDUCIAL,
    sigma_lnj=PAPER_PARAMS.sigma_lnj, lam_median=LAMBDA_MEDIAN_FIDUCIAL,
    gas_mass0=PAPER_PARAMS.gas_mass0, stars_mass0=PAPER_PARAMS.stars_mass0,
    halo_growth_rate=None,
):
    """
    Integrates the paper's own gas/star/BH reservoir ODEs (Eqs. gas,
    stars, bh_growth_general) forward along a GIVEN halo trajectory
    (`halo_mass`, shape (N, n) or (n,); `redshift`, shape (n,), shared,
    chronological -- ascending cosmic time / descending z). `halo_mass`
    may be a smooth trajectory (shoot_smooth_halo_trajectory) or a real
    merger-tree main-progenitor branch already interpolated onto this
    function's own redshift grid (interpolate_tree_onto_grid) -- "No
    additional intermediate progenitor-mass constraint is imposed"
    (Section "Merger-tree assembly histories"): the raw tree mass is
    substituted directly, unmasked.

    seed_mass : (N,) array or scalar -- black hole mass at the first
    (z_seed) step. Growth is always "hobbs_slimdisk"'s physics (Hobbs
    enclosed-mass free-fall nuclear supply + graded chi_crit
    super-Eddington cap + King 2003 feedback); chi_crit -> infinity (in
    practice, STRICT_EDDINGTON_CHI_CRIT below) reproduces "strict
    Eddington-limited growth" (Eq. strict_eddington) exactly -- see
    black_holes_growth_slimdisk.bh_growth_step_slimdisk's docstring (its
    own r_crit parameter is this same quantity, named chi_crit in the
    paper to avoid a symbol collision). chi_crit=1 does NOT do this: it
    makes the boundaries M1=A_bh/(kappa_edd*chi_crit) and
    M2=A_bh/kappa_edd coincide exactly, collapsing the genuine
    Eddington-limited regime to zero width and leaving growth equal to
    the raw, uncapped supply rate A_bh at every mass -- checked directly,
    this was a real bug in an earlier version of critical_seed_paper's
    own default (see STRICT_EDDINGTON_CHI_CRIT's comment). The fiducial
    chi_crit=10 is the default here (Section "Integration of the
    reservoir equations"), a genuine graded super-Eddington cap;
    critical-seed boundary calculations (critical_seed_paper below) pass
    chi_crit=STRICT_EDDINGTON_CHI_CRIT instead.

    a_star/epsilon: at most one may be given. epsilon=0.1 (the default
    here if neither is given) is the paper's own fiducial radiative
    efficiency, used as a free parameter in its own right -- independent
    of the a_star-derived value used specifically for the Novikov-Thorne
    spin/angular-momentum sensitivity scan (item 1a of the 2026 plots
    list).

    M_hot/phi_coldhot: the cold/hot-mode transition parameters entering
    zeta_ch (cold_hot_mode_suppression above), applied alongside UV
    suppression to the nuclear inflow below.

    sigma_lnj/lam_median: the low-spin-available-fraction parameters
    (low_spin_available_fraction above) that set what fraction of the
    galaxy-scale gas reservoir M_gas is actually available to the nuclear
    region (M_gas_nuc = f_avail*M_gas) for both the black hole's own
    supply and nuclear star formation -- both draw on this same,
    self-limiting nuclear subset and the same t_ff, not the full,
    unboundedly-growing galaxy reservoir directly (see the module
    docstring for why the latter produced a runaway).

    halo_growth_rate : optional (N, n) array, Msun/Gyr, already on this
    function's own step grid. If given, used directly as Mdot_halo for
    Mdot_b = f_b*Mdot_halo (Eq. inflow) instead of a rate derived from
    `halo_mass` itself -- zeta_UV (uv_suppression_mass_only) and zeta_ch
    (cold_hot_mode_suppression) both depend on M_halo and z only, not on
    this rate at all. When `halo_mass` is a real tree's
    main-progenitor mass, this should be grumpy_halo_growth_rate's output
    (computed from the tree's own *raw* native-grid mass history before
    interpolate_tree_onto_grid, then evaluated on this function's target
    redshift grid) -- see that function's docstring for why a rate has to
    be estimated this way at all. If omitted, falls back to
    smooth_halo_growth_rate(halo_mass, redshift) evaluated directly on
    this function's own grid, appropriate only when `halo_mass` is smooth
    by construction (e.g. shoot_smooth_halo_trajectory's output), since
    that closed-form rate has no notion of a real tree's merger jumps.

    Returns a dict: bh_mass, gas_mass, stars_mass, halo_mass, cosmic_time,
    redshift (each (N, n) except cosmic_time/redshift, shape (n,)), plus
    mdot_in/mdot_out (each (N, n)) for the outflow/inflow crossing-redshift
    diagnostic (item 4 of the 2026 plots list).
    """
    if a_star is not None and epsilon is not None:
        raise ValueError("run_reservoir_paper: pass at most one of a_star, epsilon.")

    halo_mass = np.atleast_2d(np.asarray(halo_mass, dtype=float))
    N, n = halo_mass.shape
    redshift = np.asarray(redshift, dtype=float)
    cosmic_time = time_at_z(redshift)

    if epsilon is not None:
        kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(epsilon)
    elif a_star is not None:
        kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(epsilon_from_spin(a_star))
    else:
        kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(PAPER_PARAMS.epsilon)  # paper's own fiducial epsilon

    gas_mass = np.zeros((N, n))
    stars_mass = np.zeros((N, n))
    bh_mass = np.zeros((N, n))
    mdot_in = np.zeros((N, n))
    mdot_out = np.zeros((N, n))
    A_bh_hist = np.zeros((N, n))  # nuclear BH accretion supply (pre-cap), Msun/Gyr -- diagnostic only

    gas_mass[:, 0] = gas_mass0
    stars_mass[:, 0] = stars_mass0
    bh_mass[:, 0] = seed_mass

    # halo_growth_rate feeds Mdot_b = f_b*Mdot_halo only (zeta_UV/zeta_ch
    # are both M_halo/z-only, see uv_suppression_mass_only). If the
    # caller didn't supply one, fall back to the closed-form
    # smooth_halo_growth_rate(halo_mass, z) -- only correct when
    # halo_mass is already smooth by construction (no merger jumps for
    # it to mishandle). A caller running on a real tree's mass history
    # must pass grumpy_halo_growth_rate's output explicitly; see this
    # function's own docstring.
    if halo_growth_rate is None:
        halo_growth_rate = smooth_halo_growth_rate(halo_mass, redshift)
    else:
        halo_growth_rate = np.atleast_2d(np.asarray(halo_growth_rate, dtype=float))

    for j in range(1, n):
        dt = cosmic_time[j] - cosmic_time[j - 1]
        z_mid = 0.5 * (redshift[j - 1] + redshift[j])
        hm_prev = halo_mass[:, j - 1]
        gm_prev = gas_mass[:, j - 1]
        sm_prev = stars_mass[:, j - 1]
        bh_prev = bh_mass[:, j - 1]

        # --- Eq. halo_accretion / Eq. inflow / Eq. coldhot ---
        halo_growth_rate_mid = halo_growth_rate[:, j - 1]
        mdot_b = f_b * halo_growth_rate_mid
        zeta_uv = uv_suppression_mass_only(np.full(N, z_mid), hm_prev)
        zeta_ch = cold_hot_mode_suppression(hm_prev, M_hot=M_hot, phi=phi_coldhot)
        A_acc = zeta_uv * zeta_ch * mdot_b  # Mdot_in, Msun/Gyr
        mdot_in[:, j] = A_acc

        # --- low-spin selection: only a fraction of the galaxy-scale gas
        # reservoir is available to the nuclear region at all (see module
        # docstring/low_spin_available_fraction). The angular-momentum
        # threshold is set by the *full* galaxy-scale enclosed mass
        # (M_BH+M_gas+M_star, the total baryonic reservoir the halo has
        # built up so far), not by the already-nuclear subset this
        # fraction itself produces: feeding the nuclear subset back into
        # its own threshold is an unstable fixed point (a small
        # perturbation downward shrinks the threshold, which shrinks the
        # available fraction further, collapsing to zero within a few
        # steps -- checked directly). The full reservoir grows
        # monotonically and gives a stable, self-limiting threshold instead.
        M_enc_full_prev = bh_prev + gm_prev + sm_prev
        f_avail = low_spin_available_fraction(
            M_enc_full_prev, hm_prev, z_mid, R_nuc_pc=R_nuc_pc, lam_median=lam_median, sigma_lnj=sigma_lnj
        )
        gas_nuc_prev = f_avail * gm_prev

        # --- shared nuclear clock: Eq. tff, M_enc = M_BH+M_gas_nuc+M_star ---
        M_enc_prev = bh_prev + gas_nuc_prev + sm_prev
        t_ff = bh_slimdisk.time_freefall_enclosed(M_enc_prev, R_nuc_pc=R_nuc_pc)

        # --- Eq. bh_growth_general / accretion_supply (Hobbs nuclear supply, same t_ff) ---
        A_bh = eta_acc * (gas_nuc_prev / t_ff) * compaction_boost**2
        A_bh_hist[:, j] = A_bh
        bh_new = bh_slimdisk.bh_growth_step_slimdisk(bh_prev, A_bh, kappa_edd, chi_crit, dt)
        bh_mass[:, j] = np.maximum(bh_new, 0.0)
        bh_growth_rate = (bh_mass[:, j] - bh_prev) / dt

        # --- Eq. feedback (King 2003 AGN wind, instantaneous) ---
        wind_rate = bh_slimdisk.king_agn_wind_rate(
            bh_growth_rate, halo_mass=hm_prev, redshift=z_mid, epsilon_f=epsilon_f
        )
        mdot_out[:, j] = wind_rate

        # --- Eq. stars: dStar/dt = Gas_nuc/t_star, t_star = t_ff/epsilon_sf (same nuclear reservoir/clock as the BH) ---
        sfr = gas_nuc_prev / (t_ff / epsilon_sf)

        # --- Eq. gas: dGas/dt = Mdot_in - Mdot_star - Mdot_BH - Mdot_out (full galaxy reservoir depleted, not just the nuclear subset) ---
        gas_new = gm_prev + dt * (A_acc - sfr - bh_growth_rate - wind_rate)
        gas_mass[:, j] = np.maximum(gas_new, 0.0)
        stars_mass[:, j] = sm_prev + dt * sfr

    return dict(
        bh_mass=bh_mass, gas_mass=gas_mass, stars_mass=stars_mass,
        halo_mass=halo_mass, cosmic_time=cosmic_time, redshift=redshift,
        mdot_in=mdot_in, mdot_out=mdot_out,
        A_bh=A_bh_hist, kappa_edd=float(np.asarray(kappa_edd).item()) if np.asarray(kappa_edd).size == 1 else kappa_edd,
    )


STRICT_EDDINGTON_CHI_CRIT = PAPER_PARAMS.strict_eddington_chi_crit
# bh_growth_step_slimdisk's own docstring and test suite
# (test_slimdisk_reduces_to_pznk11_two_regime_model_as_rcrit_to_infinity)
# establish that "strict Eddington-limited growth" -- dM/dt =
# min(A_bh, kappa_edd*M), no super-Eddington throttling ever engaging --
# is recovered as chi_crit -> infinity, NOT at chi_crit=1. At chi_crit=1,
# M1=A_bh/(kappa_edd*chi_crit) equals M2=A_bh/kappa_edd exactly, so the
# genuine Eddington-limited (exponential, kappa_edd-dependent) regime
# between them has zero width, and the two remaining branches both
# evaluate to A_bh identically -- i.e. dM/dt=A_bh (the raw, UNCAPPED
# supply rate) for every black hole mass, completely independent of
# kappa_edd/epsilon. Checked directly: at chi_crit=1, growth is bit-for-
# bit identical across epsilon=0.057/0.1/0.32, whereas at chi_crit=1e12 it
# correctly reduces to min(A_bh, kappa_edd*M) (verified against
# main._bh_growth_step's own two-regime solver, matching the existing
# test suite's own r_crit=1e12 convention for this limit). This was a
# real, previously undiscovered bug in this function's own former
# chi_crit=1.0 default -- every M_seed,crit computed against that default
# was actually uncapped supply-limited growth, not Eddington-limited
# growth as the paper describes and as this function's own prior
# docstring claimed.


def critical_seed_paper(
    halo_mass, redshift, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1.0, seed_mass_hi=1.0e8, n_iter=50,
    chi_crit=STRICT_EDDINGTON_CHI_CRIT, **run_kwargs,
):
    """
    Vectorised bisection for M_seed,crit against run_reservoir_paper --
    the paper-specific analogue of ashvini.critical_seed.critical_seed,
    using this module's own tree-based gas supply/zeta/shared-t_ff physics
    instead of main.run_forest()'s general-purpose model. See
    critical_seed.critical_seed's docstring for the definition
    (M_BH(z_anchor;M_seed,crit) = f_bh*M_star(z_anchor) under strict
    Eddington-limited growth) and the rationale for a straight bisection
    being well-posed here.

    chi_crit defaults to STRICT_EDDINGTON_CHI_CRIT (see that constant's
    own comment for why chi_crit=1.0, this function's former default, was
    wrong): pass an explicit finite chi_crit only to deliberately study a
    graded super-Eddington cap, not to represent "strict Eddington".

    Returns (M_seed_crit, never_reaches_target, already_above_at_lo), same
    semantics as critical_seed.critical_seed.
    """
    halo_mass = np.atleast_2d(halo_mass)
    N = halo_mass.shape[0]

    def excess(seed_mass):
        result = run_reservoir_paper(halo_mass, redshift, seed_mass, chi_crit=chi_crit, **run_kwargs)
        return result["bh_mass"][:, -1] - f_bh * result["stars_mass"][:, -1]

    lo = np.full(N, float(seed_mass_lo))
    hi = np.full(N, float(seed_mass_hi))

    excess_lo = excess(lo)
    excess_hi = excess(hi)
    already_above_at_lo = excess_lo > 0.0
    never_reaches_target = excess_hi < 0.0

    for _ in range(n_iter):
        mid = np.sqrt(lo * hi)
        excess_mid = excess(mid)
        go_higher = excess_mid < 0.0
        lo = np.where(go_higher, mid, lo)
        hi = np.where(go_higher, hi, mid)

    M_seed_crit = np.sqrt(lo * hi)
    unresolved = never_reaches_target | already_above_at_lo
    M_seed_crit = np.where(unresolved, np.nan, M_seed_crit)
    return M_seed_crit, never_reaches_target, already_above_at_lo
