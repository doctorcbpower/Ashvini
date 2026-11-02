import os
import time
import numpy as np
import h5py

from scipy.integrate import solve_ivp

from . import utils as utils
from . import supernovae_feedback as sn
from . import black_holes_growth as bh_growth
from . import black_holes_growth_slimdisk as bh_slimdisk
from . import agn_feedback as agn
from .spin import epsilon_from_spin

from .star_formation import star_formation_rate, time_freefall
from .gas_evolve import gas_inflow_rate, update_gas_reservoir
from .metallicity import evolve_gas_metals, evolve_stars_metals
from .dust import update_dust_reservoir

from .run_params import print_config
from .config import resolve
from .state import State
from .dust import supernovae_per_msun
from .metallicity import igm_metallicity

tiny = 1e-15  # small number for numerical gymnastics...

method = "LSODA"


def _effective_sn_type(cfg):
    """
    Supernova feedback type actually applied. A "delayed" wind with a zero delay is instantaneous
    feedback: the delayed lookback would otherwise point at the current step, whose star formation
    rate is not yet known, and silently switch feedback off.
    """
    sn_type, t_d = cfg.sn.type, cfg.sn.delay_time
    return "instantaneous" if (sn_type == "delayed" and t_d <= 0.0) else sn_type


def _step_feedback_type(eff_type, has_delay_history):
    """Feedback type for one step: "instantaneous" and "no" apply at every step; "delayed" only once a delay has elapsed."""
    if eff_type == "delayed" and not has_delay_history:
        return "no"
    return eff_type


def run1_scalar(halo_mass, halo_mass_rate, redshift, config=None):
    """
    Original per-halo, per-timestep solve_ivp integrator.

    Kept as a reference implementation for validating run1()/run_forest()
    below (see tests/test_run1.py::test_vectorized_matches_scalar_reference)
    -- the same role pymctrees' build_tree() plays alongside its vectorised
    build_forest_numpy()/build_forest_numba() backends. Not used by run().
    """
    cfg = resolve(config)
    cosmic_time = utils.time_at_z(redshift)  # Gyr

    n = len(cosmic_time)
    gas_mass = np.zeros(n)
    gas_metals = np.zeros(n)
    stars_mass = np.zeros(n)
    stars_metals = np.zeros(n)
    sfr = np.zeros(n)
    stellar_metallicity = np.zeros(n)
    dust_mass = np.zeros(n)
    bh_mass = np.zeros(n)

    gas_accretion_rate = gas_inflow_rate(
        redshift, halo_mass, halo_mass_rate, cfg=config
    )

    eff_sn_type = _effective_sn_type(cfg)
    growth_cap = cfg.bh.sigma_feedback.growth_cap_enabled
    delay_idx = _delay_lookback_index(cosmic_time, cfg.sn.delay_time)
    agn_delay_idx = _delay_lookback_index(cosmic_time, cfg.bh.feedback_delay_time)
    bh_growth_rate_history = np.zeros(n)

    for j in range(1, n):
        t_span = [cosmic_time[j - 1], cosmic_time[j]]

        has_delay_history = delay_idx[j] >= 0
        feedback_type = _step_feedback_type(eff_sn_type, has_delay_history)
        # The star formation rate that drives supernova feedback, for gas, metals and dust alike:
        # zero ("no"), the delayed rate ("delayed"), or the present rate ("instantaneous", frozen at the
        # start of the step for the metals and dust, as in the vectorised path).
        if feedback_type == "delayed":
            sfr_feedback = sfr[delay_idx[j]]
        elif feedback_type == "instantaneous":
            sfr_feedback = star_formation_rate(t_span[0], gas_mass[j - 1], cfg=cfg)
        else:
            sfr_feedback = 0.0

        # Black hole seeding + growth (computed before gas mass, since gas
        # mass's AGN wind term needs this step's BH growth rate)
        gas_metallicity_prev = (
            gas_metals[j - 1] / gas_mass[j - 1] if gas_mass[j - 1] > 0 else 0.0
        )
        newly_seeded, seed_mass = bh_growth.seeding_mask_and_mass(
            halo_mass[j - 1],
            redshift[j - 1],
            gas_metallicity_prev,
            bh_mass[j - 1] > 0,
            cfg=cfg,
        )
        bh_mass_y0 = float(seed_mass) if bool(newly_seeded) else bh_mass[j - 1]

        sol = solve_ivp(
            bh_growth.black_hole_growth_rate,
            t_span,
            [bh_mass_y0],
            method=method,
            args=(gas_mass[j - 1], cfg),
        )
        bh_mass[j] = max(sol.y[0, -1], 0.0)
        # Optional hard cap at PZNK11's predicted ceiling (see run_forest's
        # matching comment): dM_BH/dt is never negative in this model, so
        # clipping solve_ivp's endpoint result is exactly equivalent to a
        # "freeze once M_BH reaches the ceiling" ODE, without needing to
        # thread the cap into black_hole_growth_rate's own RHS.
        if growth_cap:
            ceiling = bh_growth.growth_ceiling(halo_mass[j - 1], redshift[j - 1], cfg)
            bh_mass[j] = min(bh_mass[j], float(ceiling))
        agn_growth_rate = (bh_mass[j] - bh_mass_y0) / (t_span[1] - t_span[0])
        bh_growth_rate_history[j] = agn_growth_rate

        has_agn_delay_history = agn_delay_idx[j] >= 0
        agn_wind_growth_rate = (
            bh_growth_rate_history[agn_delay_idx[j]] if has_agn_delay_history else 0.0
        )

        # Step-midpoint M_BH, for agn.agn_wind_mass_rate's optional
        # M-sigma self-regulation coupling (sigma_feedback_enabled) only --
        # a closed-form estimate (same _bh_growth_step used by
        # run_forest(), evaluated at dt/2 instead of dt) rather than
        # bh_mass[j] itself, so the AGN coupling doesn't depend on the
        # step's own not-yet-known end state. Cheap and unused when
        # sigma_feedback is disabled (the default).
        z_mid = utils.z_at_time(0.5 * (t_span[0] + t_span[1]))
        A_bh_mid = (cfg.bh.efficiency / bh_growth.time_freefall(z_mid)) * gas_mass[j - 1]
        kappa_edd = bh_growth.EDDINGTON_RATE_PER_UNIT_MASS * cfg.bh.eddington_multiplier
        bh_mass_mid = float(
            _bh_growth_step(bh_mass_y0, A_bh_mid, kappa_edd, 0.5 * (t_span[1] - t_span[0]))
        )
        if growth_cap:
            bh_mass_mid = min(bh_mass_mid, float(bh_growth.growth_ceiling(halo_mass[j - 1], z_mid, cfg)))

        # Update gas mass and stellar mass together: the stars formed in the step are the integral of
        # the star formation rate of the evolving gas, not of the gas mass at the start of the step.
        gas_args = (
            gas_accretion_rate[j - 1],
            halo_mass[j - 1],
            stellar_metallicity[j - 1],
            sfr_feedback,
            feedback_type,
            agn_growth_rate,       # instantaneous: BH's own accretion, a gas-mass sink
            agn_wind_growth_rate,  # possibly delayed: drives the AGN wind term
            bh_mass_mid,           # step-midpoint M_BH, for the AGN coupling switch
            config,
        )
        sol = solve_ivp(
            lambda t, y: [
                float(update_gas_reservoir(t, y[0], *gas_args)),
                float(star_formation_rate(t, y[0], cfg=cfg)),
            ],
            t_span,
            [gas_mass[j - 1], stars_mass[j - 1]],
            method=method,
        )
        gas_mass[j], stars_mass[j] = sol.y[0, -1], sol.y[1, -1]

        # Update gas metals and stellar metals together (metals are locked into stars at e_ff/t_ff * M_Z)
        metal_args = (
            gas_mass[j - 1],
            gas_accretion_rate[j - 1],
            halo_mass[j - 1],
            stellar_metallicity[j - 1],
            sfr_feedback,
            feedback_type,
            cfg,
        )
        sol = solve_ivp(
            lambda t, y: [
                float(evolve_gas_metals(t, y[0], *metal_args)),
                float(evolve_stars_metals(t, y[0], cfg=cfg)),
            ],
            t_span,
            [gas_metals[j - 1], stars_metals[j - 1]],
            method=method,
        )
        gas_metals[j], stars_metals[j] = sol.y[0, -1], sol.y[1, -1]

        # Update dust mass (driven by the same star formation rate as the supernova wind)
        sol = solve_ivp(
            update_dust_reservoir,
            t_span,
            [dust_mass[j - 1]],
            method=method,
            args=(
                gas_mass[j - 1],
                halo_mass[j - 1],
                sfr_feedback,
                stellar_metallicity[j - 1],
                cfg,
            ),
        )
        dust_mass[j] = sol.y[0, -1]

        # Enforce non-negativity (gas_mass first -- see run_forest's
        # matching comment for why sfr must be computed from the clamped
        # value, not solve_ivp's raw output)
        gas_mass[j] = max(gas_mass[j], 0.0)

        # Star formation rate at current time
        sfr[j] = star_formation_rate(cosmic_time[j], gas_mass[j], cfg=cfg)

        gas_metals[j] = max(gas_metals[j], 0.0)
        stars_mass[j] = max(stars_mass[j], 0.0)
        stars_metals[j] = max(stars_metals[j], 0.0)
        dust_mass[j] = max(dust_mass[j], 0.0)
        bh_mass[j] = max(bh_mass[j], 0.0)

        # No gas metals or dust if no gas
        if gas_mass[j] <= 0:
            gas_metals[j] = 0.0
            dust_mass[j] = 0.0

        # Stellar metallicity
        if stars_mass[j] > 0:
            stellar_metallicity[j] = stars_metals[j] / stars_mass[j]
        else:
            stellar_metallicity[j] = 0.0
            stars_metals[j] = 0.0

    return {
        "gas_mass": gas_mass,
        "stars_mass": stars_mass,
        "gas_metals": gas_metals,
        "stars_metals": stars_metals,
        "dust_mass": dust_mass,
        "bh_mass": bh_mass,
        "sfr": sfr,
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
    }


# ---------------------------------------------------------------------------
# Vectorised ODE stepping
# ---------------------------------------------------------------------------
# run1_scalar() above profiles at ~0.15-0.28s / halo for a 200-step history,
# dominated by solve_ivp call overhead (995 calls/halo, 5 per step -- see the
# code audit, Section 3). Two of the five ODEs (stars_mass, stars_metals)
# don't actually depend on the solved state within a step -- their RHS only
# depends on t (through the redshift lookup) and on quantities already fixed
# at the start of the step -- so they're pure integrals over t, replaced
# below by closed-form quadrature. The other three (gas_mass, gas_metals,
# dust_mass) are affine ODEs of the form dy/dt = forcing - decay*y, where
# forcing/decay vary slowly across a step through their redshift dependence;
# freezing them at the step midpoint admits a closed-form exponential-
# integrator update in place of solve_ivp. Because every halo in a forest
# shares the same cosmic_time/redshift grid, this collapses the entire
# per-halo Python loop into array arithmetic over all N haloes at once,
# replacing joblib+solve_ivp with vectorised NumPy -- the same strategy
# pymctrees already uses for its own tree-building hot loop.
#
# This is a genuine change in the numerical integration scheme (not just an
# interpolation swap), so run1_scalar() above is kept as a slow-but-trusted
# reference and tests/test_run1.py cross-checks the two implementations
# agree to within the accuracy expected of this approximation.
# ---------------------------------------------------------------------------


def _linear_ode_step(y0, forcing, decay, dt):
    """
    Closed-form update of dy/dt = forcing - decay*y over a step of size dt,
    assuming forcing/decay are constant across the step. Vectorised and safe
    for decay == 0 (no supernova/mass-loading term active): rather than
    dividing by decay directly, phi(x) = (1 - exp(-x)) / x is evaluated via
    a Taylor series near x = decay*dt = 0, which is well-defined there.
    """
    x = decay * dt
    exp_term = np.exp(-x)
    small = np.abs(x) < 1e-6
    phi = np.where(
        small,
        dt * (1.0 - x / 2.0 + x * x / 6.0),
        (1.0 - exp_term) / np.where(x == 0, 1.0, decay),
    )
    return y0 * exp_term + forcing * phi


def _delay_lookback_index(cosmic_time, t_d):
    """
    For each step j, the largest index k with cosmic_time[k] <= cosmic_time[j] - t_d
    (i.e. the most recent tabulated state at or before t_d in the past),
    or -1 if less than t_d has elapsed since the start of the run.

    Regression fix for a real bug: this used to be a single index offset
    (delay_counter) computed once, from however many *array steps* elapsed
    before crossing t_d early in the run, then reused as a fixed offset for
    the rest of the integration. That's only correct on a uniform
    cosmic-time grid -- trees are stepped uniformly in *redshift*, and
    dt/dz varies enormously with z (e.g. ~40x between z=25 and z=5 for a
    typical cosmology), so a fixed step-offset calibrated near z_max
    represented an effective delay tens of times longer than t_d by the
    time the integration reached low z (verified directly: 182 steps
    calibrated to 15 Myr near z=25 corresponded to ~530 Myr by z=5 for a
    z_max=25 pymctrees-adapted tree). That silently washed out the
    delayed-feedback oscillations the model is specifically built to
    produce. This recomputes the lookback per-step from actual elapsed
    time instead of a step count.
    """
    n = len(cosmic_time)
    target_t = cosmic_time - t_d
    idx = np.searchsorted(cosmic_time, target_t, side="right") - 1
    return np.clip(idx, -1, n - 1)


def _bh_growth_step(y0, A_bh, kappa_edd, dt):
    """
    Closed-form update of dM_BH/dt = min(A_bh, kappa_edd*y) over a step of
    size dt, with A_bh (gas-supply-limited rate) frozen over the step and
    kappa_edd a true constant (the Eddington-limited growth rate per unit
    mass). Unlike _linear_ode_step, this ODE is piecewise (not affine) --
    a halo below the crossover mass A_bh/kappa_edd grows exponentially
    (Eddington-limited) until it reaches it, then grows linearly at rate
    A_bh (gas-supply-limited) for the remainder of the step. Solved exactly
    (not by picking one regime for the whole step) so a step that happens
    to straddle the crossover doesn't overshoot.

    y0 == 0 (not yet seeded) is handled without special-casing: the
    Eddington branch is multiplicative, so 0 stays exactly 0.
    """
    threshold_mass = A_bh / kappa_edd  # A_bh >= 0 always, so this is safe
    eddington_limited = y0 < threshold_mass

    with np.errstate(divide="ignore", invalid="ignore"):
        y0_safe = np.where(y0 > 0, y0, 1.0)
        t_star = np.where(
            y0 > 0,
            np.log(threshold_mass / y0_safe) / kappa_edd,
            np.inf,
        )

    crosses_over = eddington_limited & (t_star < dt)
    stays_eddington = eddington_limited & ~crosses_over
    gas_limited = ~eddington_limited

    # np.select evaluates every choice array for every element regardless of
    # which branch is picked, so t_star's inf entries (uncrossed haloes)
    # would otherwise hit inf * 0 = nan under "dt - t_star" and warn, even
    # though that value is never selected. Clip to a finite stand-in first.
    t_star_safe = np.where(np.isfinite(t_star), t_star, dt)

    return np.select(
        [stays_eddington, crosses_over, gas_limited],
        [
            y0 * np.exp(kappa_edd * dt),
            threshold_mass + A_bh * (dt - t_star_safe),
            y0 + A_bh * dt,
        ],
    )


def run_forest(halo_mass, halo_mass_rate, redshift, growth_model=None,
               a_star=None, epsilon=None, r_crit=None, epsilon_f=None, eta_acc=None,
               gas_mass0=0.0, stars_mass0=0.0, config=None):
    """
    Vectorised replacement for joblib.Parallel(run1_scalar): integrates all
    N haloes in halo_mass/halo_mass_rate (shape (N, n)) simultaneously,
    looping only over the n shared redshift steps instead of over N haloes.
    See the module-level comment above for the numerical scheme.

    growth_model, a_star, epsilon, r_crit, epsilon_f, eta_acc : optional
    overrides for the "hobbs_slimdisk" branch's physics parameters, each
    defaulting to None (= use the value in the configuration: config.bh.slimdisk).
    black_holes_growth_slimdisk.py's own functions (nuclear_accretion_rate,
    eddington_rate_std_per_unit_mass, king_agn_wind_rate) already accept
    these as call-time kwargs; what was missing is that this function
    always called them with no override -- these are exactly the parameters
    the 2026 paper's sensitivity scans (spin/epsilon, r_crit, epsilon_f,
    eta_acc/f_am) need to vary run-to-run. growth_model overrides
    config.bh.growth_model; a_star/epsilon/r_crit/epsilon_f/eta_acc are only read when the
    effective growth_model is "hobbs_slimdisk" (silently unused otherwise,
    same as the existing bh.slimdisk config block).

    epsilon, if given, sets the radiative efficiency directly, bypassing
    spin.epsilon_from_spin(a_star) entirely -- the paper's own fiducial
    calculation uses epsilon=0.1 as a free parameter in its own right (see
    Section "Integration of the reservoir equations"), independent of the
    a_star=0.5 -> epsilon=0.0821 spin-derived value used specifically for
    the Novikov-Thorne spin/angular-momentum sensitivity scan (Section
    5.5) -- these are two different, both legitimate, uses of the same
    kappa_edd normalisation, not a bug in either. Passing both a_star and
    epsilon raises ValueError (ambiguous which should win).

    config : Config, optional
    All other parameters come from this configuration (default: config.DEFAULT, the one loaded from
    run_params.yaml). See config.py; nothing is read from module-level copies.

    gas_mass0, stars_mass0 : initial gas/stellar reservoir masses (Msun)
    at redshift[0] (z_seed), each defaulting to 0.0 (this function's prior,
    only behaviour). The paper's own fiducial calculation instead starts
    from M_gas(z_seed)=1e5 Msun, M_star(z_seed)=1e3 Msun (small but
    nonzero, representing whatever pre-existing reservoir a freshly
    resolved halo is assumed to already hold) -- scalars (same initial
    value for every halo) or (N,)-shaped arrays (per-halo).
    """
    if a_star is not None and epsilon is not None:
        raise ValueError("run_forest: pass at most one of a_star, epsilon (ambiguous which should set kappa_edd).")

    halo_mass = np.atleast_2d(halo_mass)
    halo_mass_rate = np.atleast_2d(halo_mass_rate)

    cfg = resolve(config)
    growth_model = growth_model if growth_model is not None else cfg.bh.growth_model
    slim = cfg.bh.slimdisk
    growth_cap = cfg.bh.sigma_feedback.growth_cap_enabled

    # Resolved once per call (not per step): these are run-level constants,
    # only actually used below when growth_model == "hobbs_slimdisk".
    _r_crit = r_crit if r_crit is not None else slim.r_crit
    _epsilon_f = epsilon_f if epsilon_f is not None else slim.epsilon_f
    _eta_acc = eta_acc if eta_acc is not None else slim.eta_acc
    if epsilon is not None:
        _kappa_edd_slimdisk = bh_slimdisk.eddington_rate_std_per_unit_mass(epsilon)
    else:
        _kappa_edd_slimdisk = bh_slimdisk.eddington_rate_std_per_unit_mass(
            epsilon_from_spin(a_star if a_star is not None else slim.a_star)
        )

    cosmic_time = utils.time_at_z(redshift)  # Gyr
    n = len(cosmic_time)
    N = halo_mass.shape[0]

    state = State.empty(N, n)
    gas_mass, gas_metals = state.gas_mass, state.gas_metals
    stars_mass, stars_metals = state.stars_mass, state.stars_metals
    sfr, stellar_metallicity = state.sfr, state.stellar_metallicity
    dust_mass, bh_mass = state.dust_mass, state.bh_mass
    gas_mass[:, 0] = gas_mass0
    stars_mass[:, 0] = stars_mass0

    IGM_metallicity = igm_metallicity(cfg)
    metallicity_yield = cfg.metals.Z_yield
    Y_d = cfg.dust.dust_yield
    Alpha, M_crit, M_swept = cfg.dust.dust_alpha, cfg.dust.m_crit, cfg.dust.m_swept
    SN_per_Msun = supernovae_per_msun(cfg)
    e_ff = cfg.sf.efficiency
    e_bh = cfg.bh.efficiency
    kappa_edd_pznk11 = bh_growth.EDDINGTON_RATE_PER_UNIT_MASS * cfg.bh.eddington_multiplier

    redshift_bcast = np.broadcast_to(redshift, (N, n))
    gas_accretion_rate = gas_inflow_rate(
        redshift_bcast, halo_mass, halo_mass_rate, cfg=config
    )

    eff_sn_type = _effective_sn_type(cfg)
    delay_idx = _delay_lookback_index(cosmic_time, cfg.sn.delay_time)
    agn_delay_idx = _delay_lookback_index(cosmic_time, cfg.bh.feedback_delay_time)
    bh_growth_rate_history = np.zeros((N, n))

    for j in range(1, n):
        dt = cosmic_time[j] - cosmic_time[j - 1]
        z_mid = utils.z_at_time(0.5 * (cosmic_time[j - 1] + cosmic_time[j]))

        gm_prev = gas_mass[:, j - 1]
        hm_prev = halo_mass[:, j - 1]
        sz_prev = stellar_metallicity[:, j - 1]
        has_gas = gm_prev > 0
        gm_prev_safe = np.where(has_gas, gm_prev, 1.0)

        has_delay_history = delay_idx[j] >= 0
        feedback_type = _step_feedback_type(eff_sn_type, has_delay_history)

        A_acc = gas_accretion_rate[:, j - 1]
        B_sf = e_ff / time_freefall(z_mid)  # scalar, shared across haloes
        ML = sn.mass_loading_factor(z_mid, hm_prev, sz_prev, cfg)  # shape (N,)

        # --- black hole seeding + growth (computed before gas mass, since
        # gas mass's AGN wind term needs this step's BH growth rate) ---
        gas_metallicity_prev = np.where(
            has_gas, gas_metals[:, j - 1] / gm_prev_safe, 0.0
        )
        newly_seeded, seed_mass = bh_growth.seeding_mask_and_mass(
            hm_prev, redshift[j - 1], gas_metallicity_prev, bh_mass[:, j - 1] > 0, cfg=cfg
        )
        bh_mass_prev = np.where(newly_seeded, seed_mass, bh_mass[:, j - 1])

        if growth_model == "hobbs_slimdisk":
            # 2026-paper alternative: enclosed-mass free-fall supply,
            # graded r_crit slim-disc cap, King (2003) energy-driven
            # feedback -- see black_holes_growth_slimdisk.py. No PZNK11
            # growth-ceiling analogue exists for this model.
            sm_prev = stars_mass[:, j - 1]
            A_bh = bh_slimdisk.nuclear_accretion_rate(
                gm_prev, bh_mass_prev, sm_prev, eta_acc=_eta_acc,
                R_nuc_pc=slim.R_nuc_pc, compaction_boost=slim.compaction_boost,
            )
            kappa_edd = _kappa_edd_slimdisk
            bh_mass_uncapped = bh_slimdisk.bh_growth_step_slimdisk(
                bh_mass_prev, A_bh, kappa_edd, _r_crit, dt
            )
            bh_mass[:, j] = np.maximum(bh_mass_uncapped, 0.0)
            agn_growth_rate = (bh_mass[:, j] - bh_mass_prev) / dt
            bh_growth_rate_history[:, j] = agn_growth_rate

            has_agn_delay_history = agn_delay_idx[j] >= 0
            agn_wind_growth_rate = (
                bh_growth_rate_history[:, agn_delay_idx[j]] if has_agn_delay_history else np.zeros(N)
            )
            agn_forcing = bh_slimdisk.king_agn_wind_rate(
                agn_wind_growth_rate, halo_mass=hm_prev, redshift=z_mid, epsilon_f=_epsilon_f
            )
        else:
            A_bh = (e_bh / bh_growth.time_freefall(z_mid)) * gm_prev
            kappa_edd = kappa_edd_pznk11
            bh_mass_uncapped = _bh_growth_step(bh_mass_prev, A_bh, kappa_edd, dt)
            # Optional hard cap at PZNK11's predicted M_sigma-relation ceiling
            # (black_holes.sigma_feedback.growth_cap_enabled, off by default --
            # see black_holes_growth.growth_ceiling). _bh_growth_step's result
            # is monotonically non-decreasing in y0 (growth rate is never
            # negative), so simply clipping to the ceiling here is exactly
            # equivalent to a "freeze once M_BH reaches the ceiling" ODE, no
            # separate crossing-time solve needed -- and correctly lets a
            # previously-capped BH resume growing if the ceiling itself rises
            # (M_sigma grows with the halo), rather than freezing permanently.
            if growth_cap:
                ceiling = bh_growth.growth_ceiling(hm_prev, z_mid, cfg)
                bh_mass_uncapped = np.minimum(bh_mass_uncapped, ceiling)
            bh_mass[:, j] = np.maximum(bh_mass_uncapped, 0.0)
            # Instantaneous: this step's own accretion, always -- the mass
            # leaves the gas reservoir when it's actually accreted, regardless
            # of whether the *wind* it powers is delayed (see below).
            agn_growth_rate = (bh_mass[:, j] - bh_mass_prev) / dt
            bh_growth_rate_history[:, j] = agn_growth_rate

            # Possibly delayed: the AGN wind mirrors SN's delayed-feedback
            # mechanism (PARAMS.bh.feedback_delay_time, default 0.0 =
            # instantaneous, in which case agn_delay_idx[j] == j always and
            # this reduces exactly to agn_growth_rate).
            has_agn_delay_history = agn_delay_idx[j] >= 0
            agn_wind_growth_rate = (
                bh_growth_rate_history[:, agn_delay_idx[j]] if has_agn_delay_history else np.zeros(N)
            )
            # Step-midpoint M_BH (same closed form as bh_mass[:, j] above, but
            # evaluated at dt/2), for agn.agn_wind_mass_rate's optional
            # M-sigma self-regulation coupling only -- see
            # gas_evolve.update_gas_reservoir's bh_mass_for_wind docstring for
            # why this must be frozen mid-step rather than left to vary with
            # the gas ODE's own state.
            bh_mass_mid = _bh_growth_step(bh_mass_prev, A_bh, kappa_edd, 0.5 * dt)
            if growth_cap:
                bh_mass_mid = np.minimum(bh_mass_mid, bh_growth.growth_ceiling(hm_prev, z_mid, cfg))
            agn_forcing = agn.agn_wind_mass_rate(
                agn_wind_growth_rate, bh_mass=bh_mass_mid, halo_mass=hm_prev, redshift=z_mid,
                **({} if config is None else {"cfg": config}),
            )

        # --- the star formation rate that drives supernova feedback, for gas, metals and dust alike:
        # zero ("no"), the present rate B_sf * M_gas ("instantaneous"), or the rate delay_time ago
        # ("delayed"). Dust follows exactly the same switch as the gas and metals. ---
        if feedback_type == "no":
            wind_sfr = np.zeros(N)
        elif feedback_type == "instantaneous":
            wind_sfr = B_sf * gm_prev
        else:  # "delayed"
            wind_sfr = sfr[:, delay_idx[j]]

        # --- gas mass: dy/dt = A_acc - B_sf*y - ML*wind_sfr - bh_accretion - AGN wind ---
        # (B_sf*y, star formation, is always self-referential, so it is always decay; only
        # "instantaneous" wind_sfr is also self-referential, the other two branches keep it as
        # forcing; bh_accretion (this step's own BH growth, a gas-mass sink -- see
        # gas_evolve.update_gas_reservoir's docstring) and the AGN wind term are forcing regardless)
        if feedback_type == "instantaneous":
            gas_forcing = A_acc - agn_growth_rate - agn_forcing
            gas_decay = B_sf * (1.0 + ML)
        else:
            gas_forcing = A_acc - ML * wind_sfr - agn_growth_rate - agn_forcing
            gas_decay = B_sf

        gas_raw = _linear_ode_step(gm_prev, gas_forcing, gas_decay, dt)
        gas_mass[:, j] = gas_raw

        # --- stars: the gas reservoir is consumed by star formation (B_sf) and, for instantaneous
        # feedback, by the wind; integrating dy/dt = F - D*y gives D * integral(y dt) = y0 + F*dt - y1
        # exactly, so the stars formed are (B_sf / D) * (y0 + F*dt - y1): no first-order lag, and
        # gas + stars + removed mass balance to round-off. A reservoir that is empty at the start
        # of the step but fills during it forms stars in that step. Where a negative-forcing step
        # would drive the raw solution below zero, the budget is clamped at zero stars.
        stars_formed = np.maximum(B_sf / gas_decay * (gm_prev + gas_forcing * dt - gas_raw), 0.0)
        stars_mass[:, j] = stars_mass[:, j - 1] + stars_formed

        # --- gas metals: dy/dt = IGM*A_acc + yield*wind_sfr - y*(B_sf + ML*wind_sfr/gas_mass) ---
        # (gas_mass here is always the frozen gm_prev, not this ode's own state; the locking term
        # B_sf*y needs no gas mass, so it also applies to a reservoir that starts the step empty)
        forcing_gm = IGM_metallicity * A_acc + metallicity_yield * wind_sfr
        decay_gm = B_sf + np.where(has_gas, ML * wind_sfr / gm_prev_safe, 0.0)
        gas_metals_raw = _linear_ode_step(gas_metals[:, j - 1], forcing_gm, decay_gm, dt)
        gas_metals[:, j] = gas_metals_raw
        # --- stellar metals: locked at B_sf * y, so the same budget identity applies ---
        metals_locked = np.maximum(
            B_sf / decay_gm * (gas_metals[:, j - 1] + forcing_gm * dt - gas_metals_raw), 0.0
        )
        stars_metals[:, j] = stars_metals[:, j - 1] + metals_locked

        # --- dust mass: dy/dt = Y_d*wind_sfr
        #        - y * wind_sfr/gas_mass * (nu_SN * M_swept * f_gas + ML) ---
        # Production, destruction and wind removal are all driven by wind_sfr (see above).
        # Destruction is the supernova rate nu_SN * wind_sfr times the dust-to-gas ratio y/gas_mass
        # times the swept-up mass (the dust in the swept-up gas), gated by f_gas.
        # ML is the same mass_loading_factor(z_mid, hm_prev, sz_prev) computed above.
        dust_loading = 1.0 - np.exp(-((gm_prev / M_crit) ** Alpha))
        dust_forcing = Y_d * wind_sfr
        dust_decay = np.where(
            has_gas, (SN_per_Msun * M_swept * dust_loading + ML) * wind_sfr / gm_prev_safe, 0.0
        )
        dust_mass[:, j] = _linear_ode_step(
            dust_mass[:, j - 1], dust_forcing, dust_decay, dt
        )

        # --- enforce non-negativity (gas_mass first: sfr below must be
        # computed from the clamped value, not the closed-form ODE
        # update's raw output, which legitimately overshoots slightly
        # negative when decay is fast relative to dt -- computing sfr from
        # the unclamped value silently produced a negative sfr that was
        # never itself clamped, later poisoning a delayed-feedback/dust
        # decay term negative and overflowing _linear_ode_step's exp(-x)
        # (verified directly: 188,702/20,000,000 (N x n_steps) raw
        # gas_mass updates go negative in a typical 10,000-halo run) ---
        gas_mass[:, j] = np.maximum(gas_mass[:, j], 0.0)

        # --- star formation rate at current (unfrozen) time/state ---
        sfr[:, j] = e_ff / time_freefall(redshift[j]) * gas_mass[:, j]

        gas_metals[:, j] = np.maximum(gas_metals[:, j], 0.0)
        stars_mass[:, j] = np.maximum(stars_mass[:, j], 0.0)
        stars_metals[:, j] = np.maximum(stars_metals[:, j], 0.0)
        dust_mass[:, j] = np.maximum(dust_mass[:, j], 0.0)

        no_gas = gas_mass[:, j] <= 0
        gas_metals[:, j] = np.where(no_gas, 0.0, gas_metals[:, j])
        dust_mass[:, j] = np.where(no_gas, 0.0, dust_mass[:, j])

        has_stars_mass = stars_mass[:, j] > 0
        stellar_metallicity[:, j] = np.where(
            has_stars_mass,
            stars_metals[:, j] / np.where(has_stars_mass, stars_mass[:, j], 1.0),
            0.0,
        )
        stars_metals[:, j] = np.where(has_stars_mass, stars_metals[:, j], 0.0)

    return {
        **state.as_dict(),
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
    }


def run1(halo_mass, halo_mass_rate, redshift, config=None):
    """
    Vectorised single-halo integrator (fast path, replaces the old
    solve_ivp-based implementation -- see run1_scalar() for the reference
    version and the module comment above for the numerical scheme).
    """
    result = run_forest(
        np.atleast_2d(halo_mass), np.atleast_2d(halo_mass_rate), redshift, config=config
    )
    return {
        key: (val[0] if key not in ("cosmic_time", "redshift") else val)
        for key, val in result.items()
    }


def run(config=None):
    cfg = resolve(config)
    print_config(cfg)

    if cfg.io.tree_source == "pymctrees":
        from . import pymctrees_adapter

        p = cfg.io.pymctrees
        if p is None:
            raise ValueError(
                "basics.tree_source is 'pymctrees' but basics.pymctrees is "
                "not configured in run_params.yaml -- see PymctreesSourceParams "
                "in run_params.py for the required fields (at minimum 'config', "
                "a path to a pymctrees YAML config)."
            )
        print(f"\nGenerating trees live via pymctrees ({p.config}) instead of reading a file...")
        halo_masses, halo_mass_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
            pymctrees_config_path=p.config, mass_bin=cfg.io.mass_bin,
            n_halos=p.n_halos, z0=p.z0, z_max=p.z_max, dz=p.dz,
            m_res=p.m_res, backend=p.backend, seed=p.seed,
        )
        # _merger_mass (per-halo, per-step discrete merger events) is
        # already folded into halo_mass_rates above -- returned separately
        # too, for any future merger-triggered physics (e.g. a starburst on
        # top of the baseline accretion-driven rate), not consumed that way
        # by run_forest() yet. Available here for callers that want it
        # directly via pymctrees_adapter.build_forest_live() instead of
        # run().
    elif cfg.io.tree_source == "file":
        halo_masses, halo_mass_rates, redshifts = utils.read_trees(
            file_path=cfg.io.tree_file, mass_bin=cfg.io.mass_bin
        )
    else:
        raise ValueError(
            f"Unknown basics.tree_source '{cfg.io.tree_source}' -- must be "
            "'file' or 'pymctrees'."
        )

    N_halos = np.shape(halo_masses)[0]
    print(f"Running {N_halos} halos (vectorised across all haloes at once)...")
    t0 = time.perf_counter()
    combined = run_forest(halo_masses, halo_mass_rates, redshifts, config=cfg)
    print(f"Done in {time.perf_counter() - t0:.2f} s")

    # Output file
    output_file = cfg.io.dir_out + f"mass_bin_{cfg.io.mass_bin}_{cfg.sn.type}.hdf5"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with h5py.File(output_file, "w") as f:
        # Save common 1D arrays at the root
        f.create_dataset("cosmic_time", data=combined["cosmic_time"])
        f.create_dataset("redshift", data=combined["redshift"])

        # Group for halo properties
        grp = f.create_group(f"mass_bin_{cfg.io.mass_bin}")
        for key, val in combined.items():
            if key not in ["cosmic_time", "redshift"]:
                grp.create_dataset(key, data=val)

    print(f"Saved outputs to {output_file}")

from .config import retire_module_values as _retire  # noqa: E402

_retire(__name__, {
    "t_d": "config.sn.delay_time", "sn_type": "config.sn.type", "UV_background": "config.reion.UVB_enabled",
    "agn_delay_time": "config.bh.feedback_delay_time", "e_ff": "config.sf.efficiency",
    "IGM_metallicity": "config.metals.Z_IGM_Zsun", "metallicity_yield": "config.metals.Z_yield",
    "Y_d": "config.dust.dust_yield", "SN_per_Msun": "config.sn.m_sn_max (or imf, m_sn_min)",
    "Alpha": "config.dust.dust_alpha", "M_crit": "config.dust.m_crit", "M_swept": "config.dust.m_swept",
})
