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
from .steppers.grid import GridStepper
from .steppers.events import EventStepper
from .components.ode import linear_ode_step as _linear_ode_step
from .components.delay import (
    effective_sn_type as _effective_sn_type,
    step_feedback_type as _step_feedback_type,
    delay_lookback_index as _delay_lookback_index,
)
from .components.blackholes import bh_growth_step as _bh_growth_step
from .dust import supernovae_per_msun
from .metallicity import igm_metallicity

tiny = 1e-15  # small number for numerical gymnastics...

method = "LSODA"


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


# The closed-form step, the delay lookback, the black hole growth step and the feedback-type helpers now
# live in ashvini/components/; these names stay for run1_scalar and older callers.


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

    stepper = GridStepper(
        config, growth_model=growth_model,
        a_star=a_star, epsilon=epsilon, r_crit=r_crit, epsilon_f=epsilon_f, eta_acc=eta_acc,
    )
    state, cosmic_time = stepper.run(
        halo_mass, halo_mass_rate, redshift, gas_mass0=gas_mass0, stars_mass0=stars_mass0
    )

    return {
        **state.as_dict(),
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
    }


def run_events(halo_mass, halo_mass_rate, redshift, config=None, sampler="smooth", seed=None,
               gas_mass0=0.0, stars_mass0=0.0, **stepper_kwargs):
    """
    The same haloes as run_forest, followed by the event integrator (steppers/events.py): individual supernovae at their own
    times. sampler is "smooth" (the expected supernovae, no randomness) or "stochastic" (seeded by `seed`).
    Returns the run_forest dictionary plus "ejected_mass" (cumulative gas removed), "returned_mass" (cumulative
    gas that fell back, zero unless events.recycling_fraction > 0) and "n_explosions" (per halo).
    """
    halo_mass = np.atleast_2d(halo_mass)
    halo_mass_rate = np.atleast_2d(halo_mass_rate)
    stepper = EventStepper(config, sampler=sampler, seed=seed, **stepper_kwargs)
    state, cosmic_time = stepper.run(halo_mass, halo_mass_rate, redshift, gas_mass0=gas_mass0, stars_mass0=stars_mass0)
    return {
        **state.as_dict(),
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
        "ejected_mass": stepper.ejected,
        "returned_mass": stepper.returned,
        "n_explosions": stepper.explosions,
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
