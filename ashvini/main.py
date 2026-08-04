import os
import time
import numpy as np
import h5py

from scipy.integrate import solve_ivp

from . import utils as utils
from . import supernovae_feedback as sn

from .star_formation import star_formation_rate, time_freefall
from .gas_evolve import gas_inflow_rate, update_gas_reservoir
from .metallicity import evolve_gas_metals, evolve_stars_metals
from .dust import update_dust_reservoir

from .run_params import PARAMS, print_config

UV_background = PARAMS.reion.UVB_enabled
t_d = PARAMS.sn.delay_time  # delay time for SNe feedback, in Gyr
sn_type = PARAMS.sn.type  # type of supernova feedback
e_ff = PARAMS.sf.efficiency
IGM_metallicity = PARAMS.metals.Z_IGM
metallicity_yield = PARAMS.metals.Z_yield
Y_d = PARAMS.dust.dust_yield
Gamma = PARAMS.dust.dust_gamma
Alpha = PARAMS.dust.dust_alpha
M_crit = PARAMS.dust.m_crit
M_swept = PARAMS.dust.m_swept


tiny = 1e-15  # small number for numerical gymnastics...

method = "LSODA"


def run1_scalar(halo_mass, halo_mass_rate, redshift):
    """
    Original per-halo, per-timestep solve_ivp integrator.

    Kept as a reference implementation for validating run1()/run_forest()
    below (see tests/test_run1.py::test_vectorized_matches_scalar_reference)
    -- the same role pymctrees' build_tree() plays alongside its vectorised
    build_forest_numpy()/build_forest_numba() backends. Not used by run().
    """
    cosmic_time = utils.time_at_z(redshift)  # Gyr
    tsn = cosmic_time[0] + t_d  # Supernova switch-on time

    n = len(cosmic_time)
    gas_mass = np.zeros(n)
    gas_metals = np.zeros(n)
    stars_mass = np.zeros(n)
    stars_metals = np.zeros(n)
    sfr = np.zeros(n)
    stellar_metallicity = np.zeros(n)
    dust_mass = np.zeros(n)

    gas_accretion_rate = gas_inflow_rate(
        redshift, halo_mass, halo_mass_rate, UV_background
    )

    delay_counter = None

    for j in range(1, n):
        t_span = [cosmic_time[j - 1], cosmic_time[j]]

        if cosmic_time[j] <= tsn:
            feedback_type = "no"
            sfr_feedback = 0.0
            delay_counter = j
        else:
            feedback_type = sn_type
            sfr_feedback = (
                sfr[j - delay_counter - 1] if delay_counter is not None else 0.0
            )

        # Update gas mass
        sol = solve_ivp(
            update_gas_reservoir,
            t_span,
            [gas_mass[j - 1]],
            method=method,
            args=(
                gas_accretion_rate[j - 1],
                halo_mass[j - 1],
                stellar_metallicity[j - 1],
                sfr_feedback,
                feedback_type,
            ),
        )
        gas_mass[j] = sol.y[0, -1]

        # Update stellar mass
        sol = solve_ivp(
            lambda t, y: [star_formation_rate(t, gas_mass[j - 1])],
            t_span,
            [stars_mass[j - 1]],
            method=method,
        )
        stars_mass[j] = sol.y[0, -1]

        # Update gas metals
        if cosmic_time[j] <= tsn:
            sfr_input = sfr[j - 1]
        else:
            sfr_input = sfr[j - 1 - delay_counter]
        sol = solve_ivp(
            evolve_gas_metals,
            t_span,
            [gas_metals[j - 1]],
            method=method,
            args=(
                gas_mass[j - 1],
                gas_accretion_rate[j - 1],
                halo_mass[j - 1],
                stellar_metallicity[j - 1],
                sfr_input,
                feedback_type,
            ),
        )
        gas_metals[j] = sol.y[0, -1]

        # Update stellar metals
        sol = solve_ivp(
            lambda t, y: [evolve_stars_metals(t, gas_metals[j - 1], gas_mass[j - 1])],
            t_span,
            [stars_metals[j - 1]],
            method=method,
        )
        stars_metals[j] = sol.y[0, -1]

        # Update dust mass
        sol = solve_ivp(
            update_dust_reservoir,
            t_span,
            [dust_mass[j - 1]],
            method=method,
            args=(
                gas_mass[j - 1],
                halo_mass[j - 1],
                sfr[j - 1 - delay_counter],
                stars_mass[j - 1 - delay_counter] - stars_mass[j - 2 - delay_counter],
                stellar_metallicity[j - 1],
            ),
        )
        dust_mass[j] = sol.y[0, -1]

        # Star formation rate at current time
        sfr[j] = star_formation_rate(cosmic_time[j], gas_mass[j])

        # Enforce non-negativity
        gas_mass[j] = max(gas_mass[j], 0.0)
        gas_metals[j] = max(gas_metals[j], 0.0)
        stars_mass[j] = max(stars_mass[j], 0.0)
        stars_metals[j] = max(stars_metals[j], 0.0)
        dust_mass[j] = max(dust_mass[j], 0.0)

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


def _sf_time_integral(cosmic_time, redshift):
    """
    I[j-1] = integral_{t[j-1]}^{t[j]} dt / t_ff(z(t)), for j = 1..n-1.

    Identical for every halo in a forest -- it depends only on the shared
    cosmic_time/redshift grid, not on halo mass -- so it's computed once via
    5-point Gauss-Legendre quadrature per step and reused by every halo,
    instead of being re-integrated per halo per step as an ODE.
    """
    nodes, weights = np.polynomial.legendre.leggauss(5)
    t0, t1 = cosmic_time[:-1], cosmic_time[1:]
    dt = t1 - t0
    t_sub = 0.5 * dt[:, None] * nodes[None, :] + 0.5 * (t0 + t1)[:, None]
    z_sub = utils.z_at_time(t_sub)
    integrand = 1.0 / time_freefall(z_sub)
    return 0.5 * dt * np.sum(integrand * weights[None, :], axis=1)


def run_forest(halo_mass, halo_mass_rate, redshift):
    """
    Vectorised replacement for joblib.Parallel(run1_scalar): integrates all
    N haloes in halo_mass/halo_mass_rate (shape (N, n)) simultaneously,
    looping only over the n shared redshift steps instead of over N haloes.
    See the module-level comment above for the numerical scheme.
    """
    halo_mass = np.atleast_2d(halo_mass)
    halo_mass_rate = np.atleast_2d(halo_mass_rate)

    cosmic_time = utils.time_at_z(redshift)  # Gyr
    n = len(cosmic_time)
    N = halo_mass.shape[0]
    tsn = cosmic_time[0] + t_d

    gas_mass = np.zeros((N, n))
    gas_metals = np.zeros((N, n))
    stars_mass = np.zeros((N, n))
    stars_metals = np.zeros((N, n))
    sfr = np.zeros((N, n))
    stellar_metallicity = np.zeros((N, n))
    dust_mass = np.zeros((N, n))

    redshift_bcast = np.broadcast_to(redshift, (N, n))
    gas_accretion_rate = gas_inflow_rate(
        redshift_bcast, halo_mass, halo_mass_rate, UV_background
    )

    sf_integral = _sf_time_integral(cosmic_time, redshift)

    delay_counter = None

    for j in range(1, n):
        dt = cosmic_time[j] - cosmic_time[j - 1]
        z_mid = utils.z_at_time(0.5 * (cosmic_time[j - 1] + cosmic_time[j]))

        gm_prev = gas_mass[:, j - 1]
        hm_prev = halo_mass[:, j - 1]
        sz_prev = stellar_metallicity[:, j - 1]
        has_gas = gm_prev > 0
        gm_prev_safe = np.where(has_gas, gm_prev, 1.0)

        if cosmic_time[j] <= tsn:
            feedback_type = "no"
            delay_counter = j
        else:
            feedback_type = sn_type

        A_acc = gas_accretion_rate[:, j - 1]
        B_sf = e_ff / time_freefall(z_mid)  # scalar, shared across haloes
        ML = sn.mass_loading_factor(z_mid, hm_prev, sz_prev)  # shape (N,)

        # --- gas mass: dy/dt = A_acc - present_sfr(y) - ML*wind_sfr ---
        # (present_sfr is always self-referential -- it's B_sf*y for THIS
        # ode -- so it's always decay; only "instantaneous" wind_sfr is
        # also self-referential, the other two branches keep it as forcing)
        if feedback_type == "no":
            wind_sfr = np.zeros(N)
            gas_forcing = A_acc
            gas_decay = B_sf
        elif feedback_type == "instantaneous":
            wind_sfr = B_sf * gm_prev  # only used below, in gas_metals forcing
            gas_forcing = A_acc
            gas_decay = B_sf * (1.0 + ML)
        else:  # "delayed"
            wind_sfr = (
                sfr[:, j - delay_counter - 1]
                if delay_counter is not None
                else np.zeros(N)
            )
            gas_forcing = A_acc - ML * wind_sfr
            gas_decay = B_sf

        gas_mass[:, j] = _linear_ode_step(gm_prev, gas_forcing, gas_decay, dt)

        # --- gas metals: dy/dt = IGM*A_acc + yield*wind_sfr
        #                          - y*(B_sf + ML*wind_sfr/gas_mass) ---
        # (gas_mass here is always the frozen gm_prev, not this ode's own
        # state, so this decay/forcing split is branch-agnostic once
        # wind_sfr is known)
        forcing_gm = IGM_metallicity * A_acc + metallicity_yield * wind_sfr
        decay_gm = np.where(has_gas, B_sf + ML * wind_sfr / gm_prev_safe, 0.0)
        gas_metals[:, j] = _linear_ode_step(
            gas_metals[:, j - 1], forcing_gm, decay_gm, dt
        )

        # --- stars mass / stars metals: pure quadrature, no y-dependence ---
        I_j = sf_integral[j - 1]
        stars_mass[:, j] = stars_mass[:, j - 1] + gm_prev * e_ff * I_j
        stars_metals[:, j] = stars_metals[:, j - 1] + np.where(
            has_gas, gas_metals[:, j - 1] * e_ff * I_j, 0.0
        )

        # --- dust mass ---
        if delay_counter is not None:
            idx1 = j - 1 - delay_counter
            idx2 = j - 2 - delay_counter
            past_sfr = sfr[:, idx1]
            past_stars_mass = stars_mass[:, idx1] - stars_mass[:, idx2]
        else:
            past_sfr = np.zeros(N)
            past_stars_mass = np.zeros(N)

        has_stars = past_stars_mass > 0
        SNe_rate = np.where(
            has_stars, Gamma * past_sfr / np.where(has_stars, past_stars_mass, 1.0), 0.0
        )
        dust_loading = 1.0 - np.exp(-((gm_prev / M_crit) ** Alpha))
        dust_forcing = Y_d * past_sfr - M_swept * SNe_rate * dust_loading
        # ML here is the same mass_loading_factor(z_mid, hm_prev, sz_prev)
        # computed above for the gas-mass/gas-metals ODEs -- dust's wind loss
        # uses the same stellar metallicity, so it's reused rather than
        # recomputed.
        dust_decay = np.where(has_gas, ML * past_sfr / gm_prev_safe, 0.0)
        dust_mass[:, j] = _linear_ode_step(
            dust_mass[:, j - 1], dust_forcing, dust_decay, dt
        )

        # --- star formation rate at current (unfrozen) time/state ---
        sfr[:, j] = e_ff / time_freefall(redshift[j]) * gas_mass[:, j]

        # --- enforce non-negativity ---
        gas_mass[:, j] = np.maximum(gas_mass[:, j], 0.0)
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
        "gas_mass": gas_mass,
        "stars_mass": stars_mass,
        "gas_metals": gas_metals,
        "stars_metals": stars_metals,
        "dust_mass": dust_mass,
        "sfr": sfr,
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
    }


def run1(halo_mass, halo_mass_rate, redshift):
    """
    Vectorised single-halo integrator (fast path, replaces the old
    solve_ivp-based implementation -- see run1_scalar() for the reference
    version and the module comment above for the numerical scheme).
    """
    result = run_forest(
        np.atleast_2d(halo_mass), np.atleast_2d(halo_mass_rate), redshift
    )
    return {
        key: (val[0] if key not in ("cosmic_time", "redshift") else val)
        for key, val in result.items()
    }


def run():
    print_config(PARAMS)

    halo_masses, halo_mass_rates, redshifts = utils.read_trees(
        file_path=PARAMS.io.tree_file, mass_bin=PARAMS.io.mass_bin
    )

    N_halos = np.shape(halo_masses)[0]
    print(f"Running {N_halos} halos (vectorised across all haloes at once)...")
    t0 = time.perf_counter()
    combined = run_forest(halo_masses, halo_mass_rates, redshifts)
    print(f"Done in {time.perf_counter() - t0:.2f} s")

    # Output file
    output_file = PARAMS.io.dir_out + f"mass_bin_{PARAMS.io.mass_bin}_{sn_type}.hdf5"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with h5py.File(output_file, "w") as f:
        # Save common 1D arrays at the root
        f.create_dataset("cosmic_time", data=combined["cosmic_time"])
        f.create_dataset("redshift", data=combined["redshift"])

        # Group for halo properties
        grp = f.create_group(f"mass_bin_{PARAMS.io.mass_bin}")
        for key, val in combined.items():
            if key not in ["cosmic_time", "redshift"]:
                grp.create_dataset(key, data=val)

    print(f"Saved outputs to {output_file}")
