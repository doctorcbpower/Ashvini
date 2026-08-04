import numpy as np

from .utils import Hubble_time, z_at_time

from .run_params import PARAMS

e_bh = PARAMS.bh.efficiency  # efficiency of black hole growth
c = 3.0e10  # speed of light in CGS (cm/s)
G = 6.674e-8  # gravitational constant in CGS
m_p = 1.67e-24  # proton mass in CGS
sigma_thomson = 6.65e-25  # Thomson cross section in CGS

Gyr_s = 3.15576e16  # seconds per Gyr


def time_freefall(redshift):
    return 0.141 * Hubble_time(redshift)


def eddington_bh_growth(M_BH):
    """
    Eddington-limited BH growth rate, in Msun/Gyr, for M_BH in Msun.

    The underlying rate is linear in mass (dM/dt = kappa * M with kappa
    computed in CGS), so mass units cancel out of the formula entirely and
    only the time unit (s -> Gyr) needs converting -- this used to return
    g/s (CGS) while being min()'d against a Msun/Gyr rate elsewhere, making
    the Eddington cap physically meaningless.
    """
    kappa_per_s = 4 * np.pi * G * m_p / sigma_thomson / c  # 1/s, mass-independent
    kappa_per_gyr = kappa_per_s * Gyr_s
    return kappa_per_gyr * np.asarray(M_BH)


EDDINGTON_RATE_PER_UNIT_MASS = float(eddington_bh_growth(1.0))  # 1/Gyr


def black_hole_growth_rate(t, M_BH, gas_mass):
    """
    dM_BH/dt (Msun/Gyr): gas-supply-limited growth capped at the Eddington
    rate. M_BH (the ODE state) is the first argument after t, matching the
    calling convention the other RHS functions in this package
    (update_gas_reservoir, evolve_gas_metals, ...) already use for
    solve_ivp.
    """
    redshift = z_at_time(t)
    growth_rate = (e_bh / time_freefall(redshift)) * gas_mass
    growth_rate = min(growth_rate, eddington_bh_growth(M_BH))
    return np.asarray(growth_rate)


def seeding_mask_and_mass(halo_mass, redshift, gas_metallicity, already_seeded):
    """
    Decide which not-yet-seeded haloes cross a BH seeding threshold this
    step, and what mass to seed them with.

    Three independently enable-able channels (run_params.yaml:
    black_holes.seeding), checked in priority order pop3 ->
    direct_collapse -> halo_mass_threshold; the first channel satisfied for
    a given halo wins. Vectorised: arguments may be scalars or same-shaped
    arrays (scalars come back as 0-d arrays).

    Returns
    -------
    triggered : bool array -- haloes newly seeded this step
    seed_mass : float array -- the mass to seed them with (0 where not
        triggered)
    """
    halo_mass = np.asarray(halo_mass, dtype=float)
    redshift = np.asarray(redshift, dtype=float)
    gas_metallicity = np.asarray(gas_metallicity, dtype=float)
    already_seeded = np.asarray(already_seeded, dtype=bool)

    candidate = ~already_seeded
    triggered = np.zeros_like(halo_mass, dtype=bool)
    seed_mass = np.zeros_like(halo_mass, dtype=float)

    seeding = PARAMS.bh.seeding

    p3 = seeding.pop3
    if p3.enabled:
        mask = (
            candidate
            & (redshift >= p3.z_min)
            & (halo_mass >= p3.M_halo_min)
            & (gas_metallicity <= p3.Z_gas_max)
        )
        seed_mass = np.where(mask, p3.M_seed, seed_mass)
        triggered = triggered | mask

    dc = seeding.direct_collapse
    if dc.enabled:
        mask = (
            candidate
            & ~triggered
            & (redshift >= dc.z_min)
            & (halo_mass >= dc.M_halo_min)
            & (gas_metallicity <= dc.Z_gas_max)
        )
        seed_mass = np.where(mask, dc.M_seed, seed_mass)
        triggered = triggered | mask

    hm = seeding.halo_mass_threshold
    if hm.enabled:
        mask = candidate & ~triggered & (halo_mass >= hm.M_halo_min)
        seed_mass = np.where(mask, hm.M_seed, seed_mass)
        triggered = triggered | mask

    return triggered, seed_mass
