import numpy as np

from .utils import Hubble_time, z_at_time

from .run_params import PARAMS

e_bh = PARAMS.bh.efficiency  # efficiency of black hole growth
eddington_multiplier = PARAMS.bh.eddington_multiplier  # f_Edd; >1 allows super-Eddington growth
c = 3.0e10  # speed of light in CGS (cm/s)
G = 6.674e-8  # gravitational constant in CGS
m_p = 1.67e-24  # proton mass in CGS
sigma_thomson = 6.65e-25  # Thomson cross section in CGS

Gyr_s = 3.15576e16  # seconds per Gyr
Msun_g = 1.98892e33  # grams per solar mass (IAU nominal)

# black_holes.sigma_feedback -- isothermal-sphere M-sigma self-regulation
# (King 2003, 2005; Power, Zubovas, Nayakshin & King 2011). Only used by
# velocity_dispersion()/m_sigma() below, and only when
# ashvini.agn_feedback's sigma_feedback_enabled is True.
_sigma_feedback_params = PARAMS.bh.sigma_feedback
f_g = _sigma_feedback_params.f_g  # baryon fraction relative to dark matter
# Electron-scattering opacity, cm^2/g. This is a genuinely different
# quantity from the `kappa_per_s`/`kappa_per_gyr` constants below (the
# mass-independent Eddington accretion rate per unit BH mass) -- both
# happen to be conventionally called "kappa" in the AGN feedback
# literature, but they enter completely different formulas (Eddington rate
# vs. M-sigma normalisation). Configurable via
# black_holes.sigma_feedback.kappa_es; None (the default) falls back to
# the physical value sigma_thomson/m_p.
kappa_es = (
    _sigma_feedback_params.kappa_es
    if _sigma_feedback_params.kappa_es is not None
    else sigma_thomson / m_p
)


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


EDDINGTON_RATE_PER_UNIT_MASS = float(eddington_bh_growth(1.0))  # 1/Gyr, f_Edd=1 (physical Eddington rate)


def velocity_dispersion(halo_mass, redshift):
    """
    Isothermal-sphere velocity dispersion sigma(M_halo, z), in cm/s (CGS).

    Derived from the isothermal-sphere virial relations
    R_v = sigma / (5*sqrt(2)*H(z)) and M_v = 2*sigma^2*R_v/G (Power,
    Zubovas, Nayakshin & King 2011, eq. 18):

        sigma(M_halo, z) = [ (5*sqrt(2)/2) * G * H(z) * M_halo ]^(1/3)

    H(z) is taken from Hubble_time(z) (= 1/H(z)), i.e. the SAME (Planck18)
    cosmology already used everywhere else in this package for cosmic
    time/redshift interpolation -- deliberately not a third cosmology
    alongside the pre-existing utils.py (Planck18) / reionization.py
    (Planck15) split documented in MODELS.md.

    halo_mass <= 0 (unformed progenitor) returns sigma = 0, the same guard
    pattern used elsewhere (supernovae_feedback.mass_loading_factor,
    reionization.uv_suppression) for undefined halo properties.
    """
    halo_mass = np.asarray(halo_mass, dtype=float)
    resolved = halo_mass > 0
    halo_mass_g = np.where(resolved, halo_mass, 1.0) * Msun_g

    H_per_s = (1.0 / Hubble_time(redshift)) / Gyr_s  # 1/Gyr -> 1/s

    sigma = (2.5 * np.sqrt(2) * G * H_per_s * halo_mass_g) ** (1.0 / 3.0)
    return np.where(resolved, sigma, 0.0)


def m_sigma(sigma):
    """
    King (2003, 2005) self-regulated BH mass scale (Power, Zubovas,
    Nayakshin & King 2011, eq. 5), where the wind's Eddington-limited
    momentum thrust balances the weight of the overlying isothermal-sphere
    gas:

        M_sigma(sigma) = (f_g * kappa_es / (pi * G^2)) * sigma^4

    sigma in cm/s (see velocity_dispersion above). Returns M_sigma in
    Msun. sigma = 0 (unformed halo) gives M_sigma = 0.
    """
    sigma = np.asarray(sigma, dtype=float)
    M_sigma_g = (f_g * kappa_es / (np.pi * G**2)) * sigma**4
    return M_sigma_g / Msun_g


def black_hole_growth_rate(t, M_BH, gas_mass):
    """
    dM_BH/dt (Msun/Gyr): gas-supply-limited growth capped at
    eddington_multiplier x the Eddington rate (eddington_multiplier=1 by
    default; >1 allows super-Eddington growth, per PARAMS.bh.
    eddington_multiplier). M_BH (the ODE state) is the first argument
    after t, matching the calling convention the other RHS functions in
    this package (update_gas_reservoir, evolve_gas_metals, ...) already
    use for solve_ivp.
    """
    redshift = z_at_time(t)
    growth_rate = (e_bh / time_freefall(redshift)) * gas_mass
    growth_rate = min(growth_rate, eddington_multiplier * eddington_bh_growth(M_BH))
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
