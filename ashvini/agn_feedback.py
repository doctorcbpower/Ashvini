# -*- coding: utf-8 -*-

import numpy as np

from .config import resolve
from . import black_holes_growth as bh_growth


def coupling_switch(bh_mass, m_sigma_mass, cfg=None):
    """
    Thrust-vs-weight coupling switch for the AGN wind, King (2003, 2005):
    AGN wind coupling to the ISM is weak (momentum-driven, effectively
    trapped near the BH) while M_BH < M_sigma, and becomes effective at
    expelling gas once M_BH >= M_sigma (see
    black_holes_growth.velocity_dispersion/m_sigma for M_sigma itself).

    Neither King (2003, 2005) nor Power, Zubovas, Nayakshin & King (2011)
    specify a particular smooth interpolation between the two regimes --
    only that a regime change happens near M_BH ~ M_sigma -- so the exact
    functional form here is a deliberate modelling choice, not something
    read off either paper. It's isolated in this one function precisely so
    it can be revisited independently of the rest of the pipeline: a
    logistic function of log10(M_BH/M_sigma), centered at M_BH = M_sigma,
    with dex width `transition_width`
    (black_holes.sigma_feedback.transition_width) -- 0 well below
    M_sigma, 1 well above it.

    Returns 0 where bh_mass <= 0 or m_sigma_mass <= 0 (no BH yet, or an
    unformed halo with M_sigma = 0) rather than evaluating the ratio.
    """
    transition_width = resolve(cfg).bh.sigma_feedback.transition_width  # dex, width of the M_BH/M_sigma switch
    bh_mass = np.asarray(bh_mass, dtype=float)
    m_sigma_mass = np.asarray(m_sigma_mass, dtype=float)
    valid = (bh_mass > 0) & (m_sigma_mass > 0)

    bh_mass_safe = np.where(valid, bh_mass, 1.0)
    m_sigma_safe = np.where(valid, m_sigma_mass, 1.0)
    log_ratio = np.log10(bh_mass_safe / m_sigma_safe)

    # Numerically stable logistic, same pattern as
    # supernovae_feedback.metallicity_function.
    switch = np.exp(-np.logaddexp(0, -log_ratio / transition_width))
    return np.where(valid, switch, 0.0)


def agn_wind_mass_rate(bh_growth_rate, bh_mass=None, halo_mass=None, redshift=None, cfg=None):
    """
    Mass-loaded AGN wind, mirroring supernovae_feedback's mass-loading
    structure: gas ejected from the ISM proportional to the BH accretion
    rate, Mdot_wind = eta_agn_eff * Mdot_BH.

    If sigma_feedback_enabled (black_holes.sigma_feedback.enabled) is
    True, eta_agn_eff replaces the constant eta_agn with
    eta_agn * coupling_switch(bh_mass, M_sigma(halo_mass, redshift)) --
    bh_mass, halo_mass, redshift are then required, and should be the
    step-midpoint-frozen values used for every other time-varying ODE
    coefficient (see main.py's run_forest/run1_scalar), not the ODE's own
    continuously-evolving state -- M_BH varies within the same step this
    term is being evaluated for, and letting eta_agn_eff depend on it at
    the step's own (not-yet-known) end would make the gas-mass ODE
    genuinely non-linear within the step, breaking the closed-form
    (_linear_ode_step) solver run_forest() relies on.

    If disabled (the default), eta_agn_eff = eta_agn identically and
    bh_mass/halo_mass/redshift are unused -- this preserves the exact
    prior constant-coupling behaviour.
    """
    bh = resolve(cfg).bh
    if not bh.sigma_feedback.enabled:
        return bh.eta_agn * bh_growth_rate

    sigma = bh_growth.velocity_dispersion(halo_mass, redshift)
    m_sigma_mass = bh_growth.m_sigma(sigma, cfg)
    eta_agn_eff = bh.eta_agn * coupling_switch(bh_mass, m_sigma_mass, cfg)
    return eta_agn_eff * bh_growth_rate

from .config import retire_module_values as _retire  # noqa: E402

_retire(__name__, {
    "eta_agn": "config.bh.eta_agn", "sigma_feedback_enabled": "config.bh.sigma_feedback.enabled",
    "transition_width": "config.bh.sigma_feedback.transition_width",
})
