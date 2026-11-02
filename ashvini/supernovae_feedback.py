# -*- coding: utf-8 -*-

import numpy as np
from . import utils as utils

from .run_params import PARAMS
from .constants import Z_SUN

epsilon_p = PARAMS.sn.epsilon_p
pi_fid = PARAMS.sn.pi_fid


def metallicity_function(stellar_metallicity, m=0.1, s=0.01, a=1, b=0.25):
    """
    Metallicity dependence of the mass loading: a smooth step from a (low metallicity) to b (high
    metallicity), centred on m with width s. stellar_metallicity is an absolute mass fraction; m and s
    are in units of Z_sun (constants.Z_SUN), so the step is centred on 0.1 Z_sun.
    """
    x = np.asarray(stellar_metallicity, dtype=float) / Z_SUN
    function_value = (
        np.exp(-np.logaddexp(0, -(x - m) / s)) * (b - a) + a
    )
    return function_value


def mass_loading_factor(redshift, halo_mass, stellar_metallicity):
    # A halo_mass of exactly 0 (e.g. a not-yet-formed progenitor -- see
    # build_trees_from_pymctrees.py's pre-formation zeroing) has no wind to
    # load: the caller (wind_mass_evolution_rate) multiplies this by
    # star_formation_rate_for_winds, which is already 0 there, so the
    # correct answer is 0 regardless. But 10**11.5/halo_mass diverges to
    # inf as halo_mass->0, and inf * (an already-correct 0) is nan, not 0
    # -- guarded here directly rather than relying on that multiplication.
    halo_mass = np.asarray(halo_mass, dtype=float)
    resolved = halo_mass > 0
    halo_mass_safe = np.where(resolved, halo_mass, 1.0)

    mass_loading_factor = (
        epsilon_p
        * pi_fid
        * ((10**11.5 / halo_mass_safe) ** (1 / 3))
        * ((9 / (1 + redshift)) ** (1 / 2))
        * metallicity_function(stellar_metallicity)
    )
    return np.where(resolved, mass_loading_factor, 0.0)


def wind_mass_evolution_rate(
    redshift, gas_mass, halo_mass, star_formation_rate_for_winds, stellar_metallicity
):
    wind_mass_rate = (
        metallicity_function(stellar_metallicity)
        * mass_loading_factor(redshift, halo_mass, stellar_metallicity)
        * star_formation_rate_for_winds
    )
    return wind_mass_rate
