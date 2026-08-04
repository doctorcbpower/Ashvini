# -*- coding: utf-8 -*-

import numpy as np
from . import utils as utils

from .run_params import PARAMS

epsilon_p = PARAMS.sn.epsilon_p
pi_fid = PARAMS.sn.pi_fid


def metallicity_function(stellar_metallicity, m=0.1, s=0.01, a=1, b=0.25):
    function_value = (
        np.exp(-np.logaddexp(0, -(stellar_metallicity - m) / s)) * (b - a) + a
    )
    return function_value


def mass_loading_factor(redshift, halo_mass, stellar_metallicity):
    mass_loading_factor = (
        epsilon_p
        * pi_fid
        * ((10**11.5 / halo_mass) ** (1 / 3))
        * ((9 / (1 + redshift)) ** (1 / 2))
        * metallicity_function(stellar_metallicity)
    )
    return mass_loading_factor


def wind_mass_evolution_rate(
    redshift, gas_mass, halo_mass, star_formation_rate_for_winds, stellar_metallicity
):
    wind_mass_rate = (
        metallicity_function(stellar_metallicity)
        * mass_loading_factor(redshift, halo_mass, stellar_metallicity)
        * star_formation_rate_for_winds
    )
    return wind_mass_rate
