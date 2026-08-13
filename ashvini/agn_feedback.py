# -*- coding: utf-8 -*-

from .run_params import PARAMS

eta_agn = PARAMS.bh.eta_agn  # AGN wind mass-loading efficiency


def agn_wind_mass_rate(bh_growth_rate):
    """
    Mass-loaded AGN wind, mirroring supernovae_feedback's mass-loading
    structure: gas ejected from the ISM proportional to the BH accretion
    rate, Mdot_wind = eta_agn * Mdot_BH.
    """
    return eta_agn * bh_growth_rate
