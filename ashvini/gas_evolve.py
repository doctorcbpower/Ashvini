import numpy as np

from . import reionization as reion
from . import supernovae_feedback as sn
from . import agn_feedback as agn
from . import utils as utils

from .utils import Omega_b, Omega_m
from .star_formation import star_formation_rate


def gas_inflow_rate(redshift, halo_mass, halo_mass_dot, UV_background=True):
    """
    Cosmological baryonic accretion rate, modulated by UV suppression.

    Args:
        redshift (float): Redshift.
        halo_mass (float): Halo mass.
        halo_mass_dot (float): Halo mass accretion rate.
        uv_suppresion (bool ?): UVB on or off.

    Returns:
        Float: The baryonic cosmological accretion rate.
    """
    gas_accretion_rate = (Omega_b / Omega_m) * halo_mass_dot
    if UV_background:
        gas_accretion_rate *= reion.uv_suppression(redshift, halo_mass, halo_mass_dot)

    return np.asarray(gas_accretion_rate)


def update_gas_reservoir(
    t,
    gas_mass,
    gas_accretion_rate,
    halo_mass,
    stellar_metallicity,
    past_sfr,
    kind="delayed",
    agn_growth_rate=0.0,
):
    """
    Eqn 1 in Menon et al 2024 with 2 and 3 substituted, plus an AGN wind
    term proportional to the BH accretion rate (agn_growth_rate).
    """

    redshift = utils.z_at_time(t)
    present_sfr = star_formation_rate(t, gas_mass=gas_mass)
    wind_sfr = past_sfr

    if kind == "no":
        wind_sfr = 0
    if kind == "instantaneous":
        wind_sfr = present_sfr

    gas_mass_evolution_rate = (
        gas_accretion_rate
        - present_sfr
        - sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity) * wind_sfr
        - agn.agn_wind_mass_rate(agn_growth_rate)
    )
    return np.asarray(gas_mass_evolution_rate)
