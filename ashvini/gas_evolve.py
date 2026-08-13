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
    bh_accretion_rate=0.0,
    agn_wind_growth_rate=0.0,
    bh_mass_for_wind=0.0,
):
    """
    Eqn 1 in Menon et al 2024 with 2 and 3 substituted, plus two BH-related
    terms:

    - bh_accretion_rate: mass flowing from the gas reservoir into the black
      hole itself (this step's own growth rate, always instantaneous --
      the accretion event and the reservoir it draws down happen at the
      same time by construction, regardless of whether the *wind* it later
      powers is delayed). This is a genuine gas-mass sink: BH growth was
      previously computed from gas_mass without ever depleting it, an
      inconsistency with the star-formation term below (which does deplete
      it) -- fixed here so BH growth and star formation compete for the
      same finite gas budget, as they physically must.
    - agn_wind_growth_rate: the BH growth rate the AGN wind term is
      computed from -- may equal bh_accretion_rate (no AGN delay
      configured, PARAMS.bh.feedback_delay_time=0) or be an earlier,
      delayed value (mirroring the SN delayed-feedback mechanism, see
      main.py's _delay_lookback_index), representing that the wind a
      given accretion episode powers isn't launched until some time later.
    - bh_mass_for_wind: the (step-midpoint-frozen) M_BH used only by
      agn.agn_wind_mass_rate's optional M-sigma self-regulation coupling
      (black_holes.sigma_feedback.enabled); ignored when that's disabled,
      the default. Frozen rather than evolving with `t` for the same
      reason agn_wind_growth_rate/halo_mass/stellar_metallicity already
      are -- see agn.agn_wind_mass_rate's docstring.
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
        - bh_accretion_rate
        - agn.agn_wind_mass_rate(
            agn_wind_growth_rate,
            bh_mass=bh_mass_for_wind,
            halo_mass=halo_mass,
            redshift=redshift,
        )
    )
    return np.asarray(gas_mass_evolution_rate)
