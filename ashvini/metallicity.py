from . import utils as utils
from . import supernovae_feedback as sn
from .config import resolve
from .constants import Z_SUN
from .star_formation import star_formation_rate, time_freefall


def igm_metallicity(cfg=None):
    """Metallicity of inflowing gas as an absolute mass fraction (the config value is in Z_sun)."""
    return resolve(cfg).metals.Z_IGM_Zsun * Z_SUN


def evolve_gas_metals(
    t,
    gas_metal_mass,
    gas_mass,
    gas_accretion_rate,
    halo_mass,
    stellar_metallicity,
    past_sfr,
    kind="delayed",
    cfg=None,
):
    """
    Gas-phase metal mass. gas_mass is the (frozen) gas mass used for the metal fraction in the
    wind-removal term; metals are locked into stars at the rate e_ff/t_ff * M_Z, which does not need
    the gas mass and so also applies to a reservoir that is empty at the start of the step.
    """
    cfg = resolve(cfg)
    redshift = utils.z_at_time(t)
    present_sfr = star_formation_rate(t, gas_mass=gas_mass, cfg=cfg)
    wind_sfr = past_sfr

    if kind == "no":
        wind_sfr = 0.0
    if kind == "instantaneous":
        wind_sfr = present_sfr

    gas_metals_rate = (
        (igm_metallicity(cfg) * gas_accretion_rate)  # Enriched gas accreting from IGM
        - (cfg.sf.efficiency / time_freefall(redshift)) * gas_metal_mass  # Removal from ISM during star formation
        + (cfg.metals.Z_yield * wind_sfr)  # Delayed enrichment of ISM by dying stars
    )
    if gas_mass > 0.0:
        gas_metals_rate -= (
            sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity, cfg)
            * (gas_metal_mass / gas_mass)
            * wind_sfr
        )  # Delayed removal from ISM by SN feedback

    return gas_metals_rate


def evolve_stars_metals(
    t,
    gas_metals,
    gas_mass=None,
    cfg=None,
):
    """Metals locked into stars: the gas-phase metal mass times the specific star formation rate e_ff/t_ff."""
    return (resolve(cfg).sf.efficiency / time_freefall(utils.z_at_time(t))) * gas_metals

from .config import retire_module_values as _retire  # noqa: E402

_retire(__name__, {"IGM_metallicity": "config.metals.Z_IGM_Zsun", "metallicity_yield": "config.metals.Z_yield"})
