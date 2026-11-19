import functools

import numpy as np

from . import supernovae_feedback as sn
from . import utils as utils
from .config import resolve
from .imf import nu_SN


@functools.lru_cache(maxsize=None)
def _nu_SN(imf, m_sn_min, m_sn_max):
    return nu_SN(imf, m_sn_min, m_sn_max)


def supernovae_per_msun(cfg=None):
    """Core-collapse supernovae per Msun of stars formed, from cfg.sn.imf, m_sn_min and m_sn_max."""
    s = resolve(cfg).sn
    return _nu_SN(s.imf, s.m_sn_min, s.m_sn_max)


def update_dust_reservoir(
    t,
    dust_mass,
    gas_mass,
    halo_mass,
    wind_sfr,
    stellar_metallicity,
    cfg=None,
):
    """
    Dust mass evolution, driven by wind_sfr: the star formation rate that drives supernova feedback
    (zero for supernova.type "no", the present rate for "instantaneous", the delayed rate for
    "delayed"), so dust, metals and gas all follow the same feedback switch.

        dM_d/dt = Y_d * SFR_wind                          (production by dying stars)
                  - nu_SN * SFR_wind * M_swept * (M_d/M_gas) * f_gas   (destruction in supernova shocks)
                  - (M_d/M_gas) * eta_w * SFR_wind        (removal by the supernova wind)

    nu_SN is the number of supernovae per Msun of stars formed (imf.nu_SN); the swept-up gas mass M_swept
    carries the dust-to-gas ratio M_d/M_gas, so the dust destroyed per supernova is M_swept * M_d/M_gas, not a
    fixed mass; f_gas = 1 - exp(-(M_gas/M_crit)^alpha) switches destruction off in gas-poor haloes.
    """
    cfg = resolve(cfg)
    d = cfg.dust
    growth_rate = d.dust_yield * wind_sfr

    if gas_mass <= 0:
        return np.asarray(growth_rate)

    dust_loading = 1 - np.exp(-((gas_mass / d.m_crit) ** d.dust_alpha))
    redshift = utils.z_at_time(t)
    ml = sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity, cfg)
    loss_per_dust = (supernovae_per_msun(cfg) * d.m_swept * dust_loading + ml) * wind_sfr / gas_mass

    return np.asarray(growth_rate - loss_per_dust * dust_mass)
