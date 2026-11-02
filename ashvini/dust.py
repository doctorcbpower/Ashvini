import numpy as np

from . import supernovae_feedback as sn
from . import utils as utils
from .imf import nu_SN

from .run_params import PARAMS

Y_d = PARAMS.dust.dust_yield
Alpha = PARAMS.dust.dust_alpha
M_crit = PARAMS.dust.m_crit
M_swept = PARAMS.dust.m_swept
# Core-collapse supernovae per Msun of stars formed, from the IMF and progenitor range in the config.
SN_per_Msun = nu_SN(PARAMS.sn.imf, PARAMS.sn.m_sn_min, PARAMS.sn.m_sn_max)


def update_dust_reservoir(
    t,
    dust_mass,
    gas_mass,
    halo_mass,
    wind_sfr,
    stellar_metallicity,
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
    growth_rate = Y_d * wind_sfr

    if gas_mass <= 0:
        return np.asarray(growth_rate)

    dust_loading = 1 - np.exp(-((gas_mass / M_crit) ** Alpha))
    redshift = utils.z_at_time(t)
    ml = sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity)
    loss_per_dust = (SN_per_Msun * M_swept * dust_loading + ml) * wind_sfr / gas_mass

    return np.asarray(growth_rate - loss_per_dust * dust_mass)
