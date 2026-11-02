import numpy as np

from .ode import linear_ode_step
from ..dust import supernovae_per_msun
from ..metallicity import igm_metallicity


class MetalsAndDust:
    """
    Metals and dust, driven by the same wind star formation rate as the gas.

    Gas metals:   dZ_g/dt = Z_IGM*A_acc + y_Z*SFR_wind - Z_g*(B_sf + eta_w*SFR_wind/M_gas)
    Stellar metals are locked at B_sf*Z_g, so the same budget identity as for the stars applies.
    Dust:         dM_d/dt = y_d*SFR_wind - M_d*(SFR_wind/M_gas)*(nu_SN*M_swept*f_gas + eta_w)
    Dust is produced by the supernovae, destroyed by their blast waves (the dust in the swept-up gas,
    nu_SN*M_swept supernova-swept mass per unit star formation, gated by f_gas) and carried off in the wind.
    Z_g/M_gas are frozen at their value at the start of the step; empty gas reservoirs have no loss term.
    """

    def __init__(self, cfg):
        self.igm_metallicity = igm_metallicity(cfg)   # mass fraction (from Z_IGM in Z_sun)
        self.metal_yield = cfg.metals.Z_yield
        self.dust_yield = cfg.dust.dust_yield
        self.alpha = cfg.dust.dust_alpha
        self.m_crit = cfg.dust.m_crit
        self.m_swept = cfg.dust.m_swept
        self.sn_per_msun = supernovae_per_msun(cfg)

    def metals_step(self, gas_metals_prev, stars_metals_prev, accretion, B_sf, ML, wind_sfr,
                    gas_prev, has_gas, gas_prev_safe, dt):
        """Returns (gas_metals_raw, stars_metals_new), both before negative values are set to zero."""
        forcing = self.igm_metallicity * accretion + self.metal_yield * wind_sfr
        decay = B_sf + np.where(has_gas, ML * wind_sfr / gas_prev_safe, 0.0)
        gas_metals_raw = linear_ode_step(gas_metals_prev, forcing, decay, dt)
        locked = np.maximum(B_sf / decay * (gas_metals_prev + forcing * dt - gas_metals_raw), 0.0)
        return gas_metals_raw, stars_metals_prev + locked

    def dust_step(self, dust_prev, ML, wind_sfr, gas_prev, has_gas, gas_prev_safe, dt):
        loading = 1.0 - np.exp(-((gas_prev / self.m_crit) ** self.alpha))
        forcing = self.dust_yield * wind_sfr
        decay = np.where(
            has_gas, (self.sn_per_msun * self.m_swept * loading + ML) * wind_sfr / gas_prev_safe, 0.0
        )
        return linear_ode_step(dust_prev, forcing, decay, dt)
