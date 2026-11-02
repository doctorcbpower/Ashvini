import numpy as np

from .ode import linear_ode_step


def gas_and_stars_step(gas_prev, accretion, B_sf, ML, wind_sfr, feedback_type, bh_accretion, agn_wind, dt):
    """
    One step of the gas reservoir and the stars formed from it:

        dM_gas/dt = A_acc - B_sf*M_gas - eta_w*SFR_wind - dM_BH/dt - AGN wind

    Star formation B_sf*M_gas is always a loss proportional to the gas; for "instantaneous" feedback the
    wind (eta_w*B_sf*M_gas) is as well, otherwise the wind is a given forcing. bh_accretion and agn_wind
    (Msun/Gyr) are given forcing losses.

    Returns (gas_raw, stars_formed, gas_forcing, gas_decay). gas_raw can be slightly negative where the
    decay is fast relative to dt; the caller sets a negative value to zero. Integrating dy/dt = F - D*y gives
    D * integral(y dt) = y0 + F*dt - y1 exactly, so the stars formed are (B_sf / D) * (y0 + F*dt - y1):
    no first-order lag, and gas + stars + removed mass balance to round-off. A reservoir empty at the start
    of the step but filling during it forms stars in that step; where a negative-forcing step would drive
    the raw solution below zero, the stars formed are set to zero.
    """
    if feedback_type == "instantaneous":
        gas_forcing = accretion - bh_accretion - agn_wind
        gas_decay = B_sf * (1.0 + ML)
    else:
        gas_forcing = accretion - ML * wind_sfr - bh_accretion - agn_wind
        gas_decay = B_sf
    gas_raw = linear_ode_step(gas_prev, gas_forcing, gas_decay, dt)
    stars_formed = np.maximum(B_sf / gas_decay * (gas_prev + gas_forcing * dt - gas_raw), 0.0)
    return gas_raw, stars_formed, gas_forcing, gas_decay
