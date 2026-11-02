"""
Ready-made parameter sets (Config overrides). A preset is a function returning a validated copy of a Config; read it
to see exactly which parameters it changes.
"""
import numpy as np

from .config import DEFAULT, override

# Stochastic_Occupation's star formation time is t_dyn = 1 / sqrt(G * 200 * rho_crit(z)) = sqrt(pi / 75) / H(z)
# (0.2047 / H); Ashvini's is 0.141 / H. SFR = eps_ff M_gas / time, so the same star formation rate needs
# Ashvini's efficiency = eps_ff * 0.141 / (sqrt(pi / 75) / H * H).
STOCHOCC_DYNAMICAL_TIME_IN_HUBBLE_TIMES = float(np.sqrt(np.pi / 75.0))
ASHVINI_FREEFALL_IN_HUBBLE_TIMES = 0.141


def ashvini_efficiency_for_stochocc(eps_ff):
    """The Ashvini star formation efficiency that gives the star formation rate of Stochastic_Occupation's eps_ff."""
    return eps_ff * ASHVINI_FREEFALL_IN_HUBBLE_TIMES / STOCHOCC_DYNAMICAL_TIME_IN_HUBBLE_TIMES


def stochocc_smooth_control(cfg=None, eps_ff=0.02178, f_p=0.7646, **changes):
    """
    The smooth limit of Stochastic_Occupation: supernovae from 8-100 Msun progenitors, each exploding after its own
    main-sequence lifetime (Georgy et al 2013, Z = 0.002), removing f_p * p_SN / v_esc of gas (momentum-driven wind),
    star formation at Stochastic_Occupation's eps_ff (0.02178 equals Ashvini's 0.015). f_p = 0.7646 is the value at
    which the momentum-driven wind equals Ashvini's calibrated mass loading at zero stellar metallicity
    (to 0.3 percent for halo masses 1e6-1e10 Msun and z = 5-15); Stochastic_Occupation's audit used 0.7625.
    """
    return override(cfg if cfg is not None else DEFAULT, **{
        "sf.efficiency": ashvini_efficiency_for_stochocc(eps_ff),
        "sn.type": "delayed", "sn.kernel": "imf", "sn.lifetime_Z": 0.002,
        "sn.m_sn_min": 8.0, "sn.m_sn_max": 100.0,
        "sn.closure": "momentum", "sn.f_p": f_p, "sn.p_SN_msun_kms": 3.0e5,
        **changes,
    })
