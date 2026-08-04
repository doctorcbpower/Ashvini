import numpy as np

from .utils import Hubble_time, z_at_time

from .run_params import PARAMS

e_bh = PARAMS.bh.efficiency  # efficiency of black hole growth
c = 3.0e10  # speed of light in CGS (cm/s)
G = 6.674e-8  # gravitational constant in CGS
m_p = 1.67e-24  # proton mass in CGS
sigma_thomson = 6.65e-25  # Thomson cross section in CGS

def time_freefall(redshift):
    return 0.141 * Hubble_time(redshift)

def black_hole_growth_rate(
    t,  # cosmic time
    gas_mass,
    M_BH
):
    redshift = z_at_time(t)

    black_hole_growth_rate = (e_bh / time_freefall(redshift)) * gas_mass
    
    black_hole_growth_rate = min(black_hole_growth_rate,eddington_bh_growth(M_BH))
    
    return np.asarray(black_hole_growth_rate)

def eddington_bh_growth(M_BH):
    return 4 * np.pi * G * M_BH * m_p/sigma_thomson/c 
