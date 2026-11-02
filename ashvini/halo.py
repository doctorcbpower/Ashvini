"""
Halo scales from the virial mass and redshift, for the supernova closures that depend on the escape speed.
Everything is a function of (M_h, z) only. CGS; halo mass in Msun. The virial radius is that of a sphere of
mean density Delta (= 200) times the critical density at z; the escape speed is sqrt(2) times the circular
speed at that radius (v_esc^2 = 2 V_c^2, a convention, not a derivation).
"""
import numpy as np

from . import utils
from .constants import G, Msun_g

MPC_CM = 3.0857e24


def hubble_s(z):
    """H(z) in 1/s."""
    return utils.h_of_z(z) * 100.0 * 1.0e5 / MPC_CM


def rho_crit(z):
    """Critical density, g/cm^3."""
    return 3.0 * hubble_s(z) ** 2 / (8.0 * np.pi * G)


def r_vir_cm(halo_mass, z, Delta=200.0):
    return (3.0 * np.asarray(halo_mass, dtype=float) * Msun_g / (4.0 * np.pi * Delta * rho_crit(z))) ** (1.0 / 3.0)


def v_circ_cms(halo_mass, z, Delta=200.0):
    m = np.asarray(halo_mass, dtype=float)
    return np.sqrt(G * m * Msun_g / r_vir_cm(m, z, Delta))


def v_esc_cms(halo_mass, z, Delta=200.0):
    return np.sqrt(2.0) * v_circ_cms(halo_mass, z, Delta)
