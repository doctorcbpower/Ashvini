"""
A disc size for the gas, and a star formation threshold that depends on it (Config section `disc`; off by default).

The gas that has settled in the galaxy is taken to form an exponential disc, surface density
Sigma(r) = Sigma_0 exp(-r / R_d), holding the whole gas mass M_gas:  Sigma_0 = M_gas / (2 pi R_d^2).

Scale length (Mo, Mao & White 1998, without the corrections for the halo concentration and for the disc's pull on the
halo):  R_d = (j_d / m_d) * lambda * R_vir / sqrt(2),  with lambda the halo spin parameter, R_vir the virial radius
(a mean density 200 times critical) and j_d / m_d the disc's share of the halo's angular momentum divided by its share
of the mass. For j_d / m_d = 1 this is the definition used by the reservoir model (reservoir_stock.py).

Star formation threshold: gas forms stars only where Sigma(r) is above a critical surface density Sigma_crit. For an
exponential disc that is every radius inside R_crit = R_d ln(Sigma_0 / Sigma_crit), and that radius holds the fraction
f = 1 - (1 + x) exp(-x) of the gas, with x = ln(Sigma_0 / Sigma_crit) = R_crit / R_d. The star formation rate is the
usual rate applied to that part of the gas only: SFR = f * eps_ff / t_ff * M_gas. If Sigma_0 <= Sigma_crit no gas forms
stars (f = 0); with Sigma_crit = 0 all of it does (f = 1). f changes smoothly as the disc grows, so star formation does
not switch on and off from step to step.
"""
import numpy as np

from .. import halo
from ..components.star_formation import FreeFallStarFormation

KPC_CM = 3.0857e21


def scale_length_kpc(halo_mass, redshift, cfg):
    """Disc scale length R_d in kpc; 0 for a halo with no mass."""
    d = cfg.disc
    halo_mass = np.asarray(halo_mass, dtype=float)
    r_vir_kpc = halo.r_vir_cm(np.where(halo_mass > 0, halo_mass, 1.0), redshift) / KPC_CM
    return np.where(halo_mass > 0, d.angular_momentum_ratio * d.spin_parameter * r_vir_kpc / np.sqrt(2.0), 0.0)


def central_surface_density(gas_mass, scale_length_kpc_):
    """Sigma_0 in Msun/pc^2 for gas_mass (Msun) in an exponential disc of scale length R_d (kpc); 0 if R_d = 0."""
    R_pc = np.asarray(scale_length_kpc_, dtype=float) * 1.0e3
    safe = np.where(R_pc > 0, R_pc, 1.0)
    return np.where(R_pc > 0, np.asarray(gas_mass, dtype=float) / (2.0 * np.pi * safe ** 2), 0.0)


def star_forming_fraction(central_sigma, critical_sigma):
    """Fraction of the gas of an exponential disc whose surface density is above critical_sigma (both in the same units)."""
    central_sigma = np.asarray(central_sigma, dtype=float)
    if critical_sigma <= 0.0:
        return np.ones_like(central_sigma)
    x = np.log(np.maximum(central_sigma / critical_sigma, 1.0))      # 0 when the centre is below the threshold
    return 1.0 - (1.0 + x) * np.exp(-x)


def diagnostics(halo_mass, gas_mass, redshift, cfg):
    """R_d (kpc), Sigma_0 (Msun/pc^2), R_crit/R_d and the star-forming fraction f, for arrays of shape (haloes, steps)."""
    R_d = scale_length_kpc(halo_mass, np.asarray(redshift), cfg)
    sigma0 = central_surface_density(gas_mass, R_d)
    crit = cfg.disc.critical_surface_density_Msun_pc2
    f = star_forming_fraction(sigma0, crit)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(sigma0 > 0, np.log(np.maximum(sigma0 / crit, 1.0)), 0.0) if crit > 0 else np.full_like(sigma0, np.inf)
    return {"scale_length_kpc": R_d, "central_surface_density": sigma0, "radius_ratio": x, "star_forming_fraction": f}


class DiscStarFormation(FreeFallStarFormation):
    """The free-fall star formation law, applied to the gas above the critical surface density only."""

    def __init__(self, cfg):
        super().__init__(cfg)
        self.cfg = cfg

    def halo_rate_coefficient(self, redshift, halo_mass, gas_mass):
        R_d = scale_length_kpc(halo_mass, redshift, self.cfg)
        f = star_forming_fraction(central_surface_density(gas_mass, R_d), self.cfg.disc.critical_surface_density_Msun_pc2)
        return self.rate_coefficient(redshift) * f
