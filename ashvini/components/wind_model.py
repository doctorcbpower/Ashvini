"""
How much gas the supernovae remove, per unit mass of stars formed (the mass loading eta_w), as a function of
redshift, halo mass and stellar metallicity. Three wind models; each has `mass_loading(z, halo_mass,
stellar_metallicity)` returning eta_w (shape of halo_mass). Choose with `supernova.wind_model`.

  ashvini    eta_w = epsilon_p * pi_fid * (1e11.5/M_h)^(1/3) * (9/(1+z))^(1/2) * f(Z_star): the calibrated
             prescription of Menon et al. It follows V_esc^-1, which is a momentum-conserving wind
             (it equals the momentum-driven wind below with f_p = 0.7625 to 0.2 percent for 1e6-1e10 Msun,
             z = 5-15, p_SN = 3e5 Msun km/s, 8-100 Msun progenitors).
  momentum   each supernova removes f_p * p_SN / v_esc;       eta_w = nu_SN * f_p * p_SN / v_esc
  energy     each supernova removes 2 eps_fb E_SN / v_esc^2;  eta_w = nu_SN * 2 eps_fb E_SN / v_esc^2

nu_SN is the number of supernova progenitors per Msun of stars formed (imf.nu_SN for the configured IMF and
progenitor range). The momentum and energy wind models give the mass an average supernova would remove if the
reservoir held enough gas; the limit at the gas actually present is applied by the time integrator that has it
(per event in the event integrator; not at all in the grid integrator, which sets a negative gas mass to zero afterwards).
"""
import numpy as np

from .. import halo
from .. import supernovae_feedback as sn
from ..constants import Msun_g
from ..imf import nu_SN


class AshviniLoading:
    def __init__(self, cfg):
        self.cfg = cfg
        self.nu_SN = nu_SN(cfg.sn.imf, cfg.sn.m_sn_min, cfg.sn.m_sn_max)

    def per_supernova(self, redshift, halo_mass, stellar_metallicity):
        """Msun removed by one supernova, for the event integrator: the mass loading divided by nu_SN."""
        return self.mass_loading(redshift, halo_mass, stellar_metallicity) / self.nu_SN

    def mass_loading(self, redshift, halo_mass, stellar_metallicity):
        return sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity, self.cfg)


class _SupernovaWind:
    def __init__(self, cfg):
        self.cfg = cfg
        self.nu_SN = nu_SN(cfg.sn.imf, cfg.sn.m_sn_min, cfg.sn.m_sn_max)

    def mass_per_supernova(self, redshift, halo_mass):
        """Msun removed by one supernova (before the limit at the gas present); 0 for an empty halo."""
        halo_mass = np.asarray(halo_mass, dtype=float)
        resolved = halo_mass > 0
        v_esc = halo.v_esc_cms(np.where(resolved, halo_mass, 1.0), redshift)
        return np.where(resolved, self._per_supernova(v_esc), 0.0)

    def mass_loading(self, redshift, halo_mass, stellar_metallicity=None):
        return self.nu_SN * self.mass_per_supernova(redshift, halo_mass)

    def per_supernova(self, redshift, halo_mass, stellar_metallicity=None):
        return self.mass_per_supernova(redshift, halo_mass)


class MomentumDrivenWind(_SupernovaWind):
    def _per_supernova(self, v_esc):
        sn_ = self.cfg.sn
        return sn_.f_p * sn_.p_SN_msun_kms / (v_esc / 1.0e5)


class EnergyDrivenWind(_SupernovaWind):
    def _per_supernova(self, v_esc, eps_fb=None):
        sn_ = self.cfg.sn
        return 2.0 * (sn_.eps_fb if eps_fb is None else eps_fb) * sn_.E_SN_erg / v_esc ** 2 / Msun_g

    def per_supernova(self, redshift, halo_mass, stellar_metallicity=None, eps_fb=None):
        """As the other wind models, with the coupling eps_fb optionally given for this call (event integrator's coupling_change_mass)."""
        halo_mass = np.asarray(halo_mass, dtype=float)
        resolved = halo_mass > 0
        v_esc = halo.v_esc_cms(np.where(resolved, halo_mass, 1.0), redshift)
        return np.where(resolved, self._per_supernova(v_esc, eps_fb), 0.0)


def make_wind_model(cfg):
    """The wind model selected by cfg.sn.wind_model."""
    return {"ashvini": AshviniLoading, "momentum": MomentumDrivenWind, "energy": EnergyDrivenWind}[cfg.sn.wind_model](cfg)
