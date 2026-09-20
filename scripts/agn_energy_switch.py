"""
ILLUSTRATIVE AGN wind for the general model (main.run_forest): momentum-driven below King's M_sigma, energy-
driven above it, blended by the existing logistic switch in M_BH / M_sigma, with BH growth capped at the
Power et al. (2011) ceiling. Installed by monkeypatching, so the package itself is unchanged.

    Mdot_wind = [ (1 - w) * L_mom + w * L_en ] * Mdot_BH,      w = agn_feedback.coupling_switch(M_BH, M_sigma)
    L_mom = f_mom * eps/(1-eps) * c / sigma                    (as in the MVM, reservoir_stock.momentum_agn_wind_rate)
    L_en  = 2 eps_f c^2 / sigma^2                              (King 2003, black_holes_growth_slimdisk.king_agn_wind_rate)

The shipped default, Mdot_wind = 0.5 Mdot_BH, is 2-4 orders of magnitude below both (docs/UV_SUPPRESSION_AUDIT.md
neighbour note in the AGN discussion). eps_f = 5e-4 is the value in run_params.yaml's slimdisk block.
"""
import contextlib

import numpy as np

from ashvini import agn_feedback as agn, black_holes_growth as bg, black_holes_growth_slimdisk as sd
from ashvini.paper_reservoir_params import PAPER_PARAMS
from ashvini.reservoir_stock import momentum_agn_wind_rate


def make_wind(epsilon_f=5e-4, f_mom=1.0, eps_rad=None, energy=True, momentum=True):
    eps_rad = PAPER_PARAMS.epsilon if eps_rad is None else eps_rad

    def wind(bh_growth_rate, bh_mass=None, halo_mass=None, redshift=None):
        # With the growth cap on, M_BH can fall when the ceiling drops between steps; a wind must not follow a
        # negative rate (it would create gas), so it is driven by accretion only.
        rate = np.maximum(np.asarray(bh_growth_rate, dtype=float), 0.0)
        sigma = bg.velocity_dispersion(halo_mass, redshift)
        w = agn.coupling_switch(bh_mass, bg.m_sigma(sigma))
        mom = momentum_agn_wind_rate(rate, halo_mass, redshift, eps_rad, f_mom) if momentum else 0.0
        en = sd.king_agn_wind_rate(rate, halo_mass, redshift, epsilon_f=epsilon_f) if energy else 0.0
        return (1.0 - w) * mom + w * en

    return wind


@contextlib.contextmanager
def installed(epsilon_f=5e-4, f_mom=1.0, cap=True, energy=True, momentum=True):
    saved = (agn.agn_wind_mass_rate, bg.growth_cap_enabled)
    agn.agn_wind_mass_rate = make_wind(epsilon_f, f_mom, energy=energy, momentum=momentum)
    bg.growth_cap_enabled = cap
    try:
        yield
    finally:
        agn.agn_wind_mass_rate, bg.growth_cap_enabled = saved
