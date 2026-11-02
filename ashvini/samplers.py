"""
Where and when the supernovae of a newly formed stellar population go off, for the event stepper.

A sampler turns "stars_formed Msun appeared at time t_form (Gyr)" into a list of Explosions: each has a time, a
progenitor mass (Msun) and a weight (1 for a whole supernova). Choose with `samplers.make_sampler(cfg, kind)`.

  StochasticSampler  A Poisson number of whole supernovae, nu_SN * stars_formed on average, with progenitor masses
                     drawn from the IMF and each exploding a main-sequence lifetime after formation. The
                     integer counts and the random masses are the stochastic part.
  SmoothSampler      The same expected supernovae spread, with fractional weights, over K fixed representative
                     masses (equal-probability quantiles of the IMF): the expected number of each mass, none
                     of the randomness. The star formation of each 1 Myr is released together, at its
                     mass-weighted mean time, so the number of events stays manageable. This is the control
                     for the stochastic runs: the 'continuous IMF' limit with the same lifetimes and closure.
"""
from dataclasses import dataclass

import numpy as np

from . import imf
from .stellar_lifetimes import t_ms_myr

MYR_GYR = 1.0e-3


@dataclass
class Explosions:
    t_gyr: np.ndarray       # explosion times, ascending
    mass_msun: np.ndarray   # progenitor masses
    weight: np.ndarray      # fraction of a whole supernova each represents


_EMPTY = Explosions(np.empty(0), np.empty(0), np.empty(0))


class StochasticSampler:
    def __init__(self, cfg, rng):
        sn = cfg.sn
        self.nu_SN = imf.nu_SN(sn.imf, sn.m_sn_min, sn.m_sn_max)
        self.m_min, self.m_max, self.slope = sn.m_sn_min, sn.m_sn_max, imf.progenitor_slope(sn.imf)
        self.lifetime_Z = sn.lifetime_Z
        self.rng = rng

    def sample(self, t_form, stars_formed):
        n = self.rng.poisson(self.nu_SN * stars_formed)
        if n == 0:
            return _EMPTY
        m = imf.sample_progenitor_masses(n, self.rng, self.m_min, self.m_max, self.slope)
        t = t_form + t_ms_myr(m, self.lifetime_Z) * MYR_GYR
        order = np.argsort(t)
        return Explosions(t[order], m[order], np.ones(n))


class SmoothSampler:
    def __init__(self, cfg, K=100, dt_emit_myr=1.0):
        sn = cfg.sn
        self.nu_SN = imf.nu_SN(sn.imf, sn.m_sn_min, sn.m_sn_max)
        self.K = K
        self.dt_emit = dt_emit_myr * MYR_GYR
        masses = imf.quantile_progenitor_masses(K, sn.m_sn_min, sn.m_sn_max, imf.progenitor_slope(sn.imf))
        lifetimes = t_ms_myr(masses, sn.lifetime_Z) * MYR_GYR
        order = np.argsort(lifetimes)
        self.masses, self.lifetimes = masses[order], lifetimes[order]
        self._mass = 0.0       # stars formed and not yet released
        self._mass_time = 0.0  # sum of mass * formation time
        self._t_last = -np.inf

    def sample(self, t_form, stars_formed):
        self._mass += stars_formed
        self._mass_time += stars_formed * t_form
        if t_form - self._t_last < self.dt_emit or self._mass == 0.0:
            return _EMPTY
        t_eff = self._mass_time / self._mass
        w = np.full(self.K, self.nu_SN * self._mass / self.K)
        self._mass = self._mass_time = 0.0
        self._t_last = t_form
        return Explosions(t_eff + self.lifetimes, self.masses, w)


def make_sampler(cfg, kind="smooth", rng=None, **kwargs):
    """kind: "smooth" (expected supernovae, no randomness) or "stochastic" (needs rng, a numpy Generator)."""
    if kind == "smooth":
        return SmoothSampler(cfg, **kwargs)
    if kind == "stochastic":
        return StochasticSampler(cfg, rng if rng is not None else np.random.default_rng())
    raise ValueError(f"sampler kind must be 'smooth' or 'stochastic'; got {kind!r}.")
