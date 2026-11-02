"""
UV-background suppression of gas accretion onto haloes (Okamoto, Gnedin & Frenk 2008).

reionization_physics.py contains the pure equilibrium prescription (M_c, mu_c, s, beta, c_omega,
each a function of its explicit physical arguments only); this module contains Ashvini's
configuration-bound wrapper around that prescription (z_rei, gamma and omega are read from the Config
passed as `cfg`, or from config.DEFAULT; beta and c_omega are derived from them on each call) plus Ashvini's own rate-dependent correction on top of it (X, epsilon, and uv_suppression's
mdot_halo/H(z) terms) -- a genuinely different, non-equilibrium physical model, not just a different
implementation of the same one.
"""
import numpy as np
import astropy.units as u

from astropy.cosmology import Planck15 as cosmo

from .config import resolve
from . import reionization_physics as _physics

H_0 = cosmo.H0  # in km / (Mpc s)
H_0 = cosmo.H0.to(u.Gyr ** (-1))  # in 1/Gyr


def _shape(cfg):
    """(z_rei, gamma, omega, beta, c_omega) for a configuration. beta and c_omega are derived, never stored."""
    r = resolve(cfg).reion
    return r.z_reion, r.gamma, r.omega, _physics.beta(r.z_reion, r.gamma), _physics.c_omega(r.omega)


def __getattr__(name):
    """Read-only access to the shape parameters of the default configuration (reion.omega, reion.beta, ...)."""
    names = ("z_rei", "gamma", "omega", "beta", "c_omega")
    if name in names:
        return _shape(None)[names.index(name)]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def s(x, y):
    """
    Suppression shape function (Okamoto, Gnedin & Frenk 2008): s(mu_c, omega) is the equilibrium
    fraction of the universal baryon accretion rate a halo retains.

    Delegates to reionization_physics.s -- see that module for the definition.
    """
    return _physics.s(x, y)


def M_c(z, cfg=None):
    """
    Characteristic mass scale for reionization (Okamoto et al. 2008): the mass at which the baryon
    fraction is suppressed by a factor of two compared to the universal value, because of the UV
    background. z_rei and gamma come from the configuration.
    """
    z_rei, gamma, *_ = _shape(cfg)
    return _physics.M_c(z, z_rei, gamma)


def mu_c(z, m_halo, cfg=None):
    """m_halo / M_c(z). Delegates to reionization_physics.mu_c."""
    z_rei, gamma, *_ = _shape(cfg)
    return _physics.mu_c(z, m_halo, z_rei, gamma)


def X(z, m_halo, cfg=None):
    _, _, omega, _, c_omega = _shape(cfg)
    M_omega = (M_c(z, cfg) / m_halo) ** (omega)
    X_val = (3 * c_omega * M_omega) / (1 + c_omega * M_omega)
    return X_val


def epsilon(z, cfg=None):
    # u = (z/beta)**gamma grows huge well within uv_suppression's z<=10 mask (about 640 at z=10 for the default
    # z_reion=7, gamma=15, and above 8000 for z_reion=6). exp(u)/(1+exp(u))**2 is written as exp(-u)/(1+exp(-u))**2, the same
    # function, which is 0 (not inf/inf = nan) when exp(u) overflows. Before this, z_reion <= 6.5 gave nan for every halo.
    _, gamma, _, beta, _ = _shape(cfg)
    with np.errstate(over="ignore"):
        u = (z / beta) ** gamma
        part2 = (gamma * (z ** (gamma - 1)) / (beta**gamma)) * np.exp(-u) / (1 + np.exp(-u)) ** 2
        epsilon = (0.63) / (1 + np.exp(u)) + part2
    return epsilon


def uv_suppression_equilibrium(z_val, m_halo, cfg=None):
    """
    Fraction of the cosmological baryon accretion a halo of mass m_halo keeps at redshift z_val: the equilibrium suppression
    s(M / M_c(z), omega) of Okamoto, Gnedin & Frenk (2008), a function of mass and redshift only (1 above z = 10, 0 for a
    halo with no mass). No dependence on the growth rate, so it does not depend on the time step of the trees.
    """
    z_val = np.asarray(z_val, dtype=float)
    m_halo = np.asarray(m_halo, dtype=float)
    suppression = np.ones_like(z_val)
    mask = z_val <= 10
    mh = m_halo[mask]
    resolved = mh > 0
    values = s(mu_c(z_val[mask], np.where(resolved, mh, 1.0), cfg), _shape(cfg)[2])
    suppression[mask] = np.where(resolved, values, 0.0)
    return suppression


def uv_suppression(z_val, m_halo, mdot_halo, cfg=None):
    z_val = np.asarray(z_val)
    m_halo = np.asarray(m_halo)
    mdot_halo = np.asarray(mdot_halo)

    suppression = np.ones_like(z_val)

    mask = z_val <= 10

    z_masked = z_val[mask]
    mh_masked = m_halo[mask]
    mdot_masked = mdot_halo[mask]

    # A halo/step with no mass and/or no accretion rate (e.g. a
    # not-yet-formed progenitor -- see build_trees_from_pymctrees.py's
    # pre-formation zeroing) has no accretion to suppress: the correct
    # answer is 0 regardless of what this function's internals would
    # otherwise evaluate to (gas_inflow_rate multiplies this output by
    # halo_mass_dot, which is already 0 there). But m_halo=0 feeds into
    # mu_c/X/s as a divide-by-zero (M_c(z)/m_halo, and s's x**(-y) with
    # x=mu_c=0), producing inf/nan that then poisons the already-correct
    # "0 *" via 0*nan=nan -- so this is guarded here directly rather than
    # relying on the caller's multiplication to save it, the same
    # principle as the mdot_halo=0 guard below.
    resolved = (mh_masked > 0) & (mdot_masked > 0)
    mh_safe = np.where(resolved, mh_masked, 1.0)
    mdot_safe = np.where(resolved, mdot_masked, 1.0)

    mu_vals = mu_c(z_masked, mh_safe, cfg)
    x_vals = X(z_masked, mh_safe, cfg)
    eps_vals = epsilon(z_masked, cfg)
    H_vals = cosmo.H(z_masked).value

    ratio_term = 2 * eps_vals * mh_safe * x_vals * (1 + z_masked) * H_vals / mdot_safe
    suppressed = s(mu_vals, _shape(cfg)[2]) * ((1 + x_vals) - ratio_term)
    suppressed = np.where(resolved, suppressed, 0.0)

    suppression[mask] = np.maximum(suppressed, 0.0)

    return suppression

from .config import retire_module_values as _retire  # noqa: E402

_retire(__name__, {
    "z_rei": "config.reion.z_reion", "gamma": "config.reion.gamma", "omega": "config.reion.omega",
    "beta": "config.reion.z_reion or gamma (beta is derived)", "c_omega": "config.reion.omega (c_omega is derived)",
})
