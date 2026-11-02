"""
Equilibrium UV-suppression physics (Okamoto, Gnedin & Frenk 2008), extracted so it can be reused
outside Ashvini (Stochastic_Occupation is the first such consumer) without pulling in `run_params`,
`run_params.yaml`, the filesystem layout, or any of Ashvini's galaxy-formation configuration.

This module contains ONLY the equilibrium filtering-mass building blocks -- the characteristic mass
`M_c(z)`, the mass ratio `mu_c`, and the suppression shape function `s`. Ashvini's own rate-dependent
correction (`X`, `epsilon`, and `uv_suppression`'s `mdot_halo`/`H(z)` terms) is NOT here: it is a
different, non-equilibrium physical model, not just a different implementation of the same one, and
stays in `ashvini/reionization.py` where it already lives, unchanged. See that module's own
docstring for why (its own audit found the rate-dependent form behaves problematically at low mass).

Every function here takes its physical parameters explicitly (`z_rei`, `gamma`, `omega`) -- none of
them are read from a global, a config object, or any Ashvini-specific state. `m_halo` is physical
Msun throughout (not Msun/h) -- Ashvini's own callers already convert from the tree's Msun/h before
calling into this module (see `ashvini/pymctrees_adapter.py`).

Mathematical definitions are copied verbatim from the pre-extraction `ashvini/reionization.py`
(module-level `beta`/`c_omega` derivation, and the bodies of `M_c`, `mu_c`, `s`); see
`tests/test_reionization.py` for the regression values pinned before this extraction and
`tests/test_reionization_physics.py` for direct tests of this module.
"""
import numpy as np


def beta(z_rei, gamma):
    """
    Derived shape parameter fixing M_c(z_rei) to the correct value implied by the reionisation
    redshift z_rei (Okamoto, Gnedin & Frenk 2008). Depends only on (z_rei, gamma).
    """
    return z_rei * (np.log(1.82 * (10 ** 3) * np.exp(-0.63 * z_rei) - 1) ** (-1 / gamma))


def c_omega(omega):
    """
    Derived normalisation constant used by the suppression shape function s(x, omega). Depends
    only on omega. (Previously duplicated inline inside s(x, y)'s own formula and, separately, as
    a module-level global in ashvini.reionization and ashvini.configure -- s() below now calls
    this directly instead of re-deriving the same sub-expression.)
    """
    return 2 ** (omega / 3) - 1


def M_c(z, z_rei, gamma):
    """
    Characteristic mass scale for reionisation (Okamoto, Gnedin & Frenk 2008), physical Msun --
    not Msun/h. The mass at which the baryon fraction is suppressed by a factor of two compared to
    the universal value, because of the UV background.
    """
    b = beta(z_rei, gamma)
    return 1.69 * (10 ** 10) * (np.exp(-0.63 * z) / (1 + np.exp((z / b) ** gamma)))


def mu_c(z, m_halo, z_rei, gamma):
    """m_halo / M_c(z, z_rei, gamma). m_halo: physical Msun."""
    return m_halo / M_c(z, z_rei, gamma)


def s(x, omega):
    """
    Suppression shape function (Okamoto, Gnedin & Frenk 2008): s(mu_c, omega) is the equilibrium
    fraction of the universal baryon accretion rate a halo retains. s(1, omega) = 0.5 exactly, for
    any omega (the definition of M_c: baryon fraction suppressed by a factor of two at m_halo=M_c).
    """
    c = c_omega(omega)
    return (1 + c * (x ** (-omega))) ** (-3 / omega)
