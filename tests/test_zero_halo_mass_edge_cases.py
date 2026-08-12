"""
Regression tests for real bugs found while validating pymctrees-generated
merger trees: a halo/step with zero mass and/or zero accretion rate (e.g. a
not-yet-formed progenitor -- see scripts/build_trees_from_pymctrees.py's
pre-formation zeroing) must produce zero suppressed accretion / zero wind
mass-loading, not NaN.

Both ashvini.reionization.uv_suppression and
ashvini.supernovae_feedback.mass_loading_factor compute a ratio whose
numerator is finite and denominator is the (possibly zero) halo mass or
accretion rate; the caller in both cases already multiplies the result by a
quantity that is itself correctly zero when the halo doesn't exist (e.g.
star_formation_rate_for_winds, halo_mass_dot), so 0 * inf/nan = nan was
silently corrupting otherwise-correct zero results before these functions
guarded the division directly.
"""

import numpy as np

from ashvini import reionization as reion
from ashvini import supernovae_feedback as sn


def test_uv_suppression_zero_halo_mass_and_rate_gives_zero_not_nan():
    z = np.array([5.0, 7.0, 9.0])
    m_halo = np.array([0.0, 0.0, 0.0])
    mdot_halo = np.array([0.0, 0.0, 0.0])

    result = reion.uv_suppression(z, m_halo, mdot_halo)

    assert not np.any(np.isnan(result))
    assert not np.any(np.isinf(result))
    assert np.all(result == 0.0)


def test_uv_suppression_mixed_zero_and_resolved_halos():
    # A batch with both not-yet-formed (zero) and resolved haloes at the
    # same step must not let the zero entries poison the resolved ones
    # (or vice versa) via shared array-wide operations. 1e8/1e9 Msun are
    # both fully suppressed at z=6 under the default reionization params
    # (a real, legitimate model result, not a bug -- checked directly),
    # so a much higher mass is used here to get a genuinely nonzero,
    # partially-suppressed resolved case to contrast against the zero one.
    z = np.array([6.0, 6.0, 6.0])
    m_halo = np.array([0.0, 1e8, 1e14])
    mdot_halo = np.array([0.0, 1e7, 1e12])

    result = reion.uv_suppression(z, m_halo, mdot_halo)

    assert not np.any(np.isnan(result))
    assert result[0] == 0.0
    assert result[1] == 0.0
    assert result[2] > 0.0


def test_uv_suppression_zero_mass_but_nonzero_rate_gives_zero():
    # An edge case distinct from the "both zero" case: zero mass alone
    # (m_halo divides M_c(z) inside X()/mu_c()) must also be guarded.
    z = np.array([6.0])
    m_halo = np.array([0.0])
    mdot_halo = np.array([1e7])

    result = reion.uv_suppression(z, m_halo, mdot_halo)
    assert not np.any(np.isnan(result))
    assert result[0] == 0.0


def test_mass_loading_factor_zero_halo_mass_gives_zero_not_nan():
    result = sn.mass_loading_factor(
        redshift=np.array([6.0, 6.0]),
        halo_mass=np.array([0.0, 1e8]),
        stellar_metallicity=np.array([0.0, 0.01]),
    )
    assert not np.any(np.isnan(result))
    assert not np.any(np.isinf(result))
    assert result[0] == 0.0
    assert result[1] > 0.0
