"""
Regression tests for ashvini.reionization, added before extracting its equilibrium physics
(M_c, mu_c, s, beta, c_omega) into a shared, PARAMS-independent module for reuse by
Stochastic_Occupation. Pins today's actual numerical behaviour so the extraction cannot silently
change any of it -- until now, nothing in the test suite asserted a numerical value for
uv_suppression, M_c, or s (only edge-case guards in test_zero_halo_mass_edge_cases.py and the
override-mechanism check in test_configure.py, which checks that `beta` changed, not that
uv_suppression's actual output changed as a result -- item 4 below closes that gap).

All values below were computed directly against the pre-extraction ashvini.reionization
(fiducial run_params.yaml: z_rei=7, gamma=15, omega=2) and must not change unless a deliberate,
documented change to the physics is intended.
"""
import numpy as np
import pytest

from ashvini import reionization as reion
from ashvini import configure


def teardown_function(_):
    # restore fiducial values so other test modules aren't affected by ordering (same pattern as
    # test_configure.py's own teardown).
    configure.set_params(z_rei=7)


# ---------------------------------------------------------------------
# 1. uv_suppression value regression
# ---------------------------------------------------------------------

def test_uv_suppression_value_regression():
    z = np.array([6.0, 6.5, 7.0, 8.0])
    m = np.array([1e8, 1e8, 1e9, 5e8])
    mdot = np.array([1e7, 1e12, 5e7, 1e7])
    result = reion.uv_suppression(z, m, mdot)
    expected = np.array([0.0, 0.78214391, 0.0, 1.0])
    assert result == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------
# 2. M_c value regression
# ---------------------------------------------------------------------

def test_M_c_value_regression():
    assert reion.M_c(0.0) == pytest.approx(8450000000.0, rel=1e-8)
    assert reion.M_c(6.5) == pytest.approx(75501002.8026284, rel=1e-8)


# ---------------------------------------------------------------------
# 3. s value regression
# ---------------------------------------------------------------------

def test_s_value_regression():
    # s(1, omega) = 0.5 exactly, for any omega -- this IS the definition of M_c ("the mass at
    # which the baryon fraction is suppressed by a factor of two", reionization.M_c's own
    # docstring): mu_c=1 means m_halo=M_c(z), and s(1,y) = (1+(2^(y/3)-1)*1)^(-3/y) = (2^(y/3))^(-3/y) = 0.5.
    assert reion.s(1.0, 2.0) == pytest.approx(0.5, rel=1e-12)
    assert reion.s(1.0, 15.0) == pytest.approx(0.5, rel=1e-12)
    assert reion.s(0.5, 2.0) == pytest.approx(0.16312095666666196, rel=1e-10)
    assert reion.s(2.0, 3.0) == pytest.approx(0.8888888888888888, rel=1e-10)


# ---------------------------------------------------------------------
# 4. Changing z_rei must change uv_suppression's *output*, not just the cached beta
#    (test_configure.py::test_z_rei_repatches_derived_beta_and_c_omega only checks the latter)
# ---------------------------------------------------------------------

def test_z_rei_change_propagates_to_uv_suppression_output():
    z, m, mdot = np.array([6.5]), np.array([1e8]), np.array([1e8])
    result_fiducial = reion.uv_suppression(z, m, mdot).copy()

    configure.set_params(z_rei=5)
    result_z5 = reion.uv_suppression(z, m, mdot).copy()

    assert result_z5 != pytest.approx(result_fiducial)
    # sanity: earlier reionisation (higher z_rei) suppresses more at fixed (z, m) here --
    # z_rei=5 is *later* than the fiducial 7, so less suppression is expected at z=6.5
    assert result_z5[0] >= result_fiducial[0]


# ---------------------------------------------------------------------
# 5. uv_suppression is bounded below at 0 but NOT bounded above at 1 -- document this explicitly
#    rather than assuming boundedness the way the equilibrium-only form (Stochastic_Occupation's
#    f_ret_reionization) is.
# ---------------------------------------------------------------------

def test_uv_suppression_has_no_upper_bound():
    # Large mdot suppresses the rate-dependent correction term, leaving s(mu_c,omega)*(1+X),
    # which can exceed 1 since X > 0 always -- np.maximum(suppressed, 0.0) is the only clip in
    # the implementation (reionization.uv_suppression's own code), floor only, no ceiling.
    result = reion.uv_suppression(np.array([6.5]), np.array([1e8]), np.array([1e14]))
    assert result[0] > 1.0
    assert result[0] == pytest.approx(1.1328384, rel=1e-5)


def test_uv_suppression_is_bounded_below_at_zero():
    # Small mdot makes the rate term dominate and negative; the result is clipped at 0, not
    # allowed to go negative.
    result = reion.uv_suppression(np.array([6.0]), np.array([1e6]), np.array([1e6]))
    assert result[0] >= 0.0


def test_uv_suppression_is_finite_for_a_late_reionisation_redshift():
    """z_reion = 6 gave nan for every halo: exp((z/beta)^gamma) overflowed to inf and inf/inf is nan."""
    import numpy as np
    from ashvini import reionization as reion
    from ashvini.config import DEFAULT, override
    for z_rei in (5.5, 6.0, 6.5, 7.0, 9.0):
        cfg = override(DEFAULT, **{"reion.z_reion": z_rei})
        z = np.linspace(0.0, 10.0, 41)
        out = reion.uv_suppression(z, 1e9 * np.ones_like(z), 1e8 * np.ones_like(z), cfg)
        assert np.all(np.isfinite(out)) and np.all(out >= 0.0), z_rei


def test_equilibrium_uv_form_is_the_mass_and_redshift_suppression_of_the_reservoir_model():
    import numpy as np
    from ashvini import reionization as reion
    from ashvini.paper_reservoir import uv_suppression_mass_only
    z = np.array([0.0, 3.0, 6.0, 9.0, 12.0, 6.0])
    m = np.array([1e9, 1e10, 1e11, 1e8, 1e9, 0.0])
    got = reion.uv_suppression_equilibrium(z, m)
    assert np.allclose(got, uv_suppression_mass_only(z, m), rtol=1e-12)
    assert got[4] == 1.0 and got[5] == 0.0                                  # above z = 10 and with no mass
    zz = np.full(50, 2.0); mm = np.logspace(8, 12, 50)
    assert np.all(np.diff(reion.uv_suppression_equilibrium(zz, mm)) >= 0)    # more massive haloes keep more of their accretion


def test_the_equilibrium_form_does_not_depend_on_the_halo_growth_rate_and_is_selectable():
    import numpy as np
    from ashvini.config import DEFAULT, override
    from ashvini.gas_evolve import gas_inflow_rate
    eq = override(DEFAULT, **{"reion.form": "equilibrium"})
    z, m = np.array([4.0]), np.array([3e9])
    a = gas_inflow_rate(z, m, np.array([1e8]), cfg=eq) / 1e8
    b = gas_inflow_rate(z, m, np.array([1e9]), cfg=eq) / 1e9
    assert np.allclose(a, b, rtol=1e-12)
    assert DEFAULT.reion.form == "rate_dependent"
    import pytest
    with pytest.raises(ValueError, match="reionization.form"):
        override(DEFAULT, **{"reion.form": "nope"})
