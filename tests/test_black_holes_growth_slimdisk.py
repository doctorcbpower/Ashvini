"""
Tests for the "hobbs_slimdisk" alternative BH growth model
(ashvini.black_holes_growth_slimdisk), added alongside the pre-existing
"pznk11_freefall" default -- see run_params.py's BlackHoleParams.growth_model
and MODELS.md.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from ashvini import black_holes_growth_slimdisk as bs
from ashvini import main
from ashvini.spin import epsilon_from_spin, isco_radius


# ---------------------------------------------------------------------------
# spin.py
# ---------------------------------------------------------------------------

def test_epsilon_from_spin_reference_values():
    # See spin.py's module docstring / the 2026 paper's Section 5.5.
    assert epsilon_from_spin(0.0) == pytest.approx(0.0572, abs=1e-4)
    assert epsilon_from_spin(0.5) == pytest.approx(0.0821, abs=1e-4)
    assert epsilon_from_spin(0.9) == pytest.approx(0.1558, abs=1e-4)
    assert epsilon_from_spin(0.998) == pytest.approx(0.3210, abs=1e-4)


def test_epsilon_from_spin_monotonic_increasing_with_spin():
    a_vals = np.linspace(-0.9, 0.998, 20)
    eps = epsilon_from_spin(a_vals)
    assert np.all(np.diff(eps) > 0)


def test_isco_radius_schwarzschild_is_6rg():
    assert isco_radius(0.0) == pytest.approx(6.0, rel=1e-6)


# ---------------------------------------------------------------------------
# nuclear_accretion_rate / time_freefall_enclosed
# ---------------------------------------------------------------------------

def test_time_freefall_enclosed_positive_and_finite_for_resolved_mass():
    t_ff = bs.time_freefall_enclosed(1e8, R_nuc_pc=100.0)
    assert np.isfinite(t_ff)
    assert t_ff > 0


def test_time_freefall_enclosed_infinite_for_empty_nucleus():
    assert np.isinf(bs.time_freefall_enclosed(0.0))
    assert np.isinf(bs.time_freefall_enclosed(-1.0))


def test_time_freefall_decreases_with_enclosed_mass():
    t_lo = bs.time_freefall_enclosed(1e6, R_nuc_pc=100.0)
    t_hi = bs.time_freefall_enclosed(1e10, R_nuc_pc=100.0)
    assert t_hi < t_lo


def test_nuclear_accretion_rate_scales_with_compaction_boost_squared():
    kwargs = dict(gas_mass=1e8, bh_mass=1e2, stars_mass=1e5, eta_acc=0.01, R_nuc_pc=100.0)
    A_boost1 = bs.nuclear_accretion_rate(**kwargs, compaction_boost=1.0)
    A_boost3 = bs.nuclear_accretion_rate(**kwargs, compaction_boost=3.0)
    assert A_boost3 == pytest.approx(9.0 * A_boost1, rel=1e-8)


def test_nuclear_accretion_rate_light_seed_dominated_by_gas_not_bh_mass():
    # For a light seed (M_BH << M_gas), enclosed mass ~ gas_mass -- the
    # whole reason the enclosed-mass estimator replaces naive Bondi
    # (Hobbs, Power, Nayakshin & King 2012): a Bondi (M_BH^2) estimate
    # would give a wildly different (much smaller) answer here.
    A_light_seed = bs.nuclear_accretion_rate(gas_mass=1e9, bh_mass=1e2, stars_mass=1e6)
    A_no_bh = bs.nuclear_accretion_rate(gas_mass=1e9, bh_mass=0.0, stars_mass=1e6)
    assert A_light_seed == pytest.approx(A_no_bh, rel=1e-6)


# ---------------------------------------------------------------------------
# bh_growth_step_slimdisk: closed form vs. reference
# ---------------------------------------------------------------------------

def test_slimdisk_reduces_to_pznk11_two_regime_model_as_rcrit_to_infinity():
    # As r_crit -> infinity, the slim-disc throttled regime (S) never
    # applies for any y0 > 0 (its upper mass boundary M1 -> 0), so this
    # must reduce exactly to main._bh_growth_step's Eddington-limited /
    # gas-supply-limited two-regime solver.
    y0 = np.array([1e2, 1e4, 1e7])
    A_bh = np.array([1e8, 1e6, 1e5])
    kappa_edd = 20.0
    dt = 0.01

    slimdisk_huge_rcrit = bs.bh_growth_step_slimdisk(y0, A_bh, kappa_edd, 1e12, dt)
    reference = main._bh_growth_step(y0, A_bh, kappa_edd, dt)
    assert slimdisk_huge_rcrit == pytest.approx(reference, rel=1e-6)


def test_slimdisk_r_crit_equal_one_is_not_strict_eddington():
    # Regression test for a real, previously-undiscovered bug (2026-09-15,
    # see docs/2026_paper_session_code_catalogue.md and the
    # ashvini_chi_crit_eddington_bug memory note): r_crit=1.0 makes the
    # closed form's two mass boundaries, M1=A_bh/(kappa_edd*r_crit) and
    # M2=A_bh/kappa_edd=M1*r_crit, coincide exactly, collapsing the genuine
    # Eddington-limited (exponential) regime to zero width. The result is
    # growth equal to the raw, uncapped supply rate A_bh at every mass,
    # completely independent of kappa_edd -- i.e. NOT strict
    # Eddington-limited growth, despite r_crit=1 having once been used
    # (wrongly) as the "strict Eddington" default in both
    # paper_reservoir.critical_seed_paper and critical_seed.critical_seed.
    # This test exists so that bug cannot silently return: it checks that
    # r_crit=1 and r_crit=1e12 give genuinely different growth whenever the
    # supply rate exceeds the Eddington rate (the regime where the bug bites).
    y0 = np.array([1e6])
    A_bh = np.array([1e10])  # supply far exceeds the Eddington rate at this mass
    kappa_edd = 20.0
    dt = 0.1

    growth_r_crit_one = bs.bh_growth_step_slimdisk(y0, A_bh, kappa_edd, 1.0, dt)
    growth_strict_eddington = bs.bh_growth_step_slimdisk(y0, A_bh, kappa_edd, 1e12, dt)

    # r_crit=1 must reduce to the raw, uncapped supply rate (the bug's
    # actual behaviour) ...
    assert growth_r_crit_one[0] == pytest.approx(y0[0] + A_bh[0] * dt, rel=1e-8)
    # ... and must NOT match genuine Eddington-limited growth, which caps
    # the rate at kappa_edd * M and so grows far more slowly here.
    assert growth_strict_eddington[0] < 0.5 * growth_r_crit_one[0]


@pytest.mark.parametrize("y0,A_bh,kappa_edd,r_crit,dt", [
    (100.0, 1e8, 20.0, 8.0, 0.01),      # starts in regime S, crosses into E within the step
    (100.0, 1e8, 20.0, 8.0, 0.5),       # starts in S, crosses S->E->G within the step
    (1e6, 1e8, 20.0, 8.0, 0.01),        # starts in regime E
    (1e9, 1e8, 20.0, 8.0, 0.01),        # starts in regime G
    (100.0, 1e5, 20.0, 3.0, 0.05),      # small supply, small r_crit
])
def test_slimdisk_closed_form_matches_numerical_solve_ivp(y0, A_bh, kappa_edd, r_crit, dt):
    def rhs(t, y):
        M = max(y[0], 1e-300)
        ratio = A_bh / (kappa_edd * M)
        cap = kappa_edd * M if ratio <= r_crit else A_bh / r_crit
        return [min(A_bh, cap)]

    sol = solve_ivp(rhs, [0.0, dt], [y0], method="LSODA", rtol=1e-11, atol=1e-8)
    numeric = sol.y[0, -1]

    closed_form = bs.bh_growth_step_slimdisk(
        np.array([y0]), np.array([A_bh]), kappa_edd, r_crit, dt
    )[0]

    assert closed_form == pytest.approx(numeric, rel=1e-4)


def test_slimdisk_growth_rate_never_negative_and_monotonic_in_time():
    y0 = np.array([1e2])
    A_bh = np.array([1e8])
    kappa_edd = 20.0
    r_crit = 8.0
    masses = [
        bs.bh_growth_step_slimdisk(y0, A_bh, kappa_edd, r_crit, dt)[0]
        for dt in [0.0, 0.001, 0.01, 0.1, 1.0]
    ]
    assert masses[0] == pytest.approx(y0[0])
    assert all(m2 >= m1 for m1, m2 in zip(masses[:-1], masses[1:]))


def test_slimdisk_no_supply_leaves_mass_unchanged():
    y0 = np.array([1e4])
    result = bs.bh_growth_step_slimdisk(y0, np.array([0.0]), 20.0, 8.0, 1.0)
    assert result[0] == pytest.approx(y0[0])


def test_slimdisk_light_seed_grows_faster_than_pznk11_default_for_same_supply():
    # The entire point of the graded r_crit cap: a light seed in a
    # supply-rich environment should grow *faster* than the standard
    # Eddington-only cap allows (super-Eddington growth), not slower.
    y0 = np.array([1e2])
    A_bh = np.array([1e8])
    kappa_edd = 20.0
    dt = 0.05

    slimdisk_result = bs.bh_growth_step_slimdisk(y0, A_bh, kappa_edd, 8.0, dt)[0]
    pznk11_result = main._bh_growth_step(y0, A_bh, kappa_edd, dt)[0]
    assert slimdisk_result > pznk11_result


# ---------------------------------------------------------------------------
# King (2003) AGN wind
# ---------------------------------------------------------------------------

def test_king_agn_wind_rate_zero_for_unformed_halo():
    rate = bs.king_agn_wind_rate(bh_growth_rate=1e4, halo_mass=0.0, redshift=10.0)
    assert rate == pytest.approx(0.0)


def test_king_agn_wind_rate_scales_linearly_with_bh_growth_rate():
    r1 = bs.king_agn_wind_rate(bh_growth_rate=1e4, halo_mass=1e11, redshift=10.0)
    r2 = bs.king_agn_wind_rate(bh_growth_rate=2e4, halo_mass=1e11, redshift=10.0)
    assert r2 == pytest.approx(2.0 * r1, rel=1e-8)


def test_king_agn_wind_rate_stronger_in_shallower_halo():
    # Mdot_wind ~ 1/sigma^2, sigma increasing with halo mass -- feedback
    # is relatively more effective (per unit BH growth) in a shallower
    # potential well.
    r_shallow = bs.king_agn_wind_rate(bh_growth_rate=1e4, halo_mass=1e9, redshift=10.0)
    r_deep = bs.king_agn_wind_rate(bh_growth_rate=1e4, halo_mass=1e13, redshift=10.0)
    assert r_shallow > r_deep


# ---------------------------------------------------------------------------
# End-to-end wiring: growth_model switch in main.run_forest
# ---------------------------------------------------------------------------

def test_run_forest_default_growth_model_unaffected_by_slimdisk_module_import():
    # Importing black_holes_growth_slimdisk (done at main.py's module
    # level now) must not change anything about the default
    # "pznk11_freefall" path's behaviour.
    assert main._DEFAULT_GROWTH_MODEL == "pznk11_freefall"

    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z)
    assert np.all(np.isfinite(result["bh_mass"]))
    assert np.all(result["bh_mass"] >= 0)


def test_run_forest_hobbs_slimdisk_switch_produces_finite_nonnegative_output():
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z, growth_model="hobbs_slimdisk")

    assert np.all(np.isfinite(result["bh_mass"]))
    assert np.all(result["bh_mass"] >= 0)
    assert np.all(np.isfinite(result["gas_mass"]))
    assert np.all(result["gas_mass"] >= 0)
    # Some seeding channel should have fired given the halo mass range above
    assert result["bh_mass"][:, -1].max() > 0


# ---------------------------------------------------------------------------
# run_forest's per-call slimdisk parameter overrides (a_star, r_crit,
# epsilon_f, eta_acc): added so paper-analysis scripts can scan these
# without reloading run_params.yaml/re-importing the module per point.
# ---------------------------------------------------------------------------

def _run_slimdisk(halo_mass, halo_mass_rate, z, **overrides):
    return main.run_forest(halo_mass, halo_mass_rate, z, growth_model="hobbs_slimdisk", **overrides)


def test_run_forest_a_star_override_changes_bh_mass_via_epsilon():
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    bh_low_spin = _run_slimdisk(halo_mass, halo_mass_rate, z, a_star=0.0)["bh_mass"][:, -1]
    bh_high_spin = _run_slimdisk(halo_mass, halo_mass_rate, z, a_star=0.998)["bh_mass"][:, -1]
    # Higher spin -> higher radiative efficiency -> lower kappa_edd (see
    # black_holes_growth_slimdisk.eddington_rate_std_per_unit_mass) -> a
    # black hole restricted to (or near) Eddington-limited growth ends up
    # less massive at fixed halo/gas history.
    assert bh_high_spin.max() < bh_low_spin.max()

    # And the a_star=0.5 default (no override) must match the fiducial
    # module-level constant exactly -- overrides shouldn't perturb the
    # existing default path.
    bh_default = main.run_forest(halo_mass, halo_mass_rate, z, growth_model="hobbs_slimdisk")["bh_mass"]
    bh_explicit_fiducial = _run_slimdisk(halo_mass, halo_mass_rate, z, a_star=bs.a_star)["bh_mass"]
    assert np.allclose(bh_default, bh_explicit_fiducial)


def test_run_forest_r_crit_override_changes_bh_mass():
    # r_crit sets both the width of the super-Eddington-capped mass range
    # (M1 = A_bh/(kappa_edd*r_crit), shrinks as r_crit grows) and the
    # capped growth rate within it (A_bh/r_crit, also shrinks as r_crit
    # grows) -- see bh_growth_step_slimdisk's docstring. These two effects
    # pull in opposite directions, so the net effect on final M_BH is not
    # guaranteed monotonic in r_crit; this only checks the override
    # actually changes the result (the wiring), not a specific direction.
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    # NB: r_crit=1.0 is NOT strict Eddington-limited growth -- it is
    # actually uncapped, supply-limited growth (see
    # test_slimdisk_r_crit_equal_one_is_not_strict_eddington above); it is
    # used here only as a convenient second r_crit value to confirm the
    # override changes the result, not as a "strict Eddington" reference.
    bh_r_crit_one = _run_slimdisk(halo_mass, halo_mass_rate, z, r_crit=1.0)["bh_mass"][:, -1]
    bh_super_eddington = _run_slimdisk(halo_mass, halo_mass_rate, z, r_crit=100.0)["bh_mass"][:, -1]
    assert not np.allclose(bh_super_eddington, bh_r_crit_one)


def test_run_forest_eta_acc_override_changes_bh_mass():
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    bh_low_eta = _run_slimdisk(halo_mass, halo_mass_rate, z, eta_acc=bs.eta_acc / 10.0)["bh_mass"][:, -1]
    bh_high_eta = _run_slimdisk(halo_mass, halo_mass_rate, z, eta_acc=bs.eta_acc * 10.0)["bh_mass"][:, -1]
    assert bh_high_eta.max() >= bh_low_eta.max()


def test_run_forest_epsilon_direct_override_matches_a_star_when_equivalent():
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    from ashvini.spin import epsilon_from_spin
    eps = epsilon_from_spin(0.9)
    bh_via_a_star = _run_slimdisk(halo_mass, halo_mass_rate, z, a_star=0.9)["bh_mass"]
    bh_via_epsilon = _run_slimdisk(halo_mass, halo_mass_rate, z, epsilon=eps)["bh_mass"]
    assert np.allclose(bh_via_a_star, bh_via_epsilon)


def test_run_forest_rejects_both_a_star_and_epsilon():
    n = 5
    z = np.linspace(5.0, 10.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e11, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    with pytest.raises(ValueError, match="at most one"):
        main.run_forest(halo_mass, halo_mass_rate, z, a_star=0.5, epsilon=0.1)


def test_run_forest_initial_reservoir_masses_default_to_zero():
    n = 20
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z)
    assert result["gas_mass"][:, 0] == 0.0
    assert result["stars_mass"][:, 0] == 0.0


def test_run_forest_nonzero_initial_reservoir_masses_take_effect():
    n = 20
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z, gas_mass0=1e5, stars_mass0=1e3)
    assert result["gas_mass"][0, 0] == 1e5
    assert result["stars_mass"][0, 0] == 1e3
    # Should propagate forward, not be reset by the first step's update.
    assert result["stars_mass"][0, 1] >= 1e3


def test_run_forest_epsilon_f_override_changes_bh_mass():
    # A stronger AGN wind removes gas directly (dGas/dt's agn_forcing
    # term), but also feeds back nonlinearly onto BH growth (a
    # gas-starved nucleus grows more slowly, which lowers agn_forcing
    # itself the following step) -- so the net effect is not guaranteed
    # to be monotonic in epsilon_f; this only checks the override
    # actually changes the result.
    #
    # Checked against bh_mass, not gas_mass: gas_mass along this
    # particular synthetic halo trajectory settles into a self-regulating
    # quasi-equilibrium (continued inflow refills whatever the wind
    # removes) whose final-step value is essentially independent of
    # epsilon_f -- confirmed directly, not an artefact of step count --
    # even though the wind measurably changes how much mass the black
    # hole and stellar component accumulate along the way. eta_acc is
    # pinned explicitly (rather than left at the module default) so this
    # test's own discriminating power doesn't silently depend on whatever
    # SlimDiskParams.eta_acc happens to default to.
    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    bh_weak_wind = _run_slimdisk(halo_mass, halo_mass_rate, z, eta_acc=0.02, epsilon_f=1e-6)["bh_mass"][:, -1]
    bh_strong_wind = _run_slimdisk(halo_mass, halo_mass_rate, z, eta_acc=0.02, epsilon_f=1e-4)["bh_mass"][:, -1]
    assert not np.allclose(bh_strong_wind, bh_weak_wind)
