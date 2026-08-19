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
    assert main.growth_model == "pznk11_freefall"

    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z)
    assert np.all(np.isfinite(result["bh_mass"]))
    assert np.all(result["bh_mass"] >= 0)


def test_run_forest_hobbs_slimdisk_switch_produces_finite_nonnegative_output(monkeypatch):
    monkeypatch.setattr(main, "growth_model", "hobbs_slimdisk")

    n = 40
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))
    result = main.run_forest(halo_mass, halo_mass_rate, z)

    assert np.all(np.isfinite(result["bh_mass"]))
    assert np.all(result["bh_mass"] >= 0)
    assert np.all(np.isfinite(result["gas_mass"]))
    assert np.all(result["gas_mass"] >= 0)
    # Some seeding channel should have fired given the halo mass range above
    assert result["bh_mass"][:, -1].max() > 0
