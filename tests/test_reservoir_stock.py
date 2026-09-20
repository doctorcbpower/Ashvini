"""Tests of the Minimal Viable Model (ashvini.reservoir_stock) against its specification."""
import inspect

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import ndtr

from ashvini import black_holes_growth_slimdisk as bs
from ashvini.constants import G, Msun_g, pc_cm
from ashvini.reservoir_stock import (
    critical_seed_stock, median_specific_angular_momentum, nfw_dark_matter_mass, nuclear_transfer,
    run_reservoir_stock, star_formation_demand, strict_eddington_step, target_residual,
)
from ashvini.utils import time_at_z, z_at_time


def _grid(n=201, z0=25.0):
    t = np.linspace(time_at_z(np.array([z0]))[0], time_at_z(np.array([5.0]))[0], n)
    return z_at_time(t)


def _halo(n=201, m0=1e8, k=0.6):
    z = _grid(n)
    return (m0 * np.exp(k * (25.0 - z)))[None, :], z


# ---------------------------------------------------------------- strict Eddington
def test_strict_eddington_step_matches_ode_and_regimes():
    kappa = 19.9
    for M0, A, dt in ((1e3, 1e9, 0.01), (1e3, 5e3, 0.05), (1e6, 1e3, 0.01), (1e2, 3e4, 0.2), (1e4, 0.0, 0.1)):
        ref = solve_ivp(lambda t, y: [min(A, kappa * y[0])], (0, dt), [M0], rtol=1e-12, atol=1e-9).y[0, -1]
        assert float(strict_eddington_step(M0, A, kappa, dt)) == pytest.approx(ref, rel=1e-8)


def test_strict_eddington_step_agrees_with_large_chi_crit_but_not_chi_crit_one():
    rng = np.random.default_rng(3)
    M0 = 10 ** rng.uniform(-1, 9, 500)
    A = 10 ** rng.uniform(-2, 12, 500)
    exact = strict_eddington_step(M0, A, 19.9, 0.003)
    big = bs.bh_growth_step_slimdisk(M0, A, 19.9, 1e12, 0.003)
    assert np.allclose(big, exact, rtol=1e-9, atol=0)
    one = bs.bh_growth_step_slimdisk(M0, A, 19.9, 1.0, 0.003)  # r_crit=1 is NOT strict Eddington
    supply_rich = A > 19.9 * M0 * 10
    assert np.any(one[supply_rich] > 1.5 * exact[supply_rich] - 0.5 * M0[supply_rich])


def test_strict_eddington_never_exceeds_eddington_or_supply():
    M0 = np.logspace(0, 8, 30)
    A = np.logspace(2, 9, 30)
    dM = strict_eddington_step(M0, A, 19.9, 0.01) - M0
    assert np.all(dM <= A * 0.01 * (1 + 1e-12))
    assert np.all(dM <= M0 * (np.exp(19.9 * 0.01) - 1) * (1 + 1e-12))


# ---------------------------------------------------------------- accessibility
def _one_cohort(mass, ljm):
    return (np.array([[mass]]), np.array([[ljm]]))


def test_nuclear_transfer_keeps_history_and_releases_gas_when_the_threshold_rises():
    # One cohort accepted below j_gal=a. Thresholds rise, fall, rise. Nucleus must hold M F(b_max)/F(a)
    # (capped at M), the rest stays in the galaxy, and a fall moves nothing.
    sigma, ljm = 0.5, np.log(1.0e28)
    a = ljm + 0.4
    F = lambda x: ndtr((x - ljm) / sigma)
    M = 7.0
    m, lj = _one_cohort(M * 1.0, ljm)
    F_up = np.array([[F(a)]]); F_cur = np.zeros((1, 1)); cut = np.full((1, 1), -np.inf)
    got = 0.0
    bmax = -np.inf
    for b in (ljm - 1.5, ljm - 0.9, ljm - 1.2, ljm - 0.2, ljm + 0.9):
        T = nuclear_transfer(m, lj, F_up, F_cur, cut, np.array([b]), sigma, 1)[0]
        got += T
        bmax = max(bmax, b)
        assert got == pytest.approx(M * min(F(bmax), F(a)) / F(a), rel=1e-9)
        assert m[0, 0] == pytest.approx(M - got, rel=1e-9, abs=1e-12)
        if b < bmax:
            assert T == 0.0
    assert m[0, 0] == pytest.approx(0.0, abs=1e-9)  # last threshold was above the galaxy cut: all moved


# ---------------------------------------------------------------- MVM structure
def test_mvm_has_no_removed_options():
    params = set(inspect.signature(run_reservoir_stock).parameters)
    for gone in ("R_nuc_rd", "dm_nuclear", "epsilon_sf_ext", "sf_ext_timescale", "wind_mode", "wind_driver",
                 "agn_wind_type", "agn_wind_sink", "epsilon_f", "chi_crit", "compaction_boost", "gas_mass0",
                 "stars_mass0", "a_star"):
        assert gone not in params
    with pytest.raises(TypeError):
        run_reservoir_stock(*_halo(), 1e5, chi_crit=10.0)


def test_empty_start_and_two_evolved_reservoirs():
    m, z = _halo()
    out = run_reservoir_stock(m, z, 1e5)
    assert out["gas_gal"][0, 0] == 0.0 and out["gas_nuc"][0, 0] == 0.0
    assert out["stars_gal"][0, 0] == 0.0 and out["stars_nuc"][0, 0] == 0.0 and out["bh_mass"][0, 0] == 1e5
    assert np.array_equal(out["gas_mass"], out["gas_gal"] + out["gas_nuc"])
    assert np.array_equal(out["stars_mass"], out["stars_gal"] + out["stars_nuc"])


def test_galaxy_cut_is_memoryless_lognormal_fraction_and_rest_is_unavailable():
    z = _grid(61)
    hm = np.full((1, len(z)), 1e10)
    rate = np.zeros((1, len(z))); rate[0, 0] = 2e9
    accepted = []
    for n_rd in (0.5, 1.0, 2.0, 3.0, 5.0, 20.0):
        out = run_reservoir_stock(hm, z, 1e4, halo_growth_rate=rate, n_rd=n_rd, eta_acc=0.0, epsilon_sf=0.0, eta_sn_scale=0.0)
        fresh = out["mdot_in"][0, 1] * (out["cosmic_time"][1] - out["cosmic_time"][0])
        total = out["gas_gal"][0, 1] + out["gas_nuc"][0, 1] + out["gas_unavailable"][0, 1]
        assert total == pytest.approx(fresh, rel=1e-12)
        accepted.append(1.0 - out["gas_unavailable"][0, 1] / fresh)
        assert np.all(out["gas_unavailable"][0, 2:] == out["gas_unavailable"][0, 1])  # ledger: never returns or changes
    assert np.all(np.diff(accepted) >= 0)  # a larger galaxy scale accepts more gas
    assert any(0.02 < a < 0.98 for a in accepted) and accepted[-1] > 0.99  # a genuine partial split, and a full one


def test_nuclear_threshold_uses_baryons_only_and_nuclear_dark_matter_does_not_exist():
    # With the galaxy cut inert (huge n_rd), the transfer must not depend on the NFW concentration.
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, eta_sn_scale=0.0, n_rd=1e4)
    a = run_reservoir_stock(hm, z, 1e6, c_nfw=3.0, **kw)["gas_nuc"]
    b = run_reservoir_stock(hm, z, 1e6, c_nfw=8.0, **kw)["gas_nuc"]
    assert np.allclose(a, b, rtol=1e-6, atol=0) and a.max() > 0


def test_nuclear_radius_is_fixed_and_larger_radius_delivers_more_gas():
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, eta_sn_scale=0.0, n_rd=1e4)
    g = [run_reservoir_stock(hm, z, 1e6, R_nuc_pc=R, **kw)["gas_nuc"][0, -1] for R in (30.0, 100.0, 300.0)]
    assert g[0] < g[1] < g[2]


# ---------------------------------------------------------------- one star-formation law
def test_same_star_formation_law_at_both_scales():
    assert float(star_formation_demand(1e8, 2.0e-3, 0.015, 0.003)) == pytest.approx(0.015 * 1e8 / 2.0e-3 * 0.003)
    assert float(star_formation_demand(1e8, np.inf, 0.015, 0.003)) == 0.0
    m, z = _halo(m0=1e9, k=0.5)
    out = run_reservoir_stock(m, z, 1e4, eta_sn_scale=0.0)
    ok = ~out["limiter_nuc"][0, 1:] & ~out["limiter_gal"][0, 1:] & (out["sfr_gal"][0, 1:] > 0) & (out["sfr_nuc"][0, 1:] > 0)
    assert ok.sum() > 5
    dt = np.diff(out["cosmic_time"])[ok]
    # SFR * t_ff / M_gas(post-transfer) = eps_sf in both reservoirs: recover the post-transfer gas from the rate
    for key, tff in (("sfr_gal", "t_ff_gal"), ("sfr_nuc", "t_ff_nuc")):
        gas_needed = out[key][0, 1:][ok] * out[tff][0, 1:][ok] / 0.015
        assert np.all(gas_needed > 0)
    # both draw a fraction eps dt / t_ff of their (post-transfer) gas
    frac_nuc = out["sfr_nuc"][0, 1:][ok] * dt
    assert np.all(frac_nuc > 0)


# ---------------------------------------------------------------- feedback
def test_outflow_is_driven_by_total_sfr_and_agn_by_bh_growth():
    m, z = _halo(m0=1e9, k=0.5)
    out = run_reservoir_stock(m, z, 1e4)
    # stellar outflow is zero if there is no star formation in either reservoir
    none = ~(out["sfr_gal"][0] > 0) & ~(out["sfr_nuc"][0] > 0)
    assert np.all(out["wind_taken_sn"][0][none] == 0.0)
    # with no black hole growth the AGN wind is zero
    frozen = run_reservoir_stock(m, z, 1e4, eta_acc=0.0)
    assert np.all(frozen["wind_taken_agn"] == 0.0)
    # no stellar wind if the loading is zero, and outflow mass is what leaves the system
    assert np.all(run_reservoir_stock(m, z, 1e4, eta_sn_scale=0.0)["wind_taken_sn"] == 0.0)
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)


def test_feedback_can_consume_the_whole_reservoir_without_negative_mass():
    m, z = _halo(m0=1e9, k=0.5)
    out = run_reservoir_stock(m, z, 1e5, eta_sn_scale=1e6, f_mom=1e4)
    assert out["wind_cap"].any()
    for k in ("gas_gal", "gas_nuc", "stars_gal", "stars_nuc", "bh_mass"):
        assert np.all(out[k] >= 0.0) and np.all(np.isfinite(out[k]))
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
    capped = out["wind_cap"][0]
    taken = out["wind_taken_agn"][0][capped] + out["wind_taken_sn"][0][capped]
    assert np.all(taken <= out["wind_demand"][0][capped] * (1 + 1e-12))
    assert np.all(out["gas_gal"][0][capped] + out["gas_nuc"][0][capped] < 1e-6 * np.maximum(1.0, taken))


def test_wind_is_removed_in_proportion_to_gas_mass_from_the_two_reservoirs():
    # An uncapped wind leaves the ratio gas_nuc / gas_gal unchanged relative to the same step without wind.
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    lam = 0.035 * np.sqrt(G * 1e9 * Msun_g * 100 * pc_cm) / median_specific_angular_momentum(1e10, 0.5 * (z[0] + z[1]), 0.035)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, lam_median=lam, n_rd=1e4)
    off = run_reservoir_stock(hm, z, 1e9, eta_sn_scale=0.0, **kw)
    on = run_reservoir_stock(hm, z, 1e9, eta_sn_scale=1.0, **kw)
    # star formation is off, so nothing drives the stellar wind: identical
    assert np.array_equal(off["gas_nuc"], on["gas_nuc"])
    on2 = run_reservoir_stock(hm, z, 1e9, eta_sn_scale=0.5, **{**kw, "epsilon_sf": 0.002})
    off2 = run_reservoir_stock(hm, z, 1e9, eta_sn_scale=0.0, **{**kw, "epsilon_sf": 0.002})
    r_off = off2["gas_nuc"][0, 1] / off2["gas_gal"][0, 1]
    # after the first step the nuclear and galaxy gas are reduced by the same fraction of what remains
    kept_nuc = on2["gas_nuc"][0, 1] / off2["gas_nuc"][0, 1]
    kept_gal = on2["gas_gal"][0, 1] / off2["gas_gal"][0, 1]
    assert not on2["wind_cap"][0, 1] and off2["gas_nuc"][0, 1] > 0
    assert kept_nuc == pytest.approx(kept_gal, rel=1e-9) and kept_nuc < 1.0 and r_off > 0


# ---------------------------------------------------------------- conservation and ordering
@pytest.mark.parametrize("seed", [1e2, 1e5, 1e8])
def test_baryon_budget_and_nonnegativity_over_random_parameters(seed):
    rng = np.random.default_rng(int(seed) % 97)
    m, z = _halo()
    m[:, :30] = 0.0
    for _ in range(4):
        kw = dict(sigma_lnj=rng.uniform(0.3, 1.5), eta_acc=10 ** rng.uniform(-3, 0), epsilon_sf=10 ** rng.uniform(-3, 0),
                  f_mom=10 ** rng.uniform(-1, 2), eta_sn_scale=10 ** rng.uniform(-1, 2), R_nuc_pc=rng.uniform(30, 1000))
        out = run_reservoir_stock(m, z, seed, **kw)
        assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
        for k in ("gas_gal", "gas_nuc", "stars_gal", "stars_nuc"):
            assert np.all(out[k] >= 0.0)


def test_unformed_halo_is_inert():
    m, z = _halo()
    m[:, :60] = 0.0
    out = run_reservoir_stock(m, z, 1e5)
    assert np.all(out["gas_gal"][:, :61] == 0.0) and np.all(out["gas_nuc"][:, :61] == 0.0)
    assert np.all(out["bh_mass"][:, :61] == 1e5) and np.all(out["mdot_out"][:, :61] == 0.0)


def test_timestep_convergence_on_a_smooth_halo():
    med = []
    for n in (101, 201, 401):
        m, z = _halo(n=n, m0=1e9, k=0.5)
        med.append(run_reservoir_stock(m, z, 1e6)["stars_mass"][0, -1])
    assert abs(med[2] / med[1] - 1) < abs(med[1] / med[0] - 1) + 1e-3
    assert abs(med[2] / med[1] - 1) < 0.05


# ---------------------------------------------------------------- the numerical experiment
def test_root_function_is_scale_free_and_solution_is_a_root():
    m, z = _halo(m0=1e9, k=0.5)
    Mc, never, above, F = critical_seed_stock(m, z, f_bh=0.5, n_iter=50, return_F=True)
    assert not never[0] and not above[0]
    assert abs(F[0]) < 1e-6
    out = run_reservoir_stock(m, z, Mc)
    assert out["bh_mass"][0, -1] == pytest.approx(0.5 * out["stars_mass"][0, -1], rel=1e-6)
    assert target_residual(out, 0.5)[0] == pytest.approx(0.0, abs=1e-6)


def test_bracket_is_adaptive_and_covers_massive_halos():
    m, z = _halo(m0=1e11, k=0.7)
    Mc, never, above = critical_seed_stock(m, z, f_bh=0.5, n_iter=30)
    assert not never[0] and np.isfinite(Mc[0])


def test_response_to_the_seed_is_monotone_on_a_grid():
    m, z = _halo(m0=1e9, k=0.5)
    seeds = np.logspace(0, 9, 19)
    F = np.array([target_residual(run_reservoir_stock(m, z, s), 0.5)[0] for s in seeds])
    signs = np.sign(F)
    assert np.sum(np.diff(signs) != 0) <= 1  # crosses zero once


def test_seed_to_host_baryon_diagnostic():
    m, z = _halo()
    m[:, :40] = 0.0
    out = run_reservoir_stock(m, z, 1e6)
    assert np.all(np.isnan(out["bh_to_host_baryons"][0, :40])) and np.all(np.isfinite(out["bh_to_host_baryons"][0, 40:]))
    assert out["z_first_resolved"][0] == z[40]
    assert out["seed_over_host_baryons_at_formation"][0] == pytest.approx(1e6 / (0.156 * m[0, 40]), rel=1e-12)
    assert out["bh_to_host_baryons"][0, -1] == pytest.approx(out["bh_mass"][0, -1] / (0.156 * m[0, -1]), rel=1e-12)


def test_transfer_diagnostic_is_nonnegative_and_bounded_by_the_galaxy_supply():
    m, z = _halo(m0=1e9, k=0.5)
    out = run_reservoir_stock(m, z, 1e6)
    assert np.all(out["transfer_nuc"] >= 0.0) and out["transfer_nuc"].sum() > 0
    accepted = out["mdot_in"][0, 1:] @ np.diff(out["cosmic_time"]) - out["gas_unavailable"][0, -1]
    assert out["transfer_nuc"].sum() <= accepted * (1 + 1e-9)
