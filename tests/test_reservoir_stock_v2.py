"""
Tests of ashvini.reservoir_stock_v2 against the frozen ashvini.reservoir_stock.

The central test is the regression: with the three switches at their legacy values the revised integrator must reproduce
the frozen one exactly, so that every difference seen with a switch on is attributable to that switch alone.
"""
import numpy as np
import pytest

from ashvini import reservoir_stock as frozen
from ashvini import reservoir_stock_v2 as v2
from ashvini.utils import time_at_z, z_at_time


def planck_grid(n, z0=25.0, z1=5.0):
    t = np.linspace(time_at_z(np.array([z0]))[0], time_at_z(np.array([z1]))[0], n)
    return z_at_time(t), t


def smooth_history(n, M0=3e10, alpha=0.7):
    z, t = planck_grid(n)
    M = np.where(M0 * np.exp(-alpha * (z - 5.0)) > 1e4, M0 * np.exp(-alpha * (z - 5.0)), 0.0)
    return z, t, M[None, :]


def jumpy_history(n, N=6, M0=3e10, seed=3):
    """Monotone, merger-like histories: smooth growth plus a few random jumps, every halo resolved from the first node."""
    rng = np.random.default_rng(seed)
    z, t = planck_grid(n)
    M = np.empty((N, n))
    for i in range(N):
        base = M0 * np.exp(-0.7 * (z - 5.0))
        jumps = np.zeros(n)
        for idx in rng.choice(np.arange(5, n - 5), 6, replace=False):
            jumps[idx:] += rng.uniform(0.05, 0.3) * base[idx]
        M[i] = np.maximum.accumulate(base + jumps)
        M[i] *= M0 / M[i, -1]
        M[i] = np.maximum(M[i], 2e4)
    return z, t, M


def rate_from_grid(M, t):
    return np.gradient(M, t, axis=1)


@pytest.mark.parametrize("builder,kw", [(smooth_history, dict()), (jumpy_history, dict())])
def test_legacy_switches_reproduce_frozen_exactly(builder, kw):
    z, t, M = builder(401, **kw)
    rate = rate_from_grid(M, t)
    for seed in (1e3, 1e6):
        for params in (dict(), dict(sigma_lnj=1.5, R_nuc_pc=300.0)):
            a = frozen.run_reservoir_stock(M, z, seed, halo_growth_rate=rate, **params)
            b = v2.run_reservoir_v2(M, z, seed, halo_growth_rate=rate, **v2.LEGACY, **params)
            for key in ("bh_mass", "stars_nuc", "stars_gal", "gas_nuc", "gas_gal", "gas_unavailable", "mdot_in", "t_ff_nuc",
                        "wind_taken_agn", "wind_taken_sn", "transfer_nuc"):
                assert np.array_equal(a[key], b[key]), key


def test_paper_cosmology_matches_astropy_and_inverts():
    from astropy.cosmology import FlatLambdaCDM
    cos = v2.PaperCosmology()
    ref = FlatLambdaCDM(H0=67.8, Om0=0.308, Tcmb0=0)          # matter + Lambda only
    z = np.array([0.0, 2.0, 5.0, 10.0, 25.0])
    assert np.allclose(cos.time_at_z(z), ref.age(z).value, rtol=1e-7)
    assert np.allclose(cos.z_at_time(cos.time_at_z(z)), z, rtol=1e-9, atol=1e-9)
    assert np.allclose(cos.H_per_s(z), (ref.H(z).value * 1e5 / (1e6 * 3.0856775814913673e18)), rtol=1e-7)


def test_uniform_grid_endpoints_and_spacing():
    cos = v2.PaperCosmology()
    z = v2.uniform_time_grid(25.0, 5.0, 801, cos)
    t = cos.time_at_z(z)
    assert z[0] == 25.0 and z[-1] == 5.0
    assert np.allclose(np.diff(t), np.diff(t)[0], rtol=1e-9)


def test_conservative_inflow_conserves_the_halo_budget_on_merger_like_histories():
    cos = v2.PaperCosmology()
    zg = v2.uniform_time_grid(25.0, 12.0, 401, cos)            # above z = 10 the UV factor is unity
    tg = cos.time_at_z(zg)
    # tree-like nodes (coarser than the grid, with mass jumps), interpolated onto the grid
    z_nodes = np.linspace(25.0, 12.0, 60)
    rng = np.random.default_rng(1)
    base = 1e9 * np.exp(-0.7 * (z_nodes - 12.0))
    steps = np.cumsum(rng.uniform(0.0, 1.0, z_nodes.size) ** 6) * 2e8
    Mn = (base + steps)[None, :] + 2e5
    Mg = v2.interpolate_tree_onto_grid_v2(Mn, z_nodes, zg, cos)
    out = v2.run_reservoir_v2(Mg, zg, 1e3, inflow="conservative", cosmology=cos)
    gain = Mg[0, -1] - Mg[0, 0]
    # zeta_ch is 1 - s(M / 4e11, 4), within 1e-5 of unity for M < 1e9; the first step uses the start-of-step mass
    assert out["cum_in"][0] == pytest.approx(0.156 * gain, rel=1e-4)
    assert out["budget_ratio"][0] == pytest.approx(1.0, rel=1e-4)
    assert abs(out["baryon_residual"][0]) < 1e-12


def test_baryon_conservation_all_switch_combinations():
    z, t, M = jumpy_history(301)
    rate = rate_from_grid(M, t)
    for dm in (False, True):
        for infl in ("rate", "conservative"):
            o = v2.run_reservoir_v2(M, z, 1e4, halo_growth_rate=rate, dm_nuclear=dm, inflow=infl, cosmology=None)
            assert np.max(np.abs(o["baryon_residual"])) < 1e-12


def test_nuclear_dark_matter_shortens_the_nuclear_free_fall_time_and_raises_delivery():
    z, t, M = smooth_history(401)
    rate = rate_from_grid(M, t)
    a = v2.run_reservoir_v2(M, z, 1e3, halo_growth_rate=rate, dm_nuclear=False, inflow="rate", cosmology=None)
    b = v2.run_reservoir_v2(M, z, 1e3, halo_growth_rate=rate, dm_nuclear=True, inflow="rate", cosmology=None)
    live = np.isfinite(a["t_ff_nuc"]) & (a["gas_nuc"] > 0)
    assert np.all(b["t_ff_nuc"][live] < a["t_ff_nuc"][live])
    assert b["transfer_nuc"].sum() > a["transfer_nuc"].sum()


def test_single_crossing_report_finds_one_crossing_on_a_smooth_history():
    cos = v2.PaperCosmology()
    zg = v2.uniform_time_grid(25.0, 5.0, 201, cos)
    t = cos.time_at_z(zg)
    M = (3e10 * np.exp(-0.7 * (zg - 5.0)))[None, :]
    M = np.where(M > 1e4, M, 0.0)
    rep = v2.single_crossing_report(M, zg, n_grid=17)
    assert rep["n_sign_changes"][0] == 1


def test_substepping_is_identity_for_one_substep_and_conserves_baryons():
    z, t, M = jumpy_history(301)
    a = v2.run_reservoir_v2(M, z, 1e4, n_sub=1)
    b = v2.run_reservoir_v2(M, z, 1e4, n_sub=4)
    assert np.max(np.abs(b["baryon_residual"])) < 1e-12
    # same physics, so the answers agree to the first-order accuracy of the explicit step
    assert np.allclose(a["stars_mass"][:, -1], b["stars_mass"][:, -1], rtol=0.15)
    c = v2.run_reservoir_v2(M, z, 1e4, n_sub=4)
    assert np.array_equal(b["bh_mass"], c["bh_mass"])


def test_substepping_converges_towards_a_refined_grid():
    """Sub-stepping the consumption part of a coarse grid must approach what a finer grid gives; the explicit step does not."""
    cos = v2.PaperCosmology()
    Mn_z = np.linspace(25.0, 5.0, 400)
    Mn = 3e10 * np.exp(-0.7 * (Mn_z - 5.0))
    Mn = np.where(Mn > 1e4, Mn, 0.0)[None, :]
    def final(n, n_sub):
        zg = v2.uniform_time_grid(25.0, 5.0, n, cos)
        Mg = v2.interpolate_tree_onto_grid_v2(Mn, Mn_z, zg, cos)
        o = v2.run_reservoir_v2(Mg, zg, 1e3, n_sub=n_sub)
        return o["stars_nuc"][0, -1] / o["stars_mass"][0, -1], o["bh_mass"][0, -1]
    ref = final(801, 16)
    err1 = abs(final(201, 1)[0] - ref[0]) / ref[0]
    err4 = abs(final(201, 4)[0] - ref[0]) / ref[0]
    assert err4 < err1


# ---------------------------------------------------------------------------------------------------------------------
# PCH08 main-branch sampler (needs a foraois checkout; skipped otherwise)
# ---------------------------------------------------------------------------------------------------------------------
import os
import pytest


@pytest.mark.skipif(not os.environ.get("FORAOIS_ROOT"), reason="FORAOIS_ROOT not set")
def test_pch08_main_branch_reproducible_and_compliant():
    """Same seed gives identical histories, the step floor is never hit (so N_upper <= target everywhere), and the grid is
    chronological, starting at zero (unformed) or a resolved mass and ending at M0 at z = 5."""
    import sys
    from pathlib import Path
    root = Path(os.environ["FORAOIS_ROOT"])
    sys.path.insert(0, str(root / "src"))
    from ashvini.pch08_main_branch import build_main_branches
    cfg = root / "config" / "menon_power_2024.yml"
    kw = dict(M0_msun=3e10, n_trees=2, z0=5.0, z_max=10.0, m_res_msun=1e4, dz_grid=0.05, seed0=7)
    hm1, zz, st = build_main_branches(cfg, **kw)
    hm2, _, _ = build_main_branches(cfg, **kw)
    assert np.array_equal(hm1, hm2)
    assert st[:, 2].sum() == 0
    assert zz[0] > zz[-1] and abs(zz[-1] - 5.0) < 1e-12
    assert np.allclose(hm1[:, -1], 3e10)
    assert np.all(np.diff(hm1, axis=1) >= -1e-6 * hm1[:, -1:])  # main branch mass never decreases forward in time
