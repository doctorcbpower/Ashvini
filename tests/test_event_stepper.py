"""M4: the event stepper, its samplers and its agreement with the grid stepper in the smooth limit."""
import numpy as np
import pytest

from ashvini import utils
from ashvini.config import DEFAULT, override
from ashvini.main import run_events, run_forest
from ashvini.samplers import SmoothSampler, StochasticSampler, make_sampler
from ashvini.steppers.events import EventStepper, exact_substep
from ashvini.trajectory import HaloHistory

NO_SEEDS = {"bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
            "bh.seeding.halo_mass_threshold.enabled": False, "reion.UVB_enabled": False}


def cfg_with(**changes):
    return override(DEFAULT, **{**NO_SEEDS, "sn.kernel": "imf", "sn.m_sn_max": 100.0, **changes})


def history(M_start, M_end, z0=20.0, z1=12.0, n=40):
    z = np.linspace(z0, z1, n)
    f = (z0 - z) / (z0 - z1)
    M = np.exp(np.log(M_start) + f * (np.log(M_end) - np.log(M_start)))
    return M[None, :], np.gradient(M, utils.time_at_z(z))[None, :], z


def test_exact_substep_is_exact_and_conserves():
    M_b, formed, offset = exact_substep(5.0, 2.0, 0.7, 1.3)
    A = 2.0 / 1.3 / 0.7
    assert M_b == pytest.approx(A + (5.0 - A) * np.exp(-0.7 * 1.3), rel=1e-13)
    assert M_b + formed == pytest.approx(5.0 + 2.0, rel=1e-14)
    assert 0 < offset < 1.3


def test_smooth_sampler_releases_the_expected_supernovae():
    cfg = cfg_with()
    s = SmoothSampler(cfg, K=50, dt_emit_myr=1.0)
    out = s.sample(0.0, 1.0e5)                       # first call releases at once
    assert out.weight.sum() == pytest.approx(s.nu_SN * 1.0e5, rel=1e-12)
    assert np.all(np.diff(out.t_gyr) >= 0) and np.all(out.mass_msun >= 8.0) and np.all(out.mass_msun <= 100.0)
    assert len(s.sample(0.0002, 1.0e4).t_gyr) == 0   # within 1 Myr: held back
    assert s.sample(0.002, 1.0e4).weight.sum() == pytest.approx(s.nu_SN * 2.0e4, rel=1e-12)


def test_stochastic_sampler_counts_are_poisson_with_the_expected_mean():
    cfg = cfg_with()
    s = StochasticSampler(cfg, np.random.default_rng(3))
    counts = np.array([len(s.sample(0.0, 1.0e5).t_gyr) for _ in range(400)])
    assert counts.mean() == pytest.approx(s.nu_SN * 1e5, rel=4 * np.sqrt(1 / 400 / (s.nu_SN * 1e5)))
    assert counts.var() == pytest.approx(counts.mean(), rel=0.3)


def test_baryons_are_conserved_with_and_without_the_gas_limit():
    for closure_changes in ({"sn.closure": "momentum", "sn.f_p": 0.76},
                            {"sn.closure": "energy", "sn.eps_fb": 1.0}):             # strong: the limit bites
        cfg = cfg_with(**closure_changes)
        Mh, Mdot, z = history(2e7, 1e8)
        for sampler in ("smooth", "stochastic"):
            r = run_events(Mh, Mdot, z, config=cfg, sampler=sampler, seed=1, gas_mass0=1e4)
            supplied = 1e4 + (utils.Omega_b / utils.Omega_m) * (Mh[0, -1] - Mh[0, 0])
            total = r["gas_mass"][0, -1] + r["stars_mass"][0, -1] + r["ejected_mass"][0, -1]
            assert total == pytest.approx(supplied, rel=1e-9)
            assert np.all(r["gas_mass"] >= 0) and r["n_explosions"][0] > 0


def test_without_feedback_everything_is_gas_or_stars():
    cfg = cfg_with(**{"sn.closure": "energy", "sn.eps_fb": 0.0})
    Mh, Mdot, z = history(1e8, 1e9)
    r = run_events(Mh, Mdot, z, config=cfg)
    assert r["ejected_mass"][0, -1] == 0.0
    supplied = (utils.Omega_b / utils.Omega_m) * (Mh[0, -1] - Mh[0, 0])
    assert r["gas_mass"][0, -1] + r["stars_mass"][0, -1] == pytest.approx(supplied, rel=1e-9)


def test_each_explosion_is_evaluated_at_its_own_time():
    cfg = cfg_with(**{"sn.closure": "momentum", "sn.f_p": 0.76})
    Mh, Mdot, z = history(2e7, 2e9, n=12)            # coarse grid, fast growth
    stepper = EventStepper(cfg)
    seen = []
    original = stepper.closure.per_supernova
    stepper.closure.per_supernova = lambda zz, M, Zs: (seen.append((zz, M)), original(zz, M, Zs))[1]
    state, _ = stepper.run(Mh, Mdot, z, gas_mass0=1e4)
    assert len(seen) > 100
    nodes_z, nodes_M = set(np.round(z, 9)), set(np.round(Mh[0], 3))
    off_node = [(zz, M) for zz, M in seen if round(zz, 9) not in nodes_z and round(M, 3) not in nodes_M]
    assert len(off_node) > 0.9 * len(seen)           # not the node values: the history at the explosion time
    hist = HaloHistory(Mh[0], z)
    assert all(M == pytest.approx(hist.mass_at_z(zz), rel=1e-12) for zz, M in seen)


def test_the_gas_limit_removes_exactly_the_gas_present():
    cfg = cfg_with(**{"sn.closure": "energy", "sn.eps_fb": 50.0})       # each explosion could remove far more than exists
    Mh, Mdot, z = history(2e6, 3e6, n=20)
    r = run_events(Mh, Mdot, z, config=cfg, gas_mass0=100.0)
    assert np.all(r["gas_mass"] >= 0.0)
    assert r["ejected_mass"][0, -1] > 0


def test_same_seed_same_result_and_smooth_is_deterministic():
    cfg = cfg_with(**{"sn.closure": "momentum", "sn.f_p": 0.76})
    Mh, Mdot, z = history(2e7, 1e8, n=20)
    a = run_events(Mh, Mdot, z, config=cfg, sampler="stochastic", seed=5)
    b = run_events(Mh, Mdot, z, config=cfg, sampler="stochastic", seed=5)
    c = run_events(Mh, Mdot, z, config=cfg, sampler="stochastic", seed=6)
    assert np.array_equal(a["stars_mass"], b["stars_mass"]) and not np.array_equal(a["stars_mass"], c["stars_mass"])
    assert np.array_equal(run_events(Mh, Mdot, z, config=cfg)["gas_mass"], run_events(Mh, Mdot, z, config=cfg)["gas_mass"])


def test_stochastic_ensemble_mean_matches_the_smooth_control_in_a_gas_rich_halo():
    """Plenty of gas, weak feedback: no gas limit is reached, so the stochastic mean is the smooth control."""
    cfg = cfg_with(**{"sn.closure": "momentum", "sn.f_p": 0.2})
    Mh, Mdot, z = history(1e8, 5e8, n=20)
    smooth = run_events(Mh, Mdot, z, config=cfg)["stars_mass"][0, -1]
    runs = np.array([run_events(Mh, Mdot, z, config=cfg, sampler="stochastic", seed=s)["stars_mass"][0, -1]
                     for s in range(24)])
    sem = runs.std(ddof=1) / np.sqrt(len(runs))
    assert abs(runs.mean() - smooth) < 4 * sem + 0.01 * smooth


def test_smooth_events_agree_with_the_grid_stepper_with_the_imf_kernel():
    cfg = cfg_with(**{"sn.closure": "momentum", "sn.f_p": 0.7646})
    z = np.linspace(20, 6, 120)
    Mh = 1e9 * np.exp(-0.5 * (z - 6))[None, :]
    Mdot = np.gradient(Mh, utils.time_at_z(z), axis=1)
    events = run_events(Mh, Mdot, z, config=cfg)
    grid = run_forest(Mh, Mdot, z, config=cfg)
    assert events["stars_mass"][0, -1] == pytest.approx(grid["stars_mass"][0, -1], rel=0.02)
    assert events["gas_mass"][0, -1] == pytest.approx(grid["gas_mass"][0, -1], rel=0.06)
    assert events["stars_metals"][0, -1] / events["stars_mass"][0, -1] == pytest.approx(
        grid["stars_metals"][0, -1] / grid["stars_mass"][0, -1], rel=0.05)
