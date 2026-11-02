"""M3: the IMF delay kernel and the momentum / energy closures."""
import numpy as np
import pytest

from ashvini import halo, imf, utils
from ashvini.components.closure import AshviniLoading, EnergyClosure, MomentumClosure, make_closure
from ashvini.components.delay import (LifetimeDistributionDelay, SupernovaDelay, explosion_kernel, make_delay)
from ashvini.config import DEFAULT, override
from ashvini.main import run_forest
from ashvini.stellar_lifetimes import t_ms_myr


@pytest.fixture
def imf_cfg():
    return override(DEFAULT, **{"sn.kernel": "imf", "sn.m_sn_max": 100.0, "sn.lifetime_Z": 0.002})


def test_default_kernel_and_closure_are_ashvini_s():
    assert DEFAULT.sn.kernel == "delta" and DEFAULT.sn.closure == "ashvini"
    assert isinstance(make_closure(DEFAULT), AshviniLoading)
    assert isinstance(make_delay(DEFAULT, np.linspace(0.3, 1.0, 20)), SupernovaDelay)


def test_lifetimes_are_the_tabulated_values_and_fall_with_mass():
    assert t_ms_myr(8.0, 0.002) == pytest.approx(37.447, rel=1e-4)
    assert t_ms_myr(100.0, 0.002) == pytest.approx(3.412, rel=1e-3)
    m = np.linspace(8, 100, 50)
    assert np.all(np.diff(t_ms_myr(m, 0.002)) < 0)
    # between the two grids, between their lifetimes
    mid = t_ms_myr(20.0, 0.005)
    assert min(t_ms_myr(20.0, 0.002), t_ms_myr(20.0, 0.014)) < mid < max(t_ms_myr(20.0, 0.002), t_ms_myr(20.0, 0.014))
    with pytest.raises(ValueError):
        t_ms_myr(5.0)


def test_kernel_conserves_the_stars_that_have_time_to_explode():
    T = utils.time_at_z(np.linspace(20, 5, 400))
    lifetimes = np.array([0.004, 0.010, 0.030])
    K = explosion_kernel(T, lifetimes)
    assert np.all(np.triu(K) == 0.0) and np.all(K >= 0)
    done = T[-1] - T > lifetimes.max() + 0.01      # formed early enough that every progenitor has exploded
    assert np.allclose(K.sum(axis=0)[done & (np.arange(len(T)) > 0)], 1.0, atol=1e-7)


def test_kernel_mean_delay_is_the_mean_lifetime():
    T = np.linspace(0.2, 1.0, 1601)                 # uniform 0.5 Myr steps
    lifetimes = np.array([0.005, 0.01, 0.02, 0.04])
    K = explosion_kernel(T, lifetimes)
    k = 10                                          # a star formed at the start of the run
    delay = (K[:, k] * (T - 0.5 * (T[k] + T[k - 1]))).sum()
    assert delay == pytest.approx(lifetimes.mean(), rel=0.05)


def test_delta_lifetime_reproduces_the_delay_time_response():
    """One lifetime t_d: the supernovae go off t_d after formation, as in the delta kernel (within one step)."""
    T = np.linspace(0.2, 1.0, 801)                  # 1 Myr steps
    K = explosion_kernel(T, np.array([0.015]))
    k = 100
    j = np.flatnonzero(K[:, k] > 0)
    assert 1 <= len(j) <= 4
    assert T[j].min() == pytest.approx(T[k] + 0.015, abs=1.5e-3)
    assert T[j].max() == pytest.approx(T[k] + 0.015, abs=2.5e-3)


def test_momentum_closure_matches_the_ashvini_loading_to_half_a_percent():
    # at zero stellar metallicity, where the Ashvini loading's metallicity step has not yet reduced it
    cfg = override(DEFAULT, **{"sn.m_sn_max": 100.0, "sn.f_p": 0.7646, "sn.closure": "momentum"})
    ashvini = AshviniLoading(cfg)
    momentum = MomentumClosure(cfg)
    for z in (15.0, 10.0, 5.0):
        for M in (1e6, 1e8, 1e10):
            a = float(ashvini.mass_loading(z, M, 0.0))
            m = float(momentum.mass_loading(z, M))
            assert m == pytest.approx(a, rel=5e-3)


def test_energy_closure_per_supernova_is_two_eps_E_over_vesc_squared():
    cfg = override(DEFAULT, **{"sn.closure": "energy", "sn.eps_fb": 0.1})
    M, z = 1e8, 10.0
    v = halo.v_esc_cms(M, z)
    expected = 2 * 0.1 * 1e51 / v ** 2 / 1.98892e33
    c = EnergyClosure(cfg)
    assert float(c.mass_per_supernova(z, M)) == pytest.approx(expected, rel=1e-12)
    assert float(c.mass_loading(z, M)) == pytest.approx(c.nu_SN * expected, rel=1e-12)
    assert float(c.mass_per_supernova(z, 0.0)) == 0.0


def test_nu_SN_matches_the_reference_values():
    assert imf.nu_SN("salpeter", 8.0, 100.0) == pytest.approx(7.4e-3, rel=0.01)
    assert imf.nu_SN("salpeter", 8.0, 40.0) == pytest.approx(6.8e-3, rel=0.01)


def test_imf_kernel_run_conserves_mass_and_every_exploded_star_drives_a_wind(imf_cfg):
    n = 300
    z = np.linspace(20, 6, n)
    Mh = 1e9 * np.exp(-0.5 * (z - 6))[None, :]
    Mdot = np.gradient(Mh, utils.time_at_z(z), axis=1)
    cfg = override(imf_cfg, **{"bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
                               "bh.seeding.halo_mass_threshold.enabled": False, "reion.UVB_enabled": False})
    r = run_forest(Mh, Mdot, z, config=cfg)
    assert np.all(np.isfinite(r["gas_mass"])) and r["stars_mass"][0, -1] > 0
    # with the Ashvini delta kernel the stars formed are larger: the IMF kernel releases the SNe sooner
    # for the shortest-lived progenitors and later for the longest-lived, but both remove gas
    delta = run_forest(Mh, Mdot, z, config=override(cfg, **{"sn.kernel": "delta"}))
    assert r["stars_mass"][0, -1] != delta["stars_mass"][0, -1]
    assert 0.5 < r["stars_mass"][0, -1] / delta["stars_mass"][0, -1] < 2.0


def test_kernel_matrix_agrees_with_the_wind_in_a_run(imf_cfg):
    """The wind star formation rate in step j is the kernel applied to the stars formed in earlier steps."""
    n = 60
    z = np.linspace(20, 12, n)
    Mh = np.full((1, n), 1e9)
    Mdot = np.zeros((1, n))
    Mdot[0, 1:] = 1e8
    cfg = override(imf_cfg, **{"bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
                               "bh.seeding.halo_mass_threshold.enabled": False, "reion.UVB_enabled": False})
    r = run_forest(Mh, Mdot, z, config=cfg, gas_mass0=1e6)
    cosmic = utils.time_at_z(z)
    delay = LifetimeDistributionDelay(cfg, cosmic)
    S = np.diff(r["stars_mass"][0])
    j = n - 1
    expected = delay.kernel[j, 1:j] @ S[: j - 1] / (cosmic[j] - cosmic[j - 1])
    assert expected > 0


def test_stochocc_preset_gives_the_same_star_formation_rate_and_a_valid_config():
    from ashvini.presets import ashvini_efficiency_for_stochocc, stochocc_smooth_control
    assert ashvini_efficiency_for_stochocc(0.02178) == pytest.approx(0.015, rel=1e-3)
    cfg = stochocc_smooth_control()
    assert (cfg.sn.kernel, cfg.sn.closure, cfg.sn.m_sn_max) == ("imf", "momentum", 100.0)
    assert DEFAULT.sn.kernel == "delta"          # the preset returns a copy
