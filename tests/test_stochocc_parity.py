"""
M3/M4 acceptance: Ashvini reproduces the Stochastic_Occupation smooth control.

Setup and the pinned numbers: tests/stochocc_parity.py (M_h(z) = 1e9 exp(-0.5 (z-5)), z = 20 -> 5, Ashvini's clock, no gate,
no UV, no black hole). The pins in tests/golden/stochocc_smooth.json are Stochastic_Occupation's own output
(WeightedQuadratureSampler, K = 200, lifetimes of Georgy+13 at Z = 0.002); `--parity` reruns it live and checks
that the pins are still what that code gives.
"""
import json

import numpy as np
import pytest

from tests import stochocc_parity as P
from ashvini.main import run_events

PIN = json.loads(P.GOLDEN.read_text()) if P.GOLDEN.exists() else None
pytestmark = pytest.mark.skipif(PIN is None, reason="tests/golden/stochocc_smooth.json not generated")

Z = (15.0, 10.0, 7.0, 5.0)


def _config(kind):
    common = {"sn.kernel": "imf", "sn.m_sn_max": 100.0, "sn.lifetime_Z": 0.002, "sf.efficiency": 0.015}
    if kind == "no_feedback":
        return P.base_config(**{**common, "sn.closure": "momentum", "sn.f_p": 0.0})
    return P.base_config(**{**common, "sn.closure": "momentum", "sn.f_p": 0.7625})


def _compare(a, case):
    pin = PIN[case]
    stars = np.array([a["stars_mass"][i] / pin[str(z)]["M_star"] for i, z in enumerate(Z)])
    gas = np.array([a["gas_mass"][i] / pin[str(z)]["M_g"] for i, z in enumerate(Z)])
    if case == "momentum":      # no chemistry in the pinned momentum case
        return stars, gas, np.ones(len(Z)), np.ones(len(Z))
    Zs = np.array([(a["stars_metals"][i] / a["stars_mass"][i]) / (pin[str(z)]["M_Zs"] / pin[str(z)]["M_star"]) for i, z in enumerate(Z)])
    Zg = np.array([(a["gas_metals"][i] / a["gas_mass"][i]) / (pin[str(z)]["M_Zg"] / pin[str(z)]["M_g"]) for i, z in enumerate(Z)])
    return stars, gas, Zs, Zg


def test_grid_stepper_without_feedback_reproduces_stars_to_a_hundredth_of_a_percent():
    """A1: feedback off. Stars agree to 1e-4 (the first-order lag in stars formed was fixed in the chemistry fixes)."""
    a, _, _ = P.ashvini(_config("no_feedback"))
    stars, gas, Zs, Zg = _compare(a, "no_feedback")
    assert np.allclose(stars, 1.0, atol=1e-4), stars
    assert np.allclose(gas, 1.0, atol=1e-3), gas


def test_grid_stepper_without_feedback_reproduces_metallicities():
    """B1: feedback off, chemistry on. The acceptance criterion was 2 percent; the agreement is 1e-4."""
    a, _, _ = P.ashvini(_config("no_feedback"))
    _, _, Zs, Zg = _compare(a, "no_feedback")
    assert np.allclose(Zs, 1.0, atol=5e-4), Zs
    assert np.allclose(Zg, 1.0, atol=5e-4), Zg


def test_grid_stepper_with_the_imf_kernel_and_momentum_closure_reproduces_the_smooth_control():
    """
    E1: the acceptance criterion was stars within 1-4 percent; the delay kernel and the closure agree to 0.06
    percent at z = 15-5 (the pin tolerance below is 0.2 percent).
    """
    a, _, _ = P.ashvini(_config("momentum"))
    stars, gas, Zs, Zg = _compare(a, "momentum")
    assert np.allclose(stars, 1.0, atol=2e-3), stars
    assert np.allclose(gas, 1.0, atol=2e-3), gas


@pytest.mark.parametrize("case", ["no_feedback", "momentum"])
def test_event_stepper_reproduces_the_smooth_control(case):
    cfg = _config(case)
    z = np.arange(20.0, 5.0 - 1e-9, -0.1)
    Mh = P.Mh_of_z(z)[None, :]
    r = run_events(Mh, np.zeros_like(Mh), z, config=cfg, sampler="smooth")
    idx = [int(np.argmin(abs(z - zz))) for zz in Z]
    a = {k: r[k][0, idx] for k in ("gas_mass", "stars_mass", "gas_metals", "stars_metals")}
    stars, gas, Zs, Zg = _compare(a, case)
    assert np.allclose(stars, 1.0, atol=3e-3), stars       # observed: 1e-3 at most
    assert np.allclose(gas, 1.0, atol=3e-3), gas
    assert np.allclose(Zs, 1.0, atol=2e-3), Zs             # observed: 1e-4 at most
    assert np.allclose(Zg, 1.0, atol=2e-3), Zg             # observed: 5e-4 at most


@pytest.mark.parity
@pytest.mark.skipif(not P.available(), reason="Stochastic_Occupation not found")
@pytest.mark.parametrize("case", ["no_feedback", "momentum"])
def test_pins_are_still_what_stochocc_gives(case):
    live = P.stochocc(zout=(10.0,), **P.CASES[case])[0]
    for key in ("M_g", "M_star", "M_Zg", "M_Zs", "M_d"):
        assert getattr(live, key) == pytest.approx(PIN[case]["10.0"][key], rel=1e-9)


@pytest.mark.skipif(not P.available(), reason="Stochastic_Occupation not found")
@pytest.mark.parametrize("closure", [{"sn.closure": "energy", "sn.eps_fb": 0.1}, {"sn.closure": "momentum", "sn.f_p": 0.7625}])
def test_stochastic_ensembles_agree_within_sampling_error(closure):
    """
    A minihalo growing from 4e5 to 2e7 Msun, where explosions reach the gas limit: the mean stars formed and gas removed
    over 200 stochastic realisations per code agree within four standard errors of the difference. (Different random
    streams, so only the statistics can agree.) Stochastic_Occupation is run live; it takes a few seconds.
    """
    Mh = lambda z: 2e7 * np.exp(-0.4 * (np.asarray(z) - 10.0))
    cfg = P.base_config(**{"sn.kernel": "imf", "sn.m_sn_max": 100.0, "sn.lifetime_Z": 0.002, "sf.efficiency": 0.015, **closure})
    z = np.arange(20.0, 10.0 - 1e-9, -0.1)
    M = Mh(z)[None, :]
    N = 200
    ash = [run_events(M, np.zeros_like(M), z, config=cfg, sampler="stochastic", seed=s) for s in range(N)]
    kw = dict(eps_fb=0.1) if "sn.eps_fb" in closure else dict(eps_fb=0.0, closure="momentum", f_p=0.7625)
    sto = [P.stochocc(zout=(10.0,), Mh=Mh, sampler="stochastic", seed=1000 + s, **kw)[0] for s in range(N)]
    for a_vals, s_vals in (([r["stars_mass"][0, -1] for r in ash], [s.M_star for s in sto]),
                           ([r["ejected_mass"][0, -1] for r in ash], [s.M_ej for s in sto])):
        a, s = np.array(a_vals), np.array(s_vals)
        se = np.hypot(a.std(ddof=1), s.std(ddof=1)) / np.sqrt(N)
        assert abs(a.mean() - s.mean()) < 4 * se, (a.mean(), s.mean(), se)
        assert a.std(ddof=1) == pytest.approx(s.std(ddof=1), rel=0.35)        # the scatter, too
