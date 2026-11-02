"""M5: the optional features of the event integrator. Each is off by default and each is checked against Stochastic_Occupation's smooth control."""
import json

import numpy as np
import pytest

from tests import stochocc_comparison as P

PIN = json.loads(P.OPTION_GOLDEN.read_text()) if P.OPTION_GOLDEN.exists() else {}


def _agree(case, rel_stars=5e-3, rel_ej=5e-3, rel_gas=5e-3):
    # the two codes differ at the 1e-3 level in the plain case (constants, Hubble time table, substep cap)
    a, s = P.ashvini_option(None if case == "plain" else case), PIN[case]
    assert a["M_star"] == pytest.approx(s["M_star"], rel=rel_stars), a
    assert a["M_ej"] == pytest.approx(s["M_ej"], rel=rel_ej), a
    assert a["M_g"] == pytest.approx(s["M_g"], rel=rel_gas, abs=100.0), a   # abs: a remainder under 1e-4 of the 3e6 Msun supplied


def test_the_plain_minihalo_agrees_with_stochocc():
    _agree("plain")


def test_options_are_all_off_by_default():
    from ashvini.config import DEFAULT
    from ashvini.run_params import EventParams
    assert DEFAULT.events == EventParams()


def test_cooling_threshold_keeps_a_small_halo_from_forming_stars_until_it_is_hot_enough():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    base = P.base_config(**P.OPTION_BASE["ash"])
    without = run_events(M, np.zeros_like(M), z, config=base)
    with_threshold = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.cooling_T_vir_K": 5.0e3}))
    first = np.argmax(with_threshold["stars_mass"][0] > 0)
    from ashvini import halo
    assert halo.T_vir_K(M[0, first], z[first]) >= 5.0e3 > halo.T_vir_K(M[0, first - 2], z[first - 2])
    assert np.all(with_threshold["stars_mass"][0, :first] == 0.0)
    assert with_threshold["gas_mass"][0, first - 1] > without["gas_mass"][0, first - 1]      # the gas waits
    # a cooling mass set at that mass lets stars form at about the same time
    mass_with_threshold = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.cooling_mass": M[0, first]}))
    assert abs(np.argmax(mass_with_threshold["stars_mass"][0] > 0) - first) <= 1     # the threshold is tested at the start of each substep
    # below reionised_below_z the post-reionisation temperature is also required
    late = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.reionised_below_z": 100.0, "events.cooling_T_vir_reionised_K": 5.0e3}))
    assert late["stars_mass"][0, -1] == pytest.approx(with_threshold["stars_mass"][0, -1], rel=1e-9)


def test_cooling_threshold_agrees_with_stochocc():
    _agree("cooling_threshold", 1.5e-2, 1.5e-2, 1.5e-2)    # the time the threshold is reached moves with the 0.1% difference in T_vir


def test_inflow_lag_delays_the_supply_and_conserves_it():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 0.0})
    plain = run_events(M, np.zeros_like(M), z, config=base)
    lagged = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.inflow_delay_in_crossing_times": 3.0}))
    # less has arrived at every time; the rest is in transit
    assert np.all(lagged["gas_mass"][0] + lagged["stars_mass"][0] <= plain["gas_mass"][0] + plain["stars_mass"][0] + 1e-6)
    assert lagged["gas_mass"][0, -1] + lagged["stars_mass"][0, -1] < plain["gas_mass"][0, -1] + plain["stars_mass"][0, -1]
    # tau -> 0 recovers the plain run
    fast = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.inflow_delay_in_crossing_times": 1e-6}))
    assert fast["stars_mass"][0, -1] == pytest.approx(plain["stars_mass"][0, -1], rel=1e-3)


def test_inflow_lag_agrees_with_stochocc():
    _agree("inflow_lag", 1.5e-2, 1.5e-2, 1.5e-2)


def test_recycling_returns_gas_conserves_baryons_and_metals():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    base = P.base_config(**P.OPTION_BASE["ash"])
    plain = run_events(M, np.zeros_like(M), z, config=base)
    rec = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.recycling_fraction": 0.5}))
    assert plain["returned_mass"][0, -1] == 0.0 and rec["returned_mass"][0, -1] > 0.0
    assert rec["returned_mass"][0, -1] <= 0.5 * rec["ejected_mass"][0, -1] * (1 + 1e-12)
    assert rec["stars_mass"][0, -1] > plain["stars_mass"][0, -1]                 # the returned gas forms more stars
    supplied = (utils.Omega_b / utils.Omega_m) * (M[0, -1] - M[0, 0])
    net_removed = rec["ejected_mass"][0, -1] - rec["returned_mass"][0, -1]
    assert rec["gas_mass"][0, -1] + rec["stars_mass"][0, -1] + net_removed == pytest.approx(supplied, rel=1e-9)
    # nothing returns if the fraction is zero, and a return of metals alone does not move the gas
    only_metals = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.recycling_fraction_Z": 1.0}))
    assert only_metals["returned_mass"][0, -1] == 0.0
    assert only_metals["gas_metals"][0, -1] > plain["gas_metals"][0, -1]


def test_recycling_agrees_with_stochocc():
    _agree("recycling", 1.5e-2, 1.5e-2, 1.5e-2)


def test_stored_energy_spends_stored_energy_on_arriving_gas_and_decays():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.constants import Msun_g
    from ashvini.main import run_events
    from ashvini.integrators.events import EventIntegrator, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)
    base = P.base_config(**P.OPTION_BASE["ash"])
    integrator = EventIntegrator(override(base, **{"events.stored_energy_loss_time_myr": 30.0}))
    hist = HaloHistory(M, z)
    t, dt = 0.5 * (hist.t_min + hist.t_max), 1.0e-3
    v = hist.v_esc_at_t(t)
    res = _Reservoirs(t, 0.0, 0.0)
    res.bank = 3.0e50                                     # enough to unbind less than the arriving gas
    decayed = res.bank * np.exp(-dt / 0.03)
    can_unbind = 2.0 * decayed / v ** 2 / Msun_g
    supply = 10.0 * can_unbind
    through = integrator.lag_and_bank(hist, res, t, dt, supply)
    assert through == pytest.approx(supply - can_unbind, rel=1e-12) and res.ejected == pytest.approx(can_unbind, rel=1e-12)
    assert res.bank == pytest.approx(0.0, abs=1e-6 * decayed)             # all of it spent
    res.bank = 3.0e50                                     # now more energy than the arriving gas needs
    res.ejected = 0.0
    little = 1e-3 * can_unbind
    assert integrator.lag_and_bank(hist, res, t, dt, little) == 0.0 and res.bank > 0.0
    # conservation in a full run, and a momentum wind_model cannot use a bank
    bank = run_events(M[None, :], np.zeros((1, len(z))), z, config=override(base, **{"events.stored_energy_loss_time_myr": 30.0}))
    supplied = (utils.Omega_b / utils.Omega_m) * (M[-1] - M[0])
    assert bank["gas_mass"][0, -1] + bank["stars_mass"][0, -1] + bank["ejected_mass"][0, -1] == pytest.approx(supplied, rel=1e-9)
    with pytest.raises(ValueError, match="energy"):
        override(base, **{"sn.wind_model": "momentum", "events.stored_energy_loss_time_myr": 30.0})


def test_stored_energy_agrees_with_stochocc():
    # the bank must matter here, or agreement would say nothing: it removes the gas and a fifth of the stars
    assert PIN["stored_energy"]["M_star"] < 0.9 * PIN["stored_energy_control"]["M_star"]
    assert PIN["stored_energy"]["M_g"] < 0.01 * PIN["stored_energy_control"]["M_g"]
    _agree("stored_energy", 1.5e-2, 1.5e-2, 1.5e-2)
    _agree("stored_energy_control", 1.5e-2, 1.5e-2, 5e-2)       # the gas is a small remainder of the supply


def test_reachable_gas_fraction_caps_what_one_explosion_can_remove():
    from ashvini.config import override
    from ashvini.integrators.events import EventIntegrator, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 50.0})       # an explosion could unbind far more than exists
    hist = HaloHistory(P.MINIHALO(z), z)
    t = 0.5 * (hist.t_min + hist.t_max)
    for fraction in (1.0, 0.3):
        integrator = EventIntegrator(override(base, **{"events.reachable_gas_fraction": fraction}))
        res = _Reservoirs(t, 1.0e5, 0.0)
        integrator.explode(hist, res, t, 20.0, 1.0)
        assert res.ejected == pytest.approx(fraction * 1.0e5, rel=1e-12)
    with pytest.raises(ValueError, match="reachable_gas_fraction"):
        override(base, **{"events.reachable_gas_fraction": 0.0})


def test_reachable_gas_fraction_agrees_with_stochocc():
    assert PIN["reachable_gas_fraction"]["M_g"] < 0.2 * PIN["stored_energy_control"]["M_g"]       # it matters: the gas is left far lower
    _agree("reachable_gas_fraction", 3e-2, 1.5e-2, 1.5e-1)    # a cap on each explosion: the gas left is a 1e-3 remainder of the supply


def test_coupling_change_uses_the_high_mass_coupling_at_the_explosion_time():
    from ashvini.config import override
    from ashvini.integrators.events import EventIntegrator, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    hist = HaloHistory(P.MINIHALO(z), z)
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 1.0})
    cfg = override(base, **{"events.coupling_change_mass": 1.0e7, "events.eps_fb_above_change_mass": 0.05})
    t_below, t_above = hist.t_min + 0.2 * (hist.t_max - hist.t_min), hist.t_max
    assert hist.mass_at_t(t_below) < 1.0e7 <= hist.mass_at_t(t_above)
    for t, eps in ((t_below, 1.0), (t_above, 0.05)):
        res = _Reservoirs(t, 1.0e12, 0.0)
        EventIntegrator(cfg).explode(hist, res, t, 20.0, 1.0)
        reference = _Reservoirs(t, 1.0e12, 0.0)
        EventIntegrator(override(base, **{"sn.eps_fb": eps})).explode(hist, reference, t, 20.0, 1.0)
        assert res.ejected == pytest.approx(reference.ejected, rel=1e-12)
    with pytest.raises(ValueError, match="energy"):
        override(base, **{"sn.wind_model": "momentum", "events.coupling_change_mass": 1.0e7})


def test_coupling_change_agrees_with_stochocc():
    assert PIN["coupling_change"]["M_star"] > 1.5 * PIN["stored_energy_control"]["M_star"]    # it matters: weaker coupling above the cut
    _agree("coupling_change", 1.5e-2, 1.5e-2, 5e-2)


def test_wind_travel_time_fixes_the_mass_at_the_explosion_and_removes_it_later():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    from ashvini.integrators.events import EventIntegrator, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    hist = HaloHistory(P.MINIHALO(z), z)
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 1.0})
    cfg = override(base, **{"events.wind_delay_in_crossing_times": 2.0})
    t = 0.5 * (hist.t_min + hist.t_max)
    integrator = EventIntegrator(cfg)
    res = _Reservoirs(t, 1.0e5, 0.0)
    integrator.explode(hist, res, t, 20.0, 1.0)
    assert res.ejected == 0.0 and res.n_explosions == 1 and len(res.queue) == 1     # nothing yet, one delivery queued
    t_delivery, _, mass, weight, ret, _, _, eject = res.queue[0]
    assert t_delivery > t and eject > 0 and ret is None
    integrator.explode(hist, res, t_delivery, mass, weight, eject=eject)
    assert res.ejected == pytest.approx(min(1.0e5, eject), rel=1e-12) and res.n_explosions == 1   # not counted twice
    # in a full run the baryons are conserved, and a longer travel time keeps more gas for star formation meanwhile
    M = P.MINIHALO(z)[None, :]
    supplied = (utils.Omega_b / utils.Omega_m) * (M[0, -1] - M[0, 0])
    r = run_events(M, np.zeros_like(M), z, config=cfg)
    assert r["gas_mass"][0, -1] + r["stars_mass"][0, -1] + r["ejected_mass"][0, -1] == pytest.approx(supplied, rel=1e-9)


def test_wind_travel_time_agrees_with_stochocc():
    assert PIN["wind_travel"]["M_star"] > 1.1 * PIN["stored_energy_control"]["M_star"]      # it matters: the gas lingers
    _agree("wind_travel", 1.5e-2, 1.5e-2, 5e-2)


def test_all_options_together_run_and_conserve_baryons():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    cfg = override(P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 0.5}), **{
        "events.cooling_T_vir_K": 3.0e3, "events.inflow_delay_in_crossing_times": 1.0, "events.recycling_fraction": 0.3,
        "events.stored_energy_loss_time_myr": 50.0, "events.reachable_gas_fraction": 0.5, "events.coupling_change_mass": 1.0e7,
        "events.eps_fb_above_change_mass": 0.1, "events.wind_delay_in_crossing_times": 1.0})
    for sampler in ("smooth", "stochastic"):
        r = run_events(M, np.zeros_like(M), z, config=cfg, sampler=sampler, seed=2)
        arrived = (utils.Omega_b / utils.Omega_m) * (M[0, -1] - M[0, 0])
        net = r["ejected_mass"][0, -1] - r["returned_mass"][0, -1]
        # gas still on its way to the reservoir is the only baryonic mass not in gas, stars or the net removed
        in_transit = arrived - (r["gas_mass"][0, -1] + r["stars_mass"][0, -1] + net)
        assert 0.0 <= in_transit < 0.5 * arrived
        assert np.all(r["gas_mass"] >= 0) and np.all(r["dust_mass"] >= 0) and np.all(np.isfinite(r["stars_metals"]))
