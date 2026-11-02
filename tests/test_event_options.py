"""M5: the optional features of the event stepper. Each is off by default and each is checked against Stochastic_Occupation's smooth control."""
import json

import numpy as np
import pytest

from tests import stochocc_parity as P

PIN = json.loads(P.OPTION_GOLDEN.read_text()) if P.OPTION_GOLDEN.exists() else {}


def _agree(case, rel_stars=5e-3, rel_ej=5e-3, rel_gas=5e-3):
    # the two codes differ at the 1e-3 level in the plain case (constants, Hubble time table, substep cap)
    a, s = P.ashvini_option(None if case == "plain" else case), PIN[case]
    assert a["M_star"] == pytest.approx(s["M_star"], rel=rel_stars), a
    assert a["M_ej"] == pytest.approx(s["M_ej"], rel=rel_ej), a
    assert a["M_g"] == pytest.approx(s["M_g"], rel=rel_gas, abs=1.0), a


def test_the_plain_minihalo_agrees_with_stochocc():
    _agree("plain")


def test_options_are_all_off_by_default():
    from ashvini.config import DEFAULT
    from ashvini.run_params import EventParams
    assert DEFAULT.events == EventParams()


def test_gate_keeps_a_small_halo_from_forming_stars_until_it_is_hot_enough():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    base = P.base_config(**P.OPTION_BASE["ash"])
    open_ = run_events(M, np.zeros_like(M), z, config=base)
    gated = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.gate_T_vir_K": 5.0e3}))
    first = np.argmax(gated["stars_mass"][0] > 0)
    from ashvini import halo
    assert halo.T_vir_K(M[0, first], z[first]) >= 5.0e3 > halo.T_vir_K(M[0, first - 2], z[first - 2])
    assert np.all(gated["stars_mass"][0, :first] == 0.0)
    assert gated["gas_mass"][0, first - 1] > open_["gas_mass"][0, first - 1]      # the gas waits
    # a mass gate set at that mass opens at about the same time
    mass_gated = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.gate_mass": M[0, first]}))
    assert abs(np.argmax(mass_gated["stars_mass"][0] > 0) - first) <= 1     # the gate is tested at substep starts
    # below gate_switch_z the post-reionisation temperature is also required
    late = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.gate_switch_z": 100.0, "events.gate_T_vir_post_K": 5.0e3}))
    assert late["stars_mass"][0, -1] == pytest.approx(gated["stars_mass"][0, -1], rel=1e-9)


def test_gate_agrees_with_stochocc():
    _agree("gate", 1.5e-2, 1.5e-2, 1.5e-2)    # the time the gate opens moves with the 0.1% difference in T_vir


def test_inflow_lag_delays_the_supply_and_conserves_it():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)[None, :]
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 0.0})
    plain = run_events(M, np.zeros_like(M), z, config=base)
    lagged = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.inflow_kappa": 3.0}))
    # less has arrived at every time; the rest is in transit
    assert np.all(lagged["gas_mass"][0] + lagged["stars_mass"][0] <= plain["gas_mass"][0] + plain["stars_mass"][0] + 1e-6)
    assert lagged["gas_mass"][0, -1] + lagged["stars_mass"][0, -1] < plain["gas_mass"][0, -1] + plain["stars_mass"][0, -1]
    # tau -> 0 recovers the plain run
    fast = run_events(M, np.zeros_like(M), z, config=override(base, **{"events.inflow_kappa": 1e-6}))
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


def test_energy_bank_spends_stored_energy_on_arriving_gas_and_decays():
    from ashvini import utils
    from ashvini.config import override
    from ashvini.constants import Msun_g
    from ashvini.main import run_events
    from ashvini.steppers.events import EventStepper, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    M = P.MINIHALO(z)
    base = P.base_config(**P.OPTION_BASE["ash"])
    stepper = EventStepper(override(base, **{"events.bank_tau_myr": 30.0}))
    hist = HaloHistory(M, z)
    t, dt = 0.5 * (hist.t_min + hist.t_max), 1.0e-3
    v = hist.v_esc_at_t(t)
    res = _Reservoirs(t, 0.0, 0.0)
    res.bank = 3.0e50                                     # enough to unbind less than the arriving gas
    decayed = res.bank * np.exp(-dt / 0.03)
    can_unbind = 2.0 * decayed / v ** 2 / Msun_g
    supply = 10.0 * can_unbind
    through = stepper.lag_and_bank(hist, res, t, dt, supply)
    assert through == pytest.approx(supply - can_unbind, rel=1e-12) and res.ejected == pytest.approx(can_unbind, rel=1e-12)
    assert res.bank == pytest.approx(0.0, abs=1e-6 * decayed)             # all of it spent
    res.bank = 3.0e50                                     # now more energy than the arriving gas needs
    res.ejected = 0.0
    little = 1e-3 * can_unbind
    assert stepper.lag_and_bank(hist, res, t, dt, little) == 0.0 and res.bank > 0.0
    # conservation in a full run, and a momentum closure cannot use a bank
    bank = run_events(M[None, :], np.zeros((1, len(z))), z, config=override(base, **{"events.bank_tau_myr": 30.0}))
    supplied = (utils.Omega_b / utils.Omega_m) * (M[-1] - M[0])
    assert bank["gas_mass"][0, -1] + bank["stars_mass"][0, -1] + bank["ejected_mass"][0, -1] == pytest.approx(supplied, rel=1e-9)
    with pytest.raises(ValueError, match="energy"):
        override(base, **{"sn.closure": "momentum", "events.bank_tau_myr": 30.0})


def test_energy_bank_agrees_with_stochocc():
    # the bank must matter here, or agreement would say nothing: it removes the gas and a fifth of the stars
    assert PIN["energy_bank"]["M_star"] < 0.9 * PIN["energy_bank_control"]["M_star"]
    assert PIN["energy_bank"]["M_g"] < 0.01 * PIN["energy_bank_control"]["M_g"]
    _agree("energy_bank", 1.5e-2, 1.5e-2, 1.5e-2)
    _agree("energy_bank_control", 1.5e-2, 1.5e-2, 5e-2)       # the gas is a small remainder of the supply


def test_access_fraction_caps_what_one_explosion_can_remove():
    from ashvini.config import override
    from ashvini.steppers.events import EventStepper, _Reservoirs
    from ashvini.trajectory import HaloHistory
    z = np.arange(20.0, 10.0 - 1e-9, -0.2)
    base = P.base_config(**{**P.OPTION_BASE["ash"], "sn.eps_fb": 50.0})       # an explosion could unbind far more than exists
    hist = HaloHistory(P.MINIHALO(z), z)
    t = 0.5 * (hist.t_min + hist.t_max)
    for fraction in (1.0, 0.3):
        stepper = EventStepper(override(base, **{"events.access_fraction": fraction}))
        res = _Reservoirs(t, 1.0e5, 0.0)
        stepper.explode(hist, res, t, 20.0, 1.0)
        assert res.ejected == pytest.approx(fraction * 1.0e5, rel=1e-12)
    with pytest.raises(ValueError, match="access_fraction"):
        override(base, **{"events.access_fraction": 0.0})


def test_access_fraction_agrees_with_stochocc():
    assert PIN["access_fraction"]["M_g"] < 0.2 * PIN["energy_bank_control"]["M_g"]       # it matters: the gas is left far lower
    _agree("access_fraction", 3e-2, 1.5e-2, 1.5e-1)    # a cap on each explosion: the gas left is a 1e-3 remainder of the supply
