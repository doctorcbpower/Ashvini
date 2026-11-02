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
    assert a["M_g"] == pytest.approx(s["M_g"], rel=rel_gas), a


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
