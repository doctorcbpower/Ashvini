"""Milestone M1: an explicit Config can be passed, is validated, and does not touch the process-wide default."""
import copy
import warnings

import numpy as np
import pytest

from ashvini import main, utils
from ashvini.config import DEFAULT, override, validate, load_config, resolve
from ashvini.state import State

warnings.filterwarnings("ignore")


def _history(n=201, M0=1e9):
    z = np.linspace(20, 5, n)
    t = utils.time_at_z(z)
    Mh = M0 * np.exp(-0.5 * (z - 5))
    return Mh[None], np.gradient(Mh, t)[None], z


def test_default_is_the_run_params_object():
    from ashvini.run_params import PARAMS
    assert DEFAULT is PARAMS and resolve(None) is DEFAULT


def test_override_returns_a_copy_and_leaves_the_default_alone():
    before = copy.deepcopy(DEFAULT)
    cfg = override(DEFAULT, **{"sf.efficiency": 0.05, "bh.seeding.pop3.enabled": False})
    assert cfg.sf.efficiency == 0.05 and cfg.bh.seeding.pop3.enabled is False
    assert DEFAULT.sf.efficiency == before.sf.efficiency
    assert DEFAULT.bh.seeding.pop3.enabled == before.bh.seeding.pop3.enabled


def test_override_rejects_unknown_parameters_and_bad_values():
    with pytest.raises(AttributeError, match="no parameter"):
        override(DEFAULT, **{"sf.efficency": 0.1})
    with pytest.raises(ValueError, match="supernova.type"):
        override(DEFAULT, **{"sn.type": "delayd"})
    with pytest.raises(ValueError):
        override(DEFAULT, **{"sn.m_sn_max": 500.0})
    with pytest.raises(ValueError, match="sigma_feedback"):
        override(DEFAULT, **{"bh.sigma_feedback.growth_cap_enabled": True})


def test_load_config_validates_the_shipped_yaml():
    cfg = load_config()
    assert cfg.sf.efficiency == DEFAULT.sf.efficiency


def test_two_configs_in_one_process_give_different_runs_without_touching_the_default():
    Mh, Md, z = _history()
    low = override(DEFAULT, **{"sf.efficiency": 0.005, "reion.UVB_enabled": False})
    high = override(DEFAULT, **{"sf.efficiency": 0.05, "reion.UVB_enabled": False})
    r_low = main.run_forest(Mh, Md, z, config=low)
    r_high = main.run_forest(Mh, Md, z, config=high)
    assert r_high["stars_mass"][0, -1] > 2 * r_low["stars_mass"][0, -1]  # feedback self-regulates, so less than 10x
    # an explicit config equal to the default reproduces a default call exactly
    r_a = main.run_forest(Mh, Md, z)
    r_b = main.run_forest(Mh, Md, z, config=copy.deepcopy(DEFAULT))
    assert r_a["stars_mass"][0, -1] == r_b["stars_mass"][0, -1]


def test_scalar_and_vectorised_integrators_accept_the_same_config():
    Mh, Md, z = _history(n=60)
    cfg = override(DEFAULT, **{"sn.type": "instantaneous", "reion.UVB_enabled": False})
    vec = main.run_forest(Mh, Md, z, config=cfg)
    scal = main.run1_scalar(Mh[0], Md[0], z, config=cfg)
    assert vec["stars_mass"][0, -1] == pytest.approx(scal["stars_mass"][-1], rel=0.05)


def test_state_dict_has_the_documented_keys_and_shapes():
    st = State.empty(3, 7)
    assert st.shape == (3, 7)
    assert set(st.as_dict()) == {"gas_mass", "gas_metals", "stars_mass", "stars_metals", "dust_mass", "bh_mass", "sfr"}
    Mh, Md, z = _history(n=30)
    res = main.run_forest(Mh, Md, z)
    assert res["gas_mass"].shape == (1, 30) and "stellar_metallicity" not in res
