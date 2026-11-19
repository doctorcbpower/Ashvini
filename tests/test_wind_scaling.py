"""The mass and redshift exponents of the Ashvini loading, and the stellar-mass wind model (docs/WIND_SCALING_AUDIT.md)."""
import numpy as np
import pytest

from ashvini import supernovae_feedback as sn_fb, utils
from ashvini.components.wind_model import StellarMassLoading, make_wind_model
from ashvini.config import DEFAULT, override
from ashvini.main import run_forest


def _old_loading(z, Mh, Zs, cfg=DEFAULT):
    """The loading before the exponents were parameters (written out)."""
    return (cfg.sn.epsilon_p * cfg.sn.pi_fid * (10**11.5 / Mh) ** (1 / 3) * (9 / (1 + z)) ** (1 / 2)
            * sn_fb.metallicity_function(Zs))


def test_default_exponents_are_the_old_ones_to_the_last_digit():
    Mh = np.logspace(6, 12, 25)
    for z in (5.0, 9.0, 15.0):
        assert np.array_equal(sn_fb.mass_loading_factor(z, Mh, 0.002), _old_loading(z, Mh, 0.002))


def test_exponents_change_the_loading_as_stated():
    cfg = override(DEFAULT, **{"sn.mass_slope": 2.0 / 3.0, "sn.redshift_slope": 0.0})
    Mh, z = 1e9, 6.0
    ratio = sn_fb.mass_loading_factor(z, Mh, 0.0, cfg) / sn_fb.mass_loading_factor(z, Mh, 0.0)
    assert ratio == pytest.approx((10**11.5 / Mh) ** (1 / 3) * ((1 + z) / 9) ** 0.5, rel=1e-12)


def _stellar(cfg_changes=None):
    cfg = override(DEFAULT, **{"sn.wind_model": "stellar_mass", **(cfg_changes or {})})
    return make_wind_model(cfg)


def test_stellar_mass_loading_is_eq34_of_kravtsov_manwadkar():
    w = _stellar()
    assert isinstance(w, StellarMassLoading)
    Ms = np.array([1e3, 1e6, 1e8, 1e9, 3e9])
    expected = 3.6 * (Ms / 1e10) ** -0.35 - 4.5
    assert np.allclose(w.mass_loading(7.0, np.full(5, 1e10), None, Ms), expected, rtol=1e-12)


def test_cutoff_is_zero_above_the_stellar_mass_where_the_formula_crosses_zero():
    w = _stellar()
    M_cut = 1e10 * (3.6 / 4.5) ** (1 / 0.35)                      # about 4.7e9 Msun
    Mh = np.full(3, 1e11)
    eta = w.mass_loading(7.0, Mh, None, np.array([0.99 * M_cut, M_cut * 1.0001, 1e11]))
    assert eta[0] > 0 and eta[1] == 0.0 and eta[2] == 0.0


def test_maximum_is_respected_and_an_empty_halo_has_no_wind():
    w = _stellar()
    assert w.mass_loading(7.0, np.array([1e8]), None, np.array([0.0]))[0] == 2000.0
    assert w.mass_loading(7.0, np.array([0.0]), None, np.array([1e5]))[0] == 0.0
    assert _stellar({"sn.eta_max": 50.0}).mass_loading(7.0, np.array([1e8]), None, np.array([10.0]))[0] == 50.0


def test_stellar_mass_wind_needs_the_stellar_mass():
    with pytest.raises(ValueError):
        _stellar().mass_loading(7.0, np.array([1e9]), None)


def _run(changes):
    n = 300
    z = np.linspace(20, 6, n)
    Mh = 1e10 * np.exp(-0.5 * (z - 6))[None, :]
    Mdot = np.gradient(Mh, utils.time_at_z(z), axis=1)
    cfg = override(DEFAULT, **{"bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
                               "bh.seeding.halo_mass_threshold.enabled": False, "reion.UVB_enabled": False, **changes})
    return run_forest(Mh, Mdot, z, config=cfg)


def test_run_with_the_stellar_mass_wind_is_finite_and_the_cutoff_matters():
    with_cut = _run({"sn.wind_model": "stellar_mass"})
    no_cut = _run({"sn.wind_model": "stellar_mass", "sn.eta_c": 0.0})
    for r in (with_cut, no_cut):
        assert np.all(np.isfinite(r["gas_mass"])) and np.all(r["gas_mass"] >= 0) and r["stars_mass"][0, -1] > 0
    assert with_cut["stars_mass"][0, -1] > no_cut["stars_mass"][0, -1]       # less wind, more stars


def test_stellar_mass_wind_is_refused_by_the_event_integrator():
    from ashvini.integrators.events import EventIntegrator
    with pytest.raises(ValueError):
        EventIntegrator(override(DEFAULT, **{"sn.wind_model": "stellar_mass"}))


def test_metallicity_step_can_be_removed_and_is_on_by_default():
    assert DEFAULT.sn.metallicity_step is True
    off = override(DEFAULT, **{"sn.metallicity_step": False})
    Mh, z = 1e9, 7.0
    poor = sn_fb.mass_loading_factor(z, Mh, 0.0, DEFAULT)
    assert sn_fb.mass_loading_factor(z, Mh, 0.0, off) == pytest.approx(poor, rel=1e-3)             # metal-poor value is the step's upper level (to the width of the step at Z = 0)
    assert sn_fb.mass_loading_factor(z, Mh, 0.02, off) == pytest.approx(sn_fb.mass_loading_factor(z, Mh, 0.0, off), rel=1e-12)   # no fall with enrichment
    assert sn_fb.mass_loading_factor(z, Mh, 0.02, DEFAULT) == pytest.approx(0.25 * poor, rel=1e-3)
