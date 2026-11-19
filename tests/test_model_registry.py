"""The knob registry and Model (docs/INTERFACE_DESIGN.md, step 1)."""
import numpy as np
import pytest

from ashvini.config import DEFAULT, override
from ashvini.model import Model, PRESETS
from ashvini.parameters import KNOBS, knob_range, parameters, _get


def test_every_knob_path_exists_and_default_is_a_number():
    for name, k in KNOBS.items():
        for sel, path in k.paths.items():
            cfg = DEFAULT if k.selector is None else override(DEFAULT, **{k.selector: sel})
            assert isinstance(_get(cfg, path), (int, float)), (name, path)


def test_default_value_lies_inside_each_knobs_range():
    for name in KNOBS:
        for preset in PRESETS:
            cfg = override(DEFAULT, **PRESETS[preset])
            if KNOBS[name].path(cfg) is None:
                continue
            lo, hi, _ = knob_range(name, cfg)
            assert lo <= KNOBS[name].default(cfg) <= hi, (name, preset)


def test_wind_loading_sets_the_field_of_the_active_wind_model():
    for preset, path in [("fiducial", "sn.epsilon_p"), ("momentum_wind", "sn.f_p"), ("energy_wind", "sn.eps_fb"), ("stellar_mass_wind", "sn.eta_norm")]:
        m = Model(preset, free=["wind_loading"])
        assert m._paths["wind_loading"] == path
        assert _get(m.config([0.7]), path) == 0.7


def test_bh_knobs_follow_the_growth_model():
    a = Model("fiducial", free=["bh_accretion_eff", "bh_feedback_coupling"])
    b = Model("slim_disc", free=["bh_accretion_eff", "bh_feedback_coupling"])
    assert list(a._paths.values()) == ["bh.efficiency", "bh.eta_agn"]
    assert list(b._paths.values()) == ["bh.slimdisk.eta_acc", "bh.slimdisk.epsilon_f"]


def test_config_at_default_theta_is_the_preset_and_does_not_touch_DEFAULT():
    before = DEFAULT.sf.efficiency
    m = Model("fiducial", free=["sf_efficiency", "z_reion"])
    cfg = m.config(None)
    assert cfg.sf.efficiency == DEFAULT.sf.efficiency and cfg.reion.z_reion == DEFAULT.reion.z_reion
    m.config([0.05, 8.0])
    assert DEFAULT.sf.efficiency == before


def test_theta_by_name_or_array_agree_and_wrong_names_are_refused():
    m = Model("fiducial", free=["sf_efficiency", "z_reion"])
    assert m.config({"z_reion": 8.0, "sf_efficiency": 0.05}).reion.z_reion == m.config([0.05, 8.0]).reion.z_reion == 8.0
    with pytest.raises(KeyError, match="missing"):
        m.config({"sf_efficiency": 0.05})
    with pytest.raises(ValueError, match="shape"):
        m.config([0.05])


def test_inapplicable_knobs_are_refused_with_the_reason():
    with pytest.raises(ValueError, match="needs sn.wind_model = 'ashvini'"):
        Model("energy_wind", free=["wind_mass_slope"])
    with pytest.raises(ValueError, match="would change nothing"):
        Model("fiducial", free=["bh_seed_mass_pop3"], fixed={"bh.seeding.pop3.enabled": False})
    Model("fiducial", free=["bh_seed_mass_pop3"])


def test_unknown_knob_suggests_the_closest():
    with pytest.raises(KeyError, match="sf_efficiency"):
        Model("fiducial", free=["sf_efficency"])
    with pytest.raises(KeyError, match="unknown preset"):
        Model("energy", free=[])


def test_fixed_and_free_cannot_overlap_and_fixed_knobs_are_applied():
    with pytest.raises(ValueError, match="both free and fixed"):
        Model("fiducial", free=["z_reion"], fixed={"z_reion": 7.5})
    m = Model("fiducial", free=["sf_efficiency"], fixed={"z_reion": 7.5, "sn.delay_time": 0.02})
    cfg = m.config(None)
    assert cfg.reion.z_reion == 7.5 and cfg.sn.delay_time == 0.02


def test_unit_cube_maps_to_the_bounds_in_log_or_linear():
    m = Model("fiducial", free=["sf_efficiency", "z_reion"])
    lo, hi = m.bounds[:, 0], m.bounds[:, 1]
    assert np.allclose(m.unit_to_theta([0, 0]), lo) and np.allclose(m.unit_to_theta([1, 1]), hi)
    mid = m.unit_to_theta([0.5, 0.5])
    assert mid[0] == pytest.approx(np.sqrt(lo[0] * hi[0])) and mid[1] == pytest.approx(0.5 * (lo[1] + hi[1]))
    assert m.contains(mid) and not m.contains([lo[0] / 2, mid[1]])


def test_bounds_can_be_overridden_and_must_be_ordered():
    assert Model("fiducial", free=["z_reion"], bounds={"z_reion": (7.0, 8.0)}).bounds.tolist() == [[7.0, 8.0]]
    with pytest.raises(ValueError, match="low < high"):
        Model("fiducial", free=["z_reion"], bounds={"z_reion": (8.0, 7.0)})


def test_run_matches_run_forest_with_the_same_config():
    from ashvini.main import run_forest
    z = np.linspace(10, 0, 40)
    mh = np.tile(np.logspace(9, 10, 40), (3, 1)) * np.array([[1.0], [2.0], [4.0]])
    rate = np.gradient(mh, axis=1) / 0.1
    m = Model("fiducial", free=["sf_efficiency"])
    a = m.run(mh, rate, z, theta=[0.04])
    b = run_forest(mh, rate, z, config=override(DEFAULT, **{"sf.efficiency": 0.04}))
    assert np.array_equal(a["stars_mass"], b["stars_mass"])


def test_parameters_table_lists_every_knob():
    text = parameters(show=False)
    assert all(n in text for n in KNOBS)
    assert "not used with the current set-up" in parameters(override(DEFAULT, **{"sn.wind_model": "energy"}), show=False)


def test_accretion_efficiency_range_reaches_a_tenth_of_a_per_cent_black_hole_to_stellar_mass_ratio():
    lo, hi, _ = knob_range("bh_accretion_eff", DEFAULT)
    assert lo <= 1e-5 and hi == 1e-2          # the fiducial model needs about 2e-5 for M_BH/M* of 0.1 per cent (docs/FIT_AUDIT.md, section 10)


def test_eddington_multiplier_knob_is_for_the_default_growth_model_only():
    m = Model("fiducial", free=["bh_eddington_multiplier"])
    assert m._paths["bh_eddington_multiplier"] == "bh.eddington_multiplier" and m.config([0.5]).bh.eddington_multiplier == 0.5
    with pytest.raises(ValueError, match="needs bh.growth_model = 'pznk11_freefall'"):
        Model("slim_disc", free=["bh_eddington_multiplier"])


def test_msigma_baryon_fraction_needs_the_growth_cap_and_the_preset_provides_it():
    with pytest.raises(ValueError, match="would change nothing"):
        Model("fiducial", free=["bh_msigma_baryon_fraction"])
    m = Model("m_sigma_cap", free=["bh_msigma_baryon_fraction"])
    cfg = m.config([0.1])
    assert cfg.bh.sigma_feedback.f_g == 0.1 and cfg.bh.sigma_feedback.enabled and cfg.bh.sigma_feedback.growth_cap_enabled
    assert m.theta0[0] == 0.16
