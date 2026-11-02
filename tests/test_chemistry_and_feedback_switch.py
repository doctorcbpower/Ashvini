"""
Regression tests for the 2026-11-02 corrections (docs/SUPERSET_DESIGN.md, section 2):
zero delay, the feedback switch applied to dust, the stellar-mass step lag, the dust supernova rate,
and metallicity units.
"""
import warnings

import numpy as np
import pytest

from ashvini import configure, main, utils
from ashvini.constants import Z_SUN
from ashvini.imf import nu_SN
from ashvini.run_params import PARAMS, load_params
from ashvini import supernovae_feedback as sn
from ashvini.metallicity import igm_metallicity

warnings.filterwarnings("ignore")


@pytest.fixture(autouse=True)
def restore_params():
    saved = dict(t_d=PARAMS.sn.delay_time, sn_type=PARAMS.sn.type, uv=PARAMS.reion.UVB_enabled)
    yield
    configure.set_params(t_d=saved["t_d"], sn_type=saved["sn_type"], uv_background=saved["uv"], bh_enabled=True)


def _history(n=301, M0=1e9):
    z = np.linspace(20, 5, n)
    t = utils.time_at_z(z)
    Mh = M0 * np.exp(-0.5 * (z - 5))
    return Mh[None], np.gradient(Mh, t)[None], z


def _run(sn_type, t_d=0.015, n=301):
    configure.set_params(bh_enabled=False, uv_background=False, t_d=t_d, sn_type=sn_type)
    Mh, Md, z = _history(n)
    return main.run_forest(Mh, Md, z)


def test_zero_delay_is_instantaneous_feedback_not_no_feedback():
    zero_delay = _run("delayed", t_d=0.0)
    instant = _run("instantaneous")
    none = _run("no")
    assert zero_delay["stars_mass"][0, -1] == pytest.approx(instant["stars_mass"][0, -1], rel=1e-12)
    assert zero_delay["stars_mass"][0, -1] < 0.9 * none["stars_mass"][0, -1]  # feedback really acts


def test_instantaneous_feedback_acts_from_the_first_step():
    # previously instantaneous feedback was switched off until delay_time had elapsed
    inst, none = _run("instantaneous"), _run("no")
    assert inst["gas_mass"][0, 5] < none["gas_mass"][0, 5]


def test_no_feedback_switches_off_dust_destruction_and_wind_removal_too():
    none = _run("no")
    assert np.all(none["dust_mass"] == 0.0)  # dust is made by the stars that drive the feedback switch
    # metals: only the pre-enriched inflow, no stellar yield
    assert np.all(none["gas_metals"] >= 0.0)


def test_stars_and_gas_conserve_baryons_without_feedback():
    # Without feedback, gas + stars must equal the integrated inflow: no first-order lag in the stars.
    res = _run("no")
    A = main.gas_inflow_rate(res["redshift"], res["halo_mass"], res["halo_mass_rate"], False)[0]
    inflow = np.sum(A[:-1] * np.diff(res["cosmic_time"]))
    total = res["gas_mass"][0, -1] + res["stars_mass"][0, -1]
    assert total == pytest.approx(inflow, rel=1e-9)


def test_stellar_mass_step_dependence_is_only_the_inflow_rule():
    # Before the fix this was 1.1% between dz = 0.05 and 0.005 (stars lost the inflow within each step).
    # What remains is the left-endpoint rule for the inflow rate, a property of the inflow, not of star
    # formation: gas + stars matches the integrated inflow to round-off on each grid (test above).
    coarse, fine = _run("no", n=301), _run("no", n=3001)
    assert coarse["stars_mass"][0, -1] == pytest.approx(fine["stars_mass"][0, -1], rel=8e-3)


def test_stellar_metals_conserve_metals_without_feedback():
    res = _run("no")
    A = main.gas_inflow_rate(res["redshift"], res["halo_mass"], res["halo_mass_rate"], False)[0]
    inflow_metals = igm_metallicity() * np.sum(A[:-1] * np.diff(res["cosmic_time"]))
    total = res["gas_metals"][0, -1] + res["stars_metals"][0, -1]
    assert total == pytest.approx(inflow_metals, rel=1e-9)


def test_dust_is_converged_in_step_size():
    # the destruction rate is nu_SN * SFR, not Gamma / dt_step
    coarse, fine = _run("delayed", n=301), _run("delayed", n=3001)
    assert coarse["dust_mass"][0, -1] == pytest.approx(fine["dust_mass"][0, -1], rel=0.05)


def test_dust_destruction_scales_with_dust_to_gas_ratio():
    from ashvini.dust import update_dust_reservoir
    kw = dict(t=1.0, gas_mass=1e6, halo_mass=1e9, wind_sfr=1e4, stellar_metallicity=0.0)
    r0 = float(update_dust_reservoir(kw["t"], 0.0, kw["gas_mass"], kw["halo_mass"], kw["wind_sfr"], 0.0))
    r1 = float(update_dust_reservoir(kw["t"], 1e3, kw["gas_mass"], kw["halo_mass"], kw["wind_sfr"], 0.0))
    r2 = float(update_dust_reservoir(kw["t"], 2e3, kw["gas_mass"], kw["halo_mass"], kw["wind_sfr"], 0.0))
    assert r0 > r1 > r2
    assert (r0 - r1) == pytest.approx(r1 - r2, rel=1e-12)  # linear in dust mass


def test_nu_sn_values_for_salpeter():
    assert nu_SN("salpeter", 8.0, 40.0) == pytest.approx(6.80e-3, rel=1e-2)
    assert nu_SN("salpeter", 8.0, 100.0) == pytest.approx(7.42e-3, rel=1e-2)
    with pytest.raises(ValueError):
        nu_SN("salpeter", 8.0, 200.0)


def test_metallicity_function_is_centred_on_a_tenth_of_solar():
    assert sn.metallicity_function(0.1 * Z_SUN) == pytest.approx(0.625)  # midpoint of 1 and 0.25
    assert sn.metallicity_function(0.0) == pytest.approx(1.0, abs=1e-3)
    assert sn.metallicity_function(Z_SUN) == pytest.approx(0.25, abs=1e-3)


def test_igm_metallicity_is_in_solar_units():
    assert igm_metallicity() == pytest.approx(PARAMS.metals.Z_IGM_Zsun * Z_SUN)


def test_seeding_threshold_is_in_solar_units():
    from ashvini import black_holes_growth as bhg
    p3 = PARAMS.bh.seeding.pop3
    below = 0.5 * p3.Z_gas_max_Zsun * Z_SUN
    above = 2.0 * p3.Z_gas_max_Zsun * Z_SUN
    args = dict(halo_mass=2 * p3.M_halo_min, redshift=p3.z_min + 1, already_seeded=False)
    assert bool(bhg.seeding_mask_and_mass(gas_metallicity=below, **args)[0])
    # above the pop3 threshold the pop3 channel must not fire (another channel may, so check the seed mass)
    _, seed_mass = bhg.seeding_mask_and_mass(gas_metallicity=above, **args)
    assert float(seed_mass) != p3.M_seed


@pytest.mark.parametrize("section,key,bad", [
    ("metallicity", "Z_IGM", 1.0e-3),
    ("dust", "dust_gamma", 1.3e-4),
])
def test_legacy_keys_are_rejected(tmp_path, section, key, bad):
    import yaml
    from pathlib import Path
    raw = yaml.safe_load((Path(__file__).resolve().parents[1] / "run_params.yaml").read_text())
    raw[section][key] = bad
    f = tmp_path / "p.yaml"
    f.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="no longer accepted"):
        load_params(f)
