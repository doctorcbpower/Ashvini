import numpy as np
from ashvini.paper_reservoir import run_reservoir_paper
from ashvini.utils import time_at_z, z_at_time


def _halo():
    t = np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], 201)
    z = z_at_time(t)
    m = 1e8 * np.exp(0.9 * (25.0 - z))[None, :]
    return m, z


def test_stellar_winds_off_matches_zero_scale():
    m, z = _halo()
    a = run_reservoir_paper(m, z, 1e5, stellar_winds=False)
    b = run_reservoir_paper(m, z, 1e5, stellar_winds=True, eta_sn_scale=0.0)
    for k in ("bh_mass", "gas_mass", "stars_mass"):
        assert np.array_equal(a[k], b[k])
    assert np.all(b["mdot_out_sn"] == 0.0)


def test_stellar_winds_remove_gas_and_reduce_stars():
    m, z = _halo()
    off = run_reservoir_paper(m, z, 1e5, stellar_winds=False)
    on = run_reservoir_paper(m, z, 1e5, stellar_winds=True)
    assert on["stars_mass"][0, -1] < off["stars_mass"][0, -1]
    assert np.all(on["mdot_out_sn"] >= 0.0) and on["mdot_out_sn"].max() > 0.0
    assert np.all(np.isfinite(on["gas_mass"]))


def test_stellar_winds_unformed_halo_is_finite():
    m, z = _halo()
    m[:, :40] = 0.0
    out = run_reservoir_paper(m, z, 1e5, stellar_winds=True)
    assert np.all(np.isfinite(out["gas_mass"])) and np.all(out["mdot_out_sn"][:, :41] == 0.0)


def test_wind_never_drains_more_than_extranuclear_gas():
    m, z = _halo()
    out = run_reservoir_paper(m, z, 1e5, stellar_winds=True, eta_sn_scale=1e6)
    assert np.all(np.isfinite(out["gas_mass"])) and np.all(out["gas_mass"] >= 0.0)
    dt = np.diff(out["cosmic_time"])
    ejected = out["mdot_out_sn"][:, 1:] * dt
    assert np.all(ejected <= out["gas_mass"][:, :-1] + 1e-6)
