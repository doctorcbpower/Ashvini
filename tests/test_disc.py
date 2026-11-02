"""The disc option: a disc size for the gas and a star formation threshold that depends on it (components/disc.py)."""
import numpy as np
import pytest

from ashvini import halo, utils
from ashvini import reservoir_stock as reservoir
from ashvini.components import disc
from ashvini.components.star_formation import FreeFallStarFormation
from ashvini.config import DEFAULT, override
from ashvini.integrators.grid import GridIntegrator
from ashvini.main import run_events, run_forest
from ashvini.paper_reservoir import OMEGA_LAMBDA_FIDUCIAL, OMEGA_M_FIDUCIAL

NO_SEEDS = {"bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
            "bh.seeding.halo_mass_threshold.enabled": False, "reion.UVB_enabled": False}


def history(M_start=2e7, M_end=2e9, z0=20.0, z1=8.0, n=120):
    z = np.linspace(z0, z1, n)
    M = np.exp(np.log(M_start) + (z0 - z) / (z0 - z1) * (np.log(M_end) - np.log(M_start)))[None, :]
    return M, np.maximum(np.gradient(M, utils.time_at_z(z), axis=1), 0.0), z


def test_disc_is_off_by_default():
    assert DEFAULT.disc.enabled is False


def test_scale_length_is_spin_times_virial_radius_over_root_two_as_in_the_reservoir_model():
    cfg = override(DEFAULT, **{"disc.spin_parameter": 0.04})
    for M, z in ((1e8, 5.0), (1e10, 5.0), (1e10, 10.0)):
        R_vir, _ = reservoir.virial_radius_velocity(M, z, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL, 67.8)
        expected_kpc = 0.04 * float(R_vir) / np.sqrt(2.0) / 3.0857e21
        assert float(disc.scale_length_kpc(M, z, cfg)) == pytest.approx(expected_kpc, rel=5e-3)   # the two cosmologies differ by 0.3 %
    doubled = override(cfg, **{"disc.angular_momentum_ratio": 2.0})
    assert float(disc.scale_length_kpc(1e9, 6.0, doubled)) == pytest.approx(2.0 * float(disc.scale_length_kpc(1e9, 6.0, cfg)), rel=1e-12)
    assert float(disc.scale_length_kpc(0.0, 6.0, cfg)) == 0.0


def test_star_forming_fraction_is_the_gas_inside_the_critical_radius():
    sigma0, crit = 40.0, 5.0
    R_d = 1.0
    r = np.linspace(0.0, 30.0, 600001)
    sigma = sigma0 * np.exp(-r / R_d)
    mass_density = 2.0 * np.pi * r * sigma
    total = np.trapezoid(mass_density, r)
    inside = np.trapezoid(np.where(sigma > crit, mass_density, 0.0), r)
    assert float(disc.star_forming_fraction(sigma0, crit)) == pytest.approx(inside / total, rel=1e-4)
    assert float(disc.star_forming_fraction(3.0, 5.0)) == 0.0          # the centre is already below the threshold
    assert float(disc.star_forming_fraction(0.0, 5.0)) == 0.0
    assert float(disc.star_forming_fraction(3.0, 0.0)) == 1.0          # no threshold
    f = disc.star_forming_fraction(np.logspace(0, 6, 50), 5.0)
    assert np.all(np.diff(f) >= 0) and f[-1] > 0.999                   # more gas, or a smaller disc: a larger fraction


def test_central_surface_density_of_an_exponential_disc():
    # 2 pi R_d^2 Sigma_0 is the total mass: 1e8 Msun in R_d = 1 kpc gives 1e8 / (2 pi 1e6 pc^2)
    assert float(disc.central_surface_density(1e8, 1.0)) == pytest.approx(1e8 / (2 * np.pi * 1e6), rel=1e-12)
    assert float(disc.central_surface_density(1e8, 0.0)) == 0.0


@pytest.mark.parametrize("stepper", ["grid", "events"])
def test_no_threshold_gives_the_model_without_a_disc(stepper):
    M, Mdot, z = history()
    base = override(DEFAULT, **NO_SEEDS)
    with_disc = override(base, **{"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 0.0})
    if stepper == "grid":
        a, b = run_forest(M, Mdot, z, config=base), run_forest(M, Mdot, z, config=with_disc)
    else:
        kw = {"sn.m_sn_max": 100.0, "sn.delay_distribution": "stellar_lifetimes"}
        a = run_events(M[:, ::4], Mdot[:, ::4], z[::4], config=override(base, **kw))
        b = run_events(M[:, ::4], Mdot[:, ::4], z[::4], config=override(with_disc, **kw))
    for key in ("gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass"):
        assert np.allclose(a[key], b[key], rtol=1e-12, atol=0.0), key


def test_a_very_large_spin_spreads_the_gas_too_thin_to_form_stars():
    M, Mdot, z = history()
    spread = override(DEFAULT, **{**NO_SEEDS, "disc.enabled": True, "disc.spin_parameter": 50.0})
    r = run_forest(M, Mdot, z, config=spread)
    assert r["stars_mass"][0, -1] == 0.0
    assert r["gas_mass"][0, -1] > 0.0
    e = run_events(M[:, ::4], Mdot[:, ::4], z[::4], config=override(spread, **{"sn.m_sn_max": 100.0, "sn.delay_distribution": "stellar_lifetimes"}))
    assert e["stars_mass"][0, -1] == 0.0


def test_without_feedback_the_threshold_can_only_reduce_the_stars_and_the_gas_waits():
    M, Mdot, z = history()
    base = override(DEFAULT, **{**NO_SEEDS, "sn.type": "no"})
    with_disc = override(base, **{"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 20.0})
    a, b = run_forest(M, Mdot, z, config=base), run_forest(M, Mdot, z, config=with_disc)
    assert b["stars_mass"][0, -1] < a["stars_mass"][0, -1]
    assert b["gas_mass"][0, -1] > a["gas_mass"][0, -1]
    supplied = a["gas_mass"][0, -1] + a["stars_mass"][0, -1]
    assert b["gas_mass"][0, -1] + b["stars_mass"][0, -1] == pytest.approx(supplied, rel=1e-6)     # nothing is lost


def test_the_star_forming_fraction_changes_smoothly_in_time():
    # a high threshold, so that the fraction passes through the whole range as the halo grows
    M, Mdot, z = history(n=300)
    cfg = override(DEFAULT, **{**NO_SEEDS, "disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 400.0})
    r = run_forest(M, Mdot, z, config=cfg)
    f = disc.diagnostics(r["halo_mass"], r["gas_mass"], z, cfg)["star_forming_fraction"][0]
    assert f[5:].min() < 0.3 and f.max() > 0.6                  # the threshold bites, and is then passed
    assert np.max(np.abs(np.diff(f[3:]))) < 0.1                 # no switching on and off from step to step (the first steps fill an empty halo)


def test_events_conserve_baryons_with_a_disc_threshold_and_form_stars_only_once_the_gas_is_dense_enough():
    M, Mdot, z = history()
    cfg = override(DEFAULT, **{**NO_SEEDS, "disc.enabled": True, "sn.m_sn_max": 100.0,
                               "sn.delay_distribution": "stellar_lifetimes", "sn.wind_model": "momentum", "sn.f_p": 0.76})
    Ms = M[:, ::2]
    r = run_events(Ms, Mdot[:, ::2], z[::2], config=cfg, gas_mass0=1e3)
    supplied = 1e3 + (utils.Omega_b / utils.Omega_m) * (Ms[0, -1] - Ms[0, 0])
    assert r["gas_mass"][0, -1] + r["stars_mass"][0, -1] + r["ejected_mass"][0, -1] == pytest.approx(supplied, rel=1e-9)
    first = int(np.argmax(r["stars_mass"][0] > 0))
    assert first > 0 and r["gas_mass"][0, first - 1] > 0              # the gas waited
    sigma0 = disc.diagnostics(r["halo_mass"], r["gas_mass"], z[::2], cfg)["central_surface_density"][0]
    assert sigma0[first - 1] < cfg.disc.critical_surface_density_Msun_pc2 or first == 1


def test_a_star_formation_law_written_for_rate_coefficient_still_works_in_the_grid_integrator():
    class Faster(FreeFallStarFormation):
        def rate_coefficient(self, redshift):
            return 2.0 * super().rate_coefficient(redshift)
    M, Mdot, z = history()
    integrator = GridIntegrator(override(DEFAULT, **NO_SEEDS))
    plain, _ = integrator.run(M, Mdot, z)
    integrator.star_formation = Faster(integrator.cfg)
    faster, _ = integrator.run(M, Mdot, z)
    assert faster.stars_mass[0, -1] > plain.stars_mass[0, -1]
