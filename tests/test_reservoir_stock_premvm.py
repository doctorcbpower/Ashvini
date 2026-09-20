import numpy as np
import pytest

from ashvini.constants import G, Msun_g, pc_cm
from ashvini.reservoir_stock_premvm import two_scale_transfer, nfw_dark_matter_mass, cohort_transfer, median_specific_angular_momentum, run_reservoir_stock
from ashvini.utils import time_at_z, z_at_time
from scipy.special import ndtr


def _grid(n=201, z0=25.0):
    t = np.linspace(time_at_z(np.array([z0]))[0], time_at_z(np.array([5.0]))[0], n)
    return z_at_time(t)


def _halo(n=201, m0=1e8, k=0.9):
    z = _grid(n)
    return (m0 * np.exp(k * (25.0 - z)))[None, :], z


def test_cohort_transfer_matches_analytic_release():
    # One cohort, a rising, then falling, then rising threshold. The mass released must be
    # M*(S(a)-S(b)) for each rise above the previous maximum, and nothing for a fall.
    sigma, ljm, M = 0.5, np.log(1.0e28), 1.0
    m = np.array([[M]]); lj = np.array([[ljm]]); surv = np.ones((1, 1)); cut = np.full((1, 1), -np.inf)
    S = lambda lnj: ndtr(-(lnj - ljm) / sigma)
    a, b, c_, d = ljm - 1.5, ljm - 0.8, ljm - 1.0, ljm - 0.1
    got = []
    for th in (a, b, c_, d):
        got.append(cohort_transfer(m, lj, surv, cut, np.array([th]), sigma, 1)[0])
    assert got[0] == pytest.approx(M * (1 - S(a)), rel=1e-12)
    assert got[1] == pytest.approx(M * (S(a) - S(b)), rel=1e-10)
    assert got[2] == 0.0
    assert got[3] == pytest.approx(M * (S(b) - S(d)), rel=1e-10)
    assert m[0, 0] == pytest.approx(M * S(d), rel=1e-10)


def test_first_step_transfer_is_lognormal_cdf_at_jthresh():
    # Heavy, frozen black hole (no accretion, no star formation): M_enc, hence j_thresh, is
    # constant, so the gas transferred from one cohort of fresh inflow is fresh * Phi(x).
    z = _grid(61)
    seed = 1e9
    rate = np.zeros((1, len(z))); rate[0, 0] = 2e6  # fresh gas only in step 1
    lam = 0.035
    jt = np.sqrt(G * seed * Msun_g * 100.0 * pc_cm)
    masses = np.logspace(9, 12, 400)
    jm = median_specific_angular_momentum(masses, 0.5 * (z[0] + z[1]), lam)
    hm = masses[np.argmin(abs(np.log(jm / jt) - np.log(1.5)))]
    m = np.full((1, len(z)), hm)
    out = run_reservoir_stock(m, z, seed, halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, epsilon_sf_ext=0.0,
                              wind_mode="none", gas_mass0=0.0, stars_mass0=0.0, lam_median=lam, sigma_lnj=0.5)
    fresh = out["mdot_in"][0, 1] * (out["cosmic_time"][1] - out["cosmic_time"][0])
    jmed = median_specific_angular_momentum(hm, 0.5 * (z[0] + z[1]), lam)
    expected = fresh * ndtr(np.log(jt / jmed) / 0.5)
    assert out["transfer"][0, 1] == pytest.approx(expected, rel=1e-6)
    assert 0.0 < expected < fresh


@pytest.mark.parametrize("mode", ["none", "extranuclear", "galaxy"])
def test_baryon_budget_is_conserved(mode):
    m, z = _halo()
    m[:, :30] = 0.0  # unformed prefix
    out = run_reservoir_stock(m, z, 1e5, wind_mode=mode)
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
    for k in ("gas_nuc", "gas_ext", "stars_nuc", "stars_ext", "bh_mass"):
        assert np.all(np.isfinite(out[k])) and np.all(out[k] >= 0.0)


def test_extreme_parameters_stay_finite_and_conserving():
    m, z = _halo()
    out = run_reservoir_stock(m, z, 1e5, wind_mode="galaxy", eta_sn_scale=1e6, epsilon_sf=5.0, eta_acc=1.0)
    assert np.all(np.isfinite(out["gas_mass"])) and np.all(out["gas_nuc"] >= 0.0) and np.all(out["gas_ext"] >= 0.0)
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
    assert out["limiter_engaged"].any()


def test_wind_none_equals_zero_loading():
    m, z = _halo()
    a = run_reservoir_stock(m, z, 1e5, wind_mode="none")
    b = run_reservoir_stock(m, z, 1e5, wind_mode="extranuclear", eta_sn_scale=0.0)
    for k in ("bh_mass", "gas_nuc", "gas_ext", "stars_nuc", "stars_ext"):
        assert np.array_equal(a[k], b[k])


def _first_step_modes():
    z = _grid(61)
    hm = np.full((1, len(z)), 1e10)
    rate = np.full((1, len(z)), 2e9)
    # tune lambda so that j_med equals j_thresh(1e9 Msun, 100 pc): about half of the fresh gas is eligible
    jt = np.sqrt(G * 1e9 * Msun_g * 100.0 * pc_cm)
    lam = 0.035 * jt / median_specific_angular_momentum(1e10, 0.5 * (z[0] + z[1]), 0.035)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, gas_mass0=0.0, stars_mass0=0.0, lam_median=lam,
              epsilon_sf_ext=0.05)
    return {mode: run_reservoir_stock(hm, z, 1e9, wind_mode=mode, **kw) for mode in ("none", "extranuclear", "galaxy")}


def test_wind_modes_act_on_the_right_gas():
    r = _first_step_modes()
    none, ext, gal = r["none"], r["extranuclear"], r["galaxy"]
    assert none["gas_nuc"][0, 1] > 0 and none["gas_ext"][0, 1] > 0  # both reservoirs populated
    # extranuclear wind: gas outside 100 pc drops, nuclear gas at step 1 is untouched
    assert ext["gas_ext"][0, 1] < none["gas_ext"][0, 1]
    assert ext["gas_nuc"][0, 1] == pytest.approx(none["gas_nuc"][0, 1], rel=1e-12)
    # galaxy wind: both drop by the same fraction (removal in proportion to gas mass)
    fn = 1.0 - gal["gas_nuc"][0, 1] / none["gas_nuc"][0, 1]
    fe = 1.0 - gal["gas_ext"][0, 1] / none["gas_ext"][0, 1]
    assert fn > 0 and fn == pytest.approx(fe, rel=1e-9)


def test_all_fresh_gas_reaches_nucleus_when_jmed_is_tiny_and_none_when_huge():
    m, z = _halo()
    lo = run_reservoir_stock(m, z, 1e5, lam_median=1e-8, wind_mode="none", epsilon_sf_ext=0.0)
    assert lo["gas_ext"].max() < 1e-6 * lo["gas_nuc"].max() + 1e-3
    hi = run_reservoir_stock(m, z, 1e5, lam_median=1e3, wind_mode="none")
    assert np.all(hi["gas_nuc"] == 0.0) and np.all(hi["bh_mass"] == 1e5)


def test_unlimited_supply_reproduces_eddington_growth():
    z = _grid(61)
    hm = np.full((1, len(z)), 1e10)
    rate = np.full((1, len(z)), 1e10)
    out = run_reservoir_stock(hm, z, 1e3, halo_growth_rate=rate, eta_acc=1e3, chi_crit=1e12, lam_median=1e-8,
                              wind_mode="none", epsilon_sf=0.0, epsilon_sf_ext=0.0, epsilon_f=0.0)  # no AGN wind: it would drain the stock
    k = 15
    expected = 1e3 * np.exp(out["kappa_edd"] * (out["cosmic_time"][k] - out["cosmic_time"][0]))
    assert out["bh_mass"][0, k] == pytest.approx(expected, rel=1e-6)


def test_jthresh_uses_nuclear_baryons_only():
    # Same halo, same fresh gas: adding a large extranuclear stellar mass must not change what is
    # transferred (only nuclear M_BH + M_gas,nuc + M_star,nuc set j_thresh).
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, epsilon_sf_ext=0.0, wind_mode="none", gas_mass0=0.0,
              n_rd=np.inf, sf_ext_timescale="hubble")  # no galaxy cut: outside stars enter no threshold
    a = run_reservoir_stock(hm, z, 1e6, stars_mass0=1e3, **kw)
    b = run_reservoir_stock(hm, z, 1e6, stars_mass0=1e10, **kw)
    assert np.array_equal(a["transfer"], b["transfer"])


@pytest.mark.parametrize("sink", ["galaxy", "nuclear"])
def test_agn_wind_sink_modes_conserve_and_stay_finite(sink):
    m, z = _halo()
    out = run_reservoir_stock(m, z, 1e5, agn_wind_sink=sink, wind_mode="galaxy")
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
    assert np.all(out["gas_nuc"] >= 0.0) and np.all(out["gas_ext"] >= 0.0)


def test_galaxy_agn_sink_lets_the_black_hole_grow_more():
    # King wind mass is ~1e4-1e5 x the mass the black hole gains: drawing it from the small
    # nuclear stock alone throttles growth; drawing it from the whole galaxy does not.
    m, z = _halo(m0=1e9, k=0.6)
    nuc = run_reservoir_stock(m, z, 1e4, agn_wind_sink="nuclear", wind_mode="none")
    gal = run_reservoir_stock(m, z, 1e4, agn_wind_sink="galaxy", wind_mode="none")
    assert gal["bh_mass"][0, -1] >= nuc["bh_mass"][0, -1]


def test_momentum_wind_loading_and_conservation():
    from ashvini.reservoir_stock_premvm import momentum_agn_wind_rate
    # loading = f_mom * eps/(1-eps) * c/sigma: 1 Msun of BH growth, sigma from the halo
    from ashvini import black_holes_growth as pz
    from ashvini.constants import c
    sig = float(pz.velocity_dispersion(1e10, 8.0))
    got = float(momentum_agn_wind_rate(1.0, 1e10, 8.0, 0.1, 1.0))
    assert got == pytest.approx(0.1 / 0.9 * c / sig, rel=1e-12)
    assert float(momentum_agn_wind_rate(1.0, 0.0, 8.0, 0.1)) == 0.0
    m, z = _halo()
    out = run_reservoir_stock(m, z, 1e5, agn_wind_type="momentum", wind_mode="galaxy")
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)


def test_momentum_wind_is_weaker_than_energy_wind_at_fiducial_coupling():
    from ashvini.reservoir_stock_premvm import momentum_agn_wind_rate
    from ashvini import black_holes_growth_slimdisk as bs
    e = float(bs.king_agn_wind_rate(1.0, halo_mass=1e10, redshift=8.0, epsilon_f=5e-4))
    p = float(momentum_agn_wind_rate(1.0, 1e10, 8.0, 0.1))
    assert p < e


def _one_cohort(mass=1.0, ljm=np.log(1.0e28)):
    return (np.array([[mass]]), np.zeros((1, 1)), np.array([[ljm]]), np.ones((1, 1)), np.ones((1, 1)),
            np.full((1, 1), -np.inf), np.full((1, 1), -np.inf))


def test_two_scale_transfer_partitions_a_cohort_exactly():
    # After any sequence of (galaxy, nuclear) thresholds the cohort splits as
    #   nucleus = M (1 - S(b)),  galaxy = M (S(b) - S(a)),  halo = M S(a),
    # with a, b the running maxima of the two thresholds.
    sigma, ljm = 0.5, np.log(1.0e28)
    h, g, lj, Sg, Sn, cg, cn = _one_cohort(1.0, ljm)
    S = lambda x: ndtr(-(x - ljm) / sigma)
    nuc = gal = 0.0
    a = b = -np.inf
    for tg, tn in ((ljm + 0.3, ljm - 1.5), (ljm + 1.0, ljm - 0.9), (ljm + 0.2, ljm - 1.2), (ljm + 1.5, ljm - 0.3)):
        Tg, Tn = two_scale_transfer(h, g, lj, Sg, Sn, cg, cn, np.array([tn]), np.array([tg]), sigma, 1)
        gal += Tg[0]; nuc += Tn[0]
        b = max(b, tn); a = max(a, tg, b)
        assert nuc == pytest.approx(1 - S(b), rel=1e-9)
        assert h[0, 0] == pytest.approx(S(a), rel=1e-9)
        assert g[0, 0] == pytest.approx(S(b) - S(a), rel=1e-6, abs=1e-12)


def test_nfw_dark_matter_mass_limits():
    assert float(nfw_dark_matter_mass(1e10, 1.0, 4.0, 0.156)) == pytest.approx(1e10 * (1 - 0.156), rel=1e-12)
    assert float(nfw_dark_matter_mass(1e10, 1e-6, 4.0, 0.156)) < 1e-6 * 1e10
    xs = np.array([0.01, 0.1, 0.5])
    assert np.all(np.diff(nfw_dark_matter_mass(1e10, xs, 4.0, 0.156)) > 0)


def test_two_scale_conserves_and_halo_gas_appears_for_a_small_galaxy_scale():
    m, z = _halo()
    kw = dict(wind_mode="galaxy", sf_ext_timescale="hubble")  # same clock at both scales: isolate the halo-gas effect
    small = run_reservoir_stock(m, z, 1e5, n_rd=0.3, **kw)
    big = run_reservoir_stock(m, z, 1e5, n_rd=1e3, **kw)
    assert np.all(np.abs(small["baryon_residual"]) < 1e-10) and np.all(np.abs(big["baryon_residual"]) < 1e-10)
    assert small["gas_halo"][0, -1] > 10 * big["gas_halo"][0, -1]
    assert small["stars_ext"][0, -1] < big["stars_ext"][0, -1]  # halo gas forms no stars


def test_legacy_single_scale_matches_two_scale_with_huge_galaxy_scale_and_hubble_clock():
    m, z = _halo()
    a = run_reservoir_stock(m, z, 1e5, n_rd=np.inf, sf_ext_timescale="hubble", wind_mode="none")
    b = run_reservoir_stock(m, z, 1e5, n_rd=1e6, sf_ext_timescale="hubble", wind_mode="none")
    for k in ("bh_mass", "gas_nuc", "stars_ext"):
        assert b[k][0, -1] == pytest.approx(a[k][0, -1], rel=1e-6, abs=1e-9)
    assert a["gas_halo"].max() == 0.0


def test_larger_nuclear_scale_lets_more_gas_reach_the_nucleus():
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, epsilon_sf_ext=0.0, wind_mode="none",
              gas_mass0=0.0, dm_nuclear=True)
    t = [run_reservoir_stock(hm, z, 1e4, R_nuc_pc=R, **kw)["transfer"].sum() for R in (30.0, 100.0, 300.0, 1000.0)]
    assert np.all(np.diff(t) > 0)


def test_dm_in_nuclear_threshold_increases_transfer():
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, epsilon_sf_ext=0.0, wind_mode="none", gas_mass0=0.0)
    a = run_reservoir_stock(hm, z, 1e4, dm_nuclear=False, **kw)["transfer"].sum()
    b = run_reservoir_stock(hm, z, 1e4, dm_nuclear=True, **kw)["transfer"].sum()
    assert b > a


def test_freefall_and_hubble_galaxy_clocks_differ_but_both_conserve():
    m, z = _halo()
    a = run_reservoir_stock(m, z, 1e5, sf_ext_timescale="freefall")
    b = run_reservoir_stock(m, z, 1e5, sf_ext_timescale="hubble")
    assert np.all(np.abs(a["baryon_residual"]) < 1e-10) and np.all(np.abs(b["baryon_residual"]) < 1e-10)
    assert a["stars_ext"][0, -1] != b["stars_ext"][0, -1]


def test_nuclear_scale_as_fraction_of_disc_scale_matches_the_equivalent_fixed_radius():
    from ashvini.paper_reservoir import virial_radius_velocity, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL
    z = _grid(61); hm = np.full((1, len(z)), 1e10); rate = np.full((1, len(z)), 2e9)
    kw = dict(halo_growth_rate=rate, eta_acc=0.0, epsilon_sf=0.0, epsilon_sf_ext=0.0, wind_mode="none",
              gas_mass0=0.0, dm_nuclear=True)
    xi, lam = 0.25, 0.035
    z_mid = 0.5 * (z[0] + z[1])
    R_vir, _ = virial_radius_velocity(1e10, z_mid, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL, 67.8)
    R_pc = xi * lam * R_vir / np.sqrt(2.0) / pc_cm
    a = run_reservoir_stock(hm, z, 1e4, R_nuc_rd=xi, **kw)["transfer"][0, 1]
    b = run_reservoir_stock(hm, z, 1e4, R_nuc_pc=R_pc, **kw)["transfer"][0, 1]
    assert a == pytest.approx(b, rel=1e-9) and a > 0
    # conservation with a moving nuclear scale
    m, zz = _halo()
    out = run_reservoir_stock(m, zz, 1e5, R_nuc_rd=xi, wind_mode="galaxy")
    assert np.all(np.abs(out["baryon_residual"]) < 1e-10)
