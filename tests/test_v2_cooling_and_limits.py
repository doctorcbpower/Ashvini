"""
Two groups of tests of ashvini.reservoir_stock_v2.

1. The cooling gate (new, default off): off means unchanged, a closed gate stops star formation and nuclear accretion without
   breaking baryon conservation, and the gate opens at the step where the virial temperature (or halo mass) first reaches the
   threshold, computed here independently.
2. Analytic limits: closed-form results that need no knowledge of the code, evaluated with astropy's constants and cosmology
   so that a slip in a unit or a constant in the model shows up as a failure. These are checks of the equations in the paper
   against the implementation, not regression tests.
"""
import numpy as np
import pytest
from astropy import constants as const, units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.special import ndtr

from ashvini import reservoir_stock_v2 as v2
from ashvini import reservoir_stock as frozen
from ashvini import black_holes_growth_slimdisk as bh_slimdisk
from ashvini import supernovae_feedback as sn
from ashvini.paper_reservoir import F_B_FIDUCIAL, LAMBDA_MEDIAN_FIDUCIAL, virial_radius_velocity
from ashvini.reservoir_stock import nfw_dark_matter_mass
from ashvini.paper_reservoir_params import PAPER_PARAMS

from test_reservoir_stock_v2 import smooth_history, jumpy_history

COSMO = FlatLambdaCDM(H0=67.8, Om0=0.308)      # the paper's halo-structure cosmology
G_CGS = const.G.cgs.value
MSUN_G = const.M_sun.cgs.value
PC_CM = const.pc.cgs.value


def t_vir_independent(M_msun, z, mu=1.22, delta=200.0):
    """T_vir = mu m_p V^2 / (2 k_B), V^2 = G M / R, M = (4 pi / 3) delta rho_crit(z) R^3, all from astropy."""
    rho = COSMO.critical_density(z).cgs.value
    M = M_msun * MSUN_G
    R = (3.0 * M / (4.0 * np.pi * delta * rho)) ** (1.0 / 3.0)
    V2 = G_CGS * M / R
    return mu * const.m_p.cgs.value * V2 / (2.0 * const.k_B.cgs.value)


# ---------------------------------------------------------------------------------------------------------------------
# 1. cooling gate
# ---------------------------------------------------------------------------------------------------------------------
KEYS = ("bh_mass", "stars_nuc", "stars_gal", "gas_nuc", "gas_gal", "gas_unavailable", "transfer_nuc")


def test_cooling_gate_always_open_or_absent_changes_nothing():
    z, t, M = smooth_history(301)
    a = v2.run_reservoir_v2(M, z, 1e3)
    b = v2.run_reservoir_v2(M, z, 1e3, cooling_T_vir_K=1e-3)       # threshold far below any halo: always open
    c = v2.run_reservoir_v2(M, z, 1e3, cooling_mass=1.0)
    for key in KEYS:
        assert np.array_equal(a[key], b[key]), key
        assert np.array_equal(a[key], c[key]), key


def test_closed_cooling_gate_stops_star_formation_and_accretion_and_keeps_baryons():
    z, t, M = smooth_history(301)
    r = v2.run_reservoir_v2(M, z, 1e3, cooling_T_vir_K=1e12)
    assert np.all(r["stars_mass"] == 0.0)
    assert np.all(r["bh_mass"] == 1e3)
    assert r["gas_mass"][0, -1] > 0.0                                # gas still arrives and accumulates
    assert np.max(np.abs(r["baryon_residual"])) < 1e-12


@pytest.mark.parametrize("kind", ["T_vir", "mass"])
def test_cooling_gate_opens_at_the_independently_computed_crossing(kind):
    z, t, M = smooth_history(801)
    T_thr = 1.0e4
    z_mid = 0.5 * (z[1:] + z[:-1])
    Tv = t_vir_independent(np.where(M[0, :-1] > 0, M[0, :-1], 1.0), z_mid)
    first_open = int(np.argmax((Tv >= T_thr) & (M[0, :-1] > 0)))      # step index k = j - 1
    assert 0 < first_open < len(z_mid) - 1
    if kind == "T_vir":
        r = v2.run_reservoir_v2(M, z, 1e3, cooling_T_vir_K=T_thr)
    else:
        r = v2.run_reservoir_v2(M, z, 1e3, cooling_mass=float(M[0, first_open]))
    sfr = r["sfr_gal"][0]
    j_first = first_open + 1                                          # output index j = k + 1
    assert np.all(sfr[:j_first] == 0.0)
    assert sfr[j_first] > 0.0


def test_cooling_gate_delays_star_formation_and_still_finds_one_critical_seed():
    z, t, M = jumpy_history(301, N=4)
    gated = v2.critical_seed_v2(M, z, n_iter=30, cooling_T_vir_K=1.0e4)[0]
    assert np.all(np.isfinite(gated))
    # Until the gate opens a gated halo has formed no stars, so for a while its stellar mass is zero where the ungated halo's is not.
    # (At z = 5 the gated halo need not have less: gas that has accumulated is not blown out early by feedback, so it can end
    # with more stars.)
    a = v2.run_reservoir_v2(M, z, 1e6)
    b = v2.run_reservoir_v2(M, z, 1e6, cooling_T_vir_K=1.0e4)
    delayed = (b["stars_mass"] == 0.0) & (a["stars_mass"] > 0.0)
    assert delayed.any()
    assert np.all(b["stars_mass"][~np.isfinite(b["stars_mass"])] == 0)


# ---------------------------------------------------------------------------------------------------------------------
# 2. analytic limits
# ---------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("M_enc,R_pc", [(1e6, 100.0), (3e8, 100.0), (1e10, 1000.0)])
def test_free_fall_time_is_the_paper_equation(M_enc, R_pc):
    """t_ff = [r^3 / (2 G M_enc)]^(1/2), Eq. (tff) of the paper."""
    expect_s = np.sqrt((R_pc * PC_CM) ** 3 / (2.0 * G_CGS * M_enc * MSUN_G))
    got_gyr = float(bh_slimdisk.time_freefall_enclosed(M_enc, R_nuc_pc=R_pc))
    assert got_gyr == pytest.approx(expect_s / (1e9 * 3.15576e7), rel=2e-3)    # Gyr_s differs from the Julian year at the 1e-3 level


def test_virial_radius_velocity_match_critical_density_definition():
    for M, zz in [(1e8, 20.0), (3e10, 10.0), (3e13, 5.0)]:
        R, V = virial_radius_velocity(M, zz, 0.308, 0.692, 67.8)
        rho = COSMO.critical_density(zz).cgs.value
        R_ref = (3.0 * M * MSUN_G / (4.0 * np.pi * 200.0 * rho)) ** (1.0 / 3.0)
        assert R == pytest.approx(R_ref, rel=2e-3)
        assert V == pytest.approx(np.sqrt(G_CGS * M * MSUN_G / R_ref), rel=1e-3)


def test_nfw_enclosed_mass_limits():
    c = PAPER_PARAMS.c_nfw
    M = 3e11
    assert nfw_dark_matter_mass(M, 1.0, c, F_B_FIDUCIAL) == pytest.approx((1.0 - F_B_FIDUCIAL) * M, rel=1e-12)
    # inner profile: m(y) = y^2 / 2 - 2 y^3 / 3 + ..., so the enclosed mass scales as x^2 at small radius
    x = 1e-4
    y = c * x
    m_c = np.log1p(c) - c / (1.0 + c)
    expect = (1.0 - F_B_FIDUCIAL) * M * (0.5 * y ** 2 - 2.0 * y ** 3 / 3.0) / m_c
    assert nfw_dark_matter_mass(M, x, c, F_B_FIDUCIAL) == pytest.approx(expect, rel=1e-5)


def test_eddington_efolding_time_is_the_salpeter_time():
    """For ample supply the black hole grows as exp(t / t_Sal), t_Sal = [eps / (1 - eps)] sigma_T c / (4 pi G m_p) ~ 50 Myr."""
    eps = 0.1
    t_edd = (const.sigma_T * const.c / (4.0 * np.pi * const.G * const.m_p)).to(u.Gyr).value
    t_sal = eps / (1.0 - eps) * t_edd
    assert t_sal == pytest.approx(0.0500, abs=5e-4)                   # the paper's "about 50 Myr"
    kappa = float(bh_slimdisk.eddington_rate_std_per_unit_mass(eps))
    # the code uses rounded cgs constants (m_p = 1.67e-24, sigma_T = 6.65e-25, c = 3e10), good to about 0.5 per cent in t_Sal
    assert 1.0 / kappa == pytest.approx(t_sal, rel=5e-3)
    M1 = frozen.strict_eddington_step(1e5, 1e30, kappa, 0.1)         # supply never binds: exactly exp(kappa dt)
    assert M1 == pytest.approx(1e5 * np.exp(kappa * 0.1), rel=1e-12)
    assert M1 == pytest.approx(1e5 * np.exp(0.1 / t_sal), rel=1e-2)


def test_supply_limited_growth_is_linear_in_the_supply():
    kappa = float(bh_slimdisk.eddington_rate_std_per_unit_mass(0.1))
    M0, A, dt = 1e8, 1e-3, 0.01                                       # supply far below kappa * M0
    assert frozen.strict_eddington_step(M0, A, kappa, dt) == pytest.approx(M0 + A * dt, rel=1e-15)


def test_star_formation_demand_is_epsilon_over_tff():
    assert frozen.star_formation_demand(1e9, 0.05, 0.015, 0.001) == pytest.approx(0.015 * 1e9 / 0.05 * 0.001, rel=1e-15)


def test_supernova_mass_loading_is_the_paper_scaling():
    """eta = 5 (10^11.5 / M)^(1/3) (9 / (1 + z))^(1/2) for the fiducial configuration (paper Section 2)."""
    for M, zz in [(3e10, 8.0), (3e11, 12.0), (3e12, 5.0)]:
        expect = 5.0 * (10 ** 11.5 / M) ** (1.0 / 3.0) * (9.0 / (1.0 + zz)) ** 0.5
        assert sn.mass_loading_factor(zz, M, 0.0) == pytest.approx(expect, rel=1e-3)


def test_agn_wind_mass_loading_is_momentum_conserving():
    """King (2003): Mdot_out v = f L / c with L = eps / (1 - eps) Mdot_BH c^2 and v = sigma, so Mdot_out / Mdot_BH = f eps c / ((1 - eps) sigma)."""
    M, zz, eps, f = 3e11, 8.0, 0.1, 1.0
    sigma = float(bh_slimdisk.pznk11.velocity_dispersion(M, zz))
    expect = f * eps / (1.0 - eps) * const.c.cgs.value / sigma
    # the code takes c = 3e10 cm/s, 0.07 per cent above the exact value
    assert frozen.momentum_agn_wind_rate(1.0, M, zz, eps, f) == pytest.approx(expect, rel=2e-3)


def test_delivered_fraction_in_one_step_is_the_lognormal_cdf():
    """From empty reservoirs, the gas delivered to the nucleus in one step is fresh * Phi((ln j_thr - ln j_med) / sigma_j),
    with j_thr = sqrt(G M_enc R_nuc) built from the dark matter inside R_nuc plus the seed, and j_med = lambda sqrt(2) V_vir R_vir."""
    M1, M2, seed, sigma_j, R_pc = 1.0e9, 1.3e9, 5.0e5, 0.5, 100.0
    cosmo = v2.PaperCosmology()
    z = np.array([12.0, 11.8])
    M = np.array([[M1, M2]])
    r = v2.run_reservoir_v2(M, z, seed, sigma_lnj=sigma_j, R_nuc_pc=R_pc, cosmology=cosmo, dm_nuclear=True)
    fresh = r["mdot_in"][0, 1] * (r["cosmic_time"][1] - r["cosmic_time"][0])
    T = r["transfer_nuc"][0, 1]
    z_mid = 0.5 * (z[0] + z[1])
    rho = COSMO.critical_density(z_mid).cgs.value
    R_vir = (3.0 * M1 * MSUN_G / (4.0 * np.pi * 200.0 * rho)) ** (1.0 / 3.0)
    V = np.sqrt(G_CGS * M1 * MSUN_G / R_vir)
    j_med = LAMBDA_MEDIAN_FIDUCIAL * np.sqrt(2.0) * V * R_vir
    c = PAPER_PARAMS.c_nfw
    m = lambda y: np.log1p(y) - y / (1.0 + y)
    M_dm = (1.0 - F_B_FIDUCIAL) * M1 * m(c * R_pc * PC_CM / R_vir) / m(c)
    j_thr = np.sqrt(G_CGS * (M_dm + seed) * MSUN_G * R_pc * PC_CM)
    expect = fresh * ndtr((np.log(j_thr) - np.log(j_med)) / sigma_j)
    assert fresh > 0 and T > 0
    assert T == pytest.approx(expect, rel=5e-3)
