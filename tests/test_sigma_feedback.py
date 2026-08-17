"""
Tests for the optional isothermal-sphere M-sigma self-regulation AGN
coupling (black_holes.sigma_feedback -- King 2003, 2005; Power, Zubovas,
Nayakshin & King 2011), replacing the constant eta_agn AGN wind coupling
with one gated by M_BH relative to M_sigma(halo_mass, redshift). See
ashvini.black_holes_growth.velocity_dispersion/m_sigma,
ashvini.agn_feedback.coupling_switch/agn_wind_mass_rate, and MODELS.md's
"AGN feedback" section.
"""

import os

import numpy as np
import pytest

from ashvini import agn_feedback as agn
from ashvini import black_holes_growth as bh_growth
from ashvini import main
from ashvini import utils

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "merger_trees_fixture.h5")


# ---------------------------------------------------------------------------
# (a) sigma(M_halo, z) and M_sigma(sigma) limits
# ---------------------------------------------------------------------------

def test_velocity_dispersion_zero_halo_mass_gives_zero():
    assert bh_growth.velocity_dispersion(0.0, 6.0) == 0.0
    result = bh_growth.velocity_dispersion(np.array([0.0, 1e12]), np.array([6.0, 6.0]))
    assert result[0] == 0.0
    assert result[1] > 0.0
    assert np.all(np.isfinite(result))


def test_velocity_dispersion_increases_with_halo_mass():
    sigma_lo = bh_growth.velocity_dispersion(1e10, 0.0)
    sigma_hi = bh_growth.velocity_dispersion(1e13, 0.0)
    assert 0.0 < sigma_lo < sigma_hi


def test_m_sigma_zero_sigma_gives_zero():
    assert bh_growth.m_sigma(0.0) == 0.0


def test_m_sigma_increases_with_sigma():
    m_lo = bh_growth.m_sigma(1e7)  # cm/s
    m_hi = bh_growth.m_sigma(2e7)
    assert 0.0 < m_lo < m_hi
    # M_sigma ~ sigma^4: doubling sigma should scale M_sigma by 2^4 = 16
    assert m_hi / m_lo == pytest.approx(16.0, rel=1e-6)


def test_m_sigma_physical_normalisation_sanity():
    # sigma = 200 km/s is the classic M-sigma-relation benchmark scale;
    # King (2003)'s f_g=0.16 normalisation gives M_sigma of a few 1e8 Msun
    # there -- a loose sanity check, not a precision calibration.
    sigma_200kms_cgs = 200e5  # cm/s
    m_sig = bh_growth.m_sigma(sigma_200kms_cgs)
    assert 1e8 < m_sig < 1e9


# ---------------------------------------------------------------------------
# (b) coupling_switch limits
# ---------------------------------------------------------------------------

def test_coupling_switch_large_ratio_approaches_one():
    switch = agn.coupling_switch(1e12, 1e8)  # M_BH/M_sigma = 1e4
    assert switch == pytest.approx(1.0, abs=1e-6)


def test_coupling_switch_small_ratio_approaches_zero():
    switch = agn.coupling_switch(1e4, 1e8)  # M_BH/M_sigma = 1e-4
    assert switch == pytest.approx(0.0, abs=1e-6)


def test_coupling_switch_centered_at_equality():
    assert agn.coupling_switch(1e8, 1e8) == pytest.approx(0.5)


def test_coupling_switch_zero_bh_or_sigma_mass_gives_zero():
    assert agn.coupling_switch(0.0, 1e8) == 0.0
    assert agn.coupling_switch(1e8, 0.0) == 0.0
    assert agn.coupling_switch(0.0, 0.0) == 0.0


def test_agn_wind_mass_rate_eff_coupling_approaches_eta_agn(monkeypatch):
    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)
    bh_growth_rate = 1e5

    # M_BH >> M_sigma: coupling switch -> 1, so eta_agn_eff -> eta_agn
    rate_saturated = agn.agn_wind_mass_rate(
        bh_growth_rate, bh_mass=1e20, halo_mass=1e12, redshift=6.0
    )
    m_sig = bh_growth.m_sigma(bh_growth.velocity_dispersion(1e12, 6.0))
    assert 1e20 / m_sig > 1e6  # sanity: genuinely deep in the saturated regime
    assert rate_saturated == pytest.approx(agn.eta_agn * bh_growth_rate, rel=1e-6)


def test_agn_wind_mass_rate_eff_coupling_approaches_zero_below_m_sigma(monkeypatch):
    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)
    bh_growth_rate = 1e5

    # M_BH << M_sigma: coupling switch -> 0
    rate_trapped = agn.agn_wind_mass_rate(
        bh_growth_rate, bh_mass=1.0, halo_mass=1e12, redshift=6.0
    )
    assert rate_trapped == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# (c) sigma_feedback.enabled=False falls back exactly to the prior behaviour
# ---------------------------------------------------------------------------

def test_disabled_sigma_feedback_matches_constant_eta_agn_function(monkeypatch):
    monkeypatch.setattr(agn, "sigma_feedback_enabled", False)
    bh_growth_rate = 1e5
    # Passing bh_mass/halo_mass/redshift must be harmless (ignored) when disabled.
    rate = agn.agn_wind_mass_rate(bh_growth_rate, bh_mass=1.0, halo_mass=1e6, redshift=20.0)
    assert rate == pytest.approx(agn.eta_agn * bh_growth_rate)


def test_disabled_sigma_feedback_reproduces_pinned_run1_baseline():
    # sigma_feedback is off by default (PARAMS.bh.sigma_feedback.enabled is
    # False in run_params.yaml) -- run1() on the pinned test fixture must
    # therefore reproduce test_run1.py's REFERENCE_FINAL_VALUES exactly,
    # confirming this feature introduces no behaviour change when disabled.
    from tests.test_run1 import REFERENCE_FINAL_VALUES

    halo_masses, halo_mass_rates, redshift = utils.read_trees(
        file_path=FIXTURE_PATH, mass_bin=1e10
    )
    result = main.run1(halo_masses[0], halo_mass_rates[0], redshift)

    for key, expected in REFERENCE_FINAL_VALUES.items():
        assert result[key][-1] == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Cross-validation: run_forest() vs run1_scalar() with sigma_feedback enabled
# ---------------------------------------------------------------------------

def test_m_sigma_matches_observed_relation_order_of_magnitude():
    # Sanity check against actual observational M-sigma calibrations, not
    # just King's own normalisation (test_m_sigma_physical_normalisation_
    # sanity above). McConnell & Ma (2013, ApJ 764, 184), all-galaxies fit:
    # M_BH = 10^8.32 * (sigma/200 km/s)^5.64 Msun. Our model instead
    # implements King (2003)'s theoretical sigma^4 scaling -- shallower
    # than the observed ~5.64 slope -- so exact agreement away from the
    # sigma=200 km/s pivot point isn't expected; this is a loose
    # order-of-magnitude check across a modest range (100-300 km/s), not a
    # precision calibration.
    for sigma_kms in (100.0, 200.0, 300.0):
        sigma_cgs = sigma_kms * 1e5  # km/s -> cm/s
        m_sig = bh_growth.m_sigma(sigma_cgs)
        m_observed = 10 ** (8.32 + 5.64 * np.log10(sigma_kms / 200.0))
        ratio = m_sig / m_observed
        assert 0.1 < ratio < 10.0, (
            f"M_sigma({sigma_kms} km/s)={m_sig:.3g} Msun is more than an "
            f"order of magnitude off the McConnell & Ma (2013) relation "
            f"({m_observed:.3g} Msun, ratio={ratio:.3g})"
        )


def test_sigma_feedback_mostly_dormant_at_group_scale_under_realistic_defaults(monkeypatch):
    # Does the mechanism ever actually engage under run_params.yaml's real
    # defaults (efficiency=0.001, eddington_multiplier=1.0), or is it
    # currently dormant for any realistic run? Checked directly against the
    # pinned test fixture (tests/fixtures/merger_trees_fixture.h5): at
    # group/cluster halo scale (1e10, 1e11 Msun mass bins), M_BH stays
    # >=~5x below M_sigma throughout every halo's history under real
    # parameters -- the mechanism is a no-op there in practice, not because
    # it's disabled but because BH growth genuinely never gets close.
    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)
    # bh_growth.e_bh / eddington_multiplier intentionally left at their
    # real PARAMS defaults (unlike the other tests here) -- this is
    # specifically checking realistic, not stress-tested, behaviour.

    fixture_path = os.path.join(os.path.dirname(__file__), "fixtures", "merger_trees_fixture.h5")
    for mass_bin in (1e10, 1e11):
        halo_masses, halo_mass_rates, redshift = utils.read_trees(fixture_path, mass_bin)
        result = main.run_forest(halo_masses, halo_mass_rates, redshift)

        sigma = bh_growth.velocity_dispersion(halo_masses, np.broadcast_to(redshift, halo_masses.shape))
        m_sigma_hist = bh_growth.m_sigma(sigma)
        ratio = np.divide(
            result["bh_mass"], m_sigma_hist, out=np.zeros_like(result["bh_mass"]),
            where=m_sigma_hist > 0,
        )
        assert np.all(np.isfinite(ratio))
        assert ratio.max() < 0.5, (
            f"mass_bin={mass_bin:.0e}: M_BH/M_sigma reached {ratio.max():.3g} under "
            f"realistic defaults -- expected the mechanism to stay dormant "
            f"(well below the M_BH=M_sigma crossing) at this halo mass scale"
        )


def test_sigma_feedback_can_engage_at_low_mass_scale_under_realistic_defaults(monkeypatch):
    # Contrast with the group-scale check above: at the low-mass end,
    # M_sigma itself is tiny (comparable to the BH seed masses themselves,
    # ~1e2-1e5 Msun -- see black_holes.seeding in run_params.yaml), so the
    # switch *can* engage under real defaults immediately upon seeding, not
    # just under the stress-test parameters used elsewhere in this file.
    # Checked directly against the fixture's 1e8 Msun mass bin, which does
    # cross M_sigma for at least one halo under real parameters -- this
    # pins that real (if narrow) behaviour rather than assuming dormancy
    # everywhere.
    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)

    fixture_path = os.path.join(os.path.dirname(__file__), "fixtures", "merger_trees_fixture.h5")
    halo_masses, halo_mass_rates, redshift = utils.read_trees(fixture_path, 1e8)
    result = main.run_forest(halo_masses, halo_mass_rates, redshift)

    sigma = bh_growth.velocity_dispersion(halo_masses, np.broadcast_to(redshift, halo_masses.shape))
    m_sigma_hist = bh_growth.m_sigma(sigma)
    ratio = np.divide(
        result["bh_mass"], m_sigma_hist, out=np.zeros_like(result["bh_mass"]),
        where=m_sigma_hist > 0,
    )
    assert np.all(np.isfinite(ratio))
    assert ratio.max() > 1.0, (
        "expected at least one halo in the 1e8 Msun fixture to cross "
        "M_BH=M_sigma under realistic defaults (M_sigma is comparable to "
        "the seed mass scale here) -- if this now fails, the low-mass "
        "M_sigma normalisation or seeding parameters have likely changed"
    )


def test_sigma_feedback_wind_alone_does_not_cap_bh_growth_near_m_sigma(monkeypatch, tmp_path):
    # The core physical question for a mechanism called "self-regulation":
    # once M_BH crosses M_sigma and the wind switches on, does M_BH
    # actually stay near M_sigma (as the full King/PZNK11 picture implies),
    # or does it keep growing past it?
    #
    # This implementation deliberately only adds a gas-removal *wind* term
    # (Mdot_wind = eta_agn_eff * Mdot_BH); it does NOT add a direct cap or
    # quenching term on M_BH's own growth (PZNK11 eq. 20-23's deviation
    # term -- explicitly out of scope, see MODELS.md). Regulation, if any,
    # can only happen indirectly: the wind depletes gas_mass, which lowers
    # the gas-supply-limited growth rate (e_bh/t_ff * gas_mass) that feeds
    # M_BH itself.
    #
    # Checked directly over a long (z=25 -> z=0.5) integration with
    # sustained cosmological gas accretion and BH growth pushed hard enough
    # to reliably cross M_sigma early: that indirect loop is NOT strong
    # enough to hold M_BH near M_sigma here -- gas resupply from
    # cosmological accretion vastly outpaces what even a maximal wind
    # (eta_agn=1.0) can remove, so M_BH keeps growing over 2 orders of
    # magnitude past M_sigma by z=0.5, not staying at order-unity. This
    # pins that real (if perhaps unintuitive) behaviour explicitly, so
    # nobody mistakes "sigma_feedback" for a growth cap it doesn't
    # implement -- see MODELS.md's "Isothermal-sphere M-sigma
    # self-regulation" section for the documented caveat.
    pytest.importorskip("pymctrees")
    from ashvini import pymctrees_adapter

    # Same synthetic-config pattern as test_pymctrees_adapter.py /
    # test_run_tree_source.py -- avoids depending on where pymctrees'
    # own example config/ directory happens to live relative to the
    # installed package.
    config_path = tmp_path / "planck_like.yml"
    config_path.write_text(
        "Run:\n"
        "  mode: camb\n"
        "  pk_kmin: 1.0e-4\n"
        "  pk_kmax: 10.0\n"
        "  pk_npoints: 500\n"
        "Cosmology:\n"
        "  H0: 67.66\n"
        "  OmegaBar: 0.048\n"
        "  OmegaM: 0.3111\n"
        "  OmegaK: 0.0\n"
        "  As: 2.1e-9\n"
        "  ns: 0.9665\n"
        "  tau_reio: 0.0561\n"
        "  mnu: 0.0\n"
        "camb:\n"
    )

    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)
    monkeypatch.setattr(agn, "eta_agn", 1.0)  # maximal possible wind coupling
    monkeypatch.setattr(bh_growth, "e_bh", 0.05)
    monkeypatch.setattr(bh_growth, "eddington_multiplier", 50.0)
    monkeypatch.setattr(main, "sn_type", "delayed")

    halo_masses, halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
        pymctrees_config_path=str(config_path), mass_bin=1e12,
        n_halos=1, z0=0.5, z_max=25.0, dz=0.02,
        m_res=100.0, backend="numpy", seed=7,
    )
    result = main.run_forest(halo_masses, halo_growth_rates, redshifts)

    sigma = bh_growth.velocity_dispersion(halo_masses[0], redshifts)
    m_sigma_hist = bh_growth.m_sigma(sigma)
    ratio = np.divide(
        result["bh_mass"][0], m_sigma_hist, out=np.zeros_like(m_sigma_hist),
        where=m_sigma_hist > 0,
    )
    assert np.any(ratio > 1.0), "expected this halo to cross M_BH=M_sigma at some point"
    assert ratio[-1] > 50.0, (
        f"final M_BH/M_sigma={ratio[-1]:.3g} -- expected the wind-only "
        f"mechanism to NOT hold M_BH near M_sigma (ratio ~ 1) under "
        f"sustained gas accretion; a value close to 1 here would mean "
        f"self-regulation is happening and this test (and the MODELS.md "
        f"caveat it pins) needs revisiting, not just updating the bound"
    )


def test_run_forest_matches_scalar_reference_with_sigma_feedback_enabled(monkeypatch):
    # Mirrors test_run1.py::test_vectorized_matches_scalar_reference, but
    # with sigma_feedback enabled and BH growth pushed to a
    # dynamically-significant scale (as in
    # test_bh_growth_feedback.py::test_run_forest_with_agn_delay_differs_from_instantaneous)
    # so the M-sigma coupling switch actually has room to matter within the
    # synthetic halo's growth history, rather than being pinned at one
    # extreme throughout.
    monkeypatch.setattr(agn, "sigma_feedback_enabled", True)
    monkeypatch.setattr(bh_growth, "e_bh", 0.3)
    monkeypatch.setattr(bh_growth, "eddington_multiplier", 500.0)

    n = 60
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, utils.time_at_z(z))

    fast = main.run_forest(halo_mass, halo_mass_rate, z)
    slow = main.run1_scalar(halo_mass, halo_mass_rate, z)

    rtol = {"bh_mass": 0.03}
    default_rtol = 0.02
    for key in ("gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass", "bh_mass", "sfr"):
        fast_val = fast[key][0, -1]
        slow_val = slow[key][-1]
        assert fast_val == pytest.approx(slow_val, rel=rtol.get(key, default_rtol)), (
            f"{key}: vectorised {fast_val!r} diverged from scalar reference "
            f"{slow_val!r} with sigma_feedback enabled"
        )

    assert np.all(np.isfinite(fast["gas_mass"]))
    assert np.all(fast["gas_mass"] >= 0)
