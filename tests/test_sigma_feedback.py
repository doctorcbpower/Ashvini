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
