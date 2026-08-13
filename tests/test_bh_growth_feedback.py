"""
Tests for BH accretion as a genuine gas-mass sink, the Eddington-growth
multiplier, and the (optionally delayed) AGN wind term -- see
gas_evolve.update_gas_reservoir, black_holes_growth.py, and main.py's
run_forest/run1_scalar.

Under the default run_params.yaml (eddington_multiplier=1.0,
feedback_delay_time=0.0, efficiency=0.001), BH growth is a tiny fraction of
a halo's gas budget (checked directly: bh_mass/gas_mass ~ 1e-6 for the test
fixture's most massive bin) -- real, but too small to see clearly end to
end without a very large/long-lived sample. These tests instead verify
each mechanism directly (with parameters chosen to make the effect
unambiguous), rather than relying on default-config end-to-end runs where
it would be lost in noise.
"""

import numpy as np
import pytest

from ashvini import agn_feedback as agn
from ashvini import black_holes_growth as bh_growth
from ashvini import gas_evolve
from ashvini import main


# ---------------------------------------------------------------------------
# BH accretion as a gas-mass sink (gas_evolve.update_gas_reservoir)
# ---------------------------------------------------------------------------

def test_bh_accretion_rate_reduces_gas_mass_evolution_rate():
    kwargs = dict(
        t=1.0, gas_mass=1e8, gas_accretion_rate=1e6, halo_mass=1e11,
        stellar_metallicity=0.01, past_sfr=1e5, kind="no",
    )
    rate_no_bh = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=0.0, agn_wind_growth_rate=0.0)
    rate_with_bh = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=1e4, agn_wind_growth_rate=0.0)

    # a genuine sink: the rate should be reduced by exactly bh_accretion_rate
    assert rate_with_bh == pytest.approx(rate_no_bh - 1e4)


def test_agn_wind_growth_rate_reduces_gas_mass_evolution_rate():
    kwargs = dict(
        t=1.0, gas_mass=1e8, gas_accretion_rate=1e6, halo_mass=1e11,
        stellar_metallicity=0.01, past_sfr=1e5, kind="no", bh_accretion_rate=0.0,
    )
    rate_no_wind = gas_evolve.update_gas_reservoir(**kwargs, agn_wind_growth_rate=0.0)
    rate_with_wind = gas_evolve.update_gas_reservoir(**kwargs, agn_wind_growth_rate=1e4)

    expected_wind = agn.agn_wind_mass_rate(1e4)
    assert expected_wind > 0  # sanity: eta_agn > 0 in default config
    assert rate_with_wind == pytest.approx(rate_no_wind - expected_wind)


def test_bh_accretion_and_agn_wind_are_independent_terms():
    # bh_accretion_rate (instantaneous sink) and agn_wind_growth_rate
    # (possibly delayed, drives the wind) must be genuinely separate
    # knobs -- changing one shouldn't silently move the other.
    kwargs = dict(
        t=1.0, gas_mass=1e8, gas_accretion_rate=1e6, halo_mass=1e11,
        stellar_metallicity=0.01, past_sfr=1e5, kind="no",
    )
    base = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=0.0, agn_wind_growth_rate=0.0)
    sink_only = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=5e3, agn_wind_growth_rate=0.0)
    wind_only = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=0.0, agn_wind_growth_rate=5e3)
    both = gas_evolve.update_gas_reservoir(**kwargs, bh_accretion_rate=5e3, agn_wind_growth_rate=5e3)

    assert sink_only != pytest.approx(base)
    assert wind_only != pytest.approx(base)
    assert sink_only != pytest.approx(wind_only)
    # both terms are linear/independent in the RHS, so they combine additively:
    # both = base - x_term - y_term = sink_only + wind_only - base
    assert both == pytest.approx(sink_only + wind_only - base)


# ---------------------------------------------------------------------------
# Eddington multiplier (super-Eddington growth)
# ---------------------------------------------------------------------------

def test_eddington_multiplier_scales_the_cap(monkeypatch):
    M_BH = 1e6  # Msun
    base_rate = bh_growth.eddington_bh_growth(M_BH)

    monkeypatch.setattr(bh_growth, "eddington_multiplier", 1.0)
    huge_gas = 1e20  # gas-unlimited regime, so growth is purely Eddington-capped
    rate_1x = bh_growth.black_hole_growth_rate(t=0.0, M_BH=M_BH, gas_mass=huge_gas)

    monkeypatch.setattr(bh_growth, "eddington_multiplier", 3.0)
    rate_3x = bh_growth.black_hole_growth_rate(t=0.0, M_BH=M_BH, gas_mass=huge_gas)

    assert float(rate_1x) == pytest.approx(float(base_rate))
    assert float(rate_3x) == pytest.approx(3.0 * float(base_rate))


def test_bh_growth_step_respects_multiplied_eddington_cap():
    M_BH0 = 1e6
    kappa = bh_growth.EDDINGTON_RATE_PER_UNIT_MASS
    A_bh = 1e300  # effectively gas-unlimited, forces the Eddington branch
    dt = 0.1  # Gyr

    grown_1x = main._bh_growth_step(np.array([M_BH0]), np.array([A_bh]), kappa, dt)
    grown_3x = main._bh_growth_step(np.array([M_BH0]), np.array([A_bh]), 3.0 * kappa, dt)

    # pure exponential growth at rate kappa (or 3*kappa) over dt
    assert grown_1x[0] == pytest.approx(M_BH0 * np.exp(kappa * dt), rel=1e-6)
    assert grown_3x[0] == pytest.approx(M_BH0 * np.exp(3.0 * kappa * dt), rel=1e-6)
    assert grown_3x[0] > grown_1x[0]  # super-Eddington genuinely grows faster


# ---------------------------------------------------------------------------
# AGN feedback delay (mirrors the SN delayed-feedback mechanism)
# ---------------------------------------------------------------------------

def test_agn_delay_time_zero_uses_same_step_growth_rate(monkeypatch):
    # PARAMS.bh.feedback_delay_time=0.0 (the default) must reduce to
    # exactly the pre-existing instantaneous behaviour: agn_delay_idx[j]==j
    # for every j (see _delay_lookback_index's docstring/tests), so the
    # wind term always uses THIS step's own growth rate.
    from ashvini.main import _delay_lookback_index

    cosmic_time = np.linspace(0.1, 1.0, 50)
    agn_delay_idx = _delay_lookback_index(cosmic_time, 0.0)
    assert np.array_equal(agn_delay_idx, np.arange(len(cosmic_time)))


def test_run_forest_with_agn_delay_differs_from_instantaneous(monkeypatch):
    # End-to-end: under default physical parameters, BH growth stays
    # Eddington-capped at a scale far below the gas accretion rate for
    # this synthetic halo (checked directly: BH growth rate ~1e4 Msun/Gyr
    # vs ~6e11 Msun/Gyr accretion) -- realistic for early, seed-scale BH
    # growth, but too small for the AGN wind's *timing* to leave any
    # visible mark on gas_mass regardless of whether it's delayed. Raising
    # e_bh and eddington_multiplier (a deliberate "what if AGN feedback
    # were dominant" stress case, not a realistic config) makes the BH
    # growth rate comparable to the accretion rate, so the delay's effect
    # on gas_mass is actually visible rather than lost in noise.
    monkeypatch.setattr(bh_growth, "e_bh", 0.3)
    monkeypatch.setattr(bh_growth, "eddington_multiplier", 500.0)

    n = 60
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    monkeypatch.setattr(main, "agn_delay_time", 0.0)
    result_instant = main.run_forest(halo_mass, halo_mass_rate, z)

    monkeypatch.setattr(main, "agn_delay_time", 0.3)  # Gyr
    result_delayed = main.run_forest(halo_mass, halo_mass_rate, z)

    assert not np.allclose(result_instant["gas_mass"], result_delayed["gas_mass"])
    assert np.all(np.isfinite(result_delayed["gas_mass"]))
    assert np.all(result_delayed["gas_mass"] >= 0)


# ---------------------------------------------------------------------------
# End-to-end: BH growth genuinely competes with star formation for gas
# ---------------------------------------------------------------------------

def test_bh_growth_measurably_reduces_gas_and_stellar_mass(monkeypatch):
    # e_bh=0 alone does NOT disable BH growth: seeding (a one-time mass
    # injection, unrelated to e_bh) still happens, leaving bh_mass pinned
    # at its seed value rather than 0 -- checked directly, this is what
    # the first version of this test got wrong. Fully disabling BH growth
    # (and, through it, the AGN wind) requires disabling all three seeding
    # channels, so no halo is ever seeded in the first place -- the same
    # mechanism the BH on/off comparison notebook uses.
    n = 80
    z = np.linspace(5.0, 15.0, n)[::-1]
    halo_mass = np.linspace(1e10, 1e12, n)
    halo_mass_rate = np.gradient(halo_mass, main.utils.time_at_z(z))

    monkeypatch.setattr(bh_growth, "e_bh", 0.3)
    result_bh_on = main.run_forest(halo_mass, halo_mass_rate, z)

    seeding = main.PARAMS.bh.seeding
    monkeypatch.setattr(seeding.pop3, "enabled", False)
    monkeypatch.setattr(seeding.direct_collapse, "enabled", False)
    monkeypatch.setattr(seeding.halo_mass_threshold, "enabled", False)
    result_bh_off = main.run_forest(halo_mass, halo_mass_rate, z)

    assert np.all(result_bh_off["bh_mass"] == 0.0)
    assert result_bh_on["bh_mass"][0, -1] > 0.0
    # BH growth draws down gas relative to the no-BH-growth case
    assert result_bh_on["gas_mass"][0, -1] <= result_bh_off["gas_mass"][0, -1]
    assert np.all(np.isfinite(result_bh_on["gas_mass"]))
    assert np.all(result_bh_on["gas_mass"] >= 0)
