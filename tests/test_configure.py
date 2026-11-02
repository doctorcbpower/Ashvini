"""
Tests for ashvini.configure.set_params -- see its module docstring: it changes values of config.DEFAULT in
place, and nothing is cached in module-level copies, so every module sees the change on its next call.
"""
import numpy as np

from ashvini import configure
from ashvini import main as amain
from ashvini import supernovae_feedback as asn
from ashvini import reionization as areion
from ashvini import reionization_physics as areion_physics
from ashvini.run_params import PARAMS


def teardown_function(_):
    # restore fiducial values so other test modules (which read these
    # globals at call time, not import time) aren't affected by ordering.
    configure.set_params(
        eps_sf=0.015, eps_fb=5, t_d=0.015, z_rei=7,
        sn_type="delayed", uv_background=True,
        bh_enabled=True, eta_agn=0.5, sigma_feedback_enabled=False,
        growth_cap_enabled=False, e_bh=0.001, eddington_multiplier=1.0,
        m_sn_max=40.0, agn_delay_time=0.0,
    )


def test_eps_sf_patches_the_config():
    from ashvini.star_formation import star_formation_rate
    configure.set_params(eps_sf=0.1)
    assert PARAMS.sf.efficiency == 0.1
    ratio = star_formation_rate(1.0, 1.0) / star_formation_rate(1.0, 1.0, cfg=_with(sf_efficiency=1.0))
    assert abs(ratio - 0.1) < 1e-12


def _with(sf_efficiency):
    from ashvini.config import override
    return override(PARAMS, **{"sf.efficiency": sf_efficiency})


def test_eps_fb_patches_supernovae_feedback():
    configure.set_params(eps_fb=7)
    assert PARAMS.sn.epsilon_p == 7
    assert float(asn.mass_loading_factor(6.0, 1e9, 0.0)) > 0
    base = float(asn.mass_loading_factor(6.0, 1e9, 0.0))
    configure.set_params(eps_fb=14)
    assert float(asn.mass_loading_factor(6.0, 1e9, 0.0)) == 2 * base


def test_t_d_and_sn_type_and_uv_patch_the_config():
    configure.set_params(t_d=0.03, sn_type="instantaneous", uv_background=False)
    assert PARAMS.sn.delay_time == 0.03
    assert PARAMS.sn.type == "instantaneous"
    assert PARAMS.reion.UVB_enabled is False


def test_z_rei_rederives_beta_and_c_omega():
    beta_fiducial = areion.beta
    configure.set_params(z_rei=10)
    assert areion.z_rei == 10
    assert areion.beta != beta_fiducial
    # beta is derived from (z_rei, gamma) on each read, never stored: it can not go stale
    assert areion.beta == areion_physics.beta(10, areion.gamma)


def test_unspecified_arguments_leave_values_unchanged():
    configure.set_params(eps_sf=0.05)
    eps_fb_before = PARAMS.sn.epsilon_p
    configure.set_params(eps_sf=0.06)
    assert PARAMS.sn.epsilon_p == eps_fb_before


def test_set_params_actually_changes_run_forest_output():
    # end-to-end check: a lower eps_sf should produce a lower SFR (and
    # hence less stellar mass) for the same halo history. halo_mass must
    # actually grow at halo_mass_rate (not just be paired with an
    # arbitrary nonzero rate) or gas_inflow_rate has nothing consistent to
    # accrete and stars_mass stays identically zero regardless of eps_sf.
    from ashvini import utils
    n = 50
    z = np.linspace(15, 5, n)  # chronological (increasing cosmic time)
    t = utils.time_at_z(z)
    halo_mass_rate = np.full((3, n), 1e8)  # Msun/Gyr
    halo_mass = 1e8 + np.cumsum(halo_mass_rate * np.gradient(t), axis=1)

    configure.set_params(eps_sf=0.015)
    result_fiducial = amain.run_forest(halo_mass, halo_mass_rate, z)

    configure.set_params(eps_sf=0.1)
    result_high_eff = amain.run_forest(halo_mass, halo_mass_rate, z)

    assert result_high_eff["stars_mass"][:, -1].sum() > result_fiducial["stars_mass"][:, -1].sum()


def test_bh_enabled_false_disables_all_seeding_channels():
    configure.set_params(bh_enabled=False)
    seeding = PARAMS.bh.seeding
    assert seeding.pop3.enabled is False
    assert seeding.direct_collapse.enabled is False
    assert seeding.halo_mass_threshold.enabled is False


def test_bh_enabled_true_restores_all_seeding_channels():
    configure.set_params(bh_enabled=False)
    configure.set_params(bh_enabled=True)
    seeding = PARAMS.bh.seeding
    assert seeding.pop3.enabled is True
    assert seeding.direct_collapse.enabled is True
    assert seeding.halo_mass_threshold.enabled is True


def test_eta_agn_patches_the_config():
    configure.set_params(eta_agn=1.0)
    assert PARAMS.bh.eta_agn == 1.0


def test_sigma_feedback_enabled_patches_the_config():
    configure.set_params(sigma_feedback_enabled=True)
    assert PARAMS.bh.sigma_feedback.enabled is True


def test_bh_enabled_false_actually_prevents_seeding_in_run_forest():
    # end-to-end check: a halo well above every seeding threshold should
    # still never grow a BH once bh_enabled=False, since
    # seeding_mask_and_mass() reads PARAMS.bh.seeding fresh every call.
    from ashvini import utils
    n = 200
    z = np.linspace(15, 5, n)
    t = utils.time_at_z(z)
    halo_mass_rate = np.full((2, n), 1e10)  # Msun/Gyr, well above every threshold quickly
    halo_mass = 1e11 + np.cumsum(halo_mass_rate * np.gradient(t), axis=1)

    configure.set_params(bh_enabled=True)
    result_on = amain.run_forest(halo_mass, halo_mass_rate, z)
    assert np.any(result_on["bh_mass"] > 0)

    configure.set_params(bh_enabled=False)
    result_off = amain.run_forest(halo_mass, halo_mass_rate, z)
    assert np.all(result_off["bh_mass"] == 0)


def test_growth_cap_enabled_patches_the_config():
    configure.set_params(growth_cap_enabled=True)
    assert PARAMS.bh.sigma_feedback.growth_cap_enabled is True


def test_e_bh_patches_the_config():
    configure.set_params(e_bh=0.05)
    assert PARAMS.bh.efficiency == 0.05


def test_m_sn_max_and_agn_delay_patch_the_config():
    configure.set_params(m_sn_max=100.0, agn_delay_time=0.2)
    assert PARAMS.sn.m_sn_max == 100.0
    assert PARAMS.bh.feedback_delay_time == 0.2


def test_e_bh_actually_changes_run_forest_bh_growth():
    # e_bh only matters when growth is gas-supply-limited, not Eddington-
    # limited (checked directly: with the default eddington_multiplier=1,
    # a freshly-seeded BH's Eddington cap is so far below the gas-supply
    # rate's implied threshold_mass that growth stays Eddington-limited --
    # and therefore e_bh-independent -- for the whole of a typical test
    # scenario). eddington_multiplier is pushed way up here so gas supply,
    # not the Eddington cap, is actually the binding constraint.
    from ashvini import utils
    n = 200
    z = np.linspace(15, 5, n)
    t = utils.time_at_z(z)
    halo_mass_rate = np.full((2, n), 1e10)
    halo_mass = 1e11 + np.cumsum(halo_mass_rate * np.gradient(t), axis=1)

    configure.set_params(e_bh=0.001, eddington_multiplier=10)
    result_low = amain.run_forest(halo_mass, halo_mass_rate, z)

    configure.set_params(e_bh=0.05, eddington_multiplier=10)
    result_high = amain.run_forest(halo_mass, halo_mass_rate, z)

    assert result_high["bh_mass"][:, -1].sum() > result_low["bh_mass"][:, -1].sum()
