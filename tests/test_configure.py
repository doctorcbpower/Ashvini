"""
Tests for ashvini.configure.set_params -- see its module docstring for why
this exists (run_params.yaml's values are cached into module-level globals
at import time, scattered across several modules, some of them duplicated).
"""
import numpy as np

from ashvini import configure
from ashvini import main as amain
from ashvini import star_formation as asf
from ashvini import supernovae_feedback as asn
from ashvini import reionization as areion
from ashvini import agn_feedback as aagn
from ashvini import black_holes_growth as abh
from ashvini.run_params import PARAMS


def teardown_function(_):
    # restore fiducial values so other test modules (which read these
    # globals at call time, not import time) aren't affected by ordering.
    configure.set_params(
        eps_sf=0.015, eps_fb=5, t_d=0.015, z_rei=7,
        sn_type="delayed", uv_background=True,
        bh_enabled=True, eta_agn=0.5, sigma_feedback_enabled=False,
        growth_cap_enabled=False, e_bh=0.001, eddington_multiplier=1.0,
    )


def test_eps_sf_patches_both_star_formation_and_main_copies():
    # main.py duplicates star_formation.e_ff as its own module-level e_ff
    # (used directly in run_forest's inlined SFR update) -- this is exactly
    # the "one of two copies not updated" bug class this test guards
    # against, given the codebase's history with silent parameter drift.
    configure.set_params(eps_sf=0.1)
    assert asf.e_ff == 0.1
    assert amain.e_ff == 0.1


def test_eps_fb_patches_supernovae_feedback():
    configure.set_params(eps_fb=7)
    assert asn.epsilon_p == 7


def test_t_d_and_sn_type_and_uv_patch_main():
    configure.set_params(t_d=0.03, sn_type="instantaneous", uv_background=False)
    assert amain.t_d == 0.03
    assert amain.sn_type == "instantaneous"
    assert amain.UV_background is False


def test_z_rei_repatches_derived_beta_and_c_omega():
    beta_fiducial = areion.beta
    configure.set_params(z_rei=10)
    assert areion.z_rei == 10
    assert areion.beta != beta_fiducial
    # sanity: recomputing beta from scratch with the same formula agrees
    gamma = areion.gamma
    expected_beta = 10 * (
        (np.log(1.82 * (10**3) * np.exp(-0.63 * 10) - 1)) ** (-1 / gamma)
    )
    assert areion.beta == expected_beta


def test_unspecified_arguments_leave_values_unchanged():
    configure.set_params(eps_sf=0.05)
    eps_fb_before = asn.epsilon_p
    configure.set_params(eps_sf=0.06)
    assert asn.epsilon_p == eps_fb_before


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


def test_eta_agn_patches_agn_feedback():
    configure.set_params(eta_agn=1.0)
    assert aagn.eta_agn == 1.0


def test_sigma_feedback_enabled_patches_agn_feedback():
    configure.set_params(sigma_feedback_enabled=True)
    assert aagn.sigma_feedback_enabled is True


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


def test_growth_cap_enabled_patches_black_holes_growth():
    configure.set_params(growth_cap_enabled=True)
    assert abh.growth_cap_enabled is True


def test_e_bh_patches_black_holes_growth():
    configure.set_params(e_bh=0.05)
    assert abh.e_bh == 0.05


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
