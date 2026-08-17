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


def teardown_function(_):
    # restore fiducial values so other test modules (which read these
    # globals at call time, not import time) aren't affected by ordering.
    configure.set_params(
        eps_sf=0.015, eps_fb=5, t_d=0.015, z_rei=7,
        sn_type="delayed", uv_background=True,
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
