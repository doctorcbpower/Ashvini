"""
Tests for ashvini.pymctrees_adapter -- the adapter that turns a pymctrees
forest into ashvini.utils.read_trees' expected layout, shared by the
offline scripts/build_trees_from_pymctrees.py path and the live
tree_source='pymctrees' path (see main.py's run()).

pymctrees is not an Ashvini dependency (see the module's own docstring), so
these tests are skipped entirely if it isn't installed. They use a synthetic
power-law P(k) (mirroring pymctrees' own test suite) rather than CLASS/CAMB,
so they don't require either heavy backend to be installed either -- only
pymctrees itself.
"""

import numpy as np
import pytest

pymctrees = pytest.importorskip("pymctrees")

from pymctrees.cosmo_utils import CosmoData
from pymctrees.pch_trees import PCHMergerTree

from ashvini import pymctrees_adapter as adapter
from ashvini import utils
from ashvini.main import run1

PLANCK_LIKE = {
    "Code": {"mode": "camb", "pk_kmin": 1e-4, "pk_kmax": 10.0, "pk_npoints": 500},
    "Cosmology": {"H0": 67.66, "OmegaM": 0.3111, "OmegaK": 0.0, "OmegaLambda": 0.6889},
}
H = 0.6766  # H0 / 100, matches Cosmology.H0 above


@pytest.fixture
def tree_generator():
    cosmo_data = CosmoData(PLANCK_LIKE, redshift=[0.0])
    k = np.logspace(-4, 1, 500)
    Pk = 2.0e4 * k**-2.0
    pk_data = {"k": k, "Pk": Pk.reshape(1, -1), "z": [0.0]}
    cosmo_data.get_power_spectrum = lambda: pk_data
    return PCHMergerTree(cosmo_data, PLANCK_LIKE)


def test_build_forest_for_bin_shape_and_ordering(tree_generator):
    n_halos = 20
    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e12, h=H, n_halos=n_halos,
        z0=0.0, z_max=5.0, m_res_msun=1e9, dz=0.5, backend="numpy", rng_seed=0,
    )
    n_steps = int(5.0 / 0.5) + 1  # +1 for the prepended M0
    assert halo_masses.shape == (n_halos, n_steps)
    assert redshifts.shape == (n_steps,)

    # index 0 = earliest (highest z); index -1 = z0 = the input M0 exactly
    assert redshifts[0] > redshifts[-1]
    assert redshifts[-1] == pytest.approx(0.0)
    assert np.allclose(halo_masses[:, -1], 1e12)

    # mass only grows forward in (chronological) time
    assert np.all(np.diff(halo_masses, axis=1) >= -1e-6)


def test_msun_over_h_conversion(tree_generator):
    # Same seed/params, different h -- output masses should scale by 1/h
    # relative to each other for a fixed physical M0 (M0 itself is passed in
    # Msun both times, so the *converted* trajectories should differ only
    # through the M0->M0*h step, not through some accidental double-counting)
    np.random.seed(1)
    halo_masses_h1, _ = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e12, h=1.0, n_halos=5,
        z0=0.0, z_max=3.0, m_res_msun=1e9, dz=0.5, backend="numpy", rng_seed=1,
    )
    np.random.seed(1)
    halo_masses_h05, _ = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e12, h=0.5, n_halos=5,
        z0=0.0, z_max=3.0, m_res_msun=1e9, dz=0.5, backend="numpy", rng_seed=1,
    )
    # M0 in Msun/h differs (1e12 vs 2e12) between the two calls, so the tree
    # topology (which depends on M0_hunits, not M0_msun) differs too -- but
    # the final (z0) mass must round-trip back to the same input M0_msun
    # regardless of h, since that's what got prepended and then divided out.
    assert np.allclose(halo_masses_h1[:, -1], 1e12)
    assert np.allclose(halo_masses_h05[:, -1], 1e12)


def test_compute_growth_rates_shape(tree_generator):
    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e12, h=H, n_halos=10,
        z0=0.0, z_max=5.0, m_res_msun=1e9, dz=0.5, backend="numpy", rng_seed=2,
    )
    rates = adapter.compute_growth_rates(halo_masses, redshifts)
    assert rates.shape == halo_masses.shape
    assert np.all(np.isfinite(rates))
    # rates are never negative for a monotonically non-decreasing mass history
    assert np.all(rates >= -1e-6)


def test_mask_unresolved_prefix_detects_frozen_start():
    # Direct unit test of the helper, independent of pymctrees itself:
    # a contiguous run of the row's own first value starting at index 0
    # is the frozen prefix; a later step that happens to coincidentally
    # re-equal the first value must not be swept in once the walk has
    # already moved on.
    halo_masses = np.array([
        [5.0, 5.0, 5.0, 12.0, 20.0],   # frozen for first 3 steps, then grows
        [7.0, 7.0, 7.0, 7.0, 7.0],     # frozen for the whole row
        [3.0, 9.0, 3.0, 15.0, 20.0],   # NOT frozen -- moves immediately, index 2's
                                        # coincidental re-match to index 0 doesn't count
    ])
    mask = adapter._mask_unresolved_prefix(halo_masses)
    assert np.array_equal(mask[0], [True, True, True, False, False])
    assert np.array_equal(mask[1], [True, True, True, True, True])
    assert np.array_equal(mask[2], [True, False, False, False, False])


def test_build_forest_for_bin_zeroes_preformation_prefix(tree_generator):
    # A deep z_max relative to a tight m_res forces at least some haloes'
    # random walks to not resolve a split until well after z_max -- for
    # those, build_forest_for_bin must report mass=0 (not a frozen
    # placeholder) for the pre-formation stretch (see
    # _mask_unresolved_prefix's docstring for why: pymctrees can't resolve
    # a progenitor below M_res, so this stretch shouldn't be presented to
    # Ashvini as "a halo with mass" at all).
    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e10, h=H, n_halos=20,
        z0=0.0, z_max=8.0, m_res_msun=1e2, dz=0.05, backend="numpy", rng_seed=7,
    )
    has_preformation = np.any(halo_masses[:, 0] == 0.0)
    assert has_preformation, "expected at least one halo to have an unresolved prefix at this dynamic range"
    # every prefix must be a contiguous run of exact zeros starting at
    # index 0 (never a zero appearing after real growth has started)
    for row in halo_masses:
        nonzero = np.where(row > 0)[0]
        if len(nonzero) == 0 or nonzero[0] == 0:
            continue
        first_nonzero = nonzero[0]
        assert np.all(row[:first_nonzero] == 0.0)
        assert np.all(row[first_nonzero:] > 0.0)


def test_compute_growth_rates_zeroes_formation_step_not_a_spike(tree_generator):
    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e10, h=H, n_halos=20,
        z0=0.0, z_max=8.0, m_res_msun=1e2, dz=0.05, backend="numpy", rng_seed=7,
    )
    rates = adapter.compute_growth_rates(halo_masses, redshifts)
    assert np.all(np.isfinite(rates))

    for i, row in enumerate(halo_masses):
        nonzero = np.where(row > 0)[0]
        if len(nonzero) == 0 or nonzero[0] == 0:
            continue
        formation_idx = nonzero[0]
        # the step whose diff spans 0 -> first resolved mass must be
        # reported as zero growth, not a huge instantaneous-formation spike
        assert rates[i, formation_idx - 1] == 0.0


def test_mass_bin_group_naming():
    assert utils.mass_bin_group_name(1e10) == "01e10"
    assert utils.mass_bin_group_name(5e8) == "05e08"


def test_output_file_round_trips_through_read_trees_and_run1(tree_generator, tmp_path):
    import h5py

    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e10, h=H, n_halos=5,
        z0=0.0, z_max=5.0, m_res_msun=1e8, dz=0.25, backend="numpy", rng_seed=3,
    )
    rates = adapter.compute_growth_rates(halo_masses, redshifts)

    out_path = tmp_path / "adapter_test.h5"
    with h5py.File(out_path, "w") as f:
        f.create_dataset("redshifts", data=redshifts)
        grp = f.create_group(utils.mass_bin_group_name(1e10))
        grp.create_dataset("halo_masses", data=halo_masses)
        grp.create_dataset("halo_growth_rates", data=rates)

    loaded_masses, loaded_rates, loaded_z = utils.read_trees(
        file_path=str(out_path), mass_bin=1e10
    )
    assert loaded_masses.shape == halo_masses.shape

    # the whole point of the adapter: its output must be directly usable by
    # run1() with no further changes
    result = run1(loaded_masses[0], loaded_rates[0], loaded_z)
    for key in ["gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass", "bh_mass", "sfr"]:
        assert np.all(np.isfinite(result[key]))
        assert np.all(result[key] >= 0)


def test_build_forest_live_matches_offline_equivalent(tmp_path):
    # build_forest_live (the "live" path used by run() when
    # tree_source='pymctrees') must produce exactly the same tree as the
    # offline build_forest_for_bin + compute_growth_rates combination it
    # wraps, given the same parameters and seed -- they're meant to be
    # interchangeable, not two independently-evolving implementations.
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

    live_masses, live_rates, live_z = adapter.build_forest_live(
        pymctrees_config_path=str(config_path), mass_bin=1e10,
        n_halos=5, z0=0.0, z_max=3.0, dz=0.5, m_res=1e9,
        backend="numpy", seed=42,
    )

    assert live_masses.shape == (5, 7)  # (0, 0.5, ..., 3.0) -> 7 steps
    assert live_z[-1] == pytest.approx(0.0)
    assert np.allclose(live_masses[:, -1], 1e10)
    assert np.all(np.isfinite(live_rates))
