"""
Tests for scripts/build_trees_from_pymctrees.py -- the adapter that turns a
pymctrees forest into ashvini.utils.read_trees' expected HDF5 layout (code
audit, Section 6).

pymctrees is not an Ashvini dependency (see the script's own docstring), so
these tests are skipped entirely if it isn't installed. They use a synthetic
power-law P(k) (mirroring pymctrees' own test suite) rather than CLASS/CAMB,
so they don't require either heavy backend to be installed either -- only
pymctrees itself.
"""

import importlib.util
import os

import numpy as np
import pytest

pymctrees = pytest.importorskip("pymctrees")

from pymctrees.cosmo_utils import CosmoData
from pymctrees.pch_trees import PCHMergerTree

from ashvini import utils
from ashvini.main import run1

SCRIPT_PATH = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "build_trees_from_pymctrees.py"
)


@pytest.fixture(scope="module")
def adapter():
    spec = importlib.util.spec_from_file_location("build_trees_from_pymctrees", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_build_forest_for_bin_shape_and_ordering(adapter, tree_generator):
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


def test_msun_over_h_conversion(adapter, tree_generator):
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


def test_compute_growth_rates_shape(adapter, tree_generator):
    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e12, h=H, n_halos=10,
        z0=0.0, z_max=5.0, m_res_msun=1e9, dz=0.5, backend="numpy", rng_seed=2,
    )
    rates = adapter.compute_growth_rates(halo_masses, redshifts)
    assert rates.shape == halo_masses.shape
    assert np.all(np.isfinite(rates))
    # rates are never negative for a monotonically non-decreasing mass history
    assert np.all(rates >= -1e-6)


def test_mass_bin_group_naming(adapter):
    assert adapter._mass_bin_group_name(1e10) == "01e10"
    assert adapter._mass_bin_group_name(5e8) == "05e08"


def test_output_file_round_trips_through_read_trees_and_run1(adapter, tree_generator, tmp_path):
    import h5py

    halo_masses, redshifts = adapter.build_forest_for_bin(
        tree_generator, M0_msun=1e10, h=H, n_halos=5,
        z0=0.0, z_max=5.0, m_res_msun=1e8, dz=0.25, backend="numpy", rng_seed=3,
    )
    rates = adapter.compute_growth_rates(halo_masses, redshifts)

    out_path = tmp_path / "adapter_test.h5"
    with h5py.File(out_path, "w") as f:
        f.create_dataset("redshifts", data=redshifts)
        grp = f.create_group(adapter._mass_bin_group_name(1e10))
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
