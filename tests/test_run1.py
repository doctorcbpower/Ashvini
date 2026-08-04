"""
Regression tests for the core ashvini pipeline.

These run against tests/fixtures/merger_trees_fixture.h5 -- a heavily
downsampled slice of the real merger-tree data (5 haloes x 3 mass bins x
301 redshift steps, ~60 KB) generated with scripts/downsample_trees.py.
See that script's docstring for how to regenerate or extend the fixture.

The point of this suite is to catch unintended behaviour changes -- e.g. from
the ODE-vectorisation work flagged in the code audit as the next efficiency
step -- not to validate the astrophysics itself. Reference values below were
captured from the current implementation and are pinned with a loose enough
tolerance (rtol=1e-6) to survive floating point / platform differences, but
tight enough to catch a real change in behaviour.
"""

import os

import numpy as np
import pytest

from ashvini import utils
from ashvini.main import run1

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "merger_trees_fixture.h5")

OUTPUT_KEYS = ["gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass", "sfr"]

# Reference "final timestep" values for halo 0 in the 1e10 mass bin, captured
# from the current implementation (see module docstring).
REFERENCE_FINAL_VALUES = {
    "gas_mass": 16427710.954957614,
    "stars_mass": 1601779.2279374518,
    "gas_metals": 65222.22581236108,
    "stars_metals": 2848.3583806431257,
    "dust_mass": 4832.189043741113,
    "sfr": 997098.4013238181,
}


@pytest.fixture(scope="module")
def tree_data():
    halo_masses, halo_mass_rates, redshift = utils.read_trees(
        file_path=FIXTURE_PATH, mass_bin=1e10
    )
    return halo_masses, halo_mass_rates, redshift


def test_read_trees_shapes(tree_data):
    halo_masses, halo_mass_rates, redshift = tree_data
    assert halo_masses.shape == halo_mass_rates.shape
    assert halo_masses.shape[1] == redshift.shape[0]
    assert halo_masses.shape[0] > 0


def test_read_trees_mass_bin_naming():
    # Regression check on the mass_bin -> HDF5 group-name convention
    # ("01e10" for 1e10, "05e08" for 5e8, ...) that utils.read_trees relies on.
    halo_masses, _, _ = utils.read_trees(file_path=FIXTURE_PATH, mass_bin=1e10)
    assert halo_masses.shape[0] == 5  # fixture was built with --n-halos 5

    with pytest.raises(KeyError):
        utils.read_trees(file_path=FIXTURE_PATH, mass_bin=1e12)  # not in fixture


def test_run1_output_shapes_and_validity(tree_data):
    halo_masses, halo_mass_rates, redshift = tree_data
    result = run1(halo_masses[0], halo_mass_rates[0], redshift)

    n = len(redshift)
    for key in OUTPUT_KEYS:
        values = result[key]
        assert values.shape == (n,), f"{key} has unexpected shape {values.shape}"
        assert np.all(np.isfinite(values)), f"{key} contains non-finite values"
        assert np.all(values >= 0), f"{key} contains negative values"

    # gas metals and dust require gas to be present
    no_gas = result["gas_mass"] <= 0
    assert np.all(result["gas_metals"][no_gas] == 0.0)
    assert np.all(result["dust_mass"][no_gas] == 0.0)


def test_run1_matches_reference_baseline(tree_data):
    halo_masses, halo_mass_rates, redshift = tree_data
    result = run1(halo_masses[0], halo_mass_rates[0], redshift)

    for key, expected in REFERENCE_FINAL_VALUES.items():
        actual = result[key][-1]
        assert actual == pytest.approx(expected, rel=1e-6), (
            f"{key} final-step value drifted from the pinned baseline "
            f"({actual!r} vs {expected!r}). If this change is intentional "
            f"(e.g. a deliberate change to the integration scheme), update "
            f"REFERENCE_FINAL_VALUES after checking the new output makes sense."
        )


def test_time_redshift_roundtrip_consistency():
    # Sanity check on the interpolation swap in utils.py (cubic interp1d ->
    # np.interp): converting z -> t -> z should be close to the identity
    # over the grid's normal operating range.
    z = np.linspace(5.0, 20.0, 50)
    t = utils.time_at_z(z)
    z_roundtrip = utils.z_at_time(t)
    assert np.allclose(z, z_roundtrip, atol=1e-3)
