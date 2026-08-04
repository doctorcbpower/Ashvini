"""
Regression tests for the core ashvini pipeline.

These run against tests/fixtures/merger_trees_fixture.h5 -- a heavily
downsampled slice of the real merger-tree data (5 haloes x 3 mass bins x
301 redshift steps, ~60 KB) generated with scripts/downsample_trees.py.
See that script's docstring for how to regenerate or extend the fixture.

The point of this suite is to catch unintended behaviour changes -- not to
validate the astrophysics itself. Reference values below were captured from
the current implementation and are pinned with a loose enough tolerance
(rtol=1e-6) to survive floating point / platform differences, but tight
enough to catch a real change in behaviour.

run1() was rewritten (see ashvini/main.py) to replace its five per-step
solve_ivp calls with a closed-form/quadrature vectorised update -- see the
module-level comment in main.py for the numerical scheme. A sixth quantity,
bh_mass, was added later (BH seeding + Eddington-limited growth + AGN
feedback on gas_mass, see ashvini/black_holes_growth.py and
ashvini/agn_feedback.py). REFERENCE_FINAL_VALUES below was recaptured from
that implementation; run1_scalar() (the old solve_ivp-based version, kept as
a reference) is cross-checked against it in
test_vectorized_matches_scalar_reference below.
"""

import os

import numpy as np
import pytest

from ashvini import utils
from ashvini.main import run1, run1_scalar
from ashvini.run_params import PARAMS

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "merger_trees_fixture.h5")

OUTPUT_KEYS = [
    "gas_mass",
    "stars_mass",
    "gas_metals",
    "stars_metals",
    "dust_mass",
    "bh_mass",
    "sfr",
]

# bh_mass carries a somewhat larger cross-validation tolerance than the other
# outputs: its piecewise (Eddington-limited vs gas-supply-limited) closed
# form is a coarser approximation than the affine closed forms used for the
# other ODEs -- see _bh_growth_step in main.py.
CROSS_VALIDATION_RTOL = {"bh_mass": 0.03}
DEFAULT_CROSS_VALIDATION_RTOL = 0.02

# Reference "final timestep" values for halo 0 in the 1e10 mass bin, captured
# from the current (vectorised) run1() implementation (see module docstring).
REFERENCE_FINAL_VALUES = {
    "gas_mass": 16460203.827404512,
    "stars_mass": 1610384.6462875996,
    "gas_metals": 65431.49820079876,
    "stars_metals": 2870.1699206575877,
    "dust_mass": 3288.939936002951,
    "bh_mass": 998.4932457898989,
    "sfr": 999070.5927788572,
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


def test_vectorized_matches_scalar_reference(tree_data):
    # Cross-check the closed-form/quadrature run1() against the old
    # solve_ivp-based run1_scalar() it replaced. Not bit-for-bit -- run1()
    # freezes redshift-dependent ODE coefficients at each step's midpoint
    # instead of letting LSODA integrate them continuously -- but final-step
    # values should agree to within a percent or so for every halo in the
    # fixture. A larger drift here would indicate a real bug in the
    # vectorised scheme, not just approximation error.
    halo_masses, halo_mass_rates, redshift = tree_data
    for i in range(halo_masses.shape[0]):
        fast = run1(halo_masses[i], halo_mass_rates[i], redshift)
        slow = run1_scalar(halo_masses[i], halo_mass_rates[i], redshift)
        for key in OUTPUT_KEYS:
            rtol = CROSS_VALIDATION_RTOL.get(key, DEFAULT_CROSS_VALIDATION_RTOL)
            assert fast[key][-1] == pytest.approx(slow[key][-1], rel=rtol), (
                f"halo {i}, {key}: vectorised final value {fast[key][-1]!r} "
                f"diverged from scalar reference {slow[key][-1]!r} by more "
                f"than the expected approximation error"
            )


def test_bh_seeding_and_growth(tree_data):
    from ashvini import black_holes_growth as bh_growth

    halo_masses, halo_mass_rates, redshift = tree_data

    valid_seed_masses = {0.0}
    for channel in (
        PARAMS.bh.seeding.pop3,
        PARAMS.bh.seeding.direct_collapse,
        PARAMS.bh.seeding.halo_mass_threshold,
    ):
        if channel.enabled:
            valid_seed_masses.add(channel.M_seed)

    for i in range(halo_masses.shape[0]):
        result = run1(halo_masses[i], halo_mass_rates[i], redshift)
        bh_mass = result["bh_mass"]

        # BH mass can only grow (seeding jumps up from 0, then Eddington-
        # limited/gas-supply-limited growth is a strictly non-negative rate)
        assert np.all(np.diff(bh_mass) >= -1e-8), f"halo {i}: bh_mass decreased"

        if np.any(bh_mass > 0):
            first_nonzero = np.argmax(bh_mass > 0)
            # the seeding jump itself should land on one of the configured
            # seed masses (to a loose tolerance -- it may have grown a
            # little within the same step it was seeded)
            seed_value = bh_mass[first_nonzero]
            assert any(
                seed_value == pytest.approx(m, rel=0.05) or seed_value >= m
                for m in valid_seed_masses
                if m > 0
            ), f"halo {i}: seed value {seed_value!r} doesn't match any configured M_seed"

        # never exceeds the Eddington rate integrated naively over the
        # whole history (a very loose sanity bound, not a tight physical one)
        assert np.all(bh_mass < 1e12), f"halo {i}: bh_mass unphysically large"


def test_time_redshift_roundtrip_consistency():
    # Sanity check on the interpolation swap in utils.py (cubic interp1d ->
    # np.interp): converting z -> t -> z should be close to the identity
    # over the grid's normal operating range.
    z = np.linspace(5.0, 20.0, 50)
    t = utils.time_at_z(z)
    z_roundtrip = utils.z_at_time(t)
    assert np.allclose(z, z_roundtrip, atol=1e-3)
