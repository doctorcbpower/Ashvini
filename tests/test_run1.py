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
    # Recaptured again after fixing a real gas-conservation gap: BH
    # accretion was computed from gas_mass (via black_hole_growth_rate)
    # but never subtracted from it -- only the *AGN wind* (eta_agn x BH
    # growth rate) affected the gas reservoir, not the accreted mass
    # itself. gas_evolve.update_gas_reservoir now takes a separate
    # bh_accretion_rate (always instantaneous, a genuine sink) alongside
    # agn_wind_growth_rate (the value the wind term uses, which can be
    # delayed -- see PARAMS.bh.feedback_delay_time / agn_delay_time in
    # main.py, mirroring the SN delayed-feedback mechanism; defaults to
    # 0.0, instantaneous, matching prior behaviour for the wind). For this
    # halo bh_mass is small relative to gas_mass throughout, so the shift
    # here is tiny (~0.003%) -- the effect is much larger for haloes where
    # BH growth is a significant fraction of the gas budget (see the BH
    # on/off comparison notebook).
    "gas_mass": 16773599.010262107,
    "stars_mass": 1328598.9040095299,
    "gas_metals": 60997.94576822442,
    "stars_metals": 2858.1035534377115,
    "dust_mass": 3008.9971514939098,
    "bh_mass": 998.4929263552812,
    "sfr": 1018092.4660432874,
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


def test_sfr_computed_from_clamped_gas_mass_not_raw_ode_output():
    # Regression test for a real bug: sfr used to be computed from
    # gas_mass's raw _linear_ode_step output *before* the "enforce
    # non-negativity" clamp ran, rather than after. The closed-form ODE
    # update legitimately overshoots slightly negative when decay is fast
    # relative to dt (checked directly against a real 10,000-halo run:
    # 188,702/20,000,000 (N x n_steps) raw gas_mass updates went negative)
    # -- gas_mass itself was still correctly clamped for every *other*
    # purpose, but the already-computed sfr silently carried the negative
    # value forward. That negative sfr later fed into a delayed-feedback
    # or dust decay term (ML*wind_sfr/gas_mass, both otherwise
    # non-negative), making the decay negative and overflowing
    # _linear_ode_step's exp(-x) a few steps later -- the source of the
    # rare NaN-halo cases seen in large parameter sweeps.
    #
    # Forces every _linear_ode_step call to return a large negative value
    # regardless of its real inputs (deterministic, no need to reverse-
    # engineer parameters that organically trigger the overshoot) and
    # checks sfr is 0 (from the clamped gas_mass), never negative.
    # Verified this test fails on the pre-fix code (sfr went to ~-1.4e11)
    # via `git stash` before adding the fix.
    from ashvini import main as amain

    orig = amain._linear_ode_step

    def force_negative(y0, forcing, decay, dt):
        return orig(y0, forcing, decay, dt) - 1e12

    amain._linear_ode_step = force_negative
    try:
        N, n = 3, 10
        redshift = np.linspace(10, 5, n)
        halo_mass = np.full((N, n), 1e10)
        halo_mass_rate = np.full((N, n), 1e9)
        result = amain.run_forest(halo_mass, halo_mass_rate, redshift)
    finally:
        amain._linear_ode_step = orig

    assert np.all(result["gas_mass"] == 0.0)
    assert np.all(result["sfr"] == 0.0)
    assert not np.any(result["sfr"] < 0)


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


def test_delay_lookback_index_uses_actual_elapsed_time():
    # Regression test for a real, pre-existing bug (confirmed present in
    # the original solas-sims/Ashvini repo too, not introduced by this
    # fork's vectorisation work): the delayed-feedback lookback used to be
    # a fixed array-index offset calibrated once early in the run, reused
    # for every later step regardless of how much cosmic_time it actually
    # spanned there. That's only correct for a uniform cosmic_time grid --
    # this directly checks _delay_lookback_index against a deliberately
    # non-uniform one (mimicking a redshift-uniform tree, where dt/dz
    # varies enormously with z) and confirms the *time gap*, not the step
    # count, stays close to t_d throughout.
    from ashvini.main import _delay_lookback_index

    # non-uniform cosmic_time: fine steps early, coarse steps late (same
    # qualitative shape as a redshift-uniform grid's z->t mapping)
    cosmic_time = np.concatenate([
        np.linspace(0.1, 0.3, 400),
        np.linspace(0.3, 1.2, 400)[1:],
    ])
    t_d = 0.015  # Gyr, matches run_params.yaml's default delay_time

    delay_idx = _delay_lookback_index(cosmic_time, t_d)

    has_history = delay_idx >= 0
    assert has_history.any() and (~has_history).any()  # exercises both branches

    actual_gap = cosmic_time[has_history] - cosmic_time[delay_idx[has_history]]
    # searchsorted picks the last tabulated point at or before the target
    # lookback time, so the true gap is always >= t_d, bounded above by
    # t_d + one local step (otherwise the *next* point would have been
    # picked instead) -- everywhere, including late in the run where steps
    # are ~20x coarser than at the start. The old fixed-index-offset bug
    # would instead have shown a gap that grew to tens of times t_d there.
    assert np.all(actual_gap >= t_d - 1e-12)
    local_step = np.diff(cosmic_time)[np.clip(delay_idx[has_history], 0, len(cosmic_time) - 2)]
    assert np.all(actual_gap <= t_d + local_step + 1e-9)


def test_time_redshift_roundtrip_consistency():
    # Sanity check on the interpolation swap in utils.py (cubic interp1d ->
    # np.interp): converting z -> t -> z should be close to the identity
    # over the grid's normal operating range.
    z = np.linspace(5.0, 20.0, 50)
    t = utils.time_at_z(z)
    z_roundtrip = utils.z_at_time(t)
    assert np.allclose(z, z_roundtrip, atol=1e-3)
