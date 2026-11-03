"""TreeSample: abundance weights and the stellar-mass selection they enable (docs/TREESAMPLE_AUDIT.md)."""
import os

import numpy as np
import pytest

from ashvini.tree_sample import (BH_RELIABLE_HALO_MASS, TreeSample, abundance_weights, draw_log_masses, halo_range_for_stellar_masses,
                                 weighted_quantile, DEFAULT_CONFIG)

SLOPE = 0.9


def power_law(m):
    """dn/dlnM = M^-0.9 (Msun): steep enough to show the Eddington bias, with a known integral."""
    return np.asarray(m, dtype=float) ** -SLOPE


def test_draws_are_stratified_reproducible_and_inside_the_range():
    a, b = draw_log_masses(9.0, 12.0, 300, seed=3), draw_log_masses(9.0, 12.0, 300, seed=3)
    assert np.array_equal(a, b) and not np.array_equal(a, draw_log_masses(9.0, 12.0, 300, seed=4))
    edges = np.linspace(9.0, 12.0, 301)
    assert np.all((a >= edges[:-1]) & (a < edges[1:]))


def test_weights_sum_to_the_number_density_in_the_range():
    lo, hi = 9.0, 12.0
    w = abundance_weights(draw_log_masses(lo, hi, 4000, 1), power_law, lo, hi)
    exact = (10 ** (-SLOPE * lo) - 10 ** (-SLOPE * hi)) / SLOPE
    assert w.sum() == pytest.approx(exact, rel=2e-3)


def _selected_mean_log_mh(log_mh_range, n, sigma, f, log_mstar, half_width, weighted, seed=5):
    lo, hi = log_mh_range
    x = draw_log_masses(lo, hi, n, seed)
    w = abundance_weights(x, power_law, lo, hi) if weighted else np.ones(n)
    log_ms = x + np.log10(f) + np.random.default_rng(seed + 1).normal(0.0, sigma, n)
    sel = abs(log_ms - log_mstar) < half_width
    return np.average(x[sel], weights=w[sel])


def test_weighted_selection_reproduces_the_eddington_bias_and_equal_weights_do_not():
    """M* = f M_h with log-normal scatter sigma, mass function M^-0.9: at fixed M* the mean log M_h is y - 0.9 ln(10) sigma^2, y = log(M*/f)."""
    sigma, f, log_ms = 0.4, 1e-2, 9.0
    y = log_ms - np.log10(f)
    expected = y - SLOPE * np.log(10) * sigma ** 2
    kw = dict(n=400_000, sigma=sigma, f=f, log_mstar=log_ms, half_width=0.03)
    weighted = _selected_mean_log_mh((9.0, 13.0), weighted=True, **kw)
    unweighted = _selected_mean_log_mh((9.0, 13.0), weighted=False, **kw)
    assert weighted == pytest.approx(expected, abs=0.02)
    assert unweighted == pytest.approx(y, abs=0.02)
    assert unweighted - weighted > 0.25


def test_the_selected_relation_does_not_depend_on_how_wide_the_halo_range_is():
    kw = dict(n=400_000, sigma=0.4, f=1e-2, log_mstar=9.0, half_width=0.03, weighted=True)
    narrow = _selected_mean_log_mh((9.5, 12.5), **kw)
    wide = _selected_mean_log_mh((8.5, 13.5), **{**kw, "n": 800_000})
    assert narrow == pytest.approx(wide, abs=0.02)


def test_a_range_that_is_too_narrow_does_bias_the_result():
    """The check the relation step needs: cutting the range inside the selection changes the answer."""
    kw = dict(n=400_000, sigma=0.4, f=1e-2, log_mstar=9.0, half_width=0.03, weighted=True)
    assert abs(_selected_mean_log_mh((10.8, 13.0), **kw) - _selected_mean_log_mh((9.0, 13.0), **kw)) > 0.1


def test_weighted_quantile():
    x = np.arange(1.0, 101.0)
    assert weighted_quantile(x, np.ones(100), 0.5) == pytest.approx(50.5, abs=0.6)
    w = np.ones(100); w[:10] = 100.0       # most of the weight on the lowest ten values
    assert weighted_quantile(x, w, 0.5) < 10.0


def test_sample_properties():
    m0 = np.array([1e10, 5e10, 1e11, 1e12])
    s = TreeSample(m0, np.array([4.0, 2.0, 1.0, 1.0]))
    assert s.n == 4 and s.number_density == 8.0
    assert list(s.bh_reliable) == [False, False, True, True] and BH_RELIABLE_HALO_MASS == 1e11
    assert s.effective_n == pytest.approx(64.0 / 22.0)
    with pytest.raises(ValueError, match="same shape"):
        TreeSample(m0, np.ones(3))
    with pytest.raises(ValueError, match="no tree histories"):
        s.run(None)


def test_halo_range_for_stellar_masses_is_ordered_and_widens_with_margin():
    lo, hi = halo_range_for_stellar_masses((8.0, 10.5), margin_dex=1.0)
    lo0, hi0 = halo_range_for_stellar_masses((8.0, 10.5), margin_dex=0.0)
    assert lo < lo0 < hi0 < hi and (lo0, hi0) == pytest.approx((10.57, 11.95), abs=0.02)
    assert halo_range_for_stellar_masses((4.0, 5.0))[0] == 8.2        # clamped to the mass function's table


def test_v1_is_z_zero_only():
    with pytest.raises(ValueError, match="z = 0 only"):
        TreeSample.build((10.0, 11.0), z=1.0)


@pytest.mark.skipif(not os.path.exists(DEFAULT_CONFIG), reason="the foraois config is not available")
def test_build_cache_and_run(tmp_path, monkeypatch):
    from ashvini.model import Model
    from ashvini import pymctrees_adapter
    kw = dict(n_trees=12, seed=7, cache_dir=str(tmp_path), z_max=10.0, dz=0.05)
    s = TreeSample.build((10.0, 12.0), **kw)
    assert s.halo_masses.shape == (12, len(s.redshifts)) and s.growth_rates.shape == s.halo_masses.shape
    assert np.allclose(s.halo_masses[:, -1], s.m0, rtol=1e-6)              # every tree ends at its own final mass
    assert np.all(np.diff(s.m0) > 0) and np.all(np.diff(s.weight) < 0)      # heavier haloes are rarer
    assert len(list(tmp_path.glob("trees_*.npz"))) == 1
    monkeypatch.setattr(pymctrees_adapter, "build_forest_masses", lambda *a, **k: pytest.fail("rebuilt instead of read from the cache"))
    again = TreeSample.build((10.0, 12.0), **kw)
    assert np.array_equal(again.halo_masses, s.halo_masses) and np.array_equal(again.weight, s.weight)
    monkeypatch.undo()
    TreeSample.build((10.0, 12.0), **{**kw, "seed": 8})
    assert len(list(tmp_path.glob("trees_*.npz"))) == 2                      # a different seed is a different sample
    out = s.run(Model("fiducial", free=["sf_efficiency"]), [0.03])
    assert out["stars_mass"].shape == s.halo_masses.shape and np.array_equal(out["weight"], s.weight)
