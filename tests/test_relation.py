"""relation(): y given x for a weighted tree sample, on fake galaxy models with analytic answers (docs/RELATION_AUDIT.md, section 7)."""
import os
import warnings

import numpy as np
import pytest

from ashvini.relation import relation, relation_from_arrays
from ashvini.tree_sample import DEFAULT_CONFIG, TreeSample, abundance_weights, draw_log_masses

SLOPE = 0.9
F = 1e-2                         # M* = F * M_h (times scatter)
BINS = np.arange(7.5, 10.6, 0.25)


def power_law(m):
    return np.asarray(m, dtype=float) ** -SLOPE


def sample(n=60_000, lo=9.0, hi=13.0, seed=2, weighted=True):
    x = draw_log_masses(lo, hi, n, seed)
    w = abundance_weights(x, power_law, lo, hi) if weighted else np.ones(n)
    return x, w


def galaxies(n=60_000, sigma=0.4, seed=2, **kw):
    log_mh, w = sample(n, seed=seed, **kw)
    eps = np.random.default_rng(seed + 10).normal(0.0, sigma, n)
    return log_mh, w, 10 ** (log_mh + np.log10(F) + eps)


def test_a_known_power_law_is_recovered_with_its_scatter():
    log_mh, w, mstar = galaxies()
    sigma_y, b, g = 0.3, 1.2, 1e-3
    y = g * mstar ** b * 10 ** np.random.default_rng(1).normal(0.0, sigma_y, len(mstar))     # independent of everything but M*
    rel = relation_from_arrays(mstar, y, w, 10 ** log_mh, BINS, n_boot=20, seed=1, log_mh_range=(9.0, 13.0))
    ok = rel.usable & (rel.ess > 500)
    assert ok.sum() >= 6
    assert np.allclose(rel.median[ok], np.log10(g) + b * rel.x_mean[ok] + 0.0, atol=0.03 + 3 * rel.median_err[ok])
    assert np.allclose(rel.scatter[ok], sigma_y, atol=0.03 + 3 * rel.scatter_err[ok])


def test_the_relation_includes_the_eddington_shift_and_equal_weights_miss_it():
    """y = M_h^1.5 exactly: log y given M* is 1.5 (log(M*/F) - 0.9 ln(10) sigma^2) for a power-law mass function."""
    sigma = 0.4
    k = 1.5
    expected = lambda x: k * (x - np.log10(F) - SLOPE * np.log(10) * sigma ** 2)
    for weighted, should_match in ((True, True), (False, False)):
        log_mh, w, mstar = galaxies(n=200_000, sigma=sigma, weighted=weighted)
        rel = relation_from_arrays(mstar, 10 ** (k * log_mh), w, 10 ** log_mh, np.arange(8.5, 9.6, 0.25), n_boot=10, min_ess=50)
        ok = rel.ess > 2000
        assert ok.sum() >= 3
        err = np.abs(rel.median[ok] - expected(rel.x_mean[ok]))
        assert (err.max() < 0.05) if should_match else (err.min() > 0.3)


def test_occupation_fraction_is_recovered_and_empty_trees_do_not_change_the_rest():
    log_mh, w, mstar = galaxies()
    lx = np.log10(mstar)
    p = np.clip((lx - 7.0) / 4.0, 0.05, 0.95)
    occupied = np.random.default_rng(3).random(len(lx)) < p
    y_all = 10 ** (lx - 2.0 + np.random.default_rng(4).normal(0, 0.2, len(lx)))
    full = relation_from_arrays(mstar, y_all, w, 10 ** log_mh, BINS, n_boot=10)
    part = relation_from_arrays(mstar, np.where(occupied, y_all, 0.0), w, 10 ** log_mh, BINS, n_boot=10)
    ok = full.ess > 500
    assert np.allclose(part.f_occ[ok], np.clip((part.x_mean[ok] - 7.0) / 4.0, 0.05, 0.95), atol=0.03)
    assert np.allclose(part.median[ok], full.median[ok], atol=0.03) and np.allclose(full.f_occ[ok], 1.0)


def test_too_few_effective_trees_gives_nan_and_a_flag():
    log_mh, w, mstar = galaxies(n=60)
    rel = relation_from_arrays(mstar, mstar, w, 10 ** log_mh, BINS, n_boot=5)
    assert rel.low_stats.any() and np.all(np.isnan(rel.median[rel.low_stats]))
    assert np.all(np.isnan(rel.f_occ[rel.ess < 20]))


def test_range_flag_fires_when_the_halo_range_is_cut_inside_the_selection_and_not_otherwise():
    for lo, expect in ((10.8, True), (9.0, False)):
        log_mh, w, mstar = galaxies(lo=lo, hi=13.0)
        rel = relation_from_arrays(mstar, mstar, w, 10 ** log_mh, np.array([8.9, 9.1]), n_boot=3, log_mh_range=(lo, 13.0))
        assert bool(rel.range_limited[0]) is expect


def test_black_hole_flag_and_warning():
    log_mh, w, mstar = galaxies(sigma=0.2)
    with pytest.warns(UserWarning, match="below 1e\\+11 Msun"):
        rel = relation_from_arrays(mstar, mstar, w, 10 ** log_mh, np.array([8.4, 8.6, 9.8, 10.0]), n_boot=3, bh=True)
    assert rel.bh_unreliable[0] and not rel.bh_unreliable[2]            # M* near 4e8 (halo 3e10) against M* near 7e9 (halo 7e11)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        relation_from_arrays(mstar, mstar, w, 10 ** log_mh, np.array([8.4, 8.6]), n_boot=3, bh=False)    # not a BH quantity: no warning


def test_bootstrap_error_scales_with_sample_size_and_is_reproducible():
    errs = []
    for n in (20_000, 80_000):
        log_mh, w, mstar = galaxies(n=n, seed=5)
        rel = relation_from_arrays(mstar, mstar ** 1.1, w, 10 ** log_mh, np.array([8.9, 9.1]), n_boot=60, seed=4)
        errs.append(rel.median_err[0])
    assert errs[0] / errs[1] == pytest.approx(2.0, rel=0.4)
    a = relation_from_arrays(mstar, mstar, w, 10 ** log_mh, np.array([8.9, 9.1]), n_boot=15, seed=9)
    b = relation_from_arrays(mstar, mstar, 7.0 * w, 10 ** log_mh, np.array([8.9, 9.1]), n_boot=15, seed=9)
    assert np.allclose(a.median, b.median) and np.allclose(a.median_err, b.median_err) and np.allclose(a.f_occ, b.f_occ)


def test_x_scale_shifts_the_bins_exactly():
    log_mh, w, mstar = galaxies(n=20_000)
    y = mstar ** 1.2
    a = relation_from_arrays(mstar, y, w, 10 ** log_mh, BINS, x_scale=2.0, n_boot=3)
    b = relation_from_arrays(mstar, y, w, 10 ** log_mh, BINS - np.log10(2.0), x_scale=1.0, n_boot=3)
    assert np.allclose(a.median, b.median, equal_nan=True) and np.allclose(a.x_mean, b.x_mean + np.log10(2.0), equal_nan=True)


def test_stars_mass_carries_the_note_and_functions_can_define_x_and_y():
    log_mh, w, mstar = galaxies(n=2000)
    result = {"stars_mass": np.stack([mstar, mstar], axis=1), "bh_mass": np.stack([mstar, 2 * mstar], axis=1), "weight": w, "m0": 10 ** log_mh}
    rel = relation(result, "bh_mass", "stars_mass", bins=np.array([8.9, 9.1]), n_boot=3, min_ess=1)
    assert "no stellar mass loss" in rel.metadata["note"] and "no stellar mass loss" in rel.table()
    ratio = relation(result, lambda r: r["bh_mass"] / r["stars_mass"], lambda r: r["stars_mass"], bins=np.array([8.9, 9.1]), n_boot=3, min_ess=1)
    assert ratio.median[0] == pytest.approx(np.log10(2.0), abs=1e-9)        # last step: 2 M* / M* = 2


@pytest.mark.skipif(not os.path.exists(DEFAULT_CONFIG), reason="the foraois config is not available")
def test_model_relation_on_a_real_tree_sample(tmp_path):
    from ashvini.model import Model
    s = TreeSample.build((10.0, 12.5), n_trees=12, seed=7, cache_dir=str(tmp_path), z_max=10.0, dz=0.05)
    bins = np.array([4.0, 7.0, 11.0])
    rel = Model("fiducial", free=["sf_efficiency"]).relation(s, [0.03], y="gas_mass", bins=bins, min_ess=1, n_boot=5)
    out = s.run(Model("fiducial", free=["sf_efficiency"]), [0.03])
    x, w = out["stars_mass"][:, -1], s.weight
    inside = (x > 0) & (np.log10(np.maximum(x, 1e-30)) >= 7.0) & (np.log10(np.maximum(x, 1e-30)) < 11.0)
    assert len(rel.median) == 2
    if inside.any():
        assert rel.f_occ[1] == pytest.approx(w[inside & (out["gas_mass"][:, -1] > 0)].sum() / w[inside].sum())
