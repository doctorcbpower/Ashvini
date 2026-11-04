"""Likelihood, prior and coverage on a closed-form model relation (docs/FIT_AUDIT.md, section 7, tests 1 to 3). No trees and no sampler."""
import numpy as np
import pytest
from scipy import stats

from ashvini.fit import Curve, Data, Problem, grid_posterior, log_likelihood, n_excluded

NODES = np.array([8.0, 9.0, 10.0, 11.0])


def curve(**kw):
    base = dict(x=NODES, median=np.array([6.0, 7.0, 8.5, 9.0]), scatter=np.array([0.1, 0.2, 0.1, 0.3]), median_err=np.array([0.0, 0.1, 0.0, 0.2]))
    base.update(kw)
    return Curve(**base)


def test_likelihood_equals_the_gaussian_terms_by_hand_including_interpolation_a_limit_and_x_errors():
    c = curve()
    d = Data(log_x=[8.5, 9.5, 10.5, 9.0], log_y=[6.8, 7.5, 8.0, 7.2], sigma_y=[0.3, 0.2, 0.4, 0.15], upper_limit=[False, False, True, False],
             sigma_x=[0.0, 0.1, 0.0, 0.0])
    s_int = 0.25
    # point 1: x = 8.5 -> median 6.5, scatter 0.15, error 0.05, slope 1.0
    v1 = 0.3 ** 2 + 0.15 ** 2 + s_int ** 2 + 0.05 ** 2
    t1 = -0.5 * ((6.8 - 6.5) ** 2 / v1 + np.log(2 * np.pi * v1))
    # point 2: x = 9.5 -> median 7.75, scatter 0.15, error 0.05, slope 1.5, x error 0.1 carried by the slope
    v2 = 0.2 ** 2 + (1.5 * 0.1) ** 2 + 0.15 ** 2 + s_int ** 2 + 0.05 ** 2
    t2 = -0.5 * ((7.5 - 7.75) ** 2 / v2 + np.log(2 * np.pi * v2))
    # point 3: upper limit at x = 10.5 -> median 8.75, scatter 0.2, error 0.1
    v3 = 0.4 ** 2 + 0.2 ** 2 + s_int ** 2 + 0.1 ** 2
    t3 = np.log(stats.norm.cdf((8.0 - 8.75) / np.sqrt(v3)))
    # point 4: exactly on a node
    v4 = 0.15 ** 2 + 0.2 ** 2 + s_int ** 2 + 0.1 ** 2
    t4 = -0.5 * ((7.2 - 7.0) ** 2 / v4 + np.log(2 * np.pi * v4))
    assert log_likelihood(d, c, s_int) == pytest.approx(t1 + t2 + t3 + t4, rel=1e-12)


def test_points_outside_the_nodes_are_left_out_and_counted_not_extrapolated():
    c = curve()
    inside = Data([9.0, 10.0], [7.0, 8.5], [0.2, 0.2])
    outside = Data([9.0, 10.0, 7.0, 11.5], [7.0, 8.5, 0.0, 99.0], [0.2, 0.2, 0.2, 0.2])
    assert n_excluded(outside, c) == 2 and n_excluded(inside, c) == 0
    assert log_likelihood(outside, c, 0.1) == log_likelihood(inside, c, 0.1)
    assert log_likelihood(Data([5.0], [1.0], [0.1]), c, 0.1) == -np.inf


def test_a_curve_from_a_relation_uses_only_unflagged_bins(monkeypatch):
    from ashvini.relation import relation_from_arrays
    rng = np.random.default_rng(3)
    n = 4000
    mh = rng.uniform(10.5, 13.0, n)
    x = 10 ** (mh - 2 + rng.normal(0, 0.2, n))
    y = 10 ** (0.9 * (np.log10(x) - 9) + 6 + rng.normal(0, 0.1, n))
    rel = relation_from_arrays(x, y, np.ones(n), 10 ** mh, np.arange(7.5, 11.5, 0.25), n_boot=10, log_mh_range=(10.5, 13.0))
    c = Curve.from_relation(rel)
    assert c.x.min() >= rel.x_mean[rel.usable].min() - 1e-12 and len(c.x) == int(np.sum(rel.usable & np.isfinite(rel.median)))
    mu, s, e, slope = c.at(c.x[len(c.x) // 2])
    assert slope == pytest.approx(0.9, abs=0.1) and s == pytest.approx(0.1, abs=0.05)


def test_data_rejects_bad_input():
    with pytest.raises(ValueError):
        Data([1.0, 2.0], [1.0], [0.1])
    with pytest.raises(ValueError):
        Data([1.0], [np.nan], [0.1])
    with pytest.raises(ValueError):
        Data([1.0], [1.0], [-0.1])


def synthetic_problem(data, fit_s_int=True):
    """Median = a + 1.5 (x - 9) at a constant model scatter 0.1: parameter a (linear, -2 to 2), plus s_int."""
    def predict(theta):
        return Curve(NODES - 0.5, theta[0] + 1.5 * (NODES - 0.5 - 9.0), np.full(4, 0.1))
    return Problem(["a"], [-2.0], [2.0], [False], predict, data, fit_s_int)


def test_the_prior_maps_log_and_linear_variables_and_respects_the_bounds():
    p = Problem(["lin", "lg"], [-2.0, 1e-4], [3.0, 1e-1], [False, True], lambda t: curve(), Data([9.0], [7.0], [0.2]))
    assert p.labels == ["lin", "lg", "s_int"]
    v = p.to_variables([0.5, 1e-3, 0.1])
    assert v == pytest.approx([0.5, -3.0, -1.0])
    theta, s = p.to_physical(v)
    assert theta == pytest.approx([0.5, 1e-3]) and s == pytest.approx(0.1)
    assert p.lower == pytest.approx([-2.0, -4.0, -2.0]) and p.upper == pytest.approx([3.0, -1.0, 0.0])
    draws = p.sample_prior(500, seed=1)
    assert all(p.contains(d) for d in draws) and draws.min(axis=0) == pytest.approx(p.lower, abs=0.2 * (p.upper - p.lower).max())
    assert p.log_prob([4.0, -3.0, -1.0]) == -np.inf and p.log_prob([0.5, -5.0, -1.0]) == -np.inf and p.log_prob([0.5, -3.0, 0.5]) == -np.inf
    assert np.isfinite(p.log_prob(v))
    with pytest.raises(ValueError):
        Problem(["lg"], [-1.0], [1.0], [True], lambda t: curve(), Data([9.0], [7.0], [0.2]))


def mock(a, s_int, rng, n=30, limit_x=None):
    """Data from the synthetic model: model scatter 0.1, intrinsic s_int, measurement error 0.3; points below a detection threshold become limits."""
    x = rng.uniform(8.0, 10.4, n)
    truth = Curve(NODES - 0.5, a + 1.5 * (NODES - 0.5 - 9.0), np.full(4, 0.1))
    mu = np.interp(x, truth.x, truth.median)
    sigma = np.full(n, 0.3)
    y = mu + rng.normal(0, np.sqrt(0.1 ** 2 + s_int ** 2), n) + rng.normal(0, sigma)
    threshold = 1.5 * (x - 9.0) + rng.normal(0, 0.3, n)                 # a detection limit set at the middle of the prior range; it does not depend on the truth
    limit = y < threshold
    return Data(x, np.where(limit, threshold, y), sigma, limit)


def test_coverage_on_a_closed_form_model_with_truths_drawn_from_the_prior():
    """150 mocks, truths drawn from the prior: the truth is in the 68% interval about 68% of the time and its rank is uniform (docs/FIT_AUDIT.md, test 3)."""
    rng = np.random.default_rng(11)
    prob0 = synthetic_problem(Data([9.0], [7.0], [0.2]))
    n_mock = 150
    inside = np.zeros((n_mock, 2), dtype=bool)
    rank = np.zeros((n_mock, 2))
    n_limits = 0
    for k in range(n_mock):
        v = prob0.sample_prior(1, seed=1000 + k)[0]
        theta, s_int = prob0.to_physical(v)
        d = mock(theta[0], s_int, rng)
        n_limits += d.upper_limit.sum()
        post = grid_posterior(synthetic_problem(d), n=(81, 41))
        for i in range(2):
            lo, hi = post.interval(i, 0.68)
            inside[k, i] = lo <= v[i] <= hi
            rank[k, i] = post.cdf_at(i, v[i])
    assert 0.2 * 30 * n_mock < n_limits < 0.8 * 30 * n_mock      # the mocks mix detections and upper limits
    sd = np.sqrt(0.68 * 0.32 / n_mock)
    for i in range(2):
        assert abs(inside[:, i].mean() - 0.68) < 3.5 * sd, (i, inside[:, i].mean())
        assert stats.kstest(rank[:, i], "uniform").pvalue > 0.01, (i, stats.kstest(rank[:, i], "uniform"))


def test_the_posterior_is_centred_on_the_truth_for_a_large_clean_sample():
    rng = np.random.default_rng(5)
    d = mock(0.3, 0.2, rng, n=400)
    post = grid_posterior(synthetic_problem(d), n=(161, 61))
    a_lo, a_hi = post.interval(0, 0.95)
    s_lo, s_hi = post.interval(1, 0.95)
    assert a_lo < 0.3 < a_hi and a_hi - a_lo < 0.2
    assert 10 ** s_lo < 0.2 < 10 ** s_hi


def test_data_entirely_outside_the_curve_is_an_error_not_a_posterior():
    """Nothing to fit: the grid says so instead of returning a number."""
    d = Data([2.0, 3.0], [1.0, 1.0], [0.1, 0.1])
    with pytest.raises(ValueError, match="likelihood is zero everywhere"):
        grid_posterior(synthetic_problem(d), n=(11, 5))


def test_grid_and_problem_are_deterministic():
    rng = np.random.default_rng(2)
    d = mock(-0.2, 0.1, rng)
    a = grid_posterior(synthetic_problem(d), n=(31, 11)).prob
    b = grid_posterior(synthetic_problem(d), n=(31, 11)).prob
    assert np.array_equal(a, b)
    p = synthetic_problem(d)
    assert np.array_equal(p.sample_prior(5, seed=4), p.sample_prior(5, seed=4))


# ---------------------------------------------------------------------------------------------------------------------------------------------
# emcee (docs/FIT_AUDIT.md, tests 4 to 7 on the closed-form model; the real-model recovery is run on OzSTAR)

from ashvini.fit import fit, problem_for_model


def test_the_fit_recovers_the_truth_and_reports_its_diagnostics():
    d = mock(0.3, 0.2, np.random.default_rng(5), n=200)
    res = fit(synthetic_problem(d), n_steps=3000, seed=3)
    assert res.converged, res.summary()
    assert res.chain.shape[2] == 2 and res.chain.shape[1] == 2 * 16 and np.all(np.isfinite(res.tau)) and 0.1 < res.acceptance < 0.9
    assert np.all(res.r_hat < 1.1)
    a = res.physical()[:, 0]
    assert np.quantile(a, 0.025) < 0.3 < np.quantile(a, 0.975)
    grid = grid_posterior(synthetic_problem(d), n=(161, 81))
    lo, hi = grid.interval(0, 0.68)
    assert np.quantile(res.samples()[:, 0], [0.16, 0.84]) == pytest.approx([lo, hi], abs=0.03)       # agrees with the exact posterior
    lo, hi = grid.interval(1, 0.68)                                                                  # (one mock's s_int can sit outside its own 95% interval)
    assert np.quantile(res.samples()[:, 1], [0.16, 0.84]) == pytest.approx([lo, hi], abs=0.15)
    assert "R-hat" in res.summary() and "tau" in res.summary()


def test_a_chain_that_is_too_short_is_reported_not_trusted():
    d = mock(0.3, 0.2, np.random.default_rng(5), n=200)
    with pytest.warns(UserWarning, match="chain too short"):
        res = fit(synthetic_problem(d), n_steps=60, seed=3)
    assert not res.converged and any("too short" in w for w in res.warnings)


def test_the_same_input_gives_the_same_chain_and_another_seed_does_not():
    d = mock(0.1, 0.1, np.random.default_rng(8), n=50)
    a = fit(synthetic_problem(d), n_steps=120, seed=4).chain
    b = fit(synthetic_problem(d), n_steps=120, seed=4).chain
    c = fit(synthetic_problem(d), n_steps=120, seed=5).chain
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_a_truth_outside_the_bounds_piles_up_against_the_bound_with_a_warning():
    d = mock(3.5, 0.1, np.random.default_rng(9), n=150)               # the prior allows -2 to 2
    with pytest.warns(UserWarning, match="piles up against the upper bound"):
        res = fit(synthetic_problem(d), n_steps=800, seed=1)
    assert np.mean(res.chain[..., 0] > 1.9) > 0.5


def test_a_model_that_never_reaches_the_data_is_an_error():
    d = Data([2.0, 3.0], [1.0, 1.0], [0.1, 0.1])
    with pytest.raises(ValueError, match="starting points"):
        fit(synthetic_problem(d), n_steps=20)


class FakeModel:
    """A stand-in for ashvini.Model: M_BH = 10^theta * M*^1 times log-normal noise, from a fixed set of 'trees'; counts the bootstrap draws asked for."""
    names = ["eff"]
    bounds = np.array([[1e-4, 1e-1]])
    log = np.array([True])

    def __init__(self, seed=0, n=3000):
        rng = np.random.default_rng(seed)
        self.mh = rng.uniform(10.5, 13.0, n)
        self.x = 10 ** (self.mh - 2.0 + rng.normal(0, 0.2, n))
        self.noise = rng.normal(0, 0.1, n)
        self.n_boot_seen = []

    def relation(self, sample, theta=None, y="bh_mass", x="stars_mass", bins=None, n_boot=200, **kw):
        from ashvini.relation import relation_from_arrays
        self.n_boot_seen.append(n_boot)
        eff = 1e-3 if theta is None else float(np.asarray(theta)[0])
        yy = eff * self.x * 10 ** self.noise
        return relation_from_arrays(self.x, yy, np.ones(len(yy)), 10 ** self.mh, np.arange(7.5, 11.5, 0.25), n_boot=n_boot, min_ess=20,
                                    log_mh_range=(10.5, 13.0))


def test_the_model_problem_holds_the_reference_errors_fixed_and_recovers_a_known_knob():
    model = FakeModel()
    truth = 4e-3
    rng = np.random.default_rng(21)
    xs = rng.uniform(8.8, 10.4, 40)
    ys = np.log10(truth) + xs + rng.normal(0, 0.25, 40)
    data = Data(xs, ys, np.full(40, 0.15))
    prob = problem_for_model(model, None, data, n_boot=30)
    assert model.n_boot_seen == [30]                                   # one bootstrap, at the reference point
    assert prob.names == ["eff"] and prob.log.tolist() == [True] and prob.lower[0] == pytest.approx(-4.0)
    lp = [prob.log_prob(np.array([v, -0.8])) for v in np.linspace(-3.5, -1.5, 41)]
    assert model.n_boot_seen[1:] == [0] * 41                           # no bootstrap inside the likelihood
    assert np.linspace(-3.5, -1.5, 41)[int(np.argmax(lp))] == pytest.approx(np.log10(truth), abs=0.15)
    res = fit(prob, n_steps=3000, seed=2)
    q = np.quantile(res.physical()[:, 0], [0.025, 0.975])
    assert q[0] < truth < q[1]


def test_forked_workers_give_the_same_chain_as_one_process():
    d = mock(0.1, 0.1, np.random.default_rng(8), n=50)
    with pytest.warns(UserWarning):
        a = fit(synthetic_problem(d), n_steps=60, seed=4, n_runs=1).chain
        b = fit(synthetic_problem(d), n_steps=60, seed=4, n_runs=1, workers=2).chain
    assert np.array_equal(a, b)
    with pytest.raises(ValueError):
        fit(synthetic_problem(d), workers=2, pool=object())
