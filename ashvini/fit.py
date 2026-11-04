"""
Observed relations and the likelihood of a model relation (docs/FIT_AUDIT.md, interface steps 4 and 5).

    data = Data(log_x, log_y, sigma_y, upper_limit=flags)          # measurements of an observed relation, in log10
    curve = Curve.from_relation(rel)                                # the model's median, scatter and errors at the centres of its unflagged bins
    ll = log_likelihood(data, curve, s_int=0.2)

Each point is Gaussian in log10 y about the model median at its x, with variance
    sigma_y^2 + (slope * sigma_x)^2 + model scatter^2 + s_int^2 + model median error^2:
the measurement error, the error in x carried into y by the local slope of the model relation, the scatter the model produces, an extra intrinsic
scatter (a fitted parameter), and the error of the model median from the finite number of trees. The median is interpolated linearly in log x between
bin centres; a point outside the range of the centres is excluded (and counted), never extrapolated. An upper limit contributes the probability that the
model lies below it.
"""
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.special import log_ndtr

S_INT_RANGE = (0.01, 1.0)          # dex; uniform in log (docs/FIT_AUDIT.md, S4)


@dataclass
class Data:
    log_x: np.ndarray
    log_y: np.ndarray                 # the measurement, or the limit where upper_limit is True
    sigma_y: np.ndarray               # dex
    upper_limit: np.ndarray = None    # bool; True where log_y is an upper limit on the value
    sigma_x: np.ndarray = None        # dex; error in log x, default zero

    def __post_init__(self):
        self.log_x, self.log_y = np.asarray(self.log_x, dtype=float), np.asarray(self.log_y, dtype=float)
        n = len(self.log_x)
        self.sigma_y = np.broadcast_to(np.asarray(self.sigma_y, dtype=float), (n,)).copy()
        self.upper_limit = np.zeros(n, dtype=bool) if self.upper_limit is None else np.broadcast_to(np.asarray(self.upper_limit, dtype=bool), (n,)).copy()
        self.sigma_x = np.zeros(n) if self.sigma_x is None else np.broadcast_to(np.asarray(self.sigma_x, dtype=float), (n,)).copy()
        if self.log_y.shape != (n,):
            raise ValueError("log_x and log_y must have the same length.")
        if not (np.all(np.isfinite(self.log_x)) and np.all(np.isfinite(self.log_y)) and np.all(np.isfinite(self.sigma_y)) and np.all(np.isfinite(self.sigma_x))):
            raise ValueError("Data contains NaN or infinite values.")
        if np.any(self.sigma_y < 0) or np.any(self.sigma_x < 0):
            raise ValueError("measurement errors must not be negative.")

    def __len__(self):
        return len(self.log_x)


@dataclass
class Curve:
    """The model's relation at nodes in log x: median of log y, the scatter about it, and the error of the median (all dex), nodes increasing."""
    x: np.ndarray
    median: np.ndarray
    scatter: np.ndarray
    median_err: np.ndarray = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.x, self.median, self.scatter = (np.asarray(a, dtype=float) for a in (self.x, self.median, self.scatter))
        self.median_err = np.zeros_like(self.x) if self.median_err is None else np.asarray(self.median_err, dtype=float)
        order = np.argsort(self.x)
        self.x, self.median, self.scatter, self.median_err = (a[order] for a in (self.x, self.median, self.scatter, self.median_err))

    @classmethod
    def from_relation(cls, rel, median_err=None):
        """The unflagged bins of a Relation. median_err replaces the bootstrap errors (the fit holds them fixed, from one reference run)."""
        ok = rel.usable & np.isfinite(rel.median) & np.isfinite(rel.scatter)
        err = rel.median_err if median_err is None else median_err
        return cls(rel.x_mean[ok], rel.median[ok], rel.scatter[ok], np.where(np.isfinite(err[ok]), err[ok], 0.0), dict(rel.metadata))

    def covers(self, log_x):
        """True where log_x lies between the first and last node."""
        if len(self.x) < 2:
            return np.zeros(np.shape(log_x), dtype=bool)
        return (np.asarray(log_x) >= self.x[0]) & (np.asarray(log_x) <= self.x[-1])

    def at(self, log_x):
        """(median, scatter, median error, local slope of the median) at log_x, by linear interpolation; defined only where `covers`."""
        log_x = np.asarray(log_x, dtype=float)
        slopes = np.diff(self.median) / np.diff(self.x)
        seg = np.clip(np.searchsorted(self.x, log_x, side="right") - 1, 0, len(slopes) - 1)
        return (np.interp(log_x, self.x, self.median), np.interp(log_x, self.x, self.scatter),
                np.interp(log_x, self.x, self.median_err), slopes[seg])


def log_likelihood(data, curve, s_int=0.0):
    """
    Log of the probability of `data` given the model relation `curve` and an intrinsic scatter s_int (dex). Points outside the nodes of the curve are
    left out (see `n_excluded`); -inf if no point is left.
    """
    use = curve.covers(data.log_x)
    if not use.any():
        return -np.inf
    mu, s_mod, err, slope = curve.at(data.log_x[use])
    var = data.sigma_y[use] ** 2 + (slope * data.sigma_x[use]) ** 2 + s_mod ** 2 + s_int ** 2 + err ** 2
    y, lim = data.log_y[use], data.upper_limit[use]
    detected = -0.5 * ((y - mu) ** 2 / var + np.log(2 * np.pi * var))
    censored = log_ndtr((y - mu) / np.sqrt(var))
    return float(np.sum(np.where(lim, censored, detected)))


def n_excluded(data, curve):
    """Number of data points outside the range of the curve's nodes, which the likelihood leaves out."""
    return int(np.sum(~curve.covers(data.log_x)))


@dataclass
class Problem:
    """
    A set of free parameters with uniform priors in their sampling variable, a function giving the model relation, and data.

    names, low, high, log: as for Model (log parameters are sampled in log10). `predict(theta)` returns a Curve for the physical values of the parameters.
    The last sampling variable is log10 of the intrinsic scatter s_int, uniform within S_INT_RANGE (or absent if fit_s_int is False, then s_int = 0).
    """
    names: list
    low: np.ndarray
    high: np.ndarray
    log: np.ndarray
    predict: callable
    data: Data
    fit_s_int: bool = True

    def __post_init__(self):
        self.low, self.high = np.asarray(self.low, dtype=float), np.asarray(self.high, dtype=float)
        self.log = np.asarray(self.log, dtype=bool)
        if np.any(self.low >= self.high) or np.any(self.log & (self.low <= 0)):
            raise ValueError("bounds must have low < high, and low > 0 for parameters sampled in log.")

    @property
    def labels(self):
        return list(self.names) + (["s_int"] if self.fit_s_int else [])

    @property
    def lower(self):
        """Lower bound of each sampling variable (log10 for log parameters, and log10 s_int last)."""
        v = np.where(self.log, np.log10(np.where(self.log, self.low, 1.0)), self.low)
        return np.append(v, np.log10(S_INT_RANGE[0])) if self.fit_s_int else v

    @property
    def upper(self):
        v = np.where(self.log, np.log10(np.where(self.log, self.high, 1.0)), self.high)
        return np.append(v, np.log10(S_INT_RANGE[1])) if self.fit_s_int else v

    def to_variables(self, physical):
        """Physical parameter values (and s_int last, if fitted) to sampling variables."""
        p = np.asarray(physical, dtype=float)
        n = len(self.names)
        v = np.where(self.log, np.log10(np.where(p[:n] > 0, p[:n], 1.0)), p[:n])
        return np.append(v, np.log10(p[n])) if self.fit_s_int else v

    def to_physical(self, v):
        """Sampling variables to (theta, s_int)."""
        v = np.asarray(v, dtype=float)
        n = len(self.names)
        theta = np.where(self.log, 10.0 ** v[:n], v[:n])
        return theta, (10.0 ** v[n] if self.fit_s_int else 0.0)

    def contains(self, v):
        v = np.asarray(v, dtype=float)
        return bool(np.all(v >= self.lower) and np.all(v <= self.upper))

    def log_prior(self, v):
        """Uniform within the bounds (0), -inf outside; the constant is dropped."""
        return 0.0 if self.contains(v) else -np.inf

    def log_prob(self, v):
        """Log posterior, up to a constant, at sampling variables v. -inf outside the bounds."""
        if not self.contains(v):
            return -np.inf
        theta, s_int = self.to_physical(v)
        return log_likelihood(self.data, self.predict(theta), s_int)

    def sample_prior(self, n, seed=0):
        """n draws of the sampling variables from the prior."""
        rng = np.random.default_rng(seed)
        return self.lower + rng.random((n, len(self.lower))) * (self.upper - self.lower)


@dataclass
class GridPosterior:
    axes: list                        # one array of sampling variables per dimension
    prob: np.ndarray                  # posterior probability per grid cell, summing to 1
    labels: list

    def marginal(self, i):
        other = tuple(j for j in range(len(self.axes)) if j != i)
        return self.prob.sum(axis=other) if other else self.prob

    def cdf_at(self, i, value):
        """Marginal posterior probability below `value` for dimension i (the rank of a truth), with the cell containing it split linearly."""
        m = self.marginal(i)
        edges = np.concatenate([[self.axes[i][0] - 0.5 * (self.axes[i][1] - self.axes[i][0])], 0.5 * (self.axes[i][1:] + self.axes[i][:-1]),
                                [self.axes[i][-1] + 0.5 * (self.axes[i][-1] - self.axes[i][-2])]])
        cum = np.concatenate([[0.0], np.cumsum(m)])
        return float(np.interp(value, edges, cum))

    def interval(self, i, level=0.68):
        """Central credible interval of dimension i (sampling variable): the points where the marginal CDF is (1-level)/2 and (1+level)/2."""
        m = self.marginal(i)
        a = self.axes[i]
        edges = np.concatenate([[a[0] - 0.5 * (a[1] - a[0])], 0.5 * (a[1:] + a[:-1]), [a[-1] + 0.5 * (a[-1] - a[-2])]])
        cum = np.concatenate([[0.0], np.cumsum(m)])
        return tuple(float(np.interp(q, cum, edges)) for q in ((1 - level) / 2, (1 + level) / 2))


def grid_posterior(problem, n=60):
    """
    The posterior on a regular grid of the sampling variables (n points per dimension; one dimension per parameter and s_int). For problems of one or two
    dimensions this is exact up to the grid and needs no sampler; docs/FIT_AUDIT.md, section 10.
    """
    n = [n] * len(problem.lower) if np.isscalar(n) else list(n)
    axes = [np.linspace(lo, hi, k) for lo, hi, k in zip(problem.lower, problem.upper, n)]
    mesh = np.meshgrid(*axes, indexing="ij")
    flat = np.stack([m.ravel() for m in mesh], axis=1)
    lp = np.array([problem.log_prob(v) for v in flat])
    lp = np.where(np.isfinite(lp), lp, -np.inf)
    if not np.isfinite(lp).any():
        raise ValueError("the likelihood is zero everywhere on the grid: no point of the prior range reaches the data.")
    p = np.exp(lp - lp.max())
    return GridPosterior(axes, (p / p.sum()).reshape([len(a) for a in axes]), problem.labels)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# Sampling (docs/FIT_AUDIT.md, N1 to N3)

PILE_UP_EDGE = 0.02          # fraction of a parameter's prior range at each end
PILE_UP_FRACTION = 0.10      # posterior mass in that edge above which the parameter is said to pile up (a flat posterior puts 2 per cent there)


@dataclass
class FitResult:
    labels: list
    chain: np.ndarray                 # (steps, walkers, ndim) in sampling variables, after the start-up steps are removed
    log_prob: np.ndarray              # (steps, walkers)
    tau: np.ndarray                   # integrated autocorrelation time per parameter, in steps (inf where it could not be estimated)
    acceptance: float
    r_hat: np.ndarray                 # split Gelman-Rubin statistic per parameter between independent ensembles (NaN with one ensemble)
    n_steps: int
    n_discard: int
    problem: Problem = None
    warnings: list = field(default_factory=list)

    @property
    def converged(self):
        return not self.warnings

    def samples(self, thin=None):
        """Flattened posterior samples in sampling variables, thinned by half the longest autocorrelation time by default."""
        step = thin or (max(1, int(0.5 * np.max(self.tau))) if np.all(np.isfinite(self.tau)) else 1)
        return self.chain[::step].reshape(-1, self.chain.shape[-1])

    def physical(self, thin=None):
        """Flattened samples in physical values (log parameters as values, s_int in dex)."""
        s = self.samples(thin)
        return np.array([np.append(*self.problem.to_physical(v)) for v in s])

    def summary(self, levels=(0.16, 0.5, 0.84)):
        """One line per parameter: the median and the 16th to 84th percentile range of the physical values, then the diagnostics."""
        p = self.physical()
        rows = []
        for i, name in enumerate(self.labels):
            q = np.quantile(p[:, i], levels)
            rows.append(f"{name}: {q[1]:.4g} (+{q[2] - q[1]:.3g} / -{q[1] - q[0]:.3g})   tau = {self.tau[i]:.0f} steps   R-hat = {self.r_hat[i]:.3f}")
        rows.append(f"{self.n_steps} steps ({self.n_discard} discarded), acceptance fraction {self.acceptance:.2f}")
        rows += [f"WARNING: {w}" for w in self.warnings] or ["no convergence or bound warnings"]
        return "\n".join(rows)

    def corner(self, truths=None, **kwargs):
        """A corner plot of the physical values (requires the `corner` package)."""
        import corner
        return corner.corner(self.physical(), labels=self.labels, truths=truths, **kwargs)


def _split_r_hat(runs):
    """Gelman-Rubin statistic from independent ensembles: each ensemble's flattened chain is split in two, giving 2 x len(runs) chains."""
    if len(runs) < 2:
        return np.full(runs[0].shape[-1], np.nan)
    chains = []
    for r in runs:
        flat = r.reshape(-1, r.shape[-1])
        h = len(flat) // 2
        chains += [flat[:h], flat[h:2 * h]]
    c = np.array(chains)                               # (m, n, ndim)
    m, n = c.shape[:2]
    within = c.var(axis=1, ddof=1).mean(axis=0)
    between = n * c.mean(axis=1).var(axis=0, ddof=1)
    var_plus = (n - 1) / n * within + between / n
    return np.sqrt(var_plus / within)


_ACTIVE_PROBLEM = None


def _active_log_prob(v):
    """The log posterior of the problem being fitted. Module level so that worker processes forked by `fit(workers=...)` can call it without pickling
    the problem (which holds the whole tree sample)."""
    return _ACTIVE_PROBLEM.log_prob(v)


def fit(problem, n_steps=2000, n_walkers=None, n_discard=None, n_runs=2, seed=0, pool=None, workers=None, min_tau_ratio=50, max_r_hat=1.1,
        min_accept=0.1, progress=False):
    """
    Sample the posterior of `problem` with emcee: `n_runs` independent ensembles of `n_walkers` (default 4 x the number of dimensions, at least 16),
    each started from different draws of the prior. Convergence is judged from the integrated autocorrelation time (the chain must be at least
    `min_tau_ratio` times longer), the acceptance fraction, and agreement between the ensembles (split R-hat below `max_r_hat`). A parameter whose
    posterior piles up against a bound of its prior is reported. All problems found are in `FitResult.warnings` and are also issued as warnings;
    nothing is hidden when a chain is too short.
    workers: number of forked processes for the walkers' likelihood evaluations (they share the problem and its tree sample without copying it).
    pool: instead, an object with a `map` method (an MPI pool, say) that can call `problem.log_prob`.
    """
    import emcee
    global _ACTIVE_PROBLEM
    if workers and pool is not None:
        raise ValueError("give workers or pool, not both.")
    log_prob = problem.log_prob
    if workers and workers > 1:
        import multiprocessing
        _ACTIVE_PROBLEM = problem
        pool = multiprocessing.get_context("fork").Pool(workers)
        log_prob = _active_log_prob
    ndim = len(problem.lower)
    n_walkers = n_walkers or max(16, 4 * ndim)
    n_walkers += n_walkers % 2                                    # emcee needs an even number
    n_discard = n_steps // 4 if n_discard is None else n_discard
    if n_steps <= n_discard:
        raise ValueError("n_steps must exceed n_discard.")
    runs, lps, tau_runs, acc = [], [], [], []
    own_pool = pool if (workers and workers > 1) else None
    try:
        for run in range(n_runs):
            start = _starting_points(problem, n_walkers, seed + 7919 * run)
            sampler = emcee.EnsembleSampler(n_walkers, ndim, log_prob, pool=pool)
            sampler.random_state = np.random.RandomState(seed + 104729 * run + 1).get_state()
            sampler.run_mcmc(start, n_steps, progress=progress)
            runs.append(sampler.get_chain(discard=n_discard))
            lps.append(sampler.get_log_prob(discard=n_discard))
            tau_runs.append(sampler.get_autocorr_time(discard=n_discard, quiet=True, tol=0))
            acc.append(float(np.mean(sampler.acceptance_fraction)))
    finally:
        if own_pool is not None:
            own_pool.close()
            own_pool.join()
    chain = np.concatenate(runs, axis=1)
    tau = np.max(tau_runs, axis=0)
    r_hat = _split_r_hat(runs)
    n_kept = n_steps - n_discard
    messages = []
    if not np.all(np.isfinite(tau)) or np.any(n_kept < min_tau_ratio * tau):
        messages.append(f"chain too short: {n_kept} kept steps against an autocorrelation time of up to {np.max(tau):.0f} steps "
                        f"(need at least {min_tau_ratio} times that). Run more steps; do not trust the intervals.")
    if np.mean(acc) < min_accept or np.mean(acc) > 0.9:
        messages.append(f"acceptance fraction {np.mean(acc):.2f} is outside {min_accept} to 0.9.")
    if n_runs > 1 and np.any(r_hat > max_r_hat):
        bad = ", ".join(f"{n} ({r:.2f})" for n, r in zip(problem.labels, r_hat) if r > max_r_hat)
        messages.append(f"the independent ensembles disagree (split R-hat above {max_r_hat}): {bad}.")
    flat = chain.reshape(-1, ndim)
    span = problem.upper - problem.lower
    for i, name in enumerate(problem.labels):
        lo_mass = np.mean(flat[:, i] < problem.lower[i] + PILE_UP_EDGE * span[i])
        hi_mass = np.mean(flat[:, i] > problem.upper[i] - PILE_UP_EDGE * span[i])
        for mass, side in ((lo_mass, "lower"), (hi_mass, "upper")):
            if mass > PILE_UP_FRACTION:
                messages.append(f"{name} piles up against the {side} bound of its prior ({100 * mass:.0f}% of the posterior within "
                                f"{100 * PILE_UP_EDGE:.0f}% of the range): the data want a value outside the allowed range, or cannot constrain it.")
    for msg in messages:
        warnings.warn("fit: " + msg, stacklevel=2)
    return FitResult(problem.labels, chain, np.concatenate(lps, axis=1), tau, float(np.mean(acc)), r_hat, n_steps, n_discard, problem, messages)


def _starting_points(problem, n, seed):
    """n draws of the prior at which the posterior is finite."""
    rng_seed, out, tries = seed, [], 0
    while len(out) < n:
        tries += 1
        if tries > 200:
            raise ValueError("could not find starting points with a finite posterior in the prior range: the model does not reach the data anywhere.")
        for v in problem.sample_prior(n, seed=rng_seed + tries):
            if np.isfinite(problem.log_prob(v)):
                out.append(v)
                if len(out) == n:
                    break
    return np.array(out)


def problem_for_model(model, sample, data, y="bh_mass", x="stars_mass", bins=None, x_scale=1.0, reference_theta=None, n_boot=100, fit_s_int=True,
                      **relation_kwargs):
    """
    A Problem whose parameters are the free knobs of `model` and whose prediction is `model.relation(sample, theta, y, x)`.
    The errors of the model medians come from the tree sample's bootstrap at `reference_theta` (default: the preset's values), once, and are held fixed
    (interpolated at each bin) while sampling, so the likelihood needs one model run and no bootstrap per evaluation (docs/FIT_AUDIT.md, N3).
    """
    if not model.names:
        raise ValueError("the model has no free knobs.")
    kw = dict(y=y, x=x, x_scale=x_scale, **relation_kwargs)
    if bins is not None:
        kw["bins"] = bins
    ref = Curve.from_relation(model.relation(sample, reference_theta, n_boot=n_boot, **kw))
    if len(ref.x) < 2:
        raise ValueError("the reference run has fewer than two unflagged bins; widen the halo range of the tree sample or the bins.")

    def predict(theta):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)                # the reference run has already warned about flagged bins
            rel = model.relation(sample, theta, n_boot=0, **kw)
        err = np.interp(rel.x_mean, ref.x, ref.median_err)
        return Curve.from_relation(rel, median_err=err)

    return Problem(list(model.names), model.bounds[:, 0], model.bounds[:, 1], model.log, predict, data, fit_s_int)
