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
