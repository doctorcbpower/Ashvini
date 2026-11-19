"""
`relation`: the distribution of y given x for a weighted tree sample (docs/RELATION_AUDIT.md, interface step 3).

    rel = model.relation(sample, theta, y="bh_mass", x="stars_mass", bins=np.arange(8.0, 11.5, 0.25))
    print(rel.table())          # per bin of x: effective number of trees, occupation fraction, median, 16th and 84th percentile, scatter, errors, flags
    rel.plot()

Bins are in log10 x. Within a bin, the statistics are weighted by the abundance of the haloes (so the mix of haloes at fixed x is right, including
the Eddington bias), and use log10 y of the trees with y > 0 only; the fraction with y > 0 is the occupation fraction. Errors come from a bootstrap over
blocks of neighbouring trees. Bins are flagged, never dropped: too few effective trees, a halo-mass range too narrow for the bin, and (for black hole
quantities) too much weight from haloes below the mass at which the black hole results are resolution-independent.
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

from .tree_sample import BH_RELIABLE_HALO_MASS, weighted_quantile

STARS_NOTE = ("stars_mass is all the mass ever formed into stars (no stellar mass loss; Salpeter IMF), while observed stellar masses are present-day "
              "(Chabrier or Kroupa). x_scale = {x_scale:g} applies no correction beyond that factor.")


@dataclass
class Relation:
    edges: np.ndarray                 # bin edges in log10 x (after x_scale)
    x_mean: np.ndarray                # weighted mean of log10 x of the trees in each bin
    ess: np.ndarray                   # Kish effective number of trees in each bin
    f_occ: np.ndarray                 # weighted fraction with y > 0
    f_occ_err: np.ndarray
    median: np.ndarray                # log10 y, occupied trees
    p16: np.ndarray
    p84: np.ndarray
    scatter: np.ndarray               # half of p84 - p16, dex
    median_err: np.ndarray
    scatter_err: np.ndarray
    low_stats: np.ndarray             # bool flags per bin
    range_limited: np.ndarray
    bh_unreliable: np.ndarray
    metadata: dict = field(default_factory=dict)

    @property
    def usable(self):
        """Bins with no flag: the ones a fit should use by default."""
        return ~(self.low_stats | self.range_limited | self.bh_unreliable)

    def table(self):
        head = f"{'log10 x':>8} {'n_eff':>7} {'f_occ':>11} {'median log y':>14} {'p16':>7} {'p84':>7} {'scatter':>13}  flags"
        rows = [head]
        for i in range(len(self.x_mean)):
            flags = ",".join(n for n, a in (("low_stats", self.low_stats), ("range_limited", self.range_limited), ("bh_unreliable", self.bh_unreliable)) if a[i])
            rows.append(f"{self.x_mean[i]:8.2f} {self.ess[i]:7.1f} {self.f_occ[i]:5.2f}+-{self.f_occ_err[i]:4.2f} "
                        f"{self.median[i]:7.2f}+-{self.median_err[i]:4.2f} {self.p16[i]:7.2f} {self.p84[i]:7.2f} {self.scatter[i]:6.2f}+-{self.scatter_err[i]:4.2f}  {flags}")
        notes = [self.metadata[k] for k in ("note",) if self.metadata.get(k)]
        return "\n".join(rows + notes)

    def plot(self, ax=None):
        """The median and the 16-84 band (open symbols where a bin is flagged), and the occupation fraction below it."""
        import matplotlib.pyplot as plt
        if ax is None:
            fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(5.5, 6), gridspec_kw={"height_ratios": [3, 1]})
        else:
            ax2 = None
        ok = self.usable
        ax.fill_between(self.x_mean, self.p16, self.p84, alpha=0.25, color="C0")
        ax.errorbar(self.x_mean[ok], self.median[ok], self.median_err[ok], fmt="o", color="C0", label="model median")
        ax.errorbar(self.x_mean[~ok], self.median[~ok], self.median_err[~ok], fmt="o", mfc="none", color="C0", label="flagged")
        ax.set_ylabel("log10 y (occupied)")
        ax.legend(frameon=False)
        if ax2 is not None:
            ax2.errorbar(self.x_mean, self.f_occ, self.f_occ_err, fmt="o-", color="C1")
            ax2.set_ylim(-0.05, 1.05)
            ax2.set_ylabel("fraction with y > 0")
            ax2.set_xlabel("log10 x")
        return ax


def _bin_statistics(bin_index, nbins, w, occupied, log_y, log_x, min_ess):
    """Per bin: (ess, mean log x, f_occ, median, p16, p84) for the given trees; NaN where there are too few effective trees."""
    out = np.full((6, nbins), np.nan)
    for b in range(nbins):
        sel = bin_index == b
        if not sel.any():
            continue
        wb = w[sel]
        ess = wb.sum() ** 2 / (wb ** 2).sum()
        out[0, b] = ess
        out[1, b] = np.average(log_x[sel], weights=wb)
        if ess < min_ess:
            continue
        occ = sel & occupied
        out[2, b] = w[occ].sum() / wb.sum()
        if occ.any():
            wo = w[occ]
            if wo.sum() ** 2 / (wo ** 2).sum() >= min_ess:
                out[3:, b] = weighted_quantile(log_y[occ], wo, [0.5, 0.16, 0.84])
    return out


def relation_from_arrays(x, y, weight, m0, bins, x_scale=1.0, n_boot=200, seed=0, min_ess=20.0, log_mh_range=None, edge_dex=0.3,
                         edge_fraction=0.02, bh=False, bh_fraction=0.10, block=10, note=""):
    """
    The relation of y given x for trees with final halo masses m0 (Msun) and abundance weights `weight`. x and y are one value per tree.
    bins: edges in log10(x * x_scale). log_mh_range: (low, high) of log10 m0 that was sampled, for the range check. bh: y is a black hole quantity.
    """
    x, y, w, m0 = (np.asarray(a, dtype=float) for a in (x, y, weight, m0))
    edges = np.asarray(bins, dtype=float)
    nb = len(edges) - 1
    has_x = x > 0
    log_x = np.log10(np.where(has_x, x * x_scale, 1.0))
    bin_index = np.where(has_x, np.digitize(log_x, edges) - 1, -1)
    bin_index = np.where((bin_index >= 0) & (bin_index < nb), bin_index, -1)
    occupied = y > 0
    log_y = np.log10(np.where(occupied, y, 1.0))
    keep = bin_index >= 0
    est = _bin_statistics(bin_index, nb, w, occupied, log_y, log_x, min_ess)

    order = np.argsort(m0)
    blocks = np.array_split(order, max(1, int(np.ceil(len(m0) / block))))
    rng = np.random.default_rng(seed)
    boot = np.full((n_boot, 6, nb), np.nan)
    for i in range(n_boot):
        idx = np.concatenate([blocks[j] for j in rng.integers(0, len(blocks), len(blocks))])
        boot[i] = _bin_statistics(bin_index[idx], nb, w[idx], occupied[idx], log_y[idx], log_x[idx], min_ess)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)       # an all-NaN bin has no error
        err = np.nanstd(boot, axis=0)
    err = np.where(np.isfinite(err), err, np.nan)

    log_m0 = np.log10(m0)
    range_limited = np.zeros(nb, dtype=bool)
    bh_unreliable = np.zeros(nb, dtype=bool)
    for b in range(nb):
        sel = bin_index == b
        if not sel.any():
            continue
        if log_mh_range is not None:
            near = (log_m0 < log_mh_range[0] + edge_dex) | (log_m0 > log_mh_range[1] - edge_dex)
            range_limited[b] = w[sel & near].sum() > edge_fraction * w[sel].sum()
        if bh:
            bh_unreliable[b] = w[sel & (m0 < BH_RELIABLE_HALO_MASS)].sum() > bh_fraction * w[sel].sum()
    low_stats = ~np.isfinite(est[2]) | ~np.isfinite(est[3])
    if bh and bh_unreliable.any():
        warnings.warn("black hole relation: bins with more than %.0f%% of their weight from haloes below %.0e Msun, where the seeding channel (and so the "
                      "black hole mass and occupation) depends on the tree resolution: log10 x = %s. See docs/TREESAMPLE_AUDIT.md, section 10."
                      % (100 * bh_fraction, BH_RELIABLE_HALO_MASS, ", ".join(f"{v:.2f}" for v in est[1][bh_unreliable])), stacklevel=3)
    meta = dict(x_scale=x_scale, n_boot=n_boot, seed=seed, min_ess=min_ess, log_mh_range=log_mh_range, edge_dex=edge_dex, edge_fraction=edge_fraction,
                bh=bh, bh_fraction=bh_fraction, block=block, note=note,
                weight_without_x=float(w[~has_x].sum() / w.sum()), weight_outside_bins=float(w[has_x & ~keep].sum() / w.sum()))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        scatter_err = np.nanstd(0.5 * (boot[:, 5] - boot[:, 4]), axis=0)
    scatter_err = np.where(np.isfinite(scatter_err), scatter_err, np.nan)
    return Relation(edges, est[1], est[0], est[2], err[2], est[3], est[4], est[5], 0.5 * (est[5] - est[4]), err[3], scatter_err,
                    low_stats, range_limited, bh_unreliable, meta)


def _quantity(result, spec):
    """One value per tree at the last step, from a name in the run's output or a function of the whole output."""
    value = spec(result) if callable(spec) else np.asarray(result[spec])
    return value[:, -1] if value.ndim == 2 else value


def relation(result, y, x="stars_mass", bins=np.arange(8.0, 11.5, 0.25), log_mh_range=None, x_scale=1.0, **kwargs):
    """
    The relation of `y` given `x` from the dictionary returned by TreeSample.run. x and y are names of its outputs ('stars_mass', 'bh_mass',
    'gas_mass', 'sfr', ...) or functions of the dictionary returning one value per tree (or an array over steps, of which the last is used).
    Other keywords go to relation_from_arrays.
    """
    note = STARS_NOTE.format(x_scale=x_scale) if x == "stars_mass" else ""
    bh = kwargs.pop("bh", y == "bh_mass")
    return relation_from_arrays(_quantity(result, x), _quantity(result, y), result["weight"], result["m0"], bins, x_scale=x_scale,
                                log_mh_range=log_mh_range, bh=bh, note=note, **kwargs)
