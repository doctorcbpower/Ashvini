"""
`TreeSample`: merger trees that form a population (docs/TREESAMPLE_AUDIT.md, interface step 2).

    sample = TreeSample.for_stellar_masses((8.0, 10.5), n_trees=1000, seed=1, cache_dir="cache")
    result = sample.run(model, theta)       # run_forest on every tree, plus sample.weight

Each tree has its own final halo mass, drawn evenly in log10 M so that every part of the range is sampled. Its weight is the abundance of haloes
of that mass (the Sheth-Tormen mass function of the foraois cosmology), so a weighted statistic over the sample is a statistic over the halo
population, and selecting on the model's own stellar mass afterwards (relation step) gives the right mix of haloes at fixed M*. Trees are built
with the resolution M_res = 1e-4 * M0, the finest the one-split-per-step condition allows, in batches by mass.

Scope (v1): z = 0, central haloes only (no environment, no subhaloes). Masses are the foraois (top-hat, Lagrangian) halo masses. A black hole
relation is reliable only for halo masses of about 1e11 Msun and above (BH_RELIABLE_HALO_MASS): below that the seeding channel that fires depends
on the tree resolution.
"""
import hashlib
import json
import os

import numpy as np

BH_RELIABLE_HALO_MASS = 1.0e11   # Msun; docs/TREESAMPLE_AUDIT.md, section 10
M_RES_FRACTION = 1.0e-4
Z_MAX, DZ, PK_KMAX = 30.0, 0.005, 3000.0   # the settings of the stored trees and of the resolution checks
DEFAULT_CONFIG = os.path.join(os.environ.get("FORAOIS_CONFIG_DIR", os.path.expanduser("~/MyCodes/foraois/config")), "planck2018_camb.yml")


def draw_log_masses(log_lo, log_hi, n, seed):
    """n values of log10 M in [log_lo, log_hi], one in each of n equal intervals (stratified, ascending), reproducible from `seed`."""
    rng = np.random.default_rng(seed)
    return log_lo + (np.arange(n) + rng.random(n)) / n * (log_hi - log_lo)


def abundance_weights(log_masses, dndlnm, log_lo, log_hi):
    """
    Weight of each tree: the number density of haloes it stands for, dn/dlnM(M) * (ln M_hi - ln M_lo) / n. The weights sum to the number density of
    haloes with masses in the range (the units of `dndlnm`). `dndlnm` is a function of M in Msun.
    """
    log_masses = np.asarray(log_masses, dtype=float)
    return dndlnm(10.0 ** log_masses) * (log_hi - log_lo) * np.log(10.0) / len(log_masses)


def weighted_quantile(x, weight, q):
    """Quantiles q (in 0..1) of x with the given weights (the inverse of the weighted empirical distribution, linearly interpolated)."""
    x, weight = np.asarray(x, dtype=float), np.asarray(weight, dtype=float)
    order = np.argsort(x)
    x, weight = x[order], weight[order]
    cum = (np.cumsum(weight) - 0.5 * weight) / weight.sum()
    return np.interp(q, cum, x)


def foraois_mass_function(config_path=DEFAULT_CONFIG, kind="sheth_tormen", pk_kmax=PK_KMAX):
    """
    dn/dlnM at z = 0 for the cosmology of a foraois config, as a function of M in Msun (comoving number density per ln M, in (Mpc/h)^-3).
    Tabulated by foraois over 1e8 to 1e15 Msun/h; asking outside that range raises.
    """
    from .pymctrees_adapter import _import_pymctrees
    cosmo_utils, _, io = _import_pymctrees()
    from foraois.mass_function_utils import MassFunctions
    params = io.get_params(config_path)
    params["Code"]["pk_kmax"] = float(pk_kmax)
    h = params["Cosmology"]["h"]
    logm, dn = MassFunctions(cosmo_utils.CosmoData(params, redshift=[0.0]), redshift=[0.0]).get_massfunc(0.0, kind)
    logm, log_dn = np.asarray(logm, dtype=float), np.log(np.asarray(dn, dtype=float).reshape(-1)[:len(logm)])

    def dndlnm(mass_msun):
        lm = np.log10(np.asarray(mass_msun, dtype=float) * h)    # foraois masses are Msun/h
        if lm.min() < logm[0] or lm.max() > logm[-1]:
            raise ValueError(f"halo masses must lie within 10^{logm[0]:g} to 10^{logm[-1]:g} Msun/h for the foraois mass function; "
                             f"got 10^{lm.min():.2f} to 10^{lm.max():.2f}.")
        return np.exp(np.interp(lm, logm, log_dn))

    return dndlnm


def moster13_mstar(halo_mass):
    """Moster et al. (2013) z = 0 stellar mass of a halo (Msun), used here only to decide how wide a halo-mass range to build, never in a result."""
    x = np.asarray(halo_mass, dtype=float) / 10 ** 11.590
    return np.asarray(halo_mass, dtype=float) * 2 * 0.0351 / (x ** -1.376 + x ** 0.608)


def halo_range_for_stellar_masses(log_mstar, margin_dex=1.0, grid=(7.0, 14.0)):
    """
    (log10 M_lo, log10 M_hi) of halo mass wide enough to contain every halo whose galaxy could have log10 M* in `log_mstar`:
    the Moster et al. (2013) relation inverted at the ends of the range, widened by `margin_dex` in halo mass on each side to cover the
    scatter and the disagreement between published relations (about 0.35 dex) and the model's own offset. It only sets the extent of the sample;
    whether it was wide enough is checked on the model's M* when a relation is made.
    """
    lm = np.linspace(*grid, 2801)
    ls = np.log10(moster13_mstar(10 ** lm))
    lo = np.interp(log_mstar[0], ls, lm) - margin_dex
    hi = np.interp(log_mstar[1], ls, lm) + margin_dex
    return float(max(lo, 8.2)), float(min(hi, 14.7))     # stay inside the mass function's table (1e8 to 10^14.9 Msun/h, about 1.5e8 to 1e15 Msun)


class TreeSample:
    """
    m0: final halo mass of each tree (Msun); weight: the abundance each stands for (see abundance_weights); halo_masses and growth_rates:
    the histories, shape (n, steps), earliest step first (as run_forest expects); redshifts: shape (steps,); settings: what built them.
    """

    def __init__(self, m0, weight, halo_masses=None, growth_rates=None, redshifts=None, settings=None):
        self.m0 = np.asarray(m0, dtype=float)
        self.weight = np.asarray(weight, dtype=float)
        self.halo_masses, self.growth_rates, self.redshifts = halo_masses, growth_rates, redshifts
        self.settings = dict(settings or {})
        if self.m0.shape != self.weight.shape:
            raise ValueError("m0 and weight must have the same shape.")

    @property
    def n(self):
        return len(self.m0)

    @property
    def number_density(self):
        """Sum of the weights: the number density of haloes in the sampled mass range (units of the mass function: (Mpc/h)^-3)."""
        return float(self.weight.sum())

    @property
    def effective_n(self):
        """Kish effective number of trees, (sum w)^2 / sum w^2: how many equally weighted trees the weighted sample is worth."""
        return float(self.weight.sum() ** 2 / (self.weight ** 2).sum())

    @property
    def bh_reliable(self):
        """True for the trees whose final halo mass is above BH_RELIABLE_HALO_MASS."""
        return self.m0 >= BH_RELIABLE_HALO_MASS

    @classmethod
    def build(cls, log_mh, n_trees=1000, seed=1, z=0.0, cache_dir=None, config_path=DEFAULT_CONFIG, mass_function="sheth_tormen",
              m_res_fraction=M_RES_FRACTION, z_max=Z_MAX, dz=DZ, pk_kmax=PK_KMAX, batch_dex=0.1):
        """
        Trees with final halo masses spread evenly in log10 M over `log_mh` = (low, high), weighted by the mass function. v1 builds z = 0 only.
        Memory is about 2 * n_trees * 6001 * 8 bytes (about 100 MB per thousand trees at the default dz).
        With `cache_dir`, an identical request (all settings, the config file's contents, the foraois version) is read back instead of rebuilt.
        """
        if z != 0.0:
            raise ValueError("v1 builds trees ending at z = 0 only (docs/INTERFACE_DESIGN.md, decision 7).")
        lo, hi = float(log_mh[0]), float(log_mh[1])
        if not lo < hi:
            raise ValueError("log_mh must be (low, high) with low < high.")
        settings = dict(log_mh=[lo, hi], n_trees=int(n_trees), seed=int(seed), z=z, mass_function=mass_function, m_res_fraction=m_res_fraction,
                        z_max=z_max, dz=dz, pk_kmax=pk_kmax, batch_dex=batch_dex, config=_file_hash(config_path), foraois=_foraois_version())
        path = None
        if cache_dir is not None:
            path = os.path.join(cache_dir, "trees_" + hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()[:16] + ".npz")
            if os.path.exists(path):
                return cls.load(path)
        from .pymctrees_adapter import build_forest_masses
        log_m = draw_log_masses(lo, hi, n_trees, seed)
        weight = abundance_weights(log_m, foraois_mass_function(config_path, mass_function, pk_kmax), lo, hi)
        hm, hg, zs, _ = build_forest_masses(config_path, 10.0 ** log_m, z, z_max, dz, m_res_fraction=m_res_fraction, batch_dex=batch_dex,
                                            seed=seed, pk_kmax_override=pk_kmax)
        sample = cls(10.0 ** log_m, weight, hm, hg, zs, settings)
        if path is not None:
            sample.save(path)
        return sample

    @classmethod
    def for_stellar_masses(cls, log_mstar, margin_dex=1.0, **kwargs):
        """A sample wide enough to hold the haloes of galaxies with log10 M* in `log_mstar` = (low, high); see halo_range_for_stellar_masses."""
        return cls.build(halo_range_for_stellar_masses(log_mstar, margin_dex), **kwargs)

    def run(self, model, theta=None, **kwargs):
        """model.run on every tree. Returns run_forest's dictionary of arrays with 'weight' and 'm0' added (one entry per tree)."""
        if self.halo_masses is None:
            raise ValueError("this TreeSample has no tree histories to run.")
        out = model.run(self.halo_masses, self.growth_rates, self.redshifts, theta=theta, **kwargs)
        out["weight"], out["m0"] = self.weight, self.m0
        return out

    def save(self, path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        tmp = path + ".tmp.npz"
        np.savez(tmp, m0=self.m0, weight=self.weight, halo_masses=self.halo_masses, growth_rates=self.growth_rates,
                 redshifts=self.redshifts, settings=json.dumps(self.settings))
        os.replace(tmp, path)

    @classmethod
    def load(cls, path):
        d = np.load(path, allow_pickle=False)
        return cls(d["m0"], d["weight"], d["halo_masses"], d["growth_rates"], d["redshifts"], json.loads(str(d["settings"])))


def _file_hash(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def _foraois_version():
    try:
        from importlib.metadata import version
        return version("foraois")
    except Exception:
        return "unknown"
