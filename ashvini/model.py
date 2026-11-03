"""
`Model`: a physical set-up plus the knobs left free (docs/INTERFACE_DESIGN.md, step 1).

    model = Model(preset="energy_wind", free=["sf_efficiency", "wind_loading", "bh_accretion_eff"], fixed={"z_reion": 7.5})
    model.names                       # the free knobs, in order
    model.bounds                      # their ranges, shape (n, 2)
    cfg = model.config({"sf_efficiency": 0.03, "wind_loading": 0.2, "bh_accretion_eff": 1e-3})
    result = model.run(halo_mass, halo_mass_rate, redshift, theta=...)

Values are always physical (never logs). A knob that would change nothing under the chosen set-up (for example `wind_mass_slope` with the
energy wind, whose slope the escape velocity fixes, or a seed mass for a disabled seeding channel) is refused with a message that says why.
"""
import difflib

import numpy as np

from .config import DEFAULT, override
from .parameters import KNOBS, knob_range, _get

PRESETS = {
    "fiducial": {},
    "momentum_wind": {"sn.wind_model": "momentum"},
    "energy_wind": {"sn.wind_model": "energy"},
    "stellar_mass_wind": {"sn.wind_model": "stellar_mass"},
    "slim_disc": {"bh.growth_model": "hobbs_slimdisk"},
}


def _knob(name):
    if name not in KNOBS:
        close = difflib.get_close_matches(name, KNOBS, n=3)
        hint = f" Did you mean {', '.join(close)}?" if close else ""
        raise KeyError(f"{name!r} is not a knob.{hint} Knobs: {', '.join(KNOBS)} (see ashvini.parameters.parameters()).")
    return KNOBS[name]


class Model:
    """
    preset: a name in PRESETS, or a Config to start from.
    free: knob names to vary; fixed: {knob name or dotted Config path: value} held at a value other than the preset's.
    bounds: {knob name: (low, high)} to override the registry range.
    """

    def __init__(self, preset="fiducial", free=(), fixed=None, bounds=None):
        base = override(DEFAULT, **PRESETS[preset]) if isinstance(preset, str) and preset in PRESETS else preset
        if isinstance(preset, str) and preset not in PRESETS:
            raise KeyError(f"unknown preset {preset!r}; presets: {', '.join(PRESETS)}.")
        fixed = dict(fixed or {})
        base = override(base, **{k: v for k, v in fixed.items() if "." in k})   # dotted paths first: they can change which path a knob uses
        knob_fixed = {k: v for k, v in fixed.items() if "." not in k}
        free = list(free)
        both = set(free) & set(knob_fixed)
        if both:
            raise ValueError(f"{sorted(both)} are both free and fixed.")
        if len(set(free)) != len(free):
            raise ValueError("a knob is listed twice in free.")
        changes = {}
        for name, value in knob_fixed.items():
            changes[self._path(base, name)] = value
        self.base = override(base, **changes)
        self.names = free
        self._paths = {name: self._path(self.base, name) for name in free}
        bounds = dict(bounds or {})
        self.bounds = np.array([bounds.get(n, knob_range(n, self.base)[:2]) for n in free], dtype=float).reshape(-1, 2)
        self.log = np.array([KNOBS[n].log for n in free], dtype=bool)
        if np.any(self.bounds[:, 0] >= self.bounds[:, 1]) or np.any(self.log & (self.bounds[:, 0] <= 0)):
            raise ValueError("bounds must have low < high, and low > 0 for knobs sampled in log.")

    @staticmethod
    def _path(cfg, name):
        k = _knob(name)
        path = k.path(cfg)
        if path is None:
            raise ValueError(f"{name} does not apply: {k.why_inactive(cfg)}. Choose a preset, or fixed={{...}} with the Config path, that uses it.")
        if k.requires and not _get(cfg, k.requires):
            raise ValueError(f"{name} would change nothing: {k.requires} is False. Enable it with fixed={{'{k.requires}': True}}.")
        return path

    @property
    def theta0(self):
        """The preset's value of each free knob, in order."""
        return np.array([_get(self.base, self._paths[n]) for n in self.names], dtype=float)

    def _as_array(self, theta):
        if theta is None:
            return self.theta0
        if isinstance(theta, dict):
            missing = set(self.names) - set(theta)
            extra = set(theta) - set(self.names)
            if missing or extra:
                raise KeyError(f"theta must give exactly the free knobs {self.names}; missing {sorted(missing)}, unexpected {sorted(extra)}.")
            return np.array([theta[n] for n in self.names], dtype=float)
        theta = np.asarray(theta, dtype=float)
        if theta.shape != (len(self.names),):
            raise ValueError(f"theta has shape {theta.shape}; expected ({len(self.names)},) for {self.names}.")
        return theta

    def contains(self, theta):
        """True if every free knob lies within its bounds."""
        t = self._as_array(theta)
        return bool(np.all((t >= self.bounds[:, 0]) & (t <= self.bounds[:, 1])))

    def unit_to_theta(self, u):
        """Map points of the unit cube to physical values (log-spaced for the knobs that span decades). u: shape (..., n)."""
        u = np.asarray(u, dtype=float)
        lo, hi = self.bounds[:, 0], self.bounds[:, 1]
        return np.where(self.log, np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))), lo + u * (hi - lo))

    def config(self, theta=None):
        """The validated Config for the free knobs at `theta` (dict by name, or array in the order of `names`; None gives the preset's values)."""
        t = self._as_array(theta)
        return override(self.base, **{self._paths[n]: float(v) for n, v in zip(self.names, t)})

    def run(self, halo_mass, halo_mass_rate, redshift, theta=None, **kwargs):
        """run_forest on the given haloes with the Config at `theta`; returns its dictionary of arrays."""
        from .main import run_forest
        return run_forest(halo_mass, halo_mass_rate, redshift, config=self.config(theta), **kwargs)

    def describe(self):
        """One line per free knob: the Config path it sets, its value in the preset and its range."""
        rows = [f"{n}: {self._paths[n]}, preset value {v:g}, range {lo:g} to {hi:g} ({'log' if lg else 'linear'})"
                for n, v, (lo, hi), lg in zip(self.names, self.theta0, self.bounds, self.log)]
        return "\n".join(rows) if rows else "no free knobs"
