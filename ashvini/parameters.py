"""
The registry of physical knobs a user may leave free (docs/INTERFACE_DESIGN.md).

One entry per knob: its physical name, the Config path it sets (which can depend on the supernova wind model or the black hole growth
model), units, a range, and a one-sentence meaning. `parameters()` prints the table. A Model (model.py) builds a Config from the values
of the knobs it leaves free.

Ranges marked "exploration" are those used in the parameter-space exploration (docs/PARAMETER_SPACE_AUDIT.md). Ranges marked "provisional"
were chosen for this registry and have not been tested: treat them as a first guess, not a finding.
"""
from dataclasses import dataclass

from .config import DEFAULT


@dataclass(frozen=True)
class Knob:
    name: str
    paths: dict          # selector value -> Config path; selector None means one path for every case
    selector: str        # Config path whose value picks the path (e.g. "sn.wind_model"), or None
    units: str
    low: float
    high: float
    log: bool            # sample in log10 (the range spans decades)
    range_source: str    # "exploration" or "provisional"
    meaning: str
    solas_name: str = ""  # the Arepo Solas parameter this maps to, if there is a direct counterpart (empty: none identified)
    requires: str = ""    # Config path that must be True for the knob to do anything (e.g. a seeding channel's "enabled")
    only_when: tuple = None  # (Config path, value): the knob does something only when the path has this value

    def path(self, cfg):
        """The Config path this knob sets under `cfg`, or None if it does not apply (the wind model or growth model that uses it is not active)."""
        if self.only_when is not None and _get(cfg, self.only_when[0]) != self.only_when[1]:
            return None
        if self.selector is None:
            return self.paths[None]
        return self.paths.get(_get(cfg, self.selector))

    def why_inactive(self, cfg):
        """A phrase saying why the knob does not apply under `cfg`."""
        if self.only_when is not None and _get(cfg, self.only_when[0]) != self.only_when[1]:
            return f"it needs {self.only_when[0]} = {self.only_when[1]!r}, which is {_get(cfg, self.only_when[0])!r}"
        return f"{self.selector} is {_get(cfg, self.selector)!r}; defined for {', '.join(map(repr, self.paths))}"

    def default(self, cfg=DEFAULT):
        p = self.path(cfg)
        return None if p is None else _get(cfg, p)


def _get(cfg, path):
    obj = cfg
    for name in path.split("."):
        obj = getattr(obj, name)
    return obj


def _one(path):
    return {None: path}


KNOBS = {k.name: k for k in [
    Knob("sf_efficiency", _one("sf.efficiency"), None, "dimensionless (fraction of gas per free-fall time)", 0.003, 0.1, True, "exploration",
         "Star formation efficiency per halo free-fall time, in the halo free-fall star formation law."),
    Knob("wind_loading", {"ashvini": "sn.epsilon_p", "momentum": "sn.f_p", "energy": "sn.eps_fb", "stellar_mass": "sn.eta_norm"}, "sn.wind_model",
         "depends on the wind model: epsilon_p and f_p are dimensionless momentum factors; eps_fb is the fraction of the supernova energy that unbinds gas; "
         "eta_norm is the mass loading at M* = 1e10 Msun",
         0.01, 25.0, True, "provisional",
         "Amplitude of the supernova wind; which Config field it sets, and its range, depend on the wind model.",
         solas_name="wind/supernova feedback efficiency (the Solas parameter depends on the wind scheme; not yet mapped)"),
    Knob("wind_mass_slope", _one("sn.mass_slope"), None, "dimensionless", 0.0, 1.0, False, "provisional",
         "Exponent of the halo-mass dependence of the Ashvini wind loading, eta ~ M_halo^-slope (1/3 is the default; 2/3 energy-driven). "
         "Only the 'ashvini' loading has free exponents: the momentum and energy wind models fix them through the escape velocity.",
         only_when=("sn.wind_model", "ashvini")),
    Knob("wind_redshift_slope", _one("sn.redshift_slope"), None, "dimensionless", 0.0, 1.5, False, "provisional",
         "Exponent of the redshift dependence of the Ashvini wind loading, eta ~ (1+z)^-slope (1/2 default; 1 energy-driven).",
         only_when=("sn.wind_model", "ashvini")),
    Knob("z_reion", _one("reion.z_reion"), None, "redshift", 6.0, 9.0, False, "exploration",
         "Redshift at which the ultraviolet background starts to suppress gas accretion onto low-mass haloes."),
    Knob("bh_seed_mass_pop3", _one("bh.seeding.pop3.M_seed"), None, "Msun", 30.0, 1000.0, True, "provisional",
         "Mass of black holes seeded by the Pop III channel (near-primordial gas at z > z_min of the channel).",
         requires="bh.seeding.pop3.enabled"),
    Knob("bh_seed_mass_direct_collapse", _one("bh.seeding.direct_collapse.M_seed"), None, "Msun", 1e4, 1e6, True, "provisional",
         "Mass of black holes seeded by direct collapse (atomic-cooling haloes of pristine gas).",
         requires="bh.seeding.direct_collapse.enabled"),
    Knob("bh_seed_mass_threshold", _one("bh.seeding.halo_mass_threshold.M_seed"), None, "Msun", 1e2, 1e5, True, "provisional",
         "Mass of black holes seeded once a halo reaches the fallback halo-mass threshold.",
         requires="bh.seeding.halo_mass_threshold.enabled"),
    Knob("bh_accretion_eff", {"pznk11_freefall": "bh.efficiency", "hobbs_slimdisk": "bh.slimdisk.eta_acc"}, "bh.growth_model",
         "pznk11_freefall: fraction of gas per free-fall time; hobbs_slimdisk: nuclear accretion efficiency", 3e-4, 1e-2, True, "exploration",
         "How fast gas is accreted onto the black hole (not the radiative efficiency).",
         solas_name="BH accretion efficiency (Bondi-like prefactor; mapping not yet checked)"),
    Knob("bh_feedback_coupling", {"pznk11_freefall": "bh.eta_agn", "hobbs_slimdisk": "bh.slimdisk.epsilon_f"}, "bh.growth_model",
         "pznk11_freefall: wind mass loaded per unit mass accreted; hobbs_slimdisk: fraction of the accretion luminosity coupled to the gas",
         5e-4, 2.0, True, "provisional",
         "Strength of the AGN wind; which Config field it sets, and its range, depend on the black hole growth model.",
         solas_name="AGN feedback coupling efficiency (mapping not yet checked)"),
]}

# per-model ranges where the union above would be misleading
_PER_MODEL_RANGE = {
    ("wind_loading", "ashvini"): (1.0, 25.0, "exploration"),
    ("wind_loading", "momentum"): (0.1, 3.0, "provisional"),
    ("wind_loading", "energy"): (0.01, 1.0, "exploration"),
    ("wind_loading", "stellar_mass"): (1.0, 10.0, "provisional"),
    ("bh_feedback_coupling", "pznk11_freefall"): (0.05, 2.0, "exploration"),
    ("bh_feedback_coupling", "hobbs_slimdisk"): (5e-5, 5e-3, "provisional"),
    ("bh_accretion_eff", "hobbs_slimdisk"): (1e-3, 2e-2, "provisional"),
}


def knob_range(name, cfg=DEFAULT):
    """(low, high, range_source) of a knob under `cfg`, using the per-model range where one is defined."""
    k = KNOBS[name]
    sel = None if k.selector is None else _get(cfg, k.selector)
    return _PER_MODEL_RANGE.get((name, sel), (k.low, k.high, k.range_source))


def parameters(cfg=DEFAULT, show=True):
    """Print (and return) the table of knobs: name, Config path under `cfg`, default, range, units and meaning."""
    lines = []
    for name, k in KNOBS.items():
        path = k.path(cfg)
        lo, hi, src = knob_range(name, cfg)
        if path is None:
            lines.append(f"{name}\n    not used with the current set-up ({k.why_inactive(cfg)})\n")
            continue
        lines.append(f"{name}  [{path}]  default {k.default(cfg)!r}  range {lo:g} to {hi:g} ({'log' if k.log else 'linear'}, {src})\n"
                     f"    units: {k.units}\n    {k.meaning}\n" + (f"    Arepo Solas: {k.solas_name}\n" if k.solas_name else ""))
    text = "\n".join(lines)
    if show:
        print(text)
    return text
