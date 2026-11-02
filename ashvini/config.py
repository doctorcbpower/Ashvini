"""
Configuration for the general galaxy-formation chain (milestone M1 of docs/SUPERSET_DESIGN.md).

`Config` is the parameter tree defined in run_params.py (sections io, sf, sn, reion, metals, dust, bh).
Physics functions take the parameters they need from a `Config` passed explicitly (`cfg=...`); a function
called without one uses `DEFAULT`, the configuration loaded from run_params.yaml. `DEFAULT` is the single
object that `configure.set_params` modifies, so code that changes parameters in place keeps working, and
a run that needs different parameters in the same process builds its own `Config` and passes it:

    cfg = override(DEFAULT, **{"sf.efficiency": 0.05, "sn.type": "instantaneous"})
    main.run_forest(halo_mass, halo_mass_rate, redshift, config=cfg)

Nothing is read from a module-level copy at import time, so there is nothing to keep in step.
"""
import copy

from .run_params import Params as Config, PARAMS, load_params

DEFAULT = PARAMS  # the same object as run_params.PARAMS

SN_TYPES = ("delayed", "instantaneous", "no")
GROWTH_MODELS = ("pznk11_freefall", "hobbs_slimdisk")
TREE_SOURCES = ("file", "pymctrees")
KERNELS = ("delta", "imf")
CLOSURES = ("ashvini", "momentum", "energy")


def resolve(cfg=None):
    """The configuration to use: `cfg` if given, otherwise the process-wide DEFAULT."""
    return DEFAULT if cfg is None else cfg


def load_config(config_file=None):
    """Load and validate a Config from a yaml file (default: run_params.yaml)."""
    cfg = load_params(config_file)
    validate(cfg)
    return cfg


def validate(cfg):
    """Raise ValueError with a message that says what to change; returns cfg."""
    if cfg.sn.type not in SN_TYPES:
        raise ValueError(f"supernova.type must be one of {SN_TYPES}; got {cfg.sn.type!r}.")
    if cfg.sn.delay_time < 0:
        raise ValueError(f"supernova.delay_time must be >= 0 Gyr; got {cfg.sn.delay_time}.")
    if cfg.sn.kernel not in KERNELS:
        raise ValueError(f"supernova.kernel must be one of {KERNELS}; got {cfg.sn.kernel!r}.")
    if cfg.sn.closure not in CLOSURES:
        raise ValueError(f"supernova.closure must be one of {CLOSURES}; got {cfg.sn.closure!r}.")
    if cfg.sn.kernel == "imf" and cfg.sn.m_sn_max > 100.0:
        raise ValueError("supernova.m_sn_max must be <= 100 Msun with kernel 'imf' (the lifetime grids stop there).")
    if cfg.events.bank_tau_myr is not None and cfg.sn.closure != "energy":
        raise ValueError("events.bank_tau_myr needs supernova.closure 'energy' (the bank stores unused explosion energy).")
    if cfg.events.eps_cut_mass is not None and cfg.sn.closure != "energy":
        raise ValueError("events.eps_cut_mass switches the energy closure's eps_fb, so it needs supernova.closure 'energy'.")
    if not 0.0 < cfg.events.access_fraction <= 1.0:
        raise ValueError(f"events.access_fraction must be in (0, 1]; got {cfg.events.access_fraction}.")
    if cfg.bh.growth_model not in GROWTH_MODELS:
        raise ValueError(f"black_holes.growth_model must be one of {GROWTH_MODELS}; got {cfg.bh.growth_model!r}.")
    if cfg.io.tree_source not in TREE_SOURCES:
        raise ValueError(f"basics.tree_source must be one of {TREE_SOURCES}; got {cfg.io.tree_source!r}.")
    if cfg.bh.sigma_feedback.growth_cap_enabled and not cfg.bh.sigma_feedback.enabled:
        raise ValueError("black_holes.sigma_feedback.growth_cap_enabled needs sigma_feedback.enabled: True.")
    from .imf import nu_SN
    nu_SN(cfg.sn.imf, cfg.sn.m_sn_min, cfg.sn.m_sn_max)  # raises on a bad IMF or progenitor range
    return cfg


def override(cfg=None, **changes):
    """
    A validated copy of `cfg` (default: DEFAULT) with dotted-path changes applied, leaving `cfg` untouched:

        override(DEFAULT, **{"sf.efficiency": 0.05, "bh.seeding.pop3.enabled": False})

    Paths follow the yaml structure as laid out in run_params.py (io, sf, sn, reion, metals, dust, bh).
    """
    new = copy.deepcopy(resolve(cfg))
    for path, value in changes.items():
        *parents, leaf = path.split(".")
        obj = new
        for name in parents:
            obj = getattr(obj, name)
        if not hasattr(obj, leaf):
            raise AttributeError(f"Config has no parameter {path!r}.")
        setattr(obj, leaf, value)
    return validate(new)


def retire_module_values(module_name, replacements):
    """
    Make assigning to a retired module-level value an error instead of a silent no-op.

    Parameters used to be copied into module-level variables at import (main.e_ff, agn_feedback.eta_agn, ...).
    They are now read from the Config at each call, so assigning to the old name would change nothing. Called
    from the module itself: `replacements` maps each retired name to the Config path that replaces it.
    """
    import sys
    import types

    module = sys.modules[module_name]

    class _Guarded(types.ModuleType):
        def __setattr__(self, name, value):
            if name in replacements:
                raise AttributeError(
                    f"{module_name}.{name} no longer exists: set {replacements[name]} instead "
                    "(configure.set_params(...), or config.override(...) and pass config=...)."
                )
            super().__setattr__(name, value)

    module.__class__ = _Guarded
