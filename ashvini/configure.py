"""
Runtime parameter overrides for parameter-sweep notebooks/scripts.

star_formation.e_ff, supernovae_feedback.epsilon_p/pi_fid, reionization's
z_rei/gamma/omega (plus beta/c_omega, derived from them), and main.py's own
UV_background/t_d/sn_type/agn_delay_time/e_ff are all read from PARAMS once,
at import time, into module-level globals -- run_params.yaml is meant to be
edited and the process restarted, not reconfigured mid-session. A parameter
sweep (e.g. reproducing Menon & Power 2024's Figures 4-7, which vary
epsilon_sf/epsilon_fb/t_d/z_rei one at a time against a fixed fiducial) needs
to change these live. set_params() does that, patching every module-level
copy of a given quantity together -- notably main.py duplicates
star_formation.e_ff as its own module-level e_ff (used directly in
run_forest's inlined closed-form star-formation update, not via
star_formation.star_formation_rate), so patching only one of the two would
silently desync the fast (run_forest) and reference (run1_scalar) paths.
"""
import numpy as np

from . import main as _main
from . import star_formation as _sf
from . import supernovae_feedback as _sn
from . import reionization as _reion
from . import agn_feedback as _agn
from . import black_holes_growth as _bh
from .run_params import PARAMS


def set_params(
    eps_sf=None,
    eps_fb=None,
    t_d=None,
    z_rei=None,
    sn_type=None,
    uv_background=None,
    bh_enabled=None,
    eta_agn=None,
    sigma_feedback_enabled=None,
    growth_cap_enabled=None,
    e_bh=None,
    eddington_multiplier=None,
):
    """
    Override one or more fiducial parameters in place, for the remainder of
    the process (until the next set_params() call). Any argument left as
    None keeps its current value. Returns nothing; call print_current() to
    inspect the resulting state.

    bh_enabled/eta_agn/sigma_feedback_enabled/growth_cap_enabled/e_bh
    control black hole growth and AGN feedback (see MODELS.md's "Black hole
    growth"/"AGN feedback" sections). growth_cap_enabled is only
    physically meaningful with sigma_feedback_enabled=True too (see
    black_holes_growth.growth_ceiling) but isn't validated against it
    here -- set both explicitly if you want the hard-cap self-regulation
    variant. Unlike eps_sf etc., PARAMS.bh.seeding.*.enabled is read
    fresh from PARAMS on every seeding_mask_and_mass() call (not cached
    into a module-level global at import time -- black_holes_growth.py
    keeps a live reference to PARAMS.bh.seeding, it never copies its
    fields out), so bh_enabled can toggle it by mutating PARAMS directly;
    eta_agn and agn_feedback.sigma_feedback_enabled *are* cached at import
    time (same pattern as eps_sf/eps_fb elsewhere in this module) and are
    patched accordingly.

    bh_enabled=False sets every seeding channel's `enabled` to False (no
    new BH is ever seeded, so growth/feedback are moot regardless of
    eta_agn/eddington_multiplier -- this is the only way to get a genuine
    "no BH" baseline, since run_forest()/run1_scalar() always compute BH
    seeding/growth/AGN terms unconditionally). bh_enabled=True restores
    all three channels to enabled -- if you'd previously disabled only
    some channels by hand via PARAMS.bh.seeding directly, this will
    re-enable all of them, not just the ones you disabled.
    """
    if eps_sf is not None:
        _sf.e_ff = eps_sf
        _main.e_ff = eps_sf  # see module docstring: duplicated in main.py

    if eps_fb is not None:
        _sn.epsilon_p = eps_fb

    if t_d is not None:
        _main.t_d = t_d

    if sn_type is not None:
        _main.sn_type = sn_type

    if uv_background is not None:
        _main.UV_background = uv_background

    if z_rei is not None:
        _reion.z_rei = z_rei
        # beta and c_omega are derived from z_rei/gamma/omega once at
        # import time (reionization.py module level) -- re-derive them
        # here with the same formulas rather than leaving them stale.
        gamma, omega = _reion.gamma, _reion.omega
        _reion.beta = z_rei * (
            (np.log(1.82 * (10**3) * np.exp(-0.63 * z_rei) - 1)) ** (-1 / gamma)
        )
        _reion.c_omega = 2 ** (omega / 3) - 1

    if bh_enabled is not None:
        seeding = PARAMS.bh.seeding
        seeding.pop3.enabled = bh_enabled
        seeding.direct_collapse.enabled = bh_enabled
        seeding.halo_mass_threshold.enabled = bh_enabled

    if eta_agn is not None:
        _agn.eta_agn = eta_agn

    if sigma_feedback_enabled is not None:
        _agn.sigma_feedback_enabled = sigma_feedback_enabled

    if growth_cap_enabled is not None:
        _bh.growth_cap_enabled = growth_cap_enabled

    if e_bh is not None:
        _bh.e_bh = e_bh

    if eddington_multiplier is not None:
        _bh.eddington_multiplier = eddington_multiplier


def print_current():
    seeding = PARAMS.bh.seeding
    print(
        f"eps_sf={_sf.e_ff} (main copy: {_main.e_ff}), eps_fb={_sn.epsilon_p}, "
        f"t_d={_main.t_d} Gyr, sn_type={_main.sn_type!r}, "
        f"UV_background={_main.UV_background}, z_rei={_reion.z_rei}\n"
        f"BH seeding enabled: pop3={seeding.pop3.enabled}, "
        f"direct_collapse={seeding.direct_collapse.enabled}, "
        f"halo_mass_threshold={seeding.halo_mass_threshold.enabled}; "
        f"eta_agn={_agn.eta_agn}, sigma_feedback_enabled={_agn.sigma_feedback_enabled}, "
        f"growth_cap_enabled={_bh.growth_cap_enabled}, e_bh={_bh.e_bh}, "
        f"eddington_multiplier={_bh.eddington_multiplier}"
    )
