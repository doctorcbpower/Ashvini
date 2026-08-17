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


def set_params(
    eps_sf=None,
    eps_fb=None,
    t_d=None,
    z_rei=None,
    sn_type=None,
    uv_background=None,
):
    """
    Override one or more fiducial parameters in place, for the remainder of
    the process (until the next set_params() call). Any argument left as
    None keeps its current value. Returns nothing; call print_current() to
    inspect the resulting state.
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


def print_current():
    print(
        f"eps_sf={_sf.e_ff} (main copy: {_main.e_ff}), eps_fb={_sn.epsilon_p}, "
        f"t_d={_main.t_d} Gyr, sn_type={_main.sn_type!r}, "
        f"UV_background={_main.UV_background}, z_rei={_reion.z_rei}"
    )
