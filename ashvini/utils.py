import numpy as np
import h5py
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u


Omega_m = cosmo.Om0
Omega_b = cosmo.Ob0
Omega_L = cosmo.Ode0


def mass_bin_group_name(mass_bin):
    """
    HDF5 group-naming convention for a mass bin, e.g. 1e10 -> "01e10",
    5e8 -> "05e08". Shared by read_trees() below and by
    ashvini.pymctrees_adapter (both the offline scripts/
    build_trees_from_pymctrees.py path and the live tree_source='pymctrees'
    path) -- kept in one place so they can't drift apart.
    """
    mass_bin = float(mass_bin)
    exponent = int(np.log10(mass_bin))
    mantissa = int(mass_bin / 10**exponent)
    return f"{mantissa:02d}e{exponent:02d}"


def read_trees(file_path, mass_bin):
    """
    Reads halo masses, halo growth rates, and redshifts for a mass bin.

    Parameters:
    - h5_file_path: str, path to the HDF5 file
    - mass_bin_float: float, e.g., 1e6, 5e7

    Returns:
    - m_halo: np.ndarray
    - halo_accretion_rate: np.ndarray
    - redshift: np.ndarray
    """
    group_name = mass_bin_group_name(mass_bin)

    with h5py.File(file_path, "r") as f:
        group = f[group_name]
        m_halo = group["halo_masses"][:]
        halo_accretion_rate = group["halo_growth_rates"][:]
        redshift = f["redshifts"][:]

    return m_halo, halo_accretion_rate, redshift


# --- Precompute time and redshift interpolation ---
# NOTE on interpolation scheme: these three lookups sit on the hot path of
# run1()'s per-timestep ODE right-hand sides (each solve_ivp call re-enters
# them at every internal solver substep). Profiling showed scipy's cubic
# interp1d objects dominate wall-clock time for a run (~50% of run1()'s
# cumulative time in a 200-step benchmark), mostly from per-call Python/
# object-dispatch overhead rather than the cubic evaluation itself. Switching
# to plain np.interp (linear, C-level, no object overhead) removes that
# dispatch cost. Linear interpolation error on this 10000-point z in [0, 50]
# grid is ~1e-5 in z (verified against the cubic result), many orders of
# magnitude below other modelling uncertainties in Ashvini, so this is a
# safe accuracy/speed trade. If higher-order accuracy is ever needed, revert
# to interp1d(..., kind="cubic") — the grid density can also be increased
# cheaply since these tables are built once at import time.
_z_vals = np.linspace(0, 50, 10000)
_t_vals = cosmo.age(_z_vals).value  # Gyr

_t_asc = _t_vals[::-1]  # ascending, required by np.interp
_z_asc = _z_vals[::-1]

_Hubble_time_vals = (1 / cosmo.H(_z_vals)).to(u.Gyr).value


def time_at_z(z):
    """Convert redshift to cosmic time (Gyr)."""
    return np.interp(z, _z_vals, _t_vals)


def z_at_time(t):
    """Convert cosmic time (Gyr) to redshift."""
    return np.interp(t, _t_asc, _z_asc)


def Hubble_time(z):
    """Return Hubble time (Gyr) at redshift z."""
    return np.interp(z, _z_vals, _Hubble_time_vals)


_h_of_z_vals = (cosmo.H(_z_vals) / (100 * u.km / u.s / u.Mpc)).to(u.dimensionless_unscaled).value


def h_of_z(z):
    """
    Dimensionless Hubble parameter h(z) = H(z) / (100 km/s/Mpc) -- e.g. for
    black_holes_growth.growth_ceiling's use of Power, Zubovas, Nayakshin &
    King (2011) eq. 21-22, which is expressed in these units (not H(z) in
    1/Gyr, which is what Hubble_time above gives the reciprocal of).
    """
    return np.interp(z, _z_vals, _h_of_z_vals)
