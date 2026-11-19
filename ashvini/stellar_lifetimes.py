"""
Main-sequence lifetimes of the stars that explode as supernovae.

The lifetime is the time for the central hydrogen to fall to a mass fraction of 1e-4 in the rotating
(v_ini / v_crit = 0.4) stellar tracks of Ekstrom et al. (2012, A&A 537, A146; Z = 0.014) and Georgy et al.
(2013, A&A 558, A103; Z = 0.002). The two grids ship with the package as ashvini/data/ms_lifetimes_*.dat
(mass in Msun, lifetime in Myr); they were extracted from the CDS tables with Stochastic_Occupation's
`stellar_lifetimes._build_lifetime_grid`. Between tabulated masses the lifetime is interpolated linearly in
log mass and log lifetime; between the two metallicities, linearly in log Z; outside them the nearest grid is
used.
"""
from pathlib import Path

import numpy as np

Z_SOLAR_GRID = 0.014   # Ekstrom et al. 2012
Z_LOW_GRID = 0.002     # Georgy et al. 2013
MASS_MIN = 8.0         # Msun, the tabulated range that the supernova progenitors must lie in
MASS_MAX = 100.0

_DATA = Path(__file__).resolve().parent / "data"
_FILES = {Z_SOLAR_GRID: "ms_lifetimes_ekstrom2012_Z0.014.dat", Z_LOW_GRID: "ms_lifetimes_georgy2013_Z0.002.dat"}
_GRIDS = {}


def _grid(Z):
    """(log10 mass, log10 lifetime/Myr) for the grid metallicity Z."""
    if Z not in _GRIDS:
        mass, t_myr = np.loadtxt(_DATA / _FILES[Z], unpack=True)
        _GRIDS[Z] = (np.log10(mass), np.log10(t_myr))
    return _GRIDS[Z]


def t_ms_myr(mass_msun, Z=Z_LOW_GRID):
    """Main-sequence lifetime in Myr for initial mass(es) in Msun (8-100) and metallicity Z (mass fraction)."""
    scalar = np.ndim(mass_msun) == 0
    mass = np.asarray(mass_msun, dtype=float)
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("Stellar masses must be finite and positive.")
    if np.any((mass < MASS_MIN) | (mass > MASS_MAX)):
        raise ValueError(f"The lifetime grids cover {MASS_MIN:g} <= M/Msun <= {MASS_MAX:g}; got masses outside this range.")
    lm = np.log10(mass)
    if Z in _FILES:
        result = 10.0 ** np.interp(lm, *_grid(Z))
    else:
        w = np.clip((np.log10(Z) - np.log10(Z_LOW_GRID)) / (np.log10(Z_SOLAR_GRID) - np.log10(Z_LOW_GRID)), 0.0, 1.0)
        lo = np.interp(lm, *_grid(Z_LOW_GRID))
        hi = np.interp(lm, *_grid(Z_SOLAR_GRID))
        result = 10.0 ** ((1.0 - w) * lo + w * hi)
    return float(result) if scalar else result
