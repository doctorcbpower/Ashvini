"""
Compare the dynamical mass accretion rate of PCH08 main branches with the N-body calibrated fit of Bera & Diemer (2026, arXiv:2606.09997).

Bera & Diemer measure the median Gamma_dyn = [log M(t) - log M(t - t_dyn)] / [log a(t) - log a(t - t_dyn)] of main-progenitor branches in
Erebos, IllustrisTNG and Thesan (z = 0-14) and fit it as a function of the peak height nu (their Eq. 1), the effective slope n_eff of the
power spectrum at the Lagrangian radius (Eq. 5) and the effective growth rate alpha_eff = -dlnD/dln(1+z) (Eq. 6):

    Gamma = A0 [ln(1 + nu)]^A1,    A_i = a_i + b_i / (3 + n_eff) + c_i alpha_eff                   (Eqs. 7, 8)

with (a0, a1, b0, b1, c0, c1) = (1.143, -0.002, 1.212, 1.830, -1.439, -0.506) (their Table 2). The fit is stated to reproduce the simulations to about
10 per cent for 0.5 < nu < 5 and 0 <= z <= 14, which spans the halos of the paper (nu = 1.5-5 at z = 5 for 3e10-3e13 Msun).
M is M200m and t_dyn = 2 R200m / V200m is the crossing time.

For a tree root at z0 (here z0 = 5, the epoch where the paper fixes the halo mass) the history before the root gives Gamma_dyn at z0 for halos
selected at fixed mass, hence at fixed nu, which is how Bera & Diemer bin their data, so the comparison is unbiased by selection.

Usage: FORAOIS_ROOT=... CORREA_CONFIG=<planck2018_camb.yml> python gen_bera_mar_benchmark.py --trees-dir DIR
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from astropy import constants as const, units as u
from astropy.cosmology import FlatLambdaCDM

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_correa_comparison as C                        # noqa: E402  (sigma(M) tables, growth factor, h)

A = dict(a0=1.143, a1=-0.002, b0=1.212, b1=1.830, c0=-1.439, c1=-0.506)
cd = C.cosmo_data
h = C.h
p = C.run_params["Cosmology"]
COSMO = FlatLambdaCDM(H0=100 * h, Om0=p["OmegaM"])      # same cosmology as the trees (foraois config)
DELTA_C = 1.686


def growth(z):
    return float(cd.get_linear_growth_and_collapse(redshift=float(z))[0])


def peak_height_neff(M_msun, z):
    """nu(M, z) and n_eff(M) with sigma from the table (sigma at z = 5 = D(5) sigma0, so nu uses the z-dependence explicitly)."""
    lg = np.log10(M_msun * h)
    s0 = cd.sigma_at_logmass(lg) / C.D_ANCHOR                      # sigma(M, z = 0)
    nu = DELTA_C / (s0 * growth(z))
    e = 0.02                                                       # dlnsigma/dlnM by central difference in log10 M
    dlns_dlnM = (np.log(cd.sigma_at_logmass(lg + e)) - np.log(cd.sigma_at_logmass(lg - e))) / (2 * e * np.log(10))
    neff = -2.0 * (dlns_dlnM * 3.0) - 3.0                          # R ~ M^(1/3): dlnsigma/dlnR = 3 dlnsigma/dlnM
    return nu, neff


def alpha_eff(z):
    eps = 1e-3
    return -(np.log(growth(z + eps)) - np.log(growth(max(z - eps, 0.0)))) / (np.log1p(z + eps) - np.log1p(max(z - eps, 0.0)))


def bera_gamma(M_msun, z):
    nu, neff = peak_height_neff(M_msun, z)
    ae = alpha_eff(z)
    A0 = A["a0"] + A["b0"] / (3 + neff) + A["c0"] * ae
    A1 = A["a1"] + A["b1"] / (3 + neff) + A["c1"] * ae
    return A0 * np.log(1 + nu) ** A1, nu, neff, ae


def t_dyn_gyr(M_msun, z):
    """Crossing time 2 R200m / V200m (independent of mass: 200 rho_m(z) sets R^3/M)."""
    rho_m = COSMO.Om0 * COSMO.critical_density0.cgs.value * (1 + z) ** 3
    M = M_msun * const.M_sun.cgs.value
    R = (3 * M / (4 * np.pi * 200 * rho_m)) ** (1 / 3)
    V = np.sqrt(const.G.cgs.value * M / R)
    return (2 * R / V) / (1e9 * 3.15576e7)


def gamma_dyn_from_tree(hm_row, zz, z0):
    """Gamma_dyn at z0 for one chronological-from-root main branch (zz ascending from z0). Returns nan if the earlier mass is unresolved."""
    M = hm_row[0]
    tdyn = t_dyn_gyr(M, z0)
    t0 = COSMO.age(z0).value
    # redshift one dynamical time earlier (a root-finder on the cosmic time grid)
    zg = np.linspace(z0, zz.max(), 4000)
    tg = COSMO.age(zg).value
    z_e = float(np.interp(t0 - tdyn, tg[::-1], zg[::-1]))
    if z_e > zz.max():
        return np.nan
    Me = float(np.exp(np.interp(z_e, zz, np.log(np.maximum(hm_row, 1e-30)))))
    if hm_row[np.searchsorted(zz, z_e)] <= 0:
        return np.nan
    return (np.log10(M) - np.log10(Me)) / (np.log10(1 / (1 + z0)) - np.log10(1 / (1 + z_e)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees-dir", required=True)
    args = ap.parse_args()
    print("halo mass at z=5   nu    n_eff alpha_eff | Bera&Diemer median Gamma | PCH08 median Gamma (16th-84th) | PCH08/Bera")
    for idx in range(13):
        M0 = 3.0e10 * 10 ** (0.25 * idx)
        t = np.load(Path(args.trees_dir) / f"trees_M{idx:02d}.npz")
        hm, zz = t["hm"], t["zz"]
        # trees are chronological in z (z = 25 ... 5 in the stored files or the reverse); order from the root outwards
        if zz[0] > zz[-1]:
            hm, zz = hm[:, ::-1], zz[::-1]
        g = np.array([gamma_dyn_from_tree(r, zz, zz[0]) for r in hm])
        gb, nu, neff, ae = bera_gamma(M0, zz[0])
        q = np.nanpercentile(g, [16, 50, 84])
        print(f"  {M0:9.2e}  {nu:5.2f} {neff:6.2f} {ae:6.3f} |  {gb:6.2f}  |  {q[1]:6.2f} ({q[0]:5.2f}-{q[2]:5.2f})  |  {q[1] / gb:5.2f}   [{int(np.isfinite(g).sum())} trees]")
