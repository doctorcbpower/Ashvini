"""
Smooth mean-accretion-history reference of Correa et al. (2015a, MNRAS 450, 1514), Section 3 summary equations:

    M(z) = M0 (1+z)^alpha exp(beta z)
    beta = -f(M0)
    alpha = [1.686 (2/pi)^0.5 dD/dz|_{z=0} + 1] f(M0)
    f(M0) = [S(M0/q) - S(M0)]^{-1/2}
    q = 4.137 z_tilde_f^{-0.9476}
    z_tilde_f = -0.0064 (log10 M0)^2 + 0.0237 log10 M0 + 1.8837

S(M) = sigma(M)^2 is the linear-theory variance at z = 0 and D(z) the linear growth factor normalised to unity at z = 0; both come
from foraois.cosmo_utils.CosmoData, so no new cosmology code is needed. M0 is the mass at z = 0.

Two points that matter when using this to benchmark halo histories selected at z = 5:
  * CosmoData here is built at the anchor redshift z = 5, so its sigma(M) is sigma(M, z=5) = D(5) sigma(M, 0). Correa's S needs the
    z = 0 value, so sigma is divided by D(5). (Without this alpha and beta come out far too large.)
  * M0 is a z = 0 mass. Matching a z = 5 mass means solving for the z = 0 mass whose mean history passes through it, which gives the
    history of the typical z = 0 halo with that z = 5 mass. A sample selected at z = 5 is biased towards early assembly, so this is a
    reference history, not the mean of the z = 5 sample, and above about 1e12 Msun at z = 5 no solution exists.

Paths: FORAOIS_ROOT (default ~/MyCodes/foraois) and CORREA_CONFIG (default <FORAOIS_ROOT>/config/planck2018_camb.yml).
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

FORAOIS_ROOT = os.environ.get("FORAOIS_ROOT", os.path.expanduser("~/MyCodes/foraois"))
FORAOIS_CONFIG = os.environ.get("CORREA_CONFIG", os.path.join(FORAOIS_ROOT, "config", "planck2018_camb.yml"))
sys.path.insert(0, os.path.join(FORAOIS_ROOT, "src"))

from foraois import cosmo_utils                       # noqa: E402
from foraois.utils import io as pymctrees_io          # noqa: E402

OUTDIR = Path(__file__).parent / "output"
Z_ANCHOR = 5.0

run_params = pymctrees_io.get_params(FORAOIS_CONFIG)
h = run_params["Cosmology"]["h"]
cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[Z_ANCHOR])
cosmo_data._prepare_sigma_grid(cosmo_data.get_power_spectrum())
D_ANCHOR = float(cosmo_data.get_linear_growth_and_collapse(redshift=Z_ANCHOR)[0])


def correa_alpha_beta(M0_msun, h=h):
    """alpha, beta, q, z_tilde_f, f(M0), dD/dz|0 for a z = 0 mass M0 [Msun]; sigma lookups are in Msun/h."""
    log10M0 = np.log10(M0_msun * h)
    z_tilde_f = -0.0064 * log10M0 ** 2 + 0.0237 * log10M0 + 1.8837
    q = 4.137 * z_tilde_f ** (-0.9476)
    S_M0 = (cosmo_data.sigma_at_logmass(log10M0) / D_ANCHOR) ** 2
    S_M0_over_q = (cosmo_data.sigma_at_logmass(log10M0 - np.log10(q)) / D_ANCHOR) ** 2
    f_M0 = 1.0 / np.sqrt(S_M0_over_q - S_M0)
    dz = 1.0e-3                                          # dD/dz at z = 0 by a one-sided difference (D(0) = 1)
    D_plus = cosmo_data.get_linear_growth_and_collapse(redshift=dz)[0]
    dD_dz0 = (D_plus - 1.0) / dz
    alpha = (1.686 * np.sqrt(2.0 / np.pi) * dD_dz0 + 1.0) * f_M0
    return float(alpha), float(-f_M0), float(q), float(z_tilde_f), float(f_M0), float(dD_dz0)


def correa_mass_history(M0_msun, z_array):
    alpha, beta, q, z_tilde_f, f_M0, dD_dz0 = correa_alpha_beta(M0_msun)
    M = M0_msun * (1.0 + np.asarray(z_array)) ** alpha * np.exp(beta * np.asarray(z_array))
    return M, dict(alpha=alpha, beta=beta, q=q, z_tilde_f=z_tilde_f, f_M0=f_M0, dD_dz0=dD_dz0)


def mass_at(M0_msun, z):
    a, b, *_ = correa_alpha_beta(M0_msun)
    return M0_msun * (1.0 + z) ** a * np.exp(b * z)


def shoot_correa_trajectory(M_target_msun, z_anchor, z_array, n_scan=160, max_decades=3.5):
    """z = 0 mass M0 for which M(z_anchor; M0) = M_target, found by scanning log M0 upwards from M_target for the first sign change
    and bisecting. Returns (nan history, nan M0) if M(z_anchor; M0) never reaches the target."""
    M0, prev = np.nan, None
    for trial in M_target_msun * np.logspace(0.3, max_decades, n_scan):
        try:
            m = mass_at(trial, z_anchor)
        except (OverflowError, ValueError, FloatingPointError):
            break
        if not np.isfinite(m):
            break
        if prev is not None and (prev[1] - M_target_msun) * (m - M_target_msun) <= 0:
            lo, hi = prev[0], trial
            for _ in range(60):
                mid = np.sqrt(lo * hi)
                if (mass_at(mid, z_anchor) - M_target_msun) * (mass_at(lo, z_anchor) - M_target_msun) <= 0:
                    hi = mid
                else:
                    lo = mid
            M0 = np.sqrt(lo * hi)
            break
        prev = (trial, m)
    if not np.isfinite(M0):
        return np.full(len(z_array), np.nan), dict(M0_correa_z0=float("nan"))
    M, params = correa_mass_history(M0, z_array)
    params["M0_correa_z0"] = M0
    return M, params


if __name__ == "__main__":
    OUTDIR.mkdir(exist_ok=True)
    z_array = np.linspace(Z_ANCHOR, 25.0, 2001)
    M_target = 5.0e10
    M, params = shoot_correa_trajectory(M_target, Z_ANCHOR, z_array)
    print("Correa (2015a) parameters for M(z=5) = %.3e Msun:" % M_target, params)
    print("M_correa(z=12) = %.4e Msun" % M[int(np.argmin(np.abs(z_array - 12.0)))])
    with open(OUTDIR / "correa_comparison_results.json", "w") as f:
        json.dump(dict(M_target=M_target, params=params, z=z_array.tolist(), M_correa=M.tolist()), f)
