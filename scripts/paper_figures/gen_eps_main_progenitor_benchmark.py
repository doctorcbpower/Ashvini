"""
Mean main-progenitor mass history of a halo of mass M0 at an arbitrary anchor redshift z0, from the excursion-set recursion of
Correa et al. (2015a, MNRAS 450, 1514), their Appendix A:

    M_{j+1} = M_j [1 - erf( (omega_{j+1} - omega_j) / sqrt(2 [S(M_j/q) - S(M_j)]) )]        (A2)

stepped back in time from the parent halo (omega = delta_c(z) = 1.686 / D(z); S(M) = sigma^2(M) at z = 0). The parameter q bounds the
main progenitor from below, M_{j+1} >= M_j / q; in Correa et al. it is fixed by the formation-time argument of their Section 2 and
fitted for a parent halo at z = 0 (q = 4.137 z_f^-0.9476). There is no such calibration for a parent at z0 = 5, so q is treated as a
range and the spread is the uncertainty of the benchmark.

Two uses
  validate   parent at z = 0: the recursion with Correa's fitted q must reproduce their closed form M0 (1+z)^alpha exp(beta z)
  benchmark  parent at z = 5 (the anchor used for the paper's halo masses): mean main-progenitor mass at z = 7, 10, 15, 20, to be
             compared with the median and mean of the PCH08 main branches (trees_M??.npz)

Usage: FORAOIS_ROOT=... CORREA_CONFIG=... python gen_eps_main_progenitor_benchmark.py [--trees-dir DIR]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gen_correa_comparison as C                     # noqa: E402  (cosmology tables, sigma, growth factor)

cd = C.cosmo_data
h = C.h


def growth(z):
    return float(cd.get_linear_growth_and_collapse(redshift=float(z))[0])


def S_of_M(M_msun):
    """Variance sigma^2(M) at z = 0, M in Msun (the table is built at z = 5 and in Msun/h)."""
    return (cd.sigma_at_logmass(np.log10(M_msun * h)) / C.D_ANCHOR) ** 2


def eps_history(M0, z0, z_out, q, dz=0.02):
    """Mean main-progenitor mass at the redshifts z_out >= z0, by iterating (A2) with step dz in z."""
    z_out = np.asarray(sorted(z_out))
    zs = np.arange(z0, z_out.max() + dz / 2, dz)
    om = np.array([1.686 / growth(z) for z in zs])
    M = np.empty_like(zs)
    M[0] = M0
    for j in range(len(zs) - 1):
        Sj = S_of_M(M[j])
        Sq = S_of_M(M[j] / q)
        if Sq <= Sj:
            M[j + 1:] = M[j]
            break
        M[j + 1] = M[j] * erfc((om[j + 1] - om[j]) / np.sqrt(2.0 * (Sq - Sj)))
    return np.interp(z_out, zs, M)


def eps_median_mmp(M0, z0, z_out):
    """Median mass of the most massive progenitor from the exact EPS conditional mass function (flat barrier, Lacey & Cole 1993).

    The mean number of progenitors above M at redshift z of a parent M0 at z0 is
        N(>M) = integral_{S0}^{S(M)} (M0 / M1) f(S1, w1 | S0, w0) (dM1 / dS1)... evaluated on a grid in S,
    with f(S1, w1 | S0, w0) = dw / sqrt(2 pi) (S1 - S0)^(-3/2) exp(-dw^2 / (2 (S1 - S0))).
    Treating the progenitors as a Poisson population, the median of the largest one solves N(>M) = ln 2.
    It ignores correlations between progenitors, so it is a reference estimate, not a tree.
    """
    w0 = 1.686 / growth(z0)
    S0 = S_of_M(M0)
    logM = np.linspace(np.log10(M0) - 8.0, np.log10(M0) - 1e-4, 4000)         # progenitor masses below M0
    Sg = np.array([S_of_M(10 ** x) for x in logM])
    out = []
    for z in z_out:
        dw = 1.686 / growth(z) - w0
        dS = Sg - S0
        ok = dS > 0
        # number density per ln M: (M0/M1) f(S1) |dS1/dlnM1|
        dSdlnM = np.abs(np.gradient(Sg, logM * np.log(10.0)))
        n_lnM = np.where(ok, (M0 / 10 ** logM) * dw / np.sqrt(2 * np.pi) * np.where(ok, dS, 1.0) ** -1.5
                         * np.exp(-dw ** 2 / (2 * np.where(ok, dS, 1.0))) * dSdlnM, 0.0)
        # N(>M) is the integral from M upwards: cumulative sum from the high-mass end
        lnM = logM * np.log(10.0)
        Ngt = np.concatenate([np.cumsum(((n_lnM[1:] + n_lnM[:-1]) * 0.5 * np.diff(lnM))[::-1])[::-1], [0.0]])
        idx = np.where(Ngt >= np.log(2.0))[0]
        out.append(10 ** logM[idx.max()] if len(idx) else np.nan)
    return np.array(out)


def validate(M0=1.0e12):
    a, b, q, *_ = C.correa_alpha_beta(M0)
    zt = [2.0, 5.0, 8.0, 10.0]
    rec = eps_history(M0, 0.0, zt, q)
    ref = np.array([M0 * (1 + z) ** a * np.exp(b * z) for z in zt])
    print(f"validation, parent {M0:.0e} Msun at z=0, q={q:.3f}: recursion / Correa closed form at z={zt}: {np.round(rec / ref, 3)}")


def benchmark(trees_dir, q_grid=(2.0, 3.0, 4.0, 6.0)):
    zt = [7.0, 10.0, 15.0, 20.0]
    print("parent at z=5: mean EPS main progenitor [Msun] for q in", q_grid, "| PCH08 median (mean) main branch | ratio PCH08 median / EPS (q=3, range over q)")
    for idx in (0, 4, 8, 12):
        M0 = 3.0e10 * 10 ** (0.25 * idx)
        eps = np.array([eps_history(M0, 5.0, zt, q) for q in q_grid])             # (n_q, n_z)
        t = np.load(Path(trees_dir) / f"trees_M{idx:02d}.npz")
        hm, zz = t["hm"], t["zz"]
        for k, z in enumerate(zt):
            j = int(np.argmin(np.abs(zz - z)))
            med, mean = np.median(hm[:, j]), np.mean(hm[:, j])
            r = med / eps[:, k]
            mm = eps_median_mmp(M0, 5.0, [z])[0]
            print(f"  M(z=5)={M0:.1e} z={z:>4.0f}: recursion {np.array2string(eps[:, k], precision=2)} | exact-EPS median MMP {mm:.2e} | PCH08 {med:.2e} ({mean:.2e}) | PCH08/MMP {med / mm:.2f}, PCH08/recursion(q=3) {r[1]:.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees-dir", default=None)
    args = ap.parse_args()
    validate()
    if args.trees_dir:
        benchmark(args.trees_dir)
