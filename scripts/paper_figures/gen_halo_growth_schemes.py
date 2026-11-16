"""
Main-progenitor growth of a halo of M0 = 3e12 Msun anchored at z = 5, in each scheme considered for the paper:

  PCH08 (compliant)      adaptive Parkinson, Cole & Helly 2008 sampler at a resolution where it keeps its splits-per-step target
  PCH08 (production)     the stored trees_M08.npz, run at m_res = 1e4 where 89 per cent of steps hit the step floor
  Zhang & Hui            fixed-step ZhangHuiMergerTree, dz = 0.05
  EPS conditional        median mass of the largest progenitor from the exact conditional mass function (Poisson), flat barrier
  EPS recursion          Correa et al. (2015a) Appendix A recursion with the progenitor-mass bound q = 2, 4, 6 (the q range is the uncertainty)
  Bera & Diemer (2026)   dlnM/dln a = Gamma(nu, n_eff, alpha_eff), the N-body calibrated fit, integrated back from z = 5
  Correa et al. (2015a)  closed form M0 (1+z)^alpha exp(beta z) for a z = 0 parent; it cannot reach 3e12 at z = 5 (its maximum is shown)

Usage: FORAOIS_ROOT=... CORREA_CONFIG=... python gen_halo_growth_schemes.py --trees-dir DIR --out growth_schemes.png
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
import gen_correa_comparison as C                       # noqa: E402
import gen_eps_main_progenitor_benchmark as E           # noqa: E402
import gen_bera_mar_benchmark as B                      # noqa: E402

M0, Z0, ZMAX = 3.0e12, 5.0, 20.0     # M0 is overwritten by --m0
ZG = np.arange(Z0, ZMAX + 1e-9, 0.25)


RES = None

def pch08(m_res, n_trees, seed):
    from ashvini.pch08_main_branch import build_main_branches
    hm, zz, st = build_main_branches(os.environ["CORREA_CONFIG"], M0, n_trees, Z0, 25.0, m_res, 0.05, seed)
    return np.asarray(hm), np.asarray(zz), int(np.asarray(st)[:, 2].sum())


def zhang_hui(n_trees, seed):
    from foraois import cosmo_utils, ZhangHuiMergerTree
    from foraois.utils import io
    from ashvini import pymctrees_adapter
    rp = io.get_params(os.environ["CORREA_CONFIG"])
    cdz = cosmo_utils.CosmoData(rp, redshift=[Z0])
    tg = ZhangHuiMergerTree(cdz, rp, model="cdm")
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, rp["Cosmology"]["h"], n_trees, z0=Z0, z_max=25.0,
                                                        m_res_msun=RES, dz=0.05, backend="numba", rng_seed=seed)
    return np.asarray(hm), np.asarray(zz)


def band(hm, zz):
    """hm chronological (z = 25 first). Returns percentiles 16/50/84 of the mass on ZG."""
    hm, zz = hm[:, ::-1], zz[::-1]                       # ascending z
    cols = [np.argmin(np.abs(zz - z)) for z in ZG]
    return np.percentile(hm[:, cols], [16, 50, 84], axis=0)


def bera_track():
    M, zc, out = M0, Z0, [M0]
    for zn in ZG[1:]:
        for sub in np.linspace(zc, zn, 6)[1:]:              # sub-steps of at most 0.05 in z
            G = B.bera_gamma(M, zc)[0]
            M *= ((1 + sub) / (1 + zc)) ** (-G)
            zc = sub
        out.append(M)
    return np.array(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--m0", type=float, default=3.0e12, help="halo mass at z = 5 [Msun]")
    ap.add_argument("--trees-dir", default=None, help="production trees (only for M0 >= 3e10)")
    ap.add_argument("--out", default="growth_schemes.png")
    ap.add_argument("--n-trees", type=int, default=60)
    a = ap.parse_args()
    M0 = a.m0
    RES = M0 / 3.0e5                                     # resolution ratio R = 3e5, the finest tested that keeps the step target

    curves = {}
    hm, zz, fl = pch08(RES, a.n_trees, 3)
    curves["PCH08, compliant (R=3e5, floor steps %d)" % fl] = band(hm, zz)
    idx = int(round(4 * np.log10(M0 / 3.0e10)))
    if a.trees_dir and abs(3.0e10 * 10 ** (0.25 * idx) / M0 - 1) < 1e-6 and 0 <= idx <= 12:
        t = np.load(Path(a.trees_dir) / f"trees_M{idx:02d}.npz")
        curves["PCH08, production trees (m_res %.0e, floor-hitting)" % float(t["m_res"])] = band(t["hm"], t["zz"])
    hz, zzz = zhang_hui(a.n_trees, 5)
    curves["Zhang & Hui (dz 0.05)"] = band(hz, zzz)
    eps = np.array([M0 if z <= Z0 else E.eps_median_mmp(M0, Z0, [z])[0] for z in ZG])
    rec = {q: E.eps_history(M0, Z0, list(ZG), q) for q in (2.0, 4.0, 6.0)}
    bera = bera_track()
    Mb, zb = None, None
    corr, cpar = C.shoot_correa_trajectory(M0, Z0, ZG)
    if np.isfinite(corr).all():
        best = cpar["M0_correa_z0"]
    else:                                               # target not reachable: use the parent that reaches the largest M(z = 5)
        best = max(10 ** np.linspace(12, 16, 60), key=lambda m: C.mass_at(m, Z0) if np.isfinite(C.mass_at(m, Z0)) else 0)
        corr = C.correa_mass_history(best, ZG)[0] if np.isfinite(C.mass_at(best, Z0)) else None

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    cols = {"PCH08, compliant": "#1f77b4", "PCH08, production": "#d62728", "Zhang": "#2ca02c"}
    palette = {"PCH08, compliant": "#1f77b4", "PCH08, production": "#d62728", "Zhang": "#2ca02c"}
    for name, p in curves.items():
        c = next(v for k, v in palette.items() if name.startswith(k))
        ax.fill_between(ZG, p[0], p[2], color=c, alpha=0.18, lw=0)
        ax.plot(ZG, p[1], color=c, lw=2, label=name)
    ax.plot(ZG, eps, "k-", lw=2, label="EPS conditional, median largest progenitor")
    ax.fill_between(ZG, rec[2.0], rec[6.0], color="0.5", alpha=0.25, lw=0, label="EPS recursion, q = 2 to 6 (q = 4 dashed)")
    ax.plot(ZG, rec[4.0], "--", color="0.35", lw=1.5)
    ax.plot(ZG, bera, color="#ff7f0e", lw=2.2, label="Bera & Diemer (2026) N-body fit, integrated")
    if corr is not None:
        ax.plot(ZG, corr, ":", color="#9467bd", lw=2, label="Correa+15 closed form (parent at z=0 chosen to reach M(z=5) closest to M0)")
    ax.set_yscale("log"); ax.set_xlim(5, ZMAX); ax.set_ylim(M0 * 1e-6, M0 * 2)
    ax.set_xlabel("redshift"); ax.set_ylabel(r"main-progenitor mass  [$M_\odot$]")
    ax.set_title(f"Growth of a {M0:.0e} $M_\\odot$ halo anchored at $z=5$")
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.tight_layout(); fig.savefig(a.out, dpi=170)

    print("median main-progenitor mass [Msun]")
    print(f"{'scheme':58s}" + "".join(f"z={z:<8g}" for z in (6, 7, 8, 10, 12, 15)))
    def row(name, f):
        print(f"{name:58s}" + "".join(f"{f(z):<10.2e}" for z in (6, 7, 8, 10, 12, 15)))
    for name, p in curves.items():
        row(name, lambda z, p=p: np.interp(z, ZG, p[1]))
    row("EPS conditional (median largest progenitor)", lambda z: np.interp(z, ZG, eps))
    for q in (2.0, 4.0, 6.0): row(f"EPS recursion q={q:g}", lambda z, q=q: np.interp(z, ZG, rec[q]))
    row("Bera & Diemer integrated", lambda z: np.interp(z, ZG, bera))
    if corr is not None:
        print(f"Correa closed form, M0(z=0)={best:.1e}: M(z=5)={C.mass_at(best, Z0):.2e}")
