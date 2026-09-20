"""Plots for C18 test 4 (2x2: stellar wind x AGN) and the phase-space diagnostic d ln R/dt vs R. R = M_BH/M*_tot (M*_tot > 1e5 Msun)."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).parent
OUT = HERE / "output"
plt.style.use(str(HERE / "mnras_science.mplstyle")); plt.rcParams["axes.grid"] = False
D = json.load(open(OUT / "c18_test4_phase_results.json"))
LAB = {3e10: r"$3\times10^{10}$", 3e11: r"$3\times10^{11}$"}
arr = lambda x: np.array(x, dtype=float)

# ---- test 4
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
combos = [("wind1_f1", r"stellar wind on, AGN on", "C0"), ("wind1_f0", r"stellar wind on, AGN off", "C2"),
          ("wind0_f1", r"stellar wind off, AGN on", "C1"), ("wind0_f0", r"stellar wind off, AGN off", "C3")]
for ax, r in zip(axes, D):
    z = arr(r["z"])
    for key, lab, c in combos:
        a = r["t4"][f"{key}_1000"]; p = arr(a["pct"]); ok = np.isfinite(p[1])
        ax.fill_between(z[ok], p[0][ok], p[2][ok], color=c, alpha=0.13, lw=0)
        for row in arr(a["ind"]): ax.plot(z, row, color=c, lw=0.3, alpha=0.45)
        ax.plot(z[ok], p[1][ok], color=c, lw=1.6, label=lab)
        b = arr(r["t4"][f"{key}_1e+07"]["pct"] if f"{key}_1e+07" in r["t4"] else r["t4"][f"{key}_10000000"]["pct"]); okb = np.isfinite(b[1])
        ax.plot(z[okb], b[1][okb], color=c, lw=1.0, ls="--")
    ax.axhline(0.5, color="k", lw=0.6); ax.set_yscale("log"); ax.set_xlim(25, 5); ax.set_ylim(1e-4, 1e2)
    ax.set_xlabel(r"$z$"); ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]] + r" ($\sigma_j{=}1.5$, $R_{\rm nuc}{=}300$ pc)", fontsize=7.2)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}$"); axes[0].legend(fontsize=6, loc="lower left")
fig.suptitle(r"solid: seed $10^3$ (median, 16--84\%, 8 trees); dashed: seed $10^7$ (median)", fontsize=7, y=1.0)
fig.tight_layout(); fig.savefig(OUT / "sat_fig6_test4.png", dpi=250); plt.close(fig)

# ---- phase space
names = [("high_f1", r"high accessibility, AGN on, wind on"), ("high_f0", r"AGN off (wind on)"), ("high_nowind", r"stellar wind off (AGN on)"), ("fid_control", r"fiducial accessibility (control)")]
cols = {"z>10": "C3", "7<z<=10": "C2", "z<=7": "C0"}
fig, axes = plt.subplots(4, 2, figsize=(7.0, 8.6), sharex=True, sharey=True)
for j, r in enumerate(D):
    for i, (nm, title) in enumerate(names):
        ax = axes[i, j]; a = r["phase"][nm]
        lr, d, z = (arr(x) for x in a["pts"])
        for band, c in cols.items():
            lo, hi = {"z>10": (10, 99), "7<z<=10": (7, 10), "z<=7": (0, 7)}[band]
            s = (z > lo) & (z <= hi)
            ax.plot(10 ** lr[s], d[s], ls="none", marker=".", ms=1.2, color=c, alpha=0.10, rasterized=True)
            rows = np.array(a["stats"][band], dtype=float); ok = np.isfinite(rows[:, 2])
            ax.fill_between(10 ** rows[ok, 0], rows[ok, 3], rows[ok, 4], color=c, alpha=0.20, lw=0)
            ax.plot(10 ** rows[ok, 0], rows[ok, 2], color=c, lw=1.5, label=band)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1.0); ax.set_ylim(-40, 40); ax.set_xlim(10 ** -2.9, 10 ** 1.5)
        for R0 in (0.1, 0.2): ax.axvline(R0, color="0.5", lw=0.5, ls=":")
        ax.text(0.03, 0.05, title, transform=ax.transAxes, fontsize=6.2, va="bottom")
        if i == 0: ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]], fontsize=7.5)
        if i == 0 and j == 0: ax.legend(fontsize=5.8, loc="upper right", title="redshift band", title_fontsize=5.8)
for ax in axes[-1]: ax.set_xlabel(r"$R=M_{\rm BH}/M_{\star,{\rm tot}}$")
for ax in axes[:, 0]: ax.set_ylabel(r"$d\ln R/dt\ [{\rm Gyr}^{-1}]$")
fig.suptitle("lines: binned median, bands: 16--84\\%, dots: raw (tree, time) samples pooled over seeds; dotted verticals: $R=0.1,0.2$", fontsize=6.5, y=1.0)
fig.tight_layout(); fig.savefig(OUT / "sat_fig7_phase.png", dpi=250); plt.close(fig)
print("saved sat_fig6_test4.png, sat_fig7_phase.png")
