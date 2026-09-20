"""
Plots 1 to 5 of the saturation tests (C18). Reads output/c18_saturation_results.json. R = M_BH/M*_tot, masked where
M*_tot < 1e5 Msun (undefined before stars exist). Thin lines are individual trees (the same 8 trees in every panel);
thick lines are the median over 100 paired trees; bands are 16 to 84 per cent.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).parent
OUT = HERE / "output"
import paper_style
paper_style.apply()
D = json.load(open(OUT / "c18_saturation_results.json"))
LAB = {3e10: r"$3\times10^{10}$", 3e11: r"$3\times10^{11}$", 3e13: r"$3\times10^{13}$"}
SEEDS = [1e2, 1e3, 1e5, 1e7]
SEEDC = paper_style.seq(4)
SL = {1e2: r"$10^2$", 1e3: r"$10^3$", 1e5: r"$10^5$", 1e7: r"$10^7$"}
arr = lambda x: np.array(x, dtype=float)


def traj(ax, z, a, c, ls="-", lw=1.6, ind=True, band=True, label=None):
    p = arr(a["pct"]); ok = np.isfinite(p[1])
    if band: ax.fill_between(z[ok], p[0][ok], p[2][ok], color=c, alpha=0.15, lw=0)
    if ind:
        for row in arr(a["ind"]): ax.plot(z, row, color=c, lw=0.35, alpha=0.5)
    ax.plot(z[ok], p[1][ok], color=c, lw=lw, ls=ls, label=label)


# Plot 1: seed memory, trajectories
fig, axes = plt.subplots(1, 3, figsize=(paper_style.FULL, 2.8), sharey=True)
for ax, r in zip(axes, D):
    z = arr(r["z"])
    for k, s in enumerate(SEEDS):
        traj(ax, z, r["t1"][f"high_{s:g}"], SEEDC[k], label=SL[s])
    ax.axhline(0.5, color="k", lw=0.6); ax.set_yscale("log"); ax.set_xlim(25, 5); ax.set_ylim(1e-4, 1e2)
    ax.set_xlabel(r"$z$"); ax.set_title(r"$M_{\rm halo}(z{=}5)=$" + LAB[r["M0"]], fontsize=8)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}$"); axes[0].legend(title=r"seed [$M_\odot$]", fontsize=7, title_fontsize=7, loc="lower left")
axes[0].text(24.5, 0.6, r"$0.5$", fontsize=7)
fig.tight_layout(); paper_style.save(fig, OUT / "sat_fig1_trajectories"); plt.close(fig)

# Plot 2: final ratio vs seed mass
fig, axes = plt.subplots(1, 3, figsize=(paper_style.FULL, 2.8), sharey=True)
for ax, r in zip(axes, D):
    for cfg, ls, col in (("high", "-", paper_style.BLUE), ("fid", "--", "0.4")):
        F = np.array([r["t1"][f"{cfg}_{s:g}"]["final"] for s in SEEDS])
        for j in range(F.shape[1]): ax.plot(SEEDS, F[:, j], color=col, lw=0.2, alpha=0.25)
        q = np.percentile(F, [16, 50, 84], axis=1)
        ax.fill_between(SEEDS, q[0], q[2], color=col, alpha=0.2, lw=0)
        ax.plot(SEEDS, q[1], color=col, lw=1.6, ls=ls, marker="o", ms=2.5)
    ax.axhline(0.5, color="k", lw=0.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-9, 3)
    ax.set_xlabel(r"$M_{\rm seed}\ [M_\odot]$"); ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]], fontsize=8)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}\,(z=5)$")
axes[0].legend([Line2D([0], [0], color=paper_style.BLUE, lw=1.6), Line2D([0], [0], color="0.4", lw=1.6, ls="--")], [r"high accessibility", r"fiducial"], fontsize=7, loc="upper left")
fig.tight_layout(); paper_style.save(fig, OUT / "sat_fig2_final_vs_seed"); plt.close(fig)

# Plot 3: final ratio vs f_mom
fig, axes = plt.subplots(1, 3, figsize=(paper_style.FULL, 2.8), sharey=True)
FM = [0.0, 0.3, 1.0, 3.0, 10.0]
for ax, r in zip(axes, D):
    for s, c, ls in ((1e3, "C1", "-"), (1e7, "C3", "--")):
        F = np.array([r["t2"][f"f{f:g}_{s:g}"]["final"] for f in FM]); q = np.percentile(F, [16, 50, 84], axis=1)
        ax.fill_between(FM, q[0], q[2], color=c, alpha=0.18, lw=0); ax.plot(FM, q[1], color=c, ls=ls, lw=1.5, marker="o", ms=2.5, label=SL[s])
    ax.axhline(0.5, color="k", lw=0.6)
    ax.set_xscale("symlog", linthresh=0.1); ax.set_yscale("log"); ax.set_ylim(5e-4, 3)
    ax.set_xticks([0, 0.3, 1, 3, 10]); ax.set_xticklabels(["0", "0.3", "1", "3", "10"])
    ax.set_xlabel(r"$f_{\rm mom}$"); ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]], fontsize=8)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}\,(z=5)$"); axes[0].legend(title=r"seed [$M_\odot$]", fontsize=7, title_fontsize=7, loc="lower left")
fig.tight_layout(); paper_style.save(fig, OUT / "sat_fig3_fmom"); plt.close(fig)

# Plot 4: trajectories for f_mom = 0, 1, 10 (seed 1e3)
fig, axes = plt.subplots(1, 3, figsize=(paper_style.FULL, 2.8), sharey=True)
for ax, r in zip(axes, D):
    z = arr(r["z"])
    for f, c in ((0.0, "C2"), (1.0, "C0"), (10.0, "C3")):
        traj(ax, z, r["t2"][f"f{f:g}_1000"], c, label=rf"$f_{{\rm mom}}={f:g}$")
    ax.axhline(0.5, color="k", lw=0.6); ax.set_yscale("log"); ax.set_xlim(25, 5); ax.set_ylim(1e-4, 1e2)
    ax.set_xlabel(r"$z$"); ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]] + r", seed $10^3$", fontsize=8)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}$"); axes[0].legend(fontsize=7, loc="lower left")
fig.tight_layout(); paper_style.save(fig, OUT / "sat_fig4_fmom_trajectories"); plt.close(fig)

# Plot 5: all gas accessible against fiducial and high accessibility
fig, axes = plt.subplots(1, 3, figsize=(paper_style.FULL, 2.8), sharey=True)
for ax, r in zip(axes, D):
    z = arr(r["z"])
    for tag, key, c in (("fiducial", "fid", "0.4"), ("high (1.5, 300 pc)", "high", "C0"), ("all accessible, 100 pc", "all_R100", "C2"), ("all accessible, 300 pc", "all_R300", "C4")):
        for s, ls in ((1e3, "-"), (1e7, "--")):
            a = (r["t1"] if key in ("fid", "high") else r["t3"])[f"{key}_{s:g}"]
            traj(ax, z, a, c, ls=ls, ind=(s == 1e3), band=(s == 1e3), lw=1.5, label=tag if s == 1e3 else None)
    ax.axhline(0.5, color="k", lw=0.6); ax.set_yscale("log"); ax.set_xlim(25, 5); ax.set_ylim(1e-9, 1e2)
    ax.set_xlabel(r"$z$"); ax.set_title(r"$M_{\rm halo}=$" + LAB[r["M0"]], fontsize=8)
axes[0].set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}$"); axes[0].legend(fontsize=7, loc="lower left")
fig.tight_layout(); paper_style.save(fig, OUT / "sat_fig5_all_accessible"); plt.close(fig)
print("saved sat_fig1..5 in", OUT)
