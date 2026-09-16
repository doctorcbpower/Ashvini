"""
Plots Figure 3 (Appendix): representative trajectories, one panel per
halo mass, with S/E/G growth-regime background shading and a cold/hot-mode
crossing marker where relevant. Reads gen_fig3_trajectories.py's output
JSON. Styled with mnras_science.mplstyle (see that file's header for why).

Vertical banding note: the background shading is not broad blocks because
the nuclear gas reservoir exhibits a known, disclosed period-2 self-
regulation oscillation (paper's own "Numerical verification" appendix) --
whenever the nuclear free-fall depletion timescale is shorter than the
fixed integration step, M_gas oscillates step-to-step rather than
declining smoothly. Because the S/E/G regime boundary depends on the
instantaneous nuclear supply rate A_bh, which swings with that
oscillation, the classification flips almost every step near the
oscillating region instead of settling into a few wide blocks. Confirmed
in the paper not to affect M_seed,crit itself (converges to a few percent
under a 10x step-count refinement).
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
OUTDIR = HERE / "output"

plt.style.use(str(HERE / "mnras_science.mplstyle"))
plt.rcParams["axes.grid"] = False

with open(OUTDIR / "fig3_trajectories_results.json") as f:
    data = json.load(f)

M_hot = data["M_hot"]
masses = ["30000000000.0", "300000000000.0", "30000000000000.0"]
titles = [
    r"$M_{\rm halo}(z=5)=3\times10^{10}\,M_\odot$",
    r"$M_{\rm halo}(z=5)=3\times10^{11}\,M_\odot$",
    r"$M_{\rm halo}(z=5)=3\times10^{13}\,M_\odot$",
]

regime_color = {"S": "#fde0dd", "E": "#deebf7", "G": "#e5f5e0"}

fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6), sharey=False)

for ax, mkey, title in zip(axes, masses, titles):
    d = data["masses"][mkey]
    z = np.array(d["redshift"])
    mask = z <= 11.0
    z = z[mask]
    halo = np.array(d["halo_mass"])[mask]
    gas = np.array(d["gas_mass"])[mask]
    stars = np.array(d["stars_mass"])[mask]
    bh = np.array(d["bh_mass"])[mask]
    regime = np.array(d["regime"])[mask]

    # shade contiguous regime blocks in the background
    change = np.where(regime[1:] != regime[:-1])[0] + 1
    edges = np.concatenate(([0], change, [len(regime) - 1]))
    for i0, i1 in zip(edges[:-1], edges[1:]):
        ax.axvspan(z[i0], z[i1], color=regime_color[regime[i0]], lw=0, zorder=0)

    ax.plot(z, halo, color="0.2", lw=1.3, label=r"$M_{\rm halo}$")
    ax.plot(z, gas, color="C0", lw=1.2, label=r"$M_{\rm gas}$")
    ax.plot(z, stars, color="C2", lw=1.2, label=r"$M_\star$")
    ax.plot(z, bh, color="C3", lw=1.4, label=r"$M_{\rm BH}$")

    if float(d["M0"]) > M_hot:
        # mark cold/hot-mode crossing: first redshift (descending) where
        # halo mass exceeds M_hot
        above = halo >= M_hot
        if above.any():
            z_cross = z[above][0]
            ax.axvline(z_cross, color="0.35", lw=0.7, ls="--", zorder=1)

    positive = np.concatenate([gas[gas > 0], stars[stars > 0], bh[bh > 0]])
    ymin = positive.min() * 0.3
    ymax = max(halo.max(), bh.max()) * 3.0

    ax.set_yscale("log")
    ax.set_xlim(11, 5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$z$")
    ax.set_title(title, fontsize=8.0)
    ax.set_ylabel(r"Mass $[M_\odot]$")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.12),
           ncol=4, fontsize=7.5, frameon=False)

fig.tight_layout()
outpath = OUTDIR / "fig_representative_trajectories.png"
fig.savefig(outpath, dpi=300)
print("saved", outpath)
