"""C18 test 5/6: R = M_BH/M*_tot(z=5) against eta_acc/eps_sf, with the phase-space zero crossing (z<=7)."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).parent; OUT = HERE / "output"
import paper_style
paper_style.apply()
D = json.load(open(OUT / "c18_eta_eps_results.json"))
LAB = {3e10: r"$3\times10^{10}$", 3e11: r"$3\times10^{11}$"}
CELL = {"nofb": "all feedback off", "fb": "AGN + stellar wind on"}
_EC = paper_style.seq(3)
ETA_C = {0.001: _EC[0], 0.005: _EC[1], 0.02: _EC[2]}
EPS_M = {0.0075: "^", 0.015: "o", 0.03: "s"}
fig, axes = plt.subplots(2, 2, figsize=(paper_style.FULL, 6.0), sharex=True, sharey=True)
x = np.logspace(-1.6, 0.6, 50)
for i, r in enumerate(D):
    for j, cell in enumerate(("nofb", "fb")):
        ax = axes[i, j]
        ax.plot(x, x, color="k", lw=0.8, ls="--"); ax.fill_between(x, 0.3 * x, x, color="0.85", lw=0, zorder=0)
        for row in [q for q in r["rows"] if q["cell"] == cell]:
            c, m = ETA_C[row["eta"]], EPS_M[row["eps"]]
            for sk, filled in (("1000", True), ("1e+07" if "1e+07" in row["seeds"] else "10000000", False)):
                R = np.array(row["seeds"][sk]["R"]); q = np.percentile(R, [16, 50, 84])
                ax.errorbar(row["ratio"], q[1], yerr=[[q[1] - q[0]], [q[2] - q[1]]], color=c, marker=m, ms=4, mfc=c if filled else "white", mec=c, lw=0.6, capsize=0)
            if row["cross_late"]:
                ax.plot(row["ratio"], 10 ** row["cross_late"][0], marker="x", color="k", ms=4, ls="none", mew=0.7)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.02, 4); ax.set_ylim(0.005, 5)
        ax.set_title(LAB[r["M0"]] + r" $M_\odot$: " + CELL[cell], fontsize=8)
        if i == 1: ax.set_xlabel(r"$\eta_{\rm acc}/\epsilon_{\rm sf}$")
        if j == 0: ax.set_ylabel(r"$M_{\rm BH}/M_{\star,{\rm tot}}\,(z=5)$")
h = [Line2D([0], [0], color=c, marker="o", ls="none", ms=4, label=rf"$\eta_{{\rm acc}}={e:g}$") for e, c in ETA_C.items()]
h += [Line2D([0], [0], color="0.3", marker=m, ls="none", ms=4, label=rf"$\epsilon_{{\rm sf}}={e:g}$") for e, m in EPS_M.items()]
h += [Line2D([0], [0], color="0.3", marker="o", ls="none", ms=4, label="seed $10^3$"), Line2D([0], [0], color="0.3", marker="o", mfc="white", ls="none", ms=4, label="seed $10^7$"),
      Line2D([0], [0], color="k", marker="x", ls="none", ms=4, label=r"phase-space zero crossing ($z\le7$)"), Line2D([0], [0], color="k", ls="--", lw=0.8, label=r"$R=\eta_{\rm acc}/\epsilon_{\rm sf}$")]
fig.tight_layout(rect=(0, 0.09, 1, 1))
fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7, frameon=False)
paper_style.save(fig, OUT / "sat_fig8_eta_eps")
print("saved sat_fig8_eta_eps.png")
