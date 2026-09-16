"""
Plots Figure 1 (critical seed-mass boundary, with seed-formation-channel
reference bands) from gen_fig1_critical_boundary.py's output JSON.
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

with open(OUTDIR / "fig1_critical_boundary_results.json") as f:
    data = json.load(f)
results = data["results"]

keys = sorted(results.keys(), key=lambda k: float(k))
m = np.array([float(k) for k in keys])
med = np.array([results[k]["median"] for k in keys])
p16 = np.array([results[k]["p16"] for k in keys])
p84 = np.array([results[k]["p84"] for k in keys])

fig, ax = plt.subplots(figsize=(6.0, 4.5))

bands = [
    ("Pop III remnants", 5e1, 5e2, "#8dd3c7"),
    ("Runaway/NSC collisions", 7e2, 6e3, "#fdb462"),
    ("Direct collapse", 8e4, 5e5, "#bc80bd"),
]
for label, lo, hi, color in bands:
    ax.axhspan(lo, hi, color=color, alpha=0.35, lw=0, zorder=0)
    ax.text(m.max() / 1.1, np.sqrt(lo * hi), label, fontsize=8.5, va="center", ha="right", color="black", zorder=5)

ax.fill_between(m, p16, p84, color="C0", alpha=0.25, label="16th--84th pct. (tree-to-tree)", zorder=3)
ax.plot(m, med, color="C0", lw=2.2, label=r"median $M_{\rm seed,crit}$", zorder=4)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$M_{\rm halo}(z=5)\ [M_\odot]$")
ax.set_ylabel(r"$M_{\rm seed,crit}\ [M_\odot]$")
ax.legend(loc="upper left", fontsize=8.5, frameon=False)
ax.set_xlim(m.min() * 0.9, m.max() * 1.1)
ax.set_ylim(3e1, 1.5e7)
fig.tight_layout()

outpath = OUTDIR / "fig_critical_boundary.png"
fig.savefig(outpath, dpi=200)
print("saved", outpath)

for target in [3e10, 5.335e10, 1e12, 1e13, 3e13]:
    idx = np.argmin(np.abs(m - target))
    print(f"M0~{target:.2e}: closest grid mass {m[idx]:.4e}, median={med[idx]:.4e}")
