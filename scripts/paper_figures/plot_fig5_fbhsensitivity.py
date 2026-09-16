"""
Plots Figure 5 (Appendix): f_BH target sensitivity. Reads
gen_fig5_fbhsensitivity.py's output JSON. Styled with mnras_science.mplstyle.
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

with open(OUTDIR / "fig5_fbhsensitivity_results.json") as f:
    data = json.load(f)

results = data["results"]


def get_series(f_bh_key):
    d = results[f_bh_key]
    keys = sorted(d.keys(), key=lambda k: float(k))
    m = np.array([float(k) for k in keys])
    med = np.array([d[k]["median"] for k in keys])
    p16 = np.array([d[k]["p16"] for k in keys])
    p84 = np.array([d[k]["p84"] for k in keys])
    return m, med, p16, p84


fig, ax = plt.subplots(figsize=(3.5, 3.0))
for f_bh_key, color, label in [
    ("0.1", "C1", r"$f_{\rm BH}=0.1$"),
    ("0.5", "C0", r"$f_{\rm BH}=0.5$ (fiducial)"),
    ("0.9", "C3", r"$f_{\rm BH}=0.9$"),
]:
    m, med, p16, p84 = get_series(f_bh_key)
    ax.fill_between(m, p16, p84, color=color, alpha=0.15, lw=0)
    ax.plot(m, med, color=color, lw=1.4, label=label)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$M_{\rm halo}(z=5)\ [M_\odot]$")
ax.set_ylabel(r"$M_{\rm seed,crit}\ [M_\odot]$")
ax.set_title(r"$f_{\rm BH}$ target sensitivity", fontsize=8.5)
ax.legend(loc="lower right", fontsize=7.0, frameon=False)
fig.tight_layout()

outpath = OUTDIR / "fig_fbh_sensitivity.png"
fig.savefig(outpath, dpi=300)
print("saved", outpath)
