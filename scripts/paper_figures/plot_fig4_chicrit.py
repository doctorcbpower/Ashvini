"""
Plots Figure 4 (Appendix): chi_crit reconciliation. Reads
gen_fig4_chicrit.py's output JSON. Styled with mnras_science.mplstyle.
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

with open(OUTDIR / "fig4_chicrit_results.json") as f:
    data = json.load(f)

seed = np.array(data["seed_masses"])
chi_crit_values = [1.0, 3.0, 10.0, 1.0e2, 1.0e12]
labels = [r"$\chi_{\rm crit}=1$", r"$\chi_{\rm crit}=3$", r"$\chi_{\rm crit}=10$",
          r"$\chi_{\rm crit}=10^2$", r"$\chi_{\rm crit}=10^{12}$"]
colors = ["C3", "C1", "C0", "C2", "0.2"]
styles = ["-", "-", "-", "--", ":"]

fig, ax = plt.subplots(figsize=(3.5, 3.0))
for chi_crit, label, color, ls in zip(chi_crit_values, labels, colors, styles):
    curve = np.array(data["curves"][str(chi_crit)])
    ax.plot(seed, curve, color=color, ls=ls, lw=1.3, label=label)
    m_crit = data["boundaries"][str(chi_crit)]
    ax.axvline(m_crit, color=color, ls=":", lw=0.7, alpha=0.6)

ax.axhline(0.5, color="0.5", lw=0.6, ls="--", alpha=0.6)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$M_{\rm seed}\ [M_\odot]$")
ax.set_ylabel(r"$f_{\rm BH} = M_{\rm BH}/M_\star\ (z=5)$")
ax.set_title(r"$\chi_{\rm crit}$ reconciliation: $M_{\rm halo}(z=5)=5\times10^{10}\,M_\odot$", fontsize=8.5)
ax.legend(loc="lower right", fontsize=6.5, frameon=False)
fig.tight_layout()

outpath = OUTDIR / "fig_chi_crit_reconciliation.png"
fig.savefig(outpath, dpi=300)
print("saved", outpath)
