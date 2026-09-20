"""
MVM Figure 2: the inverse question. M_BH(z=5)/M*_tot(z=5) for fixed seeds, fiducial (solid) and a labelled
accessibility experiment (dashed; R_nuc = 250 pc, NOT a candidate fiducial), with the target f_BH lines.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).parent
OUTDIR = HERE / "output"
import paper_style
paper_style.apply()

D = json.load(open(OUTDIR / "mvm_production_results.json"))
R = sorted(D["results"], key=lambda r: r["M0"])
M = np.array([r["M0"] for r in R])
seeds = D["meta"]["fixed_seeds"]
labels = {1e2: r"$10^{2}$ (Pop III)", 1e3: r"$10^{3}$ (runaway)", 2e5: r"$2\times10^{5}$ (direct collapse)", 1e7: r"$10^{7}$"}
med = lambda case, s: np.array([np.nanmedian(r["fixed_seed"][case][str(s)]["ratio"]) for r in R])
band = lambda case, s: np.array([np.nanpercentile(r["fixed_seed"][case][str(s)]["ratio"], [16, 84]) for r in R])

fig, ax = plt.subplots(figsize=(paper_style.COL, 3.4))
SEEDC = paper_style.seq(4)
for k, s in enumerate(seeds):
    c = SEEDC[k]
    ax.fill_between(M, band("fiducial", s)[:, 0], band("fiducial", s)[:, 1], color=c, alpha=0.15, lw=0)
    ax.plot(M, med("fiducial", s), color=c, lw=1.5, ls="-", label=labels[s])
    ax.fill_between(M, band("accessible_R250", s)[:, 0], band("accessible_R250", s)[:, 1], color=c, alpha=0.10, lw=0)
    ax.plot(M, med("accessible_R250", s), color=c, lw=1.2, ls="--")
ax.axhline(0.5, color="k", lw=0.8)
ax.text(M.min() * 1.05, 0.55, r"target $f_{\rm BH}=0.5$", fontsize=7, va="bottom")
ax.axhline(0.1, color="k", lw=0.6, ls=":")
ax.text(M.min() * 1.05, 0.11, r"$f_{\rm BH}=0.1$", fontsize=7, va="bottom")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$M_{\rm halo}(z=5)\ [M_\odot]$")
ax.set_ylabel(r"$M_{\rm BH}(z=5)\,/\,M_{\star,{\rm tot}}(z=5)$")
ax.set_ylim(1e-11, 3)
h1, l1 = ax.get_legend_handles_labels()
style = [Line2D([0], [0], color="0.3", lw=1.5, ls="-"), Line2D([0], [0], color="0.3", lw=1.2, ls="--")]
leg1 = ax.legend(h1, l1, title=r"seed [$M_\odot$]", loc="lower left", fontsize=7, title_fontsize=7)
ax.add_artist(leg1)
ax.legend(style, [r"fiducial ($R_{\rm nuc}=100$ pc)", r"accessibility test ($R_{\rm nuc}=250$ pc)"], loc="upper right", fontsize=7,
          bbox_to_anchor=(1.0, 0.86))
out = OUTDIR / "mvm_fig2_inverse.png"
paper_style.save(fig, OUTDIR / "mvm_fig2_inverse")
print("saved", out)
for case in ("fiducial", "accessible_R250"):
    print(case)
    for s in seeds:
        print(f"  seed {s:.0e}: ratio at 3e10 {med(case, s)[0]:.2e}, 3e11 {med(case, s)[4]:.2e}, 3e13 {med(case, s)[-1]:.2e};  median G at 3e10 {np.nanmedian(R[0]['fixed_seed'][case][str(s)]['G']):.3g}")
