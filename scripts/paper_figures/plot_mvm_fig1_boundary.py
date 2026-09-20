"""
MVM Figure 1: critical seed boundary M_seed,crit(M_halo) and its ratio to the target f_BH M*_tot (= 1/G_BH).
Reads output/mvm_production_results.json (gen_mvm_production.py).
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

D = json.load(open(OUTDIR / "mvm_production_results.json"))
R = sorted(D["results"], key=lambda r: r["M0"])
FBH = D["meta"]["f_bh"]
M = np.array([r["M0"] for r in R])
pct = lambda key, f=lambda x: x: np.array([np.nanpercentile(f(np.array(r[key])), [16, 50, 84]) for r in R])
crit = pct("Mcrit")
tgt = np.array([np.nanpercentile(FBH * np.array(r["Mstar"]), [16, 50, 84]) for r in R])
ratio = np.array([np.nanpercentile(np.array(r["Mcrit"]) / (FBH * np.array(r["Mstar"])), [16, 50, 84]) for r in R])

fig, (ax, bx) = plt.subplots(2, 1, figsize=(3.5, 4.5), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1], hspace=0.06))
ax.fill_between(M, crit[:, 0], crit[:, 2], color="C0", alpha=0.25, lw=0, label=r"$M_{\rm seed,crit}$, 16--84\% (trees)")
ax.plot(M, crit[:, 1], color="C0", lw=1.8, label=r"$M_{\rm seed,crit}$, median")
ax.plot(M, tgt[:, 1], color="C3", lw=1.0, ls="--", label=r"$f_{\rm BH}M_{\star,{\rm tot}}(z=5)$, median")
# the three halo masses with the deepest diagnostic coverage (convergence, sensitivities, regime map): small unlabelled markers
diag = np.array([np.argmin(abs(M - m0)) for m0 in (3e10, 3e11, 3e13)])
ax.plot(M[diag], crit[diag, 1], ls="none", marker="o", ms=3.0, mfc="white", mec="C0", mew=0.8, zorder=6)
for lab, mass in (("Pop III", 1e2), ("runaway", 1e3), ("direct collapse", 2e5)):
    ax.axhline(mass, color="0.45", lw=0.6, ls=":", zorder=0)
    ax.text(M.max() * 0.93, mass * 1.25, lab, fontsize=6.5, ha="right", va="bottom", color="0.35")
ax.axvline(4e11, color="0.6", lw=0.6, ls="-.", zorder=0)
ax.text(4.4e11, 2.0e2, r"$M_{\rm hot}$", fontsize=7, color="0.4", va="bottom")
ax.set_yscale("log")
ax.set_ylabel(r"$M_{\rm seed,crit}\ [M_\odot]$")
ax.set_ylim(5e1, 3e10)
ax.legend(loc="upper left", fontsize=6.5)

bx.fill_between(M, ratio[:, 0], ratio[:, 2], color="C0", alpha=0.25, lw=0)
bx.plot(M, ratio[:, 1], color="C0", lw=1.8)
bx.plot(M[diag], ratio[diag, 1], ls="none", marker="o", ms=3.0, mfc="white", mec="C0", mew=0.8, zorder=6)
bx.axhline(1.0, color="k", lw=0.5, ls="-")
bx.set_xscale("log")
bx.set_xlabel(r"$M_{\rm halo}(z=5)\ [M_\odot]$")
bx.set_ylabel(r"$\dfrac{M_{\rm seed,crit}}{f_{\rm BH}M_{\star,{\rm tot}}}=\dfrac{1}{G_{\rm BH}}$")
lo = min(0.8, np.nanmin(ratio) - 0.03)
bx.set_ylim(lo, 1.02)
bx.set_xlim(M.min() * 0.9, M.max() * 1.1)
out = OUTDIR / "mvm_fig1_boundary.png"
fig.savefig(out, dpi=300)
print("saved", out)
for r, c, q in zip(R, crit[:, 1], ratio[:, 1]):
    print(f"  M0={r['M0']:.2e}: median M_seed,crit={c:.3e}, ratio={q:.3f}, bracket fails={r['bracket_fail']}, max|F|={r['max_abs_F']:.0e}")
