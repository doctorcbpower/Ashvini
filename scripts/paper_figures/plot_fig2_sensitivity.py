"""
Plots Figure 2 (the three-panel eta_acc / epsilon_f / spin sensitivity
figure) from gen_fig2_sensitivity.py's output JSON.
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

with open(OUTDIR / "fig2_sensitivity_results.json") as f:
    data = json.load(f)
results = data["results"]


def get_series(name):
    d = results[name]
    keys = sorted(d.keys(), key=lambda k: float(k))
    m = np.array([float(k) for k in keys])
    med = np.array([d[k]["median"] for k in keys])
    p16 = np.array([d[k]["p16"] for k in keys])
    p84 = np.array([d[k]["p84"] for k in keys])
    return m, med, p16, p84


m_fid, med_fid, p16_fid, p84_fid = get_series("fiducial")

fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), sharey=True)

panels = [
    ("(a) Nuclear-supply efficiency $\\eta_{\\rm acc}$", [
        ("eta_acc_lo", r"$\eta_{\rm acc}=10^{-3}$", "C1"),
        ("eta_acc_hi", r"$\eta_{\rm acc}=10^{-1}$", "C2"),
    ]),
    ("(b) AGN feedback coupling $\\epsilon_f$", [
        ("epsilon_f_lo", r"$\epsilon_f=5\times10^{-5}$", "C1"),
        ("epsilon_f_hi", r"$\epsilon_f=5\times10^{-3}$", "C2"),
    ]),
    ("(c) Radiative efficiency $\\epsilon$ (spin)", [
        ("spin_chaotic", r"chaotic floor, $\epsilon=0.057$", "C1"),
        ("spin_coherent", r"coherent limit, $\epsilon=0.32$", "C2"),
    ]),
]

for ax, (title, curves) in zip(axes, panels):
    ax.fill_between(m_fid, p16_fid, p84_fid, color="C0", alpha=0.2)
    ax.plot(m_fid, med_fid, color="C0", lw=2.0, label="fiducial")
    for name, label, color in curves:
        m, med, p16, p84 = get_series(name)
        ax.plot(m, med, color=color, lw=1.6, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\rm halo}(z=5)\ [M_\odot]$")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5, frameon=False)

axes[0].set_ylabel(r"$M_{\rm seed,crit}\ [M_\odot]$")
fig.tight_layout()

outpath = OUTDIR / "fig_sensitivity_panels.png"
fig.savefig(outpath, dpi=200)
print("saved", outpath)
