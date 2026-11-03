"""Plots from the stored wind-scaling scan (scripts/wind_scaling_scan.py): $SHMR_OUT/wind_scan/*.json -> wind_scan_*.png."""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import shmr_z0_common as C
from stellar_to_halo_mass_relation import behroozi2013_mstar, moster2013_mstar

# Nadler et al. (2020, ApJ 893, 48) Fig. 6 (right): the 68 percent band of the stellar-halo mass relation inferred from the DES + PS1 satellites, as a
# function of PEAK halo mass of satellites. Read off the figure by eye (about 0.2 dex), not tabulated by the authors.
NADLER_MPEAK = np.array([4.0e7, 1.6e8, 6.4e8, 2.6e9, 1.0e10, 4.2e10, 1.0e11])
NADLER_LOG_MSTAR_LO = np.array([1.48, 2.31, 3.16, 4.27, 5.33, 6.39, 6.88])
NADLER_LOG_MSTAR_HI = np.array([2.37, 3.40, 4.65, 6.04, 6.92, 7.88, 8.34])
NADLER_LABEL = "Nadler+20 satellites, 68%"
D = os.path.join(C.OUT, "wind_scan")
runs = {d["name"]: d for d in (json.load(open(f)) for f in sorted(glob.glob(os.path.join(D, "*.json"))))}
mh = np.logspace(8, 12, 100)
ref = [(np.array([float(f(m)) for m in mh]) / mh, s, l) for f, s, l in ((moster2013_mstar, "k--", "Moster+13"), (behroozi2013_mstar, "k-", "Behroozi+13"))]
ash = [n for n in runs if n.startswith("ash")]
st = [n for n in runs if n.startswith("stellar")]


def label(n):
    if n.startswith("ash"):
        m, z = n.split("_")[1:3]
        return f"mass exp {float(m[1:]):.2f}, z exp {float(z[1:]):.2f}"
    p, c = n.split("_")[1:3]
    return f"stellar-mass, slope {float(p[1:]):.2f}, {'cutoff' if c == 'cut' else 'no cutoff'}"


# 1. z = 0 relation
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
for ax, group, title in zip(axes, (ash, st), ("Ashvini loading: halo-mass and redshift exponents", "Stellar-mass loading (Kravtsov & Manwadkar eq. 34)")):
    for y, s, l in ref:
        ax.loglog(mh, y, s, color="0.3", label=l)
    for n, c in zip(group, plt.cm.viridis(np.linspace(0, 0.9, len(group)))):
        ax.loglog(C.MASSES, runs[n]["mstar_over_mhalo_z0"], "o-", color=c, ms=4, label=label(n))
    ax.fill_between(NADLER_MPEAK, 10 ** NADLER_LOG_MSTAR_LO / NADLER_MPEAK, 10 ** NADLER_LOG_MSTAR_HI / NADLER_MPEAK, color="C3", alpha=0.25, label=NADLER_LABEL)
    ax.axvspan(3e10, 1.2e12, color="0.9", zorder=0)
    ax.set(xlabel=r"$M_{\rm halo}$ at z = 0 [Msun]", title=title, ylim=(3e-6, 0.1), xlim=(8e7, 1.2e12)); ax.legend(fontsize=7)
axes[0].set_ylabel(r"median $M_\star/M_{\rm halo}$")
fig.text(0.5, -0.01, "grey band: masses the amplitude was set on (3e10 to 1e12 Msun); everything left of it is a prediction", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(D, "wind_scan_z0_relation.png"), dpi=130, bbox_inches="tight")

# 2. the loading itself
from ashvini import supernovae_feedback as sn_fb
from ashvini.config import DEFAULT, override
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for n, c in zip(ash, plt.cm.viridis(np.linspace(0, 0.9, len(ash)))):
    r = runs[n]; cfg = override(DEFAULT, **{"sn.mass_slope": float(n.split("_")[1][1:]), "sn.redshift_slope": float(n.split("_")[2][1:]), "sn.epsilon_p": r["amplitude"]})
    axes[0].loglog(mh, sn_fb.mass_loading_factor(7.0, mh, 0.0, cfg), color=c, label=label(n))
axes[0].set(xlabel=r"$M_{\rm halo}$ [Msun]", ylabel=r"mass loading $\eta_w$ at z = 7, metal-poor", title="Ashvini loading with the fitted amplitudes"); axes[0].legend(fontsize=7)
ms = np.logspace(3, 10.5, 200)
for n, c in zip(st, plt.cm.viridis(np.linspace(0, 0.9, len(st)))):
    r = runs[n]; p = float(n.split("_")[1][1:]); cc = 4.5 if n.endswith("_cut") else 0.0
    axes[1].loglog(ms, np.minimum(2000, np.maximum(1e-3, r["amplitude"] * (ms / 1e10) ** -p - cc)), color=c, label=label(n))
axes[1].set(xlabel=r"$M_\star$ [Msun]", ylabel=r"$\eta_w$", title="Stellar-mass loading with the fitted amplitudes", ylim=(1e-2, 3e3)); axes[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(D, "wind_scan_loading.png"), dpi=130)

# 3. redshift evolution at fixed halo mass bin
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
zs = [0, 2, 5]
for ax, key in zip(axes, ("1e8-3e8", "1e9-3e9", "1e10-3e10")):
    for n in ash + st:
        y = [runs[n]["mstar_over_mhalo_by_z"][f"z{z}"][key] if z else None for z in zs[1:]]
        z0 = {"1e8-3e8": 0, "1e9-3e9": 2, "1e10-3e10": 4}[key]
        vals = [runs[n]["mstar_over_mhalo_z0"][z0]] + y
        ok = [v is not None for v in vals]
        ax.semilogy(np.array(zs)[ok], np.array(vals, dtype=float)[ok], "o-" if n.startswith("ash") else "s--", ms=4, label=label(n))
    ax.set(xlabel="redshift", title=f"M_halo(z) = {key} Msun"); ax.set_xticks(zs)
axes[0].set_ylabel(r"median $M_\star/M_{\rm halo}$"); axes[2].legend(fontsize=6, loc="center left", bbox_to_anchor=(1.01, 0.5))
fig.tight_layout(); fig.savefig(os.path.join(D, "wind_scan_redshift.png"), dpi=130, bbox_inches="tight")

# which forms fall inside the Nadler band
print("Nadler+20 band (by eye) and the median M*/M_halo of each form at z = 0:")
for m in (1e8, 10 ** 8.5, 1e9, 10 ** 9.5, 1e10):
    lo = 10 ** np.interp(np.log10(m), np.log10(NADLER_MPEAK), NADLER_LOG_MSTAR_LO) / m
    hi = 10 ** np.interp(np.log10(m), np.log10(NADLER_MPEAK), NADLER_LOG_MSTAR_HI) / m
    k = int(np.argmin(abs(np.log10(C.MASSES) - np.log10(m))))
    inside = [n for n in runs if lo <= runs[n]["mstar_over_mhalo_z0"][k] <= hi]
    print(f"  M = {C.MASSES[k]:.1e}: band {lo:.1e} to {hi:.1e}; forms inside: {len(inside)} of {len(runs)}; range of forms {min(runs[n]['mstar_over_mhalo_z0'][k] for n in runs):.1e} to {max(runs[n]['mstar_over_mhalo_z0'][k] for n in runs):.1e}")
