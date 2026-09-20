#!/usr/bin/env python3
"""Plot the z=0 SHMR from the general model and from the MVM on the same Zhang-Hui CDM trees
(data: scripts/output/shmr_mvm_vs_general_cdm_dz0.005.json, made by scripts/shmr_mvm_vs_general.py)."""
import json
import numpy as np
import matplotlib.pyplot as plt
from stellar_to_halo_mass_relation import moster2013_mstar, behroozi2013_mstar

rows = [r for r in json.load(open("scripts/output/shmr_mvm_vs_general_cdm_dz0.005.json")) if "error" not in r]
mh = np.array([r["M0"] for r in rows])
grid = np.logspace(7.8, 14.2, 200)


def series(key):
    a = np.array([[r[key][k] for k in ("p16", "med", "p84")] for r in rows])
    return np.where(a > 0, a, np.nan)


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5))
for key, color, label, ls in (("general", "#1f77b4", "general model (UV suppression, delayed SN)", "-"),
                              ("mvm", "#d62728", "MVM (high-z model, extrapolated to $z=0$)", "-")):
    s = series(key)
    ax.fill_between(mh, s[:, 0], s[:, 2], color=color, alpha=0.18, lw=0)
    ax.plot(mh, s[:, 1], "o" + ls, color=color, lw=2.6, ms=6, label=label)
    ax2.plot(mh, s[:, 1] / mh, "o" + ls, color=color, lw=2.6, ms=6)
ax.plot(mh, series("mvm_grumpy")[:, 1], "--", color="#d62728", lw=1.6, label="MVM, GRUMPY-smoothed halo growth rate")
ax2.plot(mh, series("mvm_grumpy")[:, 1] / mh, "--", color="#d62728", lw=1.6)
for f, ls, lab in ((moster2013_mstar, "--", "Moster+2013"), (behroozi2013_mstar, ":", "Behroozi+2013")):
    ax.plot(grid, f(grid), ls, color="0.3", lw=2.2, label=lab)
    ax2.plot(grid, f(grid) / grid, ls, color="0.3", lw=2.2)
for a in (ax, ax2):
    a.set_xscale("log"); a.set_yscale("log"); a.grid(alpha=0.3); a.tick_params(labelsize=12)
    a.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]", fontsize=15)
ax.set_ylim(1e-1, 3e12); ax2.set_ylim(1e-6, 1e-1)
ax.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$]", fontsize=15)
ax2.set_ylabel(r"$M_\star/M_{\rm halo}$ (median)", fontsize=15)
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig("shmr_mvm_vs_general.png", dpi=150)
print("Wrote shmr_mvm_vs_general.png")
