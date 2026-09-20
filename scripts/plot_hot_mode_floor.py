import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "scripts")
from stellar_to_halo_mass_relation import behroozi2013_mstar as B

d = json.load(open("scripts/output/hot_mode_floor_scan.json"))
m = np.array(d["masses"]); grid = np.logspace(9.8, 14.2, 100)
sel = [("none (baseline)", "tuned baseline (no hot-mode term)", "#1f77b4", "-"),
       ("M_hot=4e+11 floor=0.0", r"$M_{\rm hot}=4\times10^{11}$, floor 0 (MVM-like cutoff)", "#2ca02c", "--"),
       ("M_hot=1e+12 floor=0.1", r"$M_{\rm hot}=10^{12}$, floor 0.1", "#ff7f0e", "-"),
       ("M_hot=2e+12 floor=0.1", r"$M_{\rm hot}=2\times10^{12}$, floor 0.1", "#d62728", "-")]
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5))
for k, lab, c, ls in sel:
    ms = np.array(d["runs"][k])
    a.plot(m, ms, "o" + ls, color=c, lw=2.4, ms=5, label=lab)
    b.plot(m, ms / m, "o" + ls, color=c, lw=2.4, ms=5)
a.plot(grid, B(grid), ":", color="0.25", lw=2.4, label="Behroozi+2013")
b.plot(grid, B(grid) / grid, ":", color="0.25", lw=2.4)
for x in (a, b):
    x.set_xscale("log"); x.set_yscale("log"); x.grid(alpha=.3); x.tick_params(labelsize=11)
    x.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]", fontsize=14)
a.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$]", fontsize=14)
b.set_ylabel(r"$M_\star/M_{\rm halo}$", fontsize=14)
a.legend(fontsize=9, loc="upper left")
a.text(0.97, 0.03, "ILLUSTRATIVE: UV term, $\\epsilon_{\\rm ff}$ and hot-mode floor tuned\n(CDM, Zhang-Hui trees, 100 halos/bin)",
       transform=a.transAxes, ha="right", fontsize=8.5, bbox=dict(fc="w", ec="0.6", alpha=.85))
fig.tight_layout(); fig.savefig("shmr_hot_mode_floor.png", dpi=150); print("ok")
