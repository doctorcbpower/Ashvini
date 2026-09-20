"""
MVM Figure 3: M_BH / (f_b M_halo) against redshift for the critical seed. Thick line: representative tree (median
M_seed,crit); band: 16--84 percentiles over the trees whose host is resolved; the marker is the first resolved step.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
OUTDIR = HERE / "output"
import paper_style
paper_style.apply()

D = json.load(open(OUTDIR / "mvm_production_results.json"))
R = {round(np.log10(r["M0"]), 2): r for r in D["results"]}
MC = paper_style.seq(3)
picks = [(10.48, MC[0], r"$3\times10^{10}$"), (11.48, MC[1], r"$3\times10^{11}$"), (13.48, MC[2], r"$3\times10^{13}$")]

fig, ax = plt.subplots(figsize=(paper_style.COL, 3.0))
for key, c, lab in picks:
    r = R[key]
    z = np.array(r["z"])
    p16, p50, p84 = (np.array(r[k], dtype=float) for k in ("bh_to_host_p16", "bh_to_host_p50", "bh_to_host_p84"))
    rep = np.array(r["rep_bh_to_host"], dtype=float)
    ok = np.isfinite(p16) & np.isfinite(p84)
    ax.fill_between(z[ok], p16[ok], p84[ok], color=c, alpha=0.18, lw=0)
    ax.plot(z, rep, color=c, lw=1.6, label=r"$M_{\rm halo}(z{=}5)=$" + lab + r" $M_\odot$")
    first = np.nanargmax(np.isfinite(rep))
    ax.plot(z[first], rep[first], marker="o", ms=3.5, color=c, mec="k", mew=0.4, zorder=5)
ax.axhline(1.0, color="k", lw=0.6)
ax.text(24.7, 1.25, r"$M_{\rm BH}=f_{\rm b}M_{\rm halo}$", fontsize=7, va="bottom")
ax.set_yscale("log")
ax.set_xlim(25, 5)
ax.set_ylim(1e-4, 3e6)
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"$M_{\rm BH}\,/\,(f_{\rm b}M_{\rm halo})$")
ax.legend(loc="upper right", fontsize=7)
out = OUTDIR / "mvm_fig3_hostbaryons.png"
paper_style.save(fig, OUTDIR / "mvm_fig3_hostbaryons")
print("saved", out)
for key, c, lab in picks:
    r = R[key]
    s = np.array(r["seed_over_host_baryons_at_formation"], dtype=float)
    zf = np.array(r["z_first_resolved"], dtype=float)
    z = np.array(r["z"]); rep = np.array(r["rep_bh_to_host"], dtype=float)
    def at(zt): return rep[np.argmin(abs(z - zt))]
    print(f"M0={r['M0']:.1e}: seed/(f_b M_halo) at formation median {np.nanmedian(s):.2e} [16,84 = {np.nanpercentile(s,16):.1e}, {np.nanpercentile(s,84):.1e}], z_first median {np.nanmedian(zf):.1f}; rep tree ratio at z=15 {at(15):.2e}, z=10 {at(10):.2e}, z=5 {at(5):.2e}")
