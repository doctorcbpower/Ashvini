"""
Prints every production-derived number quoted in docs/mvm_numerical_crib_sheet.md, straight from
output/mvm_production_results.json (frozen MVM ensemble). No number in the crib sheet that comes from the
production run is typed by hand.
"""
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).parent
D = json.load(open(HERE / "output" / "mvm_production_results.json"))
meta, R = D["meta"], sorted(D["results"], key=lambda r: r["M0"])
FBH, FB = meta["f_bh"], meta["fiducial"]["f_b"]
pc = lambda x, q=(16, 50, 84): np.nanpercentile(np.array(x, dtype=float), q)

print("### META:", {k: meta[k] for k in ("date", "n_trees", "n_steps", "dz", "M_res", "f_bh", "module_sha256")})
print("\n### A. boundary, all masses: M0 | median Mcrit | 16 | 84 | 16-84 width [dex] | 1/G | G | M*/Mh | Mcrit/Mh | f_BH M* (median) | accreted/(f_b Mh) | bracket fails | max|F|")
for r in R:
    mc = pc(r["Mcrit"]); ms = pc(r["Mstar"]); rat = pc(np.array(r["Mcrit"]) / (FBH * np.array(r["Mstar"])))
    print(f"{r['M0']:.3e} | {mc[1]:.3e} | {mc[0]:.3e} | {mc[2]:.3e} | {np.log10(mc[2]/mc[0]):.3f} | {rat[1]:.4f} [{rat[0]:.4f},{rat[2]:.4f}] | {1/rat[1]:.4f} | {ms[1]/r['M0']:.3e} | {mc[1]/r['M0']:.3e} | {FBH*ms[1]:.3e} | {np.median(r['accreted_over_fb_Mhalo']):.3f} | {r['bracket_fail']} | {r['max_abs_F']:.0e}")
M = np.array([r["M0"] for r in R]); mcm = np.array([pc(r["Mcrit"])[1] for r in R])
sl = lambda i, j: np.log10(mcm[j] / mcm[i]) / np.log10(M[j] / M[i])
print(f"\n### local slopes dlogMcrit/dlogMh: 3e10->3e11 {sl(0,4):.3f};  3e11->1e12(~9.5e11) {sl(4,6):.3f};  9.5e11->3e12 {sl(6,8):.3f};  3e12->3e13 {sl(8,12):.3f};  whole range {sl(0,12):.3f}")
print("### tree scatter half-width [dex] (16-84)/2 at 3e10, 3e11, 3e13:", [f"{np.log10(pc(r['Mcrit'])[2]/pc(r['Mcrit'])[0])/2:.3f}" for r in (R[0], R[4], R[12])])
print(f"### 16-84 width range across masses [dex]: {min(np.log10(pc(r['Mcrit'])[2]/pc(r['Mcrit'])[0]) for r in R):.3f} to {max(np.log10(pc(r['Mcrit'])[2]/pc(r['Mcrit'])[0]) for r in R):.3f}")

print("\n### B. fixed seeds, achieved M_BH/M*_tot(z=5), median [16,84]; then median G")
for case in ("fiducial", "accessible_R250"):
    print(f"  -- {case}")
    for s in meta["fixed_seeds"]:
        row = []
        for idx in (0, 4, 12):
            q = pc(R[idx]["fixed_seed"][case][str(s)]["ratio"]); g = np.nanmedian(R[idx]["fixed_seed"][case][str(s)]["G"])
            row.append(f"{R[idx]['M0']:.0e}: {q[1]:.2e} [{q[0]:.1e},{q[2]:.1e}] G={g:.4g}")
        print(f"   seed {s:.0e}: " + " | ".join(row))
print("  -- accessibility test, seed 1e3, all masses: M0 | median ratio [16,84] | median G | 84th pct G | % trees with G>2 | fiducial % with G>2")
for r in R:
    a = r["fixed_seed"]["accessible_R250"]["1000.0"]; f = r["fixed_seed"]["fiducial"]["1000.0"]
    q = pc(a["ratio"]); G = np.array(a["G"])
    print(f"   {r['M0']:.2e} | {q[1]:.2e} [{q[0]:.1e},{q[2]:.1e}] | {np.median(G):.4g} | {np.percentile(G,84):.4g} | {np.mean(G>2)*100:.1f}% | {np.mean(np.array(f['G'])>2)*100:.1f}%")
best = max(np.nanpercentile(r["fixed_seed"]["accessible_R250"]["1000.0"]["ratio"], 84) for r in R)
print(f"   highest 84th-percentile ratio reached by a 1e3 seed in the accessibility test: {best:.2e}  ({np.log10(0.5/best):.2f} dex below f_BH=0.5)")
print(f"   largest median ratio (1e3, accessibility test): {max(np.nanmedian(r['fixed_seed']['accessible_R250']['1000.0']['ratio']) for r in R):.2e}")

print("\n### C. seed / host baryons: M0 | seed/(f_b Mh) at first resolved step, median [16,84] | median z_first | z where median M_BH/(f_b Mh) crosses 1 | median M_BH/(f_b Mh) at z=15,10,5 | median over trees of M_BH(z=5)/(f_b Mh)")
for r in R:
    z = np.array(r["z"]); p50 = np.array(r["bh_to_host_p50"], dtype=float)
    s = pc(r["seed_over_host_baryons_at_formation"]); zf = np.nanmedian(r["z_first_resolved"])
    ok = np.isfinite(p50); cross = "n/a"
    idx = np.where(ok[:-1] & ok[1:] & (p50[:-1] >= 1) & (p50[1:] < 1))[0]
    if len(idx): cross = f"{z[idx[0]]:.1f}"
    at = lambda zt: p50[np.argmin(abs(z - zt))]
    fin = np.median(np.array(r["MBH"]) / (FB * r["M0"]))
    print(f"   {r['M0']:.2e} | {s[1]:.2e} [{s[0]:.1e},{s[2]:.1e}] | {zf:.1f} | {cross} | {at(15):.2e}, {at(10):.2e}, {at(5):.2e} | {fin:.2e}")
print("done")
