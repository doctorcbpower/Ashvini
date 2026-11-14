"""
Merge the per-mass outputs of gen_mvm_v2_production.py into one JSON and print the table of variant effects.

Usage: python scripts/paper_figures/merge_mvm_v2_production.py [DIR]     (DIR defaults to output/mvm_v2)
Writes DIR/mvm_v2_production_results.json (same 'results' layout as a list over masses) and prints, per mass, the median
critical seed of each variant and its ratio to variant A (the frozen model) on the same trees, with the paired scatter.
"""
import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "output" / "mvm_v2"
files = sorted(d.glob("mvm_v2_M??.json"))
if not files:
    sys.exit(f"no mvm_v2_M??.json in {d}")
runs = [json.load(open(f)) for f in files]
json.dump(dict(meta=[r["meta"] for r in runs], results=[r["results"] for r in runs]), open(d / "mvm_v2_production_results.json", "w"))
names = list(runs[0]["results"])
print("median M_seed,crit [Msun] per variant, and median of the per-tree ratio to A_frozen (paired trees)")
print("%9s " % "M0" + " ".join("%22s" % n for n in names))
for r in runs:
    base = np.array(r["results"]["A_frozen"]["Mcrit"], float) if "A_frozen" in r["results"] else None
    cells = []
    for n in names:
        mc = np.array(r["results"][n]["Mcrit"], float)
        cell = "%.3e" % np.nanmedian(mc)
        if base is not None:
            cell += " (x%.3f)" % np.nanmedian(mc / base)
        cells.append("%22s" % cell)
    print("%9.2e " % r["meta"]["M0"] + " ".join(cells))
print()
print("single-crossing check (variant F): trees whose F(M_seed) does not cross zero exactly once")
tot, bad = 0, 0
for r in runs:
    sc = r["results"].get("F_full_nsub", {}).get("single_crossing")
    if sc:
        tot += len(sc["n_sign_changes"])
        bad += sc["trees_with_n_changes_ne_1"]
print(f"  {bad} of {tot} trees")
print("largest baryon-conservation residual over all variants and masses: %.1e" % max(
    v["baryon_residual_max"] for r in runs for v in r["results"].values()))
