"""
Reservoir time-step convergence of the revised MVM (variant E: dm_nuclear, conservative inflow, one cosmology) on stored trees.

The trees written by gen_mvm_v2_production.py (trees_M<idx>.npz) are interpolated onto grids of different numbers of nodes
uniform in cosmic time, and the critical seed, final stellar mass, nuclear stellar fraction and the fraction of live steps in
which outflows are limited by the available gas are compared with the production grid (801 nodes). Paired: the same trees in
every column. Nothing else changes between columns.

Usage: python v2_timestep_convergence.py TREES.npz [--n-trees 20] [--steps 401,801,801:4,801:16,1601:16] [--fixed-seed 1e3]
Each entry is NODES or NODES:N_SUB (the consumption and feedback part of each step split into N_SUB sub-steps; default 1).
Cost grows as the square of the number of nodes (one gas cohort per step): 3201 nodes takes about 16 times as long as 801.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ap = argparse.ArgumentParser()
ap.add_argument("trees")
ap.add_argument("--n-trees", type=int, default=20)
ap.add_argument("--steps", default="401,801,1601")
ap.add_argument("--fixed-seed", type=float, default=1e3, help="also report the final state of this fixed seed (Msun)")
ap.add_argument("--n-iter", type=int, default=30, help="bisection iterations (30 gives 1e-7 of the log bracket)")
a = ap.parse_args()

from ashvini import reservoir_stock_v2 as v2

d = np.load(a.trees)
hm, zz = d["hm"][: a.n_trees], d["zz"]
cos = v2.PaperCosmology()
print(f"{a.trees}: M0 = {float(d['M0']):.3e} Msun, {hm.shape[0]} trees, variant E (revised model)")
print("%6s %5s %9s %14s %14s %12s %12s %10s %12s" % ("nodes", "n_sub", "dt[Myr]", "med Mseed,crit", "ratio to 801:1", "med Mstar", "nuc/Mstar", "wind-cap", "G(seed)"))
ref = None
rows = []
for tok in a.steps.split(","):
    n, nsub = (int(x) for x in (tok.split(":") + ["1"])[:2])
    t0 = time.time()
    tz = v2.uniform_time_grid(25.0, 5.0, n, cos)
    Mh = v2.interpolate_tree_onto_grid_v2(hm, zz, tz, cos)
    Mc, never, above = v2.critical_seed_v2(Mh, tz, n_iter=a.n_iter, n_sub=nsub)
    o = v2.run_reservoir_v2(Mh, tz, Mc, n_sub=nsub)
    f = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], a.fixed_seed), n_sub=nsub)
    live = Mh[:, :-1] > 0
    wcap = ((o["wind_cap"][:, 1:] & live).sum(axis=1) / np.maximum(live.sum(axis=1), 1))
    rows.append((n, nsub, 1e3 * float(np.diff(cos.time_at_z(tz))[0]), np.nanmedian(Mc), Mc, np.median(o["stars_mass"][:, -1]),
                 np.median(o["stars_nuc"][:, -1] / o["stars_mass"][:, -1]), wcap.mean(), np.median(f["bh_mass"][:, -1] / a.fixed_seed)))
    print("  [%d nodes, n_sub %d done in %.0fs]" % (n, nsub, time.time() - t0), file=sys.stderr, flush=True)
base = next((r for r in rows if r[0] == 801 and r[1] == 1), rows[0])
for n, nsub, dtm, med, Mc, Ms, nf, wc, G in rows:
    pair = np.nanmedian(Mc / base[4])
    print("%6d %5d %9.2f %14.4e %14.4f %12.3e %12.3f %10.3f %12.3g" % (n, nsub, dtm, med, pair, Ms, nf, wc, G))
