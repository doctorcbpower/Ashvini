"""
How close do fixed seeds get to the target f_BH = 0.5 across the tested accessibility range? Direct per-tree ratios
M_BH(z=5)/M*_tot(z=5) (median [16,84]) and the fraction of trees reaching the target, plus the critical seed, for
(sigma_j, R_nuc) cells of the earlier regime grid. Branches from the frozen MVM (hash asserted). 60 trees, 801-step dt.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, generators, trees, grid
import numpy as np
assert_frozen()
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm

h, zh, _ = generators()
CELLS = [(0.5, 100.0), (0.5, 250.0), (1.0, 100.0), (1.5, 100.0), (1.0, 300.0), (1.5, 300.0)]
SEEDS = (1e3, 2e5, 1e7)
tz = grid(25.0)
print("60 trees, 801-step dt. ratio = M_BH(z=5)/M*_tot(z=5), median [16,84]; reach = fraction of trees with ratio >= 0.5 (f_BH); Mcrit = median critical seed [Msun]")
for i, M0 in enumerate((3e10, 3e11)):
    hm, zz = trees(zh, h, M0, 60, 25.0, 6000 + i)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    print(f"\nM0={M0:.0e}")
    for (s, R) in CELLS:
        kw = dict(sigma_lnj=s, R_nuc_pc=R)
        Mc, *_ = mvm.critical_seed_stock(Mh, tz, n_iter=36, halo_growth_rate=rate, **kw)
        row = []
        for sd in SEEDS:
            o = mvm.run_reservoir_stock(Mh, tz, np.full(Mh.shape[0], sd), halo_growth_rate=rate, **kw)
            r = o["bh_mass"][:, -1] / o["stars_mass"][:, -1]
            row.append(f"seed {sd:.0e}: {np.median(r):.2e} [{np.percentile(r,16):.1e},{np.percentile(r,84):.1e}] max {r.max():.2e} reach {np.mean(r>=0.5)*100:4.1f}%")
        print(f"  sigma_j={s}, R_nuc={R:5.0f}: Mcrit {np.nanmedian(Mc):.2e} | " + " | ".join(row), flush=True)
