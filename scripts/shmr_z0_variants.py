"""
Run the variants (scripts/shmr_z0_common.VARIANTS) of Ashvini's grid integrator on the stored trees of one halo mass, and
store the z = 0 state of every halo. Resume-safe. One cluster task per mass index (all variants), or name a variant.

Usage: python scripts/shmr_z0_variants.py <mass_index 0..8 | all> [variant ...]
"""
import os
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C

KEYS = ("gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass", "bh_mass", "sfr")


def run_one(variant, index):
    path = C.result_path(variant, index)
    if os.path.exists(path):
        return f"{variant} {index}: exists, skipped"
    from ashvini.config import DEFAULT, override
    from ashvini.main import run_forest
    t = np.load(C.tree_path(index))
    cfg = override(DEFAULT, **C.VARIANTS[variant])
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = run_forest(t["halo_masses"], t["halo_growth_rates"], t["redshifts"], config=cfg)
    fast = any("finer grid" in str(w.message) for w in caught)
    out = {k: r[k][:, -1] for k in KEYS}
    out["halo_mass_z0"] = r["halo_mass"][:, -1]
    out["mass_bin"] = t["mass_bin"]
    out["fast_star_formation_warning"] = fast
    out["finite"] = bool(all(np.all(np.isfinite(r[k])) for k in ("gas_mass", "stars_mass", "gas_metals")))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp.npz"
    np.savez(tmp, **out)
    os.replace(tmp, path)
    return f"{variant} {index}: {time.time() - t0:.1f} s, median M* {np.median(out['stars_mass']):.3e}" + (" [fast star formation warning]" if fast else "")


if __name__ == "__main__":
    idx = range(len(C.MASSES)) if sys.argv[1] == "all" else [int(sys.argv[1])]
    names = sys.argv[2:] or list(C.VARIANTS)
    for i in idx:
        for v in names:
            print(run_one(v, i), flush=True)
