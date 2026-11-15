"""
Nuclear-delivery experiments of the black hole paper, redone with the revised model (reservoir_stock_v2 defaults:
nuclear dark matter, conservative inflow, one cosmology, sub-stepping) on the stored production trees.

One halo mass per invocation (--mass-index 0 for 3e10, 4 for 3e11, 12 for 3e13 Msun), reading the trees that
gen_mvm_v2_production.py stored, so every cell of every table uses the same trees (as in the paper's tables).

Three blocks, each matching a table or figure of the paper:
  gmap   growth factor G = M_BH(z=5)/M_seed of a 1e3 Msun seed over a grid of sigma_j and R_nuc          (tab:gmap)
  reach  critical seed and final M_BH/M_* of fixed seeds 1e3, 2e5, 1e7 in six (sigma_j, R_nuc) cells      (tab:reach)
  seeds  M_BH/M_* against redshift and at z=5 for seeds 1e2, 1e3, 1e5, 1e7 in the fiducial and the
         high-delivery (sigma_j=1.5, R_nuc=300 pc) configurations                                         (seed-dependence section)

Usage: python gen_mvm_v2_delivery.py --mass-index I --trees-from <dir>/trees_MII.npz --outdir <dir>
Writes mvm_v2_delivery_M<II>.json.  Cost is a few minutes per mass.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

MASSES = [3.0e10 * 10 ** (0.25 * k) for k in range(13)]
SIGMA_GRID = [0.5, 0.6, 0.7, 0.8, 1.0]
RNUC_GRID = [100.0, 125.0, 150.0, 200.0, 250.0, 300.0]
REACH_CELLS = [(0.5, 100.0), (0.5, 250.0), (1.0, 100.0), (1.5, 100.0), (1.0, 300.0), (1.5, 300.0)]
REACH_SEEDS = [1e3, 2e5, 1e7]
SEED_SET = [1e2, 1e3, 1e5, 1e7]
SEED_CELLS = {"fiducial": (0.5, 100.0), "high_delivery": (1.5, 300.0)}


def pct(a, q=(16, 50, 84)):
    return [float(x) for x in np.nanpercentile(np.asarray(a, float), q)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass-index", type=int, required=True)
    ap.add_argument("--trees-from", required=True)
    ap.add_argument("--n-trees", type=int, default=60, help="trees for gmap and reach (the paper used 60 in tab:reach)")
    ap.add_argument("--n-trees-seeds", type=int, default=100, help="trees for the seed-dependence block")
    ap.add_argument("--n-steps", type=int, default=801)
    ap.add_argument("--n-sub", type=int, default=8)
    ap.add_argument("--blocks", default="gmap,reach,seeds")
    ap.add_argument("--outdir", default=str(HERE / "output" / "mvm_v2"))
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    idx = args.mass_index
    M0 = MASSES[idx]
    from ashvini import reservoir_stock_v2 as v2

    d = np.load(args.trees_from)
    hm_all, zz = d["hm"], d["zz"]
    cosmo = v2.PaperCosmology()
    tz = v2.uniform_time_grid(25.0, 5.0, args.n_steps, cosmo)
    Mh_all = v2.interpolate_tree_onto_grid_v2(hm_all, zz, tz, cosmo)
    # revised physics throughout (conservative inflow needs no growth-rate estimate)
    base = dict(dm_nuclear=True, inflow="conservative", cosmology=cosmo, n_sub=args.n_sub)
    blocks = [b.strip() for b in args.blocks.split(",")]
    res = dict(M0=M0, mass_index=idx, n_steps=args.n_steps, n_sub=args.n_sub, trees_file=str(args.trees_from),
               tree_algorithm=str(d["tree_algorithm"]) if "tree_algorithm" in d.files else "unknown")
    t0 = time.time()

    if "gmap" in blocks:
        Mh = Mh_all[:args.n_trees]
        res["gmap"] = dict(n_trees=int(Mh.shape[0]), seed=1e3, sigma=SIGMA_GRID, Rnuc=RNUC_GRID, G_p16_p50_p84=[])
        for s in SIGMA_GRID:
            row = []
            for R in RNUC_GRID:
                r = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], 1e3), sigma_lnj=s, R_nuc_pc=R, **base)
                row.append(pct(r["bh_mass"][:, -1] / 1e3))
            res["gmap"]["G_p16_p50_p84"].append(row)
        print(f"gmap done [{time.time()-t0:.0f}s]", flush=True)

    if "reach" in blocks:
        Mh = Mh_all[:args.n_trees]
        res["reach"] = dict(n_trees=int(Mh.shape[0]), cells=[])
        for s, R in REACH_CELLS:
            kw = dict(base, sigma_lnj=s, R_nuc_pc=R)
            Mc, never, above, F = v2.critical_seed_v2(Mh, tz, n_iter=40, return_F=True, **kw)
            cell = dict(sigma_lnj=s, R_nuc_pc=R, Mcrit=Mc.tolist(), bracket_fail=int(never.sum() + above.sum()), ratio={})
            for sd in REACH_SEEDS:
                r = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], sd), **kw)
                cell["ratio"][str(sd)] = (r["bh_mass"][:, -1] / r["stars_mass"][:, -1]).tolist()
            res["reach"]["cells"].append(cell)
            print(f"reach ({s}, {R:.0f}) median Mcrit {np.nanmedian(Mc):.3e} [{time.time()-t0:.0f}s]", flush=True)

    if "seeds" in blocks:
        Mh = Mh_all[:args.n_trees_seeds]
        res["seeds"] = dict(n_trees=int(Mh.shape[0]), z=tz.tolist(), cells={})
        for name, (s, R) in SEED_CELLS.items():
            res["seeds"]["cells"][name] = dict(sigma_lnj=s, R_nuc_pc=R, by_seed={})
            for sd in SEED_SET:
                r = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], sd), sigma_lnj=s, R_nuc_pc=R, **base)
                with np.errstate(all="ignore"):
                    ratio = r["bh_mass"] / r["stars_mass"]
                    good = (r["stars_mass"] > 1e5).sum(axis=0) >= 20      # plot only where M_* > 1e5 Msun, as in the paper
                    band = np.where(good[None, :], np.nanpercentile(ratio, [16, 50, 84], axis=0), np.nan)
                res["seeds"]["cells"][name]["by_seed"][str(sd)] = dict(
                    p16=band[0].tolist(), p50=band[1].tolist(), p84=band[2].tolist(),
                    final_ratio=ratio[:, -1].tolist(), G=(r["bh_mass"][:, -1] / sd).tolist(),
                    delivered_fraction=(r["transfer_nuc"].sum(axis=1) / np.maximum(r["cum_in"], 1e-300)).tolist())
            print(f"seeds {name} [{time.time()-t0:.0f}s]", flush=True)

    res["seconds"] = time.time() - t0
    with open(out / f"mvm_v2_delivery_M{idx:02d}.json", "w") as f:
        json.dump(res, f)
    print(f"DONE M0={M0:.3e} in {res['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
