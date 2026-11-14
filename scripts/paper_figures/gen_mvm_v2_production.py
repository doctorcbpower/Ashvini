"""
Production ensemble for the REVISED Minimal Viable Model (ashvini.reservoir_stock_v2), one halo mass per invocation.

Why this exists. The frozen production (gen_mvm_production.py) used (i) a halo-growth rate estimator that over-delivers baryons
on merger trees, (ii) no dark matter in the nuclear enclosed mass, and (iii) two cosmologies for time and halo structure
(see model-verification-audit-2026-11-14.md). Here the same trees are run through five variants so that the effect of each
change is measured on identical trees, and the revised model is the last one:

    A_frozen              dm_nuclear=False  inflow=rate          cosmology=Planck18 tables   (reproduces the frozen model)
    B_dm_only             dm_nuclear=True   inflow=rate          cosmology=Planck18 tables
    C_conservative_only   dm_nuclear=False  inflow=conservative  cosmology=Planck18 tables
    D_conservative_cosmo  dm_nuclear=False  inflow=conservative  cosmology=paper (one flat LCDM)
    E_full                dm_nuclear=True   inflow=conservative  cosmology=paper              (all three physics changes)
    F_full_nsub           E with the consumption and feedback part of each step split into 8 sub-steps (the revised model;
                          removes the first-order step dependence of the nuclear stellar fraction, the light-seed growth
                          factor and the fraction of steps in which outflows are gas-limited, see
                          scripts/paper_figures/diagnostics/v2_timestep_convergence.py)

For every variant: the critical seed per tree (f_BH = 0.5, bisection as in the frozen script), and at that seed the final
M_star, M_BH, growth factor, baryon-budget ratio, baryon-conservation residual and the fractions of steps limited by
gas. For A and F also the fixed-seed fates (seeds 1e2, 1e3, 2e5, 1e7 Msun) in three configurations: fiducial, the
accessibility experiment (R_nuc = 250 pc) and the high-delivery configuration used in the paper (sigma_j = 1.5, R_nuc = 300 pc).
For F also the single-crossing check of F(M_seed) on a 33-point grid for EVERY tree, and the M_BH / (f_b M_halo)
percentiles against redshift.

The trees themselves are written to trees_M<idx>.npz because the numba sampler is not seed-reproducible: that file, not the
sampler, is the record, and any later analysis can be repeated on the same trees.

Usage (from the repository root, FORAOIS_ROOT set; see docs/PRODUCTION_PROVENANCE.md and hpc/submit_mvm_v2.sbatch):
    python scripts/paper_figures/gen_mvm_v2_production.py --mass-index 0 --outdir /fred/oz009/cpower/Ashvini/mvm_v2
    options: --n-trees 240 --n-steps 801 --dz 0.05 --m-res 1e4 --crossing-grid 33 --no-fixed-seeds --trees-from FILE.npz
"""
import argparse
import datetime
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np

MASSES = [3.0e10 * 10 ** (0.25 * k) for k in range(13)]
FIXED_SEEDS = [1e2, 1e3, 2e5, 1e7]
CONFIGS = {"fiducial": dict(), "accessible_R250": dict(R_nuc_pc=250.0), "high_delivery": dict(sigma_lnj=1.5, R_nuc_pc=300.0)}
VARIANTS = {
    "A_frozen": dict(dm_nuclear=False, inflow="rate", cosmology=None),
    "B_dm_only": dict(dm_nuclear=True, inflow="rate", cosmology=None),
    "C_conservative_only": dict(dm_nuclear=False, inflow="conservative", cosmology=None),
    "D_conservative_cosmo": dict(dm_nuclear=False, inflow="conservative", cosmology="paper"),
    "E_full": dict(dm_nuclear=True, inflow="conservative", cosmology="paper"),
    "F_full_nsub": dict(dm_nuclear=True, inflow="conservative", cosmology="paper", n_sub=8),
}
FIXED_SEED_VARIANTS = ("A_frozen", "F_full_nsub")
FINAL_VARIANT = "F_full_nsub"


def git_state(path):
    try:
        head = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=30).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"], capture_output=True,
                                    text=True, timeout=30).stdout.strip())
        return dict(head=head, dirty=dirty)
    except Exception as exc:  # provenance is best effort
        return dict(head="unknown", dirty=None, error=str(exc))


def build_trees(M0, idx, args):
    """Returns (halo_masses, redshifts, tree_stats).  tree_stats is None for Zhang & Hui; for PCH08 it holds, per tree,
    the number of steps, resolved splits, and steps that hit the minimum-step floor (the last should be zero)."""
    from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG
    sys.path.insert(0, FORAOIS_SRC)
    if args.tree_algorithm == "pch08":
        # Compliant PCH08: adaptive step per branch with expected splits per step <= 0.1, main progenitor only.
        from ashvini import pch08_main_branch
        hm, zz, stats = pch08_main_branch.build_main_branches(
            FORAOIS_CONFIG, M0, args.n_trees, z0=5.0, z_max=25.0, m_res_msun=args.m_res, dz_grid=args.dz, seed0=1000 * (idx + 1),
            target_nupper=args.target_nupper)
        return np.asarray(hm), np.asarray(zz), stats
    from foraois import cosmo_utils, ZhangHuiMergerTree
    from foraois.utils import io
    from ashvini import pymctrees_adapter
    rp = io.get_params(FORAOIS_CONFIG)
    h = rp["Cosmology"]["h"]
    cd = cosmo_utils.CosmoData(rp, redshift=[5.0])
    tg = ZhangHuiMergerTree(cd, rp, model="cdm")
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(
        tg, M0, h, args.n_trees, z0=5.0, z_max=25.0, m_res_msun=args.m_res, dz=args.dz, backend="numba", rng_seed=1000 + idx)
    return np.asarray(hm), np.asarray(zz), None


def summarise_at_critical(v2, Mh, tz, Mc, kw, rate):
    """Run once at the critical seed of each tree and collect the final state and the numerical diagnostics."""
    kwr = dict(kw)
    if kw["inflow"] == "rate":
        kwr["halo_growth_rate"] = rate
    o = v2.run_reservoir_v2(Mh, tz, np.where(np.isfinite(Mc), Mc, 1.0), **kwr)
    live = (Mh[:, :-1] > 0)
    def frac(flag):
        return (flag[:, 1:] & live).sum(axis=1) / np.maximum(live.sum(axis=1), 1)
    return dict(
        Mstar=o["stars_mass"][:, -1].tolist(), MBH=o["bh_mass"][:, -1].tolist(),
        G=(o["bh_mass"][:, -1] / np.where(np.isfinite(Mc), Mc, np.nan)).tolist(),
        nuc_star_fraction=(o["stars_nuc"][:, -1] / np.maximum(o["stars_mass"][:, -1], 1e-300)).tolist(),
        budget_ratio=(o["cum_in"] / (0.156 * np.maximum(Mh[:, -1], 1e-300))).tolist(),     # cum inflow / (f_b M_halo(z=5)), as in the frozen JSON
        halo_budget_ratio=np.asarray(o["budget_ratio"]).tolist(),                         # cum inflow / (f_b * halo mass gained), NaN if rate scheme off
        baryon_residual_max=float(np.max(np.abs(o["baryon_residual"]))),
        frac_steps_limiter_nuc=frac(o["limiter_nuc"]).tolist(), frac_steps_limiter_gal=frac(o["limiter_gal"]).tolist(),
        frac_steps_wind_cap=frac(o["wind_cap"]).tolist(),
        seed_over_host_baryons_at_formation=o["seed_over_host_baryons_at_formation"].tolist(),
        z_first_resolved=o["z_first_resolved"].tolist(),
    ), o


def run_variant(v2, name, kw, hm, zz, args, fixed_seeds):
    from ashvini.paper_reservoir import grumpy_halo_growth_rate
    cosmo = v2.PaperCosmology() if kw["cosmology"] == "paper" else None
    kw = dict(kw, cosmology=cosmo)
    tz = v2.uniform_time_grid(25.0, 5.0, args.n_steps, cosmo)
    Mh = v2.interpolate_tree_onto_grid_v2(hm, zz, tz, cosmo)
    rate = grumpy_halo_growth_rate(hm, zz, tz) if kw["inflow"] == "rate" else None
    ck = dict(kw)
    if rate is not None:
        ck["halo_growth_rate"] = rate
    t0 = time.time()
    Mc, never, above, F = v2.critical_seed_v2(Mh, tz, n_iter=40, return_F=True, **ck)
    res = dict(variant=name, switches=dict(dm_nuclear=kw["dm_nuclear"], inflow=kw["inflow"], cosmology="paper" if cosmo else "planck18_tables", n_sub=kw.get("n_sub", 1)),
               Mcrit=Mc.tolist(), bracket_fail=int(never.sum() + above.sum()), max_abs_F=float(np.nanmax(np.abs(F))))
    summ, o = summarise_at_critical(v2, Mh, tz, Mc, kw, rate)
    res.update(summ)
    if name == FINAL_VARIANT and args.crossing_grid > 0:
        rep = v2.single_crossing_report(Mh, tz, n_grid=args.crossing_grid, **ck)
        res["single_crossing"] = dict(n_grid=args.crossing_grid, n_sign_changes=rep["n_sign_changes"].tolist(),
                                      trees_with_n_changes_ne_1=int((rep["n_sign_changes"] != 1).sum()),
                                      max_decrease_above_root_median=float(np.median(rep["max_decrease"])))
        ratio = o["bh_to_host_baryons"]
        with np.errstate(all="ignore"):
            good = np.isfinite(ratio).sum(axis=0) >= 20
            pct = np.where(good[None, :], np.nanpercentile(ratio, [16, 50, 84], axis=0), np.nan)
        med = np.nanmedian(Mc)
        rep_i = int(np.nanargmin(np.abs(Mc - med)))
        res.update(z=tz.tolist(), bh_to_host_p16=pct[0].tolist(), bh_to_host_p50=pct[1].tolist(), bh_to_host_p84=pct[2].tolist(),
                   rep_index=rep_i, rep_Mcrit=float(Mc[rep_i]), rep_halo=Mh[rep_i].tolist(), rep_bh=o["bh_mass"][rep_i].tolist())
    if fixed_seeds and name in FIXED_SEED_VARIANTS:
        res["fixed_seed"] = {}
        for cname, ckw in CONFIGS.items():
            res["fixed_seed"][cname] = {}
            for s in FIXED_SEEDS:
                kk = dict(kw, **ckw)
                if rate is not None:
                    kk["halo_growth_rate"] = rate
                r = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], s), **kk)
                res["fixed_seed"][cname][str(s)] = dict(ratio=(r["bh_mass"][:, -1] / r["stars_mass"][:, -1]).tolist(),
                                                        G=(r["bh_mass"][:, -1] / s).tolist(), Mstar=r["stars_mass"][:, -1].tolist(),
                                                        nuc_star_fraction=(r["stars_nuc"][:, -1] / np.maximum(r["stars_mass"][:, -1], 1e-300)).tolist(),
                                                        delivered_fraction=(r["transfer_nuc"].sum(axis=1) / np.maximum(r["cum_in"], 1e-300)).tolist())
    res["seconds"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass-index", type=int, required=True, help="0..12, M = 3e10 * 10**(0.25 k) Msun")
    ap.add_argument("--n-trees", type=int, default=240)
    ap.add_argument("--n-steps", type=int, default=801)
    ap.add_argument("--dz", type=float, default=0.05, help="zh: tree step; pch08: spacing of the output grid only (steps are adaptive)")
    ap.add_argument("--tree-algorithm", choices=["pch08", "zh"], default="pch08",
                    help="pch08: compliant adaptive Parkinson, Cole & Helly (2008) main branch (default, as in Stochastic_Occupation); "
                         "zh: Zhang & Hui fixed step at --dz (bracketing comparison only)")
    ap.add_argument("--target-nupper", type=float, default=0.1, help="PCH08 expected splits per step (PCH08 section 2.1)")
    ap.add_argument("--m-res", type=float, default=1e4)
    ap.add_argument("--crossing-grid", type=int, default=33)
    ap.add_argument("--no-fixed-seeds", action="store_true")
    ap.add_argument("--trees-from", default=None, help="reuse the trees stored in this npz instead of sampling")
    ap.add_argument("--variants", default=",".join(VARIANTS), help="comma-separated subset of the variants")
    ap.add_argument("--final-n-sub", type=int, default=8, help="sub-steps of the final variant")
    ap.add_argument("--outdir", default=str(HERE / "output" / "mvm_v2"))
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    idx = args.mass_index
    M0 = MASSES[idx]

    from ashvini import reservoir_stock_v2 as v2
    t0 = time.time()
    if args.trees_from:
        d = np.load(args.trees_from)
        hm, zz = d["hm"], d["zz"]
        print(f"trees read from {args.trees_from}", flush=True)
    else:
        hm, zz, tree_stats = build_trees(M0, idx, args)
        extra = {} if tree_stats is None else dict(tree_stats=tree_stats)
        np.savez_compressed(out / f"trees_M{idx:02d}.npz", hm=hm, zz=zz, M0=M0, dz=args.dz, m_res=args.m_res,
                            tree_algorithm=args.tree_algorithm, target_nupper=args.target_nupper, **extra)
        print(f"trees built in {time.time()-t0:.0f}s and stored", flush=True)

    VARIANTS[FINAL_VARIANT]["n_sub"] = args.final_n_sub
    results = {}
    for name in [v.strip() for v in args.variants.split(",")]:
        results[name] = run_variant(v2, name, VARIANTS[name], hm, zz, args, fixed_seeds=not args.no_fixed_seeds)
        mc = np.array(results[name]["Mcrit"], float)
        print(f"M0={M0:.3e} {name}: median M_seed,crit {np.nanmedian(mc):.4e}  bracket_fail={results[name]['bracket_fail']}"
              f"  baryon_resid={results[name]['baryon_residual_max']:.1e}  [{results[name]['seconds']:.0f}s]", flush=True)

    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    from foraois_paths import FORAOIS_ROOT
    import ashvini.reservoir_stock as frozen_mod
    meta = dict(date=datetime.datetime.now().isoformat(timespec="seconds"), M0=M0, mass_index=idx, n_trees=int(hm.shape[0]),
                n_steps=args.n_steps, dz=args.dz, m_res=args.m_res, tree_algorithm=args.tree_algorithm,
                target_nupper=args.target_nupper, z_seed=25.0, z_anchor=5.0, f_bh=0.5,
                frozen_module_sha256=sha(frozen_mod.__file__), v2_module_sha256=sha(v2.__file__), script_sha256=sha(__file__),
                ashvini_git=git_state(HERE.parents[1]), foraois_git=git_state(FORAOIS_ROOT), variants={k: VARIANTS[k] for k in results},
                configs=CONFIGS, fixed_seeds=FIXED_SEEDS, trees_file=None if args.trees_from else f"trees_M{idx:02d}.npz",
                total_seconds=time.time() - t0)
    with open(out / f"mvm_v2_M{idx:02d}.json", "w") as f:
        json.dump(dict(meta=meta, results=results), f)
    print(f"DONE M0={M0:.3e} in {time.time()-t0:.0f}s -> {out / f'mvm_v2_M{idx:02d}.json'}", flush=True)


if __name__ == "__main__":
    main()
