"""
C19: effect of the tree-builder timestep compliance on the critical seed. Frozen MVM (hash asserted; unmodified).

The production trees (dz = 0.05, M_res = 1e4 Msun) are outside the practical single-split compliance regime of the Zhang & Hui
builder (expected splits per step E of about 0.1 or less; foraois docs/PCH08_HIGH_Z_DIAGNOSTIC.md). The growth-rate estimator
(GRUMPY spline) uses the tree's own node spacing, so shrinking dz for the whole pipeline confounds tree compliance with the rate
estimator's resolution. Here the trees are built with dz = 0.05 / k and recorded ONLY at the production dz = 0.05 checkpoints
(401 nodes, z = 25 to 5), so the interpolation onto the 801-node reservoir grid and the rate estimator see the same node spacing
as production; what changes is how well the tree builder resolves the splits between checkpoints. Everything else is the
production pipeline: M_res = 1e4 Msun, z = 25 to 5, seed at the first node, f_BH = 0.5, fiducial parameters.

E is the maximum over z = 5 to 25 of foraois diagnostics.expected_eps_splits_per_step (evaluated at the anchored mass M0, which is
conservative at high z), sampled at up to 40 redshifts. E scales as 1/dz. Full compliance (E <= 0.1) is reachable at 3e10 only;
at 3e11 the smallest dz used gives E of about 0.3, and at 3e13 E stays far above 0.1.

Independent tree sets (no seed reproducibility: numba sampler); the statistic is the mean over sets of the per-set median of
M_seed,crit, as a ratio to the k = 1 (production dz) mean, with the standard error over sets. Requires a foraois checkout that
contains expected_eps_splits_per_step (the commit that added it, or later); FORAOIS_ROOT must point to it.

Usage: python c19_timestep_compliance.py   (about 10 minutes on 8 cores; prints the log to stdout, writes
output/c19_timestep_compliance_results.json)
"""
import sys, json, time, subprocess, platform, datetime
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, generators, grid
from foraois_paths import FORAOIS_ROOT
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "output" / "c19_timestep_compliance_results.json"
DZ0, M_RES, Z_MAX = 0.05, 1e4, 25.0
# mass -> (k values, number of tree sets, trees per set)
PLAN = {
    3e10: ((1, 10, 100, 1000), 6, 100),
    3e11: ((1, 10, 100, 1000, 2500), 6, 60),
    3e13: ((1, 10, 100, 1000), 6, 60),
}


def checkpoint_indices(n_nodes, k):
    """chronological node indices that lie on the production dz = 0.05 grid; the last node is z0 = 5 and is always kept"""
    return np.arange((n_nodes - 1) % k, n_nodes, k)


def worker(a):
    M0, k, s, N, seed = a
    assert_frozen()
    from ashvini import pymctrees_adapter
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    h, zh, _ = generators()
    t0 = time.time()
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(zh, M0, h, N, z0=5.0, z_max=Z_MAX, m_res_msun=M_RES, dz=DZ0 / k,
                                                        backend="numba", rng_seed=seed)
    n_fine = len(zz)
    idx = checkpoint_indices(n_fine, k)
    hm, zz = hm[:, idx], zz[idx]
    if len(zz) != 401 or abs(zz[-1] - 5.0) > 1e-9 or abs(zz[0] - 25.0) > 1e-6:
        raise SystemExit(f"checkpoint subsample failed: {len(zz)} nodes, z from {zz[0]} to {zz[-1]} (k={k})")
    tz = grid(25.0)
    Mh = interpolate_tree_onto_grid(hm, zz, tz)
    rate = grumpy_halo_growth_rate(hm, zz, tz)
    Mc, never, above = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate)
    return dict(M0=M0, k=k, set=s, N=N, seed=seed, dz=DZ0 / k, fine_nodes=int(n_fine), Mc=[float(x) for x in Mc],
                bad=int(never.sum() + above.sum()), seconds=round(time.time() - t0, 1))


def eps_table():
    from foraois import cosmo_utils, ZhangHuiMergerTree
    from foraois.utils import io
    from foraois.diagnostics import expected_eps_splits_per_step
    from foraois_paths import FORAOIS_CONFIG
    rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
    cd = cosmo_utils.CosmoData(rp, redshift=[5.0]); ZhangHuiMergerTree(cd, rp, model="cdm")  # populates the sigma grid
    E = {}
    for M0, (ks, _, _) in PLAN.items():
        for k in ks:
            dz = DZ0 / k
            zs = np.arange(5.0, 25.0, dz)
            st = max(1, len(zs) // 40)
            E[(M0, k)] = max(expected_eps_splits_per_step(cd, M0 * h, z, z + dz, M_RES * h, model="cdm") for z in zs[::st])
    return E


def git(*a):
    try:
        return subprocess.run(["git", "-C", FORAOIS_ROOT, *a], capture_output=True, text=True, check=True).stdout.strip()
    except Exception as e:  # provenance only
        return f"unavailable ({e})"


if __name__ == "__main__":
    assert_frozen()
    import numba, scipy
    print("C19: tree-builder timestep compliance vs M_seed,crit (frozen MVM, production pipeline; trees built at dz = 0.05/k, "
          "recorded at the production dz = 0.05 checkpoints)")
    print(f"date {datetime.datetime.now().isoformat(timespec='seconds')}; foraois {git('rev-parse', 'HEAD')} "
          f"(uncommitted changes: {bool(git('status', '--porcelain'))}); python {platform.python_version()}, numpy {np.__version__}, "
          f"scipy {scipy.__version__}, numba {numba.__version__}")
    E = eps_table()
    tasks = []
    for mi, (M0, (ks, nsets, N)) in enumerate(PLAN.items()):
        for ki, k in enumerate(ks):
            for s in range(nsets):
                tasks.append((M0, k, s, N, 6000 + 1000 * mi + 100 * ki + s))
    t0 = time.time()
    with Pool(3) as p:
        res = p.map(worker, tasks, chunksize=1)
    print(f"{len(tasks)} tree sets in {time.time() - t0:.0f} s; bracket failures in total: {sum(r['bad'] for r in res)}")

    summary = {}
    for M0, (ks, nsets, N) in PLAN.items():
        med = {k: np.array([np.nanmedian(r["Mc"]) for r in res if r["M0"] == M0 and r["k"] == k]) for k in ks}
        base = med[1].mean()
        print(f"\nM0 = {M0:.0e} Msun: {nsets} independent tree sets of {N} trees per setting; median M_seed,crit per set, ratio to the k = 1 mean "
              f"({base:.3e} Msun)")
        print("     k        dz     E(max z)   mean median  sd/mean   ratio to k=1   (standard error)   compliant (E <= 0.1)")
        for k in ks:
            m = med[k]; se = m.std(ddof=1) / np.sqrt(len(m)) / base
            summary[f"{M0:.0e}_{k}"] = dict(M0=M0, k=k, dz=DZ0 / k, E=E[(M0, k)], mean_median=float(m.mean()), sd_over_mean=float(m.std(ddof=1) / m.mean()),
                                            ratio=float(m.mean() / base), ratio_se=float(se), n_sets=len(m))
            print(f"  {k:5d}  {DZ0 / k:9.1e}  {E[(M0, k)]:9.3g}   {m.mean():.3e}   {m.std(ddof=1) / m.mean():6.3f}   {m.mean() / base:9.3f}       ({se:.3f})          "
                  f"{'yes' if E[(M0, k)] <= 0.1 else 'no'}")
    print("\nnotes: E is evaluated at the anchored mass M0 (conservative at high z). Only 3e10 reaches E <= 0.1; at 3e11 the smallest dz gives E about 0.3 and at 3e13 "
          "E stays far above 0.1, so the 3e13 ratios are not a compliant-limit result. Ratios below 1 mean the compliant trees give a smaller median.")
    json.dump(dict(meta=dict(date=datetime.datetime.now().isoformat(timespec="seconds"), foraois_head=git("rev-parse", "HEAD"),
                             foraois_dirty=bool(git("status", "--porcelain")), python=platform.python_version(), numpy=np.__version__,
                             numba=numba.__version__, scipy=scipy.__version__, dz0=DZ0, M_res=M_RES, plan={f"{k:.0e}": v for k, v in PLAN.items()}),
                   summary=summary, sets=[{k: v for k, v in r.items()} for r in res]), open(OUT, "w"))
    print(f"wrote {OUT}")
