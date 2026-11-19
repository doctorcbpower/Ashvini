"""
How do the black hole results depend on the tree resolution M_res / M0? (docs/TREESAMPLE_AUDIT.md, E3 and A2)

For four tree masses (1e9 to 1e12 Msun), Zhang-Hui trees are built at M_res/M0 = 1e-3 and 1e-4 (the finest the one-split-per-step condition allows
at dz = 0.005), and at 1e-4 a second time with another seed, which measures the sampling noise of the comparison. The fiducial model is run on each
and the black hole seeding is summarised: the fraction of haloes with a black hole, which seeding channel fired, the redshift of seeding, the redshift
at which the tree first resolves the halo, and the black hole and stellar masses at z = 0.

Usage: python scripts/resolution_bh_scan.py run     (resumable; one JSON per case in $SHMR_OUT/resolution_bh)
       python scripts/resolution_bh_scan.py report
"""
import glob
import json
import os
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C

warnings.filterwarnings("ignore")
MASSES = [1e9, 1e10, 1e11, 1e12]
# (M_res / M0, seed): the first two differ in resolution, the last two in seed only
CASES = [(1e-3, 101), (1e-4, 101), (1e-4, 202)]
N_HALOES = int(os.environ.get("RES_N_HALOES", "200"))
OUTDIR = os.path.join(C.OUT, "resolution_bh")


def path(mass, ratio, seed):
    return os.path.join(OUTDIR, f"M{mass:.0e}_r{ratio:.0e}_s{seed}.json")


def summarise(mass, ratio, seed, hm, hg, z):
    from ashvini.config import DEFAULT
    from ashvini.main import run_forest
    r = run_forest(hm, hg, z, config=DEFAULT)
    bh, st = r["bh_mass"], r["stars_mass"]
    has = (bh > 0).any(axis=1)
    first = np.array([np.argmax(b > 0) if (b > 0).any() else -1 for b in bh])
    seed_mass = np.array([bh[i, first[i]] if first[i] >= 0 else 0.0 for i in range(len(bh))])
    z_seed = np.where(first >= 0, z[np.maximum(first, 0)], np.nan)
    first_res = np.array([np.argmax(h > 0) for h in hm])      # chronological: index 0 is the earliest step
    pct = lambda a: [float(x) for x in np.nanpercentile(a, [16, 50, 84])] if np.isfinite(a).any() else [np.nan] * 3
    bh0, st0 = bh[:, -1], st[:, -1]
    return dict(
        mass=mass, ratio=ratio, seed=seed, n=len(hm), m_res=ratio * mass,
        frac_with_bh=float(has.mean()),
        channel_pop3=float(np.mean(has & (seed_mass < 3e2))), channel_threshold=float(np.mean(has & (seed_mass >= 3e2) & (seed_mass < 3e4))),
        channel_direct_collapse=float(np.mean(has & (seed_mass >= 3e4))),
        z_seed_pct=pct(z_seed), z_first_resolved_pct=pct(z[first_res]),
        log_mbh_pct_occupied=pct(np.log10(bh0[bh0 > 0])) if (bh0 > 0).any() else [np.nan] * 3,
        log_mstar_pct=pct(np.log10(np.maximum(st0, 1e-30))),
        median_mbh_all=float(np.median(bh0)),
    )


def run():
    from ashvini import pymctrees_adapter
    os.makedirs(OUTDIR, exist_ok=True)
    for mass in MASSES:
        for ratio, seed in CASES:
            p = path(mass, ratio, seed)
            if os.path.exists(p):
                print(p, "exists, skipped"); continue
            t0 = time.time()
            hm, hg, z, _ = pymctrees_adapter.build_forest_live(
                C.FORAOIS_CONFIG, mass, N_HALOES, C.Z0, C.Z_MAX, C.DZ, m_res=ratio * mass, seed=seed,
                algorithm=C.ALGORITHM, pk_kmax_override=C.PK_KMAX_OVERRIDE)
            out = summarise(mass, ratio, seed, hm, hg, z)
            with open(p, "w") as f:
                json.dump(out, f)
            print(f"M0 {mass:.0e} ratio {ratio:.0e} seed {seed}: {time.time() - t0:.0f} s", flush=True)


def report():
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUTDIR, "*.json")))]
    f2 = lambda v: " ".join(f"{x:6.2f}" for x in v)
    for mass in MASSES:
        print(f"\nM0 = {mass:.0e} Msun   (16/50/84 percentiles; 'seed noise' = the two M_res/M0 = 1e-4 rows)")
        for r in [r for r in rows if r["mass"] == mass]:
            tag = f"ratio {r['ratio']:.0e} seed {r['seed']}  M_res {r['m_res']:.0e}"
            print(f"  {tag}: with BH {r['frac_with_bh']:.2f}  channels pop3/thr/dc {r['channel_pop3']:.2f}/{r['channel_threshold']:.2f}/{r['channel_direct_collapse']:.2f}"
                  f"\n      z seeded {f2(r['z_seed_pct'])}   z first resolved {f2(r['z_first_resolved_pct'])}"
                  f"\n      log M_BH(z=0) occupied {f2(r['log_mbh_pct_occupied'])}   log M* {f2(r['log_mstar_pct'])}")


if __name__ == "__main__":
    {"run": run, "report": report}[sys.argv[1] if len(sys.argv) > 1 else "run"]()
