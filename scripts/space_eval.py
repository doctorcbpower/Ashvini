"""
Evaluate a shard of the parameter-space samples (scripts/space_common.py): for each sample run the grid integrator on the stored trees of all nine
halo masses and store the summary fields. Resume-safe: progress is saved every few samples.

Usage: python scripts/space_eval.py <subspace> <shard> <n_shards> [--max-samples K]
Sample s of a sub-space belongs to shard s % n_shards. Output: $SHMR_OUT/space/<subspace>_shard<shard>.npz
"""
import os
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C
import space_common as S

warnings.filterwarnings("ignore")
SAVE_EVERY = 5


def evaluate(cfg, trees):
    from ashvini.main import run_forest
    out = np.full((len(trees), len(S.FIELDS)), np.nan)
    for k, t in enumerate(trees):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = run_forest(t["halo_masses"], t["halo_growth_rates"], t["redshifts"], config=cfg)
        s = r["stars_mass"][:, -1]
        g = r["gas_mass"][:, -1]
        gz = np.where(g > 0, r["gas_metals"][:, -1] / np.where(g > 0, g, 1.0), np.nan)
        sz = np.where(s > 0, r["stars_metals"][:, -1] / np.where(s > 0, s, 1.0), np.nan)
        bad = sum(int(np.sum(~np.isfinite(r[key]))) for key in ("gas_mass", "stars_mass", "gas_metals"))
        out[k] = [np.median(s), np.percentile(s, 16), np.percentile(s, 84), np.median(g), np.nanmedian(gz) if np.any(np.isfinite(gz)) else np.nan,
                  np.nanmedian(sz) if np.any(np.isfinite(sz)) else np.nan, np.median(r["bh_mass"][:, -1]), np.mean(s > 0), bad]
    return out


def run_shard(subspace, shard, n_shards, max_samples=None):
    from ashvini.config import DEFAULT, override
    path = os.path.join(C.OUT, "space", f"{subspace}_shard{shard}.npz")
    unit = S.design(subspace)
    mine = np.arange(shard, len(unit), n_shards)
    if max_samples:
        mine = mine[:max_samples]
    results = np.full((len(mine), len(C.MASSES), len(S.FIELDS)), np.nan)
    done = np.zeros(len(mine), dtype=bool)
    if os.path.exists(path):
        d = np.load(path)
        if np.array_equal(d["samples"], mine):
            results, done = d["results"], d["done"]
    trees = [dict(np.load(C.tree_path(i))) for i in range(len(C.MASSES))]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    t0 = time.time()
    for n, s in enumerate(mine):
        if done[n]:
            continue
        cfg = override(DEFAULT, **S.config_changes(subspace, S.physical(subspace, unit[s])))
        results[n] = evaluate(cfg, trees)
        done[n] = True
        if (n + 1) % SAVE_EVERY == 0 or n == len(mine) - 1:
            tmp = path + ".tmp.npz"
            np.savez(tmp, samples=mine, results=results, done=done)
            os.replace(tmp, path)
            print(f"{subspace} shard {shard}: {done.sum()}/{len(mine)} samples, {time.time() - t0:.0f} s", flush=True)
    return path


if __name__ == "__main__":
    sub, shard, n_shards = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    k = int(sys.argv[sys.argv.index("--max-samples") + 1]) if "--max-samples" in sys.argv else None
    run_shard(sub, shard, n_shards, k)
