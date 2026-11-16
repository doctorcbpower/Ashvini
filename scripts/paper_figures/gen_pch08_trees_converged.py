"""
Build PCH08 main-branch trees for one halo-mass bin at a resolution set by the ratio R = M0 / m_res, and summarise how the
median main-progenitor mass depends on R.

Why a ratio: tests at M0 = 3e12 Msun showed that the median mass of the main progenitor at z >= 8 still changes between
m_res = 1e4 (R = 3e8) and 3e3 (R = 1e9) while 3e3 and 1e3 agree to about 1 per cent, so convergence is a statement about
the dynamic range M0 / m_res and not about m_res itself. This script lets one array task per (mass, R) test that for each mass.

Output (--outdir):
  trees_M<idx>_R<exp>.npz   hm, zz, M0, dz, m_res, R, tree_algorithm, target_nupper, tree_stats   (same layout as the production trees)
  summary_M<idx>_R<exp>.json   median and 16th/84th percentile of the main-progenitor mass at z = 5.5, 6, 7, 8, 10, 12, 15, 20

Usage: python gen_pch08_trees_converged.py --mass-index I --ratio 1e9 --n-trees 240 --outdir DIR
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG    # noqa: E402

MASSES = 3.0e10 * 10 ** (0.25 * np.arange(13))              # same grid as gen_mvm_v2_production.py
Z_REPORT = (5.5, 6, 7, 8, 10, 12, 15, 20)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass-index", type=int, required=True, help="0..12, M = 3e10 * 10**(0.25 k) Msun")
    ap.add_argument("--ratio", type=float, required=True, help="R = M0 / m_res (the earlier production trees used 3e8 at 3e12 Msun, which is too fine: 89 per cent of steps hit the step floor)")
    ap.add_argument("--n-trees", type=int, default=240)
    ap.add_argument("--dz", type=float, default=0.05, help="spacing of the output grid only; steps are adaptive")
    ap.add_argument("--target-nupper", type=float, default=0.1)
    ap.add_argument("--seed0", type=int, default=None, help="default 1000 * (idx + 1)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    sys.path.insert(0, FORAOIS_SRC)
    from ashvini import pch08_main_branch

    idx, M0 = a.mass_index, MASSES[a.mass_index]
    m_res = M0 / a.ratio
    seed0 = a.seed0 if a.seed0 is not None else 1000 * (idx + 1)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"M{idx:02d}_R{int(round(np.log10(a.ratio) * 10)):03d}"          # R = 10**(n/10), e.g. R030 -> 1e3, R095 -> ~3e9
    t0 = time.time()
    hm, zz, stats = pch08_main_branch.build_main_branches(
        FORAOIS_CONFIG, M0, a.n_trees, z0=5.0, z_max=25.0, m_res_msun=m_res, dz_grid=a.dz, seed0=seed0, target_nupper=a.target_nupper)
    hm, zz = np.asarray(hm), np.asarray(zz)
    secs = time.time() - t0
    np.savez_compressed(out / f"trees_{tag}.npz", hm=hm, zz=zz, M0=M0, dz=a.dz, m_res=m_res, R=a.ratio, tree_algorithm="pch08",
                        target_nupper=a.target_nupper, tree_stats=stats)
    summ = dict(mass_index=idx, M0=M0, ratio=a.ratio, m_res=m_res, n_trees=a.n_trees, seconds=secs,
                floor_steps=int(np.asarray(stats)[:, 2].sum()), z={})
    for z in Z_REPORT:
        col = hm[:, int(np.argmin(np.abs(zz - z)))]
        p = np.percentile(col[col > 0], [16, 50, 84]) if (col > 0).any() else [np.nan] * 3
        summ["z"][str(z)] = dict(p16=float(p[0]), median=float(p[1]), p84=float(p[2]), alive=int((col > 0).sum()))
    (out / f"summary_{tag}.json").write_text(json.dumps(summ, indent=1))
    print(f"{tag}: M0={M0:.2e} m_res={m_res:.2e} {a.n_trees} trees in {secs:.0f}s; floor steps {summ['floor_steps']}")
    for z in Z_REPORT:
        print(f"  z={z:>4}: median {summ['z'][str(z)]['median']:.3e}  alive {summ['z'][str(z)]['alive']}/{a.n_trees}")


if __name__ == "__main__":
    main()
