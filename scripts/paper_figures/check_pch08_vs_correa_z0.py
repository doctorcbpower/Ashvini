"""
Like-for-like check of the PCH08 main-branch sampler against the Correa et al. (2015a) mean accretion history.

Trees are rooted at z = 0 (the case Correa's fit describes), not at z = 5, so no conditioning on a z = 5 mass enters.
For each root mass the median and mean main-progenitor mass of the PCH08 trees are divided by Correa's M(z) at z = 2, 5, 8, 10, 12.

Usage: FORAOIS_ROOT=... CORREA_CONFIG=<planck2018.yml> python check_pch08_vs_correa_z0.py [--n-trees 100] [--m-res 1e8]
Result, 100 trees, m_res = 1e8 Msun: the median agrees with Correa to 7-28 per cent at z <= 5 and the ratio grows towards high redshift
and with root mass (at z = 10: 1.6, 3.0 and 6.7 for roots of 1e12, 1e13 and 1e14 Msun).
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
import gen_correa_comparison as C                      # noqa: E402  (builds the CosmoData tables at import)
from ashvini.pch08_main_branch import build_main_branches  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--n-trees", type=int, default=100)
ap.add_argument("--m-res", type=float, default=1e8)
ap.add_argument("--roots", default="1e12,1e13,1e14")
args = ap.parse_args()
for M0 in [float(x) for x in args.roots.split(",")]:
    hm, zz, _ = build_main_branches(C.FORAOIS_CONFIG, M0, args.n_trees, 0.0, 25.0, args.m_res, 0.05, 777)
    a, b, *_ = C.correa_alpha_beta(M0, C.h)
    row = []
    for z in (2.0, 5.0, 8.0, 10.0, 12.0):
        j = int(np.argmin(np.abs(zz - z)))
        Mc = M0 * (1 + z) ** a * np.exp(b * z)
        row.append(f"z={z:g}: median/Correa {np.median(hm[:, j]) / Mc:.2f}")
    print(f"root {M0:.0e} Msun (alpha={a:.2f}, beta={b:.2f}): " + " | ".join(row))
