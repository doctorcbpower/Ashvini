#!/usr/bin/env python3
"""
Does the M_res dependence of the z=0 stellar mass belong to the physics or to
the tree algorithm?

Runs Ashvini's default baryonic model on CDM trees from both foraois
generators (PCH08 and Zhang-Hui), at fixed halo mass, over a range of
M_res/M0. Motivated by foraois docs/PCH08_HIGH_Z_DIAGNOSTIC.md: PCH08
main progenitors at small M_res/M0 are near-deterministic and assemble
earlier than Zhang-Hui, so a trend in M_star with M_res that appears with
PCH08 but not with Zhang-Hui would be a tree-algorithm artefact.

Reports, per (algorithm, M0, M_res/M0): median and 16-84 per cent of
M_star(z=0), the median M(z=1)/M0, and the baryon-budget check
M_star + M_gas + M_BH <= f_b M_halo.

Usage: python scripts/zh_vs_pch08_mres_scan.py [--n-halos 200] [--out DIR]
"""
import argparse
import json
import os
import time

import numpy as np

from ashvini import main, pymctrees_adapter as pa
from ashvini.utils import time_at_z  # noqa: F401  (imported to fail early if cosmology is broken)

CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASSES = [1e9, 1e10, 1e11, 1e12]
MRES_FRACTIONS = [1e-2, 1e-3, 1e-4, 1e-5]
Z0, Z_MAX, DZ = 0.0, 30.0, 0.01
SEED = 42
F_B = 0.1573  # Omega_b / Omega_m, Planck 2018


def make_generators(config=CONFIG, algos=("PCH08", "ZH")):
    from foraois import CosmoData, PCHMergerTree, ZhangHuiMergerTree
    from foraois.utils import io

    run_params = io.get_params(config)
    dm_model = run_params["Code"].get("dm_model", "cdm")
    h = run_params["Cosmology"]["h"]
    cosmo = CosmoData(run_params, redshift=[Z0])
    all_gens = {
        "PCH08": lambda: PCHMergerTree(cosmo, run_params),
        "ZH": lambda: ZhangHuiMergerTree(cosmo, run_params, model=dm_model),
    }
    return h, {a: all_gens[a]() for a in algos}


def one_case(gen, h, m0, frac, n_halos):
    t0 = time.time()
    halo_masses, redshifts, smooth, merger = pa.build_forest_for_bin(
        gen, m0, h, n_halos, Z0, Z_MAX, frac * m0, DZ, "numba", SEED)
    rates = pa.compute_growth_rates(smooth, merger, redshifts)
    res = main.run_forest(halo_masses, rates, redshifts)
    t1 = time.time()

    m_star = res["stars_mass"][:, -1]
    m_gas = res["gas_mass"][:, -1] if "gas_mass" in res else np.zeros_like(m_star)
    m_bh = res["bh_mass"][:, -1] if "bh_mass" in res else np.zeros_like(m_star)
    m_halo = halo_masses[:, -1]
    ok = m_halo > 0
    i1 = int(np.argmin(np.abs(redshifts - 1.0)))
    frac_z1 = halo_masses[:, i1] / m_halo
    budget = (m_star + m_gas + m_bh)[ok] / (F_B * m_halo[ok])
    q = np.percentile(m_star[ok], [16, 50, 84])
    return {
        "M0": m0, "m_res_frac": frac, "n_resolved": int(ok.sum()), "n": n_halos,
        "mstar_p16": q[0], "mstar_med": q[1], "mstar_p84": q[2],
        "M_z1_over_M0_med": float(np.median(frac_z1[ok])),
        "budget_max": float(budget.max()), "seconds": t1 - t0,
        "keys": sorted(res.keys()),
    }


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-halos", type=int, default=200)
    ap.add_argument("--out", default="scripts/output")
    ap.add_argument("--config", default=CONFIG, help="foraois YAML config (sets the DM model)")
    ap.add_argument("--algos", nargs="*", default=["PCH08", "ZH"], choices=["PCH08", "ZH"])
    ap.add_argument("--tag", default="", help="suffix for the output JSON name")
    ap.add_argument("--masses", type=float, nargs="*", default=MASSES)
    ap.add_argument("--fracs", type=float, nargs="*", default=MRES_FRACTIONS)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    h, gens = make_generators(args.config, tuple(args.algos))
    rows = []
    for algo, gen in gens.items():
        for m0 in args.masses:
            for frac in args.fracs:
                try:
                    r = one_case(gen, h, m0, frac, args.n_halos)
                except Exception as exc:  # keep scanning, record the failure
                    r = {"M0": m0, "m_res_frac": frac, "error": repr(exc)}
                r["algo"] = algo
                rows.append(r)
                print(json.dumps({k: v for k, v in r.items() if k != "keys"}, default=float), flush=True)
                with open(os.path.join(args.out, f"zh_vs_pch08_mres_scan{args.tag}.json"), "w") as f:
                    json.dump(rows, f, indent=1, default=float)


if __name__ == "__main__":
    main_cli()
