#!/usr/bin/env python3
"""
How does the Minimal Viable Model (ashvini.reservoir_stock, the high-z black-hole model) change the
z=0 stellar-to-halo mass relation relative to Ashvini's general model (main.run_forest)?

Both models are run on the SAME Zhang-Hui CDM trees (z0 = 0, z_max = 30). Three variants per halo mass:
  general    main.run_forest with the tree's (smooth + merger)/dt growth rate
  mvm        run_reservoir_stock, same growth rate, identical trees -> isolates the physics difference
  mvm+grumpy run_reservoir_stock with the GRUMPY-smoothed rate its production ensemble uses

The MVM was built for z_seed = 25 to z = 5. Running it to z = 0 is an extrapolation: M_hot (cold-to-hot
transition, 4e11 Msun), R_nuc, sigma_lnj and the UV-suppression form are held at their frozen fiducial
values at all redshifts, there is no metallicity dependence of the stellar outflow, no recycling and no
reionisation-era gas heating beyond zeta_UV. The black hole is a 1e3 Msun seed at the first grid step;
it hardly grows and its AGN wind is negligible at this seed (checked via bh_final/seed in the output).
"""
import argparse
import json
import os
import time

import numpy as np

from ashvini import main, pymctrees_adapter as pa, reservoir_stock as mvm
from ashvini.paper_reservoir import grumpy_halo_growth_rate

CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 14, 13)
Z0, Z_MAX, DZ = 0.0, 30.0, 0.05
M_RES_FRACTION = 1e-4
SEED_BH = 1e3
F_B = 0.1573


def build_trees(n_halos, dz, m_res_frac, config):
    from foraois import CosmoData, ZhangHuiMergerTree
    from foraois.utils import io
    rp = io.get_params(config)
    cosmo = CosmoData(rp, redshift=[Z0])
    gen = ZhangHuiMergerTree(cosmo, rp, model=rp["Code"].get("dm_model", "cdm"))
    return rp["Cosmology"]["h"], gen


def run_bin(gen, h, m0, n_halos, dz, m_res_frac, seed=42):
    hm, z, smooth, merger = pa.build_forest_for_bin(gen, m0, h, n_halos, Z0, Z_MAX, m_res_frac * m0, dz, "numba", seed)
    rate = pa.compute_growth_rates(smooth, merger, z)
    out = {}
    t = time.time()
    g = main.run_forest(hm, rate, z)
    out["general"] = g["stars_mass"][:, -1]
    out["t_general"] = time.time() - t
    t = time.time()
    r = mvm.run_reservoir_stock(hm, z, SEED_BH, halo_growth_rate=rate)
    out["mvm"] = r["stars_mass"][:, -1]
    out["mvm_bh_growth"] = r["bh_mass"][:, -1] / SEED_BH
    out["mvm_gas"] = r["gas_mass"][:, -1]
    out["mvm_unavail"] = r["gas_unavailable"][:, -1]
    out["mvm_nuc_share"] = r["stars_nuc"][:, -1] / np.maximum(r["stars_mass"][:, -1], 1.0)
    out["mvm_resid"] = float(np.max(np.abs(r["baryon_residual"])))
    out["mvm_wind_capped_frac"] = float(r["wind_cap"].mean())
    out["t_mvm"] = time.time() - t
    try:
        gr = grumpy_halo_growth_rate(hm, z, z)
        out["mvm_grumpy"] = mvm.run_reservoir_stock(hm, z, SEED_BH, halo_growth_rate=gr)["stars_mass"][:, -1]
    except Exception as exc:  # keep going; recorded below
        out["mvm_grumpy_error"] = repr(exc)
    out["halo_final"] = hm[:, -1]
    return out


def summarise(m0, o):
    ok = o["halo_final"] > 0
    row = {"M0": m0, "n_resolved": int(ok.sum())}
    for k in ("general", "mvm", "mvm_grumpy"):
        if k in o:
            q = np.percentile(o[k][ok], [16, 50, 84])
            row[k] = {"p16": q[0], "med": q[1], "p84": q[2], "frac_zero": float(np.mean(o[k][ok] <= 0))}
    row["mvm_bh_growth_med"] = float(np.median(o["mvm_bh_growth"][ok]))
    row["mvm_gas_over_fbMh_med"] = float(np.median(o["mvm_gas"][ok] / (F_B * o["halo_final"][ok])))
    row["mvm_unavail_over_fbMh_med"] = float(np.median(o["mvm_unavail"][ok] / (F_B * o["halo_final"][ok])))
    row["mvm_nuc_share_med"] = float(np.median(o["mvm_nuc_share"][ok]))
    row["mvm_baryon_resid_max"] = o["mvm_resid"]
    row["mvm_wind_capped_step_frac"] = o["mvm_wind_capped_frac"]
    if "mvm_grumpy_error" in o:
        row["grumpy_error"] = o["mvm_grumpy_error"]
    return row


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-halos", type=int, default=200)
    ap.add_argument("--dz", type=float, default=DZ)
    ap.add_argument("--config", default=CONFIG)
    ap.add_argument("--tag", default="")
    ap.add_argument("--masses", type=float, nargs="*", default=list(MASS_BINS))
    ap.add_argument("--out", default="scripts/output")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    h, gen = build_trees(a.n_halos, a.dz, M_RES_FRACTION, a.config)
    rows = []
    for m0 in a.masses:
        try:
            rows.append(summarise(m0, run_bin(gen, h, m0, a.n_halos, a.dz, M_RES_FRACTION)))
        except Exception as exc:
            rows.append({"M0": m0, "error": repr(exc)})
        print(json.dumps(rows[-1], default=float), flush=True)
        with open(os.path.join(a.out, f"shmr_mvm_vs_general{a.tag}.json"), "w") as f:
            json.dump(rows, f, indent=1, default=float)


if __name__ == "__main__":
    main_cli()
