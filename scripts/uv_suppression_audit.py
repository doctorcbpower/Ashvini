#!/usr/bin/env python3
"""
Audit of the UV-suppression term in the general model (main.run_forest via gas_evolve.gas_inflow_rate).

Variants, all on identical Zhang-Hui CDM trees (z0=0, z_max=30, M_res=1e-4 M0):
  current      reionization.uv_suppression as shipped: s(mu)[(1+X) - 2 eps_code M X (1+z) H / Mdot], z <= 10
  true_eps     same form with eps = -dln M_c/dz evaluated exactly (the shipped eps_code differs for z >~ 5.5)
  mass_only    s(mu, omega) only, z <= 10 (the form the MVM uses; paper_reservoir.uv_suppression_mass_only)
  off          no suppression

Derivation: with M_gas,eq = f_b s(M/M_c) M, d M_gas,eq/dt = f_b s Mdot [(1+X) - X M (dln M_c/dt)/Mdot], and
dln M_c/dt = eps_true (1+z) H, so the correction term should read X M eps_true (1+z) H / Mdot.
Usage: python scripts/uv_suppression_audit.py [--n-halos 100]
"""
import argparse
import json
import sys

import numpy as np
from scipy.special import expit

sys.path.insert(0, "scripts")
from ashvini import main, pymctrees_adapter as pa, reionization as R  # noqa: E402
import shmr_mvm_vs_general as S  # noqa: E402
from ashvini.paper_reservoir import uv_suppression_mass_only  # noqa: E402
from stellar_to_halo_mass_relation import behroozi2013_mstar, moster2013_mstar  # noqa: E402

_ORIG = R.uv_suppression


def eps_true(z):
    u = (z / R.beta) ** R.gamma
    A = R.gamma * np.where(z > 0, z, 1.0) ** (R.gamma - 1) / R.beta ** R.gamma
    return 0.63 + A * expit(u)


def uv_true_eps(z_val, m_halo, mdot_halo):
    z_val, m_halo, mdot_halo = (np.asarray(a, dtype=float) for a in (z_val, m_halo, mdot_halo))
    out = np.ones_like(z_val)
    mask = z_val <= 10
    z, m, md = z_val[mask], m_halo[mask], mdot_halo[mask]
    ok = (m > 0) & (md > 0)
    ms, mds = np.where(ok, m, 1.0), np.where(ok, md, 1.0)
    ratio = eps_true(z) * ms * R.X(z, ms) * (1 + z) * R.cosmo.H(z).value / mds
    val = R.s(R.mu_c(z, ms), R.omega) * ((1 + R.X(z, ms)) - ratio)
    out[mask] = np.maximum(np.where(ok, val, 0.0), 0.0)
    return out


def uv_off(z_val, m_halo, mdot_halo):
    return np.ones_like(np.asarray(z_val, dtype=float))


def uv_mass_only(z_val, m_halo, mdot_halo):
    return uv_suppression_mass_only(np.asarray(z_val, dtype=float), np.asarray(m_halo, dtype=float))


VARIANTS = {"current": _ORIG, "true_eps": uv_true_eps, "mass_only": uv_mass_only, "off": uv_off}


def med_mstar(hm, rate, z, name):
    R.uv_suppression = VARIANTS[name]
    try:
        ok = hm[:, -1] > 0
        return float(np.median(main.run_forest(hm, rate, z)["stars_mass"][ok, -1]))
    finally:
        R.uv_suppression = _ORIG


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-halos", type=int, default=100)
    ap.add_argument("--config", default=S.CONFIG)
    ap.add_argument("--tag", default="_cdm")
    a = ap.parse_args()
    h, gen = S.build_trees(a.n_halos, 0.005, 1e-4, a.config)
    masses = [1e9, 3.16e9, 1e10, 3.16e10, 1e11, 3.16e11, 1e12, 3.16e12, 1e13]
    table = []
    for m0 in masses:
        hm, z, sm, mg = pa.build_forest_for_bin(gen, m0, h, a.n_halos, 0.0, 30.0, 1e-4 * m0, 0.005, "numba", 42)
        rate = pa.compute_growth_rates(sm, mg, z)
        row = {"M0": m0, "behroozi": float(behroozi2013_mstar(np.array([m0]))[0]),
               "moster": float(moster2013_mstar(np.array([m0]))[0])}
        for v in VARIANTS:
            row[v] = med_mstar(hm, rate, z, v)
        # zero-inflow fraction of the current term for resolved halos at z<=5
        R.uv_suppression = _ORIG
        zz = np.broadcast_to(z, hm.shape)
        sup = _ORIG(zz, hm, rate)
        sel = (hm > 0) & (zz <= 5) & (rate > 0)
        row["frac_steps_zero_S_current_z<=5"] = float(np.mean(sup[sel] == 0)) if sel.any() else None
        table.append(row)
        print(json.dumps(row), flush=True)
    dzt = {}
    for m0 in (1e10, 1e11):
        dzt[str(m0)] = {}
        for dz in (0.02, 0.01, 0.005):
            hm, z, sm, mg = pa.build_forest_for_bin(gen, m0, h, 60, 0.0, 30.0, 1e-4 * m0, dz, "numba", 42)
            rate = pa.compute_growth_rates(sm, mg, z)
            dzt[str(m0)][str(dz)] = {v: med_mstar(hm, rate, z, v) for v in VARIANTS}
        print(m0, json.dumps(dzt[str(m0)]), flush=True)
    json.dump({"table": table, "dz_test": dzt}, open(f"scripts/output/uv_suppression_audit{a.tag}.json", "w"), indent=1)


if __name__ == "__main__":
    main_cli()
