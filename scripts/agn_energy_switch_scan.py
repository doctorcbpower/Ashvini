#!/usr/bin/env python3
"""Test ladder for the illustrative momentum-to-energy AGN switch (scripts/agn_energy_switch.py) on top of the
tuned general model (mass-only UV, e_ff = 0.05). Same Zhang-Hui CDM trees for every variant."""
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
import agn_energy_switch as E  # noqa: E402
import shmr_mvm_vs_general as S  # noqa: E402
import uv_suppression_audit as A  # noqa: E402
from ashvini import black_holes_growth as bg, main, pymctrees_adapter as pa, reionization as R  # noqa: E402
from stellar_to_halo_mass_relation import behroozi2013_mstar  # noqa: E402

N, DZ = 100, 0.005
MASSES = [1e10, 3.16e10, 1e11, 3.16e11, 1e12, 3.16e12, 1e13, 3.16e13, 1e14]
R.uv_suppression = A.uv_mass_only
main.e_ff = 0.05
h, gen = S.build_trees(N, DZ, 1e-4, S.CONFIG)
trees = {m: pa.build_forest_for_bin(gen, m, h, N, 0.0, 30.0, 1e-4 * m, DZ, "numba", 42) for m in MASSES}
beh = np.array([float(behroozi2013_mstar(np.array([m]))[0]) for m in MASSES])

VARIANTS = {
    "shipped (eta_agn=0.5, no cap)": None,
    "AGN wind off (f_mom=0, no energy, no cap)": dict(f_mom=0.0, energy=False, cap=False),
    "cap only (shipped wind)": "cap",
    "momentum only + cap": dict(energy=False, cap=True),
    "full: momentum->energy + cap (eps_f=5e-4)": dict(epsilon_f=5e-4),
    "full, eps_f=1e-5": dict(epsilon_f=1e-5),
    "full, eps_f=1e-4": dict(epsilon_f=1e-4),
    "full, eps_f=2e-3": dict(epsilon_f=2e-3),
    "full, no cap (eps_f=5e-4)": dict(epsilon_f=5e-4, cap=False),
}


def run(spec, m):
    hm, z, sm, mg = trees[m]
    ok = hm[:, -1] > 0
    rate = pa.compute_growth_rates(sm, mg, z)
    if spec is None:
        r = main.run_forest(hm, rate, z)
    elif spec == "cap":
        bg.growth_cap_enabled = True
        try:
            r = main.run_forest(hm, rate, z)
        finally:
            bg.growth_cap_enabled = False
    else:
        with E.installed(**spec):
            r = main.run_forest(hm, rate, z)
    ms = bg.m_sigma(bg.velocity_dispersion(hm[ok, -1], 0.0))
    tot = r["stars_mass"] + r["gas_mass"] + r["bh_mass"]
    avail = 0.1573 * np.maximum.accumulate(hm, axis=1)
    budget = float(np.max(np.where(avail > 0, tot / np.where(avail > 0, avail, 1.0), 0.0)))
    return (float(np.median(r["stars_mass"][ok, -1])), float(np.median(r["bh_mass"][ok, -1])),
            float(np.median(r["bh_mass"][ok, -1] / ms)), budget)


out = {"masses": MASSES, "behroozi": beh.tolist(), "runs": {}}
print("log10(M*/Behroozi) | BH/M_sigma   at M_halo =", " ".join(f"{m:.0e}" for m in MASSES))
for name, spec in VARIANTS.items():
    res = [run(spec, m) for m in MASSES]
    out["runs"][name] = res
    lr = np.log10(np.array([r[0] for r in res]) / beh)
    print(f"{name:44s}", " ".join(f"{x:+.2f}" for x in lr), "|", " ".join(f"{r[2]:.1f}" for r in res), "| max budget ratio", f"{max(r[3] for r in res):.2f}", flush=True)
json.dump(out, open("scripts/output/agn_energy_switch_scan.json", "w"), indent=1)
