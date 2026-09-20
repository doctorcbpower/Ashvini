#!/usr/bin/env python3
"""ILLUSTRATIVE test: add a cold-to-hot inflow suppression with a residual floor to the general model.

    zeta_ch(M) = f + (1 - f) * (1 - s(M / M_hot, phi)),   phi = 4  (cold_hot_mode_suppression, as in the MVM)
    inflow = zeta_UV * zeta_ch * f_b * Mdot_halo

f = 0 is the MVM's full cutoff; f > 0 lets a fraction of accretion continue above M_hot (hot-mode cooling).
M_hot is fixed in time (a modelling limitation, see the discussion in the session notes). Baseline is the tuned
general model (mass-only UV, e_ff = 0.05, shipped AGN wind). Same Zhang-Hui CDM trees for every variant.
"""
import json
import sys

import numpy as np

sys.path.insert(0, "scripts")
import shmr_mvm_vs_general as S  # noqa: E402
import uv_suppression_audit as A  # noqa: E402
from ashvini import main, pymctrees_adapter as pa, reionization as R  # noqa: E402
from ashvini.paper_reservoir import cold_hot_mode_suppression  # noqa: E402
from stellar_to_halo_mass_relation import behroozi2013_mstar  # noqa: E402

N, DZ, PHI = 100, 0.005, 4.0
MASSES = [1e10, 3.16e10, 1e11, 3.16e11, 1e12, 3.16e12, 1e13, 3.16e13, 1e14]
main.e_ff = 0.05
h, gen = S.build_trees(N, DZ, 1e-4, S.CONFIG)
trees = {m: pa.build_forest_for_bin(gen, m, h, N, 0.0, 30.0, 1e-4 * m, DZ, "numba", 42) for m in MASSES}
beh = np.array([float(behroozi2013_mstar(np.array([m]))[0]) for m in MASSES])


def make_uv(m_hot, floor):
    def uv(z, m, mdot):
        m = np.asarray(m, dtype=float)
        zeta = floor + (1.0 - floor) * cold_hot_mode_suppression(m, M_hot=m_hot, phi=PHI)
        return A.uv_mass_only(z, m, mdot) * zeta
    return uv


def run(m_hot, floor):
    R.uv_suppression = make_uv(m_hot, floor) if m_hot else A.uv_mass_only
    out = []
    for m in MASSES:
        hm, z, sm, mg = trees[m]
        ok = hm[:, -1] > 0
        r = main.run_forest(hm, pa.compute_growth_rates(sm, mg, z), z)
        out.append(float(np.median(r["stars_mass"][ok, -1])))
    return out


res = {"masses": MASSES, "behroozi": beh.tolist(), "runs": {}}
print("log10(M*/Behroozi) at M_halo =", " ".join(f"{m:.0e}" for m in MASSES), "| rms 1e10-1e12 | rms 1e12-1e14")
cases = [("none (baseline)", None, 0.0)]
for mh in (4e11, 1e12, 2e12):
    for f in (0.0, 0.1, 0.3):
        cases.append((f"M_hot={mh:.0e} floor={f}", mh, f))
sel1 = [i for i, m in enumerate(MASSES) if 1e10 <= m <= 1e12]
sel2 = [i for i, m in enumerate(MASSES) if 1e12 <= m <= 1e14]
for name, mh, f in cases:
    med = run(mh, f)
    lr = np.log10(np.array(med) / beh)
    res["runs"][name] = med
    print(f"{name:26s}", " ".join(f"{x:+.2f}" for x in lr),
          f"| {np.sqrt(np.mean(lr[sel1]**2)):.2f} | {np.sqrt(np.mean(lr[sel2]**2)):.2f}", flush=True)
json.dump(res, open("scripts/output/hot_mode_floor_scan.json", "w"), indent=1)
