"""
Sensitivity ledger of the revised scoping model: every parameter and modelling assumption is varied one at a time about the
fiducial model, on the same stored PCH08 trees, and the response of two headline quantities is recorded.

Quantities (per halo mass, medians over the trees):
  Mcrit   critical seed mass, the seed for which M_BH = f_BH M_* at z = 5
  G       growth factor M_BH(z=5) / M_seed of a 1e3 Msun seed (the light-seed outcome)
Each row reports the change in log10 of the median relative to the baseline, so that a row is directly comparable with the
statistical error of the baseline (bootstrap over trees) and with the numerical error (n_sub, time-step).

Rows are grouped as
  numerical    n_sub 1 and 32, number of time steps
  physical     the model's free parameters, each varied over the range stated in RANGES (low, high)
  structural   assumptions that are choices of form: nuclear dark matter in the capture threshold, inflow prescription,
               cosmology tables, target definition f_BH
  cooling      the cooling gate: atomic-cooling virial-temperature gate, a molecular-cooling mass gate, and a gate that
               rises to the atomic threshold after reionisation
The ranges are the plausible range of each quantity, not a posterior: where the literature gives a spread it is the spread, otherwise
a factor of two about the fiducial value.

Usage: python gen_mvm_v2_ledger.py --mass-index I --trees-from <dir>/trees_MII.npz --outdir <dir> [--n-trees 60] [--only name,...]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
MASSES = [3.0e10 * 10 ** (0.25 * k) for k in range(13)]

# (name, group, low-kwargs, high-kwargs); a missing side is None.  Fiducial values are those of run_reservoir_v2.
RANGES = [
    ("epsilon_sf", "physical", dict(epsilon_sf=0.0075), dict(epsilon_sf=0.03)),
    ("eta_acc", "physical", dict(eta_acc=0.0025), dict(eta_acc=0.01)),
    ("R_nuc_pc", "physical", dict(R_nuc_pc=50.0), dict(R_nuc_pc=200.0)),
    ("sigma_lnj", "physical", dict(sigma_lnj=0.4), dict(sigma_lnj=0.7)),
    ("lam_median", "physical", dict(lam_median=0.025), dict(lam_median=0.045)),
    ("c_nfw", "physical", dict(c_nfw=3.0), dict(c_nfw=6.0)),
    ("M_hot", "physical", dict(M_hot=2.0e11), dict(M_hot=8.0e11)),
    ("phi_coldhot", "physical", dict(phi_coldhot=2.0), dict(phi_coldhot=8.0)),
    ("f_mom", "physical", dict(f_mom=0.5), dict(f_mom=2.0)),
    ("eta_sn_scale", "physical", dict(eta_sn_scale=0.5), dict(eta_sn_scale=2.0)),
    ("epsilon_rad", "physical", dict(epsilon=0.05), dict(epsilon=0.2)),
    ("f_b", "physical", dict(f_b=0.14), dict(f_b=0.17)),
    ("dm_nuclear_off", "structural", dict(dm_nuclear=False), None),
    ("cosmology_frozen", "structural", dict(cosmology=None), None),
    ("f_BH", "structural", dict(_f_bh=0.3), dict(_f_bh=0.7)),
    ("n_sub", "numerical", dict(n_sub=1), dict(n_sub=32)),
    ("cool_atomic_1e4K", "cooling", dict(cooling_T_vir_K=1.0e4), None),
    ("cool_H2_mass_7e5", "cooling", dict(cooling_mass=7.0e5), None),
    ("cool_T_1e3K", "cooling", dict(cooling_T_vir_K=1.0e3), None),
    ("cool_atomic_after_reion_z6", "cooling", dict(cooling_T_vir_K=1.0e3, reionised_below_z=6.0), None),
    ("cool_atomic_after_reion_z10", "cooling", dict(cooling_T_vir_K=1.0e3, reionised_below_z=10.0), None),
]


def metrics(v2, Mh, tz, kw, n_iter):
    kw = dict(kw)
    f_bh = kw.pop("_f_bh", 0.5)
    Mc, never, above = v2.critical_seed_v2(Mh, tz, f_bh=f_bh, n_iter=n_iter, **kw)
    r = v2.run_reservoir_v2(Mh, tz, np.full(Mh.shape[0], 1e3), **kw)
    return Mc, r["bh_mass"][:, -1] / 1e3, int(never.sum() + above.sum())


def med(a):
    return float(np.nanmedian(a)) if np.isfinite(a).any() else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass-index", type=int, required=True)
    ap.add_argument("--trees-from", required=True)
    ap.add_argument("--n-trees", type=int, default=40)
    ap.add_argument("--n-steps", type=int, default=801)
    ap.add_argument("--n-sub", type=int, default=4)
    ap.add_argument("--n-iter", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=400)
    ap.add_argument("--only", default="")
    ap.add_argument("--outdir", default=str(HERE / "output" / "mvm_v2"))
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    from ashvini import reservoir_stock_v2 as v2

    d = np.load(a.trees_from)
    cosmo = v2.PaperCosmology()
    tz = v2.uniform_time_grid(25.0, 5.0, a.n_steps, cosmo)
    Mh = v2.interpolate_tree_onto_grid_v2(d["hm"], d["zz"], tz, cosmo)[:a.n_trees]
    base = dict(dm_nuclear=True, inflow="conservative", cosmology=cosmo, n_sub=a.n_sub)
    only = {s.strip() for s in a.only.split(",") if s.strip()}
    t0 = time.time()

    Mc0, G0, fail0 = metrics(v2, Mh, tz, base, a.n_iter)
    rng = np.random.default_rng(12345)
    boot = []
    for _ in range(a.n_boot):                       # statistical error of the baseline: resample trees with replacement
        i = rng.integers(0, len(Mc0), len(Mc0))
        boot.append((med(Mc0[i]), med(G0[i])))
    boot = np.array(boot)
    res = dict(M0=MASSES[a.mass_index], n_trees=int(Mh.shape[0]), n_sub=a.n_sub, n_steps=a.n_steps,
               baseline=dict(Mcrit=med(Mc0), G=med(G0), bracket_fail=fail0,
                             stat_dex_Mcrit=float(np.std(np.log10(boot[:, 0]))), stat_dex_G=float(np.std(np.log10(boot[:, 1])))),
               rows=[])
    print(f"baseline M0={MASSES[a.mass_index]:.2e} Mcrit={res['baseline']['Mcrit']:.3e} G={res['baseline']['G']:.3e} "
          f"stat {res['baseline']['stat_dex_Mcrit']:.3f}/{res['baseline']['stat_dex_G']:.3f} dex [{time.time()-t0:.0f}s]", flush=True)

    for name, group, lo, hi in RANGES:
        if only and name not in only:
            continue
        row = dict(name=name, group=group)
        for side, mod in (("low", lo), ("high", hi)):
            if mod is None:
                continue
            kw = dict(base)
            kw.update(mod)
            Mc, G, fail = metrics(v2, Mh, tz, kw, a.n_iter)
            with np.errstate(all="ignore"):
                row[side] = dict(setting={k: (v if not isinstance(v, v2.PaperCosmology) else "paper") for k, v in mod.items()},
                                 Mcrit=med(Mc), G=med(G), bracket_fail=fail,
                                 dlog_Mcrit=float(np.log10(med(Mc) / res["baseline"]["Mcrit"])),
                                 dlog_G=float(np.log10(med(G) / res["baseline"]["G"])))
        res["rows"].append(row)
        print(name, {s: (round(row[s]["dlog_Mcrit"], 3), round(row[s]["dlog_G"], 3), row[s]["bracket_fail"])
                     for s in ("low", "high") if s in row}, f"[{time.time()-t0:.0f}s]", flush=True)

    res["seconds"] = time.time() - t0
    tag = f"_{'-'.join(sorted(only))}" if only else ""
    with open(out / f"mvm_v2_ledger_M{a.mass_index:02d}{tag}.json", "w") as f:
        json.dump(res, f)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
