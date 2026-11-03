"""
How the mass and redshift dependence of the supernova mass loading shapes the stellar-to-halo mass relation (docs/WIND_SCALING_AUDIT.md).

Each form of the loading has one amplitude. It is set by a secant search (in log amplitude) so that the median M* at 3e10, 1e11, 3e11 and 1e12 Msun
is as close as possible, on average in log, to the geometric mean of Moster+13 and Behroozi+13 at z = 0. What is then reported is what the form does
not control: the scatter of the residuals about that mean over the four masses (can the shape be matched?), the relation at 1e8 to 1e10 Msun, and
the redshift evolution at fixed halo mass (M*/M_halo of the main progenitors pooled over haloes, at z = 0, 2, 5).

Fixed: equilibrium UV form, current star formation law with efficiency 0.06 (the middle of the accepted range), the stored trees.
Usage: python scripts/wind_scaling_scan.py [name ...]   (no names: all)   Output: $SHMR_OUT/wind_scan/<name>.json
"""
import json
import os
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C
from stellar_to_halo_mass_relation import behroozi2013_mstar, moster2013_mstar

warnings.filterwarnings("ignore")
BASE = {"reion.form": "equilibrium", "sf.efficiency": 0.06}
CAL = [5, 6, 7, 8]                     # indices of 3e10, 1e11, 3e11, 1e12 in C.MASSES
ZS = [0.0, 2.0, 5.0]


def _forms():
    out = {}
    for ms, zs in [(1 / 3, 0.5), (1 / 3, 0.25), (1 / 3, 0.0), (0.5, 0.5), (2 / 3, 0.5), (2 / 3, 0.0), (1.0, 0.5)]:
        out[f"ash_m{ms:.2f}_z{zs:.2f}"] = ({"sn.mass_slope": ms, "sn.redshift_slope": zs}, "sn.epsilon_p", 8.0)
    for ep in (0.2, 0.35, 0.5):
        out[f"stellar_p{ep:.2f}_cut"] = ({"sn.wind_model": "stellar_mass", "sn.eta_p": ep, "sn.eta_c": 4.5}, "sn.eta_norm", 10.0)
        out[f"stellar_p{ep:.2f}_nocut"] = ({"sn.wind_model": "stellar_mass", "sn.eta_p": ep, "sn.eta_c": 0.0}, "sn.eta_norm", 3.0)
    return out


def target(m):
    return 0.5 * (np.log10(float(behroozi2013_mstar(m))) + np.log10(float(moster2013_mstar(m))))


def run_masses(cfg, trees, idx, record=False):
    from ashvini.main import run_forest
    med, rec = [], []
    for i in idx:
        t = trees[i]
        r = run_forest(t["halo_masses"], t["halo_growth_rates"], t["redshifts"], config=cfg)
        s = r["stars_mass"]
        med.append(np.median(s[:, -1]))
        if record:
            for z in ZS:
                k = int(np.argmin(abs(t["redshifts"] - z)))
                rec.append((i, z, t["halo_masses"][:, k].copy(), s[:, k].copy()))
    return np.array(med), rec


def residuals(cfg, trees):
    med, _ = run_masses(cfg, trees, CAL)
    return np.log10(np.maximum(med, 1e-30)) - np.array([target(C.MASSES[i]) for i in CAL])


def fit_amplitude(changes, key, a0, trees):
    from ashvini.config import DEFAULT, override
    def f(la):
        return residuals(override(DEFAULT, **{**BASE, **changes, key: 10 ** la}), trees)
    la = [np.log10(a0), np.log10(a0) + 0.15]
    r = [f(la[0]).mean(), f(la[1]).mean()]
    for _ in range(6):
        if abs(r[-1]) < 0.03:
            break
        slope = (r[-1] - r[-2]) / (la[-1] - la[-2])
        nxt = la[-1] - r[-1] / (slope if abs(slope) > 1e-3 else -1.0)
        la.append(float(np.clip(nxt, la[-1] - 0.7, la[-1] + 0.7)))
        r.append(f(la[-1]).mean())
    best = int(np.argmin(np.abs(r)))
    return 10 ** la[best], r[best]


def pooled(rec, z, lo, hi):
    mh = np.concatenate([h for (_, zz, h, s) in rec if zz == z])
    ms = np.concatenate([s for (_, zz, h, s) in rec if zz == z])
    sel = (mh >= lo) & (mh < hi) & (mh > 0)
    return float(np.median(ms[sel] / mh[sel])) if sel.sum() >= 30 else None


def main(names):
    from ashvini.config import DEFAULT, override
    forms = _forms()
    trees = [dict(np.load(C.tree_path(i))) for i in range(len(C.MASSES))]
    os.makedirs(os.path.join(C.OUT, "wind_scan"), exist_ok=True)
    for name in names or forms:
        path = os.path.join(C.OUT, "wind_scan", f"{name}.json")
        if os.path.exists(path):
            print(name, "exists, skipped"); continue
        changes, key, a0 = forms[name]
        t0 = time.time()
        amp, mean_res = fit_amplitude(changes, key, a0, trees)
        cfg = override(DEFAULT, **{**BASE, **changes, key: amp})
        med, rec = run_masses(cfg, trees, range(len(C.MASSES)), record=True)
        res = np.log10(np.maximum(med[CAL], 1e-30)) - np.array([target(C.MASSES[i]) for i in CAL])
        edges = {"1e8-3e8": (1e8, 3e8), "1e9-3e9": (1e9, 3e9), "1e10-3e10": (1e10, 3e10), "1e11-3e11": (1e11, 3e11)}
        out = {"name": name, "amplitude_key": key, "amplitude": amp, "mean_residual_dex": float(mean_res),
               "residual_dex_at_3e10_1e11_3e11_1e12": res.tolist(), "residual_rms_dex": float(np.sqrt(np.mean((res - res.mean()) ** 2))),
               "ratio_to_mosterbehroozi_mean_dex_all_masses": (np.log10(np.maximum(med, 1e-30)) - np.array([target(m) for m in C.MASSES])).tolist(),
               "mstar_over_mhalo_z0": (med / C.MASSES).tolist(),
               "mstar_over_mhalo_by_z": {f"z{z:g}": {k: pooled(rec, z, *v) for k, v in edges.items()} for z in ZS},
               "seconds": time.time() - t0}
        json.dump(out, open(path, "w"), indent=1)
        print(f"{name}: amplitude {amp:.3g}, residuals {np.round(res, 2)}, rms about mean {out['residual_rms_dex']:.2f} dex, {out['seconds']:.0f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
