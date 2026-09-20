"""
Production ensemble for the Minimal Viable Model (ashvini.reservoir_stock), FROZEN fiducial:
R_nuc = 100 pc, sigma_j = 0.5, eta_acc = 0.005, eps_sf = 0.015, momentum AGN wind f_mom = 1, f_BH = 0.5.
801 uniform-time steps z=25 -> 5, 240 merger trees per halo mass (dz=0.05, M_res=1e4 Msun).

For each of 13 halo masses (3e10 -> 3e13 Msun, 0.25 dex) it records
  * the critical seed per tree and, at that seed, M*_tot, M_BH, the seed/host-baryon ratio;
  * the black hole fate M_BH(z=5)/M*_tot(z=5) for fixed seeds 1e2, 1e3, 2e5, 1e7 Msun in the fiducial model
    and in a labelled accessibility experiment (R_nuc = 250 pc, sigma_j = 0.5; NOT a candidate fiducial);
  * M_BH / (f_b M_halo) against redshift (median tree, and 16/50/84 percentiles over trees).
Writes output/mvm_production_results.json. Numba tree sampling is not seed-reproducible, so the stored JSON
is the record of the ensemble.
"""
import sys
import json
import time
import hashlib
import datetime
from pathlib import Path
from multiprocessing import Pool

FORAOIS_SRC = "/Users/00075868/MyCodes/foraois/src"
FORAOIS_CONFIG = "/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"
sys.path.insert(0, FORAOIS_SRC)

import numpy as np

HERE = Path(__file__).parent
OUTDIR = HERE / "output"
OUTDIR.mkdir(exist_ok=True)

N_TREES, N_STEPS, DZ, M_RES = 240, 801, 0.05, 1e4
MASSES = [3.0e10 * 10 ** (0.25 * k) for k in range(13)]
FIXED_SEEDS = [1e2, 1e3, 2e5, 1e7]
CASES = {"fiducial": dict(), "accessible_R250": dict(R_nuc_pc=250.0)}


def worker(args):
    i, M0 = args
    from foraois import cosmo_utils, ZhangHuiMergerTree
    from foraois.utils import io
    from ashvini import pymctrees_adapter
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    from ashvini.utils import time_at_z, z_at_time

    rp = io.get_params(FORAOIS_CONFIG)
    h = rp["Cosmology"]["h"]
    cd = cosmo_utils.CosmoData(rp, redshift=[5.0])
    tg = ZhangHuiMergerTree(cd, rp, model="cdm")
    t0 = time.time()
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(
        tg, M0, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ, backend="numba", rng_seed=1000 + i)
    tz = z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], N_STEPS))
    Mh = interpolate_tree_onto_grid(hm, zz, tz)
    rate = grumpy_halo_growth_rate(hm, zz, tz)
    Mc, never, above, F = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, return_F=True)
    o = mvm.run_reservoir_stock(Mh, tz, Mc, halo_growth_rate=rate)
    cum_in = (o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(axis=1)
    med = np.nanmedian(Mc)
    rep = int(np.nanargmin(np.abs(Mc - med)))
    ratio = o["bh_to_host_baryons"]
    with np.errstate(all="ignore"):
        good = np.isfinite(ratio).sum(axis=0) >= 20
        pct = np.where(good[None, :], np.nanpercentile(ratio, [16, 50, 84], axis=0), np.nan)
    res = dict(
        M0=M0, n_trees=N_TREES, bracket_fail=int(never.sum() + above.sum()), max_abs_F=float(np.nanmax(np.abs(F))),
        Mcrit=Mc.tolist(), Mstar=o["stars_mass"][:, -1].tolist(), MBH=o["bh_mass"][:, -1].tolist(),
        seed_over_host_baryons_at_formation=o["seed_over_host_baryons_at_formation"].tolist(),
        z_first_resolved=o["z_first_resolved"].tolist(),
        accreted_over_fb_Mhalo=(cum_in / (0.156 * Mh[:, -1])).tolist(),
        z=tz.tolist(), rep_index=rep, rep_Mcrit=float(Mc[rep]),
        rep_halo=Mh[rep].tolist(), rep_bh=o["bh_mass"][rep].tolist(), rep_bh_to_host=np.where(np.isfinite(ratio[rep]), ratio[rep], np.nan).tolist(),
        bh_to_host_p16=pct[0].tolist(), bh_to_host_p50=pct[1].tolist(), bh_to_host_p84=pct[2].tolist(),
        fixed_seed={},
    )
    for name, kw in CASES.items():
        res["fixed_seed"][name] = {}
        for s in FIXED_SEEDS:
            r = mvm.run_reservoir_stock(Mh, tz, np.full(Mh.shape[0], s), halo_growth_rate=rate, **kw)
            res["fixed_seed"][name][str(s)] = dict(ratio=(r["bh_mass"][:, -1] / r["stars_mass"][:, -1]).tolist(),
                                                   G=(r["bh_mass"][:, -1] / s).tolist(), Mstar=r["stars_mass"][:, -1].tolist())
    print(f"M0={M0:.3e} done in {time.time()-t0:.0f}s: median M_seed,crit {med:.3e}", flush=True)
    return res


if __name__ == "__main__":
    t0 = time.time()
    with Pool(8) as pool:
        results = pool.map(worker, list(enumerate(MASSES)), chunksize=1)
    from ashvini.paper_reservoir_params import PAPER_PARAMS
    import ashvini.reservoir_stock as _m
    sha = hashlib.sha256(Path(_m.__file__).read_bytes()).hexdigest()
    meta = dict(
        date=datetime.datetime.now().isoformat(timespec="seconds"), n_trees=N_TREES, n_steps=N_STEPS, dz=DZ, M_res=M_RES,
        z_seed=25.0, z_anchor=5.0, f_bh=0.5, module="ashvini/reservoir_stock.py", module_sha256=sha,
        fiducial=dict(R_nuc_pc=PAPER_PARAMS.R_nuc_pc, sigma_lnj=PAPER_PARAMS.sigma_lnj, eta_acc=PAPER_PARAMS.eta_acc,
                      epsilon_sf=PAPER_PARAMS.epsilon_sf, f_mom=PAPER_PARAMS.f_mom, eta_sn_scale=PAPER_PARAMS.eta_sn_scale,
                      n_rd=PAPER_PARAMS.n_rd, c_nfw=PAPER_PARAMS.c_nfw, epsilon=PAPER_PARAMS.epsilon, f_b=PAPER_PARAMS.f_b,
                      M_hot=PAPER_PARAMS.M_hot, lam_median=PAPER_PARAMS.lam_median),
        accessibility_experiment=CASES["accessible_R250"], fixed_seeds=FIXED_SEEDS)
    with open(OUTDIR / "mvm_production_results.json", "w") as f:
        json.dump(dict(meta=meta, results=results), f)
    print(f"DONE in {time.time()-t0:.0f}s -> {OUTDIR/'mvm_production_results.json'}")
