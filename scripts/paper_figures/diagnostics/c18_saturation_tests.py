"""
C18: is the 0.1 to 0.2 M_BH/M*_tot scale a real, seed-independent outcome, and what does it depend on?
Tests 1 to 3 of the agreed ranking, on paired identical trees (100 per halo mass, 801-step dt), branching from the frozen
MVM (hash asserted; reservoir_stock.py is NOT modified).

  1. seed memory : seeds 1e2, 1e3, 1e5, 1e7 Msun at a high-accessibility configuration (sigma_j = 1.5, R_nuc = 300 pc),
                   plus the fiducial (sigma_j = 0.5, R_nuc = 100 pc) as reference.
  2. AGN strength: f_mom = 0, 0.3, 1, 3, 10 at the high-accessibility configuration (seeds 1e3 and 1e7).
  3. all gas accessible: every galaxy -> nucleus transfer succeeds (the nuclear_transfer function is replaced at run time by one
                   that moves all galaxy gas), at R_nuc = 100 and 300 pc, against the fiducial and high-accessibility models.

Writes output/c18_saturation_results.json (per-tree final ratios and downsampled trajectories, plus 8 individual trees).
"""
import sys, json, time
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, generators, trees, grid
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "output" / "c18_saturation_results.json"
MASSES = [3e10, 3e11, 3e13]
NT, STRIDE, N_IND = 100, 4, 8
HIGH = dict(sigma_lnj=1.5, R_nuc_pc=300.0)
FID = dict(sigma_lnj=0.5, R_nuc_pc=100.0)


def all_transfer(m, ln_jmed, F_up, F_cur, cut, ln_jt_nuc, sigma_lnj, n_active):
    """control: every galaxy gas parcel is accessible to the nucleus"""
    T = m[:, :n_active].sum(axis=1).copy()
    m[:, :n_active] = 0.0
    return T


def run(mvm, Mh, tz, rate, seed, all_access=False, **kw):
    orig = mvm.nuclear_transfer
    if all_access:
        mvm.nuclear_transfer = all_transfer
    try:
        o = mvm.run_reservoir_stock(Mh, tz, np.full(Mh.shape[0], float(seed)), halo_growth_rate=rate, **kw)
    finally:
        mvm.nuclear_transfer = orig
    N = Mh.shape[0]
    Ms, bh = o["stars_mass"], o["bh_mass"]
    R = np.where(Ms > 1e5, bh / np.where(Ms > 0, Ms, 1.0), np.nan)          # ratio only where stars exist (M* > 1e5 Msun)
    accepted = (o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(axis=1) - o["gas_unavailable"][:, -1]
    with np.errstate(all="ignore"):
        p = np.where((np.isfinite(R).sum(axis=0) >= 20)[None, :], np.nanpercentile(R[:, ::STRIDE], [16, 50, 84], axis=0).repeat(1, 0) if False else np.nan, np.nan)
    Rd = R[:, ::STRIDE]
    with np.errstate(all="ignore"):
        good = np.isfinite(Rd).sum(axis=0) >= 20
        pct = np.where(good[None, :], np.nanpercentile(Rd, [16, 50, 84], axis=0), np.nan)
    return dict(final=(bh[:, -1] / Ms[:, -1]).tolist(), G=(bh[:, -1] / float(seed)).tolist(), pct=pct.tolist(),
                ind=Rd[:N_IND].tolist(), deliv=(o["transfer_nuc"].sum(axis=1) / np.maximum(accepted, 1e-300)).tolist(),
                Mstar=Ms[:, -1].tolist())


def worker(a):
    i, M0 = a
    assert_frozen()
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, NT, 25.0, 7000 + i)
    tz = grid(25.0)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    res = dict(M0=M0, z=tz[::STRIDE].tolist(), t1={}, t2={}, t3={})
    for s in (1e2, 1e3, 1e5, 1e7):
        res["t1"][f"high_{s:g}"] = run(mvm, Mh, tz, rate, s, **HIGH)
        res["t1"][f"fid_{s:g}"] = run(mvm, Mh, tz, rate, s, **FID)
    for f in (0.0, 0.3, 1.0, 3.0, 10.0):
        for s in (1e3, 1e7):
            res["t2"][f"f{f:g}_{s:g}"] = run(mvm, Mh, tz, rate, s, f_mom=f, **HIGH)
    for R in (100.0, 300.0):
        for s in (1e2, 1e3, 1e5, 1e7):
            res["t3"][f"all_R{R:g}_{s:g}"] = run(mvm, Mh, tz, rate, s, all_access=True, sigma_lnj=0.5, R_nuc_pc=R)
    return res


def q(x): return np.nanpercentile(np.array(x, dtype=float), [16, 50, 84])


if __name__ == "__main__":
    assert_frozen(); t0 = time.time()
    with Pool(3) as p:
        R = p.map(worker, list(enumerate(MASSES)), chunksize=1)
    json.dump(R, open(OUT, "w"))
    print(f"saved {OUT} in {time.time()-t0:.0f}s (100 trees per mass, 801-step dt; R = M_BH/M*_tot, masked where M*_tot < 1e5 Msun)")
    for r in R:
        M0 = r["M0"]; z = np.array(r["z"])
        print(f"\n################ M0 = {M0:.0e} ################")
        print("TEST 1  seed memory: final ratio M_BH/M*(z=5), median [16,84]; high accessibility (sigma_j=1.5, R_nuc=300) and fiducial")
        for cfg in ("high", "fid"):
            for s in (1e2, 1e3, 1e5, 1e7):
                a = r["t1"][f"{cfg}_{s:g}"]; f = q(a["final"])
                print(f"   {cfg:4s} seed {s:.0e}: {f[1]:.3e} [{f[0]:.2e}, {f[2]:.2e}]   G median {np.nanmedian(a['G']):.4g}")
        # per-tree seed memory: spread of ln R across seeds, at several redshifts
        for cfg in ("high", "fid"):
            F = np.array([r["t1"][f"{cfg}_{s:g}"]["final"] for s in (1e2, 1e3, 1e5, 1e7)])       # (seeds, trees)
            spread = np.log(F.max(0) / F.min(0))
            light = np.log(np.max(F[:3], 0) / np.min(F[:3], 0))
            print(f"   {cfg:4s} per-tree max/min of final ratio across the 4 seeds: median {np.median(np.exp(spread)):.3g} [16,84 = {np.percentile(np.exp(spread),16):.3g}, {np.percentile(np.exp(spread),84):.3g}];  across the three lightest (1e2,1e3,1e5): median {np.median(np.exp(light)):.3g}")
            for zt in (15.0, 10.0, 7.0, 5.0):
                j = int(np.argmin(abs(z - zt)))
                med = np.array([r["t1"][f"{cfg}_{s:g}"]["pct"][1][j] for s in (1e2, 1e3, 1e5, 1e7)], dtype=float)
                print(f"        z={zt:4.1f}: median ratio by seed (1e2, 1e3, 1e5, 1e7) = " + ", ".join(f"{m:.2e}" for m in med) + f"   max/min = {np.nanmax(med)/np.nanmin(med):.3g}")
        print("TEST 2  AGN strength (high accessibility): final ratio median [16,84] for f_mom = 0, 0.3, 1, 3, 10")
        for s in (1e3, 1e7):
            row = []
            for f in (0.0, 0.3, 1.0, 3.0, 10.0):
                a = r["t2"][f"f{f:g}_{s:g}"]; qq = q(a["final"]); row.append(f"f={f:g}: {qq[1]:.2e} [{qq[0]:.1e},{qq[2]:.1e}]")
            print(f"   seed {s:.0e}: " + " | ".join(row))
        print("TEST 3  all gas accessible: final ratio median [16,84] and fraction of galaxy gas delivered to the nucleus")
        for s in (1e2, 1e3, 1e5, 1e7):
            row = []
            for tag, key in (("fiducial", f"fid_{s:g}"), ("high", f"high_{s:g}"), ("all R100", f"all_R100_{s:g}"), ("all R300", f"all_R300_{s:g}")):
                a = (r["t1"] if key.startswith(("fid", "high")) else r["t3"])[key]; qq = q(a["final"])
                row.append(f"{tag}: {qq[1]:.2e} [{qq[0]:.1e},{qq[2]:.1e}] (deliv {np.nanmedian(a['deliv']):.2f})")
            print(f"   seed {s:.0e}: " + " | ".join(row))
