"""
C18 test 5/6: does the no-feedback floor of R = M_BH/M*_tot scale as eta_acc / eps_sf?
Frozen MVM (hash asserted; unmodified). High delivery (sigma_j = 1.5, R_nuc = 300 pc); 3e10 and 3e11 Msun; 100 paired trees;
801-step dt. eta_acc in {0.001, 0.005, 0.02} x eps_sf in {0.0075, 0.015, 0.03}; two feedback cells: all feedback off
(f_mom = 0, eta_sn_scale = 0) and AGN + stellar wind on (fiducial coupling). Seeds 1e3 and 1e7 for the final ratio; seeds
1e2, 1e3, 1e5, 1e7 pooled for the phase-space zero crossing (z <= 7 band).
Prediction under test: R_floor ~ c * eta_acc / eps_sf (c of order 0.3 to 1), a function of the ratio only.
"""
import sys, json, time
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, generators, trees, grid
from c18_test4_phase import phase, bin_stats, crossings
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "output" / "c18_eta_eps_results.json"
MASSES = [3e10, 3e11]
ETAS, EPSS = (0.001, 0.005, 0.02), (0.0075, 0.015, 0.03)
HIGH = dict(sigma_lnj=1.5, R_nuc_pc=300.0)
CELLS = {"nofb": dict(f_mom=0.0, eta_sn_scale=0.0), "fb": dict()}


def worker(a):
    i, M0 = a
    assert_frozen()
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, 100, 25.0, 9000 + i)
    tz = grid(25.0)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    res = dict(M0=M0, rows=[])

    def run(seed, **kw):
        return mvm.run_reservoir_stock(Mh, tz, np.full(Mh.shape[0], float(seed)), halo_growth_rate=rate, **kw)

    for cell, ckw in CELLS.items():
        for eta in ETAS:
            for eps in EPSS:
                kw = dict(eta_acc=eta, epsilon_sf=eps, **HIGH, **ckw)
                row = dict(cell=cell, eta=eta, eps=eps, ratio=eta / eps, seeds={})
                for s in (1e3, 1e7):
                    o = run(s, **kw)
                    bh, Ms = o["bh_mass"], o["stars_mass"]
                    kap = o["kappa_edd"]; A = o["A_bh"][:, 1:]; M = bh[:, :-1]
                    capped = (A > 0) & (A > kap * M)
                    grow = np.diff(bh, axis=1)
                    row["seeds"][f"{s:g}"] = dict(
                        R=(bh[:, -1] / Ms[:, -1]).tolist(),
                        dBH_dNuc=((bh[:, -1] - s) / np.maximum(o["stars_nuc"][:, -1], 1e-300)).tolist(),
                        dBH_dTot=((bh[:, -1] - s) / Ms[:, -1]).tolist(),
                        nuc_share=(o["stars_nuc"][:, -1] / Ms[:, -1]).tolist(),
                        cap_steps=float(np.mean(capped.sum(1) / np.maximum((A > 0).sum(1), 1))),
                        cap_mass=float(np.median((grow * capped).sum(1) / np.maximum(grow.sum(1), 1e-300))))
                Rs, ds, zs = [], [], []
                for s in (1e2, 1e3, 1e5, 1e7):
                    R, d, z = phase(run(s, **kw), tz, cell)
                    Rs.append(R); ds.append(d); zs.append(z)
                st = bin_stats(np.concatenate(Rs), np.concatenate(ds), np.concatenate(zs))
                row["cross_late"] = crossings(st["z<=7"]); row["cross_mid"] = crossings(st["7<z<=10"])
                res["rows"].append(row)
    return res


def q(x): return np.nanpercentile(np.array(x, dtype=float), [16, 50, 84])


if __name__ == "__main__":
    assert_frozen(); t0 = time.time()
    with Pool(2) as p:
        R = p.map(worker, list(enumerate(MASSES)), chunksize=1)
    json.dump(R, open(OUT, "w"))
    print(f"saved {OUT} in {time.time()-t0:.0f}s")
    for r in R:
        for cell in ("nofb", "fb"):
            print(f"\n######## M0={r['M0']:.0e}, cell = {'ALL FEEDBACK OFF' if cell == 'nofb' else 'AGN + stellar wind on'} ########")
            print("  eta_acc  eps_sf  eta/eps | seed 1e3: R median [16,84] (R/(eta/eps)), dBH/dM*_nuc, nuc share, Edd-capped steps/mass | seed 1e7: R (R/(eta/eps)) | phase crossing z<=7 (z 7-10)")
            rows = [x for x in r["rows"] if x["cell"] == cell]
            for x in sorted(rows, key=lambda y: y["ratio"]):
                a, b = x["seeds"]["1000"], x["seeds"]["1e+07"] if "1e+07" in x["seeds"] else x["seeds"]["10000000"]
                f = q(a["R"]); g = q(b["R"])
                cr = lambda c: ", ".join(f"{10**v:.3g}" for v in c) if c else "none"
                print(f"  {x['eta']:.3f}   {x['eps']:.4f}  {x['ratio']:.3f}   | {f[1]:.3e} [{f[0]:.1e},{f[2]:.1e}] ({f[1]/x['ratio']:.2f}), {np.median(a['dBH_dNuc']):.3f}, {np.median(a['nuc_share']):.2f}, {a['cap_steps']*100:.1f}%/{a['cap_mass']*100:.1f}%"
                      f" | {g[1]:.3e} ({g[1]/x['ratio']:.2f}) | {cr(x['cross_late'])} ({cr(x['cross_mid'])})")
            # scaling fit: log R vs log(eta/eps), each seed
            for sk in ("1000", "1e+07" if "1e+07" in rows[0]["seeds"] else "10000000"):
                xs = np.log([x["ratio"] for x in rows]); ys = np.log([np.median(x["seeds"][sk]["R"]) for x in rows])
                sl, ic = np.polyfit(xs, ys, 1)
                # degeneracy: same ratio, different (eta, eps): 0.005/0.0075 vs 0.02/0.03
                d1 = [np.median(x["seeds"][sk]["R"]) for x in rows if abs(x["ratio"] - 0.6667) < 1e-3]
                print(f"  seed {sk}: log-log slope of median R against eta/eps = {sl:.2f}, intercept factor c = {np.exp(ic):.2f}; two grid points with ratio 0.667 give R = " + " and ".join(f"{v:.3e}" for v in d1))
