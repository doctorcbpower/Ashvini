"""
C18 test 4 (stellar wind on/off x AGN f_mom 0/1) and the phase-space diagnostic d ln R / dt against R, R = M_BH/M*_tot.
Frozen MVM (hash asserted, module unmodified); high delivery sigma_j = 1.5, R_nuc = 300 pc; 3e10 and 3e11 Msun only;
100 paired trees per mass, 801-step dt.

  Test 4 : seeds 1e3 and 1e7; wind on/off (eta_sn_scale = 1 / 0) x f_mom = 1 / 0.
  Phase  : d ln R/dt = d ln M_BH/dt - d ln M*_tot/dt over a +-20-step window (~26 Myr), pooled over seeds 1e2, 1e3, 1e5, 1e7,
           binned in log10 R (0.25 dex) in three redshift bands; plus a fiducial-delivery control (seeds 1e5, 1e7).
           R and its derivative are only defined where M*_tot > 1e5 Msun.
Writes output/c18_test4_phase_results.json.
"""
import sys, json, time
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, generators, trees, grid
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "output" / "c18_test4_phase_results.json"
MASSES = [3e10, 3e11]
HIGH = dict(sigma_lnj=1.5, R_nuc_pc=300.0)
FID = dict(sigma_lnj=0.5, R_nuc_pc=100.0)
STRIDE, W = 4, 20
BINS = np.arange(-3.0, 1.51, 0.25)
ZBANDS = (("z>10", 10.0, 99.0), ("7<z<=10", 7.0, 10.0), ("z<=7", 0.0, 7.0))


def phase(o, tz, seeds_tag):
    bh, Ms, t = o["bh_mass"], o["stars_mass"], o["cosmic_time"]
    ok = Ms > 1e5
    lnR = np.where(ok, np.log(np.where(ok, bh, 1.0) / np.where(ok, Ms, 1.0)), np.nan)
    dl = np.full_like(lnR, np.nan)
    dl[:, W:-W] = (lnR[:, 2 * W:] - lnR[:, :-2 * W]) / (t[2 * W:] - t[:-2 * W])[None, :]     # per Gyr
    Rv = np.exp(lnR)
    zz = np.broadcast_to(tz[None, :], lnR.shape)
    m = np.isfinite(dl) & np.isfinite(Rv)
    return Rv[m], dl[m], zz[m]


def bin_stats(R, dl, z):
    out = {}
    for name, lo, hi in ZBANDS:
        s = (z > lo) & (z <= hi)
        lr = np.log10(R[s]); d = dl[s]
        rows = []
        for a, b in zip(BINS[:-1], BINS[1:]):
            k = (lr >= a) & (lr < b)
            n = int(k.sum())
            rows.append([0.5 * (a + b), n] + ([float(np.median(d[k])), float(np.percentile(d[k], 16)), float(np.percentile(d[k], 84)), float(np.mean(d[k] < 0))] if n >= 30 else [np.nan] * 4))
        out[name] = rows
    return out


def crossings(rows, min_n=30):
    """log10 R where the median d ln R/dt goes from + (below) to - (above); interpolated"""
    x = [r for r in rows if r[1] >= min_n]
    res = []
    for a, b in zip(x[:-1], x[1:]):
        if a[2] > 0 >= b[2]:
            res.append(a[0] + (b[0] - a[0]) * a[2] / (a[2] - b[2]))
    return res


def worker(a):
    i, M0 = a
    assert_frozen()
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, 100, 25.0, 8000 + i)
    tz = grid(25.0)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    rng = np.random.default_rng(i)
    res = dict(M0=M0, z=tz[::STRIDE].tolist(), t4={}, phase={})

    def run(seed, **kw):
        return mvm.run_reservoir_stock(Mh, tz, np.full(Mh.shape[0], float(seed)), halo_growth_rate=rate, **kw)

    def traj(o):
        Ms, bh = o["stars_mass"], o["bh_mass"]
        R = np.where(Ms > 1e5, bh / np.where(Ms > 0, Ms, 1.0), np.nan)[:, ::STRIDE]
        with np.errstate(all="ignore"):
            good = np.isfinite(R).sum(axis=0) >= 20
            p = np.where(good[None, :], np.nanpercentile(R, [16, 50, 84], axis=0), np.nan)
        return p.tolist(), R[:8].tolist()

    for wind in (1.0, 0.0):
        for fm in (1.0, 0.0):
            for s in (1e3, 1e7):
                o = run(s, eta_sn_scale=wind, f_mom=fm, **HIGH)
                p, ind = traj(o)
                res["t4"][f"wind{wind:g}_f{fm:g}_{s:g}"] = dict(
                    final=(o["bh_mass"][:, -1] / o["stars_mass"][:, -1]).tolist(), Mbh=o["bh_mass"][:, -1].tolist(),
                    Mstar=o["stars_mass"][:, -1].tolist(), pct=p, ind=ind,
                    cum_out_agn=o["wind_taken_agn"].sum(axis=1).tolist(), cum_out_sn=o["wind_taken_sn"].sum(axis=1).tolist())
    cases = {"high_f1": (dict(f_mom=1.0, **HIGH), (1e2, 1e3, 1e5, 1e7)), "high_f0": (dict(f_mom=0.0, **HIGH), (1e2, 1e3, 1e5, 1e7)),
             "high_nowind": (dict(eta_sn_scale=0.0, **HIGH), (1e2, 1e3, 1e5, 1e7)), "fid_control": (dict(**FID), (1e5, 1e7))}
    for name, (kw, seeds) in cases.items():
        Rs, ds, zs = [], [], []
        for s in seeds:
            R, d, z = phase(run(s, **kw), tz, name)
            Rs.append(R); ds.append(d); zs.append(z)
        R, d, z = np.concatenate(Rs), np.concatenate(ds), np.concatenate(zs)
        st = bin_stats(R, d, z)
        pick = rng.choice(len(R), size=min(4000, len(R)), replace=False)
        res["phase"][name] = dict(stats=st, cross={k: crossings(v) for k, v in st.items()}, n=int(len(R)),
                                  pts=[np.log10(R[pick]).tolist(), d[pick].tolist(), z[pick].tolist()])
    return res


def q(x): return np.nanpercentile(np.array(x, dtype=float), [16, 50, 84])


if __name__ == "__main__":
    assert_frozen(); t0 = time.time()
    with Pool(2) as p:
        R = p.map(worker, list(enumerate(MASSES)), chunksize=1)
    json.dump(R, open(OUT, "w"))
    print(f"saved {OUT} in {time.time()-t0:.0f}s")
    for r in R:
        print(f"\n################ M0 = {r['M0']:.0e} ################\nTEST 4: final M_BH/M*_tot(z=5), median [16,84]; also median M_BH and M*_tot [Msun]; cumulative ejecta / M*_tot (AGN, stellar)")
        for s in (1e3, 1e7):
            print(f"  seed {s:.0e}")
            for wind in (1.0, 0.0):
                for fm in (1.0, 0.0):
                    a = r["t4"][f"wind{wind:g}_f{fm:g}_{s:g}"]; f = q(a["final"])
                    print(f"    stellar wind {'on ' if wind else 'off'}, f_mom={fm:g}: R = {f[1]:.3e} [{f[0]:.2e}, {f[2]:.2e}];  M_BH {np.median(a['Mbh']):.2e}, M* {np.median(a['Mstar']):.2e};  ejecta/M*: AGN {np.median(np.array(a['cum_out_agn'])/np.array(a['Mstar'])):.2f}, stellar {np.median(np.array(a['cum_out_sn'])/np.array(a['Mstar'])):.2f}")
        print("PHASE SPACE: zero crossings (+ below, - above) of the median d ln R/dt in log10 R, by redshift band; samples")
        for name, a in r["phase"].items():
            print(f"  {name:12s} n={a['n']:6d}: " + " | ".join(f"{b}: " + (", ".join(f"R={10**c:.3g}" for c in a['cross'][b]) if a['cross'][b] else "none") for b in a['cross']))
        print("  binned median [16,84] of d ln R / dt per Gyr for high_f1 (rows: log10 R centre, n)")
        for b, rows in r["phase"]["high_f1"]["stats"].items():
            print(f"   {b}: " + "; ".join(f"{x[0]:+.2f}:{x[2]:+.2f}[{x[3]:+.1f},{x[4]:+.1f}](n={x[1]})" for x in rows if x[1] >= 30))
