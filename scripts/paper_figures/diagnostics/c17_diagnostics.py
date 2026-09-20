"""
C17 diagnostics branching from the frozen MVM (695b114). The frozen module is NOT modified: the script refuses to run
unless ashvini/reservoir_stock.py has the sha256 recorded in the production JSON.

  zseed : z_seed in {15, 20, 25, 30, 35} at fixed time step (the production dt); one tree set per mass (z_max = 40), paired.
  mhot  : M_hot x0.5, x1, x2 across 7 halo masses; paired trees.
  nrd   : n_rd (galaxy scale in disc scale lengths) in {5, 10, 20}; paired trees.
  trees : Zhang & Hui (fiducial) against PCH08, 240 independent trees each, at 3e10, 3e11, 3e13.

Usage: python c17_diagnostics.py <test>
"""
import sys, json, time, hashlib
from pathlib import Path
from multiprocessing import Pool

sys.path.insert(0, "/Users/00075868/MyCodes/foraois/src")
import numpy as np

CONFIG = "/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"
PROD = Path(__file__).resolve().parents[1] / "output" / "mvm_production_results.json"
N_STEPS0 = 801


def assert_frozen():
    import ashvini.reservoir_stock as m
    sha = hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()
    frozen = json.load(open(PROD))["meta"]["module_sha256"]
    if sha != frozen:
        raise SystemExit(f"reservoir_stock.py has changed since the frozen production run ({sha[:12]} != {frozen[:12]})")


def generators():
    from foraois import cosmo_utils, ZhangHuiMergerTree, PCHMergerTree
    from foraois.utils import io
    rp = io.get_params(CONFIG)
    cd = cosmo_utils.CosmoData(rp, redshift=[5.0])
    return rp["Cosmology"]["h"], ZhangHuiMergerTree(cd, rp, model="cdm"), PCHMergerTree(cd, rp)


def grid(z_seed, n_steps0=N_STEPS0):
    """uniform-time grid from z_seed to z=5 at the production step (dt of the 801-step z=25->5 grid)"""
    from ashvini.utils import time_at_z, z_at_time
    t5, t25 = float(time_at_z(np.array([5.0]))[0]), float(time_at_z(np.array([25.0]))[0])
    dt0 = (t5 - t25) / (n_steps0 - 1)
    ts = float(time_at_z(np.array([z_seed]))[0])
    n = int(round((t5 - ts) / dt0)) + 1
    return z_at_time(np.linspace(ts, t5, n))


def trees(gen, h, M0, N, z_max, seed, dz=0.05, mres=1e4):
    from ashvini import pymctrees_adapter
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(gen, M0, h, N, z0=5.0, z_max=z_max, m_res_msun=mres, dz=dz, backend="numba", rng_seed=seed)
    return hm, zz


def solve(hm, zz, tz, **kw):
    from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
    from ashvini import reservoir_stock as mvm
    Mh = interpolate_tree_onto_grid(hm, zz, tz)
    rate = grumpy_halo_growth_rate(hm, zz, tz)
    Mc, never, above = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, **kw)
    o = mvm.run_reservoir_stock(Mh, tz, Mc, halo_growth_rate=rate, **{k: v for k, v in kw.items() if k != "f_bh"})
    cum_in = (o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(axis=1)
    return dict(Mc=Mc, Mstar=o["stars_mass"][:, -1], G=o["bh_mass"][:, -1] / Mc, Mh=Mh, tz=tz,
                acc=cum_in / (0.156 * Mh[:, -1]), bad=int(never.sum() + above.sum()))


med = np.nanmedian
def q(x): return np.nanpercentile(x, [16, 84])


# ------------------------------------------------------------------ workers
def w_zseed(a):
    i, M0 = a
    assert_frozen(); h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, 100, 40.0, 2000 + i)
    out = {}
    for zs in (15, 20, 25, 30, 35):
        out[zs] = solve(hm, zz, grid(float(zs)))
    return M0, out


def w_mhot(a):
    i, M0 = a
    assert_frozen(); h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, 100, 25.0, 3000 + i)
    out = {}
    for f in (0.5, 1.0, 2.0):
        out[f] = solve(hm, zz, grid(25.0), M_hot=4.0e11 * f)
    return M0, out


def w_nrd(a):
    i, M0 = a
    assert_frozen(); h, zh, _ = generators()
    hm, zz = trees(zh, h, M0, 100, 25.0, 4000 + i)
    out = {}
    for n in (5.0, 10.0, 20.0):
        out[n] = solve(hm, zz, grid(25.0), n_rd=n)
    return M0, out


def w_trees(a):
    i, M0, which = a
    assert_frozen(); h, zh, pch = generators()
    hm, zz = trees(zh if which == "ZH" else pch, h, M0, 240, 25.0, 5000 + i + (100 if which == "PCH" else 0))
    return M0, which, solve(hm, zz, grid(25.0))


# ------------------------------------------------------------------ reports
def rep_zseed(res):
    print("z_seed test: fixed dt (production), one tree set per mass (100 trees, z_max=40); ratios paired per tree to z_seed=25")
    print("NOTE: the MVM starts with an empty galaxy at z_seed, so a later z_seed also discards baryons the halo had already accreted.")
    for M0, out in sorted(res):
        base = out[25]["Mc"]
        print(f"\nM0={M0:.0e}   (steps at z_seed=15/20/25/30/35: {[len(out[z]['tz']) for z in (15,20,25,30,35)]})")
        print("  z_seed | median M_seed,crit | paired ratio to 25 [min,max] | halo already resolved at z_seed | M*_tot/M_halo | G_BH | bracket fails")
        for zs in (15, 20, 25, 30, 35):
            r = out[zs]; ratio = r["Mc"] / base
            print(f"   {zs:3d}   | {med(r['Mc']):.3e}        | {med(ratio):.4f} [{np.nanmin(ratio):.3f}, {np.nanmax(ratio):.3f}]      | {np.mean(r['Mh'][:,0]>0)*100:5.1f}%                        | {med(r['Mstar']/r['Mh'][:,-1]):.2e}     | {med(r['G']):.3f} | {r['bad']}")


def rep_mhot(res):
    print("M_hot test: paired 100 trees per mass, 801 steps; M_seed,crit and M*_tot/M_halo for M_hot = 2e11, 4e11 (fiducial), 8e11")
    print("  M0        | M_seed,crit [2e11 | 4e11 | 8e11]            | paired ratio (2e11, 8e11) | M*_tot/M_halo [2e11 | 4e11 | 8e11] | accreted/(f_b Mh) [2e11 | 4e11 | 8e11]")
    for M0, out in sorted(res):
        b = out[1.0]["Mc"]
        print(f"  {M0:.2e} | {med(out[0.5]['Mc']):.2e} | {med(b):.2e} | {med(out[2.0]['Mc']):.2e} | {med(out[0.5]['Mc']/b):.3f}, {med(out[2.0]['Mc']/b):.3f} | "
              f"{med(out[0.5]['Mstar']/out[0.5]['Mh'][:,-1]):.2e} | {med(out[1.0]['Mstar']/out[1.0]['Mh'][:,-1]):.2e} | {med(out[2.0]['Mstar']/out[2.0]['Mh'][:,-1]):.2e} | "
              f"{med(out[0.5]['acc']):.3f} | {med(out[1.0]['acc']):.3f} | {med(out[2.0]['acc']):.3f}")
    Ms = sorted(m for m, _ in res)
    d = {f: [med(dict(res)[m][f]["Mc"]) for m in Ms] for f in (0.5, 1.0, 2.0)}
    for f in (0.5, 1.0, 2.0):
        sl = [np.log10(d[f][j + 1] / d[f][j]) / np.log10(Ms[j + 1] / Ms[j]) for j in range(len(Ms) - 1)]
        print(f"  M_hot x{f}: local slopes between successive masses: " + ", ".join(f"{s:.2f}" for s in sl))


def rep_nrd(res):
    print("n_rd test: paired 100 trees per mass, 801 steps; n_rd = 5, 10 (fiducial), 20")
    for M0, out in sorted(res):
        b = out[10.0]["Mc"]
        print(f"  M0={M0:.0e}: median M_seed,crit " + " | ".join(f"n_rd={n:g}: {med(out[n]['Mc']):.3e} (paired ratio {med(out[n]['Mc']/b):.3f} [{np.nanmin(out[n]['Mc']/b):.3f},{np.nanmax(out[n]['Mc']/b):.3f}], M*/Mh {med(out[n]['Mstar']/out[n]['Mh'][:,-1]):.2e}, G {med(out[n]['G']):.3f})" for n in (5.0, 10.0, 20.0)))


def rep_trees(res):
    print("tree-algorithm test: Zhang & Hui (fiducial) against PCH08; 240 independent trees each, dz=0.05, M_res=1e4, 801 steps")
    d = {(M0, w): r for M0, w, r in res}
    for M0 in sorted({m for m, _ in d}):
        print(f"\nM0={M0:.0e}")
        for w in ("ZH", "PCH"):
            r = d[(M0, w)]; a, b = q(r["Mc"])
            tz = r["tz"]; at = lambda z: r["Mh"][:, int(np.argmin(abs(tz - z)))]
            print(f"  {w:3s}: median M_seed,crit {med(r['Mc']):.3e} [16,84 = {a:.2e}, {b:.2e}], 16-84 width {np.log10(b/a):.3f} dex; M*_tot/M_halo {med(r['Mstar']/r['Mh'][:,-1]):.2e}; G_BH {med(r['G']):.3f}; "
                  f"median M_halo at z=20/15/10: {med(np.where(at(20)>0, at(20), np.nan)):.2e} / {med(at(15)):.2e} / {med(at(10)):.2e}; halo resolved at z=20: {np.mean(at(20)>0)*100:.0f}%; accreted/(f_b Mh) {med(r['acc']):.3f}")
        z, p = d[(M0, "ZH")], d[(M0, "PCH")]
        print(f"  PCH08/ZH ratio of medians: M_seed,crit {med(p['Mc'])/med(z['Mc']):.3f} ({np.log10(med(p['Mc'])/med(z['Mc'])):+.3f} dex); M*_tot {med(p['Mstar'])/med(z['Mstar']):.3f}")


if __name__ == "__main__":
    test = sys.argv[1]
    assert_frozen()
    t0 = time.time()
    masses = {"zseed": [3e10, 3e11, 3e13], "nrd": [3e10, 3e11, 3e13],
              "mhot": [3e10 * 10 ** (0.25 * k) for k in range(0, 13, 2)]}
    if test == "trees":
        tasks = [(i, m, w) for i, m in enumerate([3e10, 3e11, 3e13]) for w in ("ZH", "PCH")]
        with Pool(6) as p:
            res = p.map(w_trees, tasks, chunksize=1)
        rep_trees(res)
    else:
        fn, rep = dict(zseed=(w_zseed, rep_zseed), mhot=(w_mhot, rep_mhot), nrd=(w_nrd, rep_nrd))[test]
        with Pool(min(8, len(masses[test]))) as p:
            res = p.map(fn, list(enumerate(masses[test])), chunksize=1)
        rep(res)
    print(f"\nDONE {test} in {time.time()-t0:.0f}s")
