import sys, time
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm
from ashvini.utils import time_at_z, z_at_time
M0 = float(sys.argv[1]); i0 = int(sys.argv[2])
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
cd = cosmo_utils.CosmoData(rp, redshift=[5.0]); tg = ZhangHuiMergerTree(cd, rp, model="cdm")
tz = z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], 801))
hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 20, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=i0)
Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
N = Mh.shape[0]
def diag(seed, **kw):
    o = mvm.run_reservoir_stock(Mh, tz, seed, halo_growth_rate=rate, **kw)
    kap, bh = o["kappa_edd"], o["bh_mass"]; A = o["A_bh"][:, 1:]; M = bh[:, :-1]
    active = A > 0
    ratio = np.where(active, A / (kap * M), 0.0)              # Mdot_acc / Mdot_Edd
    capped = active & (ratio > 1.0)
    grow = np.diff(bh, axis=1)
    fmass = (grow * capped).sum(1) / np.maximum(grow.sum(1), 1e-300)
    G = bh[:, -1] / bh[:, 0]
    return dict(rmax=np.median(ratio.max(1)), rmax_hi=ratio.max(), fcap=np.median(capped.sum(1) / np.maximum(active.sum(1), 1)),
                fmass=np.median(fmass), G=np.median(G), stars=np.median(o["stars_mass"][:, -1]),
                deliv=np.median(o["transfer_nuc"].sum(1) / np.maximum(((o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(1) - o["gas_unavailable"][:, -1]), 1e-300)))
def row(tag, **kw):
    Mc, *_ = mvm.critical_seed_stock(Mh, tz, n_iter=36, halo_growth_rate=rate, **kw)
    c = diag(Mc, **kw); a = diag(np.full(N, 1e3), **kw); b = diag(np.full(N, 2e5), **kw)
    print(f"{tag:26s}| deliv {c['deliv']:.1e} | crit: M_crit {np.nanmedian(Mc):9.2e} G {c['G']:8.3f} Mdot/Edd max {c['rmax']:8.2e} capped {c['fcap']*100:5.1f}% | seed 1e3: G {a['G']:9.3g} max {a['rmax']:8.2e} (worst tree {a['rmax_hi']:.1e}) capped {a['fcap']*100:5.1f}% ({a['fmass']*100:4.1f}% of mass) | seed 2e5: G {b['G']:8.3f} max {b['rmax']:8.2e} capped {b['fcap']*100:5.1f}%", flush=True)
def peak(seed, **kw):
    o = mvm.run_reservoir_stock(Mh, tz, seed, halo_growth_rate=rate, **kw)
    kap, bh = o["kappa_edd"], o["bh_mass"]; A = o["A_bh"][:, 1:]; M = bh[:, :-1]
    ratio = np.where(A > 0, A / (kap * M), 0.0); capped = (A > 0) & (ratio > 1)
    return np.median(ratio.max(1)), np.median(capped.sum(1) / np.maximum((A > 0).sum(1), 1)), np.median(bh[:, -1] / bh[:, 0])
print(f"##### M0={M0:.0e}, 20 trees, 801 steps: light-seed regime map. cell = peak Mdot_acc/Mdot_Edd (median over trees) / growth factor G, seed 1e3", flush=True)
Rs = (100.0, 125.0, 150.0, 200.0, 250.0, 300.0)
print("sigma_j \\ R_nuc[pc]  " + "".join(f"{R:>18.0f}" for R in Rs), flush=True)
for s in (0.5, 0.6, 0.7, 0.8, 1.0):
    cells = []
    for R in Rs:
        r, f, G = peak(np.full(N, 1e3), sigma_lnj=s, R_nuc_pc=R)
        cells.append(f"{r:8.2g}/{G:8.3g}")
    print(f"   {s:4.1f}            " + "".join(f"{c:>18s}" for c in cells), flush=True)
print("F(seed) sign changes over a 25-point seed grid (20 trees) in the high-delivery corners:", flush=True)
seeds = np.logspace(0, np.log10(2 * 0.5 * 0.156 * Mh[0, -1]), 25)
for tag, kw in (("fiducial", {}), ("s=1.0,R=100", dict(sigma_lnj=1.0)), ("s=0.5,R=300", dict(R_nuc_pc=300.0)), ("s=1.5,R=300", dict(sigma_lnj=1.5, R_nuc_pc=300.0)),
                ("s=1.5,R=300,eta=0.05", dict(sigma_lnj=1.5, R_nuc_pc=300.0, eta_acc=0.05)), ("s=0.5,R=100,eta=0.5", dict(eta_acc=0.5))):
    Fg = np.array([mvm.target_residual(mvm.run_reservoir_stock(Mh, tz, s_, halo_growth_rate=rate, **kw), 0.5) for s_ in seeds])
    nchg = (np.diff(np.sign(Fg), axis=0) != 0).sum(axis=0)
    print(f"   {tag:24s}: trees with 0/1/2/3+ sign changes = {[int((nchg==k).sum()) for k in (0,1,2)] + [int((nchg>=3).sum())]}", flush=True)
print("DONE", flush=True)
