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
grid = lambda n: z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], n))
def trees(N, dz=0.05, mres=1e4, seed=11):
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, N, z0=5.0, z_max=25.0, m_res_msun=1e4 if mres is None else mres, dz=dz, backend="numba", rng_seed=seed)
    return hm, zz
def msc(hm, zz, n=401, fb=0.5):
    tz = grid(n)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    Ms, ne, ab, F = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, f_bh=fb, return_F=True)
    return Ms, int(ne.sum() + ab.sum()), np.nanmax(np.abs(F))
T0 = time.time()
def log(tag, Ms, bad, Fmax, base=None):
    r = f"  per-tree ratio to base: median {np.nanmedian(Ms/base):.4f} [min {np.nanmin(Ms/base):.4f}, max {np.nanmax(Ms/base):.4f}]" if base is not None else ""
    print(f"[{time.time()-T0:5.0f}s] {tag:34s} median {np.nanmedian(Ms):.4e} [16,84]=[{np.nanpercentile(Ms,16):.2e},{np.nanpercentile(Ms,84):.2e}] bracket-fail {bad} max|F|={Fmax:.0e}{r}", flush=True)
print(f"##### MVM convergence, M0={M0:.0e}  (fixed R_nuc=100 pc, n_rd=10, c=4, sigma_j=0.5, momentum AGN f_mom=1, f_BH=0.5; base: 60 trees, 401 steps, dz=0.05, M_res=1e4)", flush=True)
hm, zz = trees(240, seed=i0)
print("== A. number of trees: disjoint sets of 60, then pooled ==", flush=True)
meds = []
for k in range(4):
    Ms, bad, Fm = msc(hm[60*k:60*(k+1)], zz); meds.append(np.nanmedian(Ms)); log(f"set {k}", Ms, bad, Fm)
    if k == 0: base60 = Ms
print(f"   std of log10(median) across the four 60-tree sets: {np.std(np.log10(meds)):.4f} dex", flush=True)
Ms, bad, Fm = msc(hm[:120], zz); log("pooled 120", Ms, bad, Fm)
Ms, bad, Fm = msc(hm, zz); log("pooled 240", Ms, bad, Fm)
h60 = hm[:60]
print("== B. time step (uniform in cosmic time, z=25->5), same 60 trees ==", flush=True)
for n in (201, 401, 801):
    Ms, bad, Fm = msc(h60, zz, n=n); log(f"n_steps={n}", Ms, bad, Fm, base=base60)
print("== C. tree resolution (new tree sets, N=60, same seed) ==", flush=True)
for dz, mres in ((0.05, 1e4), (0.1, 1e4), (0.025, 1e4), (0.05, 1e5), (0.05, 1e3)):
    a, b = trees(60, dz=dz, mres=mres, seed=i0 + 100); Ms, bad, Fm = msc(a, b)
    if (dz, mres) == (0.05, 1e4): ref = np.nanmedian(Ms)
    log(f"dz={dz}, M_res={mres:.0e}  (x{np.nanmedian(Ms)/ref:.3f} of ref)", Ms, bad, Fm)
print("DONE", flush=True)
