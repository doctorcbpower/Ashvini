import sys, time
sys.path.insert(0, "/Users/00075868/MyCodes/foraois/src")
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm
from ashvini.utils import time_at_z, z_at_time
M0 = float(sys.argv[1]); i0 = int(sys.argv[2])
rp = io.get_params("/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"); h = rp["Cosmology"]["h"]
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
print(f"##### MVM tree-resolution follow-up, M0={M0:.0e}: 240 trees per setting, two independent sets each", flush=True)
res = {}
for dz, mres in ((0.1, 1e4), (0.05, 1e4), (0.025, 1e4), (0.05, 1e5), (0.05, 1e3)):
    meds = []
    for rep in range(2):
        a, b = trees(240, dz=dz, mres=mres, seed=i0 + 200 + rep); Ms, bad, Fm = msc(a, b); meds.append(np.nanmedian(Ms))
    res[(dz, mres)] = meds
    print(f"[{time.time()-T0:5.0f}s] dz={dz:<6} M_res={mres:.0e}: medians {meds[0]:.4e}, {meds[1]:.4e}  mean {np.mean(meds):.4e}  (rep-to-rep {abs(np.log10(meds[0]/meds[1])):.4f} dex)", flush=True)
ref = np.mean(res[(0.05, 1e4)])
print("ratios to the (dz=0.05, M_res=1e4) mean:", {f"dz={k[0]},Mres={k[1]:.0e}": round(float(np.mean(v)/ref), 3) for k, v in res.items()}, flush=True)
print("DONE", flush=True)
