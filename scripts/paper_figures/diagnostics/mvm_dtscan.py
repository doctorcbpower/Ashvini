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
hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 60, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=i0)
res = {}
print(f"##### M0={M0:.0e}, same 60 trees at each time step", flush=True)
for n in (201, 401, 801, 1601):
    tz = grid(n); Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    Mm, *_ = mvm.critical_seed_stock(Mh, tz, n_iter=36, halo_growth_rate=rate)
    o = mvm.run_reservoir_stock(Mh, tz, Mm, halo_growth_rate=rate)
    zc = 0.5*(tz[1:]+tz[:-1]); dts = np.diff(o["cosmic_time"])
    early = zc > 10
    st_g = np.diff(o["stars_gal"], axis=1); st_n = np.diff(o["stars_nuc"], axis=1)
    res[n] = dict(M=Mm, ms=o["stars_mass"][:, -1], nuc=o["stars_nuc"][:, -1] / o["stars_mass"][:, -1],
                  early_frac=(st_g[:, early].sum(1) + st_n[:, early].sum(1)) / o["stars_mass"][:, -1],
                  G=o["bh_mass"][:, -1] / Mm)
    r = res[n]; b = res.get(401)
    print(f" n={n:5d}: M_seed,crit median {np.nanmedian(Mm):.4e}"
          + (f" (paired ratio to 401: {np.nanmedian(Mm/b['M']):.4f})" if b else "")
          + f"; M*_tot median {np.median(r['ms']):.3e}" + (f" (ratio {np.median(r['ms']/b['ms']):.4f})" if b else "")
          + f"; G_BH {np.nanmedian(r['G']):.3f}; nuclear share {np.median(r['nuc']):.3f}; share of M*_tot formed at z>10: {np.median(r['early_frac']):.3f}", flush=True)
