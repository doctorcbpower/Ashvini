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
tz = z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], 801))
hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 60, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=i0)
Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
def run(**kw):
    Mm, never, above, F = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, return_F=True, **kw)
    o = mvm.run_reservoir_stock(Mh, tz, Mm, halo_growth_rate=rate, **kw)
    return Mm, o, int(never.sum() + above.sum())
base, ob, _ = run()
print(f"##### MVM sensitivities, M0={M0:.0e}, 60 paired trees, 801 steps.  base median M_seed,crit = {np.nanmedian(base):.3e}  (M*_tot {np.median(ob['stars_mass'][:,-1]):.2e}, G_BH {np.nanmedian(ob['bh_mass'][:,-1]/base):.3f})", flush=True)
print(f"  {'variation':28s} {'median M_seed,crit':>18s} {'ratio to base (paired)':>26s} {'M*_tot ratio':>13s} {'G_BH':>6s}", flush=True)
CFG = [("sigma_j 0.5 -> 0.75", dict(sigma_lnj=0.75)), ("eps_sf x0.5", dict(epsilon_sf=0.0075)), ("eps_sf x2", dict(epsilon_sf=0.03)),
       ("R_nuc = 50 pc", dict(R_nuc_pc=50.0)), ("R_nuc = 200 pc", dict(R_nuc_pc=200.0)), ("f_mom = 0.3", dict(f_mom=0.3)),
       ("eta_SN x0.5", dict(eta_sn_scale=0.5)), ("eta_SN x2", dict(eta_sn_scale=2.0)), ("stellar wind off", dict(eta_sn_scale=0.0))]
for tag, kw in CFG:
    Mm, o, bad = run(**kw)
    r = Mm / base
    print(f"  {tag:28s} {np.nanmedian(Mm):18.3e}   {np.nanmedian(r):8.3f} [{np.nanmin(r):.3f}, {np.nanmax(r):.3f}]   {np.median(o['stars_mass'][:,-1]/ob['stars_mass'][:,-1]):13.3f} {np.nanmedian(o['bh_mass'][:,-1]/Mm):6.3f}" + (f"  BRACKET FAIL {bad}" if bad else ""), flush=True)
