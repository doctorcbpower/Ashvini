import sys, time
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm
from ashvini.utils import time_at_z, z_at_time
M0 = float(sys.argv[1]); i0 = int(sys.argv[2]); NS = [int(x) for x in sys.argv[3].split(",")]; NT = int(sys.argv[4])
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
cd = cosmo_utils.CosmoData(rp, redshift=[5.0]); tg = ZhangHuiMergerTree(cd, rp, model="cdm")
grid = lambda n: z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], n))
hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, NT, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=i0)
q = lambda x, p: np.nanpercentile(x, p)
T0 = time.time()
for n in NS:
    tz = grid(n); Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz); zc = 0.5 * (tz[1:] + tz[:-1])
    Mm, never, above, F = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, return_F=True)
    o = mvm.run_reservoir_stock(Mh, tz, Mm, halo_growth_rate=rate)
    ms = o["stars_mass"][:, -1]; bh = o["bh_mass"][:, -1]
    print(f"\n===== M0={M0:.0e}, {NT} trees, {n} steps  [{time.time()-T0:.0f}s] =====")
    print(f" M_seed,crit: median {np.nanmedian(Mm):.3e}  [16,84]=[{q(Mm,16):.2e}, {q(Mm,84):.2e}]  (bracket failures: never {int(never.sum())}, above {int(above.sum())}; max|F|={np.nanmax(np.abs(F)):.0e})")
    print(f" M_seed,crit / M_halo(z=5) = {np.nanmedian(Mm/Mh[:,-1]):.2e};  M*_tot/M_halo = {np.median(ms/Mh[:,-1]):.2e};  M_seed,crit/M*_tot = {np.nanmedian(Mm/ms):.3f}")
    print(f" BH growth G_BH = M_BH(z=5)/M_seed at the critical seed: median {np.nanmedian(bh/Mm):.3f} [16,84 = {q(bh/Mm,16):.3f}, {q(bh/Mm,84):.3f}], max {np.nanmax(bh/Mm):.3f}")
    print(f" halo first resolved at z = {np.nanmedian(o['z_first_resolved']):.1f};  M_seed,crit / (f_b M_halo) at that time: median {np.nanmedian(o['seed_over_host_baryons_at_formation']):.2e} [16,84 = {q(o['seed_over_host_baryons_at_formation'],16):.1e}, {q(o['seed_over_host_baryons_at_formation'],84):.1e}]")
    for zt in (20, 15, 10, 7, 5):
        j = int(np.argmin(abs(tz - zt))); r = o["bh_to_host_baryons"][:, j]
        print(f"    z={zt:2d}: M_BH / (f_b M_halo) median {np.nanmedian(r):.2e}  (fraction of trees with the BH heavier than its host's baryons: {np.mean(r[np.isfinite(r)] > 1):.2f}, host resolved in {np.isfinite(r).mean()*100:.0f}% of trees)")
    cum_in = (o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(axis=1)
    print(f" accreted baryons / (f_b M_halo(z=5)) = {np.median(cum_in/(0.156*Mh[:,-1])):.3f} (cold/hot + UV suppression);  stars/accreted {np.median(ms/cum_in):.3f}; gas left {np.median(o['gas_mass'][:,-1]/cum_in):.3f}; ejected: AGN {np.median(o['wind_taken_agn'].sum(1)/cum_in):.3f}, stellar {np.median(o['wind_taken_sn'].sum(1)/cum_in):.3f}; unavailable {np.median(o['gas_unavailable'][:,-1]/cum_in):.1e}")
    print(f" nuclear share of M*_tot {np.median(o['stars_nuc'][:,-1]/ms):.3f}")
    act = o["wind_demand"][:, 1:] > 0
    print(f" feedback-limited (demand > gas): {o['wind_cap'][:,1:][act].mean()*100:.1f}% of {act.sum()} feedback steps;  by z: " + ", ".join(f"{lo}-{hi}: {o['wind_cap'][:,1:][:, (zc<=lo)&(zc>hi)][act[:, (zc<=lo)&(zc>hi)]].mean()*100:.0f}%" for lo, hi in ((25,15),(15,10),(10,7),(7,5))))
    print(f" inflow switched off (hot-mode) in {(o['mdot_in'][:,1:]==0)[Mh[:,1:]>0].mean()*100:.1f}% of steps with an existing halo;  baryon residual max {np.abs(o['baryon_residual']).max():.0e}; limiters {int(o['limiter_nuc'].sum())}/{int(o['limiter_gal'].sum())}")
    if n == NS[-1]:
        seeds = np.logspace(0, np.log10(2 * 0.5 * 0.156 * Mh[0, -1]), 25)
        sub = slice(0, 40)
        Fg = np.array([mvm.target_residual(mvm.run_reservoir_stock(Mh[sub], tz, s, halo_growth_rate=rate[sub]), 0.5) for s in seeds])
        nchg = (np.diff(np.sign(Fg), axis=0) != 0).sum(axis=0)
        print(f" F(seed) sign changes per tree (40 trees, 25-seed grid): {np.bincount(nchg, minlength=4)[:4]} trees with 0/1/2/3")
