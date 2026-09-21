import sys
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm
from ashvini.utils import time_at_z, z_at_time
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
cd = cosmo_utils.CosmoData(rp, redshift=[5.0]); tg = ZhangHuiMergerTree(cd, rp, model="cdm")
tz = z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], 401)); zc = 0.5*(tz[1:]+tz[:-1])
for i, M0 in enumerate((3e10, 3e11)):
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 20, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=41 + i)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    Mm, *_ = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate)
    o = mvm.run_reservoir_stock(Mh, tz, Mm, halo_growth_rate=rate)
    W = o["wind_demand"][:, 1:]; cap = o["wind_cap"][:, 1:]; act = W > 0
    print(f"\nM0={M0:.0e}: steps with any feedback demand: {act.sum()};  of these the demand exceeds the available gas (reservoirs emptied) in {cap.sum()} = {cap.sum()/act.sum()*100:.1f}%")
    dom = np.where(act, o["wind_taken_agn"][:, 1:] / np.where(act, o["wind_taken_agn"][:, 1:] + o["wind_taken_sn"][:, 1:], 1), np.nan)
    print(f"   AGN share of the ejected mass per step: median {np.nanmedian(dom):.2f}; steps where AGN demand alone exceeds stellar: {int((o['wind_taken_agn'][:,1:] > o['wind_taken_sn'][:,1:])[act].sum())} of {act.sum()}")
    for lo, hi in ((25, 15), (15, 10), (10, 7), (7, 5)):
        s = (zc <= lo) & (zc > hi)
        a, c_ = act[:, s], cap[:, s]
        print(f"   z {lo:2d}->{hi:2d}: steps with feedback {a.sum():4d}, demand > gas in {c_[a].mean()*100 if a.sum() else float('nan'):5.1f}%")
    # when does the galaxy hold gas at the END of a step?  fraction of steps with exactly zero galaxy+nuclear gas after the step
    zero = ((o["gas_gal"][:, 1:] + o["gas_nuc"][:, 1:]) == 0) & (Mh[:, 1:] > 0) & (o["mdot_in"][:, 1:] > 0)
    print(f"   steps (halo accreting) that END with both reservoirs exactly empty: {zero.sum()} of {int(((Mh[:,1:]>0)&(o['mdot_in'][:,1:]>0)).sum())}")
