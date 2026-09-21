import sys
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate
from ashvini import reservoir_stock as mvm
from ashvini import reservoir_stock_premvm as pre
from ashvini.utils import time_at_z, z_at_time
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
cd = cosmo_utils.CosmoData(rp, redshift=[5.0]); tg = ZhangHuiMergerTree(cd, rp, model="cdm")
tz = z_at_time(np.linspace(time_at_z(np.array([25.0]))[0], time_at_z(np.array([5.0]))[0], 401))
q = lambda x, p: np.nanpercentile(x, p)
for i, M0 in enumerate((3e10, 3e11)):
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 20, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=31 + i)
    Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
    print(f"\n################ M0={M0:.0e}, the same 20 trees for both models ################")
    # ---------- pre-MVM reference (frozen), current defaults ----------
    Mp, *_ = pre.critical_seed_stock(Mh, tz, seed_mass_lo=1e-2, seed_mass_hi=1e10, n_iter=36, halo_growth_rate=rate)
    op = pre.run_reservoir_stock(Mh, tz, Mp, chi_crit=1e12, halo_growth_rate=rate)
    # ---------- MVM ----------
    Mm, never, above, Fm = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, return_F=True)
    o = mvm.run_reservoir_stock(Mh, tz, Mm, halo_growth_rate=rate)
    print(f"  median M_seed,crit:  pre-MVM {np.nanmedian(Mp):.3e} [16,84 = {q(Mp,16):.2e}, {q(Mp,84):.2e}]   MVM {np.nanmedian(Mm):.3e} [{q(Mm,16):.2e}, {q(Mm,84):.2e}]   (never {int(never.sum())}, above {int(above.sum())}, max|F|={np.nanmax(np.abs(Fm)):.1e})")
    print(f"  per-tree ratio MVM/pre-MVM: median {np.nanmedian(Mm/Mp):.3f}, range [{np.nanmin(Mm/Mp):.3f}, {np.nanmax(Mm/Mp):.3f}]")
    gr = o["bh_mass"][:, -1] / Mm
    print(f"  BH growth at the critical seed (MVM): median x{np.nanmedian(gr):.3f}  [16,84 = {q(gr,16):.3f}, {q(gr,84):.3f}];  pre-MVM x{np.nanmedian(op['bh_mass'][:,-1]/Mp):.3f}")
    ms = o["stars_mass"][:, -1]; print(f"  M*_tot/M_halo = {np.median(ms/Mh[:,-1]):.2e};  M_seed,crit/M*_tot = {np.median(Mm/ms):.3f}  (identity check: M_BH/(0.5 M*) - 1 = {np.nanmedian(np.abs(o['bh_mass'][:,-1]/(0.5*ms)-1)):.1e})")
    dsn, dsg = np.diff(o["stars_nuc"], axis=1), np.diff(o["stars_gal"], axis=1); tot = dsn + dsg; ok = tot > 0
    sh = np.where(ok, dsn / np.where(ok, tot, 1.0), np.nan)
    print(f"  nuclear star formation: time-integrated share of M*_tot = {np.median(o['stars_nuc'][:,-1]/ms):.3f}; per-step share median {np.nanmedian(sh):.3f}, 25/75th pct {q(sh,25):.3f}/{q(sh,75):.3f};  SFR_nuc/SFR_gal (time-integrated) = {np.median(o['stars_nuc'][:,-1]/o['stars_gal'][:,-1]):.3f}")
    alive = Mh[:, 1:] > 0
    has_gas = ((o["gas_gal"][:, :-1] + o["gas_nuc"][:, :-1]) > 0) & alive
    cap = o["wind_cap"][:, 1:] & alive
    print(f"  wind cap (feedback demand > available gas): {cap.sum()} of {has_gas.sum()} steps with gas = {cap.sum()/max(has_gas.sum(),1)*100:.1f}%;  limiter_nuc {int((o['limiter_nuc'][:,1:]&alive).sum())} steps, limiter_gal {int((o['limiter_gal'][:,1:]&alive).sum())} steps")
    cum_in = (o["mdot_in"][:, 1:] * np.diff(o["cosmic_time"])).sum(axis=1)
    agn, sne = o["wind_taken_agn"].sum(axis=1), o["wind_taken_sn"].sum(axis=1)
    print(f"  baryon budget at z=5 (fraction of total accreted mass): stars {np.median(ms/cum_in):.3f}, BH gain {np.median((o['bh_mass'][:,-1]-Mm)/cum_in):.4f}, unavailable {np.median(o['gas_unavailable'][:,-1]/cum_in):.1e}, remaining gas {np.median((o['gas_gal'][:,-1]+o['gas_nuc'][:,-1])/cum_in):.3f}, ejected by AGN {np.median(agn/cum_in):.3f}, by stellar wind {np.median(sne/cum_in):.3f}")
    print(f"  ejected mass ratio AGN : stellar = {np.median(agn/np.maximum(sne,1)):.3f};  baryon residual max {np.abs(o['baryon_residual']).max():.1e}")
    # zero seed sensitivity to nothing: F monotonic in seed?
    seeds = np.logspace(0, np.log10(2*0.5*0.156*Mh[0,-1]), 25)
    Fg = np.array([mvm.target_residual(mvm.run_reservoir_stock(Mh, tz, s, halo_growth_rate=rate), 0.5) for s in seeds])   # (seeds, trees)
    nchg = (np.diff(np.sign(Fg), axis=0) != 0).sum(axis=0)
    print(f"  F(seed) sign changes per tree over a 25-point seed grid: {np.bincount(nchg, minlength=4)[:4]} trees with 0/1/2/3 sign changes")
