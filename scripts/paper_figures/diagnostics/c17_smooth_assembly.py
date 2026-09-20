"""
C17: dependence on the halo-assembly description. The PCH08 comparison is not feasible at the required resolution (see
logs/c17_pch08_check.log); as an independent description of assembly this runs the frozen MVM on the deterministic
Fakhouri, Ma & Boylan-Kolchin (2010) mean accretion history (calibrated at z < ~2, extrapolated here), anchored to
M_halo(z=5) by shooting, and compares with the Zhang & Hui ensemble medians of the production run.
"""
import sys, json, hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from c17_diagnostics import assert_frozen, PROD
import numpy as np
assert_frozen()
from ashvini.paper_reservoir import shoot_smooth_halo_trajectory
from ashvini import reservoir_stock as mvm

D = json.load(open(PROD)); R = {round(np.log10(r["M0"]), 2): r for r in D["results"]}
print("FMBK10 mean accretion history (smooth, deterministic) against the Zhang & Hui ensemble (production medians)")
for M0, key in ((3e10, 10.48), (3e11, 11.48), (3e13, 13.48)):
    Mh, z, t, conv = shoot_smooth_halo_trajectory(M0, z_seed=25.0, z_anchor=5.0, n_steps=800)
    Mh = np.asarray(Mh, dtype=float)
    Mc, never, above, F = mvm.critical_seed_stock(Mh[None, :], z, n_iter=50, return_F=True)
    o = mvm.run_reservoir_stock(Mh[None, :], z, Mc)
    cum_in = float((o["mdot_in"][0, 1:] * np.diff(o["cosmic_time"])).sum())
    r = R[key]; zh = np.array(r["z"]); mh_rep = np.array(r["rep_halo"])
    at = lambda zz: Mh[int(np.argmin(abs(z - zz)))]
    zt = np.nanmedian(r["Mcrit"])
    print(f"M0={M0:.0e}: smooth M_halo(z=10)/M0 = {at(10)/M0:.3f}, (z=15)/M0 = {at(15)/M0:.4f}, halo mass at z_seed=25: {Mh[0]:.2e} Msun (converged={conv}); "
          f"M_seed,crit(smooth) = {float(Mc[0]):.3e} vs ZH median {zt:.3e} (ratio {float(Mc[0])/zt:.3f}); M*_tot/M_halo = {o['stars_mass'][0,-1]/M0:.2e} vs ZH {np.nanmedian(np.array(r['Mstar'])/M0):.2e}; "
          f"accreted/(f_b Mh) = {cum_in/(0.156*M0):.3f} vs ZH {np.median(r['accreted_over_fb_Mhalo']):.3f}; G_BH {float(o['bh_mass'][0,-1]/Mc[0]):.3f}; |F|={abs(float(F[0])):.0e}")
