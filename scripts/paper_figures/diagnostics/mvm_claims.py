import sys
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
hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(tg, M0, h, 60, z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=i0)
Mh = interpolate_tree_onto_grid(hm, zz, tz); rate = grumpy_halo_growth_rate(hm, zz, tz)
def crit(**kw):
    Mm, *_ = mvm.critical_seed_stock(Mh, tz, n_iter=40, halo_growth_rate=rate, **kw); return Mm
base = crit()
print(f"##### M0={M0:.0e}: base median M_seed,crit = {np.nanmedian(base):.3e}")
# ---- (a) is the Eddington limit binding? ----
def regime(seed, label, **kw):
    o = mvm.run_reservoir_stock(Mh, tz, seed, halo_growth_rate=rate, **kw)
    kap = o["kappa_edd"]; bh = o["bh_mass"]; A = o["A_bh"][:, 1:]; dt = np.diff(o["cosmic_time"])
    edd = (kap * bh[:, :-1] < A) & (A > 0)          # supply exceeds kappa*M: Eddington-limited
    grow = np.diff(bh, axis=1)
    fm = (grow * edd).sum(1) / np.maximum(grow.sum(1), 1e-300)
    efold = np.log(bh[:, -1] / bh[:, 0]) / (kap * (o["cosmic_time"][-1] - o["cosmic_time"][0]) / 1.0)
    print(f"  {label:22s}: final/seed x{np.nanmedian(bh[:,-1]/bh[:,0]):10.3f};  Eddington-limited steps {edd[grow>0].mean()*100:5.1f}% of growth steps, {np.nanmedian(fm)*100:5.1f}% of the BH mass gained;  e-folds used = {np.nanmedian(np.log(bh[:,-1]/bh[:,0])):6.3f} of {kap*(o['cosmic_time'][-1]-o['cosmic_time'][0]):.1f} available;  M_BH/M*_tot(z=5) = {np.nanmedian(bh[:,-1]/o['stars_mass'][:,-1]):.2e}")
    return o
print("(a) regime and achieved M_BH/M*_tot for fixed seeds:")
regime(base, "critical seed")
for s, lab in ((1e2, "seed 1e2 (Pop III)"), (1e3, "seed 1e3 (runaway)"), (2e5, "seed 2e5 (direct collapse)"), (1e7, "seed 1e7")):
    regime(np.full(Mh.shape[0], s), lab)
# ---- (b) the paper's headline sensitivities, which were not in the pruned MVM set ----
print("(b) paired ratios of M_seed,crit to base:")
for tag, kw in (("eta_acc x0.1", dict(eta_acc=0.0005)), ("eta_acc x10", dict(eta_acc=0.05)), ("eps_rad 0.057", dict(epsilon=0.057)), ("eps_rad 0.32", dict(epsilon=0.32)),
                ("f_BH = 0.1", dict(f_bh=0.1)), ("f_BH = 0.9", dict(f_bh=0.9))):
    Mm = crit(**kw); r = Mm / base if "f_bh" not in kw else Mm / base
    print(f"  {tag:16s}: median {np.nanmedian(Mm):.3e}, ratio {np.nanmedian(r):.3f} [{np.nanmin(r):.3f}, {np.nanmax(r):.3f}]", flush=True)
