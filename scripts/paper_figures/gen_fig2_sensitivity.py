"""
Generates the tree-ensemble data behind Figure 2 (the three-panel
eta_acc / epsilon_f / spin sensitivity figure) of the
High-Redshift-Black-Hole-Growth paper: N=500 trees at each of 15 halo
masses spanning 3e10-3e13 Msun, fiducial plus six parameter-variation
curves (eta_acc lo/hi, epsilon_f lo/hi, spin chaotic/coherent).

See gen_fig1_critical_boundary.py for the shared setup (foraois path,
config, tree generator). Writes fig2_sensitivity_results.json into ./output/.
"""
import sys
import time
import json
from pathlib import Path

FORAOIS_SRC = "/Users/00075868/MyCodes/foraois/src"
FORAOIS_CONFIG = "/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"

sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io as pymctrees_io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import interpolate_tree_onto_grid, grumpy_halo_growth_rate, critical_seed_paper
from ashvini.paper_reservoir_params import PAPER_PARAMS
from ashvini.utils import time_at_z, z_at_time

OUTDIR = Path(__file__).parent / "output"
OUTDIR.mkdir(exist_ok=True)

run_params = pymctrees_io.get_params(FORAOIS_CONFIG)
h = run_params["Cosmology"]["h"]
cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[5.0])
tree_gen = ZhangHuiMergerTree(cosmo_data, run_params, model="cdm")

N_TREES = 500
DZ = 0.05
M_RES = 1e4
MASSES = np.logspace(np.log10(3e10), np.log10(3e13), 15)

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)
EPSILON_FID = PAPER_PARAMS.epsilon

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

CURVES = {
    "fiducial":      dict(**FIDUCIAL, epsilon=EPSILON_FID),
    "eta_acc_lo":    dict(**{**FIDUCIAL, "eta_acc": 1e-3}, epsilon=EPSILON_FID),
    "eta_acc_hi":    dict(**{**FIDUCIAL, "eta_acc": 1e-1}, epsilon=EPSILON_FID),
    "epsilon_f_lo":  dict(**FIDUCIAL, epsilon=EPSILON_FID, epsilon_f=5e-5),
    "epsilon_f_hi":  dict(**FIDUCIAL, epsilon=EPSILON_FID, epsilon_f=5e-3),
    "spin_chaotic":  dict(**FIDUCIAL, epsilon=0.057),
    "spin_coherent": dict(**FIDUCIAL, epsilon=0.32),
}

results = {name: {} for name in CURVES}
t_start = time.perf_counter()

for M0 in MASSES:
    halo_masses, redshifts, smooth_accretion, merger_mass = pymctrees_adapter.build_forest_for_bin(
        tree_gen, M0, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ, backend="numba", rng_seed=1
    )
    M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
    rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

    for name, overrides in CURVES.items():
        M_crit, never, above = critical_seed_paper(
            M_h_grid, target_z, f_bh=0.5, seed_mass_lo=1e-2, seed_mass_hi=1e10, n_iter=46,
            halo_growth_rate=rate_grid, **overrides,
        )
        valid = M_crit[~np.isnan(M_crit)]
        p16, p50, p84 = (np.percentile(valid, [16, 50, 84]) if len(valid) else (np.nan, np.nan, np.nan))
        results[name][float(M0)] = dict(
            p16=float(p16), median=float(p50), p84=float(p84),
            n_valid=int(len(valid)), n_total=int(N_TREES),
            never=int(never.sum()), above_lo=int(above.sum()),
        )
    print(f"M0={M0:.3e}  fiducial_median={results['fiducial'][float(M0)]['median']:.4e}  "
          f"[elapsed {time.perf_counter()-t_start:.0f}s]", flush=True)

outpath = OUTDIR / "fig2_sensitivity_results.json"
with open(outpath, "w") as f:
    json.dump(dict(N_TREES=N_TREES, DZ=DZ, M_RES=M_RES, MASSES=MASSES.tolist(),
                    CURVES=CURVES, results=results), f, indent=2)

print("DONE. Total time:", time.perf_counter() - t_start, "->", outpath)
