"""
Generates the data behind Figure 5 (Appendix): critical seed-mass boundary
at f_BH targets 0.1, 0.5 (fiducial), and 0.9, each computed independently
from a full N=1000-tree, 25-point merger-tree ensemble spanning
M_halo(z=5)=3e10-3e13 Msun (the same grid as Figure 1). Tree ensembles are
built once per halo mass and reused across all three f_BH targets.

Writes fig5_fbhsensitivity_results.json into ./output/. Takes ~10 minutes.
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

N_TREES = 1000
DZ_FIDUCIAL = 0.05
M_RES = 1e4
MASSES = np.logspace(np.log10(3e10), np.log10(3e13), 25)
F_BH_TARGETS = [0.1, 0.5, 0.9]

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

results = {str(f): {} for f in F_BH_TARGETS}
t_start = time.perf_counter()
for M0 in MASSES:
    halo_masses, redshifts, _sm, _mm = pymctrees_adapter.build_forest_for_bin(
        tree_gen, M0, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=1,
    )
    M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
    rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

    for f_bh in F_BH_TARGETS:
        M_crit, never, above = critical_seed_paper(
            M_h_grid, target_z, f_bh=f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
            halo_growth_rate=rate_grid, **FIDUCIAL,
        )
        valid = M_crit[~np.isnan(M_crit)]
        p16, p50, p84 = (np.percentile(valid, [16, 50, 84]) if len(valid) else (np.nan, np.nan, np.nan))
        results[str(f_bh)][float(M0)] = dict(
            p16=float(p16), median=float(p50), p84=float(p84),
            n_valid=int(len(valid)), n_total=int(N_TREES),
        )
        print(f"M0={M0:.3e}  f_bh={f_bh}  median={p50:.4e}  [elapsed {time.perf_counter()-t_start:.0f}s]", flush=True)

outpath = OUTDIR / "fig5_fbhsensitivity_results.json"
with open(outpath, "w") as f:
    json.dump(dict(N_TREES=N_TREES, DZ_FIDUCIAL=DZ_FIDUCIAL, M_RES=M_RES, FIDUCIAL=FIDUCIAL, F_BH_TARGETS=F_BH_TARGETS, results=results), f, indent=2)

print("DONE.", time.perf_counter() - t_start, "->", outpath)
