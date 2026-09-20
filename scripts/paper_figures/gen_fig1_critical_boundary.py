"""
Generates the tree-ensemble data behind Figure 1 (the critical seed-mass
boundary M_seed,crit(M_halo)) of the High-Redshift-Black-Hole-Growth paper:
N=1000 independent merger trees at each of 25 halo masses spanning
3e10-3e13 Msun, dz=0.05, M_res=1e4 Msun, fiducial accretion parameters.

Requires foraois (https://github.com/doctorcbpower/foraois, or the private
pymctrees predecessor) on sys.path -- see FORAOIS_SRC below -- and this
Ashvini checkout's own venv (ashvini.paper_reservoir and friends).

Writes fig1_critical_boundary_results.json into ./output/.
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
# P(k)/sigma(M) are evaluated at the anchor redshift z0=5, and foraois normalises the collapse
# barrier to the same redshift (delta_col(z) = 1.686 D(5)/D(z)), so the trees are identical to
# those built with P(k) at z=0. Requires foraois with CosmoData.pk_redshift (2026-09-20 fix).
cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[5.0])
tree_gen = ZhangHuiMergerTree(cosmo_data, run_params, model="cdm")

N_TREES = 1000
DZ_FIDUCIAL = 0.05
M_RES = 1e4
MASSES = np.logspace(np.log10(3e10), np.log10(3e13), 25)

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

results = {}
t_start = time.perf_counter()
for M0 in MASSES:
    halo_masses, redshifts, smooth_accretion, merger_mass = pymctrees_adapter.build_forest_for_bin(
        tree_gen, M0, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=1
    )
    M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
    rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

    M_crit, never, above = critical_seed_paper(
        M_h_grid, target_z, f_bh=0.5, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
        halo_growth_rate=rate_grid, **FIDUCIAL,
    )
    valid = M_crit[~np.isnan(M_crit)]
    p16, p50, p84 = (np.percentile(valid, [16, 50, 84]) if len(valid) else (np.nan, np.nan, np.nan))
    results[float(M0)] = dict(
        p16=float(p16), median=float(p50), p84=float(p84),
        n_valid=int(len(valid)), n_total=int(N_TREES),
        never=int(never.sum()), above_lo=int(above.sum()),
    )
    print(f"M0={M0:.3e}  median={results[float(M0)]['median']:.4e}  [elapsed {time.perf_counter()-t_start:.0f}s]", flush=True)

outpath = OUTDIR / "fig1_critical_boundary_results.json"
with open(outpath, "w") as f:
    json.dump(dict(N_TREES=N_TREES, DZ_FIDUCIAL=DZ_FIDUCIAL, M_RES=M_RES, FIDUCIAL=FIDUCIAL, results=results), f, indent=2)

print("DONE.", time.perf_counter() - t_start, "->", outpath)
