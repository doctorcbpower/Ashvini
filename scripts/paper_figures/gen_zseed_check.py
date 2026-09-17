"""
Quick sensitivity check (not part of the published figure set): how much
does M_seed,crit move under the seed-formation redshift z_seed, held
fixed at 25 throughout the paper despite the true maximum available
window running from z=50 (Section 2, "approximately 23 e-folding times
... this is the maximum time interval physically available").

Single fiducial halo mass (M_halo(z=5)=5e10 Msun), N=300 trees built back
to z_max=55 (a small margin past the most extreme z_seed=50 tested, so
interpolate_tree_onto_grid never needs to extrapolate past the tree's own
range), full ensemble median -- same convention as the compaction/sigma_j
check.

Writes zseed_results.json into ./output/.
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

N_TREES = 300
DZ_FIDUCIAL = 0.05
M_RES = 1e4
M0_FIDUCIAL = 5.0e10
Z_MAX_TREE = 55.0
Z_SEED_VALUES = [50.0, 40.0, 25.0, 15.0]

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)

t0 = time.perf_counter()
halo_masses, redshifts, _sm, _mm = pymctrees_adapter.build_forest_for_bin(
    tree_gen, M0_FIDUCIAL, h, N_TREES, z0=5.0, z_max=Z_MAX_TREE, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=17,
)
print(f"trees built [elapsed {time.perf_counter()-t0:.0f}s]", flush=True)

t_anchor = float(time_at_z(np.array([5.0]))[0])

results = {"M0_fiducial": M0_FIDUCIAL, "N_TREES": N_TREES, "z_max_tree": Z_MAX_TREE, "z_seed": {}}

for z_seed in Z_SEED_VALUES:
    t_seed = float(time_at_z(np.array([z_seed]))[0])
    t_grid = np.linspace(t_seed, t_anchor, 401)
    target_z = z_at_time(t_grid)

    M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
    rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

    M_crit, never, above = critical_seed_paper(
        M_h_grid, target_z, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
        halo_growth_rate=rate_grid, **FIDUCIAL,
    )
    valid = M_crit[~np.isnan(M_crit)]
    r = dict(median=float(np.median(valid)), p16=float(np.percentile(valid, 16)),
             p84=float(np.percentile(valid, 84)), n_valid=int(len(valid)), n_total=N_TREES,
             never=int(never.sum()), above_lo=int(above.sum()))
    results["z_seed"][str(z_seed)] = r
    print(f"z_seed={z_seed:.0f}  median={r['median']:.4e}  n_valid={r['n_valid']}/{N_TREES}  [elapsed {time.perf_counter()-t0:.0f}s]", flush=True)

outpath = OUTDIR / "zseed_results.json"
with open(outpath, "w") as f:
    json.dump(results, f, indent=2)
print("DONE", outpath)
