"""
Generates the data behind Figure 4 (Appendix): chi_crit reconciliation --
f_BH = M_BH/M_star at z=5 as a function of seed mass, at the paper's
quoted fiducial halo mass (M_halo(z=5)=5e10 Msun), for
chi_crit in {1, 3, 10, 1e2, STRICT_EDDINGTON_CHI_CRIT}, each with its own
M_seed,crit boundary.

A single representative tree (closest individual M_seed,crit to the
ensemble median under strict-Eddington growth) is picked once and reused
for every chi_crit curve, so the only thing varying between curves is
chi_crit itself.

Writes fig4_chicrit_results.json into ./output/.
"""
import sys
import json
from pathlib import Path

FORAOIS_SRC = "/Users/00075868/MyCodes/foraois/src"
FORAOIS_CONFIG = "/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"

sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree
from foraois.utils import io as pymctrees_io
from ashvini import pymctrees_adapter
from ashvini.paper_reservoir import (
    interpolate_tree_onto_grid, grumpy_halo_growth_rate,
    critical_seed_paper, run_reservoir_paper,
)
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

N_TREES = 200
DZ_FIDUCIAL = 0.05
M_RES = 1e4
M0_FIDUCIAL = 5.0e10  # matches the paper's own quoted fiducial (Section 4.1)

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

halo_masses, redshifts, _sm, _mm = pymctrees_adapter.build_forest_for_bin(
    tree_gen, M0_FIDUCIAL, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=7,
)
M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

# pick the representative tree: closest individual M_seed,crit to the ensemble
# median under strict-Eddington growth (same convention used for
# gen_fig3_trajectories.py), so all chi_crit curves below share one fixed
# halo trajectory.
M_crit_strict, never, above = critical_seed_paper(
    M_h_grid, target_z, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
    halo_growth_rate=rate_grid, **FIDUCIAL,
)
valid = ~np.isnan(M_crit_strict)
median_val = np.median(M_crit_strict[valid])
idx = np.nanargmin(np.abs(M_crit_strict - median_val))

M_h_rep = M_h_grid[idx:idx + 1]
rate_rep = rate_grid[idx:idx + 1]

CHI_CRIT_VALUES = [1.0, 3.0, 10.0, 1.0e2, PAPER_PARAMS.strict_eddington_chi_crit]
SEED_MASSES = np.logspace(0, 8, 61)

results = {"M0_fiducial": M0_FIDUCIAL, "representative_tree_M_seed_crit_strict": float(M_crit_strict[idx]),
           "seed_masses": SEED_MASSES.tolist(), "curves": {}, "boundaries": {}}

for chi_crit in CHI_CRIT_VALUES:
    f_bh_z5 = []
    for m_seed in SEED_MASSES:
        out = run_reservoir_paper(
            M_h_rep, target_z, m_seed, chi_crit=chi_crit,
            halo_growth_rate=rate_rep, **FIDUCIAL,
        )
        bh_final = out["bh_mass"][0, -1]
        star_final = out["stars_mass"][0, -1]
        f_bh_z5.append(bh_final / star_final)
    results["curves"][str(chi_crit)] = f_bh_z5

    M_crit_this, never_this, above_this = critical_seed_paper(
        M_h_rep, target_z, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
        chi_crit=chi_crit, halo_growth_rate=rate_rep, **FIDUCIAL,
    )
    results["boundaries"][str(chi_crit)] = float(M_crit_this[0])
    print(f"chi_crit={chi_crit:.3e}  M_seed,crit={M_crit_this[0]:.4e}", flush=True)

outpath = OUTDIR / "fig4_chicrit_results.json"
with open(outpath, "w") as f:
    json.dump(results, f, indent=2)
print("DONE", outpath)
