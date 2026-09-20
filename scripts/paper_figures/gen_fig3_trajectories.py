"""
Generates the data behind Figure 3 (Appendix): representative M_halo,
M_gas, M_star, M_BH trajectories from z_seed=25 to z_anchor=5, one tree per
halo mass (low / mid, below M_hot / high, above M_hot), each seeded at its
own strict-Eddington critical mass and then integrated under the fiducial
graded cap (chi_crit=10). Also classifies each timestep into the S
(supply-rich/throttled) / E (Eddington-limited) / G (gas-supply-limited)
growth regime from black_holes_growth_slimdisk.bh_growth_step_slimdisk's
own M1/M2 boundaries, for the background shading.

Writes fig3_trajectories_results.json into ./output/.
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
    critical_seed_paper, run_reservoir_paper, M_HOT_FIDUCIAL,
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
MASSES = [3.0e10, 3.0e11, 3.0e13]  # low / mid (below M_hot) / high (above M_hot)

FIDUCIAL = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

output = {"M_hot": M_HOT_FIDUCIAL, "masses": {}}

for i, M0 in enumerate(MASSES):
    halo_masses, redshifts, _sm, _mm = pymctrees_adapter.build_forest_for_bin(
        tree_gen, M0, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=11 + i,
    )
    M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
    rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)

    # representative tree = closest individual M_seed,crit to the ensemble
    # median, under strict-Eddington growth (same convention as the
    # chi_crit-reconciliation figure, gen_fig4_chicrit.py).
    M_crit_strict, never, above = critical_seed_paper(
        M_h_grid, target_z, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
        halo_growth_rate=rate_grid, **FIDUCIAL,
    )
    valid = ~np.isnan(M_crit_strict)
    median_val = np.median(M_crit_strict[valid])
    idx = int(np.nanargmin(np.abs(M_crit_strict - median_val)))

    M_h_rep = M_h_grid[idx:idx + 1]
    rate_rep = rate_grid[idx:idx + 1]
    seed_mass_rep = float(M_crit_strict[idx])

    # now integrate under the FIDUCIAL graded cap (chi_crit=10), seeded at
    # this tree's own strict-Eddington critical mass.
    out = run_reservoir_paper(
        M_h_rep, target_z, seed_mass_rep, chi_crit=PAPER_PARAMS.chi_crit,
        halo_growth_rate=rate_rep, **FIDUCIAL,
    )

    kappa_edd = float(out["kappa_edd"])
    A_bh = out["A_bh"][0]
    bh_mass = out["bh_mass"][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        M1 = A_bh / (kappa_edd * PAPER_PARAMS.chi_crit)
        M2 = A_bh / kappa_edd
    regime = np.full(len(bh_mass), "S", dtype="<U1")
    regime[(bh_mass >= M1) & (bh_mass < M2)] = "E"
    regime[bh_mass >= M2] = "G"
    regime[A_bh <= 0] = "S"

    output["masses"][str(M0)] = dict(
        M0=M0, seed_mass=seed_mass_rep, kappa_edd=kappa_edd,
        redshift=target_z.tolist(),
        halo_mass=out["halo_mass"][0].tolist(),
        gas_mass=out["gas_mass"][0].tolist(),
        stars_mass=out["stars_mass"][0].tolist(),
        bh_mass=bh_mass.tolist(),
        regime=regime.tolist(),
    )
    print(f"M0={M0:.3e}  seed_mass={seed_mass_rep:.4e}  final bh={bh_mass[-1]:.4e}", flush=True)

outpath = OUTDIR / "fig3_trajectories_results.json"
with open(outpath, "w") as f:
    json.dump(output, f, indent=2)
print("DONE", outpath)
