"""
Quick sensitivity check (not part of the published figure set): how much
does M_seed,crit move under (a) the compaction enhancement Phi_hat, which
every published result in this paper holds fixed at 1 (no enhancement),
and (b) the angular-momentum scatter sigma_j, which the paper's own text
flags as adopted but not separately stress-tested?

Single fiducial halo mass (M_halo(z=5)=5e10 Msun, matching Section 4.1's
quoted fiducial and Figure 4's chi_crit check), N=300 trees, full ensemble
median (not a single representative tree, unlike Figure 3/4's diagnostic
plots) -- consistent in spirit with the eta_acc/epsilon_f stress tests.

Writes compaction_sigmaj_results.json into ./output/.
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

N_TREES = 300
DZ_FIDUCIAL = 0.05
M_RES = 1e4
M0_FIDUCIAL = 5.0e10

t_seed = float(time_at_z(np.array([25.0]))[0])
t_anchor = float(time_at_z(np.array([5.0]))[0])
t_grid = np.linspace(t_seed, t_anchor, 401)
target_z = z_at_time(t_grid)

halo_masses, redshifts, _sm, _mm = pymctrees_adapter.build_forest_for_bin(
    tree_gen, M0_FIDUCIAL, h, N_TREES, z0=5.0, z_max=25.0, m_res_msun=M_RES, dz=DZ_FIDUCIAL, backend="numba", rng_seed=13,
)
M_h_grid = interpolate_tree_onto_grid(halo_masses, redshifts, target_z)
rate_grid = grumpy_halo_growth_rate(halo_masses, redshifts, target_z)


def median_mseedcrit(**overrides):
    kwargs = dict(eta_acc=PAPER_PARAMS.eta_acc, epsilon_sf=PAPER_PARAMS.epsilon_sf, sigma_lnj=PAPER_PARAMS.sigma_lnj)
    kwargs.update(overrides)
    M_crit, never, above = critical_seed_paper(
        M_h_grid, target_z, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=1e8, n_iter=48,
        halo_growth_rate=rate_grid, **kwargs,
    )
    valid = M_crit[~np.isnan(M_crit)]
    return dict(median=float(np.median(valid)), p16=float(np.percentile(valid, 16)),
                p84=float(np.percentile(valid, 84)), n_valid=int(len(valid)), n_total=N_TREES,
                never=int(never.sum()), above_lo=int(above.sum()))


results = {"M0_fiducial": M0_FIDUCIAL, "N_TREES": N_TREES, "compaction_boost": {}, "sigma_lnj": {}}

t0 = time.perf_counter()
for phi in [1.0, 3.0, 10.0]:
    r = median_mseedcrit(compaction_boost=phi)
    results["compaction_boost"][str(phi)] = r
    print(f"compaction_boost={phi:.1f}  median={r['median']:.4e}  [elapsed {time.perf_counter()-t0:.0f}s]", flush=True)

for sj in [0.25, 0.5, 1.0]:
    r = median_mseedcrit(sigma_lnj=sj)
    results["sigma_lnj"][str(sj)] = r
    print(f"sigma_lnj={sj:.2f}  median={r['median']:.4e}  [elapsed {time.perf_counter()-t0:.0f}s]", flush=True)

outpath = OUTDIR / "compaction_sigmaj_results.json"
with open(outpath, "w") as f:
    json.dump(results, f, indent=2)
print("DONE", outpath)
