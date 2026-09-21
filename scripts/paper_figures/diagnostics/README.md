# MVM diagnostics (2026-09-20)

Scripts and logs behind the numbers in `docs/mvm_numerical_crib_sheet.md` that do **not** come from the
production ensemble (`output/mvm_production_results.json`). The merger-tree sampler (numba backend) is not
seed-reproducible, so every log is the record of its own tree ensemble; two logs of "the same" experiment
differ by sampling noise (about 1 to 2 per cent in a 60-tree median). All runs use the frozen MVM
(`ashvini/reservoir_stock.py`) at the fiducial parameters unless the script varies them, and f_BH = 0.5.

| Script | Log(s) in `logs/` | Trees / steps | What it tests |
|---|---|---|---|
| `mvm_converge.py` | `converge_60trees_*` | 60 (240 pooled), 201/401/801, dz, M_res | sampling, time step, first look at dz and M_res (60 trees) |
| `mvm_converge2.py` | `converge_dz_Mres_240trees_*` | 240 x 2 sets, 401 | dz and M_res with 240 trees per setting |
| `mvm_dtscan.py` | `timestep_paired_60trees_*` | 60 paired, 201/401/801/1601 | time-step convergence; nuclear share, early star formation, G_BH |
| `mvm_mass.py` | `mass_scaling_240trees_*` | 240, 801 (3e13 also 401) | mass scaling, host-baryon fractions, budget, feedback-limited fractions |
| `mvm_sens.py` | `sensitivities_60trees_*` | 60 paired, 801 | sigma_j, eps_sf, R_nuc, f_mom, eta_SN, wind off |
| `mvm_claims.py` | `claims_regime_and_sens_60trees_*` | 60 paired, 801 | Eddington regime for fixed seeds; eta_acc, radiative efficiency, f_BH |
| `mvm_regime.py` | `regime_grid_20trees_*` | 20, 801 | (sigma_j, R_nuc) grid and eta_acc series: peak Mdot_acc/Mdot_Edd, capped fraction, G |
| `mvm_regime2.py` | `regime_lightseed_map_20trees_*` | 20, 801 | fine light-seed map; F(seed) sign changes in accessible corners |
| `mvm_compare.py` | `premvm_vs_mvm_20trees` | 20 paired, 401 | frozen pre-MVM against MVM on identical trees |
| `mvm_cap.py` | `feedback_limited_by_z_20trees` | 20, 401 | feedback-limited fraction by redshift |
| `c17_diagnostics.py zseed\|mhot\|nrd\|trees` | `c17_zseed`, `c17_mhot`, `c17_nrd`, `c17_trees` | 100 paired (240 for trees), 801-step dt | z_seed, M_hot, n_rd, and the (invalid) PCH08 comparison |
| `c17_smooth_assembly.py` | `c17_smooth_assembly` | 1 deterministic history | frozen MVM on the Fakhouri et al. (2010) mean accretion history |
| `c17_accessibility_reach.py` | `c17_accessibility_reach` | 60, 801-step dt | fixed-seed reach of f_BH and critical-seed ratio across (sigma_j, R_nuc) cells |
| `c18_saturation_tests.py` | `c18_saturation_tests` | 100 paired, 801-step dt | seed memory, AGN strength (f_mom), all-gas-accessible control at sigma_j = 1.5, R_nuc = 300 pc |
| `c18_test4_phase.py` | `c18_test4_phase` | 100 paired, 801-step dt | stellar wind x AGN 2x2; phase-space d ln R/dt against R |
| `c18_eta_eps.py` | `c18_eta_eps` | 100 paired, 801-step dt | eta_acc x eps_sf grid, with and without feedback |
| `c17_pch08_check*.py` | `c17_pch08_check` | 40 to 60 | why PCH08 cannot be used here (near-deterministic, too-early main-progenitor history at small M_res) |

All scripts here and `gen_mvm_production.py` locate `foraois` through the environment variable `FORAOIS_ROOT` (a checkout containing `src/` and `config/`); see `docs/PRODUCTION_PROVENANCE.md` for the state used for the production ensemble and for what is not exactly reproducible.

`c17_*` diagnostics assert that `reservoir_stock.py` has the hash recorded in the production JSON, so they branch from the frozen commit.
`logs/c17_trees.log` reports the PCH08 comparison run **as it was run**; it is not a valid result (see the PCH08 evidence in the crib sheet, C17).

`c18_*` scripts assert the frozen module hash as well. The all-accessible control replaces `nuclear_transfer` at run time and restores it; the frozen file is not edited.
Data and figures: `output/c18_*.json`, `output/sat_fig1` to `sat_fig8` (`plot_c18_saturation.py`, `plot_c18_test4_phase.py`, `plot_c18_eta_eps.py`).

Caveats on individual logs:

* `premvm_vs_mvm_20trees.log`: the line "wind cap ... steps with gas" uses a wrong denominator (steps with gas at the
  start of the step) and can exceed 100 per cent. Use `feedback_limited_by_z_20trees.log`, whose denominator is
  the steps with any feedback demand, or the 801-step figures in `mass_scaling_240trees_*`.
* Diagnostics that scale with the time step (nuclear share of M_star, star formation before z = 10, G_BH - 1,
  the feedback-limited fraction) are resolution-dependent; see `timestep_paired_60trees_*`.
* `regime_*` and `converge_*` logs use fewer trees (20 or 60) than production; quote them for regimes and
  ratios, not as precision values of M_seed,crit.
