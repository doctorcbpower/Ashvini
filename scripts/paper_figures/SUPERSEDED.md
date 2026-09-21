# Superseded and exploratory material in `scripts/paper_figures/`

Written 2026-09-21. The current paper uses the Minimal Viable Model (`ashvini/reservoir_stock.py`) and the frozen production ensemble (`output/mvm_production_results.json`); see `docs/PRODUCTION_PROVENANCE.md`. Everything listed here belongs to the earlier "Differential Growth" model (`ashvini/paper_reservoir.py`, `run_reservoir_paper`, `critical_seed_paper`) or was exploratory, and is **not used for any result, number or figure of the current paper**. Nothing has been moved or deleted, so existing paths and git history are unchanged; this file is the marker. Whether a file uses the old model was checked from its imports.

The old-model generators use the older tree ensembles (1000 trees at 25 masses), the compaction factor, the chi_crit cap and the energy-driven wind, none of which exist in the MVM. Their numbers must not be quoted.

## Old-model generators and plots (tracked)

| Script | Output |
|---|---|
| `gen_fig1_critical_boundary.py`, `plot_fig1_critical_boundary.py` | `output/fig1_critical_boundary_results.json`, `output/fig_critical_boundary.png` |
| `gen_fig2_sensitivity.py`, `plot_fig2_sensitivity.py` | `output/fig2_sensitivity_results.json`, `output/fig_sensitivity_panels.png` |
| `gen_fig3_trajectories.py`, `plot_fig3_trajectories.py` | `output/fig3_trajectories_results.json`, `output/fig_representative_trajectories.png` |
| `gen_fig4_chicrit.py`, `plot_fig4_chicrit.py` | `output/fig4_chicrit_results.json` |
| `gen_fig5_fbhsensitivity.py`, `plot_fig5_fbhsensitivity.py` | `output/fig5_fbhsensitivity_results.json` |
| `gen_zseed_check.py` | `output/zseed_results.json` |
| `gen_compaction_sigmaj_check.py` | `output/compaction_sigmaj_results.json` |

These scripts reproduce the figures of the pre-MVM manuscript draft, which survive in the paper repository's git history. They still contain the machine-specific foraois path and need it edited by hand to run.

## Old-model controls and exploratory scripts (untracked)

| File | What it is |
|---|---|
| `control_fig_trajectories.py`, `paired_fig_trajectories_winds.py`, `output/control_pk0_winds_off.json`, `output/paired_pk0.json`, `output/fig_representative_trajectories_control_pk0_winds_off.png`, `output/fig_representative_trajectories_paired_pk0.png` | control runs of the old trajectory figure with the tree power spectrum evaluated at a different redshift, made while the foraois barrier normalisation was being investigated |
| `output/old_pk_z5/` | earlier copies of old-model outputs made before those controls |
| `gen_correa_comparison.py`, `output/correa_comparison_results.json` | comparison of tree main-progenitor histories with the Correa et al. (2015) smooth fit, intended to replace the Fakhouri et al. (2010) comparison in the old draft. It uses neither reservoir model. The current paper uses the Fakhouri comparison (`diagnostics/c17_smooth_assembly.py`) |

They are not tracked by git; they are listed here so that they are not mistaken for part of the frozen calculation.

## Tracked old-model outputs

`output/compaction_sigmaj_results.json`, `fig1_critical_boundary_results.json`, `fig3_trajectories_results.json`, `fig4_chicrit_results.json` and `zseed_results.json` are the committed versions. On 2026-09-20 (09:15 to 09:18, before the foraois barrier normalisation change at 11:17) the old-model generators rewrote them with different results; the tree normalisation used for that rewrite was not recorded. The rewritten files were not committed. They were restored to the committed versions, and copies of the rewritten files and their run logs were kept outside the repository. Nothing in the current paper uses any of them.

## Not superseded

`gen_mvm_production.py`, `mvm_crib_numbers.py`, `plot_mvm_fig*.py`, `plot_c18_*.py`, `paper_style.py`, `foraois_paths.py`, `mnras_science.mplstyle`, `diagnostics/` and its logs.
