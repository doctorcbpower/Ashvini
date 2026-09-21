# Provenance of the frozen production calculation

Status: record written on 2026-09-21 from the repository chronology (commit dates, file modification times, diffs and hashes). No production calculation was rerun to write it (one small scratch check, made for the audit and not stored, is mentioned in section 2 and labelled as such). Where a statement rests on file modification times and not on a logged value, it is marked *(inferred)*. The foraois commit hashes cited here are those of the foraois history after its message-only rewrite (section 2); the earlier hashes are listed there.

The calculation described here is the production ensemble behind the results of the manuscript "Critical seed masses for massive black holes in the early Universe". Model and production results are frozen.

## 1. What is frozen

| Item | Identifier |
|---|---|
| Production output | `scripts/paper_figures/output/mvm_production_results.json` |
| Its SHA-256 | `5e8ec4c44175da253f73a2f47a80226b0207005e01cdd7f96047136bca4eff05` (git blob `737e9b8c7e62cba899280b5e9d43b8a1faa6cda0`) |
| Written | 2026-09-20 15:44:29 (`meta.date` and file time agree). This is when the JSON was written, that is the end of the run. The start time and duration were not recorded |
| Model source | `ashvini/reservoir_stock.py`, SHA-256 `31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63` (recorded in the JSON `meta.module_sha256`; equal to the committed blob and to the working file) |
| Ashvini commit | `695b11467aaa68e4c7eb465e13a53e54a3135af3` (2026-09-20 15:56:21 +0800), "Add the Minimal Viable Model of the galaxy and nuclear reservoirs, its production ensemble, diagnostics and crib sheet"; on `fork/main` (`doctorcbpower/Ashvini`). This is the first commit that records the model, the production script and the output. It was made about 12 minutes after the output was written, so the run was not made from a clean checkout of it |
| Production script | `scripts/paper_figures/gen_mvm_production.py` as committed at `695b114` (SHA-256 `51b0a5aa97024f7ff3a7c78e5f8a514a81a2447210ca25f8f58efbdc69349c68`). The bytes that ran are not established (see below). A later edit replaced only its hard-coded foraois paths (section 5) |

Configuration recorded by the script and the JSON: 13 halo masses, 3e10 to 3e13 Msun at 0.25 dex; 240 Zhang and Hui trees per mass; 801 uniform-cosmic-time nodes (800 steps of about 1.3 Myr); dz = 0.05; M_res = 1e4 Msun; z = 25 to z = 5; seed inserted at the first node; f_BH = 0.5; fixed seeds 1e2, 1e3, 2e5 and 1e7 Msun in the fiducial model and in an R_nuc = 250 pc delivery experiment. The fiducial parameters are stored in `meta.fiducial` (R_nuc = 100 pc, sigma_j = 0.5, eta_acc = 0.005, eps_sf = 0.015, f_mom = 1, eta_sn_scale = 1, n_rd = 10, c_NFW = 4, epsilon = 0.1, f_b = 0.156, M_hot = 4e11, lambda = 0.035). They equal the current defaults of `ashvini/paper_reservoir_params.py`.

## 2. Which code produced it

**Ashvini.** Every module that `reservoir_stock.py` and the production script import is byte-identical between the commit `695b114` and the current HEAD, with one exception, `ashvini/pymctrees_adapter.py`. This compares commits; it is not knowledge of the state of the files when the run was made. What is known about that state:
* `ashvini/reservoir_stock.py`: established. Its SHA-256 is recorded in the JSON.
* The other imported modules (`paper_reservoir.py`, `utils.py`, `run_params.py`, `constants.py`, `reionization.py`, `spin.py`, `supernovae_feedback.py`, `black_holes_growth*.py`, `__init__.py`) have file times before the JSON was written (`paper_reservoir.py` at 2026-09-20 11:36:23, the rest earlier) *(inferred)*.
* `paper_reservoir_params.yaml` and `ashvini/paper_reservoir_params.py` were saved at 15:55:20, after the JSON was written (15:44:29) and just before the commit, so their state at run time is not known from file times. The committed change only adds fields (`stellar_winds: false`, `eta_sn_scale`, and options of the frozen pre-MVM reference); no existing value changed, and the recorded `meta.fiducial` equals the committed defaults.
* `ashvini/pymctrees_adapter.py`: no committed version changed it between `5f5f3c9` (2026-09-12) and `0df9dae` (2026-09-20 22:05), so the version at `695b114` is the `5f5f3c9` version. It was then changed in `0df9dae` and in `fc93fac` (2026-09-21 10:37), both after the run, and only in `build_forest_live` and its warning helper. The production script calls `build_forest_for_bin`, whose source is identical in every committed version from `5f5f3c9` to the current HEAD. The file's current modification time (2026-09-21) says nothing about its state at run time, and any uncommitted edits at that time cannot be recovered.

`scripts/paper_figures/gen_mvm_production.py` has a modification time 79 s after the JSON. The committed version is the only surviving record of the script; whether the last saved bytes differ from those that ran cannot be established. Its parameters are recorded in the JSON `meta`.

**foraois** (`doctorcbpower/foraois`, trees, sampler and cosmology). The production run used the foraois working tree on top of base commit `4370a4a` (`4370a4a4c59533f9a8f1f19ed8b371df82beb5a8`, 2026-09-15 11:08:30) with the collapse-barrier normalisation change already applied but not yet committed *(inferred from the file times below)*:

| Evidence | Value |
|---|---|
| `src/foraois/cosmo_utils.py` and `tests/test_cosmo_utils.py` last modified | 2026-09-20 11:17:03 and 11:16:50, before the JSON was written (15:44:29, the end of the run). The run start is not recorded, so these file times alone place the change before the end of the run, not before its start |
| Commit that recorded that change | `ec1a66f6a8b0e9b2b99315c365645407eb73b3bc`, "Normalise the collapse barrier to the redshift at which P(k) is evaluated", 2026-09-20 16:23:36, parent `4370a4a`; touches only `cosmo_utils.py` and its test |
| `cosmo_utils.py` changed after `ec1a66f`? | No |
| Other modules the tree generator executes (`zhang_hui_trees.py`, `pch_trees.py`, `transfer_functions.py`, `mass_function_utils.py`, `utils/io.py`, `utils/window_function.py`, `first_crossing*.py`) | last modified before 2026-09-20 except the docstring-only edits below |
| Config `config/menon_power_2024.yml` | unchanged since the initial release (`d9d43bb`); SHA-256 `3b3736a6b4b6d00a4b9760fb8a1b3ed5085028862fe900040bbb32bd5e671d69` |

**State of the tree code at the run.** No commit existed for it. The run was made from `4370a4a` plus an uncommitted change to `cosmo_utils.py` (and its test), which was committed afterwards as `ec1a66f`, about 40 minutes after the JSON was written. `ec1a66f` is therefore not the commit at run time. It is the recommended reproducibility pin: its content is equivalent to the run-time state for the modules the production script imports *(inferred, from the file times and diffs above)*.

**Rewrite of the foraois history (message-only).** After the first version of this record, the private foraois `main` history was rewritten (2026-09-21) solely to remove three Claude-related commit-message suffixes: a co-author trailer on one commit and an acknowledgement sentence on two others. The rewrite was verified to preserve all 28 corresponding file trees exactly, all author and committer identities, all timestamps and the parent structure; no file content changed. Only the hashes of the commits from the first rewritten commit onward changed. The hashes cited in this record are those of the rewritten history. The earlier hashes, which the Ashvini commits made before this note cite, correspond as follows:

| Earlier hash | Hash in the rewritten history | Role |
|---|---|---|
| `1ba70733cce4e5a3aa44da15a2f82f9d7cbe5c2e` | `4370a4a4c59533f9a8f1f19ed8b371df82beb5a8` | base of the run-time state |
| `ed28326488d88ebecc274c7accb75edcbd8ffe00` | `ec1a66f6a8b0e9b2b99315c365645407eb73b3bc` | barrier-fix commit; **the reproducibility pin** |
| `9d5fdae4a7afb0b59b2e74345b6ed83f8a208125` | `767398d7e1833b19b50ac8bef439d21fbd9d5f63` | first later commit |
| `74f15ad57abccf739ea9421958dbfbac2943ac0c` | `c80858d81fa2902fdbc99b8354790bd030ad9a35` | later commit |
| `7d6fa29f6f9e5d8ced41435df61218bda77cd2e1` | `d069182a728eeea41cf5d056606672c806648bff` | later commit |
| `124d4888956cd385f88f3fb860125b84d22be8df` | `8f9508731fac10dfe18e6997e44ffbc0a314372e` | later commit |
| `02157c94bad6ea95fbee6a76f53522f09d371e61` | `734c2f16cdeb8ef83f199c24be5241cfa9a8020e` | tip of `main`; the tag `v0.1.1` still pointed to the earlier hash at the time of writing |

The reproducibility pin is now `ec1a66f6a8b0e9b2b99315c365645407eb73b3bc` (formerly `ed28326`). It remains the recommended pin for the inferred run-time code state, not a claim that it was the literal production commit: no commit existed at run time, and the state was `4370a4a` plus an uncommitted change later committed as `ec1a66f`. Checks described in this record that were made before the rewrite (the AST comparison, the import closure, the scratch check and the compliance evaluation) were made on the earlier hashes; the trees are identical, so they apply unchanged to the rewritten ones. The limitations recorded here (script and adapter bytes not established, run start and environment not recorded, the stochastic Zhang and Hui sampler not seed-reproducible, and the JSON storing representative histories rather than the full ensemble) are unchanged. At the time of writing the rewritten history had not been pushed.

Unstored supporting evidence, not provenance. In a scratch computation made for the audit on 2026-09-21 (40 trees at M0 = 3e10 Msun, same script settings, not the production run, results not saved), foraois at `ec1a66f` gave a median critical seed of 9.7e7 Msun (16 to 84 per cent range 7.8e7 to 1.07e8), foraois at `4370a4a` gave 5.4e7 Msun (4.2e7 to 6.4e7), and the production JSON has 9.45e7 Msun (7.6e7 to 1.06e8) at that mass. This is consistent with the production trees having been made with the barrier change, and is one mass and a small sample.

The production script builds `CosmoData(rp, redshift=[5.0])` and `ZhangHuiMergerTree(cd, rp, model="cdm")` and calls `build_forest_for_bin(..., z0=5.0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend="numba", rng_seed=1000+i)`.

**foraois changes after the run and what they affect.**

| Commits (on the rewritten foraois `main`, now `734c2f1`) | Content | Affects the production trees? |
|---|---|---|
| `767398d`, `c80858d`, `8f95087`, `734c2f1` | documentation, reproduction scripts, diagnostics | No |
| `d069182` | figure styling (`utils/plot.py`, new `utils/paper_style.py`) | No; not imported by the tree generator (`plot.py` imports `diagnostics` lazily inside a plotting function) |
| `8f95087`, `734c2f1` | wording in the docstrings of `__init__.py`, `tree_algorithm.py` and `first_crossing.py` (changes "exact" to "binary-per-step ... first-crossing"); one changed argparse help string in `main.py` (command-line help text, a code line but not a computation); new `expected_eps_splits_per_step` and reworded docstrings in `diagnostics.py` | No. The three docstring-only files are identical to `ec1a66f` after docstrings are removed (AST comparison). `main.py` and `diagnostics.py` are not among the 15 foraois modules the production imports load, which are the same at `ec1a66f` and at the current HEAD |
| `pyproject.toml`: `d069182` (new optional `paper` extra, SciencePlots) and `8f95087` (version 0.1.0 to 0.1.1) | package metadata | No; no dependency of the package itself changed |

## 3. What is and is not exactly reproducible

* **Not exactly reproducible: the tree ensemble.** The numba Zhang and Hui sampler does not reproduce a given ensemble on rerunning. The script passes `rng_seed = 1000 + i`, which seeds only NumPy's global generator (`np.random.seed` in `build_forest_for_bin`); repeated draws of nominally identical configurations were found not to match (`scripts/paper_figures/diagnostics/README.md`). A rerun therefore gives a new ensemble. The manuscript states an observed ensemble-to-ensemble variation of about 5 per cent in the fiducial median critical seed at 3e10 Msun.
* **The JSON is the frozen record.** It stores per-tree critical seeds and end-state quantities, and percentile bands and one representative halo history per mass (`rep_halo`). It does not store every tree's halo history, so the reservoir calculation cannot be rerun on the production trees.
* **Reproducible from the stored production JSON:** the three production figures (boundary, inverse, host baryons) and the numbers in the crib sheet. `mvm_crib_numbers.py`, `plot_mvm_fig1_boundary.py`, `plot_mvm_fig2_inverse.py` and `plot_mvm_fig3_hostbaryons.py` read only that JSON. The other four manuscript figures from this work (seed trajectories, final versus seed, phase space, eta-epsilon grid) read separate diagnostic ensembles, stored in the `c18_*_results.json` files (`plot_c18_saturation.py`, `plot_c18_eta_eps.py`, `plot_c18_test4_phase.py`), not the production ensemble.
* **Deterministic given the halo histories:** `run_reservoir_stock` and `critical_seed_stock` in the frozen module. Their tests are in `tests/test_reservoir_stock.py`.
* **Diagnostics** (`scripts/paper_figures/diagnostics/`) each use their own tree ensemble; their logs are the record. `c17_diagnostics.py`, `c17_accessibility_reach.py`, `c17_smooth_assembly.py`, the three `c18_*` scripts (`c18_saturation_tests.py`, `c18_eta_eps.py`, `c18_test4_phase.py`) and `c19_timestep_compliance.py` assert the module hash (`assert_frozen`) before running the reservoir. The `c17_pch08_check*.py` scripts and the `mvm_*.py` diagnostics do not.

## 4. Known limitation of the production trees

The production timestep is outside the practical single-split compliance regime of the Zhang and Hui builder (foraois `docs/PCH08_HIGH_Z_DIAGNOSTIC.md`: expected splits per step of about 0.1 or less). Evaluated with foraois `diagnostics.expected_eps_splits_per_step(cosmo_data, M0*h, z, z+0.05, 1e4*h, model="cdm")`, using the cosmology of `menon_power_2024.yml` and `CosmoData(rp, redshift=[5.0])`, with M0 the anchored halo mass (conservative at high z), the maximum over z from 5 to 25 was about 82, 690 and 6.0e4 at M0 = 3e10, 3e11 and 3e13 Msun. These values were evaluated once, on 2026-09-21 and before the foraois history rewrite, with the commit then named `02157c9`; they were not re-evaluated after the rewrite. That commit is now `734c2f1`, whose file tree is identical to the pre-rewrite `02157c9` tree, so the evaluation applies to the same code state. It says nothing about the code used for the production run, which was not made at either commit. The values are indicative: they are not stored by a committed script or log. The manuscript states the limitation. The effect on the median critical seed was measured afterwards in diagnostic C19 (`scripts/paper_figures/diagnostics/c19_timestep_compliance.py`, log `logs/c19_timestep_compliance.log`, data `output/c19_timestep_compliance_results.json`): trees built at a smaller step and recorded at the production checkpoints give a median about 6% lower at 3e10 (compliant, E about 0.08), about 7% lower at 3e11 (E about 0.28) and about 4% lower at 3e13 (E about 60, not compliant), so the production values are probably slightly high. That is a diagnostic on separate tree ensembles; the frozen production JSON is unchanged and the other limitations of this record stand.

## 5. Machine-specific paths

`gen_mvm_production.py`, `diagnostics/mvm_*.py`, `diagnostics/c17_*.py` (and through them `c18_*.py`) previously hard-coded `/Users/00075868/MyCodes/foraois/{src,config}`. They now import `scripts/paper_figures/foraois_paths.py`, which reads the environment variable `FORAOIS_ROOT` and raises a clear error if it is unset. This changes path resolution only. The reference version of the production script is `695b114:scripts/paper_figures/gen_mvm_production.py`; the bytes that actually ran are not established (section 2).

Not changed, because they are superseded or outside the frozen calculation: the old-model generators listed in `scripts/paper_figures/SUPERSEDED.md`, and the SHMR and resolution scripts in `scripts/` (they read `planck2018_camb.yml` under the same machine-specific directory). A rerun of any of these needs the path edited by hand.

## 6. Environment

Not logged at run time. The virtual environment `.venv` present when this record was written: Python 3.11.10, numpy 2.4.6, scipy 1.17.1, astropy 8.0.1, numba 0.66.0. The test suite passes in it.

## 7. Repository locations

Ashvini: `https://github.com/doctorcbpower/Ashvini` (branch `main`; remote `fork`); its upstream is `https://github.com/Anand-JM97/Ashvini`. foraois: `https://github.com/doctorcbpower/foraois`. Manuscript: `https://github.com/doctorcbpower/High-Redshift-Black-Hole-Growth`.

## 8. Figure files in the manuscript repository

The manuscript repository's `.gitignore` excludes `*.pdf` and `*.md` (built PDFs and notes) and tracks the manuscript figures as PNG files in `figures/`. This is the existing convention, so the manuscript includes the PNG versions, at 300 dpi, and a clean checkout builds from tracked files alone (verified by building with all figure PDFs removed). Vector PDF versions of each figure are written next to the PNGs by `scripts/paper_figures/paper_style.py` (and by `figures/make_schematic_baryon_flow.py` in the manuscript repository); they are not tracked. If vector figures are wanted for submission, add `!figures/*.pdf` to that `.gitignore` and switch the `\includegraphics` extensions; nothing else depends on it. Regenerating a figure needs only the stored JSON files (the production JSON or the `c18_*_results.json` files, section 3); no simulation is rerun.
