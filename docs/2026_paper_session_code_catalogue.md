# Code catalogue: "Differential Growth" paper session (Aug 2026)

This documents every piece of code written while developing
`AshviniPapers/High_z_Black_Hole_Growth/high_z_bhs.tex`. It was all built
in an isolated scratch sandbox, **not** on top of this repository — see
"Important: this is a parallel model, not an extension" below before
merging anything.

## 1. What was built, and why

### 1.1 Standalone prototype package (`ashvini_*.py`)

Three files reimplementing a three-reservoir (gas/star/BH) galaxy model
from scratch, independent of this repo's `ashvini/` package:

- **`ashvini_cosmology.py`** — halo mass accretion rate
  (Fakhouri, Ma & Boylan-Kolchin 2010 fit, as used by Davé, Finlator &
  Oppenheimer 2012), disc dynamical time, virial velocity dispersion.
- **`ashvini_halo_model.py`** — `HaloGalaxy` class: integrates gas, star,
  and BH reservoirs forward on a fixed-step RK4 scheme. Nuclear BH
  accretion via a free-fall estimator using total enclosed mass
  (Hobbs, Power, Nayakshin & King 2012) rather than Bondi-Hoyle, with a
  compaction-driven boost (Booth & Schaye 2009 functional form, anchored
  to Lapiner, Dekel & Dubois 2021). Super-Eddington growth allowed above a
  scanned threshold `r_crit`, motivated by slim-disc photon trapping
  (Watarai et al. 2000; Madau, Haardt & Dotti 2014; Lupi et al. 2024).
  AGN feedback via King (2003) energy-driven wind. Added this session: an
  `external_halo_track` parameter letting a precomputed `(t, z, M_halo)`
  trajectory — in particular a real PCH08 merger tree main-progenitor
  branch — replace the internal smooth-fit halo growth entirely.
- **`ashvini_bh_growth.py`** — thin wrapper/driver used by early
  single-galaxy runs (`run_ashvini.py`).

### 1.2 Halo-growth methodology fix

- **`adaptive_halo_seeding.py`** — `n_hubble_volume()` (Sheth-Tormen
  cumulative halo abundance via `colossus`, used to diagnose that a fixed
  `z_start=15` initial mass was unphysical for the paper's original
  Table 1 halo masses), `z_start_adaptive()`, and
  `shoot_seed_for_z5_target()` — bisects on an early ($z=25$) seed halo
  mass so forward integration lands within ~1% of a target mass at
  $z=5$, avoiding the finite-time blow-up in the super-linear
  ($M_{\rm halo}^{1.15}$) accretion rate.
- **`anchored_grid_run.py`**, **`table1_corrected.py`** — apply the
  shooting method across a halo x seed grid; regenerate a corrected
  Table 1 (the original was built on diverged trajectories for two of
  its three halo masses).

### 1.3 Core paper result: critical seed mass boundary

- **`seed_halo_boundary.py`** — `critical_seed(M_halo0_z25, target=0.5, ...)`:
  bisects on BH seed mass (strict-Eddington only) to find
  $M_{\rm seed,crit}(M_{\rm halo})$, the seed mass below which a black
  hole cannot reach $f_{\rm BH}=0.5$ by $z=5$ without a super-Eddington
  episode. This became the primary analysis function, reused throughout.

### 1.4 Sensitivity scans

- **`sensitivity_scan.py`** — accretion-parameter sensitivity
  ($\eta_{\rm acc}$, $\Phi_{\rm max}$, $R_{\rm nuc}$): supply/Eddington
  margin plus full $r_{\rm crit}$ bisection.
- **`phenomenological_dials.py`** — radiative efficiency $\epsilon$
  (mapped to BH spin via the Novikov-Thorne ISCO relation) and an
  angular-momentum transport penalty $f_{\rm am}$ (mathematically
  equivalent to $1/\eta_{\rm acc}$), scanned via `critical_seed`.
- **`sensitivity_band.py`** — optimistic/pessimistic boundary curves for
  the paper's sensitivity-band figure; extended mid-session with
  `_v2` variants anchored to physically-motivated
  Novikov-Thorne spin values ($a_*=0$ chaotic-accretion floor vs.
  $a_*=0.9$ coherent-accretion/inefficient-transport edge) rather than
  arbitrary placeholder values.
- **`scatter_check.py`** — multi-realisation scatter from the model's
  built-in stochastic (Ornstein-Uhlenbeck) accretion-rate perturbation.

### 1.5 Merger-tree comparison (uses `pymctrees` directly, not via this repo)

- **`gen_trees.py`**, **`gen_trees_grid.py`** — generate PCH08 merger tree
  ensembles (20 trees x 10 halo masses) via `pymctrees.build_forest_numba`,
  using a `menon_power_2024.yml` cosmology config
  ($H_0=67.8$, $\Omega_m=0.308$, $\Omega_\Lambda=0.692$, $\sigma_8=0.815$,
  $n_s=0.968$ — matching Menon & Power 2024), saved as `.npz`.
- **`mergertree_scatter.py`** — `critical_seed` evaluated per individual
  tree realisation (tree-to-tree scatter: 0.63-0.70 dex).
- **`boundary_from_trees.py`** — recomputes the full boundary curve using
  tree-**median** trajectories in place of the smooth fit. Central result:
  a systematic 1.8-2.5 dex upward shift in $M_{\rm seed,crit}$, larger
  than the tree-to-tree scatter, traced to real assembly histories being
  backloaded relative to a smooth mean-accretion-rate fit (a progenitor
  bias / conditional mass function effect).
- **`band_from_trees.py`** — re-anchors the $\epsilon$/$f_{\rm am}$
  sensitivity band onto the tree-based boundary.

### 1.6 Figures

`make_paper_figures.py`, `build_workings_notebook.py`, plus inline
snippets, producing `fig_sensitivity_band.png`, `fig_sensitivity_band_v2.png`,
`fig_smooth_vs_trees.png` (now copied into
`AshviniPapers/High_z_Black_Hole_Growth/figures/`) and the earlier
`fig1_trajectory.png`, `fig2_rcrit_dial.png`, `fig3_seed_rcrit*.png`.

All of the above (source + generated `.json`/`.npz` intermediates) live
in the ephemeral session sandbox at
`ashvini-papers/papers/2026_differential_growth_overmassive_bh/` and were
never pushed anywhere persistent except the two figures above and the
paper text itself. **They will not survive this session ending unless
copied somewhere persistent** (this repo, or `AshviniPapers`) — see
recommendation below.

## 2. Important: this is a parallel model, not an extension

This repo's actual `ashvini/black_holes_growth.py` and `gas_evolve.py`
implement a **different, already-mature BH growth model**:

| | This session's prototype | This repo (`black_holes_growth.py`) |
|---|---|---|
| Nuclear accretion | Free-fall on enclosed mass $M_{\rm BH}+M_{\rm gas}+M_\star$ (Hobbs et al. 2012), boosted by a compaction factor $\hat\Phi(t)^2$ | Free-fall on $t_{\rm ff}=0.141\times$Hubble time, `black_holes.efficiency` x gas mass |
| Super-Eddington cap | Graded, continuous: full rate below $r_{\rm crit}$, $\dot M_{\rm acc}/r_{\rm crit}$ above it | Hard cap at `eddington_multiplier` x Eddington rate (constant multiplier, no state-dependent grading) |
| AGN feedback | King (2003) energy-driven, constant $\epsilon_f$ | Optional PZNK11 (Power, Zubovas, Nayakshin & King 2011) M-sigma self-regulation with a hard `growth_ceiling`, or a constant `eta_agn` wind coupling |
| Halo growth | Smooth Fakhouri+10 fit, $z=5$-anchored via a seed-mass shooting method, OR an externally supplied trajectory | Real PCH08 merger trees via `pymctrees_adapter.py`'s `build_forest_live` / offline HDF5 — **already solves exactly the problem the shooting method was invented for** |
| Seeding | Not modelled as a threshold — seed mass is a free parameter scanned over a grid | Three configurable channels (Pop III / direct collapse / halo-mass threshold), `run_params.yaml: black_holes.seeding` |

Two consequences worth flagging before any merge:

1. **The halo-growth methodology fix (Section 1.2 above) is largely
   redundant with what this repo already does correctly.** This repo
   ingests real merger trees from `pymctrees` and never integrates a
   smooth, unconditioned mean accretion-rate fit forward from an
   arbitrary $z_{\rm start}$ — the exact failure mode (finite-time
   blow-up under $M_{\rm halo}^{1.15}$ growth, extrapolated into an
   extreme-value-improbable initial mass) that motivated the shooting
   method doesn't arise here. The paper's Section 5.6 result (merger
   trees shift $M_{\rm seed,crit}$ by 1.8-2.5 dex relative to a smooth
   fit) is itself evidence *for* using this repo's tree-based approach
   over a smooth fit, not a reason to import the smooth-fit machinery.

2. **The nuclear-accretion/feedback physics is a genuine fork, not a
   bug fix.** Both models are physically motivated but make different
   choices (free-fall on total enclosed mass with a graded slim-disk cap
   and King-2003 feedback, vs. Hubble-time free-fall with a hard
   Eddington multiplier and PZNK11 M-sigma self-regulation). Merging
   the paper's prescription in means adding it as a genuinely new,
   independently selectable option (e.g. a `black_holes.growth_model:
   "hobbs_slimdisk"` vs. `"pznk11_freefall"` switch), not overwriting or
   quietly changing the existing default that this repo's own tests and
   the Menon & Power (2024) reproduction notebook already depend on.

## 3. What was actually merged (done)

Per your call ("full merge as new growth_model option"), the following
landed in this repo, none of it changing default (`growth_model:
"pznk11_freefall"`) behaviour — verified by the full existing test suite
still passing (103/103) plus new tests specifically checking the default
path is untouched:

- **`ashvini/spin.py`** — `epsilon_from_spin(a_star)` (Novikov-Thorne
  ISCO relation), `isco_radius(a_star)`. Reproduces the paper's four
  reference epsilon values exactly. Note: the prototype's `isco_radius`
  formula (ported verbatim at first) was wrong for retrograde spin
  (a*<0) — fixed here (Z1/Z2 are even in a, so the prograde/retrograde
  branch must be chosen explicitly, not by feeding a signed a* through a
  single fixed-sign formula). Doesn't affect the paper, which only uses
  a*>=0, but worth knowing if this gets reused elsewhere.
- **`ashvini/black_holes_growth_slimdisk.py`** — the full Hobbs/
  compaction/graded-r_crit/King(2003) model as a selectable alternative:
  `time_freefall_enclosed`, `nuclear_accretion_rate`,
  `eddington_rate_std_per_unit_mass`, `bh_growth_step_slimdisk` (closed-
  form three-regime solver — constant/exponential/constant — analytically
  exact, not an approximation; verified against a numerical `solve_ivp`
  reference and against the existing model's own `_bh_growth_step` in
  the `r_crit -> infinity` limit), `king_agn_wind_rate`.
- **`ashvini/seed_mass_function.py`** — the boosted-fraction-by-channel
  utility.
- **`run_params.py`/`run_params.yaml`** — `black_holes.growth_model`
  (`"pznk11_freefall"` default / `"hobbs_slimdisk"`) and a
  `black_holes.slimdisk` config block.
- **`main.py`** — `run_forest()` branches on `growth_model`; the
  `pznk11_freefall` branch is byte-identical to the pre-existing code
  (just indented under `else:`). `run1_scalar()` was **not** extended —
  it remains the `pznk11_freefall`-only reference implementation.
- **Tests**: `tests/test_black_holes_growth_slimdisk.py` (18 tests:
  spin helper, accretion-rate scaling, the closed-form solver against
  both a numerical reference and the existing model's limit, King wind
  scaling, end-to-end `run_forest` wiring for both models) and
  `tests/test_seed_mass_function.py`.
- **`MODELS.md`/`README.md`** — new sections/table rows documenting the
  above.
- **`AshviniPapers/High_z_Black_Hole_Growth/analysis/`** — every session
  script, generated `.json`/`.npz`/`.npy` intermediate, and the
  standalone `prototype_package/` (the original `ashvini_*.py` files),
  archived for reproducibility with its own README pointing back here.

Not ported (deliberately, per Section 2 above): the smooth-fit/shooting
halo-growth methodology, and the stochastic time-varying compaction
boost (now a fixed per-run multiplier, `black_holes.slimdisk.compaction_boost`).

Nothing has been committed — working tree changes only; `git status`
shows the full diff when you're ready to review and commit.

## 3a. Original recommendation (superseded by Section 3 above)

- Copy the session's scripts and intermediate results into
  `AshviniPapers/High_z_Black_Hole_Growth/analysis/` (paper-specific,
  reproducibility archive) regardless of what happens next — they are
  paper support material, not a package, and belong with the paper.
- Only add code to *this* repo (`ashvini/`) for pieces that are genuinely
  reusable beyond this one paper and that don't silently change existing
  default behaviour, e.g.:
  - A new, opt-in nuclear accretion + super-Eddington cap module
    (`ashvini/black_holes_growth_slimdisk.py` or a `growth_model` switch
    inside the existing file), implementing the Hobbs/compaction/graded-
    $r_{\rm crit}$ physics as an alternative to the current PZNK11/
    Hubble-time-free-fall default.
  - The Novikov-Thorne $\epsilon(a_*)$ helper (small, genuinely general,
    no dependency conflicts) — natural fit in `black_holes_growth.py` or
    a new `spin.py`.
  - The seed-mass-function-to-boosted-fraction conversion
    (`seed_halo_boundary.py`'s later half) — could live as a small
    post-processing utility, not core integration code.
- Do **not** port the smooth-fit shooting method
  (`adaptive_halo_seeding.py`) into this repo as a tree-generation
  alternative; this repo's `pymctrees_adapter.py` already does the
  physically better version of that job.
