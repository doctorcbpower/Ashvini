# Code catalogue: "Differential Growth" paper session (Aug 2026)

> **Historical record (2026-09-21).** This catalogue documents the earlier "Differential Growth" model (`paper_reservoir.py`) and the scripts that built the pre-MVM manuscript draft. That model and draft have been superseded by the Minimal Viable Model (`reservoir_stock.py`) and the current paper; every number, figure and parameter below refers to the old model and must not be quoted for the current paper. See `PRODUCTION_PROVENANCE.md` and `scripts/paper_figures/SUPERSEDED.md`.

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

## 4. 2026-09 session: Section 2 rewrite, critical_seed.py, paper_reservoir.py

Follow-up session, working directly against the paper draft
(`AshviniPapers/High-Redshift-Black-Hole-Growth/high-z-black-hole-growth.tex`)
rather than a standalone sandbox. Unlike Section 1-3's sandbox work, this
session's code was written directly into this repo throughout, so nothing
here is at risk of being lost the way the original prototype was.

### 4.1 What was built

- **`ashvini/main.py`**: `run_forest()` gained explicit per-call overrides
  for the "hobbs_slimdisk" branch's physics parameters (`growth_model`,
  `a_star`, `epsilon`, `r_crit`, `epsilon_f`, `eta_acc`, `gas_mass0`,
  `stars_mass0`), all defaulting to `None`/`0.0` (= prior behaviour
  unchanged). Previously these were only settable via `run_params.yaml`
  at import time, which made any kind of parameter scan require
  reloading the whole module per point. Covered by new tests in
  `tests/test_black_holes_growth_slimdisk.py`.
- **`ashvini/critical_seed.py`** (new): general-purpose vectorised
  bisection for $M_{\rm seed,crit}$ against `main.run_forest()`, one
  independent root-find per halo/tree via the existing
  `black_holes.seeding.halo_mass_threshold` channel (an `(N,)`-shaped
  `M_seed` array broadcasts against it correctly with no core-package
  change needed -- see the module docstring).
- **`ashvini/paper_reservoir.py`** (new, substantial): a *separate*
  reservoir integrator implementing the paper's own bespoke gas-supply
  physics, not `main.run_forest()`'s general-purpose model -- see its own
  module docstring for the full rationale. Implements:
  - Tree-based (or smooth Fakhouri et al. 2010) halo growth.
  - `zeta = zeta_UV * zeta_ch`: UV suppression (existing
    `reionization.uv_suppression`) times a new cold/hot-mode transition
    (`cold_hot_mode_suppression`, reusing `reionization.s` with different
    arguments, per the paper's own Eq. coldhot).
  - `low_spin_available_fraction()`: the mechanism that replaced an
    unstable/broken "feed the whole galaxy gas reservoir into a 100 pc
    free-fall time" picture (see 4.2 below) -- gas parcels retain their
    specific angular momentum, and only the low-spin tail (log-normal in
    $j$, median tied to the halo's own $\lambda$ via Bullock et al. 2001)
    can circularize within $R_{\rm nuc}$.
  - `critical_seed_paper()`: the paper-specific analogue of
    `critical_seed.critical_seed()`, against `run_reservoir_paper()`
    instead of `main.run_forest()`.

### 4.2 What was found (in rough chronological order)

1. **The paper's own Eq. halo_accretion (smooth specific baryon accretion
   rate) barely grows a halo at all between $z=25$ and $z=5$** -- checked
   directly, $M_{\rm halo}(z{=}25)/M_{\rm halo}(z{=}5)\sim0.94$ for a
   fiducial $5\times10^{10}\,M_\odot$ halo, two orders of magnitude too
   flat for real hierarchical assembly. Almost certainly a transcription
   error from an earlier drafting pass; replaced with the standard
   Fakhouri, Ma & Boylan-Kolchin (2010) fit for the smooth baseline (real
   trees are the primary halo-growth source regardless).
2. **Feeding the entire galaxy-scale $M_{\rm gas}$ directly into the
   $R_{\rm nuc}=100$ pc free-fall time produces an unphysical runaway**:
   $t_{\rm ff}\propto M_{\rm enc}^{-1/2}$ falls as the reservoir grows,
   accelerating consumption of gas that isn't actually all sitting at
   $R_{\rm nuc}$. This is what the low-spin mechanism (4.1) fixes -- but
   getting there took two wrong turns, both worth recording since they'd
   be easy to reintroduce:
   - Using the *nuclear-only* subset ($M_{\rm gas,nuc}$) to set its own
     angular-momentum threshold is an unstable fixed point -- a small
     downward fluctuation tightens the threshold, which shrinks
     $f_{\rm avail}$ further, collapsing to zero within a few steps.
     Fixed by evaluating the threshold against the *full* galaxy
     reservoir ($M_{\rm BH}+M_{\rm gas}+M_\star$) instead, which grows
     monotonically.
   - `run_reservoir_paper()` was calling `reionization.uv_suppression()`
     with a *re-derived smooth* halo growth rate
     (`smooth_halo_growth_rate(M_halo, z)`) even when `M_halo` itself came
     from a real, bursty tree. The mismatch between a tree's actual
     instantaneous growth rate and what a smooth analytic fit predicts
     for a halo of that mass/redshift drove `uv_suppression`'s internal
     formula negative and clipped it to a hard `0.0000` for stretches of
     real trajectories (not a smooth partial suppression). Fixed by
     computing the growth rate via direct finite-differencing of
     whatever `halo_mass` trajectory was actually passed in, rather than
     re-deriving it from a formula that may not match.
3. **A genuine theoretical result, not just a calibration note**: once
   black-hole growth is gas-supply-limited (rather than
   Eddington-limited), $\dot M_{\rm BH}$ and $\dot M_\star$ draw on the
   identical reservoir/clock, so $f_{\rm BH}$ asymptotes to the constant
   $\eta_{\rm acc}/\epsilon_{\rm sf}$ regardless of seed mass. This means
   $\eta_{\rm acc}/\epsilon_{\rm sf}<f_{\rm BH}$ is a *necessary condition*
   for $M_{\rm seed,crit}$ to be a meaningful, seed-dependent quantity at
   all -- above that ratio, every seed reaches the target through supply
   alone. The paper's own originally-stated fiducial values
   ($\eta_{\rm acc}=0.01$, $\epsilon_{\rm sf}=0.02$) sit *exactly* on this
   critical ratio (0.5), and $M_{\rm seed,crit}$ is a steep, near-singular
   function of the ratio as it approaches 0.5 from below -- checked
   directly (a 1% change in $\eta_{\rm acc}$ near the boundary changes
   $M_{\rm seed,crit}$ by a factor of a few). Deliberately chose *not* to
   tune into this knife-edge regime to match old reported numbers exactly
   (a referee-vulnerable choice); adopted $\eta_{\rm acc}=0.005$,
   $\epsilon_{\rm sf}=0.02$ (ratio 0.25, comfortable margin) as the
   working fiducial instead.
4. **With that fiducial choice, real `foraois` trees (not the smooth
   trajectory) at the paper's own fiducial halo mass reproduce the
   originally-reported $M_{\rm seed,crit}\approx2900\,M_\odot$ to within
   10%** ($3204\,M_\odot$, $N=1000$ trees) -- without tuning toward that
   number at all. Strong, if incomplete (see 4.3), validation that the
   reconstructed model is quantitatively consistent with whatever the
   lost original prototype actually did.
5. **Real tree-to-tree scatter in $M_{\rm seed,crit}$ is $\sim$0.05-0.07
   dex, not the paper's currently-drafted 0.6-1.0 dex claim** -- checked
   directly at $N=1000$ trees, three halo masses, and confirmed via the
   two most divergent individual trees in one ensemble (0.65 dex apart in
   mass 50 Myr after formation, 2 full redshift units apart in formation
   time) still converging to within ~15% of each other by $z=5$.
   Mechanism: every tree in a mass bin is anchored to the *identical*
   $M_{\rm halo}(z=5)$, and because accretion scales super-linearly with
   halo mass and rises toward low $z$, the late-time phase -- when every
   tree's mass is converging toward that shared endpoint -- dominates the
   *integrated* baryon budget the reservoir model cares about, far more
   than the genuinely diverse early-time phase. Real diversity exists
   (confirmed directly) but doesn't survive to $z=5$ in this model. The
   old 0.6-1.0 dex claim was very likely partly or wholly an artifact of
   the `uv_suppression` clipping bug (point 2 above) rather than genuine
   physical scatter -- this needs real revision in the paper's Section
   4.5, not just updated numbers, if it holds at the full grid.

### 4.3 What's still open

- **Only checked at 3 of the eventual 25 halo masses**, $N=1000$ trees
  each, single tree-generation seed. The tight-scatter finding (4.2,
  point 5) and the boundary shape need confirming across the full grid
  before being treated as final.
- **$M_{\rm res}$ convergence was never actually run for this specific
  calculation** (Zhang-Hui trees, $z_{\rm max}=25$, converging
  $M_{\rm seed,crit}$) -- the existing `docs/RESOLUTION_CONVERGENCE.md`
  study was for a different calculation (PCH08 trees, $z_{\rm max}=15/30$,
  converging $z=0$ $M_\star$) and doesn't automatically transfer. Currently
  using a fixed absolute $M_{\rm res}=10^4\,M_\odot$, chosen by physical
  argument (needs to resolve structure at $z\sim25$ for the halo-mass
  range considered) rather than verified by a joint $dz\times M_{\rm res}$
  grid test the way the SHMR calculation was.
- **$\sigma_{\ln j}$ (angular-momentum log-normal width) and
  $\lambda_{\rm median}$** are at literature/first-principles defaults
  ($\lambda=0.035$ from Bullock et al. 2001; $\sigma_{\ln j}=0.5$ picked,
  not derived) -- $\sigma_{\ln j}$ in particular was explicitly left as a
  free parameter to be explored, not assumed, per the paper's own Section
  2.
- **The $10^{12}\,M_\odot$/$3\times10^{13}\,M_\odot$ high-mass
  flattening** (confirmed to be $\zeta_{\rm ch}$, the cold/hot-mode
  transition, not an artifact) has not been checked against real trees at
  those masses in as much depth as the low-mass end.

## 5. 2026-09 session (cont.): full-grid $M_{\rm res}$ check; a $dz$ numerical instability found

### 5.1 $M_{\rm res}$ convergence: passed, cleanly, across the full grid

Ran the joint-grid methodology from `docs/RESOLUTION_CONVERGENCE.md`
against *this* calculation for the first time (Zhang-Hui trees,
$z_{\rm max}=25$, converging $M_{\rm seed,crit}$): $M_{\rm res}\in
\{10^3,10^4,10^5\}\,M_\odot$ (fixed absolute, not a fraction of $M_0$) at
fixed $dz=0.1$, across all 25 of the paper's halo masses ($N=100$ trees
per point, 75 points total, ~5.5 minutes). Result: **maximum spread
across the full $10^3$-$10^5\,M_\odot$ range is 0.02 dex anywhere in the
grid** -- well inside sampling noise at $N=100$, nowhere near the
fiducial ensemble's own $\sim$0.05-0.07 dex tree-to-tree scatter (Section
4, item 5). $M_{\rm res}=10^4\,M_\odot$ is safely converged with respect
to mass resolution across the entire $3\times10^{10}$-$3\times10^{13}\,
M_\odot$ range -- unlike the z=0 SHMR calculation (`RESOLUTION_CONVERGENCE.md`),
which found no ceiling even at the finest tested $M_{\rm res}$, this
quantity converges comfortably. Consistent with the mechanism identified
in Section 4, item 5: $M_{\rm seed,crit}$ is set by early-formation-time
competition (seed vs. Salpeter clock), which doesn't require resolving
anywhere near as much late-time substructure as $z=0$ $M_\star$ does.

### 5.2 $dz$: a genuine numerical instability found at the high-mass end, not yet fixed

The companion $dz$ check (same masses, fixed $M_{\rm res}=10^4\,M_\odot$)
is clean at low/fiducial mass but **non-monotonic and non-converging at
$M_{\rm halo}=3\times10^{13}\,M_\odot$**, checked directly at $N=1000$
(real statistics, not $N=100$ sampling noise) across five $dz$ values:

| $dz$ | median $M_{\rm seed,crit}$ | 16th-84th |
|---|---|---|
| 0.2 | $1.18\times10^7$ | $[1.07,1.30]\times10^7$ |
| 0.1 | $9.71\times10^6$ | $[8.97,10.58]\times10^6$ |
| 0.05 | $9.16\times10^6$ | $[8.36,9.99]\times10^6$ |
| 0.025 | $9.67\times10^6$ | $[8.61,11.01]\times10^6$ |
| 0.0125 | $1.21\times10^7$ | $[10.17,15.18]\times10^7$ |

Decreases then increases back toward the coarsest-$dz$ value, and the
16th-84th width itself roughly doubles from $dz=0.05$ to $dz=0.0125$ --
the signature of a real numerical instability, not "needs finer
resolution yet."

**Suspected cause**: `run_reservoir_paper()`'s `halo_growth_rate_mid` is
computed by finite-differencing the tree's own mass trajectory
(`(halo_mass[:,j]-hm_prev)/dt`) -- the fix, earlier this session, for a
separate bug where a *re-derived smooth* growth-rate estimate fed into
`uv_suppression()` mismatched a real tree's bursty actual rate and drove
the suppression formula negative/clipped (Section 4, item 2). A one-step
finite difference of a function with genuine discontinuities (merger
mass jumps) is itself noise-amplifying as the step shrinks: a merger
landing inside one fine step produces an enormous spike (dividing by a
tiny $\Delta t$) while adjacent steps read near-zero, rather than the
same jump being smoothed over a wider coarse step. This would hit the
most merger-active haloes/trees hardest, consistent with the high-mass
end being affected while low/fiducial mass is clean.

**Not yet fixed.** The original clipping bug is real and the trajectory-
based-rate fix was the right direction, but the naive single-step finite
difference needs replacing with something less noise-sensitive (e.g.
averaging the rate over a short window rather than a literal one-step
derivative) before $dz$ can be meaningfully convergence-tested at the
high-mass end. Do not report a fiducial $dz$ as converged for the full
grid until this is resolved -- the low-mass/fiducial-mass convergence
already demonstrated does not extend to $M_{\rm halo}\gtrsim10^{13}\,
M_\odot$.

## 6. 2026-09 session (cont.): the $dz$ instability resolved (GRUMPY rate estimator), the smooth-fit baseline abandoned

This closes out Section 5.2's "not yet fixed" status, but the resolution
is not the windowed-average fix speculated there -- that was tried next
and itself found to be a dead end, along with a smooth/merger-channel
split. Both are documented here as ruled-out approaches, not omitted, so
they are not retried.

### 6.1 Two more rate-estimator attempts, both ruled out

- **Fixed-physical-time trailing window** (`windowed_halo_growth_rate`,
  20 Myr): fixed $dz$-monotonicity, but swung the fiducial-mass
  $M_{\rm seed,crit}$ by 76x relative to the single-step-diff result,
  revealing that the earlier apparent $\sim$10% match to the paper's old
  2903 figure (Section 4, item 4) had been a coincidence of the
  since-fixed single-step-diff bug's own behaviour, not a real
  validation.
- **Smooth/merger-channel split** (`smooth_only_mass_trajectory`,
  applying $\zeta_{\rm UV}/\zeta_{\rm ch}$ only to a tree's
  `smooth_accretion` channel, bypassing suppression for merger-delivered
  mass): motivated by "a merger is an instantaneous mass transfer, not a
  rate" -- correct in spirit, but did not fix $dz$-convergence. Root
  cause traced into `foraois`: `_build_forest_flat_barrier_kernel`
  resolves at most one merger per $dz$ step, so the smooth/merger split
  itself is not $dz$-converged at fixed $M_{\rm res}$ -- finer $dz$
  keeps reclassifying more real growth from "smooth" into "merger," with
  no plateau. A `foraois`-side fix (allowing multiple resolved splits per
  step) was identified as the correct fix for *this* approach, but was
  superseded by 6.2 below before being implemented.

### 6.2 Root cause found: `uv_suppression`'s `ratio_term` needs a derivative real trees don't have

The user's question "aren't UV suppression and cold/hot-mode accretion
functions of mass, not rate?" prompted checking `reionization.uv_suppression`
directly. Its base term, $s(\mu_c,\omega)$, is the literature-standard
Okamoto et al. (2008) form (mass/$z$ only). But the codebase's own
`uv_suppression` (authored Aug 2024, `git log -S ratio_term`, no cited
derivation) adds a correction, `ratio_term`, containing $\dot M_{\rm h}$
in the denominator. Structurally this is the derivative of the
equilibrium relation $M_{\rm gas}=\bar f_b(M_{\rm h},z)M_{\rm h}$ --
literally requiring $dM_{\rm h}/dt$ to exist in the classical sense. A
real tree's raw main-progenitor mass has genuine jump discontinuities at
every resolved merger, where it does not.

The user identified the fix's actual source: Kravtsov & Manwadkar
(2022, GRUMPY)'s own published mass-accretion-history smoothing
procedure -- cubic-spline $\log M_{\rm h}$ vs. $\log t$ on the tree's
raw native-grid nodes, differentiate analytically, clip negative
$d\ln M_{\rm h}/d\ln t$ to zero (their own stated monotonicity
condition), re-spline, then evaluate
$\dot M_{\rm h}=(M_{\rm h}/t)\,d\ln M_{\rm h}/d\ln t$. Implemented as
`grumpy_halo_growth_rate()` in `paper_reservoir.py`, replacing both
`windowed_halo_growth_rate` and `smooth_only_mass_trajectory` (removed).

Because this is still exact interpolation through the tree's own raw
nodes, a merger event whose sampled step becomes small enough still
produces a locally huge derivative -- checked directly, $dz=0.0125$
produced a tree whose merger event landed in the single step immediately
preceding $z_{\rm anchor}$, giving a $>100\times$ $M_{\rm seed,crit}$
outlier. Per the user's explicit direction ("it's okay to force the user
to think" rather than add automatic suppression of bad behaviour that
doesn't require user input), **$dz$ is fixed at $dz=0.05$** as a
disclosed physical-resolution choice, analogous to $M_{\rm res}$, chosen
from an observed stable plateau ($dz=0.025$-$0.1$ agree to
$\lesssim0.1$ dex) rather than driven toward convergence.

### 6.3 Ablation: `ratio_term` turns out not to matter at all

Before committing to the GRUMPY fix's added complexity, an ablation
compared the full ratio_term-corrected $\zeta_{\rm UV}$ against the pure
mass-only Okamoto form, across the full $3\times10^{10}$-$3\times10^{13}
\,M_\odot$ grid ($N=200$): **the correction changes $M_{\rm seed,crit}$
by at most $+0.04$ dex** (fiducial mass) and $+0.00$ dex (high mass) --
indistinguishable from noise. Separately, $\zeta_{\rm ch}$ (cold/hot-mode,
purely mass-dependent, no rate at all) is the *dominant* physical lever
at high mass: forcing it to 1 swings $M_{\rm seed,crit}$ by $+2.56$ dex
(360x) at $3\times10^{13}\,M_\odot$. Conclusion: `ratio_term` was the
*only* place in the whole calculation where a rate is divided into
(`Mdot_b = f_b*halo_growth_rate` only ever multiplies by $dt$, giving a
bounded mass increment, never a divide-by-a-possibly-tiny-number) -- and
it turned out to be physically inconsequential. **Dropped** in favour of
`uv_suppression_mass_only()`, a new `paper_reservoir.py` function
implementing the plain literature form; `reionization.uv_suppression`
itself is untouched (still used by `gas_evolve.py`/`main.run_forest`'s
separate, general-purpose model, calibrated against a different smooth
trajectory where `ratio_term` was never a problem).

Full $N=1000$, 25-mass production grid re-run with the corrected
pipeline (GRUMPY rate for $\dot M_{\rm h}$, mass-only $\zeta_{\rm UV}$,
$dz=0.05$): every mass point resolves cleanly, no `never_ok` trees
anywhere, tree-to-tree scatter a smooth $0.09$-$0.18$ dex across the
whole grid.

### 6.4 The paper's old $M_{\rm seed,crit}=2903$ figure is invalid, and was never a tree-based target

The 2903 figure (paper's fiducial-mass smooth-trajectory boundary) was
traced directly in the paper's own commented-out draft history: it
predates the low-spin mechanism and the GRUMPY rate fix entirely, is
self-flagged in-draft as "not yet confirmed... not to be quoted as a
real feature," and -- decisively -- the paper's own prior tree-based
comparison (that era's model, 200 trees) already found the tree-based
fiducial answer was 5777, not 2903 (+0.30 dex shift). **2903 was never
the number a tree-based calculation was expected to reproduce**, even
historically. Do not use it to calibrate `eta_acc`/`epsilon_sf`/
`sigma_lnj` going forward; see [[ashvini-viability-theorem]] memory for
why pulling the ratio toward 2903 hits a degenerate knife-edge anyway.

### 6.5 The smooth-fit baseline itself is invalid: FMBK10-via-shooting disagrees with real trees by orders of magnitude

Attempting to reinstate *any* smooth-vs-tree comparison (needed for the
draft's "Assembly history matters" section) surfaced a much bigger
problem than $\zeta_{\rm UV}$: `shoot_smooth_halo_trajectory`'s
construction -- integrating Fakhouri, Ma & Boylan-Kolchin (2010)'s
*unconditioned* mean halo growth-rate fit forward from a shooting-chosen
$z=25$ seed to hit a target $z=5$ mass -- disagrees with the real
Zhang-Hui tree ensemble by **4-5 orders of magnitude** at high $z$. At
$z=12$, 0 of 1000 real tree realisations (fiducial mass) had even
assembled the tree's own $M_{\rm res}=10^4\,M_\odot$ floor, while the
smooth trajectory already implied $M_{\rm halo}\sim2\times10^9\,M_\odot$
there. Checked that this is not a mean/median skewness artefact: the
ensemble *mean* (including still-unformed trees as zero) is equally far
off.

Diagnosis: forward-integrating an *unconditioned* mean growth-rate fit
from a shooting-chosen seed is not the same statistical construction as
a mass accretion history genuinely *conditioned* on the target final
mass -- which is exactly what a Monte Carlo tree ensemble computes
directly. FMBK10's rate is also mildly super-linear in $M$
($\dot M\propto M^{1.1}$); the tiny shooting-chosen seed mass this
implies at $z=25$ ($e$-folding time $\sim$29 Myr, faster than the
$\sim$46 Myr Hubble time there) compounds explosively over many
$e$-foldings. Confirmed independent of this codebase: an unrelated
smooth-fit attempt in `MenonBaluPower2026`, using a different
growth-rate prescription (Krumholz & Dekel), hit the same qualitative
failure.

**Decision (explicit user instruction): drop the smooth-fit comparison
from the paper entirely, "for now."** Section 4.1's boundary and Section
4.6's "Assembly history matters" were rewritten to report only the
tree-ensemble result (median + 16th-84th scatter, $0.09$-$0.18$ dex, no
smooth baseline). `fig_smooth_vs_trees.png` was regenerated once as a
tree-only-envelope-vs-(now-invalid)-smooth-curve comparison, then its
LaTeX figure block commented out (not deleted) once the decision to drop
the comparison was made, in case a properly-conditioned baseline is
built later.

### 6.6 A better candidate exists (Correa et al. 2015) but is not ready

User supplied Correa, Wyithe, Schaye & Duffy (2015, MNRAS 450, 1514),
which derives $M(z)=M_0(1+z)^{af(M_0)}e^{-f(M_0)z}$ from Neistein et al.
(2006)'s EPS differential equation for the mass accretion history
*conditioned on final mass* -- structurally the correct object, unlike
FMBK10. Tested by directly integrating their underlying ODE (Eq. 6 of
that paper, not their closed-form fit) using `foraois`'s own
$\sigma(M)$/$D(z)$ tables (for exact consistency with the trees),
anchored at $z=5$ instead of their published $z=0$ normalisation, with
$q(M_0)$ read off their fitted $z̃_f$-$q$-$M_0$ relations evaluated at
the $z=5$ anchor mass (an approximation, since that fit was calibrated
with $M_0$ = mass at $z=0$).

Result: **dramatically better than FMBK10** -- ratios to the real tree
median of $0.9$-$1.6\times$ (not $10^4$-$10^5\times$) from $z=5.1$ down
to $z\approx7$, across the full $3\times10^{10}$-$3\times10^{13}
\,M_\odot$ grid. A first attempt (fixed-step RK4 integrating $dM/dz$
directly) produced negative masses at high $z$ from the equation's
genuine stiffness as $M\to$ small; fixed by integrating $d(\ln M)/dz$
instead (guarantees $M>0$) with `scipy.integrate.solve_ivp(method="LSODA")`
and a precomputed $D(z)$ interpolant (`get_linear_growth_and_collapse`
reruns a full trapezoid integration per call -- far too slow inside an
adaptive solver's RHS otherwise), stopping integration at the same
$M_{\rm res}=10^4\,M_\odot$ floor the trees use rather than forcing the
solver through a physically-meaningless sub-resolution singularity.

After the fix, the negative-mass bug is gone, but a **real, systematic,
mass-dependent error remains near each halo's formation epoch**
($z\gtrsim8$): over-prediction growing from $\sim5\times$ at
$3\times10^{10}\,M_\odot$ to $\sim20$-$28\times$ at
$3\times10^{13}\,M_\odot$. This trend is consistent with the $z=0$-vs-
$z=5$ anchor mismatch in the borrowed $q(M_0)$ relation, not leftover
numerical noise. **Conclusion, per explicit user direction**: this is a
genuinely promising direction (motivation for a possible follow-up paper
properly re-deriving $q$, $z̃_f$, and the constant $a$ for a general
anchor redshift, not just $z=0$), but is not ready to use as this
paper's baseline -- "using the trees in the actual calculation makes
more sense for a careful paper." No further work on this went into the
current paper; the trees-only result of Section 6.5 stands.

### 6.7 Fiducial `eta_acc`/`epsilon_sf`/`sigma_lnj` re-derived independent of 2903

With 2903 retired as a target (6.4), re-examined whether the fiducial
$\eta_{\rm acc}=0.005$, $\epsilon_{\rm sf}=0.02$, $\sigma_{\rm lnj}=0.5$
survive on their own merits. Conclusion: the *ratio* logic (viability
margin from $f_{\rm BH}=0.5$) was never target-fit and needed no change,
but $\epsilon_{\rm sf}$'s specific value had never been independently
justified at all.

User pointed to the underlying `Ashvini` model's own paper
(Menon, Balu & Power 2026, PASA, arXiv:2508.08363 -- `MenonBaluPower2026`):
its Table 1 lists $\epsilon_{\rm sf}=0.015$ as the model's own fiducial
star-formation efficiency, not $0.02$. Since $\epsilon_{\rm sf}$ is
inherited machinery (the star-formation reservoir equations are `Ashvini`'s
own, unchanged by this paper), matching its published fiducial is a much
stronger anchor than reaching for external star-formation-efficiency
literature -- **updated `epsilon_sf` default in `paper_reservoir.py`
from 0.02 to 0.015** (both the module default and every calibration
script). $\eta_{\rm acc}$ (nuclear BH accretion efficiency) and
$\sigma_{\rm lnj}$ (angular-momentum scatter) have no such precedent:
both are new mechanisms introduced by this paper, not part of the
`Ashvini` model `MenonBaluPower2026` describes, so they remain free
parameters justified only by the viability margin and (for
$\eta_{\rm acc}$) the qualitative expectation that nuclear BH accretion
is less efficient than star formation given the angular-momentum
transport bottleneck the paper's own physical picture already invokes.

New ratio: $\eta_{\rm acc}/\epsilon_{\rm sf}=0.005/0.015=0.33$ (was
$0.25$) -- still comfortable margin below the $f_{\rm BH}=0.5$
degenerate threshold (6.1's viability theorem). Full $N=1000$, 25-mass
production grid re-run with `epsilon_sf=0.015`: every point still
resolves cleanly (no `never_ok` anywhere), $M_{\rm seed,crit}$ roughly
halves at every mass point relative to the `epsilon_sf=0.02` run (as
expected -- a smaller ratio-to-threshold gap makes the boundary easier
to cross with a smaller seed), power-law slopes essentially unchanged
($\alpha\simeq0.92$ from $10^{11}$-$10^{12}\,M_\odot$, $\alpha\simeq0.04$
above), tree-to-tree scatter $0.09$-$0.19$ dex across the grid (was
$0.09$-$0.18$). All of Section 4.1's and Table `tab:trees`'s numbers,
and the abstract's scatter-range claim, updated to match.

## 7. Critical bug found and fixed: `critical_seed_paper`'s `chi_crit=1.0` default was never Eddington-limited at all

While building item 1a's spin-sensitivity band (three curves at
$\epsilon=0.057/0.1/0.32$), all three came back bit-for-bit identical --
$\epsilon$ was having literally zero effect on $M_{\rm seed,crit}$.
Traced to `critical_seed_paper`'s default `chi_crit=1.0`, in place for
this entire session (and, per its own prior docstring claiming this
"reproduces strict Eddington-limited growth exactly," likely long
before).

**The bug, verified directly** by calling
`bh_growth_step_slimdisk` with fixed supply/mass and varying only
`kappa_edd`: at `chi_crit=1`, growth equals the raw, uncapped supply
rate `A_bh` regardless of `kappa_edd` (i.e. regardless of $\epsilon$).
Structurally: the three-regime closed form has boundaries
$M_1=A_{\rm bh}/(\kappa_{\rm edd}\chi_{\rm crit})$ and
$M_2=A_{\rm bh}/\kappa_{\rm edd}=M_1\chi_{\rm crit}$; at $\chi_{\rm
crit}=1$, $M_1=M_2$ exactly, so the genuine Eddington-limited
(exponential, $\kappa_{\rm edd}$-dependent) regime between them has zero
width, and the two surviving branches (the "throttled" cap
$A_{\rm bh}/\chi_{\rm crit}$ and the "gas-supply-limited" branch
$A_{\rm bh}$) are identical -- growth is $A_{\rm bh}$ at every mass,
completely bypassing any Eddington cap. `bh_growth_step_slimdisk`'s own
docstring and test suite already state the correct limit
(`test_slimdisk_reduces_to_pznk11_two_regime_model_as_rcrit_to_infinity`,
using `r_crit=1e12` to represent "strict Eddington, no throttling"): the
bug was entirely in `critical_seed_paper`'s (and the analogous
general-purpose `critical_seed.critical_seed`'s -- same bug, same wrong
reasoning in its own docstring, **not yet fixed, flagged separately**)
choice of `chi_crit=1.0` as its own "strict Eddington" default.

**Fix**: added `STRICT_EDDINGTON_CHI_CRIT = 1e12` to `paper_reservoir.py`
(matching the test suite's own convention for this limit) and changed
`critical_seed_paper`'s default from `chi_crit=1.0` to this constant.
Both the module's own docstrings and `run_reservoir_paper`'s docstring
corrected to state the right reasoning.

**Impact, checked by rerunning the full $N=1000$, 25-mass production
grid**: moderate, not catastrophic. Medians shift up by roughly
5-20% at every mass point (e.g. fiducial $5.335\times10^{10}\,M_\odot$:
$1.56\times10^5\to1.76\times10^5\,M_\odot$), boundary stays monotonic,
every point still resolves cleanly (0 failures anywhere). Tree-to-tree
scatter widened noticeably, from $0.09$-$0.19$ dex to $0.12$-$0.31$
dex -- genuine Eddington-limited growth is exponential and legitimately
more sensitive to exactly when in cosmic time a tree's growth window
opens, so more scatter under the corrected physics is expected, not a
red flag. Power-law slopes: $\alpha\simeq0.89$ ($10^{11}$-$10^{12}
\,M_\odot$, was 0.92) and $\alpha\simeq0.09$ (above, was 0.04) --
similar shape, flatter high-mass plateau less extreme than before.
Section 4.1, Table `tab:trees`, `fig_critical_boundary.png`, and every
scatter-range mention in the abstract/discussion/conclusions updated to
match.

**Not yet done**: `critical_seed.critical_seed` (general-purpose,
`main.run_forest`-based, used outside this paper) has the identical bug
and has NOT been fixed -- flagged for a separate task, since fixing it
touches `main.py`'s own test suite and any other analysis relying on its
current (wrong) default. Also not yet done: item 1a's spin-sensitivity
band itself needs rerunning against the now-fixed `chi_crit` default (it
should show real $\epsilon$-dependence now); everything else on the
plots-needed list that calls `critical_seed_paper` (items 1a, 1b, 2, 5,
6, 7) needs the same rerun before being trusted.

## 8. Item 1a (spin-sensitivity band) rerun against the fixed chi_crit: confirms real epsilon-dependence, and a significant seed-channel finding

Rerunning `item1a_spin_band.py` against the fixed `STRICT_EDDINGTON_CHI_CRIT`
default gives three properly-separated curves (chaotic floor
$\epsilon=0.057$, fiducial $\epsilon=0.1$, Bardeen coherent limit
$\epsilon=0.32$), a factor of $\sim4$ apart at the fiducial halo mass and
holding roughly constant across the full mass range -- confirms the fix
is doing the right thing (this was bit-for-bit identical across all three
before the fix). New figure: `fig_sensitivity_band_spin.png`, replacing
the old `fig_sensitivity_band_v2.png` (smooth-fit, pre-correction
integrator, with a leftover unresolved `[CBP: ...]` editorial note baked
into its live caption -- found and removed this session, see below).

**Seed-channel implication, computed precisely** (log-normal channel
distributions: Pop III $\mu=\log_{10}150$, $\sigma=0.35$ dex;
runaway/NSC $\mu=\log_{10}2000$, $\sigma=0.3$ dex; direct collapse
$\mu=\log_{10}2\times10^5$, $\sigma=0.3$ dex -- all pre-existing paper
values, unchanged) against the corrected boundary, across all 25 halo
masses and all three spin cases:

- **Pop III remnant**: 100% require super-Eddington growth, every mass,
  every spin case. Unchanged from the (also correct) pre-fix claim.
- **Runaway/NSC collision**: 100% require super-Eddington growth, every
  mass, every spin case. **This directly contradicts** the paper's
  prior text, which called this a "transition region" depending on halo
  mass -- it never was, once Eddington limiting is genuinely enforced.
- **Direct collapse**: only escapes the super-Eddington requirement in a
  narrow low-mass window ($5.3$-$7.1\times10^{10}\,M_\odot$) at
  chaotic-floor/fiducial spin; at the coherent-accretion-limit end, the
  super-Eddington fraction never drops below 68% even at the lowest
  halo mass tested. **Also contradicts** the paper's prior claim that
  this channel lies "substantially closer to or above the boundary for
  much of the halo-mass range."

Overall: once the Eddington cap is genuinely enforced, the light-seed
problem is close to universal across all three standard formation
channels over most of the halo-mass range -- a materially different,
and more striking, result than the pre-fix draft claimed. Section 4.6
("What the boundary implies for seed formation channels"), the abstract,
Discussion \S5.2, and Conclusions item 4 all rewritten to match, per
explicit user sign-off after reviewing the exact crossing-point numbers
above.

Also removed while touching this figure: the old `fig_sensitivity_band_v2.png`
figure block had a live, unresolved `\textbf{[CBP: ...]}` editorial note
sitting in its actual caption text (found while reviewing a compiled PDF
of the paper) -- a leftover from a coauthor's own earlier editing pass,
predating this session entirely, unrelated to anything touched in
Sections 4-9 of this catalogue.

## 9. Plots-needed list fully cleared; Eq. viability and Section 4.3 fixed

Ran every remaining item on `plots-needed-2026-09-12.md` (1b, 2, 4, 5, 6,
7, 8, 9) against the fixed `chi_crit` pipeline. Full results and a running
issues list are in the session transcript; key outcomes:

- **Item 1b** (nuclear-supply band): $M_{\rm seed,crit}$ varies by
  roughly two orders of magnitude across $\eta_{\rm acc}=10^{-3}$-$10^{-1}$
  -- directly contradicts Section 4.3's old (pre-revision, self-flagged-as-
  unverified) claim that supply efficiency doesn't matter.
- **Item 2** (uncertainty budget): spin-spread panel built on a synthetic
  tree-median trajectory gives identical results for two different spin
  choices that item 1a's proper per-tree ensemble clearly separates --
  the synthetic-median-trajectory construction likely understates true
  spin-driven scatter; not yet reconciled with item 1a's more trustworthy
  per-tree approach.
- **Item 4** (feedback scan): $M_{\rm seed,crit}$-vs-$\epsilon_f$ trend
  looks legitimate (57x range, not flat); the outflow/inflow
  crossing-redshift diagnostic is broken (returns the same z regardless of
  $\epsilon_f$) -- not yet debugged.
- **Item 5** ($f_{\rm BH}$ sensitivity): clean, monotonic, no issues.
- **Item 6** ($\chi_{\rm crit}$ reconciliation): reassuring -- $\chi_{\rm
  crit}\geq10$ (the paper's own general fiducial) already gives results
  identical to $\chi_{\rm crit}=10^{12}$; only $\chi_{\rm crit}\lesssim3$
  differs. `STRICT_EDDINGTON_CHI_CRIT=1e12` wasn't an arbitrary choice --
  the S-regime genuinely stops mattering well below that value.
- **Item 8** (PCH08 cross-check): agrees with Zhang-Hui to 0.02-0.14 dex,
  inside tree-to-tree scatter -- assembly-history result isn't an
  algorithm artifact.
- **Item 9** (representative trajectories): added `A_bh`/`kappa_edd` as
  diagnostic outputs to `run_reservoir_paper`'s returned dict (additive,
  doesn't change existing behaviour) to enable S/E/G regime classification
  and cold/hot-mode-crossing annotation. Initial concern about an
  oscillating $M_\star$ curve was checked directly and is a plotting
  artifact (raw array has zero negative diffs, fully monotonic as the
  model's own equations guarantee) -- not a real bug.

**Fixed in the paper itself** (Section 2.5's Eq. viability derivation, and
Section 3's/4.3's dependent text): the "viability theorem" was itself an
artifact of the same `chi_crit=1` bug (Section 8) -- with $\chi_{\rm
crit}$ fixed, there is no simple closed-form necessary condition on
$\eta_{\rm acc}/\epsilon_{\rm sf}$ any more (checked directly up to
ratio=6.67, still well-defined). Eq. viability and its derivation were
removed and replaced with a direct-verification statement; Section 3's
$\eta_{\rm acc}$ justification updated to drop the now-false "comfortable
margin from a knife-edge" framing while keeping the still-valid
angular-momentum-bottleneck qualitative argument; Section 4.3 rewritten
(new title too) using item 1b's actual tree-based numbers, replacing the
old, self-flagged-as-unverified insensitivity claim.

**Still not propagated into the paper**: items 2, 4, 5, 6, 7, 8, 9's
figures/numbers exist but aren't yet written into the text (only 1a and
1b are integrated so far). Item 2's methodology caveat and item 4's
crossing_z bug should be resolved before those go in.

## 10. Item 4 crossing_z bug fixed; item 2 re-derived properly; everything propagated into the paper

**Item 4's crossing_z bug**: traced to the same root cause as the gas
self-regulation oscillation (Section 8 below) -- `mdot_out` is exactly
zero for the majority of trees at most individual timesteps, so the
cross-tree median used for this diagnostic flips between 0 and spikes for
reasons unrelated to `epsilon_f`. Fixed by switching to a single
representative tree (closest individual `M_seed,crit` to the ensemble
median, same convention as item 9) run on a finer `N_STEPS=4001`
reservoir grid. `crossing_z` now varies genuinely with `epsilon_f`
(range 8.4-9.5), though still noisy since a different representative
tree is selected at each `epsilon_f`; the `M_seed,crit`-vs-`epsilon_f`
trend itself (the more important panel) was always clean.

**Gas self-regulation oscillation** (found while diagnosing the above):
`run_reservoir_paper`'s explicit forward-Euler gas update exhibits a
period-2 limit cycle whenever the local depletion timescale is shorter
than the fixed reservoir step size -- pervasive (100% of trees show it
somewhere, ~10% of all steps, checked at two halo masses). User confirmed
this class of bursty self-regulation is a known, expected feature in this
model family (cf. Menon \& Power 2024, Menon, Balu \& Power 2026), albeit
via an explicit feedback delay there rather than the direct depletion
coupling here. Checked directly: `M_seed,crit` converges to a few percent
across a 10x reservoir-grid refinement (400->4000 steps) -- not extreme,
does not affect the paper's results. Disclosed in the Appendix
("Numerical verification").

**Item 2 re-derived**: the original version's "spin spread" panel used a
synthetic tree-median trajectory that erased a real spin distinction item
1a's proper per-tree ensemble shows clearly. Fixed by reusing item 1a's
own real per-tree-ensemble medians (chaotic/fiducial/coherent) directly
for the 3-point spin spread, alongside a freshly-regenerated raw
1000-tree distribution at fixed fiducial spin for the violin panel.

**Everything propagated into the paper** (compiles cleanly, 16 pages,
zero errors/undefined references):
- Item 9 (representative trajectories, with S/E/G regime shading and
  cold/hot-mode crossing annotations) -- new subsection after \S4.1.
- Item 6 (chi_crit reconciliation) -- new subsection after \S4.2; shows
  chi_crit=10 (paper's own general fiducial) is already indistinguishable
  from the formal strict-Eddington limit.
- Item 4 (feedback scan) -- new subsection after \S4.3.
- Item 8 (PCH08 cross-check) -- new paragraph within \S4.5.
- Item 5 (f_BH=0.1/0.9 sensitivity) and item 2 (uncertainty budget) -- new
  subsections after \S4.6 (seed channels), before the results summary.
- **Section 4.7 (results summary), the abstract, \S5.1 (critical seed
  mass as growth criterion), and Conclusions items 2 and 5 all rewritten**:
  these repeatedly asserted the old "gas supply/nuclear accretion
  efficiency doesn't matter, black hole remains supply-rich" claim,
  which items 1b and 4 now directly contradict (nuclear-supply efficiency
  and feedback strength are both substantial, comparably-important levers
  on the boundary, factors of tens to ~100). Replaced throughout with the
  corrected, multi-factor picture: supply efficiency, feedback, radiative
  efficiency, and assembly history are all comparably important, none
  dominates.
- Also removed, while touching the Conclusions: a second stray
  Markdown code-fence artifact (harmless to compilation in this
  position, but a leftover nonetheless) sitting just before
  `\end{enumerate}`.

## 11. `critical_seed.py`'s chi_crit bug fixed; full code audit

Fixed `ashvini/critical_seed.py`'s identical `r_crit=1.0` bug the same
way as `paper_reservoir.py`'s (added `STRICT_EDDINGTON_R_CRIT=1e12`,
corrected the default and docstring). Confirmed via `grep` that no test
suite and no other call site exists anywhere in the repo for this
function -- unlike `critical_seed_paper`, it has not been stress-tested
against real data since the fix.

**Audit performed**, all clean except one finding:
- Searched the whole `ashvini/` tree for any other `chi_crit=1`/`r_crit=1`
  default -- none found; both known instances are now fixed.
- Confirmed `main.run_forest`'s own default (`r_crit=8.0`, from
  `run_params.yaml`) was never affected by this bug at all -- the general
  Ashvini model was never in the degenerate regime.
- Full test suite: 111 passed, 0 failed, only 2 pre-existing unrelated
  warnings (M_res convergence, not from today's changes).
- `pyflakes` on both `paper_reservoir.py` and `critical_seed.py`: clean,
  no unused imports/names.
- Searched for dangling references to the abandoned
  `smooth_only_mass_trajectory`/`windowed_halo_growth_rate` functions
  (removed earlier this session): none found anywhere.
- **Found and fixed a real, separate stale-default bug**:
  `run_reservoir_paper`'s own function-signature default was still
  `eta_acc=0.01`, never updated when the paper's fiducial moved to
  `0.005` (every script this session passed it explicitly, so this never
  showed up in any result, but a bare call would have silently used the
  wrong value). Fixed to `eta_acc=0.005`. Cross-checked every other
  parameter default in that signature against the paper's own stated
  fiducial list (Section 3) -- all others already consistent.
- Reran the full $N=1000$, 25-mass production grid after all of today's
  fixes: matches the previously-reported table to $\leq1.6\%$ at every
  point (consistent with `build_forest_numba`'s already-documented
  non-reproducible parallel RNG, not a real change), and every mass
  point still resolves cleanly (0 failures anywhere). `critical_seed.py`'s
  fix, being a separate, previously-unused module, had zero effect on
  `paper_reservoir.py`'s results, as expected.
- Recompiled the paper (unaffected by today's code-only changes): clean,
  16 pages, 0 errors.

**Follow-up raised, then implemented same day**: user asked whether a central
fiducial-parameters config (following the existing `run_params.py`/
`run_params.yaml`/`PARAMS`-dataclass pattern already used for the general
model) would be worth adding for `paper_reservoir.py`, to prevent this
exact class of default-drift bug going forward. Agreed and implemented:

## 12. Central fiducial-parameters config for `paper_reservoir.py`

Added `ashvini/paper_reservoir_params.py` (a `PaperReservoirParams`
dataclass, `load_paper_reservoir_params(config_file=None)`, `print_config`,
and a module-level `PAPER_PARAMS = load_paper_reservoir_params()` singleton)
plus `paper_reservoir_params.yaml` at the project root -- the exact same
pattern as `run_params.py`/`run_params.yaml`/`PARAMS`, applied to this
paper-specific model instead of the general one.

`paper_reservoir.py`'s module-level `*_FIDUCIAL` constants
(`F_B_FIDUCIAL`, `M_HOT_FIDUCIAL`, `PHI_COLDHOT_FIDUCIAL`,
`LAMBDA_MEDIAN_FIDUCIAL`, `DELTA_VIR_FIDUCIAL`, `OMEGA_M_FIDUCIAL`,
`OMEGA_LAMBDA_FIDUCIAL`, `STRICT_EDDINGTON_CHI_CRIT`) now read from
`PAPER_PARAMS` rather than being hardcoded independently; `run_reservoir_paper`
and `critical_seed_paper`'s function-signature defaults (`epsilon_sf`,
`chi_crit`, `epsilon_f`, `eta_acc`, `R_nuc_pc`, `compaction_boost`,
`sigma_lnj`, `gas_mass0`, `stars_mass0`, `f_bh`) do the same. Every caller
can still override any of these per-call via ordinary keyword arguments
(unchanged calling convention, needed for the sensitivity scans) -- only
what "fiducial" means when no override is given now lives in one place.

Verified: every new default was checked via `inspect.signature` against
the value it replaced -- byte-identical in every case (e.g. `eta_acc`
correctly comes out as `0.005`, matching the fix from earlier today, not
the old stale `0.01`). Full test suite (111 tests) and `pyflakes` on the
new and modified files: both clean.

`ashvini/critical_seed.py`'s own `STRICT_EDDINGTON_R_CRIT` was left as a
plain module constant, not folded into this config: it is a mathematical
sentinel ("no cap ever engages"), not a tunable physical fiducial, and
that module already sources its own seeding parameters from the separate,
pre-existing general-model `run_params.PARAMS`.

**Not actioned** (flagged, out of scope for this task): `run_params.py`'s
own pre-existing `SlimDiskParams.eta_acc: float = 0.01` (the general
`main.run_forest` model's slimdisk config, separate from
`paper_reservoir.py`) is also stale relative to the paper's `0.005`
fiducial -- not yet investigated for whether it's consequential anywhere
in the general model, since no caller of `main.run_forest` with
`growth_model="hobbs_slimdisk"` has been audited for actually relying on
the bare default rather than passing `eta_acc` explicitly.

## 13. Closing the last two loose threads (2026-09-15)

Asked directly whether any code updates to `Ashvini`/`foraois` were still
needed given everything found this session. Two were:

- **Fixed** `run_params.py`'s stale `SlimDiskParams.eta_acc: float = 0.01`
  default (flagged in \S12, not actioned then) -- changed to `0.005` to
  match the paper's corrected fiducial, and updated the matching commented
  example in `run_params.yaml`. Confirmed dormant in practice (default
  `growth_model` is `"pznk11_freefall"`, and the `slimdisk:` block in
  `run_params.yaml` is fully commented out), but a live trap for the next
  `hobbs_slimdisk` caller who doesn't override `eta_acc` explicitly.
- **Added regression coverage** for the `chi_crit`/`r_crit=1.0` bug class
  (\S7, \S11), which had no test guarding against it recurring:
  `test_slimdisk_r_crit_equal_one_is_not_strict_eddington` in
  `tests/test_black_holes_growth_slimdisk.py` (checks `r_crit=1` gives the
  raw uncapped supply rate, not Eddington-capped growth), and a new
  `tests/test_critical_seed_defaults.py` asserting
  `critical_seed.critical_seed`'s and
  `paper_reservoir.critical_seed_paper`'s own default `r_crit`/`chi_crit`
  are pinned to their `STRICT_EDDINGTON_*` constants, not `1.0`.

**Found and fixed in passing**: the `eta_acc` default change exposed a
pre-existing, unrelated fragility in
`test_run_forest_epsilon_f_override_changes_gas_mass` -- it asserted on
`gas_mass`'s *final-step* value, which turns out to settle into a
self-regulating quasi-equilibrium (continued inflow refills whatever the
AGN wind removes) that is essentially independent of `epsilon_f` along
that test's specific synthetic halo trajectory, checked directly across a
wide range of `eta_acc` and step counts (not a numerical-oscillation
artefact). `bh_mass` and `stars_mass` *do* respond to `epsilon_f` along
the same trajectory. Renamed the test to
`test_run_forest_epsilon_f_override_changes_bh_mass` and pinned `eta_acc`
explicitly rather than relying on the module default. Also corrected a
factually wrong comment in `test_run_forest_r_crit_override_changes_bh_mass`
that called `r_crit=1.0` "strict Eddington" (same bug pattern as \S7,
but only in a comment -- the test's own assertion never depended on that
being true). Full suite: 114 passed (113 pre-existing + 1 renamed, +2 new
regression tests), `pyflakes` clean.

**foraois**: no changes identified or made. It was used this session
(tree generation for the production grid) but not modified, and nothing
in the corrected physics implies a bug in it.
