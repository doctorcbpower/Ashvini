# Ashvini as the main codebase: superset design (draft for agreement, 2026-11-02)

Status: proposal only. No code has been changed. Nothing below is implemented until the decisions in section 6 are agreed.

Direction (decided by the user): Ashvini is the main codebase. Stochastic_Occupation (`stochocc`) becomes a set of components and an event integrator inside Ashvini, and its analysis scripts later call Ashvini. `stochocc` is then kept unchanged as the record of the published results.

## 1. Pre-flight model audit

| # | Item | Today in Ashvini | Today in stochocc | Proposal |
|---|---|---|---|---|
| 1 | Star formation law | eps_ff M_g / (0.141 / H) | eps_ff M_g / t_dyn, t_dyn = 0.2047 / H | One `SFLaw` component. Same law up to eps_ff x 0.2047/0.141 = 1.452. Keep both as named presets. |
| 2 | Feedback timing | one lag, 15 Myr, fixed | per-progenitor lifetimes (Georgy+13) | `DelayKernel`: delta function (Ashvini) or IMF-weighted lifetime distribution (stochocc smooth limit). |
| 3 | Feedback amplitude | momentum mass loading, eta(M_h, z, Z*) | energy wind model min(M_g, 2 eps E / v_esc^2), or momentum | `WindModel` component with three members: `momentum_ml`, `momentum_event`, `energy`. Parity note shows `momentum_ml` equals `momentum_event` with f_p = 0.7625 for 8-100 Msun. |
| 4 | Discreteness | none (smooth) | Poisson-sampled supernovae, smooth control | `Sampler` option on the event integrator: `smooth`, `stochastic`. |
| 5 | Star formation cooling threshold | none | T_vir >= 1e4 K, or H2 mass cooling threshold, plus post-reionisation switch | `CoolingThreshold` component, default `always_on` to keep Ashvini's results unchanged. |
| 6 | UV suppression | rate-dependent Okamoto form (own audit flags it) | equilibrium form | `Inflow` takes a `UVSuppression` object. Both forms are available and named; equilibrium form is the shared base in `reionization_physics.py`. |
| 7 | Chemistry | gas and stellar metals, dust, with known issues | metals and dust per event | One `Chemistry` component after the known dust and Z_IGM issues are fixed. |
| 8 | Black holes, seeding, AGN | several growth laws, three seeding channels | none | Unchanged, each behind a `BlackHole` interface. |
| 9 | Reservoirs | one-zone (general chain) and MVM galaxy plus nucleus | one-zone | `Reservoirs` interface with `one_zone` and `mvm`. MVM stays unchanged and identical in every digit (decision D3). |
| 10 | Time stepping | vectorised over haloes on the tree grid | one halo at a time, events at their own times | Two time integrators on one `State` and one `Config`. |
| 11 | Extra stochocc features | none | inflow lag, recycling, stored explosion energy, f_access, M_eps_cut, wind travel time | Ported as optional components, each off by default, after the base agreement tests pass. |

Not in scope: new physics, new figures, new scans, manuscript changes.

## 2. Bugs to fix first (reproduced, 2026-11-02)

1. `supernova.delay_time: 0` with `type: delayed` gives exactly the no-feedback result. The lookback index equals the current step, and the SFR of the current step is not yet computed. It should either mean instantaneous feedback or be rejected.
2. `type: "no"` still produces dust and still removes dust with the wind; `type: "instantaneous"` still drives dust from the delayed SFR. Metals and dust must follow the same feedback switch.
3. Dust supernova rate is Gamma / dt_step (grid-dependent), with no M_d / M_g factor.
4. Stellar mass uses gas mass at the start of the step (first-order lag, 1.1% between dz = 0.05 and 0.005).
5. `Z_IGM = 1e-3` is documented in Z_sun units but used as a mass fraction.

Each fix changes results. They need agreement, a regression test that keeps the old behaviour for comparison, and a note in the change log. The unchanged MVM production run does not use the general chain, so it is unaffected.

### Status of the bug fixes (2026-11-02)

All five are fixed in the working tree (uncommitted), with tests in `tests/test_chemistry_and_feedback_switch.py` and the stored fixture baseline in `tests/test_run1.py` recaptured (old values recorded there). Also fixed on the way: `supernova.type: "instantaneous"` was switched off until `delay_time` had elapsed. Metallicity units: `Z_IGM_Zsun`, seeding `Z_gas_max_Zsun`, and the `metallicity_function` midpoint (0.1 Z_sun) are all in Z_sun; old keys are rejected with an explanatory error. New: `supernova.imf`, `m_sn_min`, `m_sn_max`. The unchanged MVM is untouched (its mass loading uses Z = 0, for which the new function is identical).

Effect on the 1e10 Msun fixture halo at the last step: gas +1.4%, stars +0.5%, gas metals -26%, stellar metals -45%, dust -13%, BH mass unchanged. Analytic 1e7-1e9 Msun halos (no BH, no UV): stars +1.2%, gas metals -55% to -78%, dust between +26% and -7%.

## 3. Architecture

- `Config`: one dataclass tree, passed explicitly. No module-level reads of `PARAMS`; `set_params` is replaced by building a `Config`. `run_params.yaml` becomes one way to build a `Config`.
- `State`: a dataclass of arrays (gas, stars, metals, dust, BH, plus the event queue for the event integrator).
- Components, each a small class with one documented method and a one-paragraph physics docstring:
  `SFLaw`, `Inflow`, `UVSuppression`, `CoolingThreshold`, `DelayDistribution`, `WindModel`, `Chemistry`, `BlackHole`, `Seeding`, `AGNFeedback`.
- Registry: `Config` names components as strings (`wind_model: "energy"`). A student adds a component by writing a class and registering it in one place.
- Time integrators: `GridIntegrator` (today's `run_forest`, vectorised) and `EventIntegrator` (ported from `stochocc.engine`). Both call the same components.
- Cosmology and constants: one place (`utils.py`, `constants.py`). The Planck15 use in `reionization.py` goes away or is made an explicit option.

## 4. Student-facing pieces

- `docs/GUIDE.md`: which of the four existing model families to use, and which to ignore (`reservoir_stock_premvm.py`, `paper_reservoir.py` integrator).
- A notebook that swaps one component at a time and plots the effect.
- Each component page lists its equation, parameters, defaults and the paper it comes from.

## 5. Order of work and tests

1. Reference-output tests: store the current outputs of the general chain and the MVM before any change.
2. Fix bugs in section 2, one commit each, each with a regression test.
3. Introduce `Config` and `State`, with no behaviour change.
4. Extract components from `run_forest`; the reference-output tests must still pass.
5. Add the `DelayDistribution` and `WindModel` options and reproduce the stochocc smooth control through `GridIntegrator` (agreement target: the existing A1, B1, E1 results in `ashvini_relation.md`).
6. Port `EventIntegrator`; target: stored stochocc results from its test suite.
7. Port the optional stochocc features one at a time.
8. Guide and notebook.

## 6. Decisions needed

- D1. DECIDED (user, 2026-11-02): default free-fall time stays 0.141/H; the Stochastic_Occupation law is a named preset (eps_ff = 0.02178). The supernova upper mass limit becomes a parameter (default 40 Msun in the general chain, 100 in the Stochastic_Occupation preset).
- D2. DECIDED (user, a main co-author, 2026-11-02): fix all five bugs. Bug 5 (user, 2026-11-02): intended value is 1e-3 Z_sun, probably muddled to a mass fraction in MBP26. Policy: every metallicity in the config carries its unit in its name (`..._Zsun` or `..._massfrac`), one `Z_SUN = 0.015` lives in `constants.py`, and the resolved values in mass fraction are printed in every run header. Same audit to be applied to the seeding `Z_gas_max` thresholds and the `metallicity_function` midpoint (m = 0.1), pending confirmation.
- D3. DECIDED (user, 2026-11-02): the MVM stays unchanged next to the framework. `reservoir_stock.py` (and the superseded `reservoir_stock_premvm.py`, `paper_reservoir.py` integrator) are not modified; the provenance hash stays valid. They may import shared constants, cosmology and reionisation physics. Rebuilding the MVM from components is a separate, later step with a round-off agreement test against the kept unchanged copy.
- D4. DECIDED (user, 2026-11-02): the Stochastic_Occupation `legacy` time scheme is dropped; only `v2` (exact linear-ODE substep, SFR-weighted event time stamp) is ported. Any early result that used `legacy` stays reproducible from the kept unchanged `stochocc` repo.
- D5. DECIDED (user, 2026-11-02): one package, with an event-integrator subpackage (layout in section 7). Stochastic_Occupation stays as a thin analysis repo that imports Ashvini, and is kept unchanged as the record of its published results once Ashvini reproduces them.

## 7. Concrete layout and milestones (proposal, for agreement before any code)

```
ashvini/
  config.py            Config dataclass tree + loaders (yaml, dict); replaces PARAMS/configure.set_params
  state.py             State dataclass (arrays), event queue for the event integrator
  constants.py  utils.py  imf.py  reionization_physics.py     (shared, already exist)
  components/
    star_formation.py  SFLaw: rate(state, t) -> dM*/dt per halo; presets "ashvini" (0.141/H), "dyn200" (0.2047/H)
    inflow.py          Inflow: mass supply per step; wraps a UVSuppression
    uv.py              UVSuppression: "none", "okamoto_eq", "ashvini_rate"
    cooling_threshold.py  CoolingThreshold: "always", "tvir", "mass", plus post-reionisation temperature
    delay.py           DelayKernel: "fixed" (Ashvini), "lifetimes" (IMF-weighted), "none"
    wind_model.py         WindModel: "momentum_ml", "momentum_event", "energy"; limit at the gas present and reachable fraction options
    chemistry.py       Chemistry: metals and dust (the corrected equations)
    blackholes.py      BlackHole: "pznk11", "slimdisk"; Seeding; AGNFeedback
    optional.py        inflow lag, recycling, energy bank, wind travel time (off by default)
  integrators/
    grid.py            GridIntegrator: today's run_forest, vectorised over haloes
    events.py          EventIntegrator: ported from stochocc.engine (v2 only)
  reservoirs/          (frozen, unchanged) reservoir_stock.py, reservoir_stock_premvm.py, paper_reservoir.py
  main.py              thin: run(), run_forest() and run1_scalar() kept as compatibility wrappers
```

Interface rule: each component is a small class with one documented method, a one-paragraph physics docstring with the equation and its source, and a `name` registered in one place (`components/__init__.py`). A student adds a component by writing the class and adding one line to the registry. No component reads a module-level global; everything comes from the `Config` it is given. Validation uses exceptions with messages, not `assert`.

Milestones, each ending with the full test suite green and the stored outputs unchanged unless stated:

- M0. Reference outputs: store general-chain outputs (several halos, all three feedback types, BH on and off, UV on and off) and the MVM output hash, as tests.
- M1. `Config` and `State` introduced; `PARAMS`-based modules read from a `Config` passed explicitly; `configure.set_params` kept as a thin wrapper. No behaviour change.
- M2. Extract components from `run_forest` (SFLaw, Inflow, UV, Wind model, Chemistry, BlackHole) behind the GridIntegrator. No behaviour change; M0 stored reference numbers hold to round-off.
- M3. `DelayDistribution` and `WindModel` options; reproduce the Stochastic_Occupation smooth control through the GridIntegrator. Acceptance: A1, B1, E1 of `ashvini_relation.md` (stars within 1-4%, metals within 2%).
- M4. EventIntegrator (v2 only). Acceptance: stochocc's stored results for the smooth sampler to round-off, and ensemble statistics for the stochastic sampler within sampling error.
- M5. Optional stochocc features, one commit each: cooling threshold, inflow lag, recycling, stored explosion energy, f_access, M_eps_cut, wind travel time.
- M6. Student guide and notebook; archive the superseded paths in a `reservoirs/` note.

Not done in this plan: new physics, new figures, changes to the MVM, a decision on rebuilding the MVM from components.

## 8. Status of milestones (2026-11-02)

- M0 done and committed (`f73cb35`): `tests/test_reference_outputs.py` stores 16 parameter cases x 3 mass bins, the scalar reference integrator, the MVM outputs and the SHA-256 of `reservoir_stock.py`; stored reference outputs in `tests/reference_outputs/`.
- M1 done in the working tree (not committed): `ashvini/config.py` (`Config` = the existing `Params` tree, `DEFAULT`, `override`, `validate`, `load_config`), `ashvini/state.py` (`State`), and every general-chain module takes `cfg=None` (explicit Config, else `config.DEFAULT`). The import-time copies (`main.e_ff`, `main.t_d`, `agn.eta_agn`, `bh_growth.e_bh`, `reion.z_rei`, ...) are gone; assigning to a retired name raises an error naming the Config path. `main.run_forest/run1/run1_scalar/run` accept `config=`. `configure.set_params` now changes `config.DEFAULT` and gained `m_sn_max` and `agn_delay_time`. `critical_seed` works on a private copy of the Config instead of mutating and restoring. The M0 stored reference outputs are unchanged by M1 (only `apply_case` was rewritten), which is the no-behaviour-change check.
- Left as they were on purpose (kept unchanged, D3): `reservoir_stock.py`, `reservoir_stock_premvm.py`, `paper_reservoir.py` and `black_holes_growth_slimdisk.py` keep their import-time defaults (`PAPER_PARAMS`, `bh.slimdisk`); the general chain passes the slim-disc parameters explicitly from the Config.

- M1 committed (`2c77c5c`).
- M2 done in the working tree (not committed): `ashvini/components/` (`ode`, `star_formation.FreeFallStarFormation`, `inflow.CosmicInflow`, `wind_model.MomentumLoading`, `delay.SupernovaDelay`, `gas.gas_and_stars_step`, `chemistry.MetalsAndDust`, `blackholes.PZNK11BlackHole` / `SlimDiscBlackHole`) and `ashvini/integrators/grid.GridIntegrator`, which holds the components as replaceable attributes. `main.run_forest` is now a thin call into `GridIntegrator.run`; `run1_scalar` remains the reference integrator. Check: on 40 random haloes x 6 cases (both growth models x three feedback types) all 42 output arrays are identical in every digit to the M1 code, and the M0 stored reference outputs are unchanged. `tests/test_components.py` swaps a component. Notes: the UV suppression lives inside the inflow component (as in `gas_evolve.gas_inflow_rate`), so there is no separate `uv.py` yet; the gas/star step is its own module (`gas.py`), which the layout above did not list. `tests/test_run1.py` now patches the ODE step in `components.gas`/`components.chemistry` instead of `main`.
- M2 committed (`a1711a1`).
- M3 done in the working tree (not committed), audit in `docs/M3_M4_AUDIT.md`. New: `config` keys `supernova.delay_distribution` ("single_delay" default | "stellar_lifetimes"), `lifetime_Z`, `n_progenitor_masses`, `wind_model` ("ashvini" default | "momentum" | "energy"), `f_p`, `p_SN_msun_kms`, `eps_fb`, `E_SN_erg`; `stellar_lifetimes.py` (Ekstrom+12 and Georgy+13 grids as 600-byte files in `ashvini/data/`, extracted once from Stochastic_Occupation's tables, nothing downloaded), `imf.py` quantile and sampling functions, `halo.py` (v_esc), `components/delay.LifetimeDistributionDelay` and `explosion_time_fractions`, `components/wind_model.{MomentumDrivenWind,EnergyDrivenWind}`, `presets.stochocc_smooth_control`. Defaults reproduce the M2 output identical in every digit (checked on 42 arrays, both growth models, three feedback types; M0 stored reference outputs unchanged).
- M4 done in the working tree (not committed): `integrators/events.py` (`EventIntegrator`, `main.run_events`), `samplers.py` (`SmoothSampler`, `StochasticSampler`), `trajectory.py`. Chemistry is in; black holes, the cooling threshold and the M5 options are not.
- Acceptance against Stochastic_Occupation's smooth control (stored reference numbers in `tests/reference_outputs/stochocc_smooth.json`, regenerated by `python -m tests.stochocc_comparison`, re-checked live with `pytest --stochocc`), z = 15, 10, 7, 5, dz = 0.005 for the grid:
  - A1 feedback off: stars within 3e-5. B1 chemistry: stellar and gas metallicity within 7e-5 and 8e-5 (grid), 5e-4 (events). E1 IMF delay distribution + momentum wind model: stars within 6e-4, gas within 6e-4 (grid); stars 8e-4, gas 1e-3 (events). The targets were 1-4 % and 2 %.
  - Not to round-off, as the M4 target said: the event integrator differs from Stochastic_Occupation at 1e-3 or below. The sources are small: Ashvini's tabulated Hubble time against an exact one, Ashvini's constants (`Msun_g`, `Gyr_s`) against Stochastic_Occupation's, the substep cap (min of 0.1 t_ff and 1 Myr against 0.1 t_dyn and 1 Myr). I did not remove them, since they are properties of the two codes' conventions, not errors.
  - Stochastic ensembles (200 per code, a minihalo growing 4e5 -> 2e7 Msun where explosions reach the gas limit; energy and momentum wind models): the mean stars formed and gas removed agree within four standard errors, the scatter within 35 %.
- Things found on the way:
  - The Ashvini mass loading equals the momentum wind model with f_p = 0.7646 only at zero stellar metallicity. The audit's 0.7625 was computed when the metallicity step was centred on 0.1 (a mass fraction); with the Z_sun fix the loading at Z_star = 1e-3 is 2.7 % lower than the momentum wind model. The momentum and energy wind models have no metallicity dependence, so a run with the IMF delay distribution and momentum wind model is Ashvini's wind at low metallicity only.
  - Stochastic_Occupation does not define the momentum wind model together with chemistry (an assertion in its engine). The event integrator allows it. The momentum reference numbers are therefore without chemistry.
  - The grid integrator's UV suppression is the rate-dependent form of `reionization.py`; the event integrator's is the equilibrium form (as in Stochastic_Occupation). With UV on, grid and event runs differ for that reason as well.
  - A first version of the delay distribution allocated 3000 x 3000 x 400 numbers (29 GB); it now uses the sorted lifetimes and their cumulative sum.
- Test time: the full suite went from 20 s to about 105 s (the dz = 0.005 grid runs, and 400 stochastic runs per ensemble test).
- Next: M5 (cooling threshold, inflow lag, recycling, stored explosion energy, reachable gas fraction, M_eps_cut, wind travel time, one commit each), then M6.
- M3, M4 committed (`85b26fb`, `84010cc`).
- M5 done in the working tree (not committed): the seven optional event-integrator features, as Config section `events` (`run_params.EventParams`), each off by default: cooling threshold for star formation, inflow lag, recycling, stored explosion energy, reachable gas fraction, change of coupling above a halo mass (`coupling_change_mass`), wind travel time. The M4 stored reference numbers are unchanged, so with the options off the time integrator is the M4 time integrator. Each option has a case on a minihalo growing 4e5 -> 2e7 Msun (z = 20 -> 10) stored from Stochastic_Occupation's smooth control; agreement with the stored reference numbers in stars formed / gas removed: plain minihalo 0.22 % / 0.02 %, cooling threshold 0.06 % / 0.05 %, inflow lag 0.4 % / 0.17 %, recycling 0.3 % / 0.001 %, stored explosion energy (eps_fb = 1) 0.02 % / 0, change of coupling above a halo mass 0.14 % / 0.04 %, wind travel time 1.1 % / 0.005 %, reachable gas fraction 2.2 % / 0.007 % (see below). Each case is checked to differ visibly from its control, so the agreement is not trivial.
  - Access fraction: a cap on each explosion's removal acts per explosion, so the answer depends on how the expected supernovae are split into events. With K = 200 representative masses as in the stored reference numbers the codes agree to 2.2 % in stars; with the Ashvini default K = 100 they differ by 26 %. The agreement uses K = 200. This applies to Stochastic_Occupation's smooth control as well and is not specific to the port.
  - Recycling, bank, wind travel and lag with chemistry: Stochastic_Occupation does not define these combinations; the event integrator does (metals and dust follow the gas). They are tested for conservation only.
  - Test time is now about 4 minutes.
- M5 committed (seven commits, `b28f284` to `d78dabe`).
- M6 done in the working tree (not committed): `docs/STUDENT_GUIDE.md` (layout, one step, recipes, how to know a change is right, cautions), `notebooks/ashvini_tour.ipynb` (executed; eight sections and six exercises; about 80 s), `docs/RESERVOIRS.md` (which reservoir modules are kept unchanged or superseded, and the rules), README sections, `tests/test_docs_examples.py` (the guide's code blocks run). The reservoir files were not moved: that would change the hash and the imports; the note is the archive.
  - Found while writing the notebook: the growth-rate column of `tests/fixtures/merger_trees_fixture.h5` integrates to about 8 % of the mass history's growth (the fixture was thinned in redshift without rescaling it). The grid integrator uses the rate column and the event integrator the mass history, so on that fixture they see different supplies. On trees from `pymctrees_adapter` (rate = (smooth + merger)/dt) the two agree. The notebook takes the rate from the mass history.
