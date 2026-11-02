# Ashvini as the main codebase: superset design (draft for agreement, 2026-11-02)

Status: proposal only. No code has been changed. Nothing below is implemented until the decisions in section 6 are agreed.

Direction (decided by the user): Ashvini is the main codebase. Stochastic_Occupation (`stochocc`) becomes a set of components and an event stepper inside Ashvini, and its analysis scripts later call Ashvini. `stochocc` is then frozen as the record of the published results.

## 1. Pre-flight model audit

| # | Item | Today in Ashvini | Today in stochocc | Proposal |
|---|---|---|---|---|
| 1 | Star formation law | eps_ff M_g / (0.141 / H) | eps_ff M_g / t_dyn, t_dyn = 0.2047 / H | One `SFLaw` component. Same law up to eps_ff x 0.2047/0.141 = 1.452. Keep both as named presets. |
| 2 | Feedback timing | one lag, 15 Myr, fixed | per-progenitor lifetimes (Georgy+13) | `DelayKernel`: delta function (Ashvini) or IMF-weighted lifetime distribution (stochocc smooth limit). |
| 3 | Feedback amplitude | momentum mass loading, eta(M_h, z, Z*) | energy closure min(M_g, 2 eps E / v_esc^2), or momentum | `Closure` component with three members: `momentum_ml`, `momentum_event`, `energy`. Parity note shows `momentum_ml` equals `momentum_event` with f_p = 0.7625 for 8-100 Msun. |
| 4 | Discreteness | none (smooth) | Poisson-sampled supernovae, smooth control | `Sampler` option on the event stepper: `smooth`, `stochastic`. |
| 5 | Star formation gate | none | T_vir >= 1e4 K, or H2 mass gate, plus post-reionisation switch | `Gate` component, default `always_on` to keep Ashvini's results unchanged. |
| 6 | UV suppression | rate-dependent Okamoto form (own audit flags it) | equilibrium form | `Inflow` takes a `UVSuppression` object. Both forms are available and named; equilibrium form is the shared base in `reionization_physics.py`. |
| 7 | Chemistry | gas and stellar metals, dust, with known issues | metals and dust per event | One `Chemistry` component after the known dust and Z_IGM issues are fixed. |
| 8 | Black holes, seeding, AGN | several growth laws, three seeding channels | none | Unchanged, each behind a `BlackHole` interface. |
| 9 | Reservoirs | one-zone (general chain) and MVM galaxy plus nucleus | one-zone | `Reservoirs` interface with `one_zone` and `mvm`. MVM stays frozen and bit-identical (decision D3). |
| 10 | Time stepping | vectorised over haloes on the tree grid | one halo at a time, events at their own times | Two steppers on one `State` and one `Config`. |
| 11 | Extra stochocc features | none | inflow lag, recycling, energy bank, f_access, M_eps_cut, wind travel time | Ported as optional components, each off by default, after the base parity tests pass. |

Not in scope: new physics, new figures, new scans, manuscript changes.

## 2. Bugs to fix first (reproduced, 2026-11-02)

1. `supernova.delay_time: 0` with `type: delayed` gives exactly the no-feedback result. The lookback index equals the current step, and the SFR of the current step is not yet computed. It should either mean instantaneous feedback or be rejected.
2. `type: "no"` still produces dust and still removes dust with the wind; `type: "instantaneous"` still drives dust from the delayed SFR. Metals and dust must follow the same feedback switch.
3. Dust supernova rate is Gamma / dt_step (grid-dependent), with no M_d / M_g factor.
4. Stellar mass uses gas mass at the start of the step (first-order lag, 1.1% between dz = 0.05 and 0.005).
5. `Z_IGM = 1e-3` is documented in Z_sun units but used as a mass fraction.

Each fix changes results. They need agreement, a regression test that pins the old behaviour for comparison, and a note in the change log. The frozen MVM production run does not use the general chain, so it is unaffected.

### Status of the bug fixes (2026-11-02)

All five are fixed in the working tree (uncommitted), with tests in `tests/test_chemistry_and_feedback_switch.py` and the pinned fixture baseline in `tests/test_run1.py` recaptured (old values recorded there). Also fixed on the way: `supernova.type: "instantaneous"` was switched off until `delay_time` had elapsed. Metallicity units: `Z_IGM_Zsun`, seeding `Z_gas_max_Zsun`, and the `metallicity_function` midpoint (0.1 Z_sun) are all in Z_sun; old keys are rejected with an explanatory error. New: `supernova.imf`, `m_sn_min`, `m_sn_max`. The frozen MVM is untouched (its mass loading uses Z = 0, for which the new function is identical).

Effect on the 1e10 Msun fixture halo at the last step: gas +1.4%, stars +0.5%, gas metals -26%, stellar metals -45%, dust -13%, BH mass unchanged. Analytic 1e7-1e9 Msun halos (no BH, no UV): stars +1.2%, gas metals -55% to -78%, dust between +26% and -7%.

## 3. Architecture

- `Config`: one dataclass tree, passed explicitly. No module-level reads of `PARAMS`; `set_params` is replaced by building a `Config`. `run_params.yaml` becomes one way to build a `Config`.
- `State`: a dataclass of arrays (gas, stars, metals, dust, BH, plus the event queue for the event stepper).
- Components, each a small class with one documented method and a one-paragraph physics docstring:
  `SFLaw`, `Inflow`, `UVSuppression`, `Gate`, `DelayKernel`, `Closure`, `Chemistry`, `BlackHole`, `Seeding`, `AGNFeedback`.
- Registry: `Config` names components as strings (`closure: "energy"`). A student adds a component by writing a class and registering it in one place.
- Steppers: `GridStepper` (today's `run_forest`, vectorised) and `EventStepper` (ported from `stochocc.engine`). Both call the same components.
- Cosmology and constants: one place (`utils.py`, `constants.py`). The Planck15 use in `reionization.py` goes away or is made an explicit option.

## 4. Student-facing pieces

- `docs/GUIDE.md`: which of the four existing model families to use, and which to ignore (`reservoir_stock_premvm.py`, `paper_reservoir.py` integrator).
- A notebook that swaps one component at a time and plots the effect.
- Each component page lists its equation, parameters, defaults and the paper it comes from.

## 5. Order of work and tests

1. Characterisation tests: pin current outputs of the general chain and the MVM before any change.
2. Fix bugs in section 2, one commit each, each with a regression test.
3. Introduce `Config` and `State`, with no behaviour change.
4. Extract components from `run_forest`; the characterisation tests must still pass.
5. Add the `DelayKernel` and `Closure` options and reproduce the stochocc smooth control through `GridStepper` (parity target: the existing A1, B1, E1 results in `ashvini_relation.md`).
6. Port `EventStepper`; target: pinned stochocc results from its test suite.
7. Port the optional stochocc features one at a time.
8. Guide and notebook.

## 6. Decisions needed

- D1. DECIDED (user, 2026-11-02): default free-fall time stays 0.141/H; the Stochastic_Occupation law is a named preset (eps_ff = 0.02178). The supernova upper mass limit becomes a parameter (default 40 Msun in the general chain, 100 in the Stochastic_Occupation preset).
- D2. DECIDED (user, a main co-author, 2026-11-02): fix all five bugs. Bug 5 (user, 2026-11-02): intended value is 1e-3 Z_sun, probably muddled to a mass fraction in MBP26. Policy: every metallicity in the config carries its unit in its name (`..._Zsun` or `..._massfrac`), one `Z_SUN = 0.015` lives in `constants.py`, and the resolved values in mass fraction are printed in every run header. Same audit to be applied to the seeding `Z_gas_max` thresholds and the `metallicity_function` midpoint (m = 0.1), pending confirmation.
- D3. DECIDED (user, 2026-11-02): the MVM stays frozen next to the framework. `reservoir_stock.py` (and the superseded `reservoir_stock_premvm.py`, `paper_reservoir.py` integrator) are not modified; the provenance hash stays valid. They may import shared constants, cosmology and reionisation physics. Rebuilding the MVM from components is a separate, later step with a round-off parity test against the frozen copy.
- D4. DECIDED (user, 2026-11-02): the Stochastic_Occupation `legacy` time scheme is dropped; only `v2` (exact linear-ODE substep, SFR-weighted event time stamp) is ported. Any early result that used `legacy` stays reproducible from the frozen `stochocc` repo.
- D5. DECIDED (user, 2026-11-02): one package, with an event-stepper subpackage (layout in section 7). Stochastic_Occupation stays as a thin analysis repo that imports Ashvini, and is frozen as the record of its published results once Ashvini reproduces them.

## 7. Concrete layout and milestones (proposal, for agreement before any code)

```
ashvini/
  config.py            Config dataclass tree + loaders (yaml, dict); replaces PARAMS/configure.set_params
  state.py             State dataclass (arrays), event queue for the event stepper
  constants.py  utils.py  imf.py  reionization_physics.py     (shared, already exist)
  components/
    star_formation.py  SFLaw: rate(state, t) -> dM*/dt per halo; presets "ashvini" (0.141/H), "dyn200" (0.2047/H)
    inflow.py          Inflow: mass supply per step; wraps a UVSuppression
    uv.py              UVSuppression: "none", "okamoto_eq", "ashvini_rate"
    gate.py            Gate: "always", "tvir", "mass", plus post-reionisation switch
    delay.py           DelayKernel: "fixed" (Ashvini), "lifetimes" (IMF-weighted), "none"
    closure.py         Closure: "momentum_ml", "momentum_event", "energy"; clip and f_access options
    chemistry.py       Chemistry: metals and dust (the corrected equations)
    blackholes.py      BlackHole: "pznk11", "slimdisk"; Seeding; AGNFeedback
    optional.py        inflow lag, recycling, energy bank, wind travel time (off by default)
  steppers/
    grid.py            GridStepper: today's run_forest, vectorised over haloes
    events.py          EventStepper: ported from stochocc.engine (v2 only)
  reservoirs/          (frozen, unchanged) reservoir_stock.py, reservoir_stock_premvm.py, paper_reservoir.py
  main.py              thin: run(), run_forest() and run1_scalar() kept as compatibility wrappers
```

Interface rule: each component is a small class with one documented method, a one-paragraph physics docstring with the equation and its source, and a `name` registered in one place (`components/__init__.py`). A student adds a component by writing the class and adding one line to the registry. No component reads a module-level global; everything comes from the `Config` it is given. Validation uses exceptions with messages, not `assert`.

Milestones, each ending with the full test suite green and the pinned outputs unchanged unless stated:

- M0. Characterisation: pin general-chain outputs (several halos, all three feedback types, BH on and off, UV on and off) and the MVM output hash, as tests.
- M1. `Config` and `State` introduced; `PARAMS`-based modules read from a `Config` passed explicitly; `configure.set_params` kept as a thin wrapper. No behaviour change.
- M2. Extract components from `run_forest` (SFLaw, Inflow, UV, Closure, Chemistry, BlackHole) behind the GridStepper. No behaviour change; M0 pins hold to round-off.
- M3. `DelayKernel` and `Closure` options; reproduce the Stochastic_Occupation smooth control through the GridStepper. Acceptance: A1, B1, E1 of `ashvini_relation.md` (stars within 1-4%, metals within 2%).
- M4. EventStepper (v2 only). Acceptance: stochocc's pinned results for the smooth sampler to round-off, and ensemble statistics for the stochastic sampler within sampling error.
- M5. Optional stochocc features, one commit each: gate, inflow lag, recycling, energy bank, f_access, M_eps_cut, wind travel time.
- M6. Student guide and notebook; archive the superseded paths in a `reservoirs/` note.

Not done in this plan: new physics, new figures, changes to the MVM, a decision on rebuilding the MVM from components.
