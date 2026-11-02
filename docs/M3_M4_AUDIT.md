# Pre-flight audit for M3 (delay distribution and wind model) and M4 (event integrator)

Written 2026-11-02 before any code. Sources read: `docs/SUPERSET_DESIGN.md`, Stochastic_Occupation `docs/ashvini_relation.md`, `engine.py`, `feedback.py`, `state.py`, `halo.py`, `imf.py` (WeightedQuadratureSampler), `trees.py` (`mass_at`). Nothing here has been run.

## 1. The question

M3: can the grid integrator reproduce the Stochastic_Occupation smooth control, so that "Ashvini with a different delay and wind model" is the same model as that control? M4: can one engine that follows individual supernovae, inside Ashvini, reproduce (a) the smooth control exactly and (b) the stochastic ensembles statistically?

## 2. What the models assume

Common to both (physical): a one-zone gas reservoir fed by f_b dM_h/dt; star formation at eps_ff M_g / (free-fall time); each massive star explodes a main-sequence lifetime after it forms; each explosion removes gas from the reservoir; metals and dust follow the same explosions.

Where they differ, and what the two codes assume:

| | Grid (Ashvini) | Events (Stochastic_Occupation) |
|---|---|---|
| Explosion time | one delay for every star (delta function) | each progenitor mass has its own lifetime (Georgy+13, Z = 0.002 default) |
| Gas removed | ML x (star formation rate one delay ago), continuous | per explosion; energy: min(M_g, 2 eps E/v_esc²), momentum: min(M_g, f_p p_SN / v_esc) |
| Reservoir limit | none in the ODE; a negative gas mass is set to zero afterwards | each explosion removes at most the gas present |
| Halo growth | rate supplied on the tree grid, constant across a step | ln M interpolated linearly in z between nodes; only growth (max(ΔM, 0)) feeds gas |
| Clock | `utils.time_at_z` (astropy Planck18) | its own `CosmicTime` (CGS integral), t_dyn = 0.2047/H |
| Supernova range | 8-40 Msun (parameter now) | 8-100 Msun |

## 3. Assumptions you specified (from the earlier decisions)

D1: keep 0.141/H; the stochocc star formation law is a preset (eps_ff = 0.02178). D2: all bugs fixed; metallicity in Z_sun. D4: only the `v2` time scheme is ported. D5: one package; time integrators/ and components/. M3 acceptance: stars within 1-4 %, metals within 2 % of the stochocc smooth control (A1, B1, E1 of `ashvini_relation.md`). M4 acceptance: stored smooth-sampler results to round-off; stochastic ensemble statistics within sampling error.

## 4. Assumptions I am inferring

1. The grid's Ashvini mass loading is a momentum wind model: `ashvini_relation.md` E1 found ML / (nu_SN p_SN / V_esc) = 0.7623-0.7636 across 1e6-1e10 Msun and z = 5-15. So `Closure` has three members: Ashvini's loading (unchanged, the default), momentum (f_p, p_SN) and energy (eps_fb, E_SN). The factor is constant to 0.2 %, not exactly, so "momentum with f_p = 0.7625" reproduces Ashvini to that accuracy, not to round-off.
2. The delay distribution on the grid is the lifetime distribution of the SN progenitors: the fraction of explosions from stars formed in step k that occur in step j, from the IMF (Salpeter 8-100 here) and the Georgy+13 lifetimes at Z = 0.002. Applied as a convolution of the star formation history. For piecewise-constant star formation this is exact; the grid remains first order in the step.
3. The single-delay distribution at 15 Myr stays the default and reproduces the current output identical in every digit (the M0 stored reference outputs and the M2 comparison are the test).
4. The event integrator uses the same components as the grid integrator for everything that is not event-specific (star formation law, inflow including the UV suppression, chemistry yields, black hole), and adds only: a queue of explosions, a sampler, the per-explosion wind model with the reservoir limit, and the exact-substep scheme `v2`.
5. For M4 agreement the stochocc halo history (ln M linear in z, max(ΔM, 0)) is an input option of the event integrator; the grid integrator keeps the supplied rate.

## 5. Where existing code is inconsistent with these

* The grid has no reservoir limit, the events do. In a halo that forms stars faster than it can hold gas (the minihaloes, the whole point of the paper) the two diverge by construction. The acceptance cases (1e9 Msun halo, from `ashvini_relation.md`) are not in that regime; the 1-4 % there comes from the delay mismatch only.
* The two clocks and the two free-fall times differ (negligibly in the clock, by the 0.2047/0.141 ratio in the law, which the preset absorbs). I have not compared the two clocks numerically.
* stochocc's smooth control is itself a quadrature (K = 100 representative progenitor masses, released in 1 Myr batches), not the exact IMF integral. The grid delay distribution built from the IMF directly is the K to infinity limit; the difference is a convergence property of the control, checked in stochocc's `control_convergence`.

## 6. Ambiguities that change results materially

A. **Reservoir limit on the grid.** (a) Leave it out of the GridStepper: the smooth-control reproduction holds where gas is plentiful and fails where gas runs out; the limit exists only in the EventStepper, where it is natural. (b) Add a limit to the GridStepper: the removal each step is capped at the gas available. This changes the stored reference Ashvini output in the gas-poor regime and makes the grid depend on the step size more. My recommendation is (a). The consequence is that "grid versus events" comparisons in minihaloes are a statement about the gas limit, not only about discreteness; the audit docs of stochocc already say this.
B. **Which clock and halo-growth rule decide the agreement tests.** Ashvini's astropy clock with the stochocc ln M interpolation as an event-integrator option, tests at 1e-4 relative; or a stochocc-compatible clock inside the agreement tests only so they hold to round-off. I recommend the first, with the second as a test fixture so the "to round-off" statement is checked once, against stochocc's own numbers, without changing Ashvini's clock.
C. **Which samplers to port in M4.** stochocc has five. The paper's control is `WeightedQuadratureSampler`; the stochastic one is `StochasticSampler`/`sample_episode`. `PacketSampler`, `DeterministicQuadratureSampler` and `SwitchOffSampler` are diagnostic variants. I recommend porting the two used in the paper and the K-quantile IMF machinery, and leaving the other three in Stochastic_Occupation.
D. **What M4 contains.** The design puts the cooling threshold, inflow lag, recycling, stored explosion energy, f_access, M_eps_cut and wind travel time in M5. Chemistry (metals, dust, the dust destruction per explosion) is not listed in M5, so I read it as part of M4; the T_vir cooling threshold is in M5, so M4's stored tests run with the cooling threshold off. Say if you want the cooling threshold in M4.
E. **Data.** The lifetime grids come from the Ekstrom+12 and Georgy+13 tables (`scripts/ekstrom_2012/` in stochocc, about the size of a few hundred kB). They would be copied into the Ashvini package as data. Only local copies are used; nothing is downloaded.

## 7. Minimum tests that separate correct from merely running

M3
* Delta delay distribution at 15 Myr, momentum-free default wind model: M0 stored reference outputs and the M2 digit-by-digit comparison unchanged.
* Kernel normalisation: the weights over all later steps sum to 1 for a star formed early enough (every progenitor explodes), and the mean delay equals the IMF-averaged lifetime.
* Mass conservation at round-off with the IMF delay distribution (gas + stars + ejected = supply).
* Wind model limits: energy wind model scaled by nu_SN reproduces `feedback.energy_response_delta_M_ej`; momentum with f_p = 0.7625 matches Ashvini's loading to 0.5 % over 1e6-1e10 Msun, z = 5-15.
* The stochocc acceptance cases A1, B1, E1: stars within 1-4 %, metals within 2 %, with the stochocc numbers in `ashvini_relation.md` as the targets.

M4
* No-feedback and no-star-formation limits against closed forms.
* Mass conservation at round-off for any sampler, with and without the reservoir limit.
* Events processed at their own times: an event inside a step uses v_esc at its time (stochocc `test_engine`).
* Smooth sampler equals stochocc's stored numbers on a synthetic trajectory.
* Stochastic ensemble mean and scatter against the stochocc ensemble at matching seeds and sample size, within the sampling error of the mean.
* With a very large number of progenitors per event the event integrator converges to the GridStepper with the IMF delay distribution (the smooth limit is shared).

## Not in scope

New physics, new figures, a decision on rebuilding the MVM, and the M5 options.
