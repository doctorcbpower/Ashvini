# Pre-flight audit: the mass and redshift dependence of the supernova mass loading

Written 2026-11-03 before any code. Sources read: Kravtsov & Manwadkar (2022, MNRAS 514, 2667), pages 1-14 (not the discussion or appendices); Pan & Kravtsov (2023), pages 1-6; `ashvini/supernovae_feedback.py`, `ashvini/components/wind_model.py`, `ashvini/integrators/grid.py`; `docs/PARAMETER_SPACE_AUDIT.md`. Nothing here has been run.

## 1. The question

Does the way the wind removes gas, as a function of halo mass, stellar mass and redshift, control (a) how steeply M*/M_halo rises with mass, (b) the stars in haloes below 1e10 Msun, and (c) how the relation evolves with redshift? GRUMPY gets its M*-M_halo relation from a loading that depends on stellar mass and has no explicit redshift dependence; Ashvini's depends on halo mass and redshift with fixed exponents. This tests which form the data and the simulations prefer. It is not a fit.

## 2. Correction to what I said earlier

I said the loading cutoff in GRUMPY "is probably what lets M*/M_halo rise" towards high mass, and that Ashvini's lack of rise might come from the fixed mass scaling. The parameter-space run shows that is too strong: the accepted parameter sets, calibrated only at 3e10 Msun and above, already lie between Moster and Behroozi at 1e10 and rise to 1e12 with the existing form. The unresolved discrepancies are below 1e10 Msun (the floor of 1e-5 against 1e-4 to 1e-3 in the extrapolated relations and in Nadler et al.) and the redshift evolution (the z = 5 offset of 3 to 20 times). Those are what these changes can address.

## 3. The models compared

| | Ashvini now | GRUMPY (Kravtsov & Manwadkar 2022, eq. 34) |
|---|---|---|
| Loading | eps_p pi_fid (1e11.5/M_h)^(1/3) (9/(1+z))^(1/2) f(Z) | max(0, 2000, eta_norm M*,10^(-eta_p) - eta_c) |
| Depends on | halo mass, redshift, metallicity (a step at 0.1 Z_sun) | stellar mass only |
| Values | eps_p = 5 default, accepted 7 to 11 | eta_norm = 3.6, eta_p = 0.35 (FIRE fit of Muratov et al. 2015), eta_c = 4.5 |
| High mass | continues to fall as M_h^(-1/3) | zero above M* of about 5e9 Msun |
| Low mass | no maximum | capped at 2000 |

Anchors: an energy-driven wind has eta proportional to M_h^(-2/3), a momentum-driven wind M_h^(-1/3) (Ashvini's is momentum-like); the energy wind model is already in the code. In GRUMPY the amplitude is eta_norm and the slope eta_p, and they show (their Fig. 3) that the M*-M_h relation moves with both.

## 4. What you specified

Look at the changes I listed: a stellar-mass-dependent loading with a cutoff, and a free mass scaling. Nothing else.

## 5. What I inferred

- That the redshift dependence should be tested too, because GRUMPY's loading has none and finds almost no evolution of M*-M_h below 1e10 Msun (their Fig. 6), while Ashvini's accepted sets evolve by 10 to 30 times at fixed mass ratio.
- That the three changes are separate steps, each with the default unchanged.

## 6. Where this may mislead

- **Stellar mass is not a mechanism.** The Muratov et al. loading is measured near 0.25 R_vir in FIRE and fitted against M*; it stands in for the depth of the potential and the surface density of star formation. In a one-zone reservoir with M* as the argument, a halo that has formed stars early is loaded differently from one of the same mass that forms them late.
- **Self-regulation through M*.** The loading falls as M* grows, so the wind weakens as the galaxy builds, which steepens the relation by construction. The explicit step in the grid integrator may need the M* from the start of the step; I would check sensitivity to the step.
- **Where the gas goes.** In GRUMPY the outflow leaves and does not return, and 44 percent of the mass of each new star is returned instantly (R = 0.44). I have not checked whether Ashvini's stellar mass counts what is formed or what stays locked, which is a 0.35 dex difference that matters against a 0.2 dex criterion.
- **The cutoff is empirical.** Muratov et al. found no measurable outflow at Milky Way mass; the cutoff is not derived. Above it Ashvini has only the AGN wind, which the parameter run found inert.
- **The comparison relation.** GRUMPY's target is a Behroozi-form relation with parameters adjusted to Nadler et al. (2020) below 1e9 Msun (log10 M1 = 11.39, log10 eps = -1.69, alpha = -2, delta = 4.34, gamma = 0.53, read from the figure caption page), and it uses the peak halo mass of satellites, not the z = 0 mass of field haloes.

## 7. Decisions that change the answer

A. **Order.** (1) Make the halo-mass exponent and the redshift exponent parameters (default 1/3 and 1/2; the output is identical at the defaults), then scan them on the stored trees. (2) Add a fourth wind model `stellar_mass` with eta_norm, eta_p, eta_c and the maximum. (3) Compare. Recommendation: do both, in this order, since (1) is cheap and tells us whether (2) is needed.
B. **Cutoff.** Include it as parameters (eta_c, with 0 meaning none), so the effect of the cutoff alone is a one-parameter comparison. Recommendation: this.
C. **Redshift dependence.** Keep it in Ashvini's form and scan its exponent from 0 to 1/2; for the stellar-mass wind model leave it out, as GRUMPY does. Recommendation: this.
D. **Third reference relation.** Add the Nadler-adjusted relation above as a target for 1e8 to 1e10 Msun, labelled as constrained by upper limits for satellites. Recommendation: yes, shown but not used in the acceptance.
E. **Mass definition.** Check which stellar mass Ashvini reports and whether instant recycling is included, and state the comparison with GRUMPY's M* on the same basis, before any scan.
F. **Scale of the runs.** The scans are on the stored trees at about 14 s per sample, which is light enough to run locally. A repeat of the full Saltelli design with the new parameters would go to OzSTAR.

## 8. Minimum checks

1. With the new parameters at their defaults, the output is identical in every digit to the stored runs.
2. With a constant loading (eta_p = 0, eta_c = 0) and no reionisation, the regulated equilibrium stellar mass matches the closed form.
3. The cutoff: with M* above the cutoff the loading is exactly zero; below it the loading equals the formula; the maximum is respected.
4. Mass balance (gas, stars, removed mass) holds to round-off for each new option.
5. The result changes by less than a few percent when the time step is halved, for the stellar-mass wind model.
6. The loading of the stellar-mass model at the Muratov values reproduces eq. 34 at tabulated points.
7. Plots are made from stored outputs, with the criteria stated in advance.

## Not in scope

Preventive feedback (the GRUMPY inflow factor eps_pr), the molecular hydrogen model, and the satellite and peak-mass comparison.

## Results of the first scan (2026-11-03; `scripts/wind_scaling_scan.py`, stored trees, equilibrium UV, efficiency 0.06)

Each form has its amplitude set so the median M* at 3e10, 1e11, 3e11 and 1e12 Msun matches the mean of Moster+13 and Behroozi+13 (tolerance 0.03 dex); one tree set. Per-form output: `scripts/output/shmr_z0/wind_scan/*.json`. The full test suite passes with the defaults unchanged (314 passed, 2 skipped).

- **Shape at 3e10 and above** (scatter of the residuals about their mean): stellar-mass loading with slope 0.5 and the cutoff 0.04 dex; Ashvini loading with mass exponent 2/3 0.05 (redshift exponent 1/2) or 0.07 (0); mass exponent 1 0.09; 1/2 0.10; the current form (1/3) 0.14 to 0.15; stellar-mass loading with the FIRE slope 0.35 0.14 with the cutoff and 0.22 without; slope 0.2 0.28 to 0.32. The cutoff improves the fit at every slope.
- **M*/M_halo at z = 0, 1e10 Msun:** 3.9e-4 (mass exponent 1) to 3.3e-3 (slope 0.2, no cutoff); the current form 1.4e-3. The steeper forms land on Moster+13 (4.6e-4), the shallow ones above Behroozi+13 (1.8e-3).
- **At 1e8 to 1e9 Msun:** 8e-6 to 5.5e-5 for every form; the wind form does not change the floor set by reionisation by more than a factor 4 to 7.
- **Redshift:** at 1e9 to 3e9 Msun and z = 5, M*/M_halo is 5.6e-4 to 3.2e-3 for all forms. Dropping the (1+z) factor of the Ashvini loading lowers z = 5 by about 1.6 times and z = 2 by about 1.3 times: part of the evolution, not all.
- **Cautions.** The slope 0.5 of the best stellar-mass fit was fitted to the same four masses (the FIRE value is 0.35); the energy-driven exponent 2/3 has a physical argument. Ashvini's M* counts all mass formed, GRUMPY's locks 56 percent (about 0.25 dex). The 1e11 to 3e11 Msun bin at z = 5 has too few haloes. The stellar-mass wind is implemented for the grid integrator only.
