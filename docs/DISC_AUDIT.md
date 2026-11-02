# Pre-flight audit: giving the gas a disc size, and a star formation threshold that depends on it

Written 2026-11-02 before any code. Sources read: the DiscModels README and function list (`compute_Rd`, `fR_fitting`, `qth.py`: `Sigma_crit`, `Q_thin`, `epicyclic_kappa`); the earlier discussion. The DiscModels code itself has not been read. Nothing here has been run.

## 1. The question

For a small set of haloes, does giving the gas a disc size, and allowing star formation only where the gas surface density is high enough, change the predicted stellar mass, metallicity and star formation history compared with treating the gas as one bulk reservoir? The point is not to make Ashvini a disc model. It is to find out whether this ingredient matters, and by how much, so that the result can be put in an ideas paper and, if it matters, passed to Arepo Solas.

## 2. Physical assumptions of the proposal

1. The gas that has settled in the galaxy forms an exponential disc, surface density `Sigma(r) = Sigma_0 exp(-r/R_d)`, with `Sigma_0 = M_gas / (2 pi R_d^2)`.
2. The scale length `R_d` follows Mo, Mao and White (1998): in the simplest form `R_d = (j_d / m_d) * lambda * R_vir / sqrt(2)`, where `lambda` is the halo spin parameter, `R_vir` the virial radius, `m_d` the disc mass as a fraction of the halo mass and `j_d` the fraction of the halo's angular momentum it holds. The fuller form multiplies by a correction that depends on the halo concentration and on `m_d` (`fR_fitting` in DiscModels), and includes the pull of the disc on the halo (adiabatic contraction).
3. Gas forms stars only where `Sigma(r)` is above a critical value `Sigma_crit`. For an exponential disc this is every radius inside `R_crit = R_d ln(Sigma_0 / Sigma_crit)`, which holds the fraction `f = 1 - (1 + x) exp(-x)` of the gas, with `x = R_crit / R_d`. The star formation rate is then the present rate applied to that gas only: `SFR = f * (eps_ff / t_ff) * M_gas`.
4. `lambda` is constant in time for a halo (the Mo, Mao and White assumption), so `R_d` grows with `R_vir`.

## 3. What you specified

A disc size for a small subset of haloes; the aim of seeing how it changes the measurements; the main benefit being a star formation threshold tied to something more than the bulk reservoir; Ashvini to keep its own niche (no calibrated or full semi-analytic model).

## 4. What I inferred

- That Mo, Mao and White is the size model, with a simple version first.
- That the threshold is a critical surface density, and that gas outside `R_crit` forms no stars.
- That `lambda` is drawn per halo from the log-normal distribution, or fixed.
- That `m_d` is taken from Ashvini's own gas plus stars, not set by hand.
- That the new code lives in Ashvini as a few lines (the formulas above), with DiscModels used only to check them.

## 5. Where this may be inconsistent

- **Equilibrium.** The size formula assumes a halo in equilibrium with fixed spin. At z > 5 haloes grow fast, so the disc size is an estimate, not a result. The sensitivity test is what tells us whether that matters.
- **Concentration.** The fuller formula needs the halo concentration, which Ashvini does not follow. The simple form does not need it; the difference between the two is part of the uncertainty to measure.
- **Two disc models in one package.** The reservoir model already uses `R_d = lambda R_vir / sqrt(2)` (the simple limit). The new component must use the same definition, so there is one disc size, not two.
- **Overlap with the cooling threshold.** Both stop star formation in small haloes. If the surface-density threshold only switches off the same haloes, it adds nothing; we must compare the two on the same haloes.
- **On/off oscillation.** A threshold on the instantaneous surface density can switch star formation on and off from step to step. Using the star-forming fraction `f` of the gas (above) is smooth and avoids this.

## 6. Choices that change the results materially

A. **Threshold type.** (a) A fixed critical surface density (a few Msun/pc^2, as in the Kennicutt relation). (b) The gas stability threshold of `qth.py`, which depends on the rotation speed through the epicyclic frequency and on the gas sound speed. (c) Molecular fraction from surface density and metallicity (as in GRUMPY). (a) has one parameter, (b) needs a rotation curve and a sound speed, (c) needs a clumping factor and a calibration at low metallicity. Recommendation: (a), then compare with (b) on the same haloes.
B. **What happens below the threshold.** (a) No star formation outside `R_crit` (the smooth star-forming fraction `f`); (b) a reduced efficiency outside. Recommendation: (a).
C. **Spin.** (a) Fixed at 0.035; (b) drawn per halo, which adds scatter and an ensemble. Recommendation: (a) for the test, (b) as a second run to see how much of the spread it adds.
D. **Where the formulas live.** (a) Re-implemented in Ashvini (a few lines) and checked against DiscModels in a test; (b) import DiscModels. Recommendation: (a): Ashvini stays self-contained and transparent, and the check keeps the two in agreement.
E. **Which haloes.** The subset: how many and which masses. Recommendation: five haloes at each of 1e8, 1e9, 1e10 and 1e11 Msun from the test trees (or live trees), at z = 5 as the end point.

## 7. Minimum tests that separate a correct implementation from one that only runs

1. With the option off, the output is identical in every digit to the current model.
2. With `Sigma_crit` set to zero (every radius forms stars), `f = 1` and the output equals the current model.
3. The star-forming fraction `f` equals the numerical integral of the exponential profile inside `R_crit`.
4. The simple disc size equals the reservoir model's definition, and the fuller form equals DiscModels `compute_Rd` for a CDM halo to within a few per cent.
5. As `lambda` becomes large, `R_d` grows, `Sigma_0` falls, and star formation stops.
6. Mass balance still holds to round-off (gas, stars and removed mass add up to the supply).
7. The result does not depend on the step size beyond what the current model does (no on/off oscillation).
8. The sensitivity table for the subset, run with and without the cooling threshold, shows what the disc threshold adds.

## Not in scope

The disc's effect on the dark matter halo, rotation curves of non-CDM models, a molecular hydrogen model, disc morphology, and stars as a separate disc component.
