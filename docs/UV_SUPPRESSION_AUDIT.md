# Audit of the UV-suppression term in the general model

Date: 2026-09-20. Scripts: `scripts/uv_suppression_audit.py` (variants on identical trees, output
`scripts/output/uv_suppression_audit_cdm.json`). Trees: Zhang-Hui, CDM, `z0 = 0`, `z_max = 30`,
`M_res = 1e-4 M0`, 100 halos per mass (60 for the dz test), `dz = 0.005` unless stated.

## The term

`reionization.uv_suppression` multiplies the cosmological baryon inflow by

```
S = max(0, s(mu, omega) [ (1 + X) - 2 eps(z) M X (1 + z) H / Mdot ]),   mu = M / M_c(z),   z <= 10
```

Differentiating an equilibrium gas mass `M_gas,eq = f_b s(M/M_c) M` gives the same form with
`eps = -d ln M_c / d z` and no factor 2, so the correction term is the rate at which the equilibrium
gas fraction falls because `M_c` grows. Inflow cannot go negative, so when that fall outpaces halo
growth the term simply stops accretion. Gas already in the halo is not removed.

## Findings

1. **Implementation of `eps(z)`.** The shipped `epsilon(z)` times 2 equals `-d ln M_c / d z` for
   `z <= 5` (checked by finite difference) but not above: the ratio is 0.99 at z = 5, 0.85 at z = 6,
   0.09 at z = 7 and 0 for z >= 8 (true value 6.9 at z = 7 and 960 at z = 10). Using the exact
   `eps` (`true_eps` in the script) changes the median `M_star(z=0)` by 5 to 25 per cent at
   `M_halo <= 3e10` and not at all above, so it is a real discrepancy with a small effect here.
   The `eps` formula in `MODELS.md` also does not match the code.
2. **Strength.** Fraction of the cosmological baryon accretion transmitted, accretion-weighted, median over halos:

   | M_halo | current, z<=5 | mass-only, z<=5 |
   |---|---|---|
   | 1e9  | 0.000 | 0.005 |
   | 1e10 | 0.000 | 0.59 |
   | 1e11 | 0.47  | 0.99 |
   | 1e12 | 0.87  | 1.00 |

   With the shipped form, halos below ~1e10 Msun receive essentially no accretion after z = 5 and
   halos at 1e11 receive about half. That is a shut-off, not a filtering-mass suppression, and it is
   what makes the low-mass SHMR (and the WDM/FDM zeros) so steep.
3. **Time-step sensitivity.** Median `M_star(z=0)` changes with the tree `dz` (0.02, 0.01, 0.005) by a
   factor 1.5 (1e10) and 2.3 (1e11) with the shipped term and by 1.03 to 1.15 with the mass-only form or
   with the term off. The ratio term divides by a per-step `Mdot` estimated from the tree, which is what
   makes it depend on the step. GRUMPY smoothing of `Mdot` does not remove this.
4. **Against the literature curves** (both extrapolated below ~1e10.5 Msun): ratio of median
   `M_star` to Behroozi+2013 is 0.01 (1e10), 0.11 (3e10), 0.58 (1e11) with the shipped term and
   1.02, 1.50, 1.19 with the mass-only form.

## What is and is not established

Established: the numbers above for these trees and this model. Not established: which form is right
physically. The mass-only form is a filtering-mass suppression of the equilibrium baryon fraction
(Okamoto et al. 2008) and is what the MVM uses; the shipped form additionally imposes that gas cannot
accrete while the equilibrium fraction is falling, which is a modelling choice this audit does not test
against simulations. The DM-model differences at low mass in the SHMR figure are controlled by this term.
