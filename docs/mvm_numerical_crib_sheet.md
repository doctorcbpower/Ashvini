# MVM numerical crib sheet

Purpose: the numerical story of the frozen Minimal Viable Model, written once, organised by claim, so that the
paper can be drafted against it. It is not text for the paper. Nothing here has been put into the manuscript.

## 0. Frozen state and conventions

* Model: `ashvini/reservoir_stock.py`, sha256 `31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63`
  (recorded in the production JSON). The exploratory model is frozen as `reservoir_stock_premvm.py`.
* Fiducial (frozen 2026-09-20): R_nuc = 100 pc, sigma_j = 0.5, eta_acc = 0.005, eps_sf = 0.015, f_mom = 1 (v_out = sigma),
  eta_SN x1, n_rd = 10, NFW c = 4 (galaxy scale only), radiative efficiency 0.1, f_b = 0.156, M_hot = 4e11 Msun,
  lambda = 0.035, z_seed = 25 -> z = 5, f_BH = 0.5. Strict Eddington (exact step). Momentum-driven AGN outflow,
  removed in proportion from the two gas reservoirs together with the stellar outflow. No recycling.
* Production: 13 halo masses, 3e10 to 3e13 Msun at 0.25 dex; 240 Zhang & Hui trees per mass; 801 uniform-time steps;
  dz = 0.05; M_res = 1e4 Msun. Data: `scripts/paper_figures/output/mvm_production_results.json`;
  numbers below marked (P) are printed by `scripts/paper_figures/mvm_crib_numbers.py`.
* Numbers marked (D) come from `scripts/paper_figures/diagnostics/logs/` (see the README there); they use their own,
  smaller tree ensembles, and the tree sampler is not seed-reproducible.
* "Median [16,84]" = median and 16th and 84th percentiles over trees. "16-84 width" = log10(p84/p16).
* "Accessibility test" = the labelled experiment R_nuc = 250 pc, sigma_j = 0.5. It is a sensitivity experiment on the same
  model, not a second physical model and not a candidate fiducial.
* Agreed framing (user, 2026-09-20): M_seed,crit is the initial BH mass required to reach a specified M_BH/M_star,tot
  at z = 5, given a particular galaxy baryon cycle and nuclear-accessibility model. It is not a measure of the Eddington
  growth time available to a seed.

## 1. Claims, quantities, results, uncertainties, and what not to say

### C1. Critical boundary

| Quantity | Result (P) | Numerical uncertainty / caveat | Do not say |
|---|---|---|---|
| M_seed,crit at M_halo(z=5) = 3e10 | 9.45e7 [7.61e7, 1.06e8] Msun | see C13 | that it is a predicted seed mass |
| at 3e11 | 1.39e9 [1.03e9, 1.61e9] Msun | " | " |
| at 3e13 | 4.85e9 [4.54e9, 5.17e9] Msun | " (dz and M_res untested here) | " |
| over 13 masses | rises monotonically 9.45e7 -> 4.85e9 Msun | | that the boundary is a power law |
| M_seed,crit / M_halo(z=5) | 3.2e-3 (3e10), 4.6e-3 (3e11), 1.6e-4 (3e13) | | |
| bracket failures; F(seed) crossing zero | 0 failures at all masses, max abs F = 2e-11; F crosses zero once in every tree tested (40 trees per mass at 3e10, 3e11, 3e13; 20 trees in every accessibility corner tested) (D) | monotonicity checked on 25-point seed grids, not proved | that the root is unique in general |

### C2. Tree-to-tree scatter

| Quantity | Result (P) | Caveat | Do not say |
|---|---|---|---|
| 16-84 width in M_seed,crit | 0.142 dex (3e10), 0.195 (3e11), 0.057 (3e13); range 0.057 to 0.195 over the 13 masses | half-widths 0.071, 0.097, 0.029 dex; width narrows above about 1e12 | that it can be compared with the old paper's 0.12 to 0.31 dex (different model and trees; definition there not re-checked) |

### C3. Mass dependence and the kink

| Quantity | Result (P) | Caveat | Do not say |
|---|---|---|---|
| local slope dlogM_seed,crit / dlogM_halo | 1.17 (3e10 -> 3e11), 0.67 (3e11 -> 9.5e11), 0.24 (9.5e11 -> 3e12), 0.09 (3e12 -> 3e13); whole range 0.57 | | that the slopes are universal |
| M_star,tot / M_halo(z=5) | 6.5e-3 (3e10), 9.7e-3 (3e11), 3.5e-4 (3e13); peak 9.7e-3 at 3e11 | | |
| mass accreted / (f_b M_halo(z=5)) | 1.16 (3e10), 1.07 (3e11), 0.84 (5.3e11), 0.17 (3e12), 0.017 (3e13) | includes the UV and cold/hot suppression | |
| origin of the kink | the imposed cold/hot inflow step at M_hot = 4e11 Msun, acting through M_star,tot | this is where the mass dependence comes from | that the kink is an emergent black-hole result; that the low 3e13 stellar mass is a physical prediction (it is the imposed zeta_ch cutting inflow) |

### C4. M_seed,crit is essentially the target: M_seed,crit = f_BH M_star,tot / G_BH

| Quantity | Result (P) | Caveat | Do not say |
|---|---|---|---|
| M_seed,crit / (f_BH M_star,tot) = 1/G_BH | 0.973 [0.971, 0.975] (3e10), 0.956 [0.953, 0.960] (3e11), 0.929 [0.922, 0.935] (3e13); falls monotonically 0.973 -> 0.929 | G_BH - 1 is resolution-dependent (C14); it shrinks as dt is refined | that the BH "grows by 3 to 8 per cent" as a converged physical quantity |
| G_BH = M_BH(z=5)/M_seed at the critical seed | 1.028 (3e10), 1.046 (3e11), 1.077 (3e13) at 801 steps | 3e10: 1.10, 1.05, 1.03, 1.02 at 201, 401, 801, 1601 steps (D) | that G_BH is a physical result; say "at most a few per cent at 801 steps" |
| f_BH M_star,tot (median) | 9.7e7, 1.46e9, 5.2e9 Msun | | |
| decomposition | M_seed,crit = f_BH M_star,tot / G_BH; M_star,tot from the baryon cycle (assembly, angular momentum, SF, feedback); G_BH from nuclear accessibility (sigma_j, R_nuc, eta_acc) | | "two things and only two set the answer" |

### C5. Regime at the fiducial: supply-limited

| Quantity | Result (D: 60 trees, 801 steps) | Caveat | Do not say |
|---|---|---|---|
| Eddington-limited steps at the critical seed | 0.0% of growth steps and 0% of the mass gained, at 3e10, 3e11, 3e13 | fiducial only | that the Eddington limit is irrelevant |
| Eddington e-folds used, of about 20.8 available | 0.028, 0.044, 0.074 (mean accretion about 0.1, 0.2, 0.4 per cent of Eddington) | | that the seed is "using the Eddington time" |
| fixed seeds 1e2 to 1e7 | Eddington-limited in about 0% of growth steps (0.4% of steps for 1e3 at 3e11, 0% of the mass) | | |
| peak Mdot_acc / Mdot_Edd (20 trees) | critical seed 0.04 (3e10), 0.07 (3e11); seed 1e3: 0.09 and 0.05, worst tree 0.4 | a maximum over time | |

### C6. Inverse question: M_BH/M_star,tot achieved by fixed seeds (fiducial)

Median [16,84] at z = 5 (P). Seeds 1e2, 1e3 (runaway-collision scale), 2e5 (direct-collapse scale) are the paper's channel
values; 1e7 is a benchmark seed mass, not a formation channel.

| Seed [Msun] | 3e10 | 3e11 | 3e13 |
|---|---|---|---|
| 1e2 | 4.0e-7 [3.6e-7, 5.1e-7] | 2.8e-8 [2.4e-8, 3.7e-8] | 7.0e-9 [6.4e-9, 7.8e-9] |
| 1e3 | 4.1e-6 [3.6e-6, 5.2e-6] | 2.8e-7 [2.4e-7, 3.8e-7] | 7.0e-8 [6.4e-8, 7.8e-8] |
| 2e5 | 8.2e-4 [7.3e-4, 1.0e-3] | 5.7e-5 [4.8e-5, 7.6e-5] | 1.5e-5 [1.3e-5, 1.7e-5] |
| 1e7 | 4.2e-2 [3.8e-2, 5.4e-2] | 2.85e-3 [2.4e-3, 3.8e-3] | 7.4e-4 [6.7e-4, 8.1e-4] |

Median growth G: 1.000 (1e2), 1.009 / 1.006 / 1.000 (1e3), 1.02 to 1.04 (2e5), 1.02 to 1.05 (1e7).
Caveat: at the fiducial the ratio is approximately M_seed/M_star,tot. Do not say: that the model predicts light seeds
cannot grow (see C7 and C8), that 1e7 is a named seed channel, or that these ratios are converged to better than the
few per cent of C13.

### C7. Accessibility test (R_nuc = 250 pc, sigma_j = 0.5): what changes and what does not

| Quantity | Result (P) | Caveat | Do not say |
|---|---|---|---|
| seed 1e3: median G | 13 (3e10), 272 (3e11), 1.6e3 (5.3e12); median 1 at 1.7e13 and 3e13 | | that R_nuc = 250 pc is a second physical model |
| seed 1e3: median ratio M_BH/M_star,tot | 5.6e-5 (3e10), 8.2e-5 (3e11), peak 1.6e-4 (5.3e12) | | that rapid growth reaches the target |
| highest 84th-percentile ratio reached by a 1e3 seed | 5.9e-4, i.e. 2.93 dex below f_BH = 0.5 | | |
| seed 1e3: fraction of trees with G > 2 | 92% (3e10), 94% (3e11), 76.7% (1.7e12), 52% (9.5e12), 41.7% (1.7e13), 40.8% (3e13); fiducial 0 to 0.8% | | that the accessibility prescription necessarily produces bimodality |
| seed 1e3 at 3e13: median vs upper percentile | median ratio 7.7e-8, 84th percentile 5.6e-4; 84th percentile G about 7,500 | | that the high-mass median drop in Figure 2 is a decline in individual BH growth: it is the median passing through a population that is becoming bimodal (the band carries the information) |
| seeds 2e5 and 1e7 | barely change (1e7: 4.4e-2 vs 4.2e-2 at 3e10) | | |
| tested link with accreted baryons | the growing fraction falls as the accreted baryon fraction falls (5.4% at 9.5e12, 1.7% at 3e13) | link not tested directly | that the cause is established |

### C8. Light-seed Eddington transition (sigma_j, R_nuc)

(D: 20 trees, 801 steps, seed 1e3.) Peak Mdot_acc/Mdot_Edd (median over trees) / growth factor G.

3e10:

| sigma_j \ R_nuc [pc] | 100 | 125 | 150 | 200 | 250 | 300 |
|---|---|---|---|---|---|---|
| 0.5 | 0.09 / 1.01 | 0.4 / 1.07 | 1.2 / 1.19 | 36 / 3.64 | 70 / 15.6 | 69 / 35.6 |
| 0.6 | 0.16 / 1.03 | 0.76 / 1.13 | 3 / 1.32 | 91 / 15.1 | 81 / 43.9 | 100 / 120 |
| 0.7 | 0.42 / 1.07 | 3.1 / 1.28 | 65 / 7.92 | 95 / 53.9 | 130 / 352 | 120 / 764 |
| 0.8 | 15 / 1.77 | 120 / 20.4 | 110 / 59.9 | 190 / 781 | 120 / 1170 | 98 / 2560 |
| 1.0 | 230 / 1480 | 200 / 2540 | 180 / 3320 | 140 / 3850 | 110 / 4070 | 90 / 4300 |

3e11:

| sigma_j \ R_nuc [pc] | 100 | 125 | 150 | 200 | 250 | 300 |
|---|---|---|---|---|---|---|
| 0.5 | 0.049 / 1.01 | 1.7 / 1.22 | 11 / 2.03 | 80 / 35.7 | 400 / 238 | 370 / 364 |
| 0.6 | 0.39 / 1.05 | 10 / 1.87 | 46 / 10.7 | 430 / 320 | 370 / 718 | 420 / 2160 |
| 0.7 | 5.1 / 1.45 | 48 / 10.1 | 430 / 304 | 400 / 1240 | 440 / 3130 | 390 / 7910 |
| 0.8 | 190 / 120 | 440 / 953 | 430 / 1800 | 470 / 5700 | 410 / 13600 | 370 / 20400 |
| 1.0 | 520 / 11600 | 510 / 22400 | 480 / 33300 | 420 / 60300 | 380 / 76800 | 340 / 86800 |

Reading (cells are 25 to 50 pc in R_nuc and 0.1 in sigma_j wide):

* Peak Mdot_acc/Mdot_Edd crosses 1 between R_nuc = 125 and 150 pc at sigma_j = 0.5 (3e10) and between 100 and 125 pc (3e11); at R_nuc = 100 pc
  it crosses between sigma_j = 0.7 and 0.8 (3e10) and between 0.6 and 0.7 (3e11).
* Substantial growth (G of 10 or more) needs R_nuc of about 200 to 250 pc (3e10) or 150 to 200 pc (3e11) at sigma_j = 0.5, or sigma_j of 0.8 to 1.0 (3e10)
  or 0.7 to 0.8 (3e11) at R_nuc = 100 pc.
* The threshold is sharp: over one grid cell the peak rises by one to two orders of magnitude while G lags the peak.

Caveats: light seed only; 20 trees; locations approximate and mass-dependent. Do not say: "a light seed is X times in R_nuc from an Eddington
regime" as a precise number; that the fiducial was chosen far from, or near, the transition; that the transition
affects the critical seed (it does not: see C9).

### C9. The critical seed stays supply-limited across the accessibility grid; eta_acc

(D: 20 trees, 801 steps.) At the critical seed, over sigma_j in {0.5, 1.0, 1.5} and R_nuc in {100, 300, 1000} pc, the peak
Mdot_acc/Mdot_Edd is at most 0.48 and capped steps are 0.0% in every cell (3e10 and 3e11); M_seed,crit falls to 1.3e7
(3e10) and 1.5e8 (3e11) at sigma_j = 1.5, R_nuc = 100 pc, with critical-seed growth G at most 1.8.

| eta_acc series (critical seed) | 3e10 | 3e11 |
|---|---|---|
| 0.005, sigma_j 0.5, R_nuc 100 (fiducial) | peak 0.04, capped 0%, G 1.03 | peak 0.07, capped 0%, G 1.05 |
| 0.05, sigma_j 0.5, R_nuc 100 | peak 0.40, capped 0%, G 1.22 | peak 0.68, capped 0%, G 1.46 |
| 0.5, sigma_j 0.5, R_nuc 100 | peak 4.3, capped 2.3%, G 4.5 | peak 7.0, capped 4.8%, G 6.8 |
| 0.05, sigma_j 1.5, R_nuc 300 | peak 4.1, capped 1.5%, G 4.0 | peak 32, capped 5.3%, G 13 |
| 0.5, sigma_j 1.5, R_nuc 300 | peak 4.4e3, capped 31%, G 6e4 (heavy-tailed) | peak 4.6e5, capped 83%, G 5e7 (heavy-tailed) |

Do not say: that the critical seed is ever Eddington-limited at the fiducial; do not put eta_acc = 0.5 in a figure or quote it as a
physical corner (half the nuclear reservoir consumed per free-fall time, beyond the paper's own eta_acc^(-1/beta)
ceiling); do not multiply a median G by a median M_seed,crit (heavy tails).

#### C9b. How far do fixed seeds get across the accessibility range? (D: `c17_accessibility_reach.log`, 60 trees, 801-step dt)

Direct per-tree M_BH(z=5)/M_star,tot(z=5), median [16,84]; "reach" = fraction of trees at or above f_BH = 0.5; Mcrit = median critical seed (f_BH = 0.5) and its
ratio to the fiducial cell (same trees in every cell).

3e10:

| sigma_j, R_nuc [pc] | Mcrit [Msun] (ratio to fiducial) | seed 1e3 | seed 2e5 | seed 1e7 |
|---|---|---|---|---|
| 0.5, 100 (fiducial) | 9.0e7 (1.00) | 4.4e-6 [3.7e-6, 5.3e-6] | 8.7e-4 [7.5e-4, 1.1e-3] | 4.5e-2 [3.9e-2, 5.5e-2] |
| 0.5, 250 (accessibility test) | 7.9e7 (0.88) | 6.4e-5 [1.7e-5, 2.8e-4] | 9.5e-4 [8.1e-4, 1.3e-3] | 4.7e-2 [4.1e-2, 5.7e-2] |
| 1.0, 100 | 3.6e7 (0.40) | 6.3e-3 [1.5e-3, 1.1e-2] | 2.8e-3 [1.5e-3, 5.2e-3] | 8.9e-2 [7.8e-2, 1.1e-1] |
| 1.5, 100 | 1.3e7 (0.14) | 7.2e-2 [5.2e-2, 9.4e-2] | 9.2e-2 [7.2e-2, 1.0e-1] | 0.42 [0.40, 0.45], reach 5.0% |
| 1.0, 300 | 2.8e7 (0.31) | 3.7e-2 [3.1e-2, 4.6e-2] | 3.9e-2 [2.9e-2, 4.9e-2] | 0.18 [0.14, 0.22] |
| 1.5, 300 | 1.4e7 (0.16) | 7.1e-2 [5.3e-2, 8.8e-2] | 8.9e-2 [7.2e-2, 1.0e-1] | 0.38 [0.34, 0.42], reach 3.3% |

3e11:

| sigma_j, R_nuc [pc] | Mcrit [Msun] (ratio to fiducial) | seed 1e3 | seed 2e5 | seed 1e7 |
|---|---|---|---|---|
| 0.5, 100 (fiducial) | 1.4e9 (1.00) | 2.9e-7 [2.3e-7, 3.7e-7] | 5.8e-5 [4.8e-5, 7.6e-5] | 2.9e-3 [2.4e-3, 3.8e-3] |
| 0.5, 250 (accessibility test) | 1.1e9 (0.82) | 6.3e-5 [3.8e-6, 2.7e-4] | 7.6e-5 [5.8e-5, 9.7e-5] | 3.1e-3 [2.5e-3, 4.0e-3] |
| 1.0, 100 | 5.3e8 (0.39) | 6.8e-3 [2.8e-3, 1.1e-2] | 1.7e-3 [5.0e-4, 6.2e-3] | 5.5e-3 [4.9e-3, 6.9e-3] |
| 1.5, 100 | 1.5e8 (0.11) | 9.2e-2 [6.9e-2, 1.2e-1] | 0.12 [8.5e-2, 0.14] | 0.17 [0.15, 0.18] |
| 1.0, 300 | 4.0e8 (0.29) | 5.4e-2 [4.6e-2, 5.8e-2] | 4.8e-2 [3.8e-2, 5.8e-2] | 4.2e-2 [3.4e-2, 5.7e-2] |
| 1.5, 300 | 1.7e8 (0.13) | 9.3e-2 [7.2e-2, 1.2e-1] | 0.11 [8.0e-2, 0.13] | 0.16 [0.13, 0.17] |

* In the accessible cells the black hole mass approaches about 0.1 to 0.2 of M_star,tot, with weak dependence on the seed (3e11, sigma_j = 1.5: 0.09, 0.12, 0.17 for
  seeds 1e3, 2e5, 1e7). At the fiducial and in the R_nuc = 250 pc test it does not (the ratio tracks the seed). No cell has more than 5 per cent of trees at f_BH = 0.5.
  This is a result of the explored models, not an established saturation scale: it may reflect the feedback prescription, the angular-momentum distribution, the nuclear radius or the star-formation law.
* The critical seed falls to 0.11 to 0.16 of its fiducial value at sigma_j = 1.5, and to 0.29 to 0.40 at sigma_j = 1.0, while remaining supply-limited (C9) and at least 1.3e7 Msun at 3e10.
* Do not say: that accessibility leaves the critical boundary unchanged; that increased accessibility "does not approach" f_BH = 0.5 (it is a factor of about 4 to 7 short for light seeds at sigma_j = 1.5); that the
  R_nuc = 250 pc test is representative of stronger accessibility. The R_nuc = 250 pc test leaves light seeds at least 2.9 dex below f_BH = 0.5 (C7).

### C10. Seed against host baryons

| Quantity | Result (P) | Caveat | Do not say |
|---|---|---|---|
| M_seed,crit / (f_b M_halo) at the first resolved step, median [16,84] | 5.1e3 [3.3e3, 1.7e4] (3e10), 4.9e4 [4.8e3, 1.3e5] (3e11), 3.8e3 [8.5e2, 2.6e4] (3e13); range 3.4e3 to 6.9e4 over the 13 masses | first resolved step = halo mass above M_res = 1e4 Msun | that the seed is assembled from the baryons in the resolved host |
| median first-resolved redshift | 21.4 (3e10), 24.0 (3e11), 24.8 (3e13) | the seed sits idle before this | |
| median M_BH / (f_b M_halo) at z = 15, 10, 5 | 10.5, 0.52, 0.021 (3e10); 33, 1.24, 0.031 (3e11); 6.7, 0.18, 0.0011 (3e13) | | |
| redshift where the median crosses 1 | 11.0 (3e10), 9.7 (3e11), 12.2 (3e13) | | |
| fraction of trees with M_BH > f_b M_halo (D: 240 trees, 801 steps) | z=15: 100%, 100%, 96%; z=10: 22%, 58%, 3%; z=7: 0% at all three (3e10, 3e11, 3e13) | | |

### C11. Dependence on f_BH

(D: 60 paired trees, 801 steps.) M_seed,crit ratio to f_BH = 0.5: f_BH = 0.1 gives 0.240, 0.233, 0.251; f_BH = 0.9 gives 1.58, 1.65,
1.59 (3e10, 3e11, 3e13). A near-linear rescaling because G_BH is close to 1: the shape of the boundary is not
changed. Do not say: that f_BH changes the shape.

### C12. Parameter sensitivities (D: 60 paired trees, 801 steps; ratio of median M_seed,crit to base)

| Variation | 3e10 | 3e11 | 3e13 |
|---|---|---|---|
| sigma_j 0.5 -> 0.75 | 0.67 | 0.69 | not run |
| eps_sf x0.5 / x2 | 0.77 / 1.20 | 0.69 / 1.37 | not run |
| R_nuc 50 / 200 pc | 1.10 / 0.89 | 1.11 / 0.87 | not run |
| f_mom 0.3 | 1.07 | 1.07 | not run |
| eta_SN x0.5 / x2 | 1.30 / 0.67 | 1.20 / 0.74 | not run |
| stellar wind off | 1.74 | 1.47 | not run |
| eta_acc x0.1 / x10 | 1.17 / 0.75 | 1.18 / 0.65 | 1.27 / 0.52 |
| radiative efficiency 0.057 / 0.32 | 1.04 / 0.90 | 1.04 / 0.91 | 1.05 / 0.87 |

In the first six rows (`sensitivities_60trees_*`) the M_seed,crit ratio equals the M_star,tot ratio and G_BH stays 1.03 to 1.08. The eta_acc and
radiative-efficiency rows (`claims_regime_and_sens_60trees_*`) recorded only M_seed,crit; M_star,tot and G_BH were not printed for them, and eta_acc x10 raises the
critical-seed G_BH to about 1.2 to 1.5 in the regime run (C9). Radiative efficiency (hence the Salpeter time) barely matters
because the Eddington limit does not bind at the critical seed. Range over the pruned set:
0.67 to 1.74. Do not say: that these ranges cover the full uncertainty of the accessibility prescription; that the
sensitivity to R_nuc of 50 to 200 pc says anything about R_nuc of order R_d or larger (the exploratory model did).

### C13. Numerical uncertainty

| Test | 3e10 | 3e11 | 3e13 |
|---|---|---|---|
| sampling: 60-tree median, std of log10 across 4 disjoint sets (D) | 0.0042 dex | 0.0078 dex | not run |
| sampling: two independent 240-tree sets, rep-to-rep (D) | 0.001 to 0.003 dex | 0.002 to 0.008 dex | not run |
| time step, ratio of median to 401 steps: 201 | 0.990 | 0.972 | 0.898 |
| time step, paired per-tree ratio to 401: 801 / 1601 | 1.005 / 1.008 | 1.013 / 1.018 | 1.045 / 1.065 |
| approximate offset of 801 below the refined limit (first-order extrapolation) | about 0.5% | about 1% | about 4% |
| dz = 0.1 / 0.025 relative to 0.05 (240 trees) | 0.985 / 1.013 | 0.988 / 1.033 | not run |
| M_res = 1e5 / 1e3 relative to 1e4 (240 trees) | 0.987 / 0.991 | 0.998 / 1.024 | not run |

* The time-step differences roughly halve with each doubling (first order). The offsets in the third row are estimates.
* dz drifts upward by 1 to 3 per cent per halving and is not flat at 0.025. dz is a resolution choice of the GRUMPY rate,
  not a parameter to be converged to zero. M_res shows no trend.
* Combined systematic on M_seed,crit is at most about 3 to 5 per cent at 3e10 and 3e11. At 3e13 the time-step offset alone is about 4 per cent and
  dz and M_res were not tested.
* Tree sampling is by the Zhang & Hui algorithm only. The alternative (PCH08) was not re-run in the MVM.
* Do not say: "converged" without these qualifications; that 801 steps is "the converged result" (it is the production choice).

### C14. Resolution-dependent diagnostics (option (a): quote only with this caveat)

(D: 60 paired trees, steps 201 / 401 / 801 / 1601.)

| Quantity | 3e10 | 3e11 | 3e13 |
|---|---|---|---|
| nuclear share of M_star,tot | 0.139 / 0.069 / 0.041 / 0.028 | 0.275 / 0.128 / 0.066 / 0.039 | 0.520 / 0.220 / 0.107 / 0.059 |
| share of M_star,tot formed at z > 10 | 0.072 / 0.035 / 0.018 / 0.009 | 0.101 / 0.052 / 0.026 / 0.013 | 0.485 / 0.218 / 0.106 / 0.061 |
| G_BH at the critical seed | 1.102 / 1.048 / 1.028 / 1.019 | 1.224 / 1.093 / 1.046 / 1.027 | 1.509 / 1.172 / 1.077 / 1.041 |
| M_star,tot ratio to 401 steps | 1.0 (401), 0.985 (801), 0.979 (1601) | 1.0, 0.969, 0.957 | 1.0, 0.961, 0.946 |

Cause: in steps where feedback empties the reservoirs, newly arrived gas gets one step of star formation before it is
ejected, so the stars formed there scale with dt. The share formed at z > 10 halves with each doubling of resolution.

Feedback-limited fraction (feedback demand exceeds the available gas in a step), at 801 steps (D: 240 trees): 38.4% (3e10),
43.7% (3e11), 26.3% (3e13) of the steps with any feedback; by redshift band 25-15 / 15-10 / 10-7 / 7-5: 100 / 99 / 32 / 0%
(3e10), 100 / 100 / 46 / 0% (3e11), 100 / 70 / 2 / 0% (3e13). At 401 steps (20 trees): 40.1% and 48.4%. These
also depend on the step. The qualitative statement that the model is feedback-limited before z of about 10 and supply-rich after z of
about 7 is stable across resolutions tested.

Do not say: nuclear share of stars, the amount of star formation before z = 10, G_BH - 1, or the feedback-limited fraction as converged
quantities. Do not describe the nuclear star formation as ordering-independent.

### C15. Pre-MVM against MVM on identical trees

(D: 20 trees, 401 steps.) Median M_seed,crit 9.35e7 -> 9.74e7 (3e10; per-tree ratio 1.034, range 1.019 to 1.047) and 1.51e9 -> 1.52e9
(3e11; ratio 1.006, range 1.005 to 1.011). The simplifications change M_seed,crit by 0.6 to 3.4 per cent. Do not say: that the
pre-MVM model is a separate result; it is a reference.

### C16. Baryon budget at z = 5 (fraction of the mass accreted by the halo; D: 240 trees, 801 steps)

| | stars | gas remaining | ejected by AGN | ejected by stellar wind |
|---|---|---|---|---|
| 3e10 | 0.038 | 0.191 | 0.231 | 0.540 |
| 3e11 | 0.060 | 0.349 | 0.187 | 0.399 |
| 3e13 | 0.133 | 0.129 | 0.345 | 0.383 |

Unavailable (failed the galaxy cut): 4e-4 to 6e-4 of the accreted mass. By step count the AGN wind dominates the ejecta,
cumulatively the stellar wind does. Do not say: that AGN feedback removes most of the gas (it removes 19 to 35 per cent at 801 steps).

### C17. Tests of the remaining structural items (branching from 695b114; the frozen model is unchanged)

Scripts and logs: `diagnostics/c17_diagnostics.py`, `c17_smooth_assembly.py`, `c17_pch08_check*.py`; `logs/c17_*.log`.
Each script refuses to run unless `reservoir_stock.py` matches the hash in the production JSON. 100 paired trees per mass, 801-step production dt.

**z_seed** (15, 20, 25, 30, 35; step fixed at the production dt; one tree set per mass with z_max = 40). Paired per-tree ratio of M_seed,crit to z_seed = 25:

| Mass | z_seed = 15 | 20 | 30 | 35 |
|---|---|---|---|---|
| 3e10 (median at 25: 9.90e7) | 1.0002 [1.000, 1.003] | 1.0000 | 1.0000 | 1.0000 [0.999, 1.001] |
| 3e11 (1.36e9) | 1.0003 [1.000, 1.006] | 1.0000 | 1.0000 | 1.0000 [0.999, 1.001] |
| 3e13 (4.85e9) | 1.0040 [1.000, 1.034] | 1.0001 [0.996, 1.003] | 1.0000 [0.998, 1.001] | 0.9998 [0.993, 1.004] |

M_seed,crit is insensitive to z_seed from 15 to 35. Reason: before the host halo exists nothing happens, and gas accreted early is tiny and
ejected. The fraction of trees whose halo already exists at z_seed is 3% (3e10), 34% (3e11), 98% (3e13) at z_seed = 25 and 100% at z_seed = 15.
Caveat: the MVM starts with an empty galaxy at z_seed, so a later z_seed also discards baryons already accreted. Do not say: that this is the
old paper's z_seed argument (the old text's reason, "the halo has assembled little mass", does not apply as stated); that z_seed = 25 is
a physical formation redshift.

**M_hot** (x0.5, x1, x2 of 4e11 Msun; same trees per mass). Median M_seed,crit [Msun] and paired ratio to the fiducial M_hot:

| M0 | M_hot = 2e11 | 4e11 | 8e11 | ratio (2e11, 8e11) | M_star,tot/M_halo (2e11, 4e11, 8e11) |
|---|---|---|---|---|---|
| 3.0e10 | 9.27e7 | 9.28e7 | 9.28e7 | 1.000, 1.000 | 6.4e-3, 6.4e-3, 6.4e-3 |
| 9.5e10 | 3.77e8 | 3.82e8 | 3.83e8 | 0.987, 1.002 | 8.2e-3, 8.3e-3, 8.4e-3 |
| 3.0e11 | 1.09e9 | 1.41e9 | 1.48e9 | 0.778, 1.049 | 7.6e-3, 9.8e-3, 1.03e-2 |
| 9.5e11 | 1.60e9 | 2.94e9 | 4.60e9 | 0.534, 1.543 | 3.6e-3, 6.6e-3, 1.03e-2 |
| 3.0e12 | 1.87e9 | 3.97e9 | 7.93e9 | 0.466, 1.997 | 1.3e-3, 2.8e-3, 5.7e-3 |
| 9.5e12 | 2.04e9 | 4.56e9 | 9.59e9 | 0.449, 2.119 | 4.6e-4, 1.0e-3, 2.2e-3 |
| 3.0e13 | 2.14e9 | 4.87e9 | 1.06e10 | 0.444, 2.177 | 1.5e-4, 3.5e-4, 7.7e-4 |

The low-mass end (3e10) does not depend on M_hot. Above about 1e12 the boundary scales almost linearly with M_hot (x0.44 to 0.47 for M_hot halved,
x2.0 to 2.2 for doubled), through the accreted baryon fraction (accreted/(f_b M_halo) at 3e13 = 0.008, 0.017, 0.034). The kink moves with M_hot:
the local log-slope drops from about 1.2 to below 0.3 at lower halo mass for smaller M_hot. Only the location of the step was varied, not its sharpness (phi = 4).
Do not say: that the high-mass boundary is independent of M_hot, that M_seed,crit(3e13) is a robust number without the qualifier "for M_hot = 4e11",
or that the hot-mode step is tested beyond its location.

**Halo-assembly description.** The alternative tree algorithm, PCH08, cannot be used at the required resolution (evidence in `logs/c17_pch08_check.log`):

* foraois PCH08 gives a near-deterministic, much earlier main-progenitor history at M_res = 1e4 Msun: for M0 = 3e10 anchored at z0 = 5, M(z=10)/M0 = 0.40
  (16-84 = 0.40, 0.40) against 0.035 [0.016, 0.082] for Zhang & Hui; the same at every dz from 0.02 to 0.2 and with the numpy and numba backends.
* At the standard anchor z0 = 0, M0 = 1e12: M(z=1)/M0 = 0.843 [0.84, 0.84] (PCH08) against 0.535 [0.36, 0.76] (Zhang & Hui).
* Cause: at M_res = 1e4 PCH08 resolves a merger in 2.3% of steps (Zhang & Hui: 41%) and books almost all accretion as smooth (merged/smooth mass 8e-6).
  PCH08 resembles Zhang & Hui only at M_res = 1e10 (about 1% of M0): M(z=1)/M0 = 0.540 [0.343, 0.662] against 0.495 [0.302, 0.667]. That cannot follow
  the assembly to z of about 20, which the MVM needs.
* The MVM run on those PCH08 trees (M_seed,crit 1.00e8, 1.68e9, 2.72e9 with zero scatter) is therefore not a result and must not be quoted.

As an independent description of assembly the frozen MVM was run on the deterministic Fakhouri, Ma & Boylan-Kolchin (2010) mean accretion history
(calibrated at z below about 2 and extrapolated here), anchored to M_halo(z=5) by shooting (`c17_smooth_assembly.py`):

| M0 | smooth M(z=10)/M0 | M_seed,crit (smooth) | Zhang & Hui median | ratio | M_star,tot/M_halo (smooth / Zhang & Hui) |
|---|---|---|---|---|---|
| 3e10 | 0.103 | 8.48e7 | 9.45e7 | 0.897 | 5.8e-3 / 6.5e-3 |
| 3e11 | 0.061 | 1.27e9 | 1.40e9 | 0.908 | 8.8e-3 / 9.7e-3 |
| 3e13 | 0.016 | 4.56e9 | 4.85e9 | 0.939 | 3.3e-4 / 3.5e-4 |

The boundary changes by 6 to 10 per cent between two assembly descriptions whose halos differ by a factor of 2.5 in mass at z = 10 (at 3e10). Do not say: that PCH08
was tested, or that the result is independent of every assembly algorithm; say that it is insensitive to the difference between the Zhang & Hui ensemble and a smooth mean history.
The statement of the old paper, "PCH08 within 0.02 to 0.14 dex", cannot be reproduced and must not be reused.

**n_rd** (galaxy scale in disc scale lengths; 5, 10 fiducial, 20; paired). M_seed,crit paired ratio to n_rd = 10 [min, max]:

| M0 | n_rd = 5 | 20 | M_star,tot/M_halo (5, 10, 20) | G_BH (5, 10, 20) |
|---|---|---|---|---|
| 3e10 | 1.112 [1.041, 1.370] | 0.792 [0.676, 0.868] | 7.3e-3, 6.5e-3, 5.1e-3 | 1.026, 1.028, 1.029 |
| 3e11 | 1.232 [1.093, 1.424] | 0.714 [0.627, 0.792] | 1.23e-2, 9.9e-3, 7.0e-3 | 1.047, 1.047, 1.046 |
| 3e13 | 1.057 [1.012, 1.151] | 0.802 [0.731, 0.848] | 3.7e-4, 3.5e-4, 2.8e-4 | 1.078, 1.078, 1.076 |

n_rd acts only through M_star,tot (G_BH does not change): it sets R_gal and hence the galaxy-scale free-fall time in the star-formation clock. Its effect (about 6 to 23 per cent
for a halving, 20 to 29 per cent for a doubling) is comparable to the other star-formation sensitivities of C12. Do not say: that n_rd is an accessibility parameter.

**Still untested in the MVM:** the UV-suppression parameters and z_reion, lambda, f_b, the NFW concentration (only c = 3 and 8 in the exploratory model, at most 4 per cent,
so it matters only at the galaxy-scale threshold), the sharpness of the cold/hot step, and the 3e13 mass for sigma_j, eps_sf, R_nuc, f_mom and eta_SN.
Any statement about them in the paper is untested by the frozen model.

## 2. Consolidated "do not say" list

1. The Eddington limit is irrelevant. (It is irrelevant to the critical boundary at the fiducial; it may be decisive for light seeds near the accessibility threshold.)
2. The model predicts that light seeds cannot grow. (At the fiducial they do not; in the accessibility test they grow by factors of 10 to 10^4.)
3. The accessibility prescription necessarily produces bimodality. (It appears in the accessibility test at high halo mass only.)
4. The R_nuc = 250 pc experiment is a second physical model, or a candidate fiducial.
5. The Figure 2 high-mass median drop is a decline in individual black hole growth.
6. Secondary star-formation or timing quantities (nuclear share, star formation at z > 10, G_BH - 1, feedback-limited fraction) are converged.
7. The 1e7 Msun seed is a seed-formation channel. It is a benchmark seed mass.
8. M_seed,crit is a prediction of a seed-formation mass. It is the initial mass required to reach the target ratio in this baryon cycle.
9. The kink in M_seed,crit or the low stellar mass at 3e13 is emergent. It follows from the imposed zeta_ch step at M_hot.
10. The light-seed threshold location is precise, or that the fiducial was placed relative to it.
11. The fiducial critical-seed result depends on proximity to the light-seed threshold. (It does not: see C9.)
12. Numerical convergence is established at 3e13 for dz and M_res, or that 801 steps is converged.
13. Anything about the Eddington time as the operative constraint on the critical seed.
14. The high-mass boundary (above about 1e12) is independent of M_hot; it scales almost linearly with it.
15. PCH08 (or any second stochastic tree algorithm) was tested. Only a smooth mean history was, and PCH08 is not usable at this resolution.
16. z_seed matters (or does not) for the reason the old paper gave; in the MVM nothing happens before the host exists.
17. n_rd is an accessibility parameter; it enters through the star-formation clock.

## 3. Old-model numbers that must not be reused

They come from the earlier three-reservoir model (nuclear-only star formation, memoryless nuclear fraction, energy-driven
wind, graded Eddington cap) and, for the tree-based ones, from trees with the P(k)-redshift error. They are not comparable to
anything above: M_seed,crit 7.3e4 (3e10), 1.6e5 (5e10), 3.9e6 (1e12), 4.8e6 (1e13); slopes 0.9 and 0.1; scatter 0.12 to 0.31 dex;
eta_acc 1.9e6 -> 2.7e4; Phi_hat factor 16; sigma_j 8.9e4 -> 4.9e5; eps_f 1.5e6 -> 2.7e4 (factor 57); f_BH 5.5e4 -> 5.4e5;
z_seed change below 2 per cent; radiative efficiency 1.6e5 -> 6.9e5; PCH08 within 0.02 to 0.14 dex. The channel comparison
(direct-collapse overlap at low halo mass) does not survive: in the MVM even a 2e5 Msun seed reaches only 8e-4 of M_star,tot at 3e10.

## 4. Provenance

Production: `gen_mvm_production.py` -> `output/mvm_production_results.json`; figures `plot_mvm_fig{1,2,3}_*.py` ->
`output/mvm_fig{1,2,3}_*.png`; numbers (P) from `mvm_crib_numbers.py`. Diagnostics (D): `diagnostics/` (README, scripts, logs).
Tests: `tests/test_reservoir_stock.py` (MVM) and `tests/test_reservoir_stock_premvm.py` (frozen reference).
The manuscript and its repository have not been modified.
