# Pre-flight audit: exploring Ashvini's parameter space for the z = 0 stellar-to-halo mass relation on OzSTAR

Written 2026-11-03 before any code or job. Based on the z = 0 results in `docs/Z0_DWARFS_AUDIT.md` and the sensitivity tables from the earlier runs. Nothing here has been run.

## 1. The question

Which parameters control the z = 0 stellar mass of haloes of each mass, how strongly, and which combinations are consistent with both published relations where they are constrained (above about 3e10 Msun)? What do those consistent combinations then predict at 1e8 to 3e9 Msun, where the published relations are extrapolations? This is the sensitivity map made quantitative, plus a calibrate-high, predict-low test. It is not a fit to the relations at low mass.

## 2. What the exploration assumes

1. The trees are fixed: Zhang-Hui, 100 haloes at each of nine masses from 1e8 to 1e12 Msun (the stored trees of the z = 0 runs), so differences between samples come from the parameters, not the trees.
2. The grid integrator with the equilibrium UV suppression (`reion.form = "equilibrium"`) is the fiducial. The rate-dependent form depends on the tree time step and stops accretion in small haloes, so it is run only as a comparison.
3. Parameters that did not change the stellar mass in the earlier sensitivity runs (dust yield and destruction, the metallicity of the inflowing gas) are left out.
4. Output of one sample: for each of the nine masses, the median, 16th and 84th percentiles of M*, and the medians of gas mass, gas metallicity, stellar metallicity and black hole mass, and the fraction of haloes with stars. A few seconds per mass, about 14 s per sample.

## 3. What you specified

Explore the parameter space on OzSTAR, with the aim of seeing how close Ashvini can get to either relation.

## 4. What I inferred

The parameter list and ranges, the sampling design, the acceptance criterion, and the four-way split by star formation law and wind model, all below.

## 5. Where this may mislead

- **The relations are extrapolations below about 3e10 Msun**, and differ by factors up to 340 at 1e8. Calibrating there would tune the model to an extrapolation, so the criterion uses only 3e10 and above.
- **Both relations agree only within a factor 2.2**, so "consistent with both" means inside the envelope of the two, widened by a tolerance. The width of the tolerance decides how many samples are accepted.
- **Fixed trees** understate the variation due to halo assembly. The scatter between haloes is kept (16th to 84th percentile), but one tree set is one realisation.
- **Degeneracies.** The mass loading amplitude and the star formation efficiency, for example, partly compensate (stars scale as the efficiency divided by the mass loading in regulated haloes). A sensitivity index measures effect, not identifiability; the accepted region will show the degeneracies directly.
- **Hard-coded values stay fixed:** the mass and redshift scaling of the Ashvini mass loading (exponents 1/3 and 1/2), the metallicity step of the mass loading, and the IMF. If the missing steepness needs a different mass scaling, that is a model change, not a parameter, and needs your agreement.

## 6. Decisions that change the answer

A. **Parameters and ranges** (log-uniform where marked). Continuous, per sub-space:

| Parameter | Range | Applies to |
|---|---|---|
| `sf.efficiency` | 0.003 to 0.1, log | current law |
| `disc.schmidt_kennicutt_normalisation` (A) | 8e-5 to 8e-4 Msun/yr/kpc², log | Schmidt-Kennicutt |
| `disc.schmidt_kennicutt_slope` | 1.0 to 1.6 | Schmidt-Kennicutt |
| `disc.angular_momentum_ratio` (sets R_d/R_vir = 0.012 x ratio) | R_d/R_vir 0.03 to 0.12 | Schmidt-Kennicutt |
| `sn.epsilon_p` | 1 to 25, log | Ashvini wind |
| `sn.eps_fb` | 0.01 to 1, log | energy-driven wind |
| `sn.delay_time` | 3 to 50 Myr | both |
| `reion.z_reion` | 6 to 9 | both |
| `reion.gamma`, `reion.omega` | 10 to 20; 1 to 3 | both |
| `bh.efficiency` | 3e-4 to 1e-2, log | both |
| `bh.eta_agn` | 0.05 to 2, log | both |

B. **Design.** Saltelli sampling for variance-based (Sobol) sensitivity indices, N = 128 base samples in each of four sub-spaces (star formation law: current or Schmidt-Kennicutt; wind: Ashvini mass loading or energy-driven), about 10 parameters each: roughly 2800 samples per sub-space, 11 000 in all, about 40 CPU hours, which is 200 array tasks of about 12 minutes. The same samples also give the accepted region (item D). Recommendation: this.
C. **Fiducial UV form:** equilibrium, with one extra sub-space on the rate-dependent form to show the difference. Recommendation: this.
D. **Acceptance criterion, set before looking.** For halo masses 3e10, 1e11, 3e11 and 1e12 Msun, the median stellar mass must lie within 0.2 dex of the interval between the two published relations (below the lower and above the upper). Accepted samples are then summarised at 1e8 to 3e9 as a band, with the fraction of samples accepted reported for each sub-space. Recommendation: this, but the 0.2 dex is your call.
E. **Outputs analysed:** the Sobol first-order and total indices of log M* at each mass, for each sub-space, with bootstrap confidence intervals; the accepted region's parameter distributions; the predicted low-mass band.
F. **Sample size of the trees:** 100 haloes per mass. A larger set (1000) would reduce the noise in the medians, at the cost of a factor 10 in time. Recommendation: 100, and a rerun of the accepted samples with 1000 haloes.

## 7. Minimum checks

1. A local pilot: 20 samples across the full ranges, to confirm the output files, the timing, that nothing is NaN or inf, and that a sample at the default parameters reproduces the stored run.
2. A repeat of the first 50 samples gives identical outputs (seeded, no randomness inside a sample).
3. The Sobol indices of a test function with known indices (Ishigami) computed by the same code.
4. A convergence check: indices from N = 64 and N = 128 agree within the bootstrap errors for the most important parameters.
5. Plots made from stored outputs only, with the criterion, the ranges and the sub-spaces as stated here.

## What running on OzSTAR needs from you

The same details as before: Slurm account, partition, module lines, scratch path. The code would be a new `hpc/submit_samples.sbatch` array and `scripts/space_*.py` (sample generation, evaluation of one shard, analysis). The trees are generated on the cluster by the existing `hpc/submit_trees.sbatch`.

## Not in scope

New physics, changes to the hard-coded mass loading scaling, subhaloes, and a fit to the dwarf end of either relation.

## Results of the run (2026-11-03, OzSTAR; 6912 samples, nine masses, 100 haloes each)

Summary and figures: `scripts/output/shmr_z0/space/` (`space_summary.txt`, `space_sobol_total.png`, `space_accepted_band.png`); the shard files (8 MB) are not in the repository. All samples finished with no non-finite values; the analysis run locally on the copied shards reproduces the cluster's summary exactly.

- **Checks.** Total Sobol indices of the dominant parameters at 1e11 Msun agree between N = 64 and N = 128 (the first half of the Sobol sequence) to 0.03 (cur_ash: efficiency 0.47 and 0.44, mass loading amplitude 0.66 and 0.68). The first-order indices are not usable at N = 128: their bootstrap half-widths (0.1 to 7) exceed the values, so only the total indices are interpreted.
- **What controls M*.** The amplitude of the supernova wind (`epsilon_p`, or `eps_fb` for the energy-driven wind) has total index 0.6 to 1.0 from 1e9 to 3e11 Msun in every sub-space. The star formation efficiency matters only above 1e11 and only for the current law (0.5 to 0.9 at 1e11 to 1e12); with the Schmidt-Kennicutt law its parameters (A, slope, disc size) have total index 0.01 to 0.1: the haloes are regulated by the wind. Reionisation (`z_rei`, `omega`) matters below 1e9 (0.1 for the equilibrium form) and, for the rate-dependent form, up to 1e10 (`z_rei` 0.3 to 0.4). `gamma`, the black hole accretion efficiency and the AGN wind coupling have total index 0.00 at every mass; the supernova delay matters (0.01 to 0.1) only below 1e9. The black holes do grow (median 1e6, 5e7 and 8e8 Msun at 1e10, 1e11 and 1e12 in cur_ash), so the AGN wind of the current model does not affect the stellar mass over this range.
- **Acceptance** (median M* within 0.2 dex of the interval between Moster+13 and Behroozi+13 at 3e10, 1e11, 3e11, 1e12): 4.3 percent (current law, Ashvini wind), 5.0 percent (current, energy wind), 14.6 percent (Schmidt-Kennicutt, Ashvini wind), 13.9 percent (Schmidt-Kennicutt, energy wind), 4.9 percent (rate-dependent UV).
- **Accepted parameters.** Wind amplitude `epsilon_p` 7.4 to 10.5 (default 5; Schmidt-Kennicutt 7.5 to 11.5) or `eps_fb` 0.18 to 0.31; efficiency 0.047 to 0.09 (default 0.015) for the current law; for Schmidt-Kennicutt A = 1e-4 to 4.6e-4 (Kennicutt: 2.5e-4), slope 1.1 to 1.5 (Kennicutt: 1.4) and R_d/R_vir 0.04 to 0.105 (GRUMPY mean: 0.06). The efficiency and the wind amplitude are only weakly correlated among the accepted samples (correlation of the logarithms 0.13).
- **Prediction below the constrained masses** (accepted samples, equilibrium UV form; median and 16-84 percent of M*/M_halo): 1e10 Msun: 1.1e-3 to 1.3e-3 (Moster 4.6e-4, Behroozi 1.8e-3); 3e9: 1.4e-4 to 2.2e-4 (Moster 9.3e-5, Behroozi 1.1e-3); 1e9: 1.9e-5 to 3.8e-5 (Moster 1.9e-5, Behroozi 7.0e-4); 1e8 to 3e8: 1e-5 to 2e-5 (Moster 8e-7 to 4e-6, Behroozi 2.7e-4 to 4.3e-4). So the model, calibrated only above 3e10 Msun, lies within a factor 3 of Moster+13 down to 1e9, between the two relations at 1e10, and far above Moster and 20 times below Behroozi at 1e8, where it reaches a floor set by stars formed before reionisation. With the rate-dependent UV form the accepted samples fall below Moster at 1e10 (1.4e-4 against 4.6e-4) and have no step-free rise.
- **Not shown by this run.** Whether either relation is right below 1e10 Msun; the effect of the tree realisation (one tree set); anything about the parameters that are fixed (the mass and redshift scaling and the metallicity step of the Ashvini mass loading).

### The accepted parameter sets at other redshifts (`scripts/space_sequence.py`; 12 accepted samples in each of the four equilibrium sub-spaces)

M*/M_halo against M_halo(z) of the main progenitors, pooled over haloes and samples (median; `scripts/output/shmr_z0/space/space_shmr_redshift.txt`, `.png`). At z = 5: about 2e-4 to 3e-4 (cur_ash, cur_en; 6e-5 to 1e-4 for sk_en, sk_ash) between 2e6 and 2e7 Msun, 4e-4 to 6e-4 at 1e8, 1.3e-3 at 6e8, 2e-3 at 2e9, 3e-3 to 4e-3 at 6e9, 4e-3 to 6e-3 at 2e10 and 4.5e-3 to 9e-3 at 5e10; the 16-84 percent range is a factor 3 to 30 at the low-mass end. At z = 5 the relation is therefore flat below about 1e8 Msun (stars formed in the first haloes, before reionisation) and about 10 to 30 times higher than at z = 0 at the same mass ratio floor (1e-5 at z = 0). Moster et al. (2013) extended to z = 5 with coefficients typed from memory (not checked against the paper) gives 2e-5 at 1e8, 1e-4 at 1e9, 4.8e-4 at 1e10 and 2.3e-3 at 1e11: the model is 20 times above that curve at 1e8, 16 at 1e9, 7 at 1e10 and 3 to 6 at 5e10, and the offset grows with redshift (about 3 times at z = 2). The model was calibrated only at z = 0, and neither relation is constrained by data at z = 5, so this is a prediction and not a test against either.
