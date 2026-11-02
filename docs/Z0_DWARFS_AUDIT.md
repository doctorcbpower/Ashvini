# Pre-flight audit: a z = 0 run on OzSTAR to compare star formation laws with local dwarfs

Written 2026-11-02 before any code or job. Sources read: Pan & Kravtsov (arXiv:2310.08636, all 14 pages), `docs/RESOLUTION_CONVERGENCE.md`, `ashvini/pymctrees_adapter.py` (tree options), the Stochastic_Occupation `hpc/README_HPC.md`, and the earlier results in this repository. Nothing here has been run.

## 1. The question

Run to z = 0, does a star formation law set by a disc (the Schmidt-Kennicutt option, with GRUMPY-like disc sizes) give a better match than the current law to the observed relations of local dwarfs: stellar mass against halo mass, gas mass against stellar mass, gas and stellar metallicity against stellar mass, and star formation rate against stellar mass (the relations of Pan & Kravtsov, Figs. 3, 4, 7)? And how much does the answer depend on the other choices (efficiency, disc size, reionisation redshift, wind model)? This is the sensitivity map applied to data, not a fit.

## 2. Assumptions the run would make

1. Each halo is one gas reservoir evolved along its main branch from z = 25-30 to z = 0, with its mass history from a tree generator (foraois).
2. The reionisation suppression (Okamoto-style, `reion.z_reion` = 7, rate-dependent form in the grid integrator) acts on every halo; no other environment.
3. Star formation: the current law (efficiency 0.015 per halo free-fall time) and the disc-based variants (`disc.star_formation_law`, disc size = `angular_momentum_ratio * spin * R_vir / sqrt(2)`; GRUMPY's mean is R_d/R_vir = 0.06, i.e. `angular_momentum_ratio` = 2.42 at spin 0.035).
4. Supernova feedback: the default delayed mass loading (15 Myr); the wind model and delay are further variants.
5. Outputs at z = 0: stellar mass, gas mass, gas metallicity, stellar metallicity, star formation rate; stored per halo.

## 3. What you specified

Run to z = 0 (to compare with local dwarfs as Pan & Kravtsov do), and on OzSTAR.

## 4. What I inferred

The comparison relations and the variants above; one halo mass range 1e8-1e12 Msun at z = 0; Zhang-Hui trees (see 5); the grid integrator only (the event integrator is too slow at 1e11 and above, and is not needed for median relations); about 100 haloes per mass.

## 5. Where this conflicts with what is already known

- **Tree choice and resolution.** The earlier work (`docs/RESOLUTION_CONVERGENCE.md`) found that the low-mass stellar mass at z = 0 grows with finer mass resolution `M_res` on PCH08 trees, an effect of the fixed time step becoming non-compliant (not of Ashvini), while Zhang-Hui trees were nearly converged in `M_res` (11 to 34 percent between `M_res/M0` = 1e-2 and 1e-5). The time step of the Zhang-Hui trees has not been checked against the expected number of splits. So the run needs a deliberate choice, and a check, not the defaults.
- **Field haloes against satellites.** Pan & Kravtsov compare Milky Way satellites (subhaloes, with tidal effects and `M_peak`) with data from satellites and the field. Ashvini evolves isolated main branches, so it corresponds to field or pre-infall haloes. The comparison with satellites is not like for like.
- **Observational data.** Only `Stochastic_Occupation/data/obs/simon2019_table1.csv` (stellar metallicities) is on this machine. The gas, SFR and halo-mass relations need data files (ALFALFA, SHIELD, Kirby et al. 2013, Berg et al. 2016, McConnachie 2012, and an abundance-matching relation).
- **GRUMPY is not reproduced.** Its star formation uses molecular gas with a fixed depletion time of 2.5 Gyr; the Schmidt-Kennicutt option does not, so the comparison tests "a disc-based rate", not GRUMPY.
- **The disc is out of equilibrium** at z > 5 and the disc size, taken from a constant spin, is an assumption.

## 6. Decisions that change the answer materially

A. **Tree generator and settings.** Zhang-Hui (nearly converged in `M_res` earlier), dz and `M_res/M0` chosen so that the splits per step are below 0.1 (the adapter warns), `z_max` 30. Recommendation: Zhang-Hui, and a convergence check on one mass at two values of dz and `M_res/M0` before the full run.
B. **Which relations and data.** At least the stellar-to-halo mass relation, gas mass against stellar mass, and gas and stellar metallicity against stellar mass. Needs data files from you (or a pointer to where you keep them).
C. **Variants.** Current law (efficiency 0.005, 0.015, 0.05); Schmidt-Kennicutt at disc sizes 0.025, 0.06 (GRUMPY) and 0.06 with 0.25 dex scatter; slope 1.0 and 1.4; reionisation redshift 6 and 8. About 10 variants.
D. **Sample.** 9 halo masses from 1e8 to 1e12 Msun at z = 0, 100 trees each, one tree set shared by all the variants.
E. **Field haloes only** for this run, with the satellite comparison left for later, or add an accretion-time cut to the halo selection.
F. **What counts as agreement.** The median relation within the observed 1-sigma band, with no tuning of parameters to the data. Stated before the run so that the comparison is a test.

## 7. Minimum checks before and after

1. A local pilot (3 masses, 5 trees, 3 variants, minutes) that checks the pipeline, the output files and the timing, and gives the cost for the full run.
2. Tree convergence at two settings of dz and `M_res/M0` (item A).
3. The default variant reproduces the stellar-to-halo mass relation of the earlier z = 0 runs (`scripts/output/shmr_dm_model_comparison_zh.json`) within the scatter.
4. Mass balance on a few haloes (gas + stars + removed mass = baryons supplied).
5. A seeded rerun of one shard gives identical output.
6. Plots made from the stored outputs only, comparing medians and 1-sigma bands with the data, with the relations and bands defined before looking.

## What running on OzSTAR needs from you

The Slurm account and partition, the Python module lines, a scratch path, and whether to clone Ashvini and foraois from GitHub on the cluster (I'd record the foraois commit, as in the Stochastic_Occupation `hpc/` directory). Compute nodes may lack internet, so the files cannot be downloaded there. Plan: a tree-generation array (one task per halo mass) writing HDF5 in Ashvini's format, then a variant array (variant by mass) running the grid integrator on those trees, both resumable and seeded, then a local plotting script.

## Not in scope

Subhalo or satellite physics, tidal stripping, a molecular-gas model, the stochastic star formation perturbation of Pan & Kravtsov, and new physics.

## Results of the first run (2026-11-03; 100 haloes at each of nine masses, Zhang-Hui trees, grid integrator)

Tables: `scripts/output/shmr_z0/summary.txt`; plot: `shmr_z0_laws.png`. Median stellar mass relative to the published relations; Moster+13 and Behroozi+13 agree within a factor 2.2 between 3e10 and 1e12 Msun and diverge below (factors 4, 12, 36 and 340 at 1e10, 3e9, 1e9 and 1e8): the relations are extrapolations below about 3e10.

- **Star formation law.** With the default (rate-dependent) UV suppression the law matters only above about 3e10 Msun. The Schmidt-Kennicutt law gives 1.1 times Moster+13 at 1e11-1e12 where the current law gives 0.2-0.3, but 2.4 times Behroozi+13 at 1e11. Below 1e10 all laws give the same plateau, M*/M_halo of about 2e-5, a factor 18-70 below the relations at 1e10.
- **Where the plateau and the step come from.** With the UV background off there is no plateau and no step (the dwarfs are then 2-3 times above Behroozi). With the default UV form 60-99 percent of the stars in haloes below 1e10 Msun form before z = 7, and almost none after z = 3.
- **Equilibrium UV form** (`reion.form = "equilibrium"`). The step disappears. With the current law (efficiency 0.015) the median is 1.1-1.6 times Behroozi+13 from 1e10 to 1e11 Msun, 0.5 and 0.33 times at 3e11 and 1e12. Stars in haloes of 3e8 Msun or less are fossils (more than 80 percent formed before z = 7), with M*/M_halo of about 1e-5 to 2e-5; at 1e9 the transition has begun (39 percent before z = 7). With the Schmidt-Kennicutt law the massive end reaches 1.3-1.7 times Moster at 1e12 but is 4.5-7 times too high at 1e10-3e10.
- **No single efficiency or law gives the shape.** M*/M_halo rises by about 30 between 1e10 and 1e12 in both published relations and by about 4 (current law) in Ashvini.
- **Found on the way:** `reionization.epsilon` gave nan for `z_reion` of about 6.5 or less (fixed, with a test).
- **Not tuned.** No parameter was adjusted to match either relation.
