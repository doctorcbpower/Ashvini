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
