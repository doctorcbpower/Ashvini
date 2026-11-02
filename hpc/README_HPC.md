# Running the z = 0 star formation law comparison on a cluster (Slurm)

What this runs: the stellar-to-halo mass relation at z = 0 for 10 variants of Ashvini's star formation (the current law at three efficiencies,
the Schmidt-Kennicutt disc law at several disc sizes and slopes, two reionisation redshifts), on Zhang-Hui trees at nine halo masses from 1e8 to
1e12 Msun, compared with the Moster et al. (2013) and Behroozi et al. (2013) relations. See `docs/Z0_DWARFS_AUDIT.md` for the question, the
assumptions and the checks. I have not run this on OzSTAR; the Slurm account, partition and Python module lines are placeholders marked TODO.

## Cost (measured on a laptop, 2026-11-03)
- Trees: about 35 s per halo mass, almost independent of the number of haloes (the time is the set-up, not the trees), so 1000 haloes per mass costs about the same as 5.
- Variants: seconds per variant and mass for 100 haloes; roughly linear in the number of haloes.
- So the default size (100 haloes, 9 masses, 10 variants) is about 15-20 minutes on one core; the cluster is for larger samples (`SHMR_N_HALOES=1000`, more masses, more variants).
- Memory: about 1 GB per task at 100 haloes (a tree array is haloes x 6001 x 8 bytes).

## Steps
1. On the cluster, clone Ashvini and foraois (the runs here used foraois commit `94ddf3a`, see `hpc/foraois_commit.txt`), or unpack the tarballs from `bash hpc/make_package.sh`.
2. Edit the `module` lines in `hpc/setup_env.sh` and the `--account`/`--partition` lines in the two `.sbatch` files, then `bash hpc/setup_env.sh <path-to-foraois>`.
3. Submit: `bash hpc/submit_all.sh` (trees, then the variants once the trees are done). Or by hand:
   `sbatch --array=0-8 hpc/submit_trees.sbatch` then `sbatch --array=0-8 --dependency=afterok:<jobid> hpc/submit_variants.sbatch`.
4. Both steps are resumable: finished files are skipped, so a failed task can be resubmitted with the same command (`--array=3,7`).
5. Copy `$SHMR_OUT` (default `scripts/output/shmr_z0`) back, then `python scripts/shmr_z0_plot.py` locally.

## Things to know
- Compute nodes may lack internet: nothing here downloads anything (the foraois configs are read from `FORAOIS_CONFIG_DIR`).
- Everything is seeded (tree seed 42 + mass index), so a rerun on the same software gives identical trees.
- Environment variables: `SHMR_OUT`, `FORAOIS_CONFIG_DIR`, `SHMR_N_HALOES` (see `scripts/shmr_z0_common.py`).
