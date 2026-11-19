#!/bin/bash
# Parameter-space exploration: generate the trees if they are missing, then evaluate the 6912 samples in 200 tasks.
# Run from the repository root after hpc/setup_env.sh. Then copy $SHMR_OUT/space back and run scripts/space_analyse.py locally.
set -euo pipefail
mkdir -p logs
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
if ls "$SHMR_OUT"/trees/trees_8.npz >/dev/null 2>&1; then
  dep=""
else
  tree_job=$(sbatch --parsable --array=0-8 hpc/submit_trees.sbatch)
  echo "trees: job $tree_job"
  dep="--dependency=afterok:$tree_job"
fi
sbatch --array=0-199 $dep hpc/submit_samples.sbatch
