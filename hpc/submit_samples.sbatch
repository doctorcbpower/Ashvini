#!/bin/bash
#SBATCH --job-name=space_samples
#SBATCH --account=TODO            # TODO: your Slurm account
#SBATCH --partition=TODO          # TODO: a CPU partition
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/space_samples_%A_%a.out
# Usage (from the repository root): mkdir -p logs; sbatch --array=0-199 hpc/submit_samples.sbatch
# Task k evaluates shard k of every sub-space (samples k, k + 200, ... of each); about 35 samples of 15 s each. Resume-safe: progress is
# saved every 5 samples and a resubmitted task continues where it stopped. Needs the stored trees of hpc/submit_trees.sbatch.
set -euo pipefail
# module load python/3.11        # TODO: same module as in setup_env.sh
source venv/bin/activate
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
N_SHARDS="${N_SHARDS:-200}"
cd scripts
for sub in cur_ash cur_en sk_ash sk_en cur_ash_rate; do
  python space_eval.py "$sub" "${SLURM_ARRAY_TASK_ID:?run with --array=0-199}" "$N_SHARDS"
done
