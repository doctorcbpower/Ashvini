#!/bin/bash
# Parameter-space exploration: task k evaluates shard k of every sub-space (about 35 samples of 15 s each). Resume-safe: progress is saved
# every 5 samples and a resubmitted task continues where it stopped. Needs the stored trees (hpc/submit_trees.sbatch). Usage, from the repository root:
#   mkdir -p logs; sbatch --array=0-199 hpc/submit_samples.sbatch
#SBATCH --job-name=space_samples
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=logs/space_samples_%A_%a.out
## #SBATCH --account=<your-project-account>
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"; source "${VENV:-venv}/bin/activate"        # your own venv: sbatch --export=ALL,VENV=/path/to/venv,...
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
N_SHARDS="${N_SHARDS:-200}"
cd scripts
for sub in cur_ash cur_en sk_ash sk_en cur_ash_rate; do
  python space_eval.py "$sub" "${SLURM_ARRAY_TASK_ID:?run with --array=0-199}" "$N_SHARDS"
done
