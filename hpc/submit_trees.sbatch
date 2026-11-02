#!/bin/bash
#SBATCH --job-name=shmr_trees
#SBATCH --account=TODO            # TODO: your Slurm account
#SBATCH --partition=TODO          # TODO: a CPU partition
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/shmr_trees_%A_%a.out
# Usage (from the repository root): mkdir -p logs; sbatch --array=0-8 hpc/submit_trees.sbatch
# One task per halo mass index. Resume-safe.
set -euo pipefail
# module load python/3.11        # TODO: same module as in setup_env.sh
source venv/bin/activate
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FORAOIS_CONFIG_DIR="${FORAOIS_CONFIG_DIR:-$PWD/../foraois/config}"
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_trees.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
