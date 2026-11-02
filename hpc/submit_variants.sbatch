#!/bin/bash
#SBATCH --job-name=shmr_variants
#SBATCH --account=TODO            # TODO: your Slurm account
#SBATCH --partition=TODO          # TODO: a CPU partition
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/shmr_variants_%A_%a.out
# Usage (from the repository root): sbatch --array=0-8 [--dependency=afterok:<trees jobid>] hpc/submit_variants.sbatch
# One task per halo mass index, running every variant. Resume-safe.
set -euo pipefail
# module load python/3.11        # TODO: same module as in setup_env.sh
source venv/bin/activate
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_variants.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
