#!/bin/bash
# Run every variant of the z = 0 star formation comparison, one task per halo mass (one core, a few minutes each). Resume-safe. Usage:
#   sbatch --array=0-8 [--dependency=afterok:<trees jobid>] hpc/submit_variants.sbatch
#SBATCH --job-name=shmr_variants
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/shmr_variants_%A_%a.out
## #SBATCH --account=<your-project-account>
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"; source "${VENV:-venv}/bin/activate"        # your own venv: sbatch --export=ALL,VENV=/path/to/venv,...
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_variants.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
