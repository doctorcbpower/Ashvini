#!/bin/bash
# Generate the merger trees, one task per halo mass (light: one core, about a minute each). Usage, from the repository root:
#   mkdir -p logs; sbatch --array=0-8 hpc/submit_trees.sbatch
#SBATCH --job-name=shmr_trees
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/shmr_trees_%A_%a.out
## #SBATCH --account=<your-project-account>
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"; source "${VENV:-venv}/bin/activate"        # your own venv: sbatch --export=ALL,VENV=/path/to/venv,...
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FORAOIS_CONFIG_DIR="${FORAOIS_CONFIG_DIR:-$PWD/../foraois/config}"
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_trees.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
