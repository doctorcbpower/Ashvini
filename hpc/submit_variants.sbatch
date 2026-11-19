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
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# Environment: VENV=/path/to/venv (a venv with bin/activate), else ./venv if it exists, else the environment active when you submitted
# (sbatch exports it with --export=ALL; this is the route for a conda environment).
if [ -n "${VENV:-}" ]; then source "$VENV/bin/activate"; elif [ -f venv/bin/activate ]; then source venv/bin/activate; fi
echo "python: $(command -v python) ($(python --version 2>&1))"
python - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment: " + ", ".join(missing) if missing else 0)
PY
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_variants.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
