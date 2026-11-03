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
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# Environment: VENV=/path/to/venv (a venv with bin/activate), else ./venv if it exists, else the environment active when you submitted
# (sbatch exports it with --export=ALL; this is the route for a conda environment).
if [ -n "${VENV:-}" ]; then source "$VENV/bin/activate"; elif [ -f venv/bin/activate ]; then source venv/bin/activate; fi
echo "python: $(command -v python) ($(python --version 2>&1))"
python - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","numba","camb","foraois"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment: " + ", ".join(missing) if missing else 0)
PY
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FORAOIS_CONFIG_DIR="${FORAOIS_CONFIG_DIR:-$PWD/../foraois/config}"
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
python scripts/shmr_z0_trees.py "${SLURM_ARRAY_TASK_ID:?run with --array=0-8}"
