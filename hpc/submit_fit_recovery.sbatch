#!/bin/bash
# Mock recovery of the fit on the real Model (docs/FIT_AUDIT.md, section 14): task k fits case k (0-2: mock from the same trees, 3-5: from an independent
# tree realisation), 32 forked workers each. Needs the foraois configuration (the first task builds the trees and stores them in $TREE_CACHE; start
# task 0 alone first, or run `python scripts/fit_recovery.py run 0 --steps 20` once, so that the other tasks find them). Usage, from the repository root:
#   mkdir -p logs; sbatch --array=0-5 hpc/submit_fit_recovery.sbatch;  then:  python scripts/fit_recovery.py report
#SBATCH --job-name=fit_recovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/fit_recovery_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
if [ -n "${VENV:-}" ]; then source "$VENV/bin/activate"; elif [ -f venv/bin/activate ]; then source venv/bin/activate; fi
echo "python: $(command -v python) ($(python --version 2>&1))"
python - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","emcee","foraois"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment: " + ", ".join(missing) if missing else 0)
PY
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}"
export TREE_CACHE="${TREE_CACHE:-$SHMR_OUT/tree_cache}"
python scripts/fit_recovery.py run "${SLURM_ARRAY_TASK_ID:?run with --array=0-5}" --steps "${FIT_STEPS:-4000}"
