#!/bin/bash
# Mock recovery of the fit on the real Model (docs/FIT_AUDIT.md, section 14): task k fits case k (0-2: mock from the same trees, 3-5: from an independent
# tree realisation), 32 forked workers each. Needs the foraois configuration. Each task first builds the two tree samples if they are not in $TREE_CACHE
# (about a minute each; tasks that start together may all build them, which only wastes that minute), then runs the fit. Do not run the fit on a login
# node. The chains are checkpointed under $SHMR_OUT/fit_recovery/checkpoints: a task that reaches the time limit is resubmitted with the same command
# and continues where it stopped. Usage, from the repository root:
#   mkdir -p /fred/oz009/cpower/Ashvini/logs; sbatch --array=0-5 hpc/submit_fit_recovery.sbatch;  then:  python scripts/fit_recovery.py report
#SBATCH --job-name=fit_recovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=/fred/oz009/cpower/Ashvini/logs/fit_recovery_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# Environment: VENV=/path/to/venv, else the venv that was active when you submitted (VIRTUAL_ENV), else ./venv. It is activated again here, because the
# batch shell may load modules after your environment arrives and put the module's python first; the venv's python is then called by full path.
VENV="${VENV:-${VIRTUAL_ENV:-}}"
if [ -z "$VENV" ] && [ -f venv/bin/activate ]; then VENV="$PWD/venv"; fi
if [ -n "$VENV" ]; then source "$VENV/bin/activate"; PY="$VENV/bin/python"; else PY="$(command -v python)"; fi
echo "python: $PY ($($PY --version 2>&1)); VENV='${VENV}'"
$PY - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","emcee","foraois"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment (%s): " % sys.executable + ", ".join(missing) if missing else 0)
PY
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export SHMR_OUT="${SHMR_OUT:-/fred/oz009/cpower/Ashvini}"
export TREE_CACHE="${TREE_CACHE:-$SHMR_OUT/tree_cache}"
$PY scripts/fit_recovery.py trees
$PY scripts/fit_recovery.py run "${SLURM_ARRAY_TASK_ID:?run with --array=0-5}" --steps "${FIT_STEPS:-4000}"
