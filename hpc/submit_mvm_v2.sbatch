#!/bin/bash
# Revised MVM production (docs: model-verification-audit-2026-11-14.md). One array task per halo mass, 13 masses:
#   mkdir -p /fred/oz009/cpower/Ashvini/logs /fred/oz009/cpower/Ashvini/mvm_v2
#   export FORAOIS_ROOT=/path/to/foraois   (v0.1.2 or later; its git state is recorded in the output)
#   sbatch --array=0-12 hpc/submit_mvm_v2.sbatch
# then, after all tasks finish:  python scripts/paper_figures/merge_mvm_v2_production.py /fred/oz009/cpower/Ashvini/mvm_v2
# Trees are sampled once per task and stored in trees_M<idx>.npz next to the JSON; re-run an analysis on the same trees with
#   python scripts/paper_figures/gen_mvm_v2_production.py --mass-index I --trees-from <dir>/trees_MII.npz --outdir <other dir>
# Each task is single-threaded numpy; request time generously for the high masses until a first run gives the timings.
#SBATCH --job-name=mvm_v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/fred/oz009/cpower/Ashvini/logs/mvm_v2_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
VENV="${VENV:-${VIRTUAL_ENV:-}}"
if [ -z "$VENV" ] && [ -f venv/bin/activate ]; then VENV="$PWD/venv"; fi
if [ -n "$VENV" ]; then source "$VENV/bin/activate"; PY="$VENV/bin/python"; else PY="$(command -v python)"; fi
echo "python: $PY ($($PY --version 2>&1)); VENV='${VENV}'; FORAOIS_ROOT='${FORAOIS_ROOT:?set FORAOIS_ROOT}'"
$PY - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","numba","foraois"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment (%s): " % sys.executable + ", ".join(missing) if missing else 0)
PY
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OUT="${MVM_V2_OUT:-/fred/oz009/cpower/Ashvini/mvm_v2}"
$PY -m pytest -q -p no:cacheprovider tests/test_reservoir_stock_v2.py
$PY scripts/paper_figures/gen_mvm_v2_production.py --mass-index "${SLURM_ARRAY_TASK_ID:?run with --array=0-12}" --outdir "$OUT"
