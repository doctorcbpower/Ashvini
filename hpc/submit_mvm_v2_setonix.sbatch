#!/bin/bash
# Revised MVM production on Setonix (Pawsey), project pawsey1164. One array task per halo mass (13 masses, shared-node
# 'work' partition, 4 cores each). See docs/MVM_V2.md. DRAFT, not yet run on Setonix: the first array task is the test.
#
# Usage, from the repository root:
#   mkdir -p logs
#   export FORAOIS_ROOT=$MYSOFTWARE/foraois        # v0.1.2 or later; use your checkout
#   sbatch --array=0-12 hpc/submit_mvm_v2_setonix.sbatch
# Defaults: VENV=~/software/myenv_3.11, data in $MYSCRATCH/Ashvini/mvm_v2 (override with MVM_V2_OUT), logs in ./logs.
# Setonix scratch is purged, so copy trees_M*.npz and the JSON files to permanent storage once the run is finished.
# Then:  python scripts/paper_figures/merge_mvm_v2_production.py $MYSCRATCH/Ashvini/mvm_v2
#SBATCH --job-name=mvm_v2
#SBATCH --account=pawsey1164
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mvm_v2_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# If the venv was built on a Pawsey python module, load that module here first, for example: module load python/3.11.6
VENV="${VENV:-$HOME/software/myenv_3.11}"
source "$VENV/bin/activate"
PY="$VENV/bin/python"
echo "python: $PY ($($PY --version 2>&1)); FORAOIS_ROOT='${FORAOIS_ROOT:?set FORAOIS_ROOT}'; task ${SLURM_ARRAY_TASK_ID:?run with --array=0-12}"
$PY - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","numba","foraois","pytest"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment (%s): " % sys.executable + ", ".join(missing) if missing else 0)
PY
export PYTHONPATH="$PWD:${FORAOIS_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_cache_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
OUT="${MVM_V2_OUT:-${MYSCRATCH:?set MYSCRATCH or MVM_V2_OUT}/Ashvini/mvm_v2}"
mkdir -p "$OUT"
$PY -m pytest -q -p no:cacheprovider tests/test_reservoir_stock_v2.py
$PY scripts/paper_figures/gen_mvm_v2_production.py --mass-index "$SLURM_ARRAY_TASK_ID" --outdir "$OUT"
