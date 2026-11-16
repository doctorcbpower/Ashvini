#!/bin/bash
# Production PCH08 trees for all 13 mass bins at a resolution R = M0 / m_res chosen from the convergence run (floor_steps must be 0).
#   S=/scratch/pawsey1164/cpower/Ashvini
#   sbatch --array=0-12 --export=ALL,RATIO=3e9,OUT=$S/trees_R3e9 --output=$S/logs/trees_%A_%a.out hpc/submit_pch08_trees_production_setonix.sbatch
# Output trees_M<idx>_R<..>.npz are read by gen_mvm_v2_ledger.py / gen_mvm_v2_production.py through --trees-from (same layout).
#SBATCH --job-name=pch08_trees
#SBATCH --account=pawsey1164
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/trees_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
VENV="${VENV:-$HOME/software/myenv_3.11}"
FORAOIS_ROOT="${FORAOIS_ROOT:-$HOME/software/foraois}"
source "$VENV/bin/activate"
export FORAOIS_ROOT PYTHONPATH="$PWD:${FORAOIS_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
IDX=${SLURM_ARRAY_TASK_ID:?run with --array=0-12}
RATIO="${RATIO:?set RATIO, e.g. 3e9}"
OUT="${OUT:-${MYSCRATCH:-/scratch/pawsey1164/$USER}/Ashvini/trees}"
mkdir -p "$OUT" logs
python scripts/paper_figures/gen_pch08_trees_converged.py --mass-index "$IDX" --ratio "$RATIO" --n-trees "${N_TREES:-240}" --outdir "$OUT"
