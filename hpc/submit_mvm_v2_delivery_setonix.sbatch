#!/bin/bash
# Delivery experiments (growth-factor map, reach table, seed dependence) with the revised model on the stored production
# trees. One array task per halo mass: --array=0-2 runs the mass indices 0, 4 and 12 (3e10, 3e11 and 3e13 Msun).
# Usage, from the repository root (same environment notes as submit_mvm_v2_setonix.sbatch):
#   S=/scratch/pawsey1164/cpower/Ashvini
#   sbatch --array=0-2 --export=ALL,MVM_V2_TREES_DIR=$S,MVM_V2_OUT=$S/delivery --output=$S/logs/delivery_%A_%a.out hpc/submit_mvm_v2_delivery_setonix.sbatch
#SBATCH --job-name=mvm_v2_deliv
#SBATCH --account=pawsey1164
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/delivery_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
VENV="${VENV:-$HOME/software/myenv_3.11}"
FORAOIS_ROOT="${FORAOIS_ROOT:-$HOME/software/foraois}"
source "$VENV/bin/activate"
export FORAOIS_ROOT PYTHONPATH="$PWD:${FORAOIS_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MASS_INDICES=(0 4 12)
IDX=${MASS_INDICES[${SLURM_ARRAY_TASK_ID:?run with --array=0-2}]}
TREES="${MVM_V2_TREES_DIR:?set MVM_V2_TREES_DIR to the directory holding trees_M*.npz}/trees_M$(printf %02d "$IDX").npz"
OUT="${MVM_V2_OUT:-${MYSCRATCH:-/scratch/pawsey1164/$USER}/Ashvini/delivery}"
mkdir -p "$OUT"
echo "mass index $IDX, trees $TREES, output $OUT"
python scripts/paper_figures/gen_mvm_v2_delivery.py --mass-index "$IDX" --trees-from "$TREES" --outdir "$OUT"
