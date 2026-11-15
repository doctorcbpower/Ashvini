#!/bin/bash
# Sensitivity ledger of the revised model (one-at-a-time variation of every parameter and assumption, cooling gates
# included) on the stored PCH08 trees. One array task per halo mass: --array=0-3 runs mass indices 0, 4, 8, 12
# (3e10, 3e11, 3e12, 3e13 Msun). Same environment notes as submit_mvm_v2_setonix.sbatch.
#   S=/scratch/pawsey1164/cpower/Ashvini
#   sbatch --array=0-3 --export=ALL,MVM_V2_TREES_DIR=$S,MVM_V2_OUT=$S/ledger --output=$S/logs/ledger_%A_%a.out hpc/submit_mvm_v2_ledger_setonix.sbatch
# Optional: MVM_V2_ARGS="--n-trees 60 --n-sub 8" for the final run (about 4x the cost of the defaults below).
#SBATCH --job-name=mvm_v2_ledger
#SBATCH --account=pawsey1164
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ledger_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
VENV="${VENV:-$HOME/software/myenv_3.11}"
FORAOIS_ROOT="${FORAOIS_ROOT:-$HOME/software/foraois}"
source "$VENV/bin/activate"
export FORAOIS_ROOT PYTHONPATH="$PWD:${FORAOIS_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MASS_INDICES=(0 4 8 12)
IDX=${MASS_INDICES[${SLURM_ARRAY_TASK_ID:?run with --array=0-3}]}
TREES="${MVM_V2_TREES_DIR:?set MVM_V2_TREES_DIR to the directory holding trees_M*.npz}/trees_M$(printf %02d "$IDX").npz"
OUT="${MVM_V2_OUT:-${MYSCRATCH:-/scratch/pawsey1164/$USER}/Ashvini/ledger}"
mkdir -p "$OUT"
echo "mass index $IDX, trees $TREES, output $OUT"
python scripts/paper_figures/gen_mvm_v2_ledger.py --mass-index "$IDX" --trees-from "$TREES" --outdir "$OUT" \
    --n-trees 40 --n-sub 4 --n-iter 30 ${MVM_V2_ARGS:-}
