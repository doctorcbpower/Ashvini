#!/bin/bash
# Submit the parameter-space exploration on OzSTAR with all output under /fred/oz009/cpower/Ashvini
# (trees/, space/, logs/ and, after the plot job, space_summary.txt and the figures). Run from the repository root, after hpc/setup_env.sh:
#   FORAOIS_CONFIG_DIR=/path/to/foraois/config bash hpc/submit_space_ozstar.sh
# Optional: ACCOUNT (default oz009, taken from the /fred/oz009 path: check it is your project), VENV (default venv), SHMR_OUT.
set -euo pipefail
export SHMR_OUT="${SHMR_OUT:-/fred/oz009/cpower/Ashvini}"
ACCOUNT="${ACCOUNT:-oz009}"
: "${FORAOIS_CONFIG_DIR:?set FORAOIS_CONFIG_DIR to the foraois config directory (the one holding planck2018_camb.yml)}"
export FORAOIS_CONFIG_DIR
mkdir -p "$SHMR_OUT/logs" "$SHMR_OUT/trees" "$SHMR_OUT/space"
COMMON=(--account="$ACCOUNT" --export=ALL)
LOGS="$SHMR_OUT/logs"

if [ -f "$SHMR_OUT/trees/trees_8.npz" ] && [ -f "$SHMR_OUT/trees/trees_0.npz" ]; then
  dep=()
  echo "trees found in $SHMR_OUT/trees: skipping tree generation"
else
  tree_job=$(sbatch --parsable "${COMMON[@]}" --array=0-8 --output="$LOGS/shmr_trees_%A_%a.out" hpc/submit_trees.sbatch)
  echo "trees: job $tree_job"
  dep=(--dependency=afterok:"$tree_job")
fi
sample_job=$(sbatch --parsable "${COMMON[@]}" "${dep[@]}" --array=0-199 --output="$LOGS/space_samples_%A_%a.out" hpc/submit_samples.sbatch)
echo "samples: job $sample_job (200 tasks)"
plot_job=$(sbatch --parsable "${COMMON[@]}" --dependency=afterok:"$sample_job" --output="$LOGS/plot_%j.out" hpc/submit_plot.sbatch)
echo "plot: job $plot_job (runs when all samples are done; its figures and tables go to $SHMR_OUT/space)"
