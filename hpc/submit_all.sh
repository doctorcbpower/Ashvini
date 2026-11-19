#!/bin/bash
# Submit the trees, then the variants after the trees finish. Run from the repository root after hpc/setup_env.sh.
set -euo pipefail
mkdir -p logs
tree_job=$(sbatch --parsable --array=0-8 hpc/submit_trees.sbatch)
echo "trees: job $tree_job"
sbatch --array=0-8 --dependency=afterok:"$tree_job" hpc/submit_variants.sbatch
