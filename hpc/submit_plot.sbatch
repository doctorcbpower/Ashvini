#!/bin/bash
# Make the figures and summary tables from the stored results (light: one core, a few minutes). Usage:
#   sbatch hpc/submit_plot.sbatch   (or run the lines below on the login node)
#SBATCH --job-name=shmr_plot
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/plot_%j.out
## #SBATCH --account=<your-project-account>
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"; source "${VENV:-venv}/bin/activate"        # your own venv: sbatch --export=ALL,VENV=/path/to/venv,...
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}" MPLBACKEND=Agg
python scripts/shmr_z0_plot.py          # the z = 0 comparison of the star formation laws (needs the variant results)
python scripts/space_analyse.py         # Sobol indices, accepted samples, predicted low-mass band (needs all the parameter-space shards)
