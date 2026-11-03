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
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# Environment: VENV=/path/to/venv (a venv with bin/activate), else ./venv if it exists, else the environment active when you submitted
# (sbatch exports it with --export=ALL; this is the route for a conda environment).
if [ -n "${VENV:-}" ]; then source "$VENV/bin/activate"; elif [ -f venv/bin/activate ]; then source venv/bin/activate; fi
echo "python: $(command -v python) ($(python --version 2>&1))"
python - <<PY
import importlib.util, sys
missing = [m for m in ["numpy","scipy","astropy","h5py","yaml","matplotlib"] if importlib.util.find_spec(m) is None]
sys.exit("missing Python packages in this environment: " + ", ".join(missing) if missing else 0)
PY
export SHMR_OUT="${SHMR_OUT:-$PWD/scripts/output/shmr_z0}" MPLBACKEND=Agg
python scripts/shmr_z0_plot.py          # the z = 0 comparison of the star formation laws (needs the variant results)
python scripts/space_analyse.py         # Sobol indices, accepted samples, predicted low-mass band (needs all the parameter-space shards)
