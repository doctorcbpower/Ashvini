#!/bin/bash
# One-off environment setup on the cluster login node. Usage: bash hpc/setup_env.sh <path-to-foraois>
# EDIT the module lines for your cluster (check `module avail python`). Versions used locally: Python 3.11, numpy 2, scipy, numba, camb, PyYAML, h5py, astropy.
set -euo pipefail
FORAOIS_DIR="${1:?usage: bash hpc/setup_env.sh <path-to-foraois>}"
# module purge
# module load python/3.11      # TODO: replace with the cluster's Python >= 3.11 module
python3 -m venv "$PWD/venv"
source "$PWD/venv/bin/activate"
pip install --upgrade pip
pip install numpy scipy numba matplotlib pyyaml camb h5py astropy pytest emcee corner
pip install -e "$FORAOIS_DIR"
pip install -e "$PWD"
python - <<'PY'
import numpy, scipy, numba, camb, foraois, ashvini
print("numpy", numpy.__version__, "scipy", scipy.__version__, "numba", numba.__version__, "camb", camb.__version__)
PY
echo "environment ready: source $PWD/venv/bin/activate"
