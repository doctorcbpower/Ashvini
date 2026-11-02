"""
Shared settings for the z = 0 stellar-to-halo mass relation runs that compare star formation laws
(scripts/shmr_z0_trees.py, shmr_z0_variants.py, shmr_z0_plot.py; hpc/ for the cluster).

Environment variables (all optional):
  SHMR_OUT            output directory (default scripts/output/shmr_z0)
  FORAOIS_CONFIG_DIR  directory of the foraois cluster of planck2018_camb.yml (default ~/MyCodes/foraois/config)
  SHMR_N_HALOES       haloes per mass (default 100)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

OUT = os.environ.get("SHMR_OUT", os.path.join(HERE, "output", "shmr_z0"))
FORAOIS_CONFIG_DIR = os.environ.get("FORAOIS_CONFIG_DIR", os.path.expanduser("~/MyCodes/foraois/config"))
FORAOIS_CONFIG = os.path.join(FORAOIS_CONFIG_DIR, "planck2018_camb.yml")
N_HALOES = int(os.environ.get("SHMR_N_HALOES", "100"))

# halo masses at z = 0 (Msun)
MASSES = np.logspace(8, 12, 9)
# tree settings: Zhang-Hui trees, finest verified settings of docs/RESOLUTION_CONVERGENCE.md and shmr_dm_model_comparison.py
Z0, Z_MAX, DZ = 0.0, 30.0, 0.005
M_RES_FRACTION = 1e-4
PK_KMAX_OVERRIDE = 3000.0
ALGORITHM = "zh"
SEED = 42

TINY = 0.0247  # fiducial R_d / R_vir = 0.035 / sqrt(2); see components/disc.py


def disc_size(chi):
    """disc.angular_momentum_ratio that gives R_d / R_vir = chi at the fixed spin 0.035."""
    return chi * np.sqrt(2.0) / 0.035


SK = {"disc.enabled": True, "disc.star_formation_law": "schmidt_kennicutt"}
VARIANTS = {
    "current_eps0.015": {},
    "current_eps0.005": {"sf.efficiency": 0.005},
    "current_eps0.05": {"sf.efficiency": 0.05},
    "sk_chi0.025": dict(SK),
    "sk_chi0.06": {**SK, "disc.angular_momentum_ratio": disc_size(0.06)},
    "sk_chi0.06_small": {**SK, "disc.angular_momentum_ratio": disc_size(0.06 * 10 ** -0.25)},
    "sk_chi0.06_large": {**SK, "disc.angular_momentum_ratio": disc_size(0.06 * 10 ** 0.25)},
    "sk_chi0.06_slope1": {**SK, "disc.angular_momentum_ratio": disc_size(0.06), "disc.schmidt_kennicutt_slope": 1.0},
    "current_zrei6": {"reion.z_reion": 6.0},
    "current_zrei8": {"reion.z_reion": 8.0},
    "current_noUV": {"reion.UVB_enabled": False},
    "sk_chi0.06_noUV": {**SK, "disc.angular_momentum_ratio": disc_size(0.06), "reion.UVB_enabled": False},
}


def tree_path(index):
    return os.path.join(OUT, "trees", f"trees_{index}.npz")


def result_path(variant, index):
    return os.path.join(OUT, "results", f"{variant}_{index}.npz")
