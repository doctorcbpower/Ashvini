"""
Generate the merger trees for one halo mass of the z = 0 stellar-to-halo mass relation runs, and store them.
Resume-safe (skips an existing file). One cluster task per mass index.

Usage: python scripts/shmr_z0_trees.py <mass_index 0..8 | all>
"""
import os
import subprocess
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C

warnings.filterwarnings("ignore")


def main(index):
    path = C.tree_path(index)
    if os.path.exists(path):
        print(f"mass {index}: trees exist, skipping")
        return
    from ashvini import pymctrees_adapter
    mass = C.MASSES[index]
    t0 = time.time()
    hm, hg, z, merger = pymctrees_adapter.build_forest_live(
        C.FORAOIS_CONFIG, mass, C.N_HALOES, C.Z0, C.Z_MAX, C.DZ, m_res=C.M_RES_FRACTION * mass,
        seed=C.SEED + index, algorithm=C.ALGORITHM, pk_kmax_override=C.PK_KMAX_OVERRIDE)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp.npz"
    np.savez(tmp, halo_masses=hm, halo_growth_rates=hg, redshifts=z, merger_mass=merger, mass_bin=mass,
             settings=str(dict(algorithm=C.ALGORITHM, z_max=C.Z_MAX, dz=C.DZ, m_res_fraction=C.M_RES_FRACTION,
                               pk_kmax=C.PK_KMAX_OVERRIDE, seed=C.SEED + index, n_haloes=C.N_HALOES)))
    os.replace(tmp, path)
    print(f"mass {index} ({mass:.2e} Msun): {hm.shape[0]} haloes x {hm.shape[1]} steps in {time.time() - t0:.0f} s -> {path}")


if __name__ == "__main__":
    arg = sys.argv[1]
    for i in (range(len(C.MASSES)) if arg == "all" else [int(arg)]):
        main(i)
