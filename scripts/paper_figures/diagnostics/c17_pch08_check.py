import sys
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree, PCHMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
def stats(gen, label, M0=3e10, z0=5.0, backend="numba", N=60, seed=3):
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(gen, M0, h, N, z0=z0, z_max=25.0, m_res_msun=1e4, dz=0.05, backend=backend, rng_seed=seed)
    at = lambda z: hm[:, int(np.argmin(abs(zz - z)))]
    m10, m15 = at(10.0), at(15.0)
    print(f"  {label:44s} M(z=10)/M0 median {np.median(m10)/M0:.3g} [16,84 = {np.percentile(m10,16)/M0:.2g}, {np.percentile(m10,84)/M0:.2g}]  M(z=15)/M0 {np.median(m15)/M0:.3g}  identical trees? {bool(np.all(hm == hm[0]))}")
for zpk in (5.0, 0.0):
    cd = cosmo_utils.CosmoData(rp, redshift=[zpk])
    print(f"CosmoData(redshift=[{zpk}]), trees anchored at z0=5, M0=3e10")
    stats(ZhangHuiMergerTree(cd, rp, model="cdm"), "ZH  numba")
    stats(PCHMergerTree(cd, rp), "PCH08 numba", backend="numba")
    stats(PCHMergerTree(cd, rp), "PCH08 numpy", backend="numpy", N=30)
print("Correa+15-type expectation for a 3e10 halo at z=5: M(z=10)/M0 of order 0.05-0.2 (the ZH value is the fiducial)")
