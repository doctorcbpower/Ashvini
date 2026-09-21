import sys
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); from foraois_paths import FORAOIS_SRC, FORAOIS_CONFIG; sys.path.insert(0, FORAOIS_SRC)
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree, PCHMergerTree
from foraois.utils import io
from ashvini import pymctrees_adapter
rp = io.get_params(FORAOIS_CONFIG); h = rp["Cosmology"]["h"]
cd0 = cosmo_utils.CosmoData(rp, redshift=[0.0]); cd5 = cosmo_utils.CosmoData(rp, redshift=[5.0])
def st(gen, M0, z0, zt, N=60, dz=0.05, mres=1e4, seed=3):
    hm, zz, *_ = pymctrees_adapter.build_forest_for_bin(gen, M0, h, N, z0=z0, z_max=zt + 1.0, m_res_msun=mres, dz=dz, backend="numba", rng_seed=seed)
    m = hm[:, int(np.argmin(abs(zz - zt)))] / M0
    return f"{np.median(m):.3g} [{np.percentile(m,16):.2g},{np.percentile(m,84):.2g}]"
print("A. standard anchor z0=0, M0=1e12 (P(k) at z=0): M(z)/M0 median [16,84]")
for zt in (1.0, 2.0):
    print(f"   z={zt}:  ZH {st(ZhangHuiMergerTree(cd0, rp, model='cdm'), 1e12, 0.0, zt)}   PCH08 {st(PCHMergerTree(cd0, rp), 1e12, 0.0, zt)}")
print("B. anchor z0=5 (P(k) at z=5), M(z=10)/M0 for several halo masses")
for M0 in (3e10, 3e11, 3e12, 3e13):
    print(f"   M0={M0:.0e}: ZH {st(ZhangHuiMergerTree(cd5, rp, model='cdm'), M0, 5.0, 10.0)}   PCH08 {st(PCHMergerTree(cd5, rp), M0, 5.0, 10.0)}")
print("C. PCH08, M0=3e10, z0=5: dependence on dz and M_res")
for dz in (0.02, 0.05, 0.1, 0.2):
    print(f"   dz={dz}: {st(PCHMergerTree(cd5, rp), 3e10, 5.0, 10.0, dz=dz)}")
for mres in (1e3, 1e4, 1e6):
    print(f"   M_res={mres:.0e}: {st(PCHMergerTree(cd5, rp), 3e10, 5.0, 10.0, mres=mres)}")
