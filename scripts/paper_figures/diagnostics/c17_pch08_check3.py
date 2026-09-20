import sys
sys.path.insert(0, "/Users/00075868/MyCodes/foraois/src")
import numpy as np
from foraois import cosmo_utils, ZhangHuiMergerTree, PCHMergerTree
from foraois.utils import io
rp = io.get_params("/Users/00075868/MyCodes/foraois/config/menon_power_2024.yml"); h = rp["Cosmology"]["h"]
cd0 = cosmo_utils.CosmoData(rp, redshift=[0.0])
for name, gen in (("PCH08", PCHMergerTree(cd0, rp)), ("ZH", ZhangHuiMergerTree(cd0, rp, model="cdm"))):
    for M_res in (1e4, 1e8, 1e10):
        out = gen.build_forest_numba(M0_array=np.full(40, 1e12 * h), z0=0.0, z_max=3.0, M_res=M_res * h, dz=0.05)
        mh, zs, sm, mg = out[:4] if len(out) == 4 else (out[0], out[1], out[3], out[4])
        m1 = mh[:, np.argmin(abs(zs - 1.0)) - 1] / (1e12 * h)
        print(f"{name:5s} M_res={M_res:.0e}: fraction of steps with a resolved merger {np.mean(mg > 0)*100:6.2f}%;  total merged mass / total smooth mass {mg.sum()/max(sm.sum(),1e-300):.3g};  M(z=1)/M0 median {np.median(m1):.3f} [{np.percentile(m1,16):.3f},{np.percentile(m1,84):.3f}]")
