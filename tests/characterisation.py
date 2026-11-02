"""
Characterisation cases (milestone M0 of docs/SUPERSET_DESIGN.md): a fixed set of runs of the general chain
and of the frozen MVM whose outputs are pinned in tests/golden/*.json. They were captured from the code
as it stood on 2026-11-02 (after the chemistry and units corrections) and must not change during the
refactor into Config / State / components (M1-M2) except by a deliberate, documented recapture.

`apply_case` is the only place that knows how to set parameters for a case. It uses the mechanisms that
exist today (module-level values); milestone M1 replaces its body with Config overrides. The golden
numbers do not move when that happens, which is the point.
"""
import os
import hashlib
import json

import numpy as np

from ashvini import main, utils
from ashvini import star_formation as sf
from ashvini import supernovae_feedback as sn
from ashvini import reionization as reion
from ashvini import reionization_physics as reion_physics
from ashvini import agn_feedback as agn
from ashvini import black_holes_growth as bhg
from ashvini import dust as dust_mod
from ashvini.imf import nu_SN
from ashvini.run_params import PARAMS

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "fixtures", "merger_trees_fixture.h5")
GOLDEN_DIR = os.path.join(HERE, "golden")
MASS_BINS = (1e8, 1e10, 1e11)
KEYS = ("gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass", "bh_mass", "sfr")
MVM_SHA256 = "31c701dd8f298d4b7bfc2bfb4d74f90fbe77106c817670a1432563b86d255f63"

# --- the cases -------------------------------------------------------------------------------------------
CASES = {}
for _sn in ("delayed", "instantaneous", "no"):
    for _uv in (True, False):
        for _bh in (True, False):
            CASES[f"sn-{_sn}_uv-{_uv}_bh-{_bh}"] = dict(sn_type=_sn, uv=_uv, bh=_bh)
CASES["sigma-feedback"] = dict(sigma=True, eta_agn=1.0, e_bh=0.05, edd=50.0)
CASES["sigma-growth-cap"] = dict(sigma=True, cap=True, eta_agn=1.0, e_bh=0.05, edd=50.0)
CASES["slimdisk"] = dict(growth_model="hobbs_slimdisk")
CASES["agn-delay-0.3Gyr"] = dict(agn_delay=0.3, e_bh=0.3, edd=500.0)
CASES["zero-delay"] = dict(t_d=0.0)
CASES["strong-feedback"] = dict(eps_sf=0.05, eps_p=3.0)
CASES["z_rei-10"] = dict(z_rei=10.0)
CASES["m_sn_max-100"] = dict(m_sn_max=100.0)


def apply_case(mp, case):
    """Set parameters for one case through `mp` (a pytest MonkeyPatch), so everything is restored after."""
    if "sn_type" in case: mp.setattr(main, "sn_type", case["sn_type"])
    if "t_d" in case: mp.setattr(main, "t_d", case["t_d"])
    if "uv" in case: mp.setattr(main, "UV_background", case["uv"])
    if "agn_delay" in case: mp.setattr(main, "agn_delay_time", case["agn_delay"])
    if "eps_sf" in case:
        mp.setattr(sf, "e_ff", case["eps_sf"]); mp.setattr(main, "e_ff", case["eps_sf"])
    if "eps_p" in case: mp.setattr(sn, "epsilon_p", case["eps_p"])
    if "eta_agn" in case: mp.setattr(agn, "eta_agn", case["eta_agn"])
    if "e_bh" in case: mp.setattr(bhg, "e_bh", case["e_bh"])
    if "edd" in case: mp.setattr(bhg, "eddington_multiplier", case["edd"])
    if case.get("sigma"): mp.setattr(agn, "sigma_feedback_enabled", True)
    if case.get("cap"): mp.setattr(bhg, "growth_cap_enabled", True)
    if "z_rei" in case:
        z = case["z_rei"]
        mp.setattr(reion, "z_rei", z)
        mp.setattr(reion, "beta", reion_physics.beta(z, reion.gamma))
        mp.setattr(reion, "c_omega", reion_physics.c_omega(reion.omega))
    if "m_sn_max" in case:
        nu = nu_SN(PARAMS.sn.imf, PARAMS.sn.m_sn_min, case["m_sn_max"])
        mp.setattr(main, "SN_per_Msun", nu); mp.setattr(dust_mod, "SN_per_Msun", nu)
    if case.get("bh") is False:
        for ch in (PARAMS.bh.seeding.pop3, PARAMS.bh.seeding.direct_collapse, PARAMS.bh.seeding.halo_mass_threshold):
            mp.setattr(ch, "enabled", False)


def read_bin(mass_bin):
    return utils.read_trees(FIXTURE, mass_bin)


def run_case(case):
    """Run the general chain for every fixture mass bin; returns {bin: {key: {"final": [...], "sum": x}}}."""
    out = {}
    kwargs = {"growth_model": case["growth_model"]} if "growth_model" in case else {}
    for mb in MASS_BINS:
        halo_mass, rate, z = read_bin(mb)
        res = main.run_forest(halo_mass, rate, z, **kwargs)
        out[f"{mb:g}"] = {k: {"final": res[k][:, -1].tolist(), "sum": float(np.sum(res[k]))} for k in KEYS}
    return out


def run_scalar_reference():
    halo_mass, rate, z = read_bin(1e10)
    res = main.run1_scalar(halo_mass[0], rate[0], z)
    return {k: {"final": [float(res[k][-1])], "sum": float(np.sum(res[k]))} for k in KEYS}


# --- the frozen MVM -------------------------------------------------------------------------------------
MVM_KEYS = ("bh_mass", "gas_gal", "gas_nuc", "gas_unavailable", "stars_gal", "stars_nuc", "mdot_in", "mdot_out", "A_bh")


def run_mvm():
    from ashvini.reservoir_stock import run_reservoir_stock, critical_seed_stock
    halo_mass, _, z = read_bin(1e10)
    res = run_reservoir_stock(halo_mass, z, 1.0e3)
    out = {k: {"final": res[k][:, -1].tolist(), "sum": float(np.sum(res[k]))} for k in MVM_KEYS}
    out["baryon_residual"] = {"final": res["baryon_residual"].tolist(), "sum": 0.0}
    crit, never, above = critical_seed_stock(halo_mass, z, n_iter=25)
    out["critical_seed"] = {"final": [None if np.isnan(c) else float(c) for c in crit], "sum": 0.0}
    return out


def mvm_sha256():
    path = os.path.join(os.path.dirname(HERE), "ashvini", "reservoir_stock.py")
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def regenerate():
    """Recapture all golden files from the current code. Only for a deliberate, documented change."""
    from _pytest.monkeypatch import MonkeyPatch
    general = {}
    for name, case in CASES.items():
        with MonkeyPatch.context() as mp:
            apply_case(mp, case)
            general[name] = run_case(case)
        print("captured", name)
    general["scalar-reference"] = run_scalar_reference()
    json.dump(general, open(os.path.join(GOLDEN_DIR, "general_chain.json"), "w"), indent=0)
    json.dump(run_mvm(), open(os.path.join(GOLDEN_DIR, "mvm.json"), "w"), indent=0)
    print("done")


if __name__ == "__main__":
    regenerate()
