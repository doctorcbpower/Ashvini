"""
Drive Stochastic_Occupation's smooth control and Ashvini on one halo history (the setup of Stochastic_Occupation
docs/ashvini_relation.md, parity audit): M_h(z) = 1e9 exp(-0.5 (z - 5)) Msun, z = 20 -> 5, Ashvini's clock and
baryon fraction, no star formation gate, no UV suppression, no black hole, no gas at z = 20.

Needs Stochastic_Occupation's source on the path: set STOCHOCC_SRC, or keep the repository at
~/MyCodes/Stochastic_Occupation. `available()` says whether it can be imported.
"""
import os
import sys
from pathlib import Path

import numpy as np

from ashvini import config, utils
from ashvini.constants import Z_SUN
from ashvini.config import override
from ashvini.main import run_forest

GYR = 3.15576e16
ZGRID = (15.0, 10.0, 7.0, 5.0)
FB = utils.Omega_b / utils.Omega_m
# Ashvini's e_ff / (0.141 t_H) equals stochocc's eps_ff / (0.2047 t_H) for this eps_ff (the 0.2047 is the value of
# the stochocc trajectory used in the parity audit)
EPS_FF_STOCHOCC = 0.015 * 0.2047 / 0.141


def _src():
    return Path(os.environ.get("STOCHOCC_SRC", Path.home() / "MyCodes" / "Stochastic_Occupation" / "src"))


def available():
    if not (_src() / "stochocc").is_dir():
        return False
    if str(_src()) not in sys.path:
        sys.path.insert(0, str(_src()))
    try:
        import stochocc.engine  # noqa: F401
    except Exception:
        return False
    return True


def Mh_of_z(z):
    return 1e9 * np.exp(-0.5 * (np.asarray(z) - 5.0))


def ashvini(cfg, dz=0.005, zout=ZGRID):
    z = np.arange(20.0, 5.0 - 1e-9, -dz)
    t = utils.time_at_z(z)
    Mh = Mh_of_z(z)
    Mdot = np.zeros_like(Mh)
    Mdot[:-1] = np.diff(Mh) / np.diff(t)
    Mdot[-1] = Mdot[-2]
    r = run_forest(Mh[None], Mdot[None], z, config=cfg)
    idx = [int(np.argmin(abs(z - zz))) for zz in zout]
    return {k: r[k][0, idx] for k in ("gas_mass", "stars_mass", "gas_metals", "stars_metals", "dust_mass")}, r, z


def base_config(**changes):
    """Ashvini with the features the parity setup switches off: no UV suppression, no black hole seeding."""
    c = {"reion.UVB_enabled": False,
         "bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
         "bh.seeding.halo_mass_threshold.enabled": False,
         "metals.Z_yield": 0.06, "dust.dust_yield": 0.004}
    c.update(changes)
    return override(config.DEFAULT, **c)


def stochocc(zout=ZGRID, eps_fb=0.0, chem=False, z_igm=1.0e-3 * Z_SUN, Mh=None, sampler="smooth", seed=0, **extra):
    from stochocc import halo as H
    from stochocc.engine import Params, evolve_branch
    from stochocc.imf import WeightedQuadratureSampler
    cosmo = H.Cosmology(H0_kms_mpc=67.66, OmegaM=utils.Omega_m, OmegaLambda=1 - utils.Omega_m, h=0.6766, fb=FB)

    class Traj:
        Delta = 200.0
        def __init__(self, z0=20.0, z1=5.0):
            self.cosmo = cosmo
            self.t_min = float(utils.time_at_z(z0)) * GYR
            self.t_max = float(utils.time_at_z(z1)) * GYR
        def z_at_t(self, t): return float(utils.z_at_time(t / GYR))
        def mass_msun_at_t(self, t): return float((Mh or Mh_of_z)(self.z_at_t(t)))
        def v_esc_at_t(self, t): return float(H.v_esc_cms(self.mass_msun_at_t(t), self.z_at_t(t), cosmo, 200.0))
        def T_vir_at_t(self, t): return float(H.T_vir(self.mass_msun_at_t(t), self.z_at_t(t), cosmo, 200.0))
        def t_dyn_at_t(self, t): return float(0.2047 / (utils.cosmo.H(self.z_at_t(t)).to("1/s").value))

    tr = Traj()
    p = Params(f_b=FB, eps_ff=EPS_FF_STOCHOCC, eps_fb=eps_fb, time_scheme="v2", sf_dt_max_myr=1.0,
               lifetime_Z=0.002, chemistry=chem, Z_IGM=z_igm, Y_Z=0.06, Y_d=0.004, **{"T_vir_thresh": 0.0, **extra})
    out = []
    for zz in zout:
        smp = (WeightedQuadratureSampler(p.nu_SN, 200, 1.0, p.imf_slope, p.lifetime_scale, p.lifetime_Z)
               if sampler == "smooth" else None)
        t1 = float(utils.time_at_z(zz)) * GYR
        out.append(evolve_branch(tr, p, np.random.default_rng(seed), tr.t_min, t1, 0.2 * (t1 - tr.t_min),
                                 M_g0=0.0, record_history=False, sampler=smp)[0])
    return out


# --- pinned Stochastic_Occupation numbers -------------------------------------------------------------------------
CASES = {
    # feedback off, chemistry on (A1, B1)
    "no_feedback": dict(eps_fb=0.0, chem=True),
    # momentum closure, f_p = 0.7625 (E1); Stochastic_Occupation does not define the momentum closure with chemistry
    "momentum": dict(eps_fb=0.0, chem=False, closure="momentum", f_p=0.7625),
}
GOLDEN = Path(__file__).parent / "golden" / "stochocc_smooth.json"


def _run_one(args):
    case, zz = args
    st = stochocc(zout=(zz,), **CASES[case])[0]
    return case, zz, {k: getattr(st, k) for k in ("M_g", "M_star", "M_Zg", "M_Zs", "M_d", "M_ej")}


def regenerate(processes=8):
    """Run Stochastic_Occupation's smooth control for every case and redshift and write the pins (minutes)."""
    import json
    from multiprocessing import Pool
    assert available(), "Stochastic_Occupation not found: set STOCHOCC_SRC"
    jobs = [(c, zz) for c in CASES for zz in ZGRID]
    with Pool(processes) as pool:
        results = pool.map(_run_one, jobs)
    out = {c: {} for c in CASES}
    for case, zz, vals in results:
        out[case][str(zz)] = vals
    out["_meta"] = {"setup": "see tests/stochocc_parity.py", "stochocc_source": str(_src()),
                    "sampler": "WeightedQuadratureSampler(K=200, dt_emit=1 Myr)", "lifetime_Z": 0.002}
    GOLDEN.write_text(json.dumps(out, indent=1))
    return out


# --- the optional event-stepper features (milestone M5): a minihalo growing from 4e5 to 2e7 Msun between z = 20 and 10 ---
MINIHALO = lambda z: 2.0e7 * np.exp(-0.4 * (np.asarray(z) - 10.0))
OPTION_BASE = {"ash": {"sn.kernel": "imf", "sn.m_sn_max": 100.0, "sn.lifetime_Z": 0.002, "sf.efficiency": 0.015,
                       "sn.closure": "energy", "sn.eps_fb": 0.1},
               "sto": {"eps_fb": 0.1}}
# case name -> the Ashvini Config changes and the Stochastic_Occupation Params that mean the same thing
OPTION_CASES = {
    "inflow_lag": {"ash": {"events.inflow_kappa": 3.0}, "sto": {"inflow_kappa": 3.0}},
    "recycling": {"ash": {"events.recycling_fraction": 0.5, "events.recycling_kappa": 1.0},
                  "sto": {"rec_frac": 0.5, "rec_kappa": 1.0}},
    # strong coupling, so that explosions regularly find less gas than their energy could unbind; the control has no bank
    "energy_bank_control": {"ash": {"sn.eps_fb": 1.0}, "sto": {"eps_fb": 1.0}},
    "energy_bank": {"ash": {"sn.eps_fb": 1.0, "events.bank_tau_myr": 100.0}, "sto": {"eps_fb": 1.0, "bank_tau_myr": 100.0}},
    "access_fraction": {"ash": {"sn.eps_fb": 1.0, "events.access_fraction": 0.05}, "sto": {"eps_fb": 1.0, "f_access": 0.05}},
    "gate": {"ash": {"events.gate_T_vir_K": 5.0e3}, "sto": {"T_vir_thresh": 5.0e3}},
}
OPTION_GOLDEN = Path(__file__).parent / "golden" / "stochocc_options.json"


def ashvini_option(case, sampler="smooth", seed=None):
    """Event stepper on the minihalo with the case's options; returns the final-state masses."""
    from ashvini.main import run_events
    z = np.arange(20.0, 10.0 - 1e-9, -0.1)
    M = MINIHALO(z)[None, :]
    cfg = base_config(**{**OPTION_BASE["ash"], **(OPTION_CASES[case]["ash"] if case else {})})
    # K = 200 representative masses, as in the pinned Stochastic_Occupation control: a cap on what one explosion removes
    # acts per explosion, so the answer depends on how the supernovae are split
    kw = {"sampler_kwargs": {"K": 200}} if sampler == "smooth" else {}
    r = run_events(M, np.zeros_like(M), z, config=cfg, sampler=sampler, seed=seed, **kw)
    return {"M_star": r["stars_mass"][0, -1], "M_g": r["gas_mass"][0, -1], "M_ej": r["ejected_mass"][0, -1],
            "M_ret": r["returned_mass"][0, -1]}


def stochocc_option(case, sampler="smooth", seed=0):
    st = stochocc(zout=(10.0,), Mh=MINIHALO, sampler=sampler, seed=seed,
                  **{**OPTION_BASE["sto"], **(OPTION_CASES[case]["sto"] if case else {})})[0]
    return {"M_star": st.M_star, "M_g": st.M_g, "M_ej": st.M_ej, "M_ret": st.M_ret}


def regenerate_options(cases=None, processes=8):
    """Run Stochastic_Occupation's smooth control for the option cases (and the plain minihalo, key 'plain') and update the pins."""
    import json
    from multiprocessing import Pool
    assert available(), "Stochastic_Occupation not found: set STOCHOCC_SRC"
    names = ["plain"] + list(cases or OPTION_CASES)
    with Pool(min(processes, len(names))) as pool:
        results = pool.map(_option_one, names)
    out = json.loads(OPTION_GOLDEN.read_text()) if OPTION_GOLDEN.exists() else {}
    out.update(dict(results))
    OPTION_GOLDEN.write_text(json.dumps(out, indent=1))
    return out


def _option_one(name):
    return name, stochocc_option(None if name == "plain" else name)


if __name__ == "__main__":
    import sys
    if "options" in sys.argv:
        regenerate_options([a for a in sys.argv[1:] if a != "options"] or None)
    else:
        regenerate()
