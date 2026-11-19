"""
Parameter-space exploration of the z = 0 stellar-to-halo mass relation (docs/PARAMETER_SPACE_AUDIT.md).

Sub-spaces: a star formation law (the current free-fall law, or the Schmidt-Kennicutt disc law) times a supernova wind model
(Ashvini's mass loading, or the energy-driven wind), all with the equilibrium UV suppression, plus the current law with Ashvini's
mass loading and the rate-dependent UV suppression for comparison. Each sub-space is sampled with the Saltelli scheme
(A, B and one matrix AB_i per parameter), N * (D + 2) samples, from which Sobol first-order and total indices follow.

The design is regenerated deterministically from the sub-space name and N, so nothing but the evaluated outputs is stored.
"""
import numpy as np
from scipy.stats import qmc

import shmr_z0_common as C

# name -> (low, high, log?, the Config path it sets)
PARAMETERS = {
    "eps": (0.003, 0.1, True, "sf.efficiency"),
    "A": (8e-5, 8e-4, True, "disc.schmidt_kennicutt_normalisation"),
    "N": (1.0, 1.6, False, "disc.schmidt_kennicutt_slope"),
    "chi": (0.03, 0.12, False, "disc.angular_momentum_ratio"),      # R_d / R_vir; converted by C.disc_size
    "epsilon_p": (1.0, 25.0, True, "sn.epsilon_p"),
    "eps_fb": (0.01, 1.0, True, "sn.eps_fb"),
    "delay": (0.003, 0.05, False, "sn.delay_time"),
    "z_rei": (6.0, 9.0, False, "reion.z_reion"),
    "gamma": (10.0, 20.0, False, "reion.gamma"),
    "omega": (1.0, 3.0, False, "reion.omega"),
    "bh_eff": (3e-4, 1e-2, True, "bh.efficiency"),
    "eta_agn": (0.05, 2.0, True, "bh.eta_agn"),
}
COMMON = ["delay", "z_rei", "gamma", "omega", "bh_eff", "eta_agn"]
SK_BASE = {"disc.enabled": True, "disc.star_formation_law": "schmidt_kennicutt"}

# sub-space -> (fixed Config overrides, the parameters that vary)
SUBSPACES = {
    "cur_ash": ({"reion.form": "equilibrium"}, ["eps", "epsilon_p"] + COMMON),
    "cur_en": ({"reion.form": "equilibrium", "sn.wind_model": "energy"}, ["eps", "eps_fb"] + COMMON),
    "sk_ash": ({"reion.form": "equilibrium", **SK_BASE}, ["A", "N", "chi", "epsilon_p"] + COMMON),
    "sk_en": ({"reion.form": "equilibrium", "sn.wind_model": "energy", **SK_BASE}, ["A", "N", "chi", "eps_fb"] + COMMON),
    "cur_ash_rate": ({"reion.form": "rate_dependent"}, ["eps", "epsilon_p"] + COMMON),
}
N_BASE = 128
# fields stored per sample and mass
FIELDS = ("stars_median", "stars_p16", "stars_p84", "gas_median", "gas_Z_median", "stars_Z_median", "bh_median", "frac_with_stars", "n_nonfinite")


def n_samples(subspace, N=N_BASE):
    return N * (len(SUBSPACES[subspace][1]) + 2)


def design(subspace, N=N_BASE):
    """Unit-cube samples, shape (N * (D + 2), D): the matrices A, B, AB_1 ... AB_D stacked (Saltelli; first-order and total indices)."""
    D = len(SUBSPACES[subspace][1])
    seed = 1000 + list(SUBSPACES).index(subspace)
    u = qmc.Sobol(d=2 * D, scramble=True, seed=seed).random(N)
    A, B = u[:, :D], u[:, D:]
    blocks = [A, B]
    for i in range(D):
        AB = A.copy()
        AB[:, i] = B[:, i]
        blocks.append(AB)
    return np.vstack(blocks)


def physical(subspace, unit_row):
    """{parameter name: value} for one row of the unit cube."""
    out = {}
    for name, u in zip(SUBSPACES[subspace][1], unit_row):
        lo, hi, log, _ = PARAMETERS[name]
        out[name] = float(np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))) if log else lo + u * (hi - lo))
    return out


def config_changes(subspace, values):
    """Config overrides for one sample of a sub-space."""
    ch = dict(SUBSPACES[subspace][0])
    for name, v in values.items():
        path = PARAMETERS[name][3]
        ch[path] = C.disc_size(v) if name == "chi" else v
    return ch


# --- Sobol indices from the Saltelli matrices ---

def sobol_indices(f, N, D):
    """
    First-order (Saltelli 2010) and total (Jansen 1999) Sobol indices of a model output f, with the samples stacked as in
    `design`: f[:N] for A, f[N:2N] for B, then D blocks of N for AB_i. Returns (S1, ST), each of length D.
    """
    f = np.asarray(f, dtype=float)
    fA, fB = f[:N], f[N:2 * N]
    var = np.var(np.concatenate([fA, fB]), ddof=1)
    S1, ST = np.empty(D), np.empty(D)
    for i in range(D):
        fAB = f[(2 + i) * N:(3 + i) * N]
        S1[i] = np.mean(fB * (fAB - fA)) / var
        ST[i] = 0.5 * np.mean((fA - fAB) ** 2) / var
    return S1, ST


def bootstrap_indices(f, N, D, n_boot=200, seed=0):
    """Bootstrap (resampling the N base rows) 95 percent half-widths of S1 and ST."""
    f = np.asarray(f, dtype=float)
    rng = np.random.default_rng(seed)
    blocks = f.reshape(D + 2, N)
    s1s, sts = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, N, N)
        s1, st = sobol_indices(blocks[:, idx].reshape(-1), N, D)
        s1s.append(s1); sts.append(st)
    s1s, sts = np.array(s1s), np.array(sts)
    return 0.5 * (np.percentile(s1s, 97.5, axis=0) - np.percentile(s1s, 2.5, axis=0)), \
        0.5 * (np.percentile(sts, 97.5, axis=0) - np.percentile(sts, 2.5, axis=0))
