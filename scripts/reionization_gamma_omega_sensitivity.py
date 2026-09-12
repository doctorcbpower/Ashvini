#!/usr/bin/env python3
"""
Sensitivity of the z=0 stellar-to-halo mass relation (SHMR) to the
Okamoto, Gao & Theuns (2008) reionization-suppression model's gamma/omega
shape parameters (ashvini/reionization.py), given that turning UVB off
entirely (see stellar_to_halo_mass_relation.py's on/off comparison) removed
essentially all of the low-mass suppression and the sharp M_halo~1e10-1e11
kink in the default SHMR -- this checks whether a softer (not just
on/off) reionization-suppression shape brings the low-mass end closer to
the Behroozi+2013/Moster+2013 literature curves instead of just removing
the effect outright.

gamma sets how sharply M_c(z) (the characteristic suppression mass)
transitions across z_reion; omega sets how sharply the suppression itself
turns on in halo mass at fixed z (the s(mu, omega) step function in
uv_suppression). Both, plus the beta/c_omega constants derived from them,
are captured as ashvini.reionization module-level globals at import time
(not re-read from PARAMS per call), so they must be patched directly on
that module, not via PARAMS.reion.gamma/omega -- same pattern as
main.UV_background/main.sn_type.

Defaults from run_params.yaml: gamma=15, omega=2.
"""
import time as _time

import numpy as np
import matplotlib.pyplot as plt

from ashvini import main, pymctrees_adapter, reionization as reion
from ashvini.run_params import PARAMS

from stellar_to_halo_mass_relation import moster2013_mstar, behroozi2013_mstar

PYMCTREES_CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 14, 13)
N_HALOS = 1000
Z0, Z_MAX, DZ = 0.0, 15.0, 0.05
M_RES = None
SEED = 42
BACKEND = "numba"

Z_REI = PARAMS.reion.z_reion  # fixed at the run_params.yaml value throughout; only gamma/omega vary here


def set_reion_shape(gamma, omega):
    """Patch ashvini.reionization's module-level gamma/omega (and the
    beta/c_omega constants derived from them) in place, returning the
    previous values for restoration. See module docstring for why PARAMS
    alone isn't enough."""
    old = (reion.gamma, reion.omega, reion.beta, reion.c_omega)
    reion.gamma = gamma
    reion.omega = omega
    reion.beta = Z_REI * ((np.log(1.82 * (10**3) * np.exp(-0.63 * Z_REI) - 1)) ** (-1 / gamma))
    reion.c_omega = 2 ** (omega / 3) - 1
    return old


def restore_reion_shape(old):
    reion.gamma, reion.omega, reion.beta, reion.c_omega = old


def run_ensemble(gamma, omega):
    old = set_reion_shape(gamma, omega)
    per_bin_median = []
    try:
        for mass_bin in MASS_BINS:
            t0 = _time.time()
            halo_masses, halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
                pymctrees_config_path=PYMCTREES_CONFIG, mass_bin=mass_bin,
                n_halos=N_HALOS, z0=Z0, z_max=Z_MAX, dz=DZ,
                m_res=M_RES, backend=BACKEND, seed=SEED,
            )
            result = main.run_forest(halo_masses, halo_growth_rates, redshifts)
            t1 = _time.time()

            m_halo_final = halo_masses[:, -1]
            m_star_final = result["stars_mass"][:, -1]
            resolved = m_halo_final > 0
            n_formed = int(np.sum(resolved))
            median_star = np.median(m_star_final[resolved]) if n_formed else np.nan
            print(f"[gamma={gamma},omega={omega}] M_halo(z=0)={mass_bin:.2e}: "
                  f"median M_star={median_star:.3e} ({t1-t0:.1f}s)", flush=True)
            per_bin_median.append((mass_bin, median_star))
    finally:
        restore_reion_shape(old)
    return np.array(per_bin_median)


def main_():
    configs = [
        (15, 2, "default (gamma=15, omega=2)", "#1f77b4"),
        (15, 1, "gamma=15, omega=1 (softer mass cutoff)", "#9467bd"),
        (15, 4, "gamma=15, omega=4 (sharper mass cutoff)", "#8c564b"),
        (5, 2, "gamma=5, omega=2 (softer z transition)", "#ff7f0e"),
        (30, 2, "gamma=30, omega=2 (sharper z transition)", "#e377c2"),
    ]

    results = {}
    for gamma, omega, label, color in configs:
        results[label] = (run_ensemble(gamma, omega), color)

    mhalo_grid = np.logspace(7.8, 14.2, 200)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax, ax2 = axes

    for label, (median, color) in results.items():
        ax.plot(median[:, 0], median[:, 1], "o-", color=color, lw=1.8, ms=5, label=label)
        ax2.plot(median[:, 0], median[:, 1] / median[:, 0], "o-", color=color, lw=1.8, ms=5)

    ax.plot(mhalo_grid, moster2013_mstar(mhalo_grid), "--", color="0.3", lw=1.5, label="Moster+2013 (literature)")
    ax.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid), ":", color="0.3", lw=1.8, label="Behroozi+2013 (literature)")
    ax2.plot(mhalo_grid, moster2013_mstar(mhalo_grid) / mhalo_grid, "--", color="0.3", lw=1.5)
    ax2.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid) / mhalo_grid, ":", color="0.3", lw=1.8)

    for a in (ax, ax2):
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]")
        a.grid(alpha=0.3)
    ax.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$] (per-bin median)")
    ax2.set_ylabel(r"$M_\star/M_{\rm halo}$ at $z=0$ (per-bin median)")
    ax.legend(fontsize=7)

    fig.suptitle("Ashvini SHMR sensitivity to reionization-suppression shape (Okamoto+2008 gamma/omega)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("reionization_gamma_omega_sensitivity.png", dpi=150)
    print("Wrote reionization_gamma_omega_sensitivity.png")


if __name__ == "__main__":
    main_()
