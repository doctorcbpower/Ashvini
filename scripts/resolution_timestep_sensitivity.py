#!/usr/bin/env python3
"""
Checks whether the low-mass z=0 SHMR (M_halo ~ 1e8-1e12 Msun, where the
reionization on/off test showed the suppression/kink is dominated by
UVB physics -- see stellar_to_halo_mass_relation.py) is also sensitive to
the tree-growth timestep (dz) or mass resolution (M_res = m_res_fraction *
mass_bin), independent of reionization. Reionization stays at its
run_params.yaml default (on) throughout -- this isolates numerical
resolution/timestepping choices, not the reionization-shape question
already covered by reionization_gamma_omega_sensitivity.py.

Motivation for checking timestepping specifically: dz is fixed in
REDSHIFT, but dz maps to very different amounts of cosmic TIME dt at
different z (dt/dz shrinks quickly with increasing z), and Ashvini's own
physics (SF free-fall time, SN delay_time=0.015 Gyr, AGN feedback_delay_time)
operates on a cosmic-time clock, not a redshift clock -- so a dz coarse
enough to under-resolve those timescales could plausibly show up as a
mass-dependent artifact if different mass bins spend different amounts of
their assembly history at different z.
"""
import time as _time

import numpy as np
import matplotlib.pyplot as plt

from ashvini import main, pymctrees_adapter

from stellar_to_halo_mass_relation import moster2013_mstar, behroozi2013_mstar

PYMCTREES_CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 12, 9)  # focus on the low-mass/transition region
N_HALOS = 1000
Z0, Z_MAX = 0.0, 15.0
SEED = 42
BACKEND = "numba"


def run_ensemble(dz, m_res_fraction):
    per_bin_median = []
    for mass_bin in MASS_BINS:
        m_res = m_res_fraction * mass_bin
        t0 = _time.time()
        halo_masses, halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
            pymctrees_config_path=PYMCTREES_CONFIG, mass_bin=mass_bin,
            n_halos=N_HALOS, z0=Z0, z_max=Z_MAX, dz=dz,
            m_res=m_res, backend=BACKEND, seed=SEED,
        )
        result = main.run_forest(halo_masses, halo_growth_rates, redshifts)
        t1 = _time.time()

        m_halo_final = halo_masses[:, -1]
        m_star_final = result["stars_mass"][:, -1]
        resolved = m_halo_final > 0
        n_formed = int(np.sum(resolved))
        median_star = np.median(m_star_final[resolved]) if n_formed else np.nan
        n_steps = int(round((Z_MAX - Z0) / dz))
        print(f"[dz={dz},m_res_frac={m_res_fraction}] M_halo={mass_bin:.2e}: "
              f"{n_steps} steps, median M_star={median_star:.3e} ({t1-t0:.1f}s)", flush=True)
        per_bin_median.append((mass_bin, median_star))
    return np.array(per_bin_median)


def main_():
    configs = [
        (0.05, 1e-3, "default (dz=0.05, M_res=1e-3*M0)", "#1f77b4"),
        (0.01, 1e-3, "dz=0.01 (5x finer timestep)", "#d62728"),
        (0.005, 1e-3, "dz=0.005 (10x finer timestep)", "#9467bd"),
        (0.05, 1e-4, "M_res=1e-4*M0 (10x finer resolution)", "#2ca02c"),
    ]

    results = {}
    for dz, m_res_fraction, label, color in configs:
        results[label] = (run_ensemble(dz, m_res_fraction), color)

    mhalo_grid = np.logspace(7.8, 12.2, 100)

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

    fig.suptitle("Ashvini low-mass SHMR: sensitivity to timestep (dz) and mass resolution (M_res), reionization ON (default)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("resolution_timestep_sensitivity.png", dpi=150)
    print("Wrote resolution_timestep_sensitivity.png")


if __name__ == "__main__":
    main_()
