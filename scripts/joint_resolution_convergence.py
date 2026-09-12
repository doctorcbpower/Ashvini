#!/usr/bin/env python3
"""
Joint dz x M_res convergence test for the low-mass z=0 SHMR, following up
on resolution_timestep_sensitivity.py's finding that dz and M_res are each
independently unconverged at run_params.yaml's defaults (dz=0.05,
M_res=1e-3*M0), pulling the low-mass stellar mass in OPPOSITE directions
(finer dz -> lower M_star, by ~7x from dz=0.05->0.01; finer M_res -> higher
M_star, by ~7.6x from 1e-3->1e-4) -- so a joint grid is needed to find
where both are simultaneously fine enough, rather than trusting either
one-at-a-time test.

Reionization stays at its run_params.yaml default (on) throughout.
"""
import time as _time

import numpy as np
import matplotlib.pyplot as plt

from ashvini import main, pymctrees_adapter

from stellar_to_halo_mass_relation import moster2013_mstar, behroozi2013_mstar

PYMCTREES_CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 12, 9)
N_HALOS = 1000
Z0, Z_MAX = 0.0, 30.0  # raised from 15 -- see docs/RESOLUTION_CONVERGENCE.md: z_max=15 was pinning
                       # formation redshift at the finest M_res tested, confounding the convergence test
SEED = 42
BACKEND = "numba"

DZ_VALUES = [0.05, 0.01, 0.005]
MRES_FRAC_VALUES = [1e-3, 1e-4, 1e-5]


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
        print(f"[dz={dz},m_res_frac={m_res_fraction}] M_halo={mass_bin:.2e}: "
              f"median M_star={median_star:.3e} ({t1-t0:.1f}s)", flush=True)
        per_bin_median.append((mass_bin, median_star))
    return np.array(per_bin_median)


def main_():
    grid = {}
    for dz in DZ_VALUES:
        for m_res_frac in MRES_FRAC_VALUES:
            grid[(dz, m_res_frac)] = run_ensemble(dz, m_res_frac)

    # Convergence check: relative change in the lowest-mass bin's median
    # M_star as each axis is refined one step further, holding the other
    # axis fixed at its finest tested value -- if both are small, the
    # finest-tested corner is a safe converged answer; if not, the grid
    # needs to extend further in whichever axis is still moving.
    finest_dz, finest_mres = min(DZ_VALUES), min(MRES_FRAC_VALUES)

    def lowmass_star(dz, m_res_frac):
        arr = grid[(dz, m_res_frac)]
        return arr[0, 1]  # first (lowest) mass bin

    print("\n--- Convergence summary (lowest mass bin, M_halo={:.1e}) ---".format(MASS_BINS[0]))
    for dz in DZ_VALUES:
        print(f"dz={dz}, M_res_frac={finest_mres} (finest): M_star={lowmass_star(dz, finest_mres):.3e}")
    for m_res_frac in MRES_FRAC_VALUES:
        print(f"dz={finest_dz} (finest), M_res_frac={m_res_frac}: M_star={lowmass_star(finest_dz, m_res_frac):.3e}")

    mhalo_grid = np.logspace(7.8, 12.2, 100)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax, ax2 = axes
    cmap_dz = {0.05: "o-", 0.01: "s--", 0.005: "^:"}
    colors_mres = {1e-3: "#1f77b4", 1e-4: "#2ca02c", 1e-5: "#d62728"}

    for (dz, m_res_frac), median in grid.items():
        style = cmap_dz[dz]
        color = colors_mres[m_res_frac]
        label = f"dz={dz}, M_res=1e{int(np.log10(m_res_frac))}*M0"
        ax.plot(median[:, 0], median[:, 1], style, color=color, lw=1.5, ms=5, label=label)
        ax2.plot(median[:, 0], median[:, 1] / median[:, 0], style, color=color, lw=1.5, ms=5)

    ax.plot(mhalo_grid, moster2013_mstar(mhalo_grid), "--", color="0.3", lw=1.5, label="Moster+2013")
    ax.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid), ":", color="0.3", lw=1.8, label="Behroozi+2013")
    ax2.plot(mhalo_grid, moster2013_mstar(mhalo_grid) / mhalo_grid, "--", color="0.3", lw=1.5)
    ax2.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid) / mhalo_grid, ":", color="0.3", lw=1.8)

    for a in (ax, ax2):
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]")
        a.grid(alpha=0.3)
    ax.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$] (per-bin median)")
    ax2.set_ylabel(r"$M_\star/M_{\rm halo}$ at $z=0$ (per-bin median)")
    ax.legend(fontsize=6, ncol=2)

    fig.suptitle("Ashvini low-mass SHMR: joint dz x M_res convergence grid, reionization ON (default)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig("joint_resolution_convergence.png", dpi=150)
    print("Wrote joint_resolution_convergence.png")


if __name__ == "__main__":
    main_()
