#!/usr/bin/env python3
"""
Pulls each halo's formation redshift (the first/highest-z step where its
tracked mass exceeds 0 -- i.e. the first time PCH08/foraois resolves a
progenitor above M_res) out of the same dz x M_res grid used in
joint_resolution_convergence.py, to test a specific hypothesis about why
M_star kept growing by ~7-8x per decade of M_res refinement with no sign
of turning over (see docs/RESOLUTION_CONVERGENCE.md):

Does refining M_res mostly reveal MORE, comparably-early small-scale
structure (in which case formation redshift should stay roughly flat as
M_res is refined, since you're resolving more contemporaneous structure at
similar epochs) -- or does it mostly push each halo's own EFFECTIVE
FORMATION TIME earlier (in which case median formation redshift should
keep climbing with no sign of leveling off, mirroring M_star's own
behaviour, since a halo's "pre-formation prefix" before it first exceeds
M_res is zeroed out entirely by build_forest_for_bin/_mask_unresolved_prefix
rather than credited as smooth accretion)?

halo_masses has shape (n_halos, n_steps+1), index 0 = earliest (highest z)
per pymctrees_adapter.py's own chronological convention. A halo's
formation index is the first column where its mass is > 0; formation
redshift is redshifts[] at that index. A halo that never resolves (mass
0 at every step, i.e. mass_bin's own M_res is too coarse for it to form
at all within z_max) is excluded, same resolved-mask convention used
throughout the SHMR scripts.
"""
import warnings

import numpy as np
import matplotlib.pyplot as plt

from ashvini import pymctrees_adapter

# This script deliberately spans the already-documented unconverged
# m_res_fraction=1e-3 default (see docs/RESOLUTION_CONVERGENCE.md) as one
# point on the grid, so build_forest_live's convergence warning would
# otherwise fire on every one of those calls -- suppress it here, not by
# weakening the warning itself.
warnings.filterwarnings("ignore", message="build_forest_live: m_res/mass_bin")

PYMCTREES_CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 12, 9)
N_HALOS = 1000
Z0, Z_MAX = 0.0, 30.0  # raised from 15 -- the earlier z_max=15 run pinned formation redshift
                       # at the ceiling for M_res_frac<=1e-4, confounding the convergence test
SEED = 42
BACKEND = "numba"

DZ_VALUES = [0.05, 0.01, 0.005]
MRES_FRAC_VALUES = [1e-3, 1e-4, 1e-5]


def formation_redshifts(dz, m_res_fraction, mass_bin):
    m_res = m_res_fraction * mass_bin
    halo_masses, _halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
        pymctrees_config_path=PYMCTREES_CONFIG, mass_bin=mass_bin,
        n_halos=N_HALOS, z0=Z0, z_max=Z_MAX, dz=dz,
        m_res=m_res, backend=BACKEND, seed=SEED,
    )
    n_halos, n_steps = halo_masses.shape
    formed_mask = halo_masses > 0
    ever_formed = formed_mask.any(axis=1)
    first_idx = np.argmax(formed_mask, axis=1)  # first True per row; garbage where never formed

    z_form = redshifts[first_idx[ever_formed]]
    return z_form, int(np.sum(ever_formed))


def main_():
    print("--- Median formation redshift vs. M_res, at each dz (lowest 3 mass bins) ---")
    low_mass_bins = MASS_BINS[:3]

    fig, axes = plt.subplots(1, len(low_mass_bins), figsize=(15, 4.5), sharey=True)

    results = {}  # (dz, m_res_frac, mass_bin) -> z_form array
    for mass_bin in low_mass_bins:
        for dz in DZ_VALUES:
            for m_res_frac in MRES_FRAC_VALUES:
                z_form, n_formed = formation_redshifts(dz, m_res_frac, mass_bin)
                results[(dz, m_res_frac, mass_bin)] = z_form
                print(f"M_halo={mass_bin:.2e}, dz={dz}, M_res_frac={m_res_frac}: "
                      f"{n_formed}/{N_HALOS} formed, median z_form={np.median(z_form):.2f}, "
                      f"90th pctile z_form={np.percentile(z_form, 90):.2f}, max z_form={z_form.max():.2f}")

    colors_mres = {1e-3: "#1f77b4", 1e-4: "#2ca02c", 1e-5: "#d62728"}
    styles_dz = {0.05: "o-", 0.01: "s--", 0.005: "^:"}

    for ax, mass_bin in zip(axes, low_mass_bins):
        for dz in DZ_VALUES:
            medians = [np.median(results[(dz, m_res_frac, mass_bin)]) for m_res_frac in MRES_FRAC_VALUES]
            ax.plot(MRES_FRAC_VALUES, medians, styles_dz[dz], color="0.2", alpha=0.4 + 0.3 * DZ_VALUES.index(dz),
                    label=f"dz={dz}")
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel(r"$M_{\rm res}/M_0$")
        ax.set_title(f"$M_2={mass_bin:.1e}\\,M_\\odot$")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("median formation redshift")
    axes[0].legend(fontsize=8)

    fig.suptitle("Median halo formation redshift vs. mass resolution (does it plateau or keep climbing?)")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig("formation_redshift_vs_resolution.png", dpi=150)
    print("Wrote formation_redshift_vs_resolution.png")


if __name__ == "__main__":
    main_()
