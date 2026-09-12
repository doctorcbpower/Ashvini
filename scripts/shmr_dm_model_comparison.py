#!/usr/bin/env python3
"""
Stellar-to-halo mass relation (SHMR) at z=0 for CDM, WDM, and FDM, using
Ashvini's baryonic model on foraois-generated PCH08 trees.

SIDM is NOT included: foraois's own collapse/sidm.py raises
NotImplementedError unconditionally (no collapse barrier has been derived
or even scoped for SIDM -- see foraois's ROADMAP.md and sidm.py's module
docstring) -- there is no placeholder to fall back on the way FDM has one,
so attempting it would not produce a real result.

Resolution/timestep settings follow docs/RESOLUTION_CONVERGENCE.md's
findings: dz has no independent effect once M_res is adequately resolved,
so dz=0.01 is used (cheaper than 0.005, empirically equivalent once M_res
is fine); m_res_fraction=1e-5 is the finest value verified NOT to hit the
z_max ceiling artifact at z_max=30 (still not proven fully converged --
see that doc's follow-up item 1 -- so treat absolute low-mass M_star
values here as the best currently-available estimate, not a final
answer). z_max=30 (not the notebooks' usual 15) specifically to avoid the
formation-redshift pinning artifact documented there.
"""
import time as _time

import numpy as np
import matplotlib.pyplot as plt

from ashvini import main, pymctrees_adapter

from stellar_to_halo_mass_relation import moster2013_mstar, behroozi2013_mstar

FORAOIS_CONFIG_DIR = "/Users/00075868/MyCodes/foraois/config"
DM_MODELS = {
    "CDM": f"{FORAOIS_CONFIG_DIR}/planck2018_camb.yml",
    "WDM (3 keV thermal relic)": f"{FORAOIS_CONFIG_DIR}/planck2018_wdm.yml",
    "FDM (m=1e-22 eV)": f"{FORAOIS_CONFIG_DIR}/planck2018_fdm.yml",
}
COLORS = {
    "CDM": "#1f77b4",
    "WDM (3 keV thermal relic)": "#ff7f0e",
    "FDM (m=1e-22 eV)": "#2ca02c",
}

MASS_BINS = np.logspace(8, 14, 13)
N_HALOS = 1000
Z0, Z_MAX, DZ = 0.0, 30.0, 0.01
M_RES_FRACTION = 1e-5
SEED = 42
BACKEND = "numba"


def run_ensemble(config_path):
    per_bin_median = []
    halo_pts, star_pts = [], []
    for mass_bin in MASS_BINS:
        m_res = M_RES_FRACTION * mass_bin
        t0 = _time.time()
        try:
            halo_masses, halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
                pymctrees_config_path=config_path, mass_bin=mass_bin,
                n_halos=N_HALOS, z0=Z0, z_max=Z_MAX, dz=DZ,
                m_res=m_res, backend=BACKEND, seed=SEED,
            )
        except Exception as exc:
            print(f"[{config_path}] M_halo={mass_bin:.2e}: FAILED ({exc!r})", flush=True)
            per_bin_median.append((mass_bin, np.nan))
            continue
        result = main.run_forest(halo_masses, halo_growth_rates, redshifts)
        t1 = _time.time()

        m_halo_final = halo_masses[:, -1]
        m_star_final = result["stars_mass"][:, -1]
        resolved = m_halo_final > 0
        n_formed = int(np.sum(resolved))
        median_star = np.median(m_star_final[resolved]) if n_formed else np.nan
        print(f"[{config_path}] M_halo={mass_bin:.2e}: {n_formed}/{N_HALOS} resolved, "
              f"median M_star={median_star:.3e} ({t1-t0:.1f}s)", flush=True)
        per_bin_median.append((mass_bin, median_star))
        halo_pts.append(m_halo_final[resolved])
        star_pts.append(m_star_final[resolved])

    halo_pts = np.concatenate(halo_pts) if halo_pts else np.array([])
    star_pts = np.concatenate(star_pts) if star_pts else np.array([])
    return np.array(per_bin_median), halo_pts, star_pts


def main_():
    import warnings
    warnings.filterwarnings("ignore", message="build_forest_live: m_res/mass_bin")

    results = {}
    for label, config_path in DM_MODELS.items():
        print(f"\n=== {label} ===")
        median, halo_pts, star_pts = run_ensemble(config_path)
        results[label] = (median, halo_pts, star_pts)

    mhalo_grid = np.logspace(7.8, 14.2, 200)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax, ax2 = axes

    for label, (median, halo_pts, star_pts) in results.items():
        color = COLORS[label]
        ax.scatter(halo_pts, star_pts, s=6, alpha=0.08, color=color)
        ax.plot(median[:, 0], median[:, 1], "o-", color=color, lw=2.6, ms=6, label=label)
        ax2.plot(median[:, 0], median[:, 1] / median[:, 0], "o-", color=color, lw=2.6, ms=6)

    ax.plot(mhalo_grid, moster2013_mstar(mhalo_grid), "--", color="0.3", lw=2.2, label="Moster+2013 (literature)")
    ax.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid), ":", color="0.3", lw=2.4, label="Behroozi+2013 (literature)")
    ax2.plot(mhalo_grid, moster2013_mstar(mhalo_grid) / mhalo_grid, "--", color="0.3", lw=2.2)
    ax2.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid) / mhalo_grid, ":", color="0.3", lw=2.4)

    for a in (ax, ax2):
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]", fontsize=15)
        a.tick_params(labelsize=12)
        a.grid(alpha=0.3)
    ax.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$]", fontsize=15)
    ax2.set_ylabel(r"$M_\star/M_{\rm halo}$ at $z=0$ (per-bin median)", fontsize=15)
    ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig("shmr_dm_model_comparison.png", dpi=150)
    print("Wrote shmr_dm_model_comparison.png")


if __name__ == "__main__":
    main_()
