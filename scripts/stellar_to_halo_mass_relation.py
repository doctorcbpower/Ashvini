#!/usr/bin/env python3
"""
Stellar-to-halo mass relation (SHMR) at z=0, using Ashvini's baryonic
model run on pymctrees-generated (foraois) PCH08 merger trees.

Each mass_bin below is the halo mass M_halo(z=0) an entire ensemble of
Monte Carlo trees is grown to exactly (this is how build_forest_live's
mass_bin parameter works -- see pymctrees_adapter.py's docstring): the
scatter in the resulting M_star at z=0 across n_halos realizations of the
same mass_bin comes entirely from each halo's individual Monte Carlo
accretion/merger history, not from any scatter in M_halo itself. This
matches the same per-bin-ensemble methodology already used for Figure 8
(stellar fraction vs halo mass) in notebooks/menon_power_2024_all_figures.ipynb,
just evaluated at z=0 with M_star directly rather than M_star/(M_star+M_gas)
at z=5/7/10.

All Ashvini physics parameters (star formation efficiency, delay time,
reionization, dust, black holes, ...) are left at run_params.yaml's
defaults except where explicitly toggled (see run_bh_on/run_bh_off below);
only the dark-matter/cosmology config (pymctrees_config_path) and
tree-generation settings below are otherwise explicit.

Moster et al. (2013) and Behroozi, Wechsler & Conroy (2013) z=0
double-power-law SHMR fits are overlaid for sanity-check comparison only --
Ashvini/foraois are not tuned/calibrated to reproduce either, so agreement
or disagreement is informative, not a pass/fail test.

BH growth off = all three seeding channels (pop3, direct_collapse,
halo_mass_threshold) disabled, so no halo ever seeds a black hole -- same
approach as notebooks/pymctrees_ashvini_demo.ipynb's own bh_on/bh_off
comparison cells, not a separate/new switch.
"""
import time as _time

import numpy as np
import matplotlib.pyplot as plt

from ashvini import main, pymctrees_adapter
from ashvini.run_params import PARAMS

PYMCTREES_CONFIG = "/Users/00075868/MyCodes/foraois/config/planck2018_camb.yml"
MASS_BINS = np.logspace(8, 14, 13)  # Msun, at z0=0
N_HALOS = 1000
Z0, Z_MAX, DZ = 0.0, 15.0, 0.05
M_RES = None  # defaults to 1e-3 * mass_bin
SEED = 42
BACKEND = "numba"


def moster2013_mstar(mhalo):
    """Moster, Naab & White (2013) double power-law SHMR fit, z=0
    coefficients (their Table 1) -- literature reference curve only, not
    a target Ashvini is calibrated against."""
    M1 = 10 ** 11.590
    N = 0.0351
    beta = 1.376
    gamma = 0.608
    x = mhalo / M1
    return mhalo * (2 * N) / (x ** (-beta) + x ** gamma)


def behroozi2013_mstar(mhalo):
    """Behroozi, Wechsler & Conroy (2013) SHMR fit, z=0 best-fit
    parameters (their Table 2) -- literature reference curve only, not a
    target Ashvini is calibrated against. Eq. 3-4 of that paper:
    log10(Ms(Mh)) = log10(eps*M1) + f(log10(Mh/M1)) - f(0), with
    f(x) = -log10(10^(alpha x)+1) + delta*(log10(1+exp(x)))^gamma / (1+exp(10^-x))."""
    log10_M1 = 11.514
    log10_eps = -1.777
    alpha = -1.412
    delta = 3.508
    gamma = 0.316

    def f(x):
        return -np.log10(10 ** (alpha * x) + 1.0) + delta * (np.log10(1.0 + np.exp(x))) ** gamma / (1.0 + np.exp(10 ** (-x)))

    log10_mhalo = np.log10(mhalo)
    log10_mstar = (log10_eps + log10_M1) + f(log10_mhalo - log10_M1) - f(0.0)
    return 10 ** log10_mstar


def run_ensemble(disable_bh_seeding=False, disable_uvb=False):
    orig_enabled = {
        "pop3": PARAMS.bh.seeding.pop3.enabled,
        "direct_collapse": PARAMS.bh.seeding.direct_collapse.enabled,
        "halo_mass_threshold": PARAMS.bh.seeding.halo_mass_threshold.enabled,
    }
    if disable_bh_seeding:
        PARAMS.bh.seeding.pop3.enabled = False
        PARAMS.bh.seeding.direct_collapse.enabled = False
        PARAMS.bh.seeding.halo_mass_threshold.enabled = False

    # main.UV_background is a module-level snapshot of PARAMS.reion.UVB_enabled
    # taken at import time (same pattern as main.sn_type -- see
    # notebooks/pymctrees_ashvini_demo.ipynb's own sn_type toggling), so
    # PARAMS.reion.UVB_enabled itself must be patched via main.UV_background,
    # not by assigning to PARAMS.
    orig_uvb = main.UV_background
    if disable_uvb:
        main.UV_background = False

    tag = f"bh_seeding={'off' if disable_bh_seeding else 'on'},uvb={'off' if disable_uvb else 'on'}"
    halo_mass_pts, stellar_mass_pts, per_bin_median = [], [], []
    try:
        for mass_bin in MASS_BINS:
            t0 = _time.time()
            halo_masses, halo_growth_rates, redshifts, _merger_mass = pymctrees_adapter.build_forest_live(
                pymctrees_config_path=PYMCTREES_CONFIG, mass_bin=mass_bin,
                n_halos=N_HALOS, z0=Z0, z_max=Z_MAX, dz=DZ,
                m_res=M_RES, backend=BACKEND, seed=SEED,
            )
            t1 = _time.time()
            result = main.run_forest(halo_masses, halo_growth_rates, redshifts)
            t2 = _time.time()

            m_halo_final = halo_masses[:, -1]
            m_star_final = result["stars_mass"][:, -1]
            n_seeded = int(np.sum(result["bh_mass"][:, -1] > 0))

            resolved = m_halo_final > 0
            n_formed = int(np.sum(resolved))
            print(f"[{tag}] "
                  f"M_halo(z=0)={mass_bin:.2e} Msun: {n_formed}/{N_HALOS} resolved, {n_seeded} BH-seeded, "
                  f"median M_star={np.median(m_star_final[resolved]) if n_formed else float('nan'):.3e} Msun "
                  f"(tree {t1-t0:.1f}s, ashvini {t2-t1:.1f}s)", flush=True)

            halo_mass_pts.append(m_halo_final[resolved])
            stellar_mass_pts.append(m_star_final[resolved])
            per_bin_median.append((mass_bin, np.median(m_star_final[resolved]) if n_formed else np.nan))
    finally:
        PARAMS.bh.seeding.pop3.enabled = orig_enabled["pop3"]
        PARAMS.bh.seeding.direct_collapse.enabled = orig_enabled["direct_collapse"]
        PARAMS.bh.seeding.halo_mass_threshold.enabled = orig_enabled["halo_mass_threshold"]
        main.UV_background = orig_uvb

    return (
        np.concatenate(halo_mass_pts),
        np.concatenate(stellar_mass_pts),
        np.array(per_bin_median),
    )


def main_():
    halo_pts_on, star_pts_on, median_on = run_ensemble(disable_bh_seeding=False, disable_uvb=False)
    _, _, median_bh_off = run_ensemble(disable_bh_seeding=True, disable_uvb=False)
    _, _, median_uvb_off = run_ensemble(disable_bh_seeding=False, disable_uvb=True)

    mhalo_grid = np.logspace(7.8, 14.2, 200)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.scatter(halo_pts_on, star_pts_on, s=8, alpha=0.15, color="#1f77b4")
    ax.plot(median_on[:, 0], median_on[:, 1], "o-", color="#1f77b4", lw=1.8, ms=6, label="Ashvini, default (BH+UVB on)")
    ax.plot(median_bh_off[:, 0], median_bh_off[:, 1], "o-", color="#d62728", lw=1.8, ms=6, label="Ashvini, BH growth OFF")
    ax.plot(median_uvb_off[:, 0], median_uvb_off[:, 1], "o-", color="#2ca02c", lw=1.8, ms=6, label="Ashvini, UVB/reionization OFF")
    ax.plot(mhalo_grid, moster2013_mstar(mhalo_grid), "--", color="0.3", lw=1.5, label="Moster+2013 (z=0, literature)")
    ax.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid), ":", color="0.3", lw=1.8, label="Behroozi+2013 (z=0, literature)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]")
    ax.set_ylabel(r"$M_\star(z=0)$ [M$_\odot$]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax2 = axes[1]
    ax2.plot(median_on[:, 0], median_on[:, 1] / median_on[:, 0], "o-", color="#1f77b4", lw=1.8, ms=6)
    ax2.plot(median_bh_off[:, 0], median_bh_off[:, 1] / median_bh_off[:, 0], "o-", color="#d62728", lw=1.8, ms=6)
    ax2.plot(median_uvb_off[:, 0], median_uvb_off[:, 1] / median_uvb_off[:, 0], "o-", color="#2ca02c", lw=1.8, ms=6)
    ax2.plot(mhalo_grid, moster2013_mstar(mhalo_grid) / mhalo_grid, "--", color="0.3", lw=1.5)
    ax2.plot(mhalo_grid, behroozi2013_mstar(mhalo_grid) / mhalo_grid, ":", color="0.3", lw=1.8)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$M_{\rm halo}(z=0)$ [M$_\odot$]")
    ax2.set_ylabel(r"$M_\star/M_{\rm halo}$ at $z=0$ (per-bin median)")
    ax2.grid(alpha=0.3)

    fig.suptitle(f"Ashvini SHMR, z=0 (n_halos={N_HALOS}/bin, CDM): BH growth and reionization on/off, vs. literature")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("stellar_to_halo_mass_relation.png", dpi=150)
    print("Wrote stellar_to_halo_mass_relation.png")


if __name__ == "__main__":
    main_()
