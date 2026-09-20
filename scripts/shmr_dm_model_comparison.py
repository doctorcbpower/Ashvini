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
Z0, Z_MAX, DZ = 0.0, 30.0, 0.005
M_RES_FRACTION = 1e-4
TUNED_E_FF = 0.05
TUNED_M_HOT, TUNED_HOT_FLOOR, HOT_PHI = 2e12, 0.1, 4.0
ALGORITHM = "zh"  # Zhang-Hui trees; see docs/RESOLUTION_CONVERGENCE.md (correction) and scripts/zh_vs_pch08_mres_scan.py
SEED = 42
BACKEND = "numba"


SUMMARY = []


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
                m_res=m_res, backend=BACKEND, seed=SEED, algorithm=ALGORITHM,
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
        frac_star = float(np.mean(m_star_final[resolved] > 0)) if n_formed else np.nan
        SUMMARY.append({"M_halo": float(mass_bin), "n_resolved": n_formed, "n": N_HALOS,
                        "median": float(median_star), "p16": float(np.percentile(m_star_final[resolved], 16)) if n_formed else np.nan,
                        "p84": float(np.percentile(m_star_final[resolved], 84)) if n_formed else np.nan,
                        "frac_with_stars": frac_star, "config": config_path})
        print(f"[{config_path}] M_halo={mass_bin:.2e}: {n_formed}/{N_HALOS} resolved, "
              f"median M_star={median_star:.3e} ({t1-t0:.1f}s)", flush=True)
        per_bin_median.append((mass_bin, median_star))
        halo_pts.append(m_halo_final[resolved])
        star_pts.append(m_star_final[resolved])

    halo_pts = np.concatenate(halo_pts) if halo_pts else np.array([])
    star_pts = np.concatenate(star_pts) if star_pts else np.array([])
    return np.array(per_bin_median), halo_pts, star_pts


def apply_variant(variant):
    """'standard' leaves Ashvini's defaults. 'tuned' is ILLUSTRATIVE ONLY. It (1) swaps the UV-suppression term for
    the mass-only Okamoto+2008 form (as the MVM uses; docs/UV_SUPPRESSION_AUDIT.md), (2) raises the star-formation
    efficiency e_ff from 0.015 to 0.05, and (3) multiplies the inflow by a cold-to-hot suppression with a residual
    floor, zeta_ch = f + (1 - f)(1 - s(M/M_hot, phi)), M_hot = 2e12 Msun, f = 0.1, phi = 4. All three were chosen by
    eye against Behroozi+2013 (scripts/output/eff_epsp_scan.json, hot_mode_floor_scan.json). epsilon_p, the SN
    loading and the AGN wind are unchanged."""
    if variant == "tuned":
        from ashvini import reionization, main as ash_main
        from ashvini.paper_reservoir import cold_hot_mode_suppression, uv_suppression_mass_only

        def uv_with_hot_mode(z, m, mdot):
            m = np.asarray(m, dtype=float)
            zeta = TUNED_HOT_FLOOR + (1.0 - TUNED_HOT_FLOOR) * cold_hot_mode_suppression(m, M_hot=TUNED_M_HOT, phi=HOT_PHI)
            return uv_suppression_mass_only(np.asarray(z, dtype=float), m) * zeta

        reionization.uv_suppression = uv_with_hot_mode
        ash_main.e_ff = TUNED_E_FF


def main_():
    import argparse
    import warnings
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["standard", "tuned"], default="standard")
    variant = ap.parse_args().variant
    apply_variant(variant)
    suffix = "" if variant == "standard" else "_tuned"
    warnings.filterwarnings("ignore", message="build_forest_live: m_res/mass_bin")

    results = {}
    for label, config_path in DM_MODELS.items():
        print(f"\n=== {label} ===")
        median, halo_pts, star_pts = run_ensemble(config_path)
        results[label] = (median, halo_pts, star_pts)

    import json
    with open(f"scripts/output/shmr_dm_model_comparison_zh{suffix}.json", "w") as f:
        json.dump(SUMMARY, f, indent=1)
    np.savez(f"scripts/output/shmr_dm_model_comparison_zh{suffix}.npz",
             **{f"{k}|{name}": v for k, arrs in results.items() for name, v in zip(("median", "halo", "star"), arrs)})
    make_figure(results, variant, suffix)


def make_figure(results, variant, suffix):
    """Full-width (7.1 in) figure in the SciencePlots style shared with the paper's other figures."""
    import json
    import sys
    sys.path.insert(0, "scripts/paper_figures")
    import paper_style as ps
    ps.apply()
    colors = {"CDM": ps.BLUE, "WDM (3 keV thermal relic)": ps.RED, "FDM (m=1e-22 eV)": ps.GREEN}
    names = {"CDM": "CDM", "WDM (3 keV thermal relic)": "WDM, 3 keV", "FDM (m=1e-22 eV)": r"FDM, $10^{-22}$ eV"}
    mhalo_grid = np.logspace(7.8, 14.2, 200)

    if variant == "tuned":  # every halo forms stars here, so a stars-fraction panel would carry no information
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(ps.FULL, 2.9))
        ax3 = None
    else:
        fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(ps.FULL, 2.9), gridspec_kw={"width_ratios": [1, 1, 0.7]})

    for label, (median, halo_pts, star_pts) in results.items():
        keep = star_pts > 0
        ax.scatter(halo_pts[keep], star_pts[keep], s=1.5, alpha=0.08, color=colors[label], lw=0, rasterized=True)
        m = np.where(median[:, 1] > 0, median[:, 1], np.nan)  # a zero median is not drawn
        ax.plot(median[:, 0], m, "o-", color=colors[label], lw=1.5, ms=3, label=names[label])
        ax2.plot(median[:, 0], m / median[:, 0], "o-", color=colors[label], lw=1.5, ms=3)
        if ax3 is not None:
            frac = [x["frac_with_stars"] for x in SUMMARY if x["config"] == DM_MODELS[label]]
            ax3.plot(median[:, 0], frac, "o-", color=colors[label], lw=1.5, ms=3)

    for f, ls, lab in ((moster2013_mstar, "--", "Moster+2013"), (behroozi2013_mstar, ":", "Behroozi+2013")):
        ax.plot(mhalo_grid, f(mhalo_grid), ls, color="0.3", lw=1.3, label=lab)
        ax2.plot(mhalo_grid, f(mhalo_grid) / mhalo_grid, ls, color="0.3", lw=1.3)

    if variant == "tuned":
        try:
            std = [x for x in json.load(open("scripts/output/shmr_dm_model_comparison_zh.json")) if x["config"] == DM_MODELS["CDM"]]
            sm = np.array([[x["M_halo"], x["median"]] for x in std])
            sm[sm[:, 1] <= 0, 1] = np.nan
            ax.plot(sm[:, 0], sm[:, 1], "-", color=ps.BLUE, lw=0.8, alpha=0.5, label="CDM, untuned")
            ax2.plot(sm[:, 0], sm[:, 1] / sm[:, 0], "-", color=ps.BLUE, lw=0.8, alpha=0.5)
        except FileNotFoundError:
            pass
        ax.text(0.97, 0.03, "illustrative: UV term, $\\epsilon_{\\rm ff}$ and hot-mode floor tuned", transform=ax.transAxes,
                ha="right", fontsize=5.5, color="0.25")

    for a in [x for x in (ax, ax2, ax3) if x is not None]:
        a.set_xscale("log")
        a.set_xlabel(r"$M_{\rm halo}(z=0)\ [{\rm M}_\odot]$")
    for a in (ax, ax2):
        a.set_yscale("log")
    ax.set_ylim(1e-1, 3e12)
    ax2.set_ylim(1e-6, 1e-1)
    ax.set_ylabel(r"$M_\star(z=0)\ [{\rm M}_\odot]$")
    ax2.set_ylabel(r"$M_\star/M_{\rm halo}$")
    if ax3 is not None:
        ax3.set_ylabel(r"fraction of halos with $M_\star>0$")
        ax3.set_ylim(-0.02, 1.02)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    ps.save(fig, f"shmr_dm_model_comparison{suffix}")
    print(f"Wrote shmr_dm_model_comparison{suffix}.pdf/.png")


if __name__ == "__main__":
    main_()
