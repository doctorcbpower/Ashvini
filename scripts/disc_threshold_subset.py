"""
Does giving the gas a disc, with a star formation threshold on its surface density, change the predictions?

For a small subset of haloes (the first three of each mass bin of the test trees) run the event integrator and the grid
integrator with and without the disc option, over a range of its parameters, and compare with a cooling threshold
(the halo must reach a virial temperature before gas forms stars) on the same haloes. Output: for each variant the
median over the haloes of (stellar mass, gas metallicity, redshift of the first star) relative to the plain model, and
the lowest and highest of the haloes. Writes scripts/output/disc_threshold_subset.csv.

The growth-rate column of the test trees does not match their mass history (they were thinned in redshift), so the
supply is taken from the mass history for both integrators.

Usage: python scripts/disc_threshold_subset.py
"""
import csv
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
warnings.filterwarnings("ignore")

from ashvini import utils  # noqa: E402
from ashvini.components import disc  # noqa: E402
from ashvini.config import DEFAULT, override  # noqa: E402
from ashvini.main import run_events, run_forest  # noqa: E402

TREES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures", "merger_trees_fixture.h5")
BINS = (1e8, 1e10, 1e11)
N_HALOES = 3
EVENT_BINS = (1e8, 1e10)                 # the event integrator is slow for the largest haloes (every supernova is followed)
BASE = {"reion.UVB_enabled": False, "bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
        "bh.seeding.halo_mass_threshold.enabled": False,
        "sn.m_sn_max": 100.0, "sn.delay_distribution": "stellar_lifetimes", "sn.wind_model": "momentum", "sn.f_p": 0.7646}
VARIANTS = {
    "plain": {},
    "disc, Sigma_crit = 1": {"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 1.0},
    "disc, Sigma_crit = 5 (default)": {"disc.enabled": True},
    "disc, Sigma_crit = 50": {"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 50.0},
    "disc, Sigma_crit = 500": {"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 500.0},
    "disc, spin 0.02": {"disc.enabled": True, "disc.spin_parameter": 0.02},
    "disc, spin 0.07": {"disc.enabled": True, "disc.spin_parameter": 0.07},
    "disc, j_d/m_d = 0.6": {"disc.enabled": True, "disc.angular_momentum_ratio": 0.6},
    "cooling threshold 3000 K": {"events.cooling_T_vir_K": 3.0e3},
    "cooling threshold 10000 K": {"events.cooling_T_vir_K": 1.0e4},
    "disc 50 + cooling 3000 K": {"disc.enabled": True, "disc.critical_surface_density_Msun_pc2": 50.0, "events.cooling_T_vir_K": 3.0e3},
}


def measures(r, z):
    s = r["stars_mass"][0]
    g = r["gas_mass"][0, -1]
    first = int(np.argmax(s > 0)) if np.any(s > 0) else None
    return (s[-1], r["gas_metals"][0, -1] / g if g > 0 else np.nan, z[first] if first is not None else np.nan)


def main():
    rows = []
    for mb in BINS:
        hm, _, z = utils.read_trees(TREES, mb)
        for stepper in ("events", "grid"):
            if stepper == "events" and mb not in EVENT_BINS:
                continue
            results = {v: [] for v in VARIANTS}
            for i in range(N_HALOES):
                M = hm[i:i + 1]
                Mdot = np.maximum(np.gradient(M, utils.time_at_z(z), axis=1), 0.0)
                for name, change in VARIANTS.items():
                    if stepper == "grid" and "cooling" in name:
                        continue                         # the cooling threshold is an option of the event integrator
                    cfg = override(DEFAULT, **{**BASE, **change})
                    r = (run_events(M, Mdot, z, config=cfg) if stepper == "events" else run_forest(M, Mdot, z, config=cfg))
                    results[name].append(measures(r, z))
                print(f"  done {stepper} M_h={mb:.0e} halo {i}", flush=True)
            plain = np.array(results["plain"])
            for name, vals in results.items():
                if not vals:
                    continue
                v = np.array(vals)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio_s = np.where(plain[:, 0] > 0, v[:, 0] / plain[:, 0], np.nan)
                    ratio_z = np.where(plain[:, 1] > 0, v[:, 1] / plain[:, 1], np.nan)
                dz_first = v[:, 2] - plain[:, 2]
                rows.append({"halo_mass": mb, "integrator": stepper, "variant": name,
                             "Mstar_ratio_median": np.nanmedian(ratio_s), "Mstar_ratio_min": np.nanmin(ratio_s), "Mstar_ratio_max": np.nanmax(ratio_s),
                             "Zgas_ratio_median": np.nanmedian(ratio_z), "Zgas_ratio_min": np.nanmin(ratio_z), "Zgas_ratio_max": np.nanmax(ratio_z),
                             "first_star_dz_median": np.nanmedian(dz_first), "n_halos_without_stars": int(np.sum(~(v[:, 0] > 0)))})
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "output"), exist_ok=True)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "disc_threshold_subset.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"{'M_h':>6} {'integr.':>7} {'variant':34s} {'M* ratio (min-max)':>24s} {'Zgas ratio (min-max)':>24s} {'dz first star':>13s}")
    for r in rows:
        print(f"{r['halo_mass']:6.0e} {r['integrator']:>7} {r['variant']:34s} "
              f"{r['Mstar_ratio_median']:7.2f} ({r['Mstar_ratio_min']:5.2f}-{r['Mstar_ratio_max']:5.2f})  "
              f"{r['Zgas_ratio_median']:7.2f} ({r['Zgas_ratio_min']:5.2f}-{r['Zgas_ratio_max']:5.2f})  {r['first_star_dz_median']:+8.2f}")


if __name__ == "__main__":
    main()
