"""
Does a star formation law set by a disc differ from the current one, which takes a fraction of the gas per halo free-fall
time? And can the current law be made consistent with the disc law by changing its efficiency?

Compared on the same haloes (the first three of each mass bin of the test trees), with the grid integrator:
  current            SFR = eps * M_gas / (0.141 / H(z)), eps = sf.efficiency (0.015)
  Schmidt-Kennicutt  disc.star_formation_law = "schmidt_kennicutt": Sigma_SFR = A Sigma_gas^N over the exponential disc
Then, for the Schmidt-Kennicutt law, (1) the efficiency per halo free-fall time that it amounts to, and how much that
changes with time and mass, and (2) the constant efficiency that gives the same final stellar mass, and what else that
gets right.

A law of SFR = eps_d M_gas / t_disc (disc dynamical time) is not compared: for a disc whose size is a fixed fraction of the
virial radius it is exactly the current law with eps = eps_d * 0.141 * 10 sqrt(2) / (j_d/m_d * spin) (0.855 for eps_d = 0.015).
See components/disc.py and tests/test_disc.py.

Usage: python scripts/disc_star_formation_laws.py
"""
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
warnings.filterwarnings("ignore")

from ashvini import utils  # noqa: E402
from ashvini.config import DEFAULT, override  # noqa: E402
from ashvini.star_formation import time_freefall  # noqa: E402
from ashvini.integrators.grid import GridIntegrator  # noqa: E402

TREES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures", "merger_trees_fixture.h5")
BASE = {"reion.UVB_enabled": False, "bh.seeding.pop3.enabled": False, "bh.seeding.direct_collapse.enabled": False,
        "bh.seeding.halo_mass_threshold.enabled": False}
SK = {"disc.enabled": True, "disc.star_formation_law": "schmidt_kennicutt"}
LAWS = [("current (eps 0.015)", {}),
        ("current, eps 0.005", {"sf.efficiency": 0.005}),
        ("current, eps 0.05", {"sf.efficiency": 0.05}),
        ("current + threshold 5 Msun/pc^2", {"disc.enabled": True}),
        ("Schmidt-Kennicutt", SK),
        ("Schmidt-Kennicutt, A / 3", {**SK, "disc.schmidt_kennicutt_normalisation": 2.5e-4 / 3}),
        ("Schmidt-Kennicutt, A x 3", {**SK, "disc.schmidt_kennicutt_normalisation": 7.5e-4})]


def run(M, Mdot, z, **change):
    cfg = override(DEFAULT, **{**BASE, **change})
    return GridIntegrator(cfg).run(M, Mdot, z)[0], cfg


def summary(s):
    g, sfr = s.gas_mass[0], s.sfr[0]
    ok = (sfr > 0) & (g > 0)
    return dict(stars=s.stars_mass[0, -1], gas=g[-1], Z=s.gas_metals[0, -1] / g[-1] if g[-1] > 0 else np.nan,
                dep=np.median(g[ok] / sfr[ok]) if ok.any() else np.nan)


def best_constant_efficiency(M, Mdot, z, target):
    grid = np.logspace(-3, 0.5, 60)
    stars = np.array([run(M, Mdot, z, **{"sf.efficiency": e})[0].stars_mass[0, -1] for e in grid])
    return grid[np.argmin(np.abs(np.log(stars / target)))]


def main():
    print("1. Stellar mass, gas mass, gas metallicity at z = 5 and median depletion time, relative to the current law (median of 3 haloes)")
    print(f"{'M_h':>6} {'law':34s} {'M*':>6} {'Mgas':>6} {'Zgas':>6} {'t_dep (Gyr)':>12}")
    series = {}
    for mb in (1e8, 1e10, 1e11):
        hm, _, z = utils.read_trees(TREES, mb)
        rows = {name: [] for name, _ in LAWS}
        for i in range(3):
            M = hm[i:i + 1]; Mdot = np.maximum(np.gradient(M, utils.time_at_z(z), axis=1), 0.0)
            for name, ch in LAWS:
                s, cfg = run(M, Mdot, z, **ch)
                rows[name].append(summary(s))
                if name == "Schmidt-Kennicutt":
                    series[(mb, i)] = (M, Mdot, z, s, cfg)
        ref = rows[LAWS[0][0]]
        for name, _ in LAWS:
            r = rows[name]
            rel = lambda k: np.nanmedian([r[i][k] / ref[i][k] if ref[i][k] > 0 else np.nan for i in range(3)])
            print(f"{mb:6.0e} {name:34s} {rel('stars'):6.2f} {rel('gas'):6.2f} {rel('Z'):6.2f} {np.nanmedian([x['dep'] for x in r]):12.2f}")
        print()

    print("2. The efficiency per halo free-fall time that the Schmidt-Kennicutt law amounts to (5th, 50th, 95th percentile over time;")
    print("   and at z = 20, 14, 10, 8, 6 for the first halo), and what a single constant efficiency does")
    for mb in (1e8, 1e10, 1e11):
        eff_all, fits = [], []
        for i in range(3):
            M, Mdot, z, s, cfg = series[(mb, i)]
            ok = (s.gas_mass[0] > 0) & (s.sfr[0] > 0)
            eff = np.where(ok, s.sfr[0] / np.where(s.gas_mass[0] > 0, s.gas_mass[0], 1.0) * time_freefall(z), np.nan)
            eff_all.append(eff)
            e = best_constant_efficiency(M, Mdot, z, s.stars_mass[0, -1])
            f, _ = run(M, Mdot, z, **{"sf.efficiency": e})
            at = [int(np.argmin(abs(z - zz))) for zz in (20, 14, 10, 8, 6)]
            fits.append((e, f.gas_mass[0, -1] / s.gas_mass[0, -1], (f.gas_metals[0, -1] / f.gas_mass[0, -1]) / (s.gas_metals[0, -1] / s.gas_mass[0, -1]),
                         [f.sfr[0][k] / s.sfr[0][k] if s.sfr[0][k] > 0 else np.nan for k in at]))
        e0 = eff_all[0]
        z0 = series[(mb, 0)][2]
        at = [int(np.argmin(abs(z0 - zz))) for zz in (20, 14, 10, 8, 6)]
        print(f"\n  M_h = {mb:.0e}: implied eps  5/50/95 %: " + " | ".join(f"{np.nanpercentile(e, 5):.3f}/{np.nanmedian(e):.3f}/{np.nanpercentile(e, 95):.3f}" for e in eff_all))
        print(f"     halo 0 at z = 20, 14, 10, 8, 6: " + "  ".join(f"{e0[k]:.3f}" for k in at))
        print(f"     constant eps with the same final M*: " + ", ".join(f"{x[0]:.3f}" for x in fits) + "   -> gas at z=5 / SK: " + ", ".join(f"{x[1]:.2f}" for x in fits)
              + "   Zgas / SK: " + ", ".join(f"{x[2]:.2f}" for x in fits))
        print(f"     SFR / SK at z = 20, 14, 10, 8, 6 (halo 0): " + "  ".join(f"{x:.2f}" for x in fits[0][3]))


if __name__ == "__main__":
    main()
