"""
Plot the z = 0 stellar-to-halo mass relation of each variant against the Moster et al. (2013) and Behroozi et al. (2013) relations,
from the stored results, and write a table of the median stellar mass relative to Behroozi et al. (2013) at each halo mass.

Both published relations are fits to galaxies above about 1e10 Msun in halo mass and are extrapolated below: treat the
comparison at 1e8-1e10 Msun as a comparison with an extrapolation. Ashvini's haloes are isolated main branches (no
satellites). Nothing has been tuned to these relations.

Usage: python scripts/shmr_z0_plot.py
"""
import os

import numpy as np

import shmr_z0_common as C
from stellar_to_halo_mass_relation import behroozi2013_mstar, moster2013_mstar


def load(variant):
    rows = []
    for i in range(len(C.MASSES)):
        p = C.result_path(variant, i)
        if os.path.exists(p):
            d = np.load(p)
            rows.append((float(d["mass_bin"]), d["halo_mass_z0"], d["stars_mass"], d["gas_mass"], d["gas_metals"], bool(d["fast_star_formation_warning"])))
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    variants = [v for v in C.VARIANTS if load(v)]
    fig, axes = plt.subplots(1, 4, figsize=(22, 5), sharey=True)
    groups = [("star formation efficiency of the current law", ["current_eps0.005", "current_eps0.015", "current_eps0.05"]),
              ("Schmidt-Kennicutt law, disc size", ["sk_chi0.025", "sk_chi0.06", "sk_chi0.06_small", "sk_chi0.06_large", "sk_chi0.06_slope1"]),
              ("equilibrium UV suppression (reion.form)", ["eq_eps0.005", "eq_eps0.015", "eq_eps0.05", "eq_sk_chi0.06", "eq_sk_chi0.06_A/3", "eq_zrei6", "eq_zrei8"]),
              ("reionisation (current law unless stated)", ["current_zrei6", "current_eps0.015", "current_zrei8", "current_noUV", "sk_chi0.06_noUV"])]
    mh = np.logspace(8, 12, 100)
    lines = []
    header = "variant".ljust(24) + "".join(f"{m:10.1e}" for m in C.MASSES)
    for ax, (title, names) in zip(axes, groups):
        ax.loglog(mh, behroozi2013_mstar(mh) / mh, "k-", lw=2, label="Behroozi+13")
        ax.loglog(mh, moster2013_mstar(mh) / mh, "k--", lw=2, label="Moster+13")
        for v in names:
            rows = load(v)
            if not rows:
                continue
            x = [np.median(r[1]) for r in rows]
            y = [np.median(r[2]) / np.median(r[1]) if np.median(r[2]) > 0 else np.nan for r in rows]
            lo = [np.percentile(r[2] / r[1], 16) for r in rows]; hi = [np.percentile(r[2] / r[1], 84) for r in rows]
            line, = ax.loglog(x, y, "o-", label=v)
            ax.fill_between(x, np.maximum(lo, 1e-9), np.maximum(hi, 1e-9), alpha=0.15, color=line.get_color())
        ax.set(xlabel=r"$M_{\rm halo}(z=0)\ [{\rm M}_\odot]$", title=title, ylim=(1e-6, 0.1))
        ax.legend(fontsize=7)
    axes[0].set_ylabel(r"$M_\star/M_{\rm halo}$")
    fig.tight_layout()
    os.makedirs(C.OUT, exist_ok=True)
    fig.savefig(os.path.join(C.OUT, "shmr_z0_laws.png"), dpi=130)
    for name, relation in (("Behroozi et al. (2013)", behroozi2013_mstar), ("Moster et al. (2013)", moster2013_mstar)):
        lines += ["", f"median stellar mass relative to {name}", header]
        lines.append("published M* (Msun)".ljust(24) + "".join(f"{float(relation(m)):10.1e}" for m in C.MASSES))
        for v in variants:
            rows = load(v)
            text = v.ljust(24)
            for m in C.MASSES:
                r = next((rr for rr in rows if abs(np.log10(rr[0] / m)) < 1e-6), None)
                if r is None:
                    text += "       n/a"
                    continue
                med = np.median(r[2])
                text += f"{med / float(relation(np.median(r[1]))):10.2f}" if med > 0 else "      zero"
            lines.append(text)
    with open(os.path.join(C.OUT, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
