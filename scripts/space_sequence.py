"""
The stellar-to-halo mass relation at several redshifts for the accepted parameter sets of the exploration (docs/PARAMETER_SPACE_AUDIT.md).

For each equilibrium-UV sub-space, K accepted samples (chosen at random, seeded) are rerun on the stored trees and the halo and stellar masses of every
halo are stored at z = 0, 1, 2, 3, 5, 7, 10. The relation at redshift z is then M*/M_halo against M_halo(z) of the main progenitors, pooled over the
nine tree masses, the haloes and the samples.

Moster et al. (2013) at z > 0: the redshift dependence below uses coefficients typed from memory of their Table 1 and NOT checked against the paper
(the z = 0 values are the ones in stellar_to_halo_mass_relation.py). Both published relations are heavily extrapolated at z = 5.

Usage: python scripts/space_sequence.py [K per sub-space, default 12]
"""
import os
import sys
import time
import warnings

import numpy as np

import shmr_z0_common as C
import space_common as S
import space_analyse as A

warnings.filterwarnings("ignore")
SUBS = ["cur_ash", "cur_en", "sk_ash", "sk_en"]
ZS = [0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]


def moster_z(mh, z):
    """Moster, Naab & White (2013) with the redshift evolution z/(z+1) of M1, N, beta, gamma. UNVERIFIED coefficients (see module docstring)."""
    f = z / (z + 1.0)
    M1 = 10 ** (11.590 + 1.195 * f)
    N = 0.0351 - 0.0247 * f
    beta = 1.376 - 0.826 * f
    gamma = 0.608 + 0.329 * f
    x = mh / M1
    return mh * 2 * N / (x ** (-beta) + x ** gamma)


def run_selection(K, seed=7):
    from ashvini.config import DEFAULT, override
    from ashvini.main import run_forest
    trees = [dict(np.load(C.tree_path(i))) for i in range(len(C.MASSES))]
    rng = np.random.default_rng(seed)
    for sub in SUBS:
        path = os.path.join(C.OUT, "space", f"sequence_{sub}.npz")
        if os.path.exists(path):
            print(sub, "exists, skipped"); continue
        res = A.load_subspace(sub)
        acc = np.flatnonzero(A.accepted_mask(res))
        pick = np.sort(rng.choice(acc, size=min(K, len(acc)), replace=False))
        unit = S.design(sub)
        stars = np.zeros((len(pick), len(C.MASSES), len(ZS), len(trees[0]["halo_masses"])))
        halo = np.zeros_like(stars)
        t0 = time.time()
        for n, s in enumerate(pick):
            cfg = override(DEFAULT, **S.config_changes(sub, S.physical(sub, unit[s])))
            for i, t in enumerate(trees):
                r = run_forest(t["halo_masses"], t["halo_growth_rates"], t["redshifts"], config=cfg)
                for j, z in enumerate(ZS):
                    k = int(np.argmin(abs(t["redshifts"] - z)))
                    stars[n, i, j] = r["stars_mass"][:, k]
                    halo[n, i, j] = t["halo_masses"][:, k]
            print(f"{sub}: {n + 1}/{len(pick)} samples, {time.time() - t0:.0f} s", flush=True)
        np.savez(path, samples=pick, zs=ZS, stars=stars, halo=halo)


def pooled_relation(sub, zi, edges):
    d = np.load(os.path.join(C.OUT, "space", f"sequence_{sub}.npz"))
    mh = d["halo"][:, :, zi, :].reshape(-1)
    ms = d["stars"][:, :, zi, :].reshape(-1)
    ok = mh > 0
    mh, ms = mh[ok], ms[ok]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (mh >= lo) & (mh < hi)
        if sel.sum() < 30:
            out.append((np.sqrt(lo * hi), np.nan, np.nan, np.nan, int(sel.sum()))); continue
        r = ms[sel] / mh[sel]
        out.append((np.exp(np.mean(np.log(mh[sel]))), *np.percentile(r, [16, 50, 84]), int(sel.sum())))
    return out


def main(K):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    run_selection(K)
    edges = 10 ** np.arange(5.5, 12.01, 0.5)
    show = [0.0, 2.0, 5.0]
    fig, axes = plt.subplots(1, len(show), figsize=(6 * len(show), 5), sharey=True)
    lines = ["SHMR of the accepted parameter sets (equilibrium UV form), pooled over haloes and samples; M*/M_halo (16/50/84 percent) in bins of M_halo(z)",
             "Moster+13 for z > 0 uses coefficients typed from memory and not checked against the paper"]
    for ax, z in zip(axes, show):
        zi = ZS.index(z)
        mhgrid = np.logspace(5.5, 12, 100)
        ax.loglog(mhgrid, moster_z(mhgrid, z) / mhgrid, "k--", label=f"Moster+13 at z={z:g} (z > 0 coefficients unchecked)")
        if z == 0:
            from stellar_to_halo_mass_relation import behroozi2013_mstar
            ax.loglog(mhgrid, behroozi2013_mstar(mhgrid) / mhgrid, "k-", label="Behroozi+13 at z=0")
        for sub, color in zip(SUBS, ["C0", "C1", "C2", "C3"]):
            rel = np.array([r[:4] for r in pooled_relation(sub, zi, edges)])
            ok = np.isfinite(rel[:, 2]) & (rel[:, 2] > 0)
            ax.loglog(rel[ok, 0], rel[ok, 2], "o-", color=color, label=sub)
            ax.fill_between(rel[ok, 0], rel[ok, 1], rel[ok, 3], color=color, alpha=0.12)
            if z in (5.0, 0.0):
                lines.append(f"z = {z:g}, {sub}: " + "; ".join(f"M_h {r[0]:.1e}: {r[2]:.1e} ({r[1]:.1e}-{r[3]:.1e}, n={r[4]})" for r in pooled_relation(sub, zi, edges) if np.isfinite(r[2])))
        if z == 5.0:
            lines.append("  Moster+13 (unchecked z coefficients) at those masses: " + ", ".join(f"{m:.1e}: {moster_z(m, 5.0) / m:.1e}" for m in 10 ** np.arange(6, 11.6, 1.0)))
        ax.set(xlabel=r"$M_{\rm halo}(z)$ [Msun]", title=f"z = {z:g}", ylim=(1e-6, 0.1)); ax.legend(fontsize=7)
    axes[0].set_ylabel(r"$M_\star/M_{\rm halo}$")
    fig.tight_layout()
    fig.savefig(os.path.join(C.OUT, "space", "space_shmr_redshift.png"), dpi=120)
    open(os.path.join(C.OUT, "space", "space_shmr_redshift.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
