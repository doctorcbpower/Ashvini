"""
Analyse the parameter-space samples (space_eval.py) as set out in docs/PARAMETER_SPACE_AUDIT.md:
  1. Sobol first-order and total indices of log10 M*(z = 0) at each halo mass, with bootstrap 95 percent half-widths;
  2. the samples whose median stellar mass is consistent with both published relations at the masses where they are constrained
     (3e10, 1e11, 3e11, 1e12 Msun): within TOLERANCE_DEX of the interval between Moster et al. (2013) and Behroozi et al. (2013);
  3. what those accepted samples predict at 1e8 - 3e9 Msun, where the relations are extrapolations.
Usage: python scripts/space_analyse.py
"""
import os

import numpy as np

import shmr_z0_common as C
import space_common as S
from stellar_to_halo_mass_relation import behroozi2013_mstar, moster2013_mstar

TOLERANCE_DEX = 0.2
CONSTRAINED = [5, 6, 7, 8]          # 3e10, 1e11, 3e11, 1e12 Msun
LOW = [0, 1, 2, 3, 4]               # 1e8 ... 1e10 Msun
PLOT = os.path.join(C.OUT, "space")


def load_subspace(sub):
    """Assembled results, shape (n_samples, n_masses, n_fields); samples not yet evaluated are NaN."""
    n = S.n_samples(sub)
    out = np.full((n, len(C.MASSES), len(S.FIELDS)), np.nan)
    for name in sorted(os.listdir(PLOT)) if os.path.isdir(PLOT) else []:
        if name.startswith(sub + "_shard") and name.endswith(".npz") and not name.endswith(".tmp.npz"):
            d = np.load(os.path.join(PLOT, name))
            ok = d["done"]
            out[d["samples"][ok]] = d["results"][ok]
    return out


def published(mass):
    b, m = float(behroozi2013_mstar(mass)), float(moster2013_mstar(mass))
    return min(b, m), max(b, m)


def accepted_mask(res):
    ok = np.ones(len(res), dtype=bool)
    for k in CONSTRAINED:
        lo, hi = published(C.MASSES[k])
        logm = np.log10(np.maximum(res[:, k, 0], 1.0))
        ok &= (logm >= np.log10(lo) - TOLERANCE_DEX) & (logm <= np.log10(hi) + TOLERANCE_DEX)
    return ok


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lines = [f"Parameter-space exploration; tolerance {TOLERANCE_DEX} dex around the interval between Moster+13 and Behroozi+13 at "
             + ", ".join(f"{C.MASSES[k]:.0e}" for k in CONSTRAINED) + " Msun", ""]
    fig_s, axs = plt.subplots(1, len(S.SUBSPACES), figsize=(5 * len(S.SUBSPACES), 4.5))
    fig_p, axp = plt.subplots(1, len(S.SUBSPACES), figsize=(5 * len(S.SUBSPACES), 4.5), sharey=True)
    for ax_s, ax_p, sub in zip(np.atleast_1d(axs), np.atleast_1d(axp), S.SUBSPACES):
        names = S.SUBSPACES[sub][1]
        res = load_subspace(sub)
        done = np.isfinite(res[:, 0, 0])
        lines.append(f"== {sub}: {done.sum()} of {len(res)} samples evaluated")
        if done.sum() < len(res):
            lines += ["   incomplete: indices and acceptance are only computed for complete sub-spaces", ""]
            continue
        bad = int(np.nansum(res[:, :, 8]))
        if bad:
            lines.append(f"   WARNING: {bad} non-finite values in the evolved arrays")
        N, D = S.N_BASE, len(names)
        ST = np.zeros((D, len(C.MASSES)))
        for k in range(len(C.MASSES)):
            f = np.log10(np.maximum(np.nan_to_num(res[:, k, 0], nan=1.0), 1.0))
            s1, st = S.sobol_indices(f, N, D)
            w1, wt = S.bootstrap_indices(f, N, D)
            ST[:, k] = st
            if k in (2, 4, 6, 8):
                lines.append(f"   M_h = {C.MASSES[k]:.0e}: " + "; ".join(f"{n} S1 {a:+.2f}+-{b:.2f} ST {c:.2f}+-{d:.2f}"
                                                                         for n, a, b, c, d in zip(names, s1, w1, st, wt)))
        im = ax_s.imshow(ST, aspect="auto", vmin=0, vmax=1, cmap="viridis")
        ax_s.set_xticks(range(len(C.MASSES)), [f"{m:.0e}" for m in C.MASSES], rotation=60, fontsize=7)
        ax_s.set_yticks(range(D), names, fontsize=8); ax_s.set_title(f"{sub}: total Sobol index of log M*", fontsize=9)
        ok = accepted_mask(res)
        lines.append(f"   accepted: {ok.sum()} of {len(res)} ({100 * ok.mean():.1f} percent)")
        mh = np.logspace(8, 12, 100)
        ax_p.loglog(mh, np.array([float(behroozi2013_mstar(m)) for m in mh]) / mh, "k-", label="Behroozi+13")
        ax_p.loglog(mh, np.array([float(moster2013_mstar(m)) for m in mh]) / mh, "k--", label="Moster+13")
        if ok.any():
            band = np.array([np.percentile(res[ok, k, 0], [16, 50, 84]) for k in range(len(C.MASSES))]) / C.MASSES[:, None]
            ax_p.fill_between(C.MASSES, np.maximum(band[:, 0], 1e-9), np.maximum(band[:, 2], 1e-9), alpha=0.25, color="C0")
            ax_p.plot(C.MASSES, np.maximum(band[:, 1], 1e-9), "o-", color="C0", label="accepted samples")
            for k in LOW:
                lo, hi = published(C.MASSES[k])
                lines.append(f"   predicted M*/M_halo at {C.MASSES[k]:.0e}: 16/50/84 percent {band[k, 0]:.1e} / {band[k, 1]:.1e} / {band[k, 2]:.1e}"
                             f"  (published: Moster {float(moster2013_mstar(C.MASSES[k])) / C.MASSES[k]:.1e}, Behroozi {float(behroozi2013_mstar(C.MASSES[k])) / C.MASSES[k]:.1e})")
            for n_, vals in zip(names, np.array([[S.physical(sub, u)[n] for n in names] for u in S.design(sub)[ok]]).T):
                lines.append(f"   accepted {n_}: 16/50/84 percent {np.percentile(vals, 16):.3g} / {np.percentile(vals, 50):.3g} / {np.percentile(vals, 84):.3g}")
        ax_p.set_title(f"{sub}", fontsize=9); ax_p.set_xlabel(r"$M_{\rm halo}$ [Msun]"); ax_p.set_ylim(1e-6, 0.1); ax_p.legend(fontsize=7)
        lines.append("")
    os.makedirs(PLOT, exist_ok=True)
    fig_s.tight_layout(); fig_s.savefig(os.path.join(PLOT, "space_sobol_total.png"), dpi=120)
    fig_p.tight_layout(); fig_p.savefig(os.path.join(PLOT, "space_accepted_band.png"), dpi=120)
    open(os.path.join(PLOT, "space_summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
