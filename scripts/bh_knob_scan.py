"""
One-at-a-time scan of the knobs against M_BH at fixed M* (docs/FIT_AUDIT.md, section 8, step 1).

Each knob is set to the low end, the preset value and the high end of its registry range, everything else at the fiducial preset, on the same tree sample.
For the bins of M* that are not flagged, prints the change in the median log10 M_BH and in log10 M* (via the relation itself) and the occupation fraction.
Usage: python scripts/bh_knob_scan.py [cache_dir]
"""
import sys
import warnings

import numpy as np

import ashvini as av
from ashvini.parameters import KNOBS, knob_range

warnings.filterwarnings("ignore")
CACHE = sys.argv[1] if len(sys.argv) > 1 else "cache"
BINS = np.arange(8.75, 10.3, 0.25)
KNOBS_SCANNED = ["bh_accretion_eff", "bh_feedback_coupling", "bh_seed_mass_pop3", "bh_seed_mass_direct_collapse", "bh_seed_mass_threshold",
                 "sf_efficiency", "wind_loading", "z_reion"]


def run(sample, name=None, value=None):
    model = av.Model("fiducial", free=[name] if name else [])
    return model.relation(sample, [value] if name else None, y="bh_mass", x="stars_mass", bins=BINS, n_boot=0, min_ess=10)


if __name__ == "__main__":
    sample = av.TreeSample.for_stellar_masses((8.0, 10.5), n_trees=1000, seed=1, cache_dir=CACHE)
    base = run(sample)
    ok = base.usable & np.isfinite(base.median)
    print("bins (log M*):", np.round(base.x_mean[ok], 2))
    print("fiducial log M_BH:", np.round(base.median[ok], 2), " f_occ:", np.round(base.f_occ[ok], 2))
    for name in KNOBS_SCANNED:
        k = KNOBS[name]
        lo, hi, _ = knob_range(name)
        default = k.default(av.Model("fiducial").base)
        for label, v in (("low", lo), ("high", hi)):
            r = run(sample, name, v)
            good = ok & np.isfinite(r.median)
            d_bh = r.median[good] - base.median[good]
            print(f"{name:30s} {label:4s} {v:9.3g} (default {default:9.3g}): d log M_BH at fixed M* = {np.round(d_bh, 2)}  d scatter = {np.round((r.scatter - base.scatter)[good], 2)}"
                  f"  f_occ min {np.nanmin(r.f_occ[good]):.2f}")
    # how far down must the accretion efficiency go to reach M_BH/M* of 0.1 per cent?
    print("\nbh_accretion_eff beyond the registry range (log M_BH - log M* at the middle bin):")
    mid = int(np.argmin(abs(base.x_mean - 9.5)))
    for v in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        r = run(sample, "bh_accretion_eff", v)
        print(f"  {v:8.0e}: {r.median[mid] - r.x_mean[mid]:6.2f}  (f_occ {r.f_occ[mid]:.2f})")

    # other BH options (fixed values, not registry knobs), each against the fiducial preset
    print("\nother options: slope of log M_BH against log M* over the unflagged bins, offset at log M* = 9.5, scatter")

    def summary(label, model):
        r = model.relation(sample, None, y="bh_mass", x="stars_mass", bins=BINS, n_boot=0, min_ess=10)
        good = r.usable & np.isfinite(r.median)
        if not good.any():
            print(f"  {label:55s} no bin has enough trees (M* is far from the fiducial range)"); return
        x, y = r.x_mean[good], r.median[good]
        print(f"  {label:55s} slope {np.polyfit(x, y, 1)[0]:5.2f}  offset {np.interp(9.5, x, y) - 9.5:6.2f}  scatter {np.nanmedian(r.scatter[good]):.2f}")

    summary("fiducial", av.Model("fiducial"))
    for path, v in (("bh.eddington_multiplier", 0.3), ("bh.eddington_multiplier", 3.0), ("bh.feedback_delay_time", 0.3), ("bh.sigma_feedback.enabled", True)):
        summary(f"{path} = {v}", av.Model("fiducial", fixed={path: v}))
    summary("sigma_feedback enabled + growth cap", av.Model("fiducial", fixed={"bh.sigma_feedback.enabled": True, "bh.sigma_feedback.growth_cap_enabled": True}))
    summary("slim_disc preset (hobbs_slimdisk at its defaults)", av.Model("slim_disc"))

    # the knobs added after the first scan: the M-sigma growth cap's baryon fraction and the Eddington multiplier
    print("\nM-sigma cap and Eddington multiplier (slope, log M_BH/M* at log M* = 9.5, scatter):")

    def summary_theta(label, model, theta):
        r = model.relation(sample, theta, y="bh_mass", x="stars_mass", bins=BINS, n_boot=0, min_ess=10)
        good = r.usable & np.isfinite(r.median)
        x, y = r.x_mean[good], r.median[good]
        print(f"  {label:48s} slope {np.polyfit(x, y, 1)[0]:5.2f}  {np.interp(9.5, x, y) - 9.5:6.2f}  scatter {np.nanmedian(r.scatter[good]):.2f}")

    for f in (0.05, 0.16, 0.3):
        summary_theta(f"m_sigma_cap, f_g = {f}", av.Model("m_sigma_cap", free=["bh_msigma_baryon_fraction"]), [f])
    for e in (1e-5, 1e-3):
        summary_theta(f"m_sigma_cap, accretion eff = {e:g}", av.Model("m_sigma_cap", free=["bh_accretion_eff"]), [e])
    for e in (0.1, 0.3, 0.5, 1.0, 10.0):
        summary_theta(f"fiducial, Eddington multiplier = {e:g}", av.Model("fiducial", free=["bh_eddington_multiplier"]), [e])
    summary_theta("fiducial, accretion eff = 2e-5", av.Model("fiducial", free=["bh_accretion_eff"]), [2e-5])
