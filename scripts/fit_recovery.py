"""
Mock recovery on the real Model (docs/FIT_AUDIT.md, section 7 test 4 and section 14): can `fit` recover known parameters from the model's own relation?

    python scripts/fit_recovery.py run CASE [--n-trees 1000] [--steps 4000] [--workers 8]     # CASE 0 to 5; writes <out>/fit_recovery/case_<CASE>.json
    python scripts/fit_recovery.py trees [--n-trees 1000]                                     # only builds and stores the two tree samples (no fit)
    python scripts/fit_recovery.py report                                                      # table of the finished cases

Output goes to $SHMR_OUT, else /fred/oz009/cpower/Ashvini if that directory exists (OzSTAR), else scripts/output/shmr_z0; trees are cached in
<out>/tree_cache.
Free knobs: bh_accretion_eff (the one a M_BH-M* relation constrains, with the star formation efficiency fixed at the preset) and bh_feedback_coupling
(eta_agn, which section 10 found to change nothing, so it should come back as its prior). Plus the intrinsic scatter s_int.
Cases 0-2: truths in the preset's neighbourhood, at the preset, and near the lower bound of the accretion efficiency. Cases 3-5: the same truths with
the mock generated from an independent tree realisation (seed 2) while the fit uses seed 1 (design choice A4).
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np
from scipy import stats

from ashvini.fit import Curve, Data, fit, problem_for_model
from ashvini.model import Model
from ashvini.tree_sample import TreeSample

CLUSTER_OUT = "/fred/oz009/cpower/Ashvini"          # OzSTAR project storage: used when it exists, unless SHMR_OUT says otherwise
LOCAL_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "shmr_z0")
BASE = os.environ.get("SHMR_OUT") or (CLUSTER_OUT if os.path.isdir(CLUSTER_OUT) else LOCAL_OUT)
OUT = os.path.join(BASE, "fit_recovery")
FREE = ["bh_accretion_eff", "bh_feedback_coupling"]
TRUTHS = [dict(bh_accretion_eff=2e-4, bh_feedback_coupling=0.5), dict(bh_accretion_eff=1e-3, bh_feedback_coupling=0.1),
          dict(bh_accretion_eff=3e-6, bh_feedback_coupling=1.0)]
S_INT_TRUE, SIGMA_Y, N_POINTS = 0.2, 0.3, 40            # docs/FIT_AUDIT.md, I2
LOG_MSTAR = (8.0, 10.5)


def make_mock(model, sample, truth, seed):
    """40 points uniform in log M* over the unflagged bins of the model's relation at `truth`, with measurement error 0.3 dex and intrinsic scatter 0.2."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        rel = model.relation(sample, truth, n_boot=0)
    curve = Curve.from_relation(rel, median_err=np.zeros_like(rel.median))
    if len(curve.x) < 2:
        raise SystemExit(f"the relation at the truth has {len(curve.x)} unflagged bins: use more trees (--n-trees; 1000 gives about 6 bins).")
    rng = np.random.default_rng(seed)
    x = rng.uniform(curve.x[0], curve.x[-1], N_POINTS)
    mu, scatter, _, _ = curve.at(x)
    y = mu + rng.normal(0, np.sqrt(scatter ** 2 + S_INT_TRUE ** 2), N_POINTS) + rng.normal(0, SIGMA_Y, N_POINTS)
    return Data(x, y, np.full(N_POINTS, SIGMA_Y)), curve


def run(case, n_trees, steps, workers, cache):
    truth = TRUTHS[case % 3]
    independent = case >= 3
    t0 = time.time()
    kw = dict(n_trees=n_trees, cache_dir=cache)
    sample = TreeSample.for_stellar_masses(LOG_MSTAR, seed=1, **kw)
    sample_mock = TreeSample.for_stellar_masses(LOG_MSTAR, seed=2, **kw) if independent else sample
    model = Model("fiducial", free=FREE)
    data, true_curve = make_mock(model, sample_mock, truth, seed=100 + case)
    problem = problem_for_model(model, sample, data, n_boot=100)
    truth_phys = [truth[n] for n in FREE] + [S_INT_TRUE]
    t1 = time.time()
    res = fit(problem, n_steps=steps, workers=workers, seed=case)
    v_true = problem.to_variables(truth_phys)
    flat = res.chain.reshape(-1, res.chain.shape[-1])
    out = dict(case=case, truth=dict(zip(problem.labels, truth_phys)), mock_from="independent trees" if independent else "same trees", n_trees=n_trees,
               steps=steps, workers=workers, seconds_setup=t1 - t0, seconds_fit=time.time() - t1, acceptance=res.acceptance, tau=res.tau.tolist(),
               r_hat=res.r_hat.tolist(), converged=res.converged, warnings=res.warnings, summary=res.summary(), params={})
    for i, name in enumerate(problem.labels):
        lo, hi = problem.lower[i], problem.upper[i]
        q = np.quantile(flat[:, i], [0.0015, 0.025, 0.16, 0.5, 0.84, 0.975, 0.9985])
        u = (flat[:, i] - lo) / (hi - lo)
        out["params"][name] = dict(sampling_variable_quantiles=q.tolist(), truth_variable=float(v_true[i]), in_68=bool(q[2] <= v_true[i] <= q[4]),
                                   in_95=bool(q[1] <= v_true[i] <= q[5]), ks_uniform_p=float(stats.kstest(u[:: max(1, len(u) // 2000)], "uniform").pvalue),
                                   prior_fraction_in_central_half=float(np.mean((u > 0.25) & (u < 0.75))))
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(os.path.join(OUT, f"case_{case}_chain.npz"), chain=res.chain, log_prob=res.log_prob, labels=problem.labels)
    with open(os.path.join(OUT, f"case_{case}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(res.summary())
    return out


def build_trees(n_trees, cache):
    """Build and store the tree samples that every case uses (seed 1 for the fit, seed 2 for the independent mocks); no model runs."""
    for seed in (1, 2):
        s = TreeSample.for_stellar_masses(LOG_MSTAR, seed=seed, n_trees=n_trees, cache_dir=cache)
        print(f"seed {seed}: {s.n} trees, effective number {s.effective_n:.0f}")


def report():
    rows = []
    for p in sorted(glob.glob(os.path.join(OUT, "case_*.json"))):
        r = json.load(open(p))
        rows.append(f"case {r['case']} ({r['mock_from']}), {r['steps']} steps, {r['seconds_fit'] / 60:.0f} min, converged={r['converged']}")
        for name, d in r["params"].items():
            rows.append(f"    {name:22s} truth {r['truth'][name]:9.3g}  in 68%: {d['in_68']!s:5}  in 95%: {d['in_95']!s:5}  "
                        f"flat-posterior KS p {d['ks_uniform_p']:.3f}  mass in central half of range {d['prior_fraction_in_central_half']:.2f}")
        rows += [f"    WARNING: {w}" for w in r["warnings"]]
    print("\n".join(rows) or f"no finished cases in {OUT}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("case", type=int, choices=range(6))
    r.add_argument("--n-trees", type=int, default=1000)
    r.add_argument("--steps", type=int, default=4000)
    r.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    r.add_argument("--cache", default=os.environ.get("TREE_CACHE", os.path.join(BASE, "tree_cache")))
    t = sub.add_parser("trees")
    t.add_argument("--n-trees", type=int, default=1000)
    t.add_argument("--cache", default=r.get_default("cache"))
    sub.add_parser("report")
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.case, a.n_trees, a.steps, a.workers, a.cache)
    elif a.cmd == "trees":
        build_trees(a.n_trees, a.cache)
    else:
        report()
