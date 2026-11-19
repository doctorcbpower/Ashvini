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


def summarise(chain, tau, lower, upper, truth_variable):
    """
    Per parameter, from a chain of shape (steps, walkers, ndim) in sampling variables: the median and 68 and 95 per cent intervals, the posterior
    probability below the truth (its rank: uniform over many mocks if the fit is calibrated), whether the truth is inside each interval, and how far
    the posterior is from flat across the prior range (the largest departure of a tenth of the range from 0.1). Samples are thinned by the longest
    autocorrelation time before the deciles are counted, so that the numbers do not rest on strongly correlated samples.
    """
    step = max(1, int(np.ceil(np.max(tau))))
    thin = chain[::step].reshape(-1, chain.shape[-1])
    flat = chain.reshape(-1, chain.shape[-1])
    out = {}
    for i in range(chain.shape[-1]):
        q = np.quantile(flat[:, i], [0.025, 0.16, 0.5, 0.84, 0.975])
        deciles = np.histogram((thin[:, i] - lower[i]) / (upper[i] - lower[i]), bins=10, range=(0, 1))[0] / len(thin)
        rank = float(np.mean(flat[:, i] < truth_variable[i]))
        out[i] = dict(sampling_variable_quantiles=q.tolist(), truth_variable=float(truth_variable[i]), rank_of_truth=rank,
                      in_68=bool(q[1] <= truth_variable[i] <= q[3]), in_95=bool(q[0] <= truth_variable[i] <= q[4]),
                      width_68_dex=float(q[3] - q[1]), decile_mass=deciles.tolist(), max_decile_departure=float(np.max(np.abs(deciles - 0.1))),
                      thinned_samples=int(len(thin)))
    return out


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
    checkpoint = os.path.join(OUT, "checkpoints", f"case_{case}_n{n_trees}")     # a task that stops at the time limit resumes from here when resubmitted
    res = fit(problem, n_steps=steps, workers=workers, seed=case, checkpoint=checkpoint, log_every=50)
    v_true = problem.to_variables(truth_phys)
    out = dict(case=case, truth=dict(zip(problem.labels, truth_phys)), mock_from="independent trees" if independent else "same trees", n_trees=n_trees,
               steps=steps, workers=workers, seconds_setup=t1 - t0, seconds_fit_this_run=time.time() - t1, acceptance=res.acceptance, tau=res.tau.tolist(),
               r_hat=res.r_hat.tolist(), converged=res.converged, warnings=res.warnings, summary=res.summary(), params={})
    for i, d in summarise(res.chain, res.tau, problem.lower, problem.upper, v_true).items():
        out["params"][problem.labels[i]] = d
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
    """A table of the finished cases. The statistics are recomputed from the stored chains, so cases run with an older version of this script are included."""
    from ashvini.fit import S_INT_RANGE
    model = Model("fiducial", free=FREE)
    lower = np.append(np.log10(model.bounds[:, 0]), np.log10(S_INT_RANGE[0]))        # all free knobs here are log parameters
    upper = np.append(np.log10(model.bounds[:, 1]), np.log10(S_INT_RANGE[1]))
    rows = []
    for p in sorted(glob.glob(os.path.join(OUT, "case_*.json"))):
        r = json.load(open(p))
        rows.append(f"case {r['case']} ({r['mock_from']}), {r['steps']} steps, {r['seconds_fit_this_run'] / 60:.0f} min this run, converged={r['converged']}, "
                    f"tau up to {max(r['tau']):.0f} steps")
        chain_file = p.replace(".json", "_chain.npz")
        names = list(r["params"])
        truth_v = np.array([r["params"][n]["truth_variable"] for n in names])
        stats_ = summarise(np.load(chain_file)["chain"], r["tau"], lower, upper, truth_v) if os.path.exists(chain_file) else None
        for i, name in enumerate(names):
            d = stats_[i] if stats_ else r["params"][name]
            rows.append(f"    {name:22s} truth {r['truth'][name]:9.3g}  rank {d.get('rank_of_truth', float('nan')):.2f}  in 68%: {d['in_68']!s:5}  in 95%: {d['in_95']!s:5}  "
                        f"68% width {d.get('width_68_dex', float('nan')):.2f} dex  largest departure from flat {d.get('max_decile_departure', float('nan')):.3f}")
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
