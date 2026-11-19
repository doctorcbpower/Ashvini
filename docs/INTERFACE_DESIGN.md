# Ashvini: a scaling-relation interface (design sketch, nothing implemented)

Goal: going from "I want the BH mass–stellar mass relation at z = 0 and I have data" to a fitted set of
parameters should take a few lines, and every line should be about the physics, not the plumbing.

Status: sketch for comment. Names are proposals. Nothing here exists in the code except `run_forest`.

## 1. The workflow in five lines

```python
import ashvini as av

sample = av.TreeSample.for_stellar_masses(logMstar=(8.5, 11.5), z=0.0, n_haloes=2000, seed=1)   # step 1
model  = av.Model(preset="fiducial", free=["sf_efficiency", "wind_loading", "bh_accretion_eff", "bh_feedback_coupling"])  # step 2
rel    = model.relation(sample, y="bh_mass", x="stellar_mass", z=0.0)                           # step 3
data   = av.Data.from_table("reines_volonteri_2015.csv", x="logMstar", y="logMBH", yerr="err", upper_limits="is_limit")  # step 4
chain  = av.fit(model, sample, data, nwalkers=32, nsteps=3000)                                  # step 5
chain.corner(); chain.overlay(data)
```

Each step is a separate object, so you can stop after any of them (for example, after step 3, plot the model
relation with no data).

## 2. The five pieces

### 2.1 `TreeSample`: which haloes
Builds the set of merger trees whose galaxies end up spanning the stellar-mass range at the target redshift.

- **The scientific problem.** The sample is selected on stellar mass, but trees come labelled by halo mass.
  Selecting needs an assumed stellar-to-halo mass relation *including its scatter*. If the assumed relation is
  tighter than the model's, the fit sees a biased sample (haloes that the model would place in the range
  are missing, and vice versa).
- **Proposed treatment.** Draw halo masses from the halo mass function over a wide range (so that the range is
  covered *by the model*, not by the assumption), run, then select on the model's own M* afterwards and weight
  by the halo mass function. The assumed relation then only decides how far to extend the halo mass range,
  not the answer. This costs more haloes than selecting up front; the sample size is the price.
- Trees are cached to disk with a hash of (halo mass range, z, resolution, seed), so the same call twice
  loads, not rebuilds.
- Option to reuse a fixed realisation (deterministic likelihood, which emcee needs) or to draw a fresh one.

### 2.2 `Model`: the physics knobs
A thin, documented layer over `override(DEFAULT, ...)`.

- A **parameter registry**: one table, one line per knob, giving the physical name, config field(s) it sets,
  units, default, a sensible range, and a one-sentence description. `av.parameters()` prints it.
- **Presets** (`fiducial`, `momentum_wind`, `energy_wind`, `grumpy_like`) fix everything not listed in `free`.
- Parameters that map to Arepo Solas sub-grid parameters carry a `solas_name` entry in the registry, so the
  bridge to the simulations is explicit and visible.
- Free parameters are sampled in log where they span decades; the registry says which.
- Unknown names raise an error that lists the closest valid names.

### 2.3 `Model.relation(...)`: the model's prediction
Runs the sample and reduces to a relation.

- Returns median, 16th and 84th percentiles of y in bins of x (or, with `scatter="samples"`, the individual
  points), plus the number of haloes per bin and the fraction of haloes with no black hole.
- **Zero-BH haloes are reported, never dropped silently.** Whether they enter the median is a choice the user
  makes with `occupied_only=True/False`, because it changes the answer a great deal at low mass.
- Quantities available: anything in the `run_forest` output (`stars_mass`, `bh_mass`, `gas_mass`, `sfr`, ...),
  at any redshift, plus ratios (`sSFR`, `bh_mass/stellar_mass`).
- Returns an object that plots itself: `rel.plot()`.

### 2.4 `Data`: the observations
- `Data.from_table(...)` reads a CSV or a pandas frame; columns are named by the user.
- Carries uncertainties on both axes where available, upper limits, and (optionally) a selection note.
- **Selection of the observed sample is the data's job to declare.** A BH sample selected on BH detectability
  is biased high in M_BH at fixed M*; the likelihood needs to know. See 4.

### 2.5 `fit(...)`: the sampler
- Default sampler emcee; the log-probability is built from the model relation and the data, and is
  accessible (`av.logprob(theta, model, sample, data)`) so any other sampler can be used.
- Priors come from the registry ranges (uniform in log or linear as flagged), overridable.
- Returns a `Chain` with `corner()`, `overlay(data)` (posterior draws of the relation over the data),
  `table()` (median and 68% interval), and `save()/load()`.
- **Built-in recovery test**: `av.fit(..., mock=True)` generates mock data from the model at known
  parameters, fits it, and reports whether the truth is recovered. This is the first thing to run on any new
  relation, and it takes the same two lines.

## 3. Things that must be true for this to be pleasant

1. **Speed.** Currently about 10 s (estimated, not timed) per evaluation of 900 haloes. For interactive
   use we need either fewer haloes during exploration (`n_haloes=200`, with the noise it implies), or an
   emulator. Proposal: `fit(..., emulate=True)` trains a Gaussian process on a short space-filling design
   of runs, then samples the cheap surrogate; it reports the emulator's own error against held-out runs.
2. **Reproducibility.** One `seed` argument controls trees, sampling and walkers. A run records the full
   resolved configuration, git hash, and registry version with its results.
3. **No environment variables.** Output location is an argument with a default (`./ashvini_output`).
4. **OzSTAR.** `fit(..., backend="slurm")` or a CLI wrapper writes the job script and collects chains; the
   interactive form is the same call with the default local backend. Only the small tests run locally.
5. **Honest defaults.** If a parameter is unconstrained by the data, the summary says so (posterior ≈
   prior), and does not report a best-fit number as though it meant something.

## 4. Decisions

All seven were agreed as recommended (3 Nov 2026). The "My recommendation" column is now the decision.

| # | Question | My recommendation |
|---|----------|-------------------|
| 1 | Sample selection: select on model M* afterwards (weighted) or on an assumed relation up front? | Afterwards and weighted. More haloes, but the fit does not inherit an assumed SHMR. |
| 2 | Selection effects in the data (BH detectability, dwarf AGN completeness): model them in the likelihood, or restrict to a mass range where they are small? | Start with a restricted range and a declared caveat. A selection model is a second project. |
| 3 | Zero-BH haloes: part of the relation or reported separately? | Separately, with an explicit occupation fraction; the fit can then use both. |
| 4 | Scatter: fit the median only, or the median and scatter? | Both. The scatter from assembly history alone (0.3–0.5 dex at low mass) is a real constraint and an easy way to learn that the model needs more intrinsic scatter. |
| 5 | Which knobs are in the registry for v1? | SF efficiency, wind loading amplitude and exponents, reionisation redshift, BH seed mass, BH accretion efficiency, BH feedback coupling. Everything else stays as a preset. |
| 6 | Emulator in v1 or after a first direct emcee run? | After. Get the direct version working and timed first. |
| 7 | Redshift: z = 0 only for v1? | Yes. Higher z needs a verified z-dependent SHMR for the tree selection and observed relations with their own selection. |

## 5. What this does not do

- It does not make the model calibrated. It makes testing an idea and measuring how constrained it is
  quick, which is the stated purpose.
- It does not fix known weaknesses: the z = 5 stellar-mass offset, the assembly-only scatter being larger than
  published limits at low mass, the metallicity step as a hidden assumption. A fit will absorb these into
  the parameters, so posteriors should be read with that in mind.
- The BH parameters had zero Sobol index for M* in the earlier exploration. For M_BH they may be
  degenerate with each other (accretion efficiency vs feedback coupling vs seed mass). A one-at-a-time check
  of which parameters move M_BH at fixed M* is the first science step, before any sampling.

## 6. Suggested order, if agreed

1. Parameter registry and `Model` (no new physics; documents and validates what exists).
2. `TreeSample` with the afterwards-selection and caching.
3. `relation()` and plotting.
4. Mock-recovery test using steps 1–3.
5. `Data` and `fit` with emcee.
6. Emulator and the OzSTAR backend, once timings are known.

Each step is testable on its own, and steps 1–3 already give a useful tool without any fitting.
