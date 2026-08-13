# **Ashvini** galaxy formation and evolution model

Repo for the **Ashvini** galaxy formation and evolution model, a lightweight Python package for simulating galaxy formation and evolution in a cosmological context.

The model works by balancing the mass fluxes between the baryonic components of galaxies -- gas, stars, dust, and (optionally) a central black hole -- along dark matter halo merger histories. Given a halo's merger history, **Ashvini** models:

1. Star formation rate
1. Gas mass and gas-phase metallicity
1. Stellar mass and stellar metallicity
1. Dust mass
1. Black hole mass (seeding + Eddington-limited growth + AGN feedback)

## Installation

Clone the repository:

```
git clone https://github.com/Anand-JM97/Ashvini.git
cd Ashvini
pip install -r requirements.txt
pip install -e .
```

Required packages (all pip-installable): `numpy`, `scipy`, `astropy`, `h5py`, `pyyaml`.

For running the test suite, also install the dev requirements:

```
pip install -r requirements-dev.txt
```

## Usage

Configure a run by editing `run_params.yaml` (star formation efficiency, supernova feedback type/timing, reionization, metal yields, dust, black hole growth/seeding/AGN feedback -- see the comments in the file for each section), then run:

```
python3 run.py
```

This integrates every halo in the input merger-tree file and writes the results to an HDF5 file under `data/outputs/`.

### Input halo merger trees

**Ashvini** needs halo merger histories as input. You can download a merger tree file [here](https://drive.google.com/file/d/1eAiONNCOHSAw829n3zbR1izsIU6JNCO4/view?usp=sharing). This consists of merger histories of 100 haloes in each mass bin, where the bins are linearly distributed within $10^6\leq M_{\rm h}/M_\odot \leq 10^{11}$ at $z=5$.

Merger trees can also be generated directly with [pymctrees](https://github.com/doctorcbpower/pymctrees) and adapted to Ashvini's expected HDF5 layout with `scripts/build_trees_from_pymctrees.py`:

```
pip install -e /path/to/pymctrees[camb]   # pymctrees is not an Ashvini dependency
python scripts/build_trees_from_pymctrees.py \
    /path/to/pymctrees/config/planck2018.yml \
    data/inputs/merger_trees_pymctrees.h5 \
    --mass-bins 1e6 1e7 1e8 1e9 1e10 1e11 --n-halos 100 --z0 5.0 --z-max 20.0
```

See the script's docstring for the unit/ordering conversions it handles (pymctrees works in Msun/h and grows trees backward from z0; Ashvini expects plain Msun in forward chronological order).

Alternatively, `run()` can generate trees live -- no intermediate file -- by setting `basics.tree_source: pymctrees` in `run_params.yaml` instead of `file` (see the commented-out example in `run_params.yaml` for the full `basics.pymctrees` block: config path, halo count, z0/z_max/dz, mass resolution, backend, seed). Tree generation is fast enough that there's no caching -- every run regenerates the forest. For a forest you want to reuse across many `run_params.yaml` sweeps, build one offline instead with the script above and point `tree_file` at it.

## Package structure

| Module | Responsibility |
|---|---|
| `main.py` | Drives the integration: `run1()` (vectorised, fast path), `run1_scalar()` (solve_ivp reference), `run_forest()` (all haloes at once), `run()` (CLI entry point, reads or generates trees, writes HDF5) |
| `run_params.py` / `run_params.yaml` | Dataclass-based config loader |
| `pymctrees_adapter.py` | pymctrees <-> Ashvini tree-format conversion (unit/ordering conversions, pre-formation-prefix handling), shared by `scripts/build_trees_from_pymctrees.py` (offline) and `run()`'s live `tree_source: pymctrees` path |
| `gas_evolve.py` | Cosmological gas accretion and the gas-mass reservoir ODE |
| `star_formation.py` | Star formation rate |
| `metallicity.py` | Gas-phase and stellar metal enrichment |
| `dust.py` | Dust mass evolution |
| `supernovae_feedback.py` | Mass-loaded supernova winds |
| `black_holes_growth.py` | BH seeding (three configurable channels) and Eddington-limited growth (`black_holes.eddington_multiplier` allows super-Eddington) |
| `agn_feedback.py` | AGN-driven gas wind, proportional to BH accretion rate |
| `reionization.py` | UV background suppression of gas accretion |
| `utils.py` | Merger-tree I/O, cosmic time/redshift interpolation |

## Interactive exploration

`notebooks/pymctrees_ashvini_demo.ipynb` demonstrates the live pymctrees hook end to end:

1. Generates trees for a few z=5 mass bins and runs Ashvini under delayed vs. instantaneous supernova feedback (default parameters otherwise) -- the bursty, oscillatory gas-mass signature delayed feedback is meant to produce (and its absence under instantaneous feedback) is directly visible, most strongly at the low-mass end.
2. Verifies BH growth and feedback directly: seeding events, the Eddington-limited growth cap (`black_holes.eddington_multiplier`), and a sanity check that growth never exceeds it.
3. Samples a much wider, denser mass range (24 log-spaced points, 1e7-1e11 Msun) and compares stellar/gas mass growth with vs. without BH growth and feedback (seeding disabled entirely) -- the same competitive-gas-budget mechanism `tests/test_bh_growth_feedback.py` checks as a unit test.
4. Shows how to systematically swap the underlying dark matter model (CDM/WDM/FDM/FDM+sharp-k, by pointing `PYMCTREES_CONFIG` at a different pymctrees config file), with a worked CDM-vs-FDM comparison.

Open it with `jupyter notebook notebooks/pymctrees_ashvini_demo.ipynb` after `pip install -e /path/to/pymctrees[camb] jupyter`.

### Numerical scheme

`run1()`/`run_forest()` integrate the five baryonic ODEs (gas mass, gas metals, stellar mass, stellar metals, dust mass) and the BH-mass ODE using closed-form/quadrature updates rather than `scipy.integrate.solve_ivp`, and vectorise the integration across every halo in a forest simultaneously (they all share the same redshift grid). `run1_scalar()` retains the original per-halo, per-timestep `solve_ivp` integration as a trusted reference for validating the fast path -- see the module-level comments in `main.py` for the numerical details, and `tests/test_run1.py` for the cross-validation.

## Testing

```
pip install -r requirements-dev.txt
pytest tests/
```

The suite runs against a small downsampled merger-tree fixture (`tests/fixtures/merger_trees_fixture.h5`, regenerated with `scripts/downsample_trees.py`) and checks output shapes/validity, a pinned reference baseline, agreement between the vectorised and reference integrators, and BH seeding/growth invariants. `test_run_params.py` and `test_run_tree_source.py` cover the `tree_source: pymctrees` config parsing (including a regression test for `run_params.yaml`'s bare-exponent-number YAML gotcha, e.g. `1e10` parsing as a string) and a full live-generation `run()` smoke test; both are skipped if pymctrees isn't installed. `test_bh_growth_feedback.py` covers BH accretion as a genuine gas-mass sink (competing with star formation for the same budget), the Eddington-growth multiplier, and the (optionally delayed) AGN wind term.

## Citation

For more information about the model, please refer to (and kindly cite!) the following publications:

1. Menon & Power 2024, **On bursty star formation during cosmological reionisation – how does it influence the baryon mass content of dark matter halos?**, _Publications of the Astronomical Society of Australia_, 41, id.e049, 11 pp. [DOI: 10.1017/pasa.2024.39](https://ui.adsabs.harvard.edu/abs/2024PASA...41...49M/abstract)
1. Menon, Balu & Power 2025, **On bursty star formation during cosmological reionization -- influence on the metal and dust content of low-mass galaxies**, submitted to _Publications of the Astronomical Society of Australia_, [arXiv link](https://arxiv.org/abs/2508.08363)
