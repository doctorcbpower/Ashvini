# Ashvini: a guide for people learning and changing the model

This guide is for someone who has run a galaxy model before but has not seen this one. Read it with `notebooks/ashvini_tour.ipynb` open: the notebook runs each idea below. `MODELS.md` has the equations; this guide says where they live and how to change them.

## 1. What the model does

A dark matter halo grows along a merger tree. Gas falls in at the cosmic baryon fraction of the halo's growth rate (less after reionisation). The gas forms stars at a fixed fraction per free-fall time. Some of the massive stars explode as supernovae after their lifetimes, and the explosions blow gas out of the galaxy. Metals and dust are made by the explosions and carried off in the wind. A black hole can seed and grow, and its wind removes gas too. All of it is one reservoir of gas per halo: the model has no radial structure.

Two ways of advancing this in time share the same physics:

* the **grid integrator** (`GridStepper`) advances every halo of a forest together, step by step on the tree's redshift grid, with supernova feedback as a smooth rate;
* the **event integrator** (`EventStepper`) follows one halo at a time, through individual supernovae at the moment each explodes.

## 2. Five minutes

```python
from ashvini import utils
from ashvini.config import DEFAULT, override
from ashvini.main import run_forest, run_events

halo_mass, rate, z = utils.read_trees("tests/fixtures/merger_trees_fixture.h5", 1e10)
result = run_forest(halo_mass, rate, z)                        # default parameters
other = run_forest(halo_mass, rate, z, config=override(DEFAULT, **{"sf.efficiency": 0.03}))
print(result["stars_mass"][:, -1], other["stars_mass"][:, -1])  # stars at the last redshift, per halo
```

`result` is a dictionary of arrays of shape (number of haloes, number of steps): `gas_mass`, `stars_mass`, `gas_metals`, `stars_metals`, `dust_mass`, `bh_mass`, `sfr`, plus the inputs and the cosmic time. Masses are in Msun, time in Gyr.

Note: the growth-rate column of the test fixture does not integrate to its mass history (the fixture was thinned in redshift). Fine for checking that code runs; for science, take the rate from a consistent tree, or from the masses.

## 3. How the code is laid out

| Layer | Where | What it is |
|---|---|---|
| Parameters | `config.py`, `run_params.py`, `run_params.yaml` | One `Config` tree (sections `sf`, `sn`, `reion`, `metals`, `dust`, `bh`, `events`). `DEFAULT` is the file's values; `override(cfg, **{"sn.type": "no"})` gives a changed copy. |
| State | `state.py` | The arrays that evolve, one per quantity, shape (haloes, steps). |
| Components | `components/` | One piece of physics each: `star_formation`, `inflow`, `delay` (how long after formation stars explode), `closure` (the wind model: how much gas each supernova removes), `gas` (the gas and stellar mass balance), `chemistry` (metals, dust), `blackholes`, `ode` (the exact solution of a linear equation over a step). |
| Time integrators | `steppers/grid.py`, `steppers/events.py` | Apply the components in order and advance in time. |
| Supporting | `imf.py`, `stellar_lifetimes.py`, `halo.py`, `samplers.py`, `trajectory.py`, `presets.py` | The IMF and supernova rate, the lifetime tables, halo scales (virial radius, escape speed, T_vir), the random or expected supernovae of a stellar population, the halo history for the event integrator, ready-made parameter sets. |
| Frozen | `reservoir_stock*.py`, `paper_reservoir*.py` | Models of the 2026 papers. Do not edit them; see `docs/RESERVOIRS.md`. |

Rule of thumb: **physics goes in a component, parameters in the Config, the order of operations in the time integrator**. If you find yourself editing the integrator to change an equation, the equation probably wants to be a component.

## 4. One step of the grid integrator

For each step from `j-1` to `j` (`steppers/grid.py`, `GridStepper.run`), with every rate held at its value at the start of the step and each reservoir advanced by the exact solution of its linear ODE:

1. the black hole is seeded if the halo qualifies, grows, and launches a wind (`components/blackholes.py`);
2. the star formation that drives the supernova wind is chosen: none, the present rate, or the rate some time ago (`components/delay.py`);
3. the gas is advanced, and the stars formed follow from mass balance: for `dy/dt = F - D y` the time integral of `D y` over the step is `y0 + F dt - y1` (`components/gas.py`), so gas, stars and removed mass add up to the supply to round-off;
4. metals and dust are advanced with the same wind rate (`components/chemistry.py`);
5. any negative mass is set to zero, and the star formation rate and stellar metallicity are recorded.

The event integrator (`steppers/events.py`) replaces steps 2 to 5 with a list of pending explosions: stars form in short substeps; `samplers.py` turns the stars formed into supernovae with explosion times; each supernova removes gas at the moment it explodes.

## 5. Recipes

**Change a parameter.** `override(DEFAULT, **{"dotted.path": value})`. The path follows `run_params.yaml`. A wrong name raises an error.

**Choose the physics from the Config.** `sn.type` (delayed, instantaneous, no); `sn.kernel`, the distribution of delays between a star forming and exploding (`delta`: every star explodes `delay_time` after forming; `imf`: each massive star after its own lifetime); `sn.closure`, the wind model, i.e. how much gas each supernova removes (`ashvini`: the calibrated mass loading; `momentum`: a momentum-driven wind; `energy`: an energy-driven wind); `bh.growth_model`. `presets.stochocc_smooth_control()` sets up the Stochastic_Occupation smooth limit in one call.

**Replace a component.**

```python
from ashvini.steppers.grid import GridStepper
from ashvini.components.star_formation import FreeFallStarFormation

class MyLaw(FreeFallStarFormation):
    def rate_coefficient(self, redshift):          # 1/Gyr, the same for all haloes
        return 2.0 * super().rate_coefficient(redshift)

stepper = GridStepper()
stepper.star_formation = MyLaw(DEFAULT)
state, cosmic_time = stepper.run(halo_mass, rate, z)
```

Look at the component you are replacing: its methods are the interface. The integrator calls only those.

**Add a parameter.** (1) Add a field to the right dataclass in `run_params.py` with a default that keeps the old behaviour; (2) read it from `cfg` where it is used, never from a module variable; (3) check it in `config.validate` if some values are wrong; (4) add it, commented, to `run_params.yaml`; (5) write a test that the default changes nothing and that the new value does what you expect; (6) add it to the tables in `MODELS.md`. Names carry their units (`Z_IGM_Zsun`, `bank_tau_myr`).

**Try the event options.** The `events` section adds: a cooling threshold (no star formation until the halo is hot enough to cool), a delay in the gas reaching the galaxy, gas falling back after it was blown out, explosion energy kept for removing gas that arrives later, a limit on the fraction of the gas one supernova can reach, a different energy coupling above a halo mass, and a wind travel time. All are off by default. `MODELS.md` has the table.

## 6. How to know a change is right

The tests are the specification you can run (`pytest tests/`, about five minutes; `pytest --parity` also reruns Stochastic_Occupation, if it is next to this repository).

* **Reference outputs** (`tests/test_characterisation.py`, `tests/golden/`): the outputs of the general chain for 16 parameter sets are stored, and a rearrangement of the code must reproduce them. If you change the model on purpose, run `python -m tests.characterisation`, read what changed, and say in the commit message why the numbers moved.
* **Identical arrays**: the M2 and M3 rearrangements were checked to give identical arrays, not just close ones.
* **Mass balance**: gas + stars + removed mass equals the supply to round-off (`tests/test_event_stepper.py`, `test_event_options.py`).
* **Limiting cases**: no feedback gives gas + stars = supply; zero efficiency gives no stars.
* **Agreement with Stochastic_Occupation**: both integrators reproduce its smooth (non-random) supernova model (`tests/test_stochocc_parity.py`). `docs/SUPERSET_DESIGN.md` section 8 has the numbers.

A run that completes and gives a pretty figure is not evidence of a correct model. Before you write physics, write down what it should do in a case you can solve by hand, and test that.

## 7. Conventions

* Units: Msun, Gyr, cm and s inside the wind models; metallicity as a mass fraction inside the code, as a multiple of Z_sun (0.015) in the Config.
* Rates are per Gyr; `sfr` in the output is Msun per Gyr.
* The grid steps are on the tree's redshift grid; the answer depends on the step size (first order). The event integrator's substeps are at most 1 Myr and 0.1 free-fall times.

## 8. What to be careful about

* The grid integrator's UV suppression of gas accretion depends on the halo growth rate (`reionization.py`); the event integrator's depends on mass and redshift only. They differ when the UV background is on.
* The momentum-driven and energy-driven winds have no dependence on metallicity; Ashvini's own mass loading does (`supernovae_feedback.py`). They agree only at low stellar metallicity.
* If one supernova can reach only part of the gas (`events.access_fraction`), the answer depends on how the expected supernovae are divided into separate explosions (`sampler_kwargs={"K": ...}`).
* The event integrator has no black holes yet.
* The tree's growth-rate column and its mass history must be consistent (section 2).

## 9. Where next

Work through the exercises at the end of the notebook. Then pick one equation in `MODELS.md`, find the component that implements it, and change it.
