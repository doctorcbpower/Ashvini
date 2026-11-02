"""M2: the components behind GridStepper can be used on their own and replaced one at a time."""
import numpy as np
import pytest

from ashvini import main
from ashvini.components.ode import linear_ode_step
from ashvini.components.star_formation import FreeFallStarFormation
from ashvini.config import DEFAULT
from ashvini.steppers.grid import GridStepper


@pytest.fixture
def forest():
    n = 40
    redshift = np.linspace(20, 6, n)
    halo_mass = 1e10 * np.exp(-np.linspace(3, 0, n))[None, :] * np.array([[1.0], [3.0]])
    rate = np.gradient(halo_mass, axis=1) * 5.0
    return halo_mass, rate, redshift


def test_run_forest_is_the_default_stepper(forest):
    halo_mass, rate, redshift = forest
    result = main.run_forest(halo_mass, rate, redshift)
    state, _ = GridStepper().run(halo_mass, rate, redshift)
    for name, values in state.as_dict().items():
        assert np.array_equal(result[name], values)


def test_replacing_the_star_formation_law_changes_the_stars_not_the_inflow(forest):
    halo_mass, rate, redshift = forest

    class DoubleRate(FreeFallStarFormation):
        def rate_coefficient(self, redshift):
            return 2.0 * super().rate_coefficient(redshift)

    base, _ = GridStepper().run(halo_mass, rate, redshift)
    stepper = GridStepper()
    stepper.star_formation = DoubleRate(DEFAULT)
    faster, _ = stepper.run(halo_mass, rate, redshift)
    assert np.all(faster.stars_mass[:, -1] > base.stars_mass[:, -1])
    assert np.all(faster.gas_mass[:, -1] < base.gas_mass[:, -1])


def test_linear_ode_step_matches_the_analytic_solution():
    y0, F, D, dt = 3.0, 2.0, 0.7, 1.3
    exact = F / D + (y0 - F / D) * np.exp(-D * dt)
    assert linear_ode_step(np.array(y0), F, D, dt) == pytest.approx(exact, rel=1e-14)
    assert linear_ode_step(np.array(y0), F, 0.0, dt) == pytest.approx(y0 + F * dt, rel=1e-14)
