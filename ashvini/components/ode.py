import numpy as np


def linear_ode_step(y0, forcing, decay, dt):
    """
    Closed-form update of dy/dt = forcing - decay*y over a step of size dt,
    assuming forcing/decay are constant across the step. Vectorised and safe
    for decay == 0 (no supernova/mass-loading term active): rather than
    dividing by decay directly, phi(x) = (1 - exp(-x)) / x is evaluated via
    a Taylor series near x = decay*dt = 0, which is well-defined there.
    """
    x = decay * dt
    exp_term = np.exp(-x)
    small = np.abs(x) < 1e-6
    phi = np.where(
        small,
        dt * (1.0 - x / 2.0 + x * x / 6.0),
        (1.0 - exp_term) / np.where(x == 0, 1.0, decay),
    )
    return y0 * exp_term + forcing * phi
