import numpy as np

from .. import utils
from ..config import resolve
from ..state import State
from ..components.blackholes import make_black_hole
from ..components.chemistry import MetalsAndDust
from ..components.closure import make_closure
from ..components.delay import make_delay
from ..components.gas import gas_and_stars_step
from ..components.inflow import CosmicInflow
from ..components.star_formation import FreeFallStarFormation


class GridStepper:
    """
    Integrates N haloes at once over the shared redshift grid of a merger tree, one step at a time. In
    each step the components act in this order: black hole (seed, grow, wind), supernova delay (which
    star formation drives the wind), gas and stars, metals, dust. Within a step the rates are frozen at the
    start of the step and each reservoir is advanced by the closed-form solution of its linear ODE.

    Every component is an attribute and can be replaced before `run` with any object that has the same
    methods, e.g.

        stepper = GridStepper(cfg)
        stepper.star_formation = MyStarFormation(cfg)

    Components that need the time grid (the supernova delay and the black hole) are built in `run`, once
    the redshift grid is known; pass `black_hole_factory` or `delay_factory` to replace those.
    """

    def __init__(self, config=None, growth_model=None, **slimdisk_overrides):
        self.config = config  # the explicit Config, or None for config.DEFAULT
        self.cfg = resolve(config)
        self.growth_model = growth_model if growth_model is not None else self.cfg.bh.growth_model
        self.slimdisk_overrides = slimdisk_overrides
        self.star_formation = FreeFallStarFormation(self.cfg)
        self.inflow = CosmicInflow(config)
        self.closure = make_closure(self.cfg)
        self.chemistry = MetalsAndDust(self.cfg)
        self.delay_factory = lambda cosmic_time: make_delay(self.cfg, cosmic_time)
        self.black_hole_factory = lambda cosmic_time, n_halos: make_black_hole(
            self.growth_model, self.cfg, self.config, cosmic_time, n_halos, **self.slimdisk_overrides
        )

    def run(self, halo_mass, halo_mass_rate, redshift, gas_mass0=0.0, stars_mass0=0.0):
        """Returns the State (every array of shape (N, n)) and the cosmic time in Gyr."""
        cosmic_time = utils.time_at_z(redshift)  # Gyr
        n = len(cosmic_time)
        N = halo_mass.shape[0]

        state = State.empty(N, n)
        gas_mass, gas_metals = state.gas_mass, state.gas_metals
        stars_mass, stars_metals = state.stars_mass, state.stars_metals
        sfr, stellar_metallicity = state.sfr, state.stellar_metallicity
        dust_mass, bh_mass = state.dust_mass, state.bh_mass
        gas_mass[:, 0] = gas_mass0
        stars_mass[:, 0] = stars_mass0

        redshift_bcast = np.broadcast_to(redshift, (N, n))
        gas_accretion_rate = self.inflow.rate(redshift_bcast, halo_mass, halo_mass_rate)
        delay = self.delay_factory(cosmic_time)
        black_hole = self.black_hole_factory(cosmic_time, N)

        for j in range(1, n):
            dt = cosmic_time[j] - cosmic_time[j - 1]
            z_mid = utils.z_at_time(0.5 * (cosmic_time[j - 1] + cosmic_time[j]))

            gm_prev = gas_mass[:, j - 1]
            hm_prev = halo_mass[:, j - 1]
            has_gas = gm_prev > 0
            gm_prev_safe = np.where(has_gas, gm_prev, 1.0)

            feedback_type = delay.feedback_type(j)

            A_acc = gas_accretion_rate[:, j - 1]
            B_sf = self.star_formation.rate_coefficient(z_mid)  # scalar, shared across haloes
            ML = self.closure.mass_loading(z_mid, hm_prev, stellar_metallicity[:, j - 1])  # shape (N,)

            # --- black hole first: the gas equation needs this step's accretion and wind ---
            gas_metallicity_prev = np.where(has_gas, gas_metals[:, j - 1] / gm_prev_safe, 0.0)
            bh_mass[:, j], bh_accretion, agn_wind = black_hole.step(
                j, dt, redshift[j - 1], z_mid, hm_prev, gm_prev, stars_mass[:, j - 1],
                gas_metallicity_prev, bh_mass[:, j - 1],
            )

            wind_sfr = delay.wind_sfr(j, feedback_type, B_sf, gm_prev, state)

            gas_raw, stars_formed, _, _ = gas_and_stars_step(
                gm_prev, A_acc, B_sf, ML, wind_sfr, feedback_type, bh_accretion, agn_wind, dt
            )
            gas_mass[:, j] = gas_raw
            stars_mass[:, j] = stars_mass[:, j - 1] + stars_formed

            gas_metals[:, j], stars_metals[:, j] = self.chemistry.metals_step(
                gas_metals[:, j - 1], stars_metals[:, j - 1], A_acc, B_sf, ML, wind_sfr,
                gm_prev, has_gas, gm_prev_safe, dt,
            )
            dust_mass[:, j] = self.chemistry.dust_step(
                dust_mass[:, j - 1], ML, wind_sfr, gm_prev, has_gas, gm_prev_safe, dt
            )

            # --- enforce non-negativity. Negative gas mass is set to zero first: the star formation rate below must be
            # computed from the value after negative masses are set to zero, not the closed-form update's raw output, which legitimately
            # overshoots slightly negative when the decay is fast relative to dt. ---
            gas_mass[:, j] = np.maximum(gas_mass[:, j], 0.0)
            sfr[:, j] = self.star_formation.rate_coefficient(redshift[j]) * gas_mass[:, j]

            gas_metals[:, j] = np.maximum(gas_metals[:, j], 0.0)
            stars_mass[:, j] = np.maximum(stars_mass[:, j], 0.0)
            stars_metals[:, j] = np.maximum(stars_metals[:, j], 0.0)
            dust_mass[:, j] = np.maximum(dust_mass[:, j], 0.0)

            no_gas = gas_mass[:, j] <= 0
            gas_metals[:, j] = np.where(no_gas, 0.0, gas_metals[:, j])
            dust_mass[:, j] = np.where(no_gas, 0.0, dust_mass[:, j])

            has_stars_mass = stars_mass[:, j] > 0
            stellar_metallicity[:, j] = np.where(
                has_stars_mass,
                stars_metals[:, j] / np.where(has_stars_mass, stars_mass[:, j], 1.0),
                0.0,
            )
            stars_metals[:, j] = np.where(has_stars_mass, stars_metals[:, j], 0.0)

        return state, cosmic_time

