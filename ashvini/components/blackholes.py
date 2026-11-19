import numpy as np

from .. import agn_feedback as agn
from .. import black_holes_growth as bh_growth
from .. import black_holes_growth_slimdisk as bh_slimdisk
from ..spin import epsilon_from_spin
from .delay import delay_lookback_index


def bh_growth_step(y0, A_bh, kappa_edd, dt):
    """
    Closed-form update of dM_BH/dt = min(A_bh, kappa_edd*y) over a step of
    size dt, with A_bh (gas-supply-limited rate) frozen over the step and
    kappa_edd a true constant (the Eddington-limited growth rate per unit
    mass). Unlike linear_ode_step, this ODE is piecewise (not affine) --
    a halo below the crossover mass A_bh/kappa_edd grows exponentially
    (Eddington-limited) until it reaches it, then grows linearly at rate
    A_bh (gas-supply-limited) for the remainder of the step. Solved exactly
    (not by picking one regime for the whole step) so a step that happens
    to straddle the crossover doesn't overshoot.

    y0 == 0 (not yet seeded) is handled without special-casing: the
    Eddington branch is multiplicative, so 0 stays exactly 0.
    """
    threshold_mass = A_bh / kappa_edd  # A_bh >= 0 always, so this is safe
    eddington_limited = y0 < threshold_mass

    with np.errstate(divide="ignore", invalid="ignore"):
        y0_safe = np.where(y0 > 0, y0, 1.0)
        t_star = np.where(
            y0 > 0,
            np.log(threshold_mass / y0_safe) / kappa_edd,
            np.inf,
        )

    crosses_over = eddington_limited & (t_star < dt)
    stays_eddington = eddington_limited & ~crosses_over
    gas_limited = ~eddington_limited

    # np.select evaluates every choice array for every element regardless of
    # which branch is picked, so t_star's inf entries (uncrossed haloes)
    # would otherwise hit inf * 0 = nan under "dt - t_star" and warn, even
    # though that value is never selected. Clip to a finite stand-in first.
    t_star_safe = np.where(np.isfinite(t_star), t_star, dt)

    return np.select(
        [stays_eddington, crosses_over, gas_limited],
        [
            y0 * np.exp(kappa_edd * dt),
            threshold_mass + A_bh * (dt - t_star_safe),
            y0 + A_bh * dt,
        ],
    )


class _BlackHole:
    """
    Seeding, growth and AGN wind of the central black hole, one step at a time. A subclass sets the growth
    law in `grow`; `step` does the seeding, then the growth, then the (possibly delayed) wind. The AGN wind
    mirrors the supernova delay: cfg.bh.feedback_delay_time (0 = instantaneous).
    """

    def __init__(self, cfg, config, cosmic_time, n_halos):
        self.cfg = cfg
        self.config = config  # the explicit Config, or None
        self.growth_cap = cfg.bh.sigma_feedback.growth_cap_enabled
        self.agn_delay_idx = delay_lookback_index(cosmic_time, cfg.bh.feedback_delay_time)
        self.growth_rate_history = np.zeros((n_halos, len(cosmic_time)))
        self.n_halos = n_halos

    def seed(self, hm_prev, z_prev, gas_metallicity_prev, bh_mass_prev_col):
        newly_seeded, seed_mass = bh_growth.seeding_mask_and_mass(
            hm_prev, z_prev, gas_metallicity_prev, bh_mass_prev_col > 0, cfg=self.cfg
        )
        return np.where(newly_seeded, seed_mass, bh_mass_prev_col)

    def delayed_growth_rate(self, j):
        if self.agn_delay_idx[j] >= 0:
            return self.growth_rate_history[:, self.agn_delay_idx[j]]
        return np.zeros(self.n_halos)

    def step(self, j, dt, z_prev, z_mid, hm_prev, gm_prev, sm_prev, gas_metallicity_prev, bh_mass_prev_col):
        """Returns (bh_mass_new, bh_accretion_rate, agn_wind_rate), the rates in Msun/Gyr."""
        bh_prev = self.seed(hm_prev, z_prev, gas_metallicity_prev, bh_mass_prev_col)
        bh_new, mid = self.grow(dt, z_mid, hm_prev, gm_prev, sm_prev, bh_prev)
        bh_new = np.maximum(bh_new, 0.0)
        # Instantaneous: this step's own accretion, always -- the mass leaves the gas reservoir when
        # it is actually accreted, regardless of whether the *wind* it powers is delayed.
        growth_rate = (bh_new - bh_prev) / dt
        self.growth_rate_history[:, j] = growth_rate
        wind = self.wind(self.delayed_growth_rate(j), mid, z_mid, hm_prev)
        return bh_new, growth_rate, wind


class PZNK11BlackHole(_BlackHole):
    """
    Gas-supply-limited growth with an Eddington cap (Pacucci et al 2015 / PZNK11): dM_BH/dt =
    min(eps_BH/t_ff * M_gas, kappa_edd*M_BH), optionally capped at the M-sigma ceiling
    (cfg.bh.sigma_feedback.growth_cap_enabled), with the wind from agn_feedback.agn_wind_mass_rate.
    """

    def __init__(self, cfg, config, cosmic_time, n_halos):
        super().__init__(cfg, config, cosmic_time, n_halos)
        self.efficiency = cfg.bh.efficiency
        self.kappa_edd = bh_growth.EDDINGTON_RATE_PER_UNIT_MASS * cfg.bh.eddington_multiplier

    def grow(self, dt, z_mid, hm_prev, gm_prev, sm_prev, bh_prev):
        A_bh = (self.efficiency / bh_growth.time_freefall(z_mid)) * gm_prev
        grown = bh_growth_step(bh_prev, A_bh, self.kappa_edd, dt)
        # The growth step is monotonically non-decreasing in y0, so clipping to the ceiling is exactly a
        # "freeze once M_BH reaches the ceiling" ODE, and lets a capped BH resume if the ceiling rises.
        if self.growth_cap:
            grown = np.minimum(grown, bh_growth.growth_ceiling(hm_prev, z_mid, self.cfg))
        # Step-midpoint M_BH, for the wind's optional M-sigma self-regulation coupling only (see
        # gas_evolve.update_gas_reservoir's bh_mass_for_wind docstring).
        mid = bh_growth_step(bh_prev, A_bh, self.kappa_edd, 0.5 * dt)
        if self.growth_cap:
            mid = np.minimum(mid, bh_growth.growth_ceiling(hm_prev, z_mid, self.cfg))
        return grown, mid

    def wind(self, growth_rate, bh_mass_mid, z_mid, hm_prev):
        return agn.agn_wind_mass_rate(
            growth_rate, bh_mass=bh_mass_mid, halo_mass=hm_prev, redshift=z_mid,
            **({} if self.config is None else {"cfg": self.config}),
        )


class SlimDiscBlackHole(_BlackHole):
    """
    2026-paper alternative: enclosed-mass free-fall supply, graded r_crit slim-disc cap and King (2003)
    energy-driven feedback (black_holes_growth_slimdisk.py). No PZNK11 growth-ceiling analogue exists.
    Per-call overrides of a_star/epsilon/r_crit/epsilon_f/eta_acc take priority over cfg.bh.slimdisk.
    """

    def __init__(self, cfg, config, cosmic_time, n_halos, a_star=None, epsilon=None,
                 r_crit=None, epsilon_f=None, eta_acc=None):
        super().__init__(cfg, config, cosmic_time, n_halos)
        slim = cfg.bh.slimdisk
        self.slim = slim
        self.r_crit = r_crit if r_crit is not None else slim.r_crit
        self.epsilon_f = epsilon_f if epsilon_f is not None else slim.epsilon_f
        self.eta_acc = eta_acc if eta_acc is not None else slim.eta_acc
        if epsilon is not None:
            self.kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(epsilon)
        else:
            self.kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(
                epsilon_from_spin(a_star if a_star is not None else slim.a_star)
            )

    def grow(self, dt, z_mid, hm_prev, gm_prev, sm_prev, bh_prev):
        A_bh = bh_slimdisk.nuclear_accretion_rate(
            gm_prev, bh_prev, sm_prev, eta_acc=self.eta_acc,
            R_nuc_pc=self.slim.R_nuc_pc, compaction_boost=self.slim.compaction_boost,
        )
        return bh_slimdisk.bh_growth_step_slimdisk(bh_prev, A_bh, self.kappa_edd, self.r_crit, dt), None

    def wind(self, growth_rate, bh_mass_mid, z_mid, hm_prev):
        return bh_slimdisk.king_agn_wind_rate(
            growth_rate, halo_mass=hm_prev, redshift=z_mid, epsilon_f=self.epsilon_f
        )


def make_black_hole(growth_model, cfg, config, cosmic_time, n_halos, **slimdisk_overrides):
    """The black hole component for cfg.bh.growth_model (or the given override)."""
    if growth_model == "hobbs_slimdisk":
        return SlimDiscBlackHole(cfg, config, cosmic_time, n_halos, **slimdisk_overrides)
    return PZNK11BlackHole(cfg, config, cosmic_time, n_halos)
