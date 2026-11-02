"""
The event stepper (milestone M4 of docs/SUPERSET_DESIGN.md): follows individual supernovae. Between events a halo's gas
reservoir is advanced exactly (a linear ODE solved in closed form); each supernova is processed at its own
explosion time, with the halo's escape speed at that time, and removes at most the gas then present. It is the
event-driven counterpart of GridStepper for one halo at a time, and shares the star formation law, the
Ashvini mass-loading closure, and the metal and dust yields with it (Config sections sf, sn, metals, dust).

What it does:
  supply      f_b * max(M_h(t_b) - M_h(t_a), 0) over each substep, reduced by the UV suppression
              (equilibrium form of Okamoto et al 2008, as in Stochastic_Occupation; the grid stepper uses the
              rate-dependent form of ashvini.reionization)
  star formation   dM_g/dt = S - B_sf M_g exactly over each substep (B_sf from the star formation law, evaluated
              at the start of the substep); the stars formed are M_a + S dt - M_b exactly, and are handed to the
              sampler, which returns the explosions they will produce
  explosion   closure.per_supernova(z, M_h, Z_star) * weight, limited at the gas present
  chemistry   metals follow the same exact substep; each explosion adds its yield of metals and dust, destroys
              dust (weight * m_swept * M_d/M_g * f_gas, at most the dust present), and carries off metals and
              dust in proportion to the gas it removes
Not in this stepper: black holes and AGN winds (bh_mass is zero). Optional, all off by default (Config section `events`):
the star formation gate.

Times are in Gyr, masses in Msun.
"""
import heapq
import itertools

import numpy as np

from .. import halo as halo_scales
from .. import utils
from ..config import resolve
from ..constants import Gyr_s
from ..reionization_physics import mu_c, s as uv_shape
from ..components.closure import make_closure
from ..components.star_formation import FreeFallStarFormation
from ..imf import nu_SN
from ..samplers import make_sampler
from ..state import State
from ..trajectory import HaloHistory

_GL_X, _GL_W = np.polynomial.legendre.leggauss(16)
MYR_GYR = 1.0e-3


def exact_substep(M_a, supply, k, dt):
    """
    Exact solution of dM/dt = S - k M over dt for a constant supply rate S = supply / dt. Returns (M_b, formed,
    offset): the final mass, the mass consumed by the loss term k M (= M_a + supply - M_b exactly), and the
    consumption-weighted mean time of that loss measured from the start (16-point Gauss-Legendre quadrature of
    the exact rate).
    """
    x = k * dt
    A = supply / dt / k
    B = M_a - A
    em1 = np.expm1(-x)
    M_b = A + B * (1.0 + em1)
    formed = A * x - B * em1
    nodes = 0.5 * dt * (_GL_X + 1.0)
    rate = k * (A + B * np.exp(-k * nodes))
    num, den = float(np.sum(_GL_W * nodes * rate)), float(np.sum(_GL_W * rate))
    return M_b, formed, (num / den if den > 0.0 else 0.5 * dt)


class _Reservoirs:
    """Mutable masses of one halo's baryons and the queue of pending explosions."""

    def __init__(self, t, gas, stars):
        self.t = t
        self.gas, self.stars, self.ejected, self.returned = gas, stars, 0.0, 0.0
        self.gas_metals = self.stars_metals = self.dust = 0.0
        self.metals_ejected = self.dust_ejected = self.dust_destroyed = 0.0
        self.queue, self._count = [], itertools.count()
        self.n_explosions = 0

    def push(self, t, mass, weight):
        heapq.heappush(self.queue, (t, next(self._count), mass, weight))

    def next_time(self):
        return self.queue[0][0] if self.queue else None


class EventStepper:
    def __init__(self, config=None, sampler="smooth", seed=None, sampler_kwargs=None,
                 sf_dt_max_myr=1.0, step_cap_factor=0.1):
        self.config = config
        self.cfg = cfg = resolve(config)
        self.sampler_kind, self.seed = sampler, seed
        self.sampler_kwargs = sampler_kwargs or {}
        self.sf_dt_max = sf_dt_max_myr * MYR_GYR
        self.step_cap_factor = step_cap_factor
        self.star_formation = FreeFallStarFormation(cfg)
        self.closure = make_closure(cfg)
        self.nu_SN = nu_SN(cfg.sn.imf, cfg.sn.m_sn_min, cfg.sn.m_sn_max)
        self.f_b = utils.Omega_b / utils.Omega_m
        self.Z_IGM = cfg.metals.Z_IGM_Zsun * _z_sun()
        self.yield_Z_per_sn = cfg.metals.Z_yield / self.nu_SN       # Msun of metals per whole supernova
        self.yield_dust_per_sn = cfg.dust.dust_yield / self.nu_SN
        self.gate = self.star_formation_gate                        # replace to try another gate: f(history, t) -> bool

    # --- pieces -------------------------------------------------------------------------------------

    def star_formation_gate(self, history, t):
        """
        Whether gas can form stars at time t: the halo has reached the mass (events.gate_mass) or virial temperature
        (events.gate_T_vir_K) at which it cools, and, below events.gate_switch_z, also events.gate_T_vir_post_K (the
        temperature below which reionisation keeps gas from cooling). No gate configured: always open.
        """
        ev = self.cfg.events
        if ev.gate_mass is None and ev.gate_T_vir_K is None and ev.gate_switch_z is None:
            return True
        M, z = history.mass_at_t(t), history.z_at_t(t)
        if ev.gate_mass is not None:
            open_ = M >= ev.gate_mass
        elif ev.gate_T_vir_K is not None:
            open_ = bool(halo_scales.T_vir_K(M, z) >= ev.gate_T_vir_K)
        else:
            open_ = True
        if ev.gate_switch_z is not None and z < ev.gate_switch_z:
            open_ = open_ and bool(halo_scales.T_vir_K(M, z) >= ev.gate_T_vir_post_K)
        return open_

    def supply(self, history, t_a, t_b):
        """Gas mass entering the reservoir during [t_a, t_b]."""
        M_a, M_b = history.mass_at_t(t_a), history.mass_at_t(t_b)
        f_ret = 1.0
        cfg = self.cfg
        if cfg.reion.UVB_enabled and M_a > 0:
            f_ret = float(uv_shape(mu_c(history.z_at_t(t_a), M_a, cfg.reion.z_reion, cfg.reion.gamma), cfg.reion.omega))
        return self.f_b * f_ret * max(M_b - M_a, 0.0)

    def advance(self, history, res, rng_sampler, t_a, t_b):
        """Advance from t_a to t_b with no explosion in between (see module docstring)."""
        dt = t_b - t_a
        if dt <= 0.0:
            res.t = t_b
            return
        supply = self.supply(history, t_a, t_b)
        if not self.gate(history, t_a) or res.gas + supply <= 0.0:
            res.gas += supply
            res.gas_metals += self.Z_IGM * supply
            res.t = t_b
            return
        k = float(self.star_formation.rate_coefficient(history.z_at_t(t_a)))
        gas_b, formed, offset = exact_substep(res.gas, supply, k, dt)
        res.gas = gas_b
        res.stars += formed
        metals_b, locked, _ = exact_substep(res.gas_metals, self.Z_IGM * supply, k, dt)
        res.gas_metals = metals_b
        res.stars_metals += locked
        if formed > 0.0:
            explosions = rng_sampler.sample(t_a + offset, formed)
            for t_sn, m, w in zip(explosions.t_gyr, explosions.mass_msun, explosions.weight):
                res.push(max(float(t_sn), t_b), float(m), float(w))
        res.t = t_b

    def explode(self, history, res, t_sn, mass, weight):
        """One explosion at its own time (see module docstring)."""
        cfg = self.cfg
        z = history.z_at_t(t_sn)
        M_h = history.mass_at_t(t_sn)
        dZ = weight * self.yield_Z_per_sn
        dD = weight * self.yield_dust_per_sn
        res.gas_metals += dZ
        res.dust += dD
        if res.gas > 0.0 and res.dust > 0.0:
            f_gas = 1.0 - np.exp(-((res.gas / cfg.dust.m_crit) ** cfg.dust.dust_alpha))
            destroyed = min(res.dust, weight * cfg.dust.m_swept * (res.dust / res.gas) * f_gas)
            res.dust -= destroyed
            res.dust_destroyed += destroyed
        Z_star = res.stars_metals / res.stars if res.stars > 0 else 0.0
        wanted = weight * float(self.closure.per_supernova(z, M_h, Z_star))
        removed = min(max(res.gas, 0.0), wanted)
        if removed > 0.0 and res.gas > 0.0:
            fraction = min(removed / res.gas, 1.0)
            lost_Z, lost_d = res.gas_metals * fraction, res.dust * fraction
            res.gas_metals -= lost_Z
            res.dust -= lost_d
            res.metals_ejected += lost_Z
            res.dust_ejected += lost_d
        res.gas -= removed
        res.ejected += removed
        res.n_explosions += 1
        res.t = t_sn

    def evolve(self, history, res, rng_sampler, t_end):
        """Evolve a halo's reservoirs to t_end, processing every explosion at its own time."""
        while res.t < t_end - 1.0e-9 * (t_end - history.t_min + 1.0e-30):
            t_ff = 0.141 * utils.Hubble_time(history.z_at_t(res.t))   # Gyr
            step = min(self.step_cap_factor * t_ff, self.sf_dt_max)
            t_next = min(res.t + step, t_end)
            t_evt = res.next_time()
            if t_evt is not None and t_evt < t_next:
                self.advance(history, res, rng_sampler, res.t, t_evt)
                _, _, mass, weight = heapq.heappop(res.queue)
                self.explode(history, res, t_evt, mass, weight)
            else:
                self.advance(history, res, rng_sampler, res.t, t_next)

    # --- the run ------------------------------------------------------------------------------------

    def run(self, halo_mass, halo_mass_rate, redshift, gas_mass0=0.0, stars_mass0=0.0):
        """
        Same call as GridStepper.run; returns (State, cosmic_time_Gyr). halo_mass_rate is not used: the supply is
        taken from the mass history itself, which is interpolated between the nodes (trajectory.HaloHistory).
        The state is recorded at every node of the redshift grid. Haloes are evolved one after another.
        """
        halo_mass = np.atleast_2d(halo_mass)
        redshift = np.asarray(redshift, dtype=float)
        cosmic_time = utils.time_at_z(redshift)
        N, n = halo_mass.shape
        state = State.empty(N, n)
        self.explosions = np.zeros(N, dtype=int)
        self.ejected = np.zeros((N, n))
        self.returned = np.zeros((N, n))
        rng = np.random.default_rng(self.seed)
        for i in range(N):
            history = HaloHistory(halo_mass[i], redshift)
            sampler = make_sampler(self.cfg, self.sampler_kind, rng, **self.sampler_kwargs)
            res = _Reservoirs(cosmic_time[0], gas_mass0, stars_mass0)
            res.gas_metals = self.Z_IGM * gas_mass0
            for j in range(n):
                if j > 0:
                    self.evolve(history, res, sampler, cosmic_time[j])
                state.gas_mass[i, j] = res.gas
                state.stars_mass[i, j] = res.stars
                state.gas_metals[i, j] = res.gas_metals
                state.stars_metals[i, j] = res.stars_metals
                state.dust_mass[i, j] = res.dust
                state.stellar_metallicity[i, j] = res.stars_metals / res.stars if res.stars > 0 else 0.0
                self.ejected[i, j] = res.ejected
                self.returned[i, j] = res.returned
            state.sfr[i, 1:] = np.diff(state.stars_mass[i]) / np.diff(cosmic_time)
            self.explosions[i] = res.n_explosions
        return state, cosmic_time


def _z_sun():
    from ..constants import Z_SUN
    return Z_SUN
