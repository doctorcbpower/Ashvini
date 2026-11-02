import numpy as np


def effective_sn_type(cfg):
    """
    Supernova feedback type actually applied. A "delayed" wind with a zero delay is instantaneous
    feedback: the delayed lookback would otherwise point at the current step, whose star formation
    rate is not yet known, and silently switch feedback off.
    """
    sn_type, t_d = cfg.sn.type, cfg.sn.delay_time
    return "instantaneous" if (sn_type == "delayed" and t_d <= 0.0) else sn_type


def step_feedback_type(eff_type, has_delay_history):
    """Feedback type for one step: "instantaneous" and "no" apply at every step; "delayed" only once a delay has elapsed."""
    if eff_type == "delayed" and not has_delay_history:
        return "no"
    return eff_type


def delay_lookback_index(cosmic_time, t_d):
    """
    For each step j, the largest index k with cosmic_time[k] <= cosmic_time[j] - t_d
    (i.e. the most recent tabulated state at or before t_d in the past),
    or -1 if less than t_d has elapsed since the start of the run.

    Regression fix for a real bug: this used to be a single index offset
    (delay_counter) computed once, from however many *array steps* elapsed
    before crossing t_d early in the run, then reused as a fixed offset for
    the rest of the integration. That's only correct on a uniform
    cosmic-time grid -- trees are stepped uniformly in *redshift*, and
    dt/dz varies enormously with z (e.g. ~40x between z=25 and z=5 for a
    typical cosmology), so a fixed step-offset calibrated near z_max
    represented an effective delay tens of times longer than t_d by the
    time the integration reached low z (verified directly: 182 steps
    calibrated to 15 Myr near z=25 corresponded to ~530 Myr by z=5 for a
    z_max=25 pymctrees-adapted tree). That silently washed out the
    delayed-feedback oscillations the model is specifically built to
    produce. This recomputes the lookback per-step from actual elapsed
    time instead of a step count.
    """
    n = len(cosmic_time)
    target_t = cosmic_time - t_d
    idx = np.searchsorted(cosmic_time, target_t, side="right") - 1
    return np.clip(idx, -1, n - 1)


class SupernovaDelay:
    """
    Which star formation drives the supernova wind in each step (the same star formation rate drives the
    gas, metal and dust losses): zero ("no"), the present rate B_sf * M_gas ("instantaneous"), or the rate
    one delay time ago ("delayed"; no wind until a delay has elapsed). Here every star explodes exactly
    t_d after it forms, a single delay (supernova.kernel = "delta").
    """

    def __init__(self, cfg, cosmic_time):
        self.eff_type = effective_sn_type(cfg)
        self.delay_idx = delay_lookback_index(cosmic_time, cfg.sn.delay_time)

    def feedback_type(self, j):
        return step_feedback_type(self.eff_type, self.delay_idx[j] >= 0)

    def wind_sfr(self, j, feedback_type, B_sf, gas_prev, state):
        """Star formation rate (Msun/Gyr) that drives the wind in step j; state.sfr is the (N, n) history."""
        if feedback_type == "no":
            return np.zeros(len(gas_prev))
        if feedback_type == "instantaneous":
            return B_sf * gas_prev
        return state.sfr[:, self.delay_idx[j]]  # "delayed"


def lifetime_distribution_integral(tau_gyr, lifetimes_gyr):
    """
    G(tau) = integral from 0 to tau of F, where F(s) is the fraction of supernova progenitors whose lifetime
    is at most s (equal-probability progenitor masses, so F is the fraction of the given lifetimes <= s).
    G(tau) = mean(max(tau - t_i, 0)) exactly, and 0 for tau <= 0; evaluated from the sorted lifetimes and their
    cumulative sum, so the cost is that of a search, not of a sum over lifetimes for every tau.
    """
    t = np.sort(np.asarray(lifetimes_gyr, dtype=float))
    csum = np.concatenate([[0.0], np.cumsum(t)])
    tau = np.asarray(tau_gyr, dtype=float)
    k = np.searchsorted(t, tau, side="right")      # lifetimes shorter than or equal to tau
    return (k * tau - csum[k]) / len(t)


def explosion_kernel(cosmic_time, lifetimes_gyr):
    """
    K[j, k]: the fraction of the stars formed (uniformly in time) in step k, between cosmic_time[k-1] and
    cosmic_time[k], whose supernovae go off in step j. With F the fraction exploded by age s and G its
    integral, K[j, k] = [G(T_j - T_{k-1}) - G(T_j - T_k) - G(T_{j-1} - T_{k-1}) + G(T_{j-1} - T_k)] / dT_k.
    Exact for a star formation rate that is constant within each step.

    Explosions that fall in the step of formation itself (progenitors with a lifetime shorter than the step) go
    off at the start of the next step: the star formation of step j is only known at its end. Each row k of K
    then sums to 1 minus the fraction not yet exploded at the end of the run.
    """
    T = np.asarray(cosmic_time, dtype=float)
    n = len(T)
    Tj, Tj1 = T[:, None], np.concatenate([[T[0]], T[:-1]])[:, None]
    Tk, Tk1 = T[None, :], np.concatenate([[T[0]], T[:-1]])[None, :]
    dTk = Tk - Tk1
    G = lambda x: lifetime_distribution_integral(x, lifetimes_gyr)
    with np.errstate(divide="ignore", invalid="ignore"):
        K = (G(Tj - Tk1) - G(Tj - Tk) - G(Tj1 - Tk1) + G(Tj1 - Tk)) / dTk
    K[:, 0] = 0.0                       # node 0 is an initial condition, not a step
    K[0, :] = 0.0
    K = np.tril(K)                      # a star cannot explode before it forms
    K[K < 1.0e-9] = 0.0                 # exactly cancelling terms leave round-off of this size
    diagonal = np.diag(K).copy()
    K[np.arange(1, n), np.arange(0, n - 1)] += diagonal[:-1]   # same-step explosions move to the next step
    np.fill_diagonal(K, 0.0)
    return K


class LifetimeDistributionDelay:
    """
    Supernova delay distribution from the stellar lifetimes (supernova.kernel = "imf"): each core-collapse progenitor
    explodes after its own main-sequence lifetime t_MS(m, Z), the masses following the IMF between m_sn_min
    and m_sn_max (supernova.n_progenitor_masses equal-probability bins). The stars formed in step k, S_k, are
    spread over the later steps with the weights of explosion_kernel, and the star formation rate that drives
    the wind in step j is the sum of S_k K[j, k] / dt_j over earlier steps. It reduces to the delta
    function as the lifetimes approach a single value; the mean delay is the IMF-averaged lifetime.

    The wind is "no" or "instantaneous" exactly as in SupernovaDelay when supernova.type says so.
    """

    def __init__(self, cfg, cosmic_time):
        from ..imf import quantile_progenitor_masses, progenitor_slope
        from ..stellar_lifetimes import t_ms_myr
        self.eff_type = cfg.sn.type
        masses = quantile_progenitor_masses(
            cfg.sn.n_progenitor_masses, cfg.sn.m_sn_min, cfg.sn.m_sn_max, progenitor_slope(cfg.sn.imf)
        )
        self.lifetimes_gyr = t_ms_myr(masses, cfg.sn.lifetime_Z) * 1.0e-3
        self.cosmic_time = np.asarray(cosmic_time, dtype=float)
        self.kernel = explosion_kernel(self.cosmic_time, self.lifetimes_gyr)

    def feedback_type(self, j):
        return self.eff_type

    def wind_sfr(self, j, feedback_type, B_sf, gas_prev, state):
        if feedback_type == "no":
            return np.zeros(len(gas_prev))
        if feedback_type == "instantaneous":
            return B_sf * gas_prev
        dt = self.cosmic_time[j] - self.cosmic_time[j - 1]
        formed = np.diff(state.stars_mass[:, :j], axis=1)          # stars formed in steps 1..j-1
        return formed @ self.kernel[j, 1:j] / dt


def make_delay(cfg, cosmic_time):
    """The supernova delay selected by cfg.sn.kernel."""
    return (LifetimeDistributionDelay if cfg.sn.kernel == "imf" else SupernovaDelay)(cfg, cosmic_time)
