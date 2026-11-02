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
    one delay time ago ("delayed"; no wind until a delay has elapsed). In the Ashvini model every star
    explodes exactly t_d after it forms, a delta-function delay kernel.
    """

    def __init__(self, cfg, cosmic_time):
        self.eff_type = effective_sn_type(cfg)
        self.delay_idx = delay_lookback_index(cosmic_time, cfg.sn.delay_time)

    def feedback_type(self, j):
        return step_feedback_type(self.eff_type, self.delay_idx[j] >= 0)

    def wind_sfr(self, j, feedback_type, B_sf, gas_prev, sfr):
        """Star formation rate (Msun/Gyr) that drives the wind in step j; sfr is the (N, n) history."""
        if feedback_type == "no":
            return np.zeros(len(gas_prev))
        if feedback_type == "instantaneous":
            return B_sf * gas_prev
        return sfr[:, self.delay_idx[j]]  # "delayed"
