"""
Alternative BH nuclear accretion + super-Eddington cap + AGN feedback
model ("hobbs_slimdisk"), selected via black_holes.growth_model:
"hobbs_slimdisk" in run_params.yaml (default remains "pznk11_freefall",
i.e. the pre-existing black_holes_growth.py path -- unchanged).

Ported and adapted from the standalone prototype developed for the 2026
"Differential Growth" paper (see
docs/2026_paper_session_code_catalogue.md). Three physical differences
from the existing black_holes_growth.py/agn_feedback.py default:

1. Nuclear accretion supply is a free-fall estimate on the *total
   enclosed mass* M_enc = M_BH + M_gas + M_star at a fixed nuclear radius
   R_nuc (Hobbs, Power, Nayakshin & King 2012), not a fixed
   0.141 x Hubble-time free-fall time. This matters specifically for a
   light seed, where M_BH does not dominate M_enc and a Bondi-type
   M_BH^2 estimate would be far too small (see nuclear_accretion_rate's
   docstring).
2. The super-Eddington cap is graded rather than a hard multiplier: full
   supply-limited growth is allowed up to r_crit x the standard Eddington
   rate; beyond that the effective rate saturates at Mdot_acc / r_crit
   (Watarai et al. 2000; Madau, Haardt & Dotti 2014; Lupi et al. 2024
   slim-disc picture), rather than a single constant f_Edd applied at
   every mass.
3. AGN feedback is King (2003) energy-driven (Mdot_wind = 2*epsilon_f*
   c^2/sigma^2 * Mdot_BH) rather than a constant eta_agn mass-loading
   factor or the PZNK11 M-sigma switch.

Radiative efficiency epsilon (which sets the standard Eddington rate
normalisation, kappa_edd below) is exposed as a genuine free parameter
here, via ashvini.spin.epsilon_from_spin(a_star) -- unlike
black_holes_growth.eddington_bh_growth, which has no epsilon dependence
at all (see that module; its kappa_per_s bakes in a specific, implicit
convention). See MODELS.md for the full equations and the paper's
Section 5.5 for why epsilon (and, through it, BH spin) turns out to be
the dominant source of uncertainty in this growth model.

Compaction boost (Booth & Schaye 2009 functional form, anchored to
Lapiner, Dekel & Dubois 2021): the paper's own standalone prototype used
a stochastic, time-varying Phi_hat(t) (an Ornstein-Uhlenbeck process) to
represent bursty compaction. That is NOT reproduced here -- it doesn't
fit the closed-form-per-step integration scheme run_forest() relies on
(Phi_hat would need to vary *within* a step for its stochastic character
to mean anything, breaking the "coefficients frozen over a step" pattern
every other RHS term in this package already uses). Here
`compaction_boost` is a single configurable multiplier
(black_holes.slimdisk.compaction_boost in run_params.yaml, default 1.0 =
no boost), constant for a given run -- representing a fixed assumed
compaction state, not a time-varying process. The full stochastic
version remains paper-analysis-only; see
docs/2026_paper_session_code_catalogue.md.
"""

import numpy as np

from .run_params import PARAMS
from .spin import epsilon_from_spin
from . import black_holes_growth as pznk11  # reuse velocity_dispersion() for the King wind term

from .constants import G, c, m_p, sigma_thomson, Gyr_s, Msun_g, pc_cm

_slim = PARAMS.bh.slimdisk

eta_acc = _slim.eta_acc
R_nuc_pc = _slim.R_nuc_pc
compaction_boost = _slim.compaction_boost
r_crit = _slim.r_crit
epsilon_f = _slim.epsilon_f
a_star = _slim.a_star
epsilon_rad = epsilon_from_spin(a_star)


def time_freefall_enclosed(M_enc_msun, R_nuc_pc=R_nuc_pc):
    """
    Free-fall time (Gyr) at fixed nuclear radius R_nuc (pc), on the total
    enclosed mass M_enc = M_BH + M_gas + M_star (Msun):

        t_ff = sqrt(R_nuc^3 / (2 * G * M_enc))

    Replaces naive Bondi-Hoyle (Mdot ~ M_BH^2), which is the wrong
    estimator whenever the BH's own mass does not dominate the local
    enclosed mass -- exactly the light-seed case this growth model exists
    to handle (Hobbs, Power, Nayakshin & King 2012). M_enc <= 0
    (unformed/empty nucleus) returns +inf (zero accretion), guarded the
    same way other undefined-quantity cases are handled elsewhere in this
    package (e.g. black_holes_growth.velocity_dispersion).
    """
    M_enc_msun = np.asarray(M_enc_msun, dtype=float)
    resolved = M_enc_msun > 0
    M_enc_g = np.where(resolved, M_enc_msun, 1.0) * Msun_g
    R_nuc_cm = R_nuc_pc * pc_cm

    t_ff_s = np.sqrt(R_nuc_cm**3 / (2.0 * G * M_enc_g))
    t_ff_gyr = t_ff_s / Gyr_s
    return np.where(resolved, t_ff_gyr, np.inf)


def nuclear_accretion_rate(gas_mass, bh_mass, stars_mass,
                            eta_acc=eta_acc, R_nuc_pc=R_nuc_pc,
                            compaction_boost=compaction_boost):
    """
    Mdot_acc (Msun/Gyr) = eta_acc * (gas_mass / t_ff(M_enc)) * compaction_boost^2,
    the Booth & Schaye (2009)-style density-boosted free-fall supply rate
    (beta=2 fixed, matching the paper; compaction_boost plays the role of
    their Phi_hat, capped externally by the caller if a Phi_max ceiling is
    wanted -- not enforced here since compaction_boost is a fixed run
    parameter, not a state variable that could exceed a ceiling on its
    own).
    """
    M_enc = np.asarray(bh_mass, dtype=float) + np.asarray(gas_mass, dtype=float) + np.asarray(stars_mass, dtype=float)
    t_ff = time_freefall_enclosed(M_enc, R_nuc_pc)
    with np.errstate(divide="ignore", invalid="ignore"):
        Mdot_ff = np.where(np.isfinite(t_ff) & (t_ff > 0), np.asarray(gas_mass, dtype=float) / np.where(t_ff > 0, t_ff, 1.0), 0.0)
    return eta_acc * Mdot_ff * compaction_boost**2


def eddington_rate_std_per_unit_mass(epsilon=epsilon_rad):
    """
    kappa_edd (1/Gyr): standard Eddington-limited growth rate per unit
    M_BH, kappa_edd = (1-epsilon)/epsilon * 4*pi*G*m_p/(sigma_T*c) --
    i.e. 1/t_Sal with t_Sal = [epsilon/(1-epsilon)] * t_Edd,
    t_Edd = sigma_T*c/(4*pi*G*m_p) = 0.450 Gyr -- matching the paper's
    Section 2.1 convention exactly (and NOT the same convention as
    black_holes_growth.eddington_bh_growth, which has no epsilon
    dependence at all -- these are two independently configurable growth
    models, not required to agree on this choice).
    """
    kappa_per_s = 4.0 * np.pi * G * m_p / (sigma_thomson * c)
    kappa_per_gyr = kappa_per_s * Gyr_s
    return ((1.0 - epsilon) / epsilon) * kappa_per_gyr


EDDINGTON_RATE_PER_UNIT_MASS_FIDUCIAL = float(eddington_rate_std_per_unit_mass(epsilon_rad))


def bh_growth_step_slimdisk(y0, A_bh, kappa_edd, r_crit, dt):
    """
    Closed-form update of

        dM/dt = min(A_bh, cap(M)),
        cap(M) = kappa_edd*M            if A_bh/(kappa_edd*M) <= r_crit
                 A_bh / r_crit          otherwise

    over a step of size dt, with A_bh (nuclear supply, e.g. from
    nuclear_accretion_rate) and kappa_edd frozen constant over the step
    (same "coefficients frozen at step start/midpoint" pattern
    main._bh_growth_step already uses for the pznk11_freefall model).

    Three regimes in M, boundaries M1 = A_bh/(kappa_edd*r_crit) and
    M2 = A_bh/kappa_edd = M1*r_crit (r_crit >= 1 assumed: it is a
    threshold *ratio* of supply to standard Eddington beyond which the
    slim-disc cap engages, so M1 <= M2):

        M < M1  : "S" (supply-rich, slim-disc-throttled) -- growth is
                  CONSTANT at A_bh/r_crit (cap saturates the graded
                  Eq.~3 ceiling; a genuinely super-Eddington, capped rate)
        M1<=M<M2: "E" (standard-Eddington-limited) -- growth is
                  EXPONENTIAL at kappa_edd*M
        M >= M2 : "G" (gas-supply-limited) -- growth is CONSTANT at A_bh

    Continuous at both boundaries by construction (checked in
    tests/test_black_holes_growth_slimdisk.py). As r_crit -> infinity,
    M1 = A_bh/(kappa_edd*r_crit) -> 0, so regime S (which only applies
    below M1) vanishes for any y0 > 0 and this reduces exactly to
    main._bh_growth_step's two-regime solver (Eddington-limited vs.
    gas-supply-limited only -- no super-Eddington cap ever engages) --
    verified directly in the test suite as the primary correctness check
    on this closed form, alongside a direct comparison to a numerical
    solve_ivp reference at a realistic finite r_crit.

    Vectorised (y0, A_bh, kappa_edd, r_crit may be scalars or
    broadcastable arrays); A_bh <= 0 (no supply, e.g. gas_mass=0) returns
    y0 unchanged rather than evaluating M1/M2 (which are undefined at
    A_bh=0).
    """
    y0 = np.asarray(y0, dtype=float)
    A_bh = np.asarray(A_bh, dtype=float)
    kappa_edd = np.asarray(kappa_edd, dtype=float)
    r_crit = np.asarray(r_crit, dtype=float)

    no_supply = A_bh <= 0
    A_bh_safe = np.where(no_supply, 1.0, A_bh)
    kappa_safe = np.where(kappa_edd > 0, kappa_edd, 1.0)

    M1 = A_bh_safe / (kappa_safe * r_crit)
    M2 = M1 * r_crit  # = A_bh_safe / kappa_safe

    rate_S = A_bh_safe / r_crit  # constant growth rate while M < M1

    # --- case: starts at/above M2 (regime G only) ---
    in_G0 = y0 >= M2
    final_G0 = y0 + A_bh_safe * dt

    # --- case: starts in [M1, M2) (regime E, possibly crossing into G) ---
    in_E0 = (y0 >= M1) & (y0 < M2)
    y0_E_safe = np.where(y0 > 0, y0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        t2_from_E0 = np.where(y0 > 0, np.log(M2 / y0_E_safe) / kappa_safe, np.inf)
    E0_only = in_E0 & (t2_from_E0 >= dt)
    E0_then_G = in_E0 & (t2_from_E0 < dt)
    with np.errstate(over="ignore"):
        final_E0_only = y0 * np.exp(np.minimum(kappa_safe * dt, 700.0))
    t2_from_E0_safe = np.where(np.isfinite(t2_from_E0), t2_from_E0, dt)
    final_E0_then_G = M2 + A_bh_safe * (dt - t2_from_E0_safe)

    # --- case: starts in regime S (y0 < M1) ---
    in_S0 = y0 < M1
    t1 = (M1 - y0) / rate_S  # time to reach M1 at the constant rate_S
    S0_only = in_S0 & (t1 >= dt)
    final_S0_only = y0 + rate_S * dt

    # crosses into E at t1; time spent in E to reach M2 is ln(r_crit)/kappa
    with np.errstate(divide="ignore", invalid="ignore"):
        t_E_span = np.where(r_crit > 0, np.log(np.maximum(r_crit, 1.0)) / kappa_safe, 0.0)
    t2_from_S0 = t1 + t_E_span
    S0_then_E = in_S0 & (t1 < dt) & (t2_from_S0 >= dt)
    with np.errstate(over="ignore"):
        final_S0_then_E = M1 * np.exp(np.minimum(kappa_safe * (dt - t1), 700.0))

    S0_then_E_then_G = in_S0 & (t1 < dt) & (t2_from_S0 < dt)
    t2_from_S0_safe = np.where(np.isfinite(t2_from_S0), t2_from_S0, dt)
    final_S0_then_E_then_G = M2 + A_bh_safe * (dt - t2_from_S0_safe)

    result = np.select(
        [in_G0, E0_only, E0_then_G, S0_only, S0_then_E, S0_then_E_then_G],
        [final_G0, final_E0_only, final_E0_then_G, final_S0_only, final_S0_then_E, final_S0_then_E_then_G],
        default=y0,
    )
    return np.where(no_supply, y0, result)


def king_agn_wind_rate(bh_growth_rate, halo_mass, redshift, epsilon_f=epsilon_f):
    """
    King (2003) energy-driven AGN wind: Mdot_wind = (2*epsilon_f*c^2 /
    sigma^2) * Mdot_BH, sigma the isothermal-sphere velocity dispersion
    (reusing black_holes_growth.velocity_dispersion, CGS) -- the
    self-regulation mechanism underlying the M_BH-sigma relation
    (Power, Zubovas, Nayakshin & King 2011), used here as a direct,
    constant-epsilon_f coupling rather than the existing package's
    optional M_BH/M_sigma logistic switch (agn_feedback.coupling_switch).
    sigma <= 0 (unformed halo) returns 0.
    """
    sigma = pznk11.velocity_dispersion(halo_mass, redshift)  # cm/s
    bh_growth_rate = np.asarray(bh_growth_rate, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    valid = sigma > 0
    sigma_safe = np.where(valid, sigma, 1.0)
    wind_rate = (2.0 * epsilon_f * c**2 / sigma_safe**2) * bh_growth_rate
    return np.where(valid, wind_rate, 0.0)
