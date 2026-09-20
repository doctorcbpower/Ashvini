"""
FROZEN PRE-MVM REFERENCE (2026-09-20). Do not modify. Superseded by ashvini.reservoir_stock (the
Minimal Viable Model); kept only so that the MVM can be compared with the exploratory model on
identical trees. It contains features that are NOT part of the MVM (moving R_nuc, nuclear dark
matter, three gas layers, two star-formation clocks, energy-driven AGN wind, several wind geometries,
the graded Eddington cap, initial conditions).

Reservoir with tracked gas stocks at two scales: a nuclear region and the galaxy.

Model specification (agreed with the user, 2026-09-20; scientific assumptions are
listed here separately from implementation choices).

Scientific assumptions
----------------------
1. Fresh gas enters the halo at Mdot_in = zeta_UV zeta_ch f_b Mdot_halo. Each parcel
   keeps its specific angular momentum j, drawn from a log-normal with median
   j_med = lambda sqrt(2) V_vir R_vir (evaluated when the gas arrives) and width
   sigma_lnj.
2. The same angular-momentum criterion applies at two scales R: a parcel can settle
   inside R if j <= j_thresh(R) = sqrt(G M_enc(R) R).
     * galaxy scale, R_gal = n_rd R_d with R_d = lambda R_vir / sqrt(2) (the disc scale
       length that goes with j_med), M_enc(R_gal) = M_DM,NFW(<R_gal) + M_BH + all gas
       in the galaxy + all stars;
     * nuclear scale, R_nuc (a free parameter, fiducial 100 pc), M_enc(R_nuc) =
       M_BH + M_gas,nuc + M_star,nuc, plus the NFW dark matter inside R_nuc only if
       ``dm_nuclear`` is set.
   The nuclear cut lies inside the galaxy cut. Gas with j > j_thresh(R_gal) stays in the
   halo (no star formation, not subject to galactic winds); gas that fails the nuclear
   test stays in the galaxy. Rejected gas becomes eligible later if j_thresh rises.
   Transfer is instantaneous once eligible.
3. Nuclear gas is a stock: it is fed by the transfer and drained by nuclear star
   formation (eps_sf M_gas,nuc / t_ff), black hole growth (eta_acc M_gas,nuc / t_ff,
   subject to the Eddington cap) and, if the AGN wind is drawn from the nucleus, by that
   wind. t_ff is the free-fall time on M_enc at R_nuc.
4. Galaxy gas (inside R_gal, outside R_nuc) forms stars on eps_ff / t_ff, with
   t_ff = t_ff(R_gal) computed on M_enc(R_gal) ('freefall') or 0.141 t_Hubble(z)
   ('hubble', Menon & Power 2024 Eq. 2). Stars formed there are outside R_nuc and do not
   enter M_enc(R_nuc); they do enter M_enc(R_gal).
5. Stellar winds: 'none'; 'extranuclear' (remove galaxy gas outside R_nuc only);
   'galaxy' (remove galaxy gas everywhere inside R_gal, nuclear included, in proportion
   to gas mass, with no extra resistance assigned to nuclear gas). Mass loading is
   Ashvini's supernovae_feedback.mass_loading_factor in the metal-poor limit, driven by
   the total ('total') or the nuclear ('nuclear') SFR.
6. AGN wind, either 'momentum' (Mdot_out v_out = f_mom L/c, L = eps/(1-eps) Mdot_BH c^2,
   v_out = sigma; King 2003 radiative limit) or 'energy' (King 2003/2005 adiabatic
   limit, Mdot_out = 2 eps_f c^2/sigma^2 Mdot_BH). It is drawn from the galaxy gas
   (extranuclear plus nuclear, in proportion to gas mass) or from the nuclear stock only.
7. The M_BH / M_star target uses the total stellar mass (nuclear + galaxy).
8. The seeding epoch z_seed (start of the integration) is unchanged.

Implementation choices (not scientific assumptions)
---------------------------------------------------
* Each step's inflow is a cohort with two running cuts in j (one per scale). A cohort's
  halo part is the log-normal above the galaxy cut, its galaxy part lies between the
  nuclear and galaxy cuts, so all transfers are exact (no j-binning error in the tail).
  Star formation and winds remove galaxy gas from all cohorts in proportion, which
  leaves each cohort's conditional j distribution unchanged.
* Explicit steps: transfers first (using M_enc at the start of the step), then nuclear
  processes on the post-transfer stock, then galaxy star formation, then the AGN wind,
  then the stellar wind. If the nuclear sinks would exceed the nuclear stock in a step
  they are scaled down together (``limiter_engaged``); baryon mass is conserved exactly.
* The NFW concentration is a parameter (not calibrated at z > 5 or for M < 1e8 Msun).
* Gas present at z_seed (gas_mass0) is deposited as fresh gas at the first step at which
  the halo exists, since j_med is undefined before that.
* ``n_rd = inf`` removes the galaxy cut (all fresh gas enters the galaxy at once).
"""

import numpy as np
from scipy.special import ndtr

from . import black_holes_growth_slimdisk as bh_slimdisk
from . import supernovae_feedback as sn
from .constants import G, Msun_g, pc_cm, c
from .paper_reservoir import (
    F_B_FIDUCIAL, M_HOT_FIDUCIAL, PHI_COLDHOT_FIDUCIAL, LAMBDA_MEDIAN_FIDUCIAL,
    OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL,
    virial_radius_velocity, uv_suppression_mass_only, cold_hot_mode_suppression,
    smooth_halo_growth_rate, critical_seed_paper,
)
from .paper_reservoir_params import PAPER_PARAMS
from .spin import epsilon_from_spin
from .utils import time_at_z, Hubble_time

WIND_MODES = ("none", "extranuclear", "galaxy")
WIND_DRIVERS = ("total", "nuclear")
AGN_WIND_SINKS = ("galaxy", "nuclear")
AGN_WIND_TYPES = ("energy", "momentum")
SF_EXT_TIMESCALES = ("freefall", "hubble")


def momentum_agn_wind_rate(bh_growth_rate, halo_mass, redshift, epsilon, f_mom=1.0):
    """
    Momentum-driven AGN wind (King 2003): the outflow carries the photon momentum,
    Mdot_out v_out = f_mom L / c, with L = [epsilon/(1-epsilon)] Mdot_BH c^2 and
    v_out = sigma, the isothermal velocity dispersion used by the energy-driven wind.
    The mass loading is f_mom [epsilon/(1-epsilon)] c / sigma (a few hundred to a few
    thousand at high z, against 1e4 to 1e5 for the energy-driven wind at epsilon_f=5e-4).
    Returns 0 for an unformed halo.
    """
    sigma = np.asarray(bh_slimdisk.pznk11.velocity_dispersion(halo_mass, redshift), dtype=float)  # cm/s
    valid = sigma > 0
    load = f_mom * epsilon / (1.0 - epsilon) * c / np.where(valid, sigma, 1.0)
    return np.where(valid, load * np.asarray(bh_growth_rate, dtype=float), 0.0)


def median_specific_angular_momentum(halo_mass, z, lam_median=LAMBDA_MEDIAN_FIDUCIAL, H0_km_s_mpc=67.8):
    """j_med = lambda sqrt(2) V_vir R_vir (cm^2/s); 0 where there is no halo."""
    halo_mass = np.asarray(halo_mass, dtype=float)
    has = halo_mass > 0
    R, V = virial_radius_velocity(np.where(has, halo_mass, 1.0), z, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL, H0_km_s_mpc)
    return np.where(has, lam_median * np.sqrt(2.0) * V * R, 0.0)


def nfw_dark_matter_mass(halo_mass, x, c_nfw=PAPER_PARAMS.c_nfw, f_b=F_B_FIDUCIAL):
    """Dark matter mass inside r = x R_vir for an NFW halo of concentration c_nfw."""
    m = lambda y: np.log1p(y) - y / (1.0 + y)
    return (1.0 - f_b) * np.asarray(halo_mass, dtype=float) * m(c_nfw * x) / m(c_nfw)


def _ln_jthresh(M_enc_msun, R_cm):
    """ln sqrt(G M_enc R) (cgs); -inf where M_enc <= 0."""
    M = np.maximum(np.asarray(M_enc_msun, dtype=float), 0.0)
    with np.errstate(divide="ignore"):
        out = 0.5 * (np.log(G) + np.log(np.where(M > 0, M, 1.0) * Msun_g) + np.log(R_cm))
    return np.where(M > 0, out, -np.inf)


def two_scale_transfer(h, g, ln_jmed, S_gal, S_nuc, cut_gal, cut_nuc, ln_jt_nuc, ln_jt_gal, sigma_lnj, n_active):
    """
    Advance the two running cuts of the first ``n_active`` cohorts and move the gas
    that the new cuts release: halo -> galaxy (j below the galaxy cut) and galaxy ->
    nucleus (j below the nuclear cut).

    Cohort k has j ~ LogNormal(ln_jmed[k], sigma_lnj). Its halo part ``h[k]`` has
    j > exp(cut_gal[k]); its galaxy part ``g[k]`` has exp(cut_nuc[k]) < j < exp(cut_gal[k]);
    the rest has moved to the nucleus. ``S_gal``, ``S_nuc`` are P(j > cut) for the
    unconditioned distribution. Cuts never decrease and cut_gal >= cut_nuc. ``h``, ``g``
    and the cut arrays are updated in place. Returns (mass moved into the galaxy,
    mass moved into the nucleus) per halo.
    """
    N = h.shape[0]
    sl = slice(0, n_active)
    alive = (h[:, sl] + g[:, sl]) > 0
    b_new = np.maximum(cut_nuc[:, sl], ln_jt_nuc[:, None])
    a_new = np.maximum(np.maximum(cut_gal[:, sl], ln_jt_gal[:, None]), b_new)
    T_gal = np.zeros(N)
    T_nuc = np.zeros(N)

    ii, cc = np.nonzero(alive & (a_new > cut_gal[:, sl]))
    if ii.size:
        a = a_new[ii, cc]
        S_new = ndtr(-(a - ln_jmed[ii, cc]) / sigma_lnj)
        S_old = S_gal[ii, cc]
        frac = np.clip(np.where(S_old > 0, (S_old - S_new) / np.where(S_old > 0, S_old, 1.0), 0.0), 0.0, 1.0)
        moved = h[ii, cc] * frac
        h[ii, cc] -= moved
        g[ii, cc] += moved
        S_gal[ii, cc] = S_new
        cut_gal[ii, cc] = a
        np.add.at(T_gal, ii, moved)

    ii, cc = np.nonzero(alive & (b_new > cut_nuc[:, sl]))
    if ii.size:
        b = b_new[ii, cc]
        S_new = ndtr(-(b - ln_jmed[ii, cc]) / sigma_lnj)
        S_old = S_nuc[ii, cc]
        denom = S_old - S_gal[ii, cc]  # probability of the galaxy part of the cohort
        frac = np.clip(np.where(denom > 0, (S_old - S_new) / np.where(denom > 0, denom, 1.0), 1.0), 0.0, 1.0)
        moved = g[ii, cc] * frac
        g[ii, cc] -= moved
        S_nuc[ii, cc] = S_new
        cut_nuc[ii, cc] = b
        np.add.at(T_nuc, ii, moved)
    return T_gal, T_nuc


def cohort_transfer(m, ln_jmed, surv, ln_cut, ln_jt, sigma_lnj, n_active):
    """
    Single-scale special case of ``two_scale_transfer``: move to the nucleus the gas in
    the first ``n_active`` cohorts whose j lies below exp(ln_jt) and above the cohort's
    current cut. ``m``, ``surv`` and ``ln_cut`` are updated in place; returns the mass
    moved per halo.
    """
    N = m.shape[0]
    T = np.zeros(N)
    need = (ln_jt[:, None] > ln_cut[:, :n_active]) & (m[:, :n_active] > 0)
    if not need.any():
        return T
    ii, cc = np.nonzero(need)
    x = (ln_jt[ii] - ln_jmed[ii, cc]) / sigma_lnj
    S_new = ndtr(-x)
    S_old = surv[ii, cc]
    frac = np.clip(np.where(S_old > 0, (S_old - S_new) / np.where(S_old > 0, S_old, 1.0), 0.0), 0.0, 1.0)
    moved = m[ii, cc] * frac
    m[ii, cc] -= moved
    surv[ii, cc] = S_new
    ln_cut[ii, cc] = ln_jt[ii]
    np.add.at(T, ii, moved)
    return T


def run_reservoir_stock(
    halo_mass, redshift, seed_mass, f_b=F_B_FIDUCIAL,
    epsilon_sf=PAPER_PARAMS.epsilon_sf, epsilon_sf_ext=PAPER_PARAMS.epsilon_sf_ext,
    a_star=None, epsilon=None,
    chi_crit=PAPER_PARAMS.chi_crit, epsilon_f=PAPER_PARAMS.epsilon_f, eta_acc=PAPER_PARAMS.eta_acc,
    R_nuc_pc=PAPER_PARAMS.R_nuc_pc, compaction_boost=PAPER_PARAMS.compaction_boost,
    M_hot=M_HOT_FIDUCIAL, phi_coldhot=PHI_COLDHOT_FIDUCIAL,
    sigma_lnj=PAPER_PARAMS.sigma_lnj, lam_median=LAMBDA_MEDIAN_FIDUCIAL,
    gas_mass0=PAPER_PARAMS.gas_mass0, stars_mass0=PAPER_PARAMS.stars_mass0,
    halo_growth_rate=None,
    wind_mode=PAPER_PARAMS.wind_mode, wind_driver=PAPER_PARAMS.wind_driver,
    eta_sn_scale=PAPER_PARAMS.eta_sn_scale, agn_wind_sink=PAPER_PARAMS.agn_wind_sink,
    agn_wind_type=PAPER_PARAMS.agn_wind_type, f_mom=PAPER_PARAMS.f_mom,
    n_rd=PAPER_PARAMS.n_rd, c_nfw=PAPER_PARAMS.c_nfw, dm_nuclear=PAPER_PARAMS.dm_nuclear,
    sf_ext_timescale=PAPER_PARAMS.sf_ext_timescale, R_nuc_rd=PAPER_PARAMS.R_nuc_rd,
):
    """
    Integrate the two-scale stock reservoir (see module docstring) along a given halo
    trajectory. Arguments follow ``paper_reservoir.run_reservoir_paper``; the additional
    ones are:

    epsilon_sf_ext, sf_ext_timescale : galaxy-scale star-formation efficiency and clock.
    wind_mode, wind_driver, eta_sn_scale : stellar winds.
    agn_wind_type, f_mom, epsilon_f, agn_wind_sink : AGN wind (type and coupling; which gas).
    R_nuc_pc : nuclear scale (free parameter), or R_nuc_rd : the nuclear scale as a fraction of
        the disc scale length R_d = lam R_vir/sqrt2 of the current halo (overrides R_nuc_pc;
        the scale then follows the halo, an approximation for gas already in the nucleus).
        n_rd : galaxy scale in disc scale lengths
        (inf removes the galaxy cut). c_nfw : NFW concentration. dm_nuclear : include the
        dark matter inside R_nuc in M_enc(R_nuc).

    Returns a dict with the ``run_reservoir_paper`` keys (``gas_mass`` and ``stars_mass``
    are totals) plus ``gas_nuc``, ``gas_ext`` (galaxy gas outside R_nuc), ``gas_halo``
    (gas with j above the galaxy cut), ``gas_pending``, ``stars_nuc``, ``stars_ext``,
    ``transfer`` (mass moved into the nucleus per step), ``transfer_gal`` (into the
    galaxy), ``limiter_engaged`` and ``baryon_residual`` (should be ~1e-12).
    """
    if a_star is not None and epsilon is not None:
        raise ValueError("run_reservoir_stock: pass at most one of a_star, epsilon.")
    for name, val, allowed in (("wind_mode", wind_mode, WIND_MODES), ("wind_driver", wind_driver, WIND_DRIVERS),
                               ("agn_wind_sink", agn_wind_sink, AGN_WIND_SINKS),
                               ("agn_wind_type", agn_wind_type, AGN_WIND_TYPES),
                               ("sf_ext_timescale", sf_ext_timescale, SF_EXT_TIMESCALES)):
        if val not in allowed:
            raise ValueError(f"{name} must be one of {allowed}; got {val!r}.")
    two_scale = bool(np.isfinite(n_rd))
    if not two_scale and sf_ext_timescale == "freefall":
        raise ValueError("sf_ext_timescale='freefall' needs a finite n_rd; use 'hubble' with n_rd=inf.")

    halo_mass = np.atleast_2d(np.asarray(halo_mass, dtype=float))
    N, n = halo_mass.shape
    redshift = np.asarray(redshift, dtype=float)
    cosmic_time = time_at_z(redshift)

    if epsilon is not None:
        eps_rad = epsilon
    elif a_star is not None:
        eps_rad = epsilon_from_spin(a_star)
    else:
        eps_rad = PAPER_PARAMS.epsilon
    kappa_edd = bh_slimdisk.eddington_rate_std_per_unit_mass(eps_rad)

    if halo_growth_rate is None:
        halo_growth_rate = smooth_halo_growth_rate(halo_mass, redshift)
    else:
        halo_growth_rate = np.atleast_2d(np.asarray(halo_growth_rate, dtype=float))

    seed_mass = np.broadcast_to(np.asarray(seed_mass, dtype=float), (N,))

    def zeros():
        return np.zeros((N, n))

    bh, gas_nuc, gas_ext, gas_halo, gas_pend = zeros(), zeros(), zeros(), zeros(), zeros()
    stars_nuc, stars_ext = zeros(), zeros()
    mdot_in, mdot_out, mdot_out_sn = zeros(), zeros(), zeros()
    A_bh_hist, transfer_hist, transfer_gal_hist = zeros(), zeros(), zeros()
    limiter = np.zeros((N, n), dtype=bool)

    bh[:, 0] = seed_mass
    gas_pend[:, 0] = gas_mass0
    stars_ext[:, 0] = stars_mass0

    # cohorts (one per step): halo part h, galaxy part g, ln j_med, survival and ln cut per scale
    h, g = zeros(), zeros()
    ln_jmed = zeros()
    S_gal, S_nuc = np.ones((N, n)), np.ones((N, n))
    cut_gal, cut_nuc = np.full((N, n), -np.inf), np.full((N, n), -np.inf)

    cum_in = np.zeros(N)
    cum_out = np.zeros(N)
    for j in range(1, n):
        dt = cosmic_time[j] - cosmic_time[j - 1]
        z_mid = 0.5 * (redshift[j - 1] + redshift[j])
        c_ = j - 1
        hm = halo_mass[:, c_]
        has_halo = hm > 0
        gn, sn_prev, bh_prev = gas_nuc[:, c_], stars_nuc[:, c_], bh[:, c_]
        ext_prev, se_prev = gas_ext[:, c_], stars_ext[:, c_]
        pend = gas_pend[:, c_].copy()

        # --- fresh inflow (Eq. inflow, Eq. coldhot) ---
        A_acc = uv_suppression_mass_only(np.full(N, z_mid), hm) * cold_hot_mode_suppression(
            hm, M_hot=M_hot, phi=phi_coldhot
        ) * f_b * halo_growth_rate[:, c_]
        mdot_in[:, j] = A_acc
        fresh = A_acc * dt
        cum_in += fresh
        release = np.where(has_halo, pend, 0.0)
        pend -= release
        gas_pend[:, j] = pend
        jmed = median_specific_angular_momentum(hm, z_mid, lam_median)
        h[:, c_] = fresh + release
        ln_jmed[:, c_] = np.log(np.where(jmed > 0, jmed, 1.0))

        # --- thresholds at the two scales (state at the start of the step) ---
        R_vir_cm, _ = virial_radius_velocity(np.where(has_halo, hm, 1.0), z_mid, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL, 67.8)
        if R_nuc_rd is None:
            R_nuc_cm = R_nuc_pc * pc_cm * np.ones(N)
        else:
            R_nuc_cm = R_nuc_rd * lam_median * R_vir_cm / np.sqrt(2.0)
        M_dm_nuc = nfw_dark_matter_mass(hm, R_nuc_cm / R_vir_cm, c_nfw, f_b) if dm_nuclear else 0.0
        ln_jt_nuc = np.where(has_halo, _ln_jthresh(bh_prev + gn + sn_prev + M_dm_nuc, R_nuc_cm), -np.inf)
        if two_scale:
            R_gal_cm = n_rd * lam_median * R_vir_cm / np.sqrt(2.0)
            M_dm_gal = nfw_dark_matter_mass(hm, R_gal_cm / R_vir_cm, c_nfw, f_b)
            M_enc_gal0 = M_dm_gal + bh_prev + gn + ext_prev + sn_prev + se_prev
            ln_jt_gal = np.where(has_halo, _ln_jthresh(M_enc_gal0, R_gal_cm), -np.inf)
        else:
            ln_jt_gal = np.where(has_halo, np.inf, -np.inf)
        ln_jt_gal = np.maximum(ln_jt_gal, ln_jt_nuc)

        T_gal, T = two_scale_transfer(h, g, ln_jmed, S_gal, S_nuc, cut_gal, cut_nuc, ln_jt_nuc, ln_jt_gal, sigma_lnj, c_ + 1)
        transfer_hist[:, j] = T
        transfer_gal_hist[:, j] = T_gal

        # --- nuclear processes on the post-transfer stock ---
        g1 = gn + T
        t_ff = bh_slimdisk.time_freefall_enclosed(bh_prev + g1 + sn_prev + M_dm_nuc, R_nuc_pc=R_nuc_cm / pc_cm)
        with np.errstate(divide="ignore", invalid="ignore"):
            inv_tff = np.where(np.isfinite(t_ff) & (t_ff > 0), 1.0 / t_ff, 0.0)
        A_bh = eta_acc * g1 * inv_tff * compaction_boost**2
        A_bh_hist[:, j] = A_bh
        bh_full = np.maximum(bh_slimdisk.bh_growth_step_slimdisk(bh_prev, A_bh, kappa_edd, chi_crit, dt), 0.0)
        dbh_full = bh_full - bh_prev
        if agn_wind_type == "energy":
            agn_full = bh_slimdisk.king_agn_wind_rate(dbh_full / dt, halo_mass=hm, redshift=z_mid, epsilon_f=epsilon_f) * dt
        else:
            agn_full = momentum_agn_wind_rate(dbh_full / dt, hm, z_mid, eps_rad, f_mom) * dt
        sf_full = epsilon_sf * g1 * inv_tff * dt
        nuclear_agn = agn_wind_sink == "nuclear"
        demand = dbh_full + (agn_full if nuclear_agn else 0.0) + sf_full
        phi = np.where(demand > g1, g1 / np.where(demand > 0, demand, 1.0), 1.0)
        limiter[:, j] = phi < 1.0
        dbh = dbh_full * phi
        sf_n = sf_full * phi
        agn = agn_full * phi  # the wind is linear in the growth actually achieved
        bh[:, j] = bh_prev + dbh
        g2 = g1 - dbh - (agn if nuclear_agn else 0.0) - sf_n
        stars_nuc[:, j] = sn_prev + sf_n
        if nuclear_agn:
            cum_out += agn

        # --- galaxy star formation ---
        ext = g[:, :c_ + 1].sum(axis=1)
        if sf_ext_timescale == "hubble":
            t_sf = 0.141 * Hubble_time(z_mid)
        else:
            M_enc_gal1 = M_dm_gal + bh_prev + g2 + ext + sn_prev + se_prev
            t_sf = bh_slimdisk.time_freefall_enclosed(M_enc_gal1, R_nuc_pc=R_gal_cm / pc_cm)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_sf = np.where(np.isfinite(t_sf) & (t_sf > 0), 1.0 - np.exp(-epsilon_sf_ext * dt / t_sf), 0.0)
        f_sf = np.broadcast_to(f_sf, (N,)).astype(float)
        sf_e = f_sf * ext
        g[:, :c_ + 1] *= (1.0 - f_sf)[:, None]
        ext = ext - sf_e
        stars_ext[:, j] = se_prev + sf_e

        # --- AGN wind drawn from the whole galaxy, in proportion to gas mass ---
        if not nuclear_agn:
            pool = ext + g2
            take = np.minimum(agn, pool)
            keep = np.where(pool > 0, 1.0 - take / np.where(pool > 0, pool, 1.0), 1.0)
            g[:, :c_ + 1] *= keep[:, None]
            ext = ext * keep
            g2 = g2 * keep
            cum_out += take
            agn = take

        # --- stellar winds ---
        sn_lost = np.zeros(N)
        if wind_mode != "none":
            drv = (sf_n + sf_e) if wind_driver == "total" else sf_n
            W = eta_sn_scale * sn.mass_loading_factor(z_mid, hm, 0.0) * drv
            pool = ext if wind_mode == "extranuclear" else ext + g2
            take = np.minimum(W, pool)
            keep = np.where(pool > 0, 1.0 - take / np.where(pool > 0, pool, 1.0), 1.0)
            g[:, :c_ + 1] *= keep[:, None]
            ext = ext * keep
            if wind_mode == "galaxy":
                g2 = g2 * keep
            sn_lost = take
        cum_out += sn_lost
        mdot_out[:, j] = (agn + sn_lost) / dt
        mdot_out_sn[:, j] = sn_lost / dt

        gas_nuc[:, j] = np.maximum(g2, 0.0)
        gas_ext[:, j] = ext
        gas_halo[:, j] = h[:, :c_ + 1].sum(axis=1)

    gas_total = gas_nuc + gas_ext + gas_halo + gas_pend
    stars_total = stars_nuc + stars_ext
    baryon_in = gas_mass0 + stars_mass0 + cum_in
    baryon_now = gas_total[:, -1] + stars_total[:, -1] + (bh[:, -1] - bh[:, 0]) + cum_out
    residual = (baryon_now - baryon_in) / np.maximum(baryon_in, 1.0)

    return dict(
        bh_mass=bh, gas_mass=gas_total, stars_mass=stars_total,
        halo_mass=halo_mass, cosmic_time=cosmic_time, redshift=redshift,
        mdot_in=mdot_in, mdot_out=mdot_out, mdot_out_sn=mdot_out_sn, A_bh=A_bh_hist,
        kappa_edd=float(np.asarray(kappa_edd).item()) if np.asarray(kappa_edd).size == 1 else kappa_edd,
        gas_nuc=gas_nuc, gas_ext=gas_ext, gas_halo=gas_halo, gas_pending=gas_pend,
        stars_nuc=stars_nuc, stars_ext=stars_ext,
        transfer=transfer_hist, transfer_gal=transfer_gal_hist, limiter_engaged=limiter, baryon_residual=residual,
    )


def critical_seed_stock(halo_mass, redshift, **kwargs):
    """M_seed,crit for the stock reservoir: same bisection as ``critical_seed_paper``."""
    return critical_seed_paper(halo_mass, redshift, runner=run_reservoir_stock, **kwargs)
