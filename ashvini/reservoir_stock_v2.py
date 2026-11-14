"""
Revised Minimal Viable Model (MVM v2): ``ashvini.reservoir_stock`` with three switchable changes.

The frozen model in ``ashvini/reservoir_stock.py`` (SHA-256 31c701dd...) is NOT modified. This module copies its
integration loop line by line and adds three options, each of which removes an approximation identified in
``model-verification-audit-2026-11-14.md`` (High-Redshift-Black-Hole-Growth paper folder):

1. ``dm_nuclear`` : include the NFW dark matter inside R_nuc in the enclosed mass used for the nuclear capture threshold
   and the nuclear free-fall time, M_enc(R_nuc) = M_DM,NFW(<R_nuc) + M_BH + M_gas,nuc + M_star,nuc. The frozen model
   omits it although it exceeds the nuclear baryonic mass of a light seed by orders of magnitude, and the paper's own
   test showed that this changes the light-seed regime. Same NFW profile (concentration c_nfw, dark fraction 1 - f_b)
   as the galaxy-scale term.
2. ``inflow``     : "rate" is the frozen scheme (GRUMPY spline rate multiplied by the step). "conservative" delivers in each
   step the mass the halo actually gained on the integration grid, so that the cumulative baryon inflow equals
   f_b times the halo mass gained (times the suppression factors) by construction. The rate scheme was measured to
   over-deliver by a median of 17 per cent at 3e10 Msun on merger trees (it is exact on smooth histories).
3. ``cosmology``  : None keeps the frozen mixture (cosmic time and H(z) in sigma_h from the astropy Planck18 tables, halo
   structure from H0 = 67.8, Omega_m = 0.308). A ``PaperCosmology`` instance uses that one flat LCDM cosmology for the
   time axis, the tree interpolation, the virial quantities and sigma_h. The two differ by 4.4 Myr in the duration of
   the z = 25 to 5 run, about 0.09 e-folds of Salpeter growth.

With ``dm_nuclear=False, inflow="rate", cosmology=None`` (``LEGACY``) the results are identical to the frozen model
(tests/test_reservoir_stock_v2.py asserts this to machine precision). Everything else, including the order of operations
within a step and all parameter defaults, is unchanged; the docstring of ``reservoir_stock.py`` describes it.

The public defaults of ``run_reservoir_v2`` and ``critical_seed_v2`` are the REVISED physics (all three changes on).
"""

import numpy as np
from scipy.special import ndtr

from . import black_holes_growth_slimdisk as bh_slimdisk
from . import supernovae_feedback as sn
from .constants import G, Msun_g, pc_cm, c, Gyr_s
from .paper_reservoir import (
    F_B_FIDUCIAL, M_HOT_FIDUCIAL, PHI_COLDHOT_FIDUCIAL, LAMBDA_MEDIAN_FIDUCIAL,
    OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL,
    virial_radius_velocity, uv_suppression_mass_only, cold_hot_mode_suppression,
)
from .paper_reservoir_params import PAPER_PARAMS
from .utils import time_at_z as _time_at_z_planck18
from .reservoir_stock import (
    strict_eddington_step, momentum_agn_wind_rate, median_specific_angular_momentum, nfw_dark_matter_mass,
    ln_jthresh, star_formation_demand, nuclear_transfer, _settle, target_residual,
)

LEGACY = dict(dm_nuclear=False, inflow="rate", cosmology=None)


# ---------------------------------------------------------------------------------------------------------------
# One flat LCDM cosmology for the whole calculation
# ---------------------------------------------------------------------------------------------------------------
class PaperCosmology:
    """
    Flat matter + Lambda cosmology (no radiation) with closed-form cosmic time.

    t(z) = 2 / (3 H0 sqrt(OL)) * arcsinh( sqrt(OL / Om) (1 + z)^(-3/2) ), which inverts to
    1 + z = [ sqrt(Om / OL) sinh(3 H0 sqrt(OL) t / 2) ]^(-2/3).
    Defaults are the cosmology of the merger trees and the virial quantities (Menon & Power 2024).
    """

    def __init__(self, omega_m=OMEGA_M_FIDUCIAL, omega_lambda=OMEGA_LAMBDA_FIDUCIAL, H0_km_s_mpc=67.8):
        if abs(omega_m + omega_lambda - 1.0) > 1e-9:
            raise ValueError("PaperCosmology is flat: omega_m + omega_lambda must equal 1")
        self.omega_m, self.omega_lambda, self.H0_km_s_mpc = float(omega_m), float(omega_lambda), float(H0_km_s_mpc)
        self.H0_per_s = self.H0_km_s_mpc * 1.0e5 / (1.0e6 * pc_cm)       # km/s/Mpc -> 1/s
        self.H0_per_Gyr = self.H0_per_s * Gyr_s

    def H_per_s(self, z):
        return self.H0_per_s * np.sqrt(self.omega_m * (1.0 + np.asarray(z, dtype=float)) ** 3 + self.omega_lambda)

    def time_at_z(self, z):
        """cosmic time in Gyr"""
        x = np.sqrt(self.omega_lambda / self.omega_m) * (1.0 + np.asarray(z, dtype=float)) ** -1.5
        return 2.0 / (3.0 * self.H0_per_Gyr * np.sqrt(self.omega_lambda)) * np.arcsinh(x)

    def z_at_time(self, t_gyr):
        arg = 1.5 * self.H0_per_Gyr * np.sqrt(self.omega_lambda) * np.asarray(t_gyr, dtype=float)
        return (np.sqrt(self.omega_m / self.omega_lambda) * np.sinh(arg)) ** (-2.0 / 3.0) - 1.0

    def __repr__(self):
        return f"PaperCosmology(omega_m={self.omega_m}, omega_lambda={self.omega_lambda}, H0={self.H0_km_s_mpc})"


def uniform_time_grid(z_start=25.0, z_end=5.0, n_steps=801, cosmology=None):
    """Redshifts of n_steps nodes uniform in cosmic time from z_start to z_end, in the given cosmology (None: Planck18 tables,
    as in the frozen production script)."""
    from .utils import z_at_time as _z_at_time_planck18
    if cosmology is None:
        t = np.linspace(_time_at_z_planck18(np.array([z_start]))[0], _time_at_z_planck18(np.array([z_end]))[0], n_steps)
        return _z_at_time_planck18(t)
    t = np.linspace(cosmology.time_at_z(z_start), cosmology.time_at_z(z_end), n_steps)
    z = cosmology.z_at_time(t)
    z[0], z[-1] = z_start, z_end          # remove round-off at the end points
    return z


def interpolate_tree_onto_grid_v2(tree_mass_msun, tree_redshift, target_redshift, cosmology=None):
    """As ``paper_reservoir.interpolate_tree_onto_grid`` (linear in cosmic time) but in the chosen cosmology."""
    tree_mass_msun = np.atleast_2d(tree_mass_msun)
    tfun = _time_at_z_planck18 if cosmology is None else cosmology.time_at_z
    t_tree, t_target = tfun(np.asarray(tree_redshift, dtype=float)), tfun(np.asarray(target_redshift, dtype=float))
    return np.array([np.interp(t_target, t_tree, row) for row in tree_mass_msun])


def momentum_agn_wind_rate_cosmo(bh_growth_rate, halo_mass, redshift, epsilon, cosmology, f_mom=1.0):
    """Same momentum-driven wind as ``reservoir_stock.momentum_agn_wind_rate`` with H(z) from ``cosmology``:
    sigma = [(5 sqrt(2) / 2) G H(z) M_halo]^(1/3), mass loading f_mom [eps / (1 - eps)] c / sigma."""
    halo_mass = np.asarray(halo_mass, dtype=float)
    resolved = halo_mass > 0
    M_g = np.where(resolved, halo_mass, 1.0) * Msun_g
    sigma = (2.5 * np.sqrt(2) * G * cosmology.H_per_s(redshift) * M_g) ** (1.0 / 3.0)
    sigma = np.where(resolved, sigma, 0.0)
    valid = sigma > 0
    load = f_mom * epsilon / (1.0 - epsilon) * c / np.where(valid, sigma, 1.0)
    return np.where(valid, load * np.asarray(bh_growth_rate, dtype=float), 0.0)


# ---------------------------------------------------------------------------------------------------------------
# The integration (copy of reservoir_stock.run_reservoir_stock with the three switches)
# ---------------------------------------------------------------------------------------------------------------
def run_reservoir_v2(
    halo_mass, redshift, seed_mass, f_b=F_B_FIDUCIAL,
    epsilon_sf=PAPER_PARAMS.epsilon_sf, epsilon=None, eta_acc=PAPER_PARAMS.eta_acc,
    R_nuc_pc=PAPER_PARAMS.R_nuc_pc, n_rd=PAPER_PARAMS.n_rd, c_nfw=PAPER_PARAMS.c_nfw,
    sigma_lnj=PAPER_PARAMS.sigma_lnj, lam_median=LAMBDA_MEDIAN_FIDUCIAL,
    M_hot=M_HOT_FIDUCIAL, phi_coldhot=PHI_COLDHOT_FIDUCIAL,
    f_mom=PAPER_PARAMS.f_mom, eta_sn_scale=PAPER_PARAMS.eta_sn_scale,
    halo_growth_rate=None,
    dm_nuclear=True, inflow="conservative", cosmology="paper", n_sub=1,
):
    """
    Integrate the revised MVM along a halo trajectory. Arguments and returned dictionary are those of
    ``reservoir_stock.run_reservoir_stock``, with:

    dm_nuclear : bool. Add the NFW dark matter inside R_nuc to M_enc(R_nuc) (threshold and nuclear free-fall time).
    inflow : "conservative" (halo mass gained on the grid; ``halo_growth_rate`` is ignored) or "rate" (frozen scheme,
        ``halo_growth_rate`` in Msun/Gyr required).
    cosmology : "paper" (a default ``PaperCosmology``), a ``PaperCosmology`` instance, or None (frozen Planck18 tables).
        ``redshift`` must already be the grid appropriate to this choice (``uniform_time_grid``).
    n_sub : int. 1 (default) is the frozen explicit step. n_sub > 1 splits the consumption and feedback part of each step
        into n_sub sub-steps (see the comment in the loop); inflow and the nuclear transfer remain once per step.

    Extra returned keys: ``cum_in`` (cumulative baryon inflow, Msun), ``halo_gain`` (sum over steps of the positive halo
    mass increments, Msun), ``halo_decrease_clipped`` (sum of the halo mass decreases that the conservative scheme does not
    deliver as negative inflow), ``budget_ratio`` = cum_in / (f_b halo_gain), at most 1 and equal to the mass-weighted
    mean suppression when the scheme is conservative.
    """
    if cosmology == "paper":
        cosmology = PaperCosmology()
    if inflow not in ("conservative", "rate"):
        raise ValueError("inflow must be 'conservative' or 'rate'")
    n_sub = int(n_sub)
    if n_sub < 1:
        raise ValueError("n_sub must be >= 1")
    halo_mass = np.atleast_2d(np.asarray(halo_mass, dtype=float))
    N, n = halo_mass.shape
    redshift = np.asarray(redshift, dtype=float)
    cosmic_time = _time_at_z_planck18(redshift) if cosmology is None else cosmology.time_at_z(redshift)
    om = OMEGA_M_FIDUCIAL if cosmology is None else cosmology.omega_m
    ol = OMEGA_LAMBDA_FIDUCIAL if cosmology is None else cosmology.omega_lambda
    H0 = 67.8 if cosmology is None else cosmology.H0_km_s_mpc
    eps_rad = PAPER_PARAMS.epsilon if epsilon is None else epsilon
    kappa_edd = float(bh_slimdisk.eddington_rate_std_per_unit_mass(eps_rad))

    if inflow == "rate":
        if halo_growth_rate is None:
            raise ValueError("inflow='rate' needs halo_growth_rate (Msun/Gyr), e.g. from paper_reservoir.grumpy_halo_growth_rate")
        halo_growth_rate = np.atleast_2d(np.asarray(halo_growth_rate, dtype=float))
    seed_mass = np.broadcast_to(np.asarray(seed_mass, dtype=float), (N,))

    def zeros():
        return np.zeros((N, n))

    bh, gas_gal, gas_nuc, gas_un = zeros(), zeros(), zeros(), zeros()
    stars_gal, stars_nuc = zeros(), zeros()
    mdot_in, mdot_out, mdot_out_sn, A_hist = zeros(), zeros(), zeros(), zeros()
    sfr_gal, sfr_nuc, tff_gal_h, tff_nuc_h = zeros(), zeros(), zeros(), zeros()
    wind_demand, w_agn, w_sn = zeros(), zeros(), zeros()
    transfer_nuc = zeros()
    lim_nuc = np.zeros((N, n), dtype=bool)
    lim_gal = np.zeros((N, n), dtype=bool)
    wcap = np.zeros((N, n), dtype=bool)
    bh[:, 0] = seed_mass

    m = zeros()
    ln_jmed = zeros()
    F_up = zeros()
    F_cur = zeros()
    cut = np.full((N, n), -np.inf)
    cum_in = np.zeros(N)
    cum_out = np.zeros(N)
    halo_gain = np.zeros(N)
    halo_clip = np.zeros(N)
    R_nuc_cm = R_nuc_pc * pc_cm

    for j in range(1, n):
        dt = cosmic_time[j] - cosmic_time[j - 1]
        z_mid = 0.5 * (redshift[j - 1] + redshift[j])
        k = j - 1
        hm = halo_mass[:, k]
        has_halo = hm > 0
        gn0, gg0 = gas_nuc[:, k], gas_gal[:, k]
        sn0, sg0, bh0 = stars_nuc[:, k], stars_gal[:, k], bh[:, k]

        # --- 1. inflow, galaxy cut (memoryless), nuclear transfer (with history) ---
        if inflow == "rate":
            rate_k = halo_growth_rate[:, k]
        else:
            dM = halo_mass[:, j] - hm                      # mass the halo gains in this step on the integration grid
            halo_gain += np.maximum(dM, 0.0)
            halo_clip += np.maximum(-dM, 0.0)
            rate_k = np.maximum(dM, 0.0) / dt              # Msun/Gyr; times dt below this is exactly the mass gained
        A_in = uv_suppression_mass_only(np.full(N, z_mid), hm) * cold_hot_mode_suppression(
            hm, M_hot=M_hot, phi=phi_coldhot) * f_b * rate_k
        mdot_in[:, j] = A_in
        fresh = A_in * dt
        cum_in += fresh
        R_vir_cm, _ = virial_radius_velocity(np.where(has_halo, hm, 1.0), z_mid, om, ol, H0)
        R_gal_cm = n_rd * lam_median * R_vir_cm / np.sqrt(2.0)
        M_dm_gal = nfw_dark_matter_mass(hm, R_gal_cm / R_vir_cm, c_nfw, f_b)
        M_dm_nuc = nfw_dark_matter_mass(hm, R_nuc_cm / R_vir_cm, c_nfw, f_b) if dm_nuclear else 0.0
        jmed = median_specific_angular_momentum(hm, z_mid, lam_median)
        ln_jmed[:, k] = np.log(np.where(jmed > 0, jmed, 1.0))
        a_gal = np.where(has_halo, ln_jthresh(M_dm_gal + bh0 + gn0 + gg0 + sn0 + sg0, R_gal_cm), -np.inf)
        f_gal = np.where(has_halo, ndtr((a_gal - ln_jmed[:, k]) / sigma_lnj), 0.0)
        m[:, k] = fresh * f_gal
        F_up[:, k] = f_gal
        gas_un[:, j] = gas_un[:, k] + fresh * (1.0 - f_gal)
        ln_jt_nuc = np.where(has_halo, ln_jthresh(M_dm_nuc + bh0 + gn0 + sn0, R_nuc_cm), -np.inf)
        T = nuclear_transfer(m, ln_jmed, F_up, F_cur, cut, ln_jt_nuc, sigma_lnj, k + 1)
        transfer_nuc[:, j] = T
        G_gal1 = m[:, :k + 1].sum(axis=1)
        G_nuc1 = gn0 + T

        if n_sub == 1:
            # --- 2. demands from the post-transfer state ---
            tff_nuc = bh_slimdisk.time_freefall_enclosed(M_dm_nuc + bh0 + G_nuc1 + sn0, R_nuc_pc=R_nuc_pc)
            tff_gal = bh_slimdisk.time_freefall_enclosed(M_dm_gal + bh0 + G_gal1 + G_nuc1 + sn0 + sg0, R_nuc_pc=R_gal_cm / pc_cm)
            tff_nuc_h[:, j], tff_gal_h[:, j] = tff_nuc, tff_gal
            with np.errstate(divide="ignore", invalid="ignore"):
                inv_nuc = np.where(np.isfinite(tff_nuc) & (tff_nuc > 0), 1.0 / tff_nuc, 0.0)
            A_bh = eta_acc * G_nuc1 * inv_nuc
            A_hist[:, j] = A_bh
            dbh_full = strict_eddington_step(bh0, A_bh, kappa_edd, dt) - bh0
            sf_n_full = star_formation_demand(G_nuc1, tff_nuc, epsilon_sf, dt)
            sf_g_full = star_formation_demand(G_gal1, tff_gal, epsilon_sf, dt)

            # --- 3. primary consumption, proportional sharing if a reservoir cannot meet its demands ---
            dem_n = dbh_full + sf_n_full
            phi_n = np.where(dem_n > G_nuc1, G_nuc1 / np.where(dem_n > 0, dem_n, 1.0), 1.0)
            phi_g = np.where(sf_g_full > G_gal1, G_gal1 / np.where(sf_g_full > 0, sf_g_full, 1.0), 1.0)
            lim_nuc[:, j], lim_gal[:, j] = phi_n < 1.0, phi_g < 1.0
            dbh, sf_n, sf_g = dbh_full * phi_n, sf_n_full * phi_n, sf_g_full * phi_g
            bh[:, j] = bh0 + dbh
            stars_nuc[:, j] = sn0 + sf_n
            stars_gal[:, j] = sg0 + sf_g
            sfr_nuc[:, j], sfr_gal[:, j] = sf_n / dt, sf_g / dt
            G_nuc2 = _settle(G_nuc1 - dbh - sf_n, G_nuc1, "nuclear gas")
            G_gal2 = _settle(G_gal1 - sf_g, G_gal1, "galaxy gas")
            m[:, :k + 1] *= np.where(G_gal1 > 0, G_gal2 / np.where(G_gal1 > 0, G_gal1, 1.0), 1.0)[:, None]

            # --- 4. feedback driven by what was achieved, removed in proportion to gas mass, capped at the gas ---
            if cosmology is None:
                W_agn = momentum_agn_wind_rate(dbh / dt, hm, z_mid, eps_rad, f_mom) * dt
            else:
                W_agn = momentum_agn_wind_rate_cosmo(dbh / dt, hm, z_mid, eps_rad, cosmology, f_mom) * dt
            W_sn = eta_sn_scale * sn.mass_loading_factor(z_mid, hm, 0.0) * (sf_n + sf_g)
            W = W_agn + W_sn
            pool = G_gal2 + G_nuc2
            take = np.minimum(W, pool)
            keep = np.where(pool > 0, 1.0 - take / np.where(pool > 0, pool, 1.0), 1.0)
            m[:, :k + 1] *= keep[:, None]
            share = np.where(W > 0, take / np.where(W > 0, W, 1.0), 0.0)
            wcap[:, j] = W > pool
            wind_demand[:, j], w_agn[:, j], w_sn[:, j] = W, W_agn * share, W_sn * share
            cum_out += take
            mdot_out[:, j] = take / dt
            mdot_out_sn[:, j] = W_sn * share / dt

            # --- 5. advance ---
            gas_nuc[:, j] = G_nuc2 * keep
            gas_gal[:, j] = G_gal2 * keep
        else:
            # --- steps 2-5 with n_sub sub-steps ---------------------------------------------------------------
            # The explicit consumption and the capped, proportional feedback removal are first-order accurate; where the
            # feedback depletes a reservoir within one step (the "wind_cap" steps) the result depends on the step. Here the
            # step is split into n_sub equal sub-steps, with the free-fall times, the Eddington step, star formation and both
            # outflows recomputed from the updated state in each. Inflow and the nuclear transfer are done once per step
            # (above). All gas cohorts of a halo are scaled by the same factor, so they are rescaled once at the end.
            dts = dt / n_sub
            Gg, Gn, bh_c, sn_c, sg_c = G_gal1.copy(), G_nuc1.copy(), bh0.copy(), sn0.copy(), sg0.copy()
            dbh_t, sfn_t, sfg_t = np.zeros(N), np.zeros(N), np.zeros(N)
            Wd_t, Wa_t, Ws_t, take_t = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
            A_t, tn_first, tg_first = np.zeros(N), None, None
            any_ln, any_lg, any_wc = np.zeros(N, dtype=bool), np.zeros(N, dtype=bool), np.zeros(N, dtype=bool)
            for _s in range(n_sub):
                tff_nuc = bh_slimdisk.time_freefall_enclosed(M_dm_nuc + bh_c + Gn + sn_c, R_nuc_pc=R_nuc_pc)
                tff_gal = bh_slimdisk.time_freefall_enclosed(M_dm_gal + bh_c + Gg + Gn + sn_c + sg_c, R_nuc_pc=R_gal_cm / pc_cm)
                if tn_first is None:
                    tn_first, tg_first = tff_nuc, tff_gal
                with np.errstate(divide="ignore", invalid="ignore"):
                    inv_nuc = np.where(np.isfinite(tff_nuc) & (tff_nuc > 0), 1.0 / tff_nuc, 0.0)
                A_bh = eta_acc * Gn * inv_nuc
                A_t += A_bh / n_sub
                dbh_full = strict_eddington_step(bh_c, A_bh, kappa_edd, dts) - bh_c
                sf_n_full = star_formation_demand(Gn, tff_nuc, epsilon_sf, dts)
                sf_g_full = star_formation_demand(Gg, tff_gal, epsilon_sf, dts)
                dem_n = dbh_full + sf_n_full
                phi_n = np.where(dem_n > Gn, Gn / np.where(dem_n > 0, dem_n, 1.0), 1.0)
                phi_g = np.where(sf_g_full > Gg, Gg / np.where(sf_g_full > 0, sf_g_full, 1.0), 1.0)
                any_ln |= phi_n < 1.0
                any_lg |= phi_g < 1.0
                dbh, sf_n, sf_g = dbh_full * phi_n, sf_n_full * phi_n, sf_g_full * phi_g
                bh_c, sn_c, sg_c = bh_c + dbh, sn_c + sf_n, sg_c + sf_g
                Gn2 = _settle(Gn - dbh - sf_n, Gn, "nuclear gas")
                Gg2 = _settle(Gg - sf_g, Gg, "galaxy gas")
                if cosmology is None:
                    W_agn = momentum_agn_wind_rate(dbh / dts, hm, z_mid, eps_rad, f_mom) * dts
                else:
                    W_agn = momentum_agn_wind_rate_cosmo(dbh / dts, hm, z_mid, eps_rad, cosmology, f_mom) * dts
                W_sn = eta_sn_scale * sn.mass_loading_factor(z_mid, hm, 0.0) * (sf_n + sf_g)
                W = W_agn + W_sn
                pool = Gg2 + Gn2
                take = np.minimum(W, pool)
                keep = np.where(pool > 0, 1.0 - take / np.where(pool > 0, pool, 1.0), 1.0)
                share = np.where(W > 0, take / np.where(W > 0, W, 1.0), 0.0)
                any_wc |= W > pool
                dbh_t += dbh; sfn_t += sf_n; sfg_t += sf_g
                Wd_t += W; Wa_t += W_agn * share; Ws_t += W_sn * share; take_t += take
                Gn, Gg = Gn2 * keep, Gg2 * keep
            tff_nuc_h[:, j], tff_gal_h[:, j] = tn_first, tg_first
            A_hist[:, j] = A_t
            lim_nuc[:, j], lim_gal[:, j], wcap[:, j] = any_ln, any_lg, any_wc
            bh[:, j] = bh_c
            stars_nuc[:, j], stars_gal[:, j] = sn_c, sg_c
            sfr_nuc[:, j], sfr_gal[:, j] = sfn_t / dt, sfg_t / dt
            m[:, :k + 1] *= np.where(G_gal1 > 0, Gg / np.where(G_gal1 > 0, G_gal1, 1.0), 1.0)[:, None]
            wind_demand[:, j], w_agn[:, j], w_sn[:, j] = Wd_t, Wa_t, Ws_t
            cum_out += take_t
            mdot_out[:, j] = take_t / dt
            mdot_out_sn[:, j] = Ws_t / dt
            gas_nuc[:, j] = Gn
            gas_gal[:, j] = Gg

    baryon_now = gas_gal[:, -1] + gas_nuc[:, -1] + gas_un[:, -1] + stars_gal[:, -1] + stars_nuc[:, -1] \
        + (bh[:, -1] - bh[:, 0]) + cum_out
    residual = (baryon_now - cum_in) / np.maximum(cum_in, 1.0)

    resolved = halo_mass > 0
    first = np.argmax(resolved, axis=1)
    formed = resolved.any(axis=1)
    host_baryons = f_b * halo_mass
    with np.errstate(divide="ignore", invalid="ignore"):
        bh_to_host = np.where(resolved, bh / host_baryons, np.nan)
        seed_at_formation = np.where(formed, bh[:, 0] / host_baryons[np.arange(N), first], np.nan)
        budget_ratio = np.where(halo_gain > 0, cum_in / (f_b * halo_gain), np.nan)
    z_first = np.where(formed, redshift[first], np.nan)

    return dict(
        transfer_nuc=transfer_nuc, bh_to_host_baryons=bh_to_host, seed_over_host_baryons_at_formation=seed_at_formation, z_first_resolved=z_first,
        bh_mass=bh, gas_gal=gas_gal, gas_nuc=gas_nuc, gas_unavailable=gas_un, stars_gal=stars_gal, stars_nuc=stars_nuc,
        gas_mass=gas_gal + gas_nuc, stars_mass=stars_gal + stars_nuc,
        halo_mass=halo_mass, cosmic_time=cosmic_time, redshift=redshift,
        mdot_in=mdot_in, mdot_out=mdot_out, mdot_out_sn=mdot_out_sn, A_bh=A_hist,
        sfr_gal=sfr_gal, sfr_nuc=sfr_nuc, t_ff_gal=tff_gal_h, t_ff_nuc=tff_nuc_h, kappa_edd=kappa_edd,
        limiter_nuc=lim_nuc, limiter_gal=lim_gal, wind_cap=wcap, wind_demand=wind_demand,
        wind_taken_agn=w_agn, wind_taken_sn=w_sn, baryon_residual=residual,
        cum_in=cum_in, halo_gain=halo_gain, halo_decrease_clipped=halo_clip, budget_ratio=budget_ratio,
    )


# ---------------------------------------------------------------------------------------------------------------
# Critical seed and the single-crossing check
# ---------------------------------------------------------------------------------------------------------------
def critical_seed_v2(halo_mass, redshift, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=None,
                     n_iter=40, return_F=False, **run_kwargs):
    """
    M_seed,crit by bisection in ln M_seed on F = M_BH / (f_BH M_star,tot) - 1 at the last grid point; identical logic and
    bracket to ``reservoir_stock.critical_seed_stock``. Returns (M_seed_crit, never_reaches_target, already_above_at_lo)
    (plus F at the solution if return_F). Assumes F crosses zero once: use ``single_crossing_report`` to test that.
    """
    halo_mass = np.atleast_2d(halo_mass)
    N = halo_mass.shape[0]
    lo = np.full(N, float(seed_mass_lo))
    hi = (2.0 * f_bh * F_B_FIDUCIAL * halo_mass[:, -1]) if seed_mass_hi is None else np.full(N, float(seed_mass_hi))

    def F(seed):
        return target_residual(run_reservoir_v2(halo_mass, redshift, seed, **run_kwargs), f_bh)

    above, never = F(lo) > 0.0, F(hi) < 0.0
    for _ in range(n_iter):
        mid = np.sqrt(lo * hi)
        neg = F(mid) < 0.0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    M = np.where(above | never, np.nan, np.sqrt(lo * hi))
    if return_F:
        return M, never, above, F(np.where(np.isnan(M), lo, M))
    return M, never, above


def single_crossing_report(halo_mass, redshift, f_bh=PAPER_PARAMS.f_bh, n_grid=33, seed_mass_lo=1e-2, **run_kwargs):
    """
    Test, for every tree, that F(M_seed) crosses zero exactly once, on a log-spaced seed grid between the bisection bracket
    limits. Returns per-tree arrays: n_sign_changes (the quantity that matters for the bisection), max_decrease (largest fall
    of F between adjacent grid seeds; informational: F can fall at seeds far above the root, where the black hole's own
    feedback removes gas, without affecting the root), the seeds and F. A grid cannot exclude a crossing pair narrower than
    its spacing, so use n_grid >= 33 (about 0.9 dex per interval over 1e-2 to 1e12 Msun) and treat the result as a check,
    not a proof.
    """
    halo_mass = np.atleast_2d(halo_mass)
    N = halo_mass.shape[0]
    hi = 2.0 * f_bh * F_B_FIDUCIAL * halo_mass[:, -1]
    frac = np.linspace(0.0, 1.0, n_grid)
    seeds = np.exp(np.log(seed_mass_lo) + frac[None, :] * (np.log(hi) - np.log(seed_mass_lo))[:, None])
    Fm = np.empty((N, n_grid))
    for i in range(n_grid):
        Fm[:, i] = target_residual(run_reservoir_v2(halo_mass, redshift, seeds[:, i], **run_kwargs), f_bh)
    sign = np.sign(Fm)
    n_changes = (sign[:, 1:] != sign[:, :-1]).sum(axis=1)
    max_decrease = np.maximum(0.0, (Fm[:, :-1] - Fm[:, 1:])).max(axis=1)
    return dict(n_sign_changes=n_changes, max_decrease=max_decrease, seeds=seeds, F=Fm)
