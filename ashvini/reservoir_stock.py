"""
Minimal Viable Model (MVM) of the galaxy and nuclear gas reservoirs.

Question: for a halo of specified M_halo(z=5), what is the minimum black hole seed mass at
z_seed = 25 that reaches M_BH = f_BH M_star,tot by z = 5 when black hole growth is restricted to
Eddington-limited gas accretion? The model isolates the interaction of (1) stochastic halo assembly,
(2) galaxy-scale gas supply and retention, (3) angular-momentum-dependent delivery of gas to a nuclear
reservoir, (4) competition for that gas, and (5) the finite Eddington growth time. It is not a model of
the ISM, of a nuclear star cluster, or of sub-pc fuelling. The pre-MVM exploratory model is frozen in
``ashvini.reservoir_stock_premvm``.

Scientific content
------------------
* Halo: a Zhang & Hui main-progenitor tree, interpolated in cosmic time, with the GRUMPY rate.
* Inflow: Mdot_in = zeta_UV zeta_ch f_b Mdot_halo.
* Angular momentum: each step's inflow is a cohort with j ~ LogNormal(ln j_med, sigma_j),
  j_med = lambda sqrt(2) V_vir R_vir. Specific angular momentum is conserved.
* Galaxy cut (memoryless): the fraction of a cohort with j <= j_thresh(R_gal) enters the galaxy gas
  reservoir; the rest is unavailable (booked in a ledger, no dynamics, never returns).
  j_thresh(R) = sqrt(G M_enc(R) R), R_gal = n_rd R_d, R_d = lambda R_vir / sqrt(2),
  M_enc(R_gal) = M_DM,NFW(<R_gal) + M_BH + all gas + all stars.
* Nuclear transfer (with history): galaxy gas with j <= j_thresh(R_nuc) moves to the nuclear reservoir,
  M_enc(R_nuc) = M_BH + M_gas,nuc + M_star,nuc (no dark matter; R_nuc is fixed). Each cohort keeps the
  running maximum of the nuclear threshold, so gas that was inaccessible becomes accessible if
  j_thresh(R_nuc) later rises above its j.
* Two evolved gas reservoirs only: M_gas,gal and M_gas,nuc.
* One star-formation law at both scales, Mdot_star = eps_sf M_gas / t_ff(R), t_ff = sqrt(R^3 / 2 G M_enc(R)).
* One unresolved galaxy-scale stellar outflow, Mdot_out,star = eta(M_halo, z) Mdot_star,tot (metal-poor
  limit, instantaneous), removed in proportion to gas mass from the two reservoirs.
* Black hole: Mdot_acc = eta_acc M_gas,nuc / t_ff,nuc; Mdot_BH = min(Mdot_acc, kappa_edd M_BH) (exact
  strict-Eddington step; no compaction factor; no super-Eddington growth).
* AGN outflow: momentum-driven, Mdot_out v_out = f_mom L / c, L = [eps/(1-eps)] Mdot_BH c^2, v_out = sigma;
  removed from the two reservoirs in proportion to gas mass together with the stellar outflow.
* No recycling: outflow mass leaves the system.
* Excluded: BH mergers and dynamical seeding, torques or angular-momentum loss, disc dynamics, spin,
  metallicity evolution, delayed feedback, clumps, multiphase ISM, wind recycling, nuclear dark matter,
  moving R_nuc, a separate nuclear feedback model.
* Initial conditions: empty galaxy, the seed black hole only, at z_seed.

Numerical convention (one explicit step, in causal order)
---------------------------------------------------------
1. Inflow and transfers. Fresh inflow is split by the galaxy cut; the nuclear transfer is evaluated
   with thresholds computed from the state at the start of the step (the last established state).
2. Demands. From the resulting post-transfer state, t_ff at both scales, the black hole step, the
   nuclear and galaxy star formation are evaluated.
3. Primary consumption. Black hole growth and nuclear star formation draw on the nuclear gas, galaxy
   star formation on the galaxy gas. If a reservoir's demands exceed its gas they are scaled down together
   (proportional sharing; counted in ``limiter_nuc`` and ``limiter_gal``).
4. Feedback. The stellar and AGN outflows are driven by the star formation and black hole growth
   actually achieved in step 3. Their combined demand is removed from the gas that remains, in proportion
   to gas mass, and is capped at that gas: feedback can consume the whole reservoir, which is a
   consequence of the prescription (counted in ``wind_cap``); a negative mass is an error.
5. Advance.
"""

import numpy as np
from scipy.special import ndtr

from . import black_holes_growth_slimdisk as bh_slimdisk
from . import supernovae_feedback as sn
from .constants import G, Msun_g, pc_cm, c
from .paper_reservoir import (
    F_B_FIDUCIAL, M_HOT_FIDUCIAL, PHI_COLDHOT_FIDUCIAL, LAMBDA_MEDIAN_FIDUCIAL,
    OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL,
    virial_radius_velocity, uv_suppression_mass_only, cold_hot_mode_suppression, smooth_halo_growth_rate,
)
from .paper_reservoir_params import PAPER_PARAMS
from .utils import time_at_z


def strict_eddington_step(M0, A, kappa, dt):
    """
    Exact update of dM/dt = min(A, kappa M) over dt for constant A and kappa: strict Eddington-limited
    growth with supply limit A. Growth is exponential while kappa M < A and linear at rate A afterwards.
    Note that ``bh_growth_step_slimdisk(..., r_crit=1)`` is NOT this: with r_crit = 1 it returns the raw supply.
    """
    M0 = np.asarray(M0, dtype=float)
    A = np.asarray(A, dtype=float)
    has_supply = A > 0
    A_s = np.where(has_supply, A, 1.0)
    supply_limited = kappa * M0 >= A_s
    positive = M0 > 0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        t_star = np.where(positive, np.log(A_s / (kappa * np.where(positive, M0, 1.0))) / kappa, np.inf)
        grown = np.where(t_star >= dt, M0 * np.exp(kappa * dt), A_s / kappa + A_s * (dt - t_star))
    out = np.where(supply_limited, M0 + A_s * dt, np.where(positive, grown, M0))
    return np.where(has_supply, out, M0)


def momentum_agn_wind_rate(bh_growth_rate, halo_mass, redshift, epsilon, f_mom=1.0):
    """
    Momentum-driven AGN outflow (King 2003): Mdot_out v_out = f_mom L / c with
    L = [epsilon/(1-epsilon)] Mdot_BH c^2 and v_out = sigma (isothermal), so the mass loading is
    f_mom [epsilon/(1-epsilon)] c / sigma. Zero for an unformed halo.
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


def ln_jthresh(M_enc_msun, R_cm):
    """ln sqrt(G M_enc R) (cgs); -inf where M_enc <= 0."""
    M = np.maximum(np.asarray(M_enc_msun, dtype=float), 0.0)
    with np.errstate(divide="ignore"):
        out = 0.5 * (np.log(G) + np.log(np.where(M > 0, M, 1.0) * Msun_g) + np.log(R_cm))
    return np.where(M > 0, out, -np.inf)


def star_formation_demand(gas, t_ff, epsilon, dt):
    """Mass a reservoir would turn into stars in dt, epsilon M_gas / t_ff (0 where t_ff is infinite)."""
    t_ff = np.asarray(t_ff, dtype=float)
    ok = np.isfinite(t_ff) & (t_ff > 0)
    return np.where(ok, epsilon * np.asarray(gas, dtype=float) / np.where(ok, t_ff, 1.0) * dt, 0.0)


def nuclear_transfer(m, ln_jmed, F_up, F_cur, cut, ln_jt_nuc, sigma_lnj, n_active):
    """
    Move to the nucleus the galaxy gas whose j lies below the nuclear threshold, keeping each cohort's
    history.

    Cohort k arrived with j ~ LogNormal(ln_jmed[k], sigma_lnj) truncated at j < j_gal (the galaxy cut then),
    so its gas has j in (exp(cut[k]), exp(ln a_k)) with CDF values (F_cur[k], F_up[k]) for the untruncated
    distribution. Raising the nuclear cut to b moves the fraction (F(b) - F_cur) / (F_up - F_cur) of what is
    left; a cut never decreases, so a later fall in the threshold moves nothing and a later rise releases
    previously inaccessible gas. ``m``, ``F_cur`` and ``cut`` are updated in place; returns the mass moved
    per halo.
    """
    N = m.shape[0]
    T = np.zeros(N)
    sl = slice(0, n_active)
    need = (ln_jt_nuc[:, None] > cut[:, sl]) & (m[:, sl] > 0)
    ii, cc = np.nonzero(need)
    if ii.size == 0:
        return T
    b = ln_jt_nuc[ii]
    up = F_up[ii, cc]
    old = F_cur[ii, cc]
    new = np.minimum(ndtr((b - ln_jmed[ii, cc]) / sigma_lnj), up)
    denom = up - old
    frac = np.clip(np.where(denom > 0, (new - old) / np.where(denom > 0, denom, 1.0), 1.0), 0.0, 1.0)
    moved = m[ii, cc] * frac
    m[ii, cc] -= moved
    F_cur[ii, cc] = np.maximum(new, old)
    cut[ii, cc] = b
    np.add.at(T, ii, moved)
    return T


def _settle(x, ref, name):
    """Zero round-off-level negatives; anything larger is a bug and raises."""
    if np.any(x < -1e-9 * np.maximum(np.abs(ref), 1.0)):
        raise FloatingPointError(f"negative {name}: min {np.min(x):.3e}")
    return np.maximum(x, 0.0)


def run_reservoir_stock(
    halo_mass, redshift, seed_mass, f_b=F_B_FIDUCIAL,
    epsilon_sf=PAPER_PARAMS.epsilon_sf, epsilon=None, eta_acc=PAPER_PARAMS.eta_acc,
    R_nuc_pc=PAPER_PARAMS.R_nuc_pc, n_rd=PAPER_PARAMS.n_rd, c_nfw=PAPER_PARAMS.c_nfw,
    sigma_lnj=PAPER_PARAMS.sigma_lnj, lam_median=LAMBDA_MEDIAN_FIDUCIAL,
    M_hot=M_HOT_FIDUCIAL, phi_coldhot=PHI_COLDHOT_FIDUCIAL,
    f_mom=PAPER_PARAMS.f_mom, eta_sn_scale=PAPER_PARAMS.eta_sn_scale,
    halo_growth_rate=None,
):
    """
    Integrate the MVM (see the module docstring) along a given halo trajectory.

    halo_mass : (N, n) or (n,) Msun on the shared chronological grid ``redshift`` (n,); the first step is
        z_seed. seed_mass : black hole mass at the first step. halo_growth_rate : (N, n) Msun/Gyr, e.g. from
        ``paper_reservoir.grumpy_halo_growth_rate`` (default: smooth closed form, valid only for smooth
        trajectories).

    Returns a dict: bh_mass, gas_gal, gas_nuc, gas_unavailable (ledger), stars_gal, stars_nuc, gas_mass
    (= gal + nuc), stars_mass (= gal + nuc), halo_mass, cosmic_time, redshift, mdot_in, mdot_out (AGN + stellar),
    mdot_out_sn, A_bh, sfr_gal, sfr_nuc, t_ff_gal, t_ff_nuc (Gyr), kappa_edd, plus per-step diagnostics
    limiter_nuc, limiter_gal, wind_cap (booleans), wind_demand, wind_taken_agn, wind_taken_sn (Msun) and
    baryon_residual (relative, ~1e-15). Output-only diagnostics: transfer_nuc (mass moved to the nuclear reservoir per step), bh_to_host_baryons (M_BH / f_b M_halo per step, NaN
    before the halo exists), seed_over_host_baryons_at_formation (seed / f_b M_halo at the first resolved step)
    and z_first_resolved.
    """
    halo_mass = np.atleast_2d(np.asarray(halo_mass, dtype=float))
    N, n = halo_mass.shape
    redshift = np.asarray(redshift, dtype=float)
    cosmic_time = time_at_z(redshift)
    eps_rad = PAPER_PARAMS.epsilon if epsilon is None else epsilon
    kappa_edd = float(bh_slimdisk.eddington_rate_std_per_unit_mass(eps_rad))

    if halo_growth_rate is None:
        halo_growth_rate = smooth_halo_growth_rate(halo_mass, redshift)
    else:
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

    # galaxy-gas cohorts: mass left, ln j_med at arrival, CDF at the galaxy cut then, CDF and ln of the nuclear cut
    m = zeros()
    ln_jmed = zeros()
    F_up = zeros()
    F_cur = zeros()
    cut = np.full((N, n), -np.inf)
    cum_in = np.zeros(N)
    cum_out = np.zeros(N)
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
        A_in = uv_suppression_mass_only(np.full(N, z_mid), hm) * cold_hot_mode_suppression(
            hm, M_hot=M_hot, phi=phi_coldhot) * f_b * halo_growth_rate[:, k]
        mdot_in[:, j] = A_in
        fresh = A_in * dt
        cum_in += fresh
        R_vir_cm, _ = virial_radius_velocity(np.where(has_halo, hm, 1.0), z_mid, OMEGA_M_FIDUCIAL, OMEGA_LAMBDA_FIDUCIAL, 67.8)
        R_gal_cm = n_rd * lam_median * R_vir_cm / np.sqrt(2.0)
        M_dm_gal = nfw_dark_matter_mass(hm, R_gal_cm / R_vir_cm, c_nfw, f_b)
        jmed = median_specific_angular_momentum(hm, z_mid, lam_median)
        ln_jmed[:, k] = np.log(np.where(jmed > 0, jmed, 1.0))
        a_gal = np.where(has_halo, ln_jthresh(M_dm_gal + bh0 + gn0 + gg0 + sn0 + sg0, R_gal_cm), -np.inf)
        f_gal = np.where(has_halo, ndtr((a_gal - ln_jmed[:, k]) / sigma_lnj), 0.0)
        m[:, k] = fresh * f_gal
        F_up[:, k] = f_gal
        gas_un[:, j] = gas_un[:, k] + fresh * (1.0 - f_gal)
        ln_jt_nuc = np.where(has_halo, ln_jthresh(bh0 + gn0 + sn0, R_nuc_cm), -np.inf)
        T = nuclear_transfer(m, ln_jmed, F_up, F_cur, cut, ln_jt_nuc, sigma_lnj, k + 1)
        transfer_nuc[:, j] = T
        G_gal1 = m[:, :k + 1].sum(axis=1)
        G_nuc1 = gn0 + T

        # --- 2. demands from the post-transfer state ---
        tff_nuc = bh_slimdisk.time_freefall_enclosed(bh0 + G_nuc1 + sn0, R_nuc_pc=R_nuc_pc)
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
        W_agn = momentum_agn_wind_rate(dbh / dt, hm, z_mid, eps_rad, f_mom) * dt
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

    baryon_now = gas_gal[:, -1] + gas_nuc[:, -1] + gas_un[:, -1] + stars_gal[:, -1] + stars_nuc[:, -1] \
        + (bh[:, -1] - bh[:, 0]) + cum_out
    residual = (baryon_now - cum_in) / np.maximum(cum_in, 1.0)

    # output-only diagnostics: how heavy is the black hole compared with the baryons its host has accreted?
    resolved = halo_mass > 0
    first = np.argmax(resolved, axis=1)
    formed = resolved.any(axis=1)
    host_baryons = f_b * halo_mass
    with np.errstate(divide="ignore", invalid="ignore"):
        bh_to_host = np.where(resolved, bh / host_baryons, np.nan)
        seed_at_formation = np.where(formed, bh[:, 0] / host_baryons[np.arange(N), first], np.nan)
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
    )


def target_residual(result, f_bh):
    """F = M_BH(z=5) / (f_BH M_star,tot(z=5)) - 1, per halo."""
    return result["bh_mass"][:, -1] / (f_bh * result["stars_mass"][:, -1]) - 1.0


def critical_seed_stock(halo_mass, redshift, f_bh=PAPER_PARAMS.f_bh, seed_mass_lo=1e-2, seed_mass_hi=None,
                        n_iter=40, return_F=False, **run_kwargs):
    """
    M_seed,crit by bisection (in ln M_seed) on F(M_seed) = M_BH / (f_BH M_star,tot) - 1 at z = 5.

    The default upper bracket, 2 f_BH f_b M_halo(z=5), always has F > 0: M_star,tot cannot exceed the baryon
    supply f_b M_halo(z=5), and M_BH >= M_seed. Returns (M_seed_crit, never_reaches_target, already_above_at_lo)
    with the same semantics as ``paper_reservoir.critical_seed_paper``; NaN where the bracket fails. With
    return_F=True also returns F at the solution. Bisection assumes F crosses zero once, from below; check
    with ``target_residual`` on a seed grid if in doubt.
    """
    halo_mass = np.atleast_2d(halo_mass)
    N = halo_mass.shape[0]
    lo = np.full(N, float(seed_mass_lo))
    hi = (2.0 * f_bh * F_B_FIDUCIAL * halo_mass[:, -1]) if seed_mass_hi is None else np.full(N, float(seed_mass_hi))

    def F(seed):
        return target_residual(run_reservoir_stock(halo_mass, redshift, seed, **run_kwargs), f_bh)

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
