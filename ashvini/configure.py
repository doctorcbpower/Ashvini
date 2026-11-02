"""
Runtime parameter overrides for parameter-sweep notebooks/scripts.

The parameters live in one place, `config.DEFAULT` (the Config loaded from run_params.yaml); nothing is
copied into module-level variables, so changing a value here is seen by every module on its next call.
`set_params` changes values of `config.DEFAULT` in place, for the remainder of the process (until the next
call); it is the interactive wrapper around what a script should do instead when it can: build a Config with
`config.override(...)` and pass it to `main.run_forest(..., config=cfg)`, which leaves `DEFAULT` untouched
and allows several parameter sets in one process.

A parameter sweep (e.g. reproducing Menon & Power 2024's Figures 4-7, which vary
epsilon_sf/epsilon_fb/t_d/z_rei one at a time against a fixed fiducial) is the typical use.
"""
from .config import DEFAULT
from .imf import nu_SN


def set_params(
    eps_sf=None,
    eps_fb=None,
    t_d=None,
    z_rei=None,
    sn_type=None,
    uv_background=None,
    bh_enabled=None,
    eta_agn=None,
    sigma_feedback_enabled=None,
    growth_cap_enabled=None,
    e_bh=None,
    eddington_multiplier=None,
    m_sn_max=None,
    agn_delay_time=None,
):
    """
    Override one or more fiducial parameters in place, for the remainder of
    the process (until the next set_params() call). Any argument left as
    None keeps its current value. Returns nothing; call print_current() to
    inspect the resulting state.

    eps_sf / eps_fb are star_formation.efficiency and supernova.epsilon_p; t_d is supernova.delay_time (Gyr);
    z_rei is reionization.z_reion; m_sn_max is the upper progenitor mass for supernovae (it changes only the
    supernova rate in the dust destruction term); agn_delay_time is black_holes.feedback_delay_time (Gyr).

    bh_enabled/eta_agn/sigma_feedback_enabled/growth_cap_enabled/e_bh
    control black hole growth and AGN feedback (see MODELS.md's "Black hole
    growth"/"AGN feedback" sections). growth_cap_enabled is only
    physically meaningful with sigma_feedback_enabled=True too (see
    black_holes_growth.growth_ceiling) but isn't validated against it
    here -- set both explicitly if you want the hard-cap self-regulation
    variant.

    bh_enabled=False sets every seeding channel's `enabled` to False (no
    new BH is ever seeded, so growth/feedback are moot regardless of
    eta_agn/eddington_multiplier -- this is the only way to get a genuine
    "no BH" baseline, since run_forest()/run1_scalar() always compute BH
    seeding/growth/AGN terms unconditionally). bh_enabled=True restores
    all three channels to enabled -- if you'd previously disabled only
    some channels by hand via the seeding config directly, this will
    re-enable all of them, not just the ones you disabled.
    """
    cfg = DEFAULT
    if eps_sf is not None:
        cfg.sf.efficiency = eps_sf
    if eps_fb is not None:
        cfg.sn.epsilon_p = eps_fb
    if t_d is not None:
        cfg.sn.delay_time = t_d
    if sn_type is not None:
        cfg.sn.type = sn_type
    if uv_background is not None:
        cfg.reion.UVB_enabled = uv_background
    if z_rei is not None:
        cfg.reion.z_reion = z_rei
    if bh_enabled is not None:
        seeding = cfg.bh.seeding
        seeding.pop3.enabled = bh_enabled
        seeding.direct_collapse.enabled = bh_enabled
        seeding.halo_mass_threshold.enabled = bh_enabled
    if eta_agn is not None:
        cfg.bh.eta_agn = eta_agn
    if sigma_feedback_enabled is not None:
        cfg.bh.sigma_feedback.enabled = sigma_feedback_enabled
    if growth_cap_enabled is not None:
        cfg.bh.sigma_feedback.growth_cap_enabled = growth_cap_enabled
    if e_bh is not None:
        cfg.bh.efficiency = e_bh
    if eddington_multiplier is not None:
        cfg.bh.eddington_multiplier = eddington_multiplier
    if m_sn_max is not None:
        nu_SN(cfg.sn.imf, cfg.sn.m_sn_min, m_sn_max)  # validates the range
        cfg.sn.m_sn_max = m_sn_max
    if agn_delay_time is not None:
        cfg.bh.feedback_delay_time = agn_delay_time


def print_current():
    cfg = DEFAULT
    seeding = cfg.bh.seeding
    print(
        f"eps_sf={cfg.sf.efficiency}, eps_fb={cfg.sn.epsilon_p}, "
        f"t_d={cfg.sn.delay_time} Gyr, sn_type={cfg.sn.type!r}, "
        f"UV_background={cfg.reion.UVB_enabled}, z_rei={cfg.reion.z_reion}, m_sn_max={cfg.sn.m_sn_max}\n"
        f"BH seeding enabled: pop3={seeding.pop3.enabled}, "
        f"direct_collapse={seeding.direct_collapse.enabled}, "
        f"halo_mass_threshold={seeding.halo_mass_threshold.enabled}; "
        f"eta_agn={cfg.bh.eta_agn}, sigma_feedback_enabled={cfg.bh.sigma_feedback.enabled}, "
        f"growth_cap_enabled={cfg.bh.sigma_feedback.growth_cap_enabled}, e_bh={cfg.bh.efficiency}, "
        f"eddington_multiplier={cfg.bh.eddington_multiplier}, agn_delay_time={cfg.bh.feedback_delay_time}"
    )
