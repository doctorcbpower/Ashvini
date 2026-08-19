"""
Convert a critical-seed-mass boundary M_seed,crit(M_halo) (e.g. from
ashvini.black_holes_growth_slimdisk / the "hobbs_slimdisk" growth model,
or the standalone seed_halo_boundary.py analysis script in
AshviniPapers/High_z_Black_Hole_Growth/analysis/) into a predicted
fraction of a given seed-formation channel that requires a
super-Eddington episode, given a log-normal-in-log10(mass) population
model for that channel's seed mass.

Illustrative, not a rigorously derived population synthesis -- see
docs/2026_paper_session_code_catalogue.md and the paper's Section 5.3/
Section 6.5 caveat. Ported from the standalone prototype
(seed_halo_boundary.py's later half); this is the one piece of that
script judged small and general enough to live in the package rather
than only in the paper's analysis archive.
"""

import numpy as np
from scipy.stats import norm


# mu = log10(median seed mass / Msun), sigma = dex scatter. See
# Inayoshi, Visbal & Haiman (2020) for the review this is anchored to;
# DevecchiVolonteri2009/LodatoNatarajan2006 for the runaway-collision /
# direct-collapse channel numbers specifically. All three need a
# from-scratch literature check before being treated as more than
# illustrative -- see the paper's caveats.
DEFAULT_CHANNELS = {
    "pop3": dict(mu=np.log10(150.0), sigma=0.35),
    "runaway_collision": dict(mu=np.log10(2000.0), sigma=0.30),
    "direct_collapse": dict(mu=np.log10(2.0e5), sigma=0.30),
}


def fraction_requiring_boost(M_seed_crit, mu, sigma):
    """
    Fraction of a log-normal-in-log10(mass) seed population below a
    critical seed mass M_seed_crit (Msun) -- i.e. the fraction requiring
    a super-Eddington episode to reach the boundary's target f_BH by the
    boundary's anchor redshift.

    mu, sigma : median log10(seed mass/Msun) and dex scatter of the
        channel's seed-mass distribution.
    """
    M_seed_crit = np.asarray(M_seed_crit, dtype=float)
    with np.errstate(divide="ignore"):
        log_M_crit = np.where(M_seed_crit > 0, np.log10(np.where(M_seed_crit > 0, M_seed_crit, 1.0)), -np.inf)
    return norm.cdf(log_M_crit, loc=mu, scale=sigma)


def boosted_fractions_by_channel(M_seed_crit_of_halo_mass, channels=None):
    """
    M_seed_crit_of_halo_mass : dict or array-like mapping halo mass ->
        M_seed_crit (Msun), e.g. Table~1 of the 2026 paper.
    channels : dict of {name: dict(mu=..., sigma=...)}, default
        DEFAULT_CHANNELS above.

    Returns dict {channel_name: {halo_mass: fraction_requiring_boost}}.
    """
    channels = channels or DEFAULT_CHANNELS
    if hasattr(M_seed_crit_of_halo_mass, "items"):
        halo_masses = list(M_seed_crit_of_halo_mass.keys())
        M_crit_vals = np.array([M_seed_crit_of_halo_mass[m] for m in halo_masses])
    else:
        halo_masses = list(range(len(M_seed_crit_of_halo_mass)))
        M_crit_vals = np.asarray(M_seed_crit_of_halo_mass, dtype=float)

    out = {}
    for name, params in channels.items():
        fracs = fraction_requiring_boost(M_crit_vals, params["mu"], params["sigma"])
        out[name] = dict(zip(halo_masses, fracs))
    return out
