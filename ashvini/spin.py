"""
Black hole spin -> radiative efficiency, via the Novikov & Thorne (1973)
innermost-stable-circular-orbit (ISCO) binding energy.

Used by black_holes_growth_slimdisk.py (the "hobbs_slimdisk" growth
model, see MODELS.md) to convert a spin parameter a* into the radiative
efficiency epsilon that sets the standard Eddington rate
(t_Sal = [epsilon/(1-epsilon)] * t_Edd). Reference values (see the
2026 differential-growth paper, Section 5.5):

    a* = 0     (King & Pringle 2006 chaotic-accretion floor) -> epsilon = 0.0572
    a* = 0.5   (fiducial)                                    -> epsilon = 0.0821
    a* = 0.9   (coherent-accretion-motivated pessimistic edge)-> epsilon = 0.1558
    a* = 0.998 (maximal prograde, Thorne 1974 bound)          -> epsilon = 0.3210
"""

import numpy as np


def isco_radius(a_star):
    """
    ISCO radius r_isco(a*) in units of gravitational radius r_g = GM/c^2
    (Bardeen, Press & Teukolsky 1972, eq. 2.21; Novikov & Thorne 1973).

    a_star in [-1, 1]: a_star >= 0 is a prograde (co-rotating) orbit at
    spin magnitude a_star; a_star < 0 is a retrograde (counter-rotating)
    orbit at spin magnitude |a_star|. Z1 and Z2 below are even functions
    of a (unchanged by a -> -a), so the prograde/retrograde distinction
    comes entirely from which root is taken -- "-" for prograde, "+" for
    retrograde -- not from the sign of a_star inside Z1/Z2 itself; this is
    handled explicitly (`sign` below) rather than by naively passing a
    signed a_star through a single fixed-sign formula, which would
    silently return the *prograde* value at |a_star| for a retrograde
    input (checked directly against known limits: retrograde extremal,
    a_star=-1, must give r_isco=9, not 1):

        Z1 = 1 + (1-a^2)^(1/3) * [(1+a)^(1/3) + (1-a)^(1/3)],  a = |a_star|
        Z2 = sqrt(3*a^2 + Z1^2)
        r_isco/r_g = 3 + Z2 + sign*sqrt((3-Z1)*(3+Z1+2*Z2))
        sign = -1 (prograde, a_star >= 0) or +1 (retrograde, a_star < 0)
    """
    a_star = np.asarray(a_star, dtype=float)
    a = np.abs(a_star)
    sign = np.where(a_star >= 0, -1.0, 1.0)

    Z1 = 1.0 + (1.0 - a**2) ** (1.0 / 3.0) * (
        (1.0 + a) ** (1.0 / 3.0) + (1.0 - a) ** (1.0 / 3.0)
    )
    Z2 = np.sqrt(3.0 * a**2 + Z1**2)
    return 3.0 + Z2 + sign * np.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2))


def epsilon_from_spin(a_star):
    """
    Radiative efficiency epsilon(a*) = 1 - E_isco, the specific binding
    energy per unit rest-mass energy at the ISCO (Novikov & Thorne 1973):

        E_isco = (1 - 2/(3*r_isco))^(1/2)

    with r_isco in units of r_g (isco_radius above). This is the standard
    "thin, radiatively efficient disc" definition of epsilon; it is used
    here only as a mapping from spin to a *value* of epsilon (fed into the
    Eddington-rate normalisation, black_holes_growth_slimdisk.eddington_rate_std)
    -- not as a claim that accretion in the slim/super-Eddington regime is
    itself thin-disc-like (see MODELS.md).
    """
    r_isco = isco_radius(a_star)
    E_isco = np.sqrt(1.0 - 2.0 / (3.0 * r_isco))
    return 1.0 - E_isco
