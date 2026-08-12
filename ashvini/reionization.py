import numpy as np
import astropy.units as u

from astropy.cosmology import Planck15 as cosmo

from .run_params import PARAMS

z_rei = PARAMS.reion.z_reion
gamma = PARAMS.reion.gamma
omega = PARAMS.reion.omega


H_0 = cosmo.H0  # in km / (Mpc s)
H_0 = cosmo.H0.to(u.Gyr ** (-1))  # in 1/Gyr


beta = z_rei * ((np.log(1.82 * (10**3) * np.exp(-0.63 * z_rei) - 1)) ** (-1 / gamma))
c_omega = 2 ** (omega / 3) - 1


def s(x, y):
    """
    A step function used in the expression for accretion suppression due to UV background.
    Args:
        x, y (float): Two parameters representing any physical quantity.

    Returns:
        The value of the function s(x,y) based on the argument of the parameters.
    """
    s_val = (1 + (2 ** (y / 3) - 1) * (x ** (-y))) ** (-3 / y)
    return s_val


def M_c(z):
    """
    Characteristic mass scale for reionization (Okamoto et al. 2008).
    Args:
        z (float): Parameter for redshift.

    Returns:
        Float: The characteristic mass scale at which the baryon fraction is suppressed by a factor of two,
        compared to the universal value, because of background UV.
    """
    M_val = 1.69 * (10**10) * (np.exp(-0.63 * z) / (1 + np.exp((z / beta) ** gamma)))
    return M_val


def mu_c(z, m_halo):
    mu = m_halo / M_c(z)
    return mu


def X(z, m_halo):
    M_omega = (M_c(z) / m_halo) ** (omega)
    X_val = (3 * c_omega * M_omega) / (1 + c_omega * M_omega)
    return X_val


def epsilon(z):
    # (z/beta)**gamma grows huge well within uv_suppression's z<=10 mask
    # (e.g. ~640 at z=10 for the default z_reion=7, gamma=15) -- exp of
    # that is still finite, but its square overflows float64, correctly
    # driving part2 -> 0 in the limit. np.errstate suppresses the resulting
    # (harmless, verified against a finite/non-nan result) RuntimeWarning
    # rather than leaving every real run of this function noisy.
    with np.errstate(over="ignore"):
        part2 = (
            (gamma * (z ** (gamma - 1)) / (beta**gamma))
            * (np.exp((z / beta) ** gamma))
            / ((1 + np.exp((z / beta) ** gamma)) ** 2)
        )
        epsilon = (0.63) / (1 + np.exp((z / beta) ** gamma)) + part2
    return epsilon


def uv_suppression(z_val, m_halo, mdot_halo):
    z_val = np.asarray(z_val)
    m_halo = np.asarray(m_halo)
    mdot_halo = np.asarray(mdot_halo)

    suppression = np.ones_like(z_val)

    mask = z_val <= 10

    z_masked = z_val[mask]
    mh_masked = m_halo[mask]
    mdot_masked = mdot_halo[mask]

    # A halo/step with no mass and/or no accretion rate (e.g. a
    # not-yet-formed progenitor -- see build_trees_from_pymctrees.py's
    # pre-formation zeroing) has no accretion to suppress: the correct
    # answer is 0 regardless of what this function's internals would
    # otherwise evaluate to (gas_inflow_rate multiplies this output by
    # halo_mass_dot, which is already 0 there). But m_halo=0 feeds into
    # mu_c/X/s as a divide-by-zero (M_c(z)/m_halo, and s's x**(-y) with
    # x=mu_c=0), producing inf/nan that then poisons the already-correct
    # "0 *" via 0*nan=nan -- so this is guarded here directly rather than
    # relying on the caller's multiplication to save it, the same
    # principle as the mdot_halo=0 guard below.
    resolved = (mh_masked > 0) & (mdot_masked > 0)
    mh_safe = np.where(resolved, mh_masked, 1.0)
    mdot_safe = np.where(resolved, mdot_masked, 1.0)

    mu_vals = mu_c(z_masked, mh_safe)
    x_vals = X(z_masked, mh_safe)
    eps_vals = epsilon(z_masked)
    H_vals = cosmo.H(z_masked).value

    ratio_term = 2 * eps_vals * mh_safe * x_vals * (1 + z_masked) * H_vals / mdot_safe
    suppressed = s(mu_vals, omega) * ((1 + x_vals) - ratio_term)
    suppressed = np.where(resolved, suppressed, 0.0)

    suppression[mask] = np.maximum(suppressed, 0.0)

    return suppression
