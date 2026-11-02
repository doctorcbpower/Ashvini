"""
Initial mass function and the number of core-collapse supernovae per solar mass of stars formed.

Only one thing is needed from the IMF in the general chain: nu_SN, the number of stars in
[m_sn_min, m_sn_max] per Msun of stars formed, used for the supernova rate that destroys dust.
The mass-loading factor of the stellar wind (supernovae_feedback.py) is a separate calibrated
prescription and does not depend on these limits.

The normalisation integral runs over 0.1-100 Msun; m_sn_max may be smaller than 100 Msun (stars above
it are assumed not to explode as supernovae) but the IMF itself is not truncated there.
"""
import numpy as np
from scipy.integrate import quad

M_IMF_LOW = 0.1
M_IMF_HIGH = 100.0


def _salpeter(m):
    return m ** -2.35


def _kroupa(m):
    """Kroupa (2001): slope 1.3 for 0.08-0.5 Msun, 2.3 above, continuous at 0.5 Msun."""
    return m ** -1.3 if m < 0.5 else 0.5 * m ** -2.3


def _chabrier(m):
    """Chabrier (2003): lognormal (m_c = 0.079, sigma = 0.69) below 1 Msun, slope 2.3 above."""
    ln = lambda x: np.exp(-(np.log10(x) - np.log10(0.079)) ** 2 / (2 * 0.69 ** 2)) / x
    return ln(m) if m < 1.0 else ln(1.0) * m ** -2.3


IMF_DNDM = {"salpeter": _salpeter, "kroupa": _kroupa, "chabrier": _chabrier}


def nu_SN(imf="salpeter", m_sn_min=8.0, m_sn_max=40.0):
    """Supernova progenitors per Msun of stars formed (1/Msun)."""
    if imf not in IMF_DNDM:
        raise ValueError(f"Unknown IMF '{imf}' -- must be one of {sorted(IMF_DNDM)}.")
    if not (M_IMF_LOW < m_sn_min < m_sn_max <= M_IMF_HIGH):
        raise ValueError(f"Need {M_IMF_LOW} < m_sn_min < m_sn_max <= {M_IMF_HIGH}; got {m_sn_min}, {m_sn_max}.")
    f = IMF_DNDM[imf]
    breaks = [b for b in (0.5, 1.0) if M_IMF_LOW < b < M_IMF_HIGH]
    mass, _ = quad(lambda m: m * f(m), M_IMF_LOW, M_IMF_HIGH, points=breaks, limit=200)
    number, _ = quad(f, m_sn_min, m_sn_max)
    return number / mass
