import numpy as np

from .utils import Hubble_time, z_at_time
from .config import resolve


def time_freefall(redshift):
    """Free-fall time 0.141 t_Hubble(z) in Gyr, shared by star formation and black hole growth."""
    return 0.141 * Hubble_time(redshift)


def star_formation_rate(
    t,  # cosmic time
    gas_mass,
    cfg=None,
):
    """
    Equation 2 in Menon et al 2024: SFR = eps_ff / t_ff * M_gas, eps_ff = cfg.sf.efficiency.
    """

    redshift = z_at_time(t)

    star_formation_rate = (resolve(cfg).sf.efficiency / time_freefall(redshift)) * gas_mass
    return np.asarray(star_formation_rate)

from .config import retire_module_values as _retire  # noqa: E402

_retire(__name__, {"e_ff": "config.sf.efficiency"})
