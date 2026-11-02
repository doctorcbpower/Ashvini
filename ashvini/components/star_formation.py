from ..star_formation import time_freefall


class FreeFallStarFormation:
    """
    Star formation at a fixed fraction of the gas per free-fall time: SFR = B_sf(z) * M_gas, with
    B_sf = eps_ff / t_ff(z) and t_ff = 0.141 t_Hubble(z) (Menon et al 2024, Eq. 2). eps_ff is
    cfg.sf.efficiency. B_sf is the same for every halo at a given redshift.

    To try another law, write a class with the same `rate_coefficient(redshift)` and hand it to the
    stepper (see steppers/grid.py).
    """

    def __init__(self, cfg):
        self.efficiency = cfg.sf.efficiency

    def rate_coefficient(self, redshift):
        """B_sf in 1/Gyr."""
        return self.efficiency / time_freefall(redshift)
