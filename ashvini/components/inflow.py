from ..gas_evolve import gas_inflow_rate


class CosmicInflow:
    """
    Gas accretion rate on to the galaxy: the halo mass growth rate times the baryon fraction, reduced by
    the photoionisation (UV) suppression after reionisation (Okamoto et al 2008; cfg.reionization).
    See gas_evolve.gas_inflow_rate for the equation. The UV suppression lives inside this function, so
    replacing the inflow replaces it too.
    """

    def __init__(self, config=None):
        self.config = config  # None = config.DEFAULT

    def rate(self, redshift, halo_mass, halo_mass_rate):
        """Msun/Gyr, same shape as halo_mass."""
        return gas_inflow_rate(redshift, halo_mass, halo_mass_rate, cfg=self.config)
