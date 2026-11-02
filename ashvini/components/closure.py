from .. import supernovae_feedback as sn


class MomentumLoading:
    """
    Mass loading of the supernova wind: gas blown out per unit mass of stars formed, from momentum
    conservation, eta_w(z, M_halo, Z_star) (supernovae_feedback.mass_loading_factor; the stellar
    metallicity enters through the supernova energy-per-mass function).
    """

    def __init__(self, cfg):
        self.cfg = cfg

    def mass_loading(self, redshift, halo_mass, stellar_metallicity):
        return sn.mass_loading_factor(redshift, halo_mass, stellar_metallicity, self.cfg)
