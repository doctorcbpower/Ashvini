"""
The history of one halo as a function of time, for the event integrator: M_h(t), z(t), and the escape speed.
Built from the arrays a merger tree gives on its redshift grid. The mass between two nodes is interpolated
linearly in ln M against z (linearly in M if either node has no mass, e.g. before the halo forms); outside the
tabulated redshift range it is zero. Times are in Gyr, on Ashvini's clock (utils.time_at_z).
"""
import numpy as np

from . import halo, utils


class HaloHistory:
    def __init__(self, halo_mass, redshift):
        order = np.argsort(redshift)                    # ascending z for interpolation
        self.z = np.asarray(redshift, dtype=float)[order]
        self.M = np.asarray(halo_mass, dtype=float)[order]
        self.t_min = float(utils.time_at_z(self.z[-1]))
        self.t_max = float(utils.time_at_z(self.z[0]))

    def z_at_t(self, t):
        return float(utils.z_at_time(t))

    def mass_at_z(self, z):
        z_nodes, M = self.z, self.M
        if z < z_nodes[0] - 1e-12 or z > z_nodes[-1] + 1e-12:
            return 0.0
        k = min(max(np.searchsorted(z_nodes, z, side="right") - 1, 0), len(z_nodes) - 2)
        if z_nodes[k + 1] == z_nodes[k]:
            return float(M[k])
        f = (z - z_nodes[k]) / (z_nodes[k + 1] - z_nodes[k])
        if M[k] > 0 and M[k + 1] > 0:
            return float(np.exp((1 - f) * np.log(M[k]) + f * np.log(M[k + 1])))
        return float((1 - f) * M[k] + f * M[k + 1])

    def mass_at_t(self, t):
        return self.mass_at_z(self.z_at_t(t))

    def v_esc_at_t(self, t):
        z = self.z_at_t(t)
        return float(halo.v_esc_cms(self.mass_at_z(z), z))
