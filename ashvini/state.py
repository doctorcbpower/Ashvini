"""
State of the general chain (milestone M1 of docs/SUPERSET_DESIGN.md): the time series of every evolved
quantity for N haloes on the shared redshift grid, in one place. Each array has shape (N, n_steps).
The steppers write into it; `as_dict` gives the dictionary that run_forest has always returned.
"""
from dataclasses import dataclass, field, fields

import numpy as np


@dataclass
class State:
    gas_mass: np.ndarray            # Msun, gas in the galaxy
    gas_metals: np.ndarray          # Msun, metals in the gas
    stars_mass: np.ndarray          # Msun, stars formed
    stars_metals: np.ndarray        # Msun, metals locked in stars
    dust_mass: np.ndarray           # Msun, dust in the gas
    bh_mass: np.ndarray             # Msun, central black hole
    sfr: np.ndarray                 # Msun/Gyr, star formation rate at each node
    stellar_metallicity: np.ndarray  # mass fraction, stars_metals / stars_mass

    @classmethod
    def empty(cls, n_halos, n_steps):
        return cls(**{f.name: np.zeros((n_halos, n_steps)) for f in fields(cls)})

    @property
    def shape(self):
        return self.gas_mass.shape

    def as_dict(self):
        """The evolved quantities by name, as returned in the run_forest result (the stellar metallicity is internal)."""
        return {f.name: getattr(self, f.name) for f in fields(self) if f.name != "stellar_metallicity"}
