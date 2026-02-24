import os
import numpy as np
import h5py

from scipy.integrate import solve_ivp

from tqdm import tqdm

from joblib import Parallel, delayed
from .utils import tqdm_joblib

from . import utils as utils

from .star_formation import star_formation_rate
from .gas_evolve import gas_inflow_rate, update_gas_reservoir
from .metallicity import evolve_gas_metals, evolve_stars_metals
from .dust import update_dust_reservoir

from .run_params import PARAMS, print_config

UV_background = PARAMS.reion.UVB_enabled
t_d = PARAMS.sn.delay_time  # delay time for SNe feedback, in Gyr
sn_type = PARAMS.sn.type  # type of supernova feedback


tiny = 1e-15  # small number for numerical gymnastics...

method = "LSODA"

def combined_rhs(
    t, 
    y,
    gas_accretion_rate, 
    halo_mass, 
    stellar_metallicity,
    sfr_feedback, 
    sfr_metals_feedback, 
    past_stars_mass,
    feedback_type,
    ):
    
    """
    Single combined RHS for all five state variables:
        y = [gas_mass, stars_mass, gas_metals, stars_metals, dust_mass]
    """
    
    gas_mass, stars_mass, gas_metals, stars_metals, dust_mass = y

    d_gas = update_gas_reservoir(
        t, 
        gas_mass, 
        gas_accretion_rate, 
        halo_mass,
        stars_metals, 
        sfr_feedback, 
        feedback_type,
    )

    d_stars = star_formation_rate(t, gas_mass)

    d_gas_met = evolve_gas_metals(
        t, 
        gas_metals, 
        gas_mass, 
        gas_accretion_rate, 
        halo_mass,
        stars_metals, 
        sfr_metals_feedback, 
        feedback_type,
    )

    d_stars_met = evolve_stars_metals(t, gas_metals, gas_mass)

    d_dust = update_dust_reservoir(
        t, 
        dust_mass, 
        gas_mass, 
        halo_mass,
        sfr_feedback, 
        past_stars_mass, 
        stars_metals,
    )
    return [d_gas, d_stars, d_gas_met, d_stars_met, d_dust]

def run1(halo_mass, halo_mass_rate, redshift):
    cosmic_time = utils.time_at_z(redshift)  # Gyr
    tsn = cosmic_time[0] + t_d  # Supernova switch-on time

    n = len(cosmic_time)
    gas_mass = np.zeros(n)
    gas_metals = np.zeros(n)
    stars_mass = np.zeros(n)
    stars_metals = np.zeros(n)
    sfr = np.zeros(n)
    stellar_metallicity = np.zeros(n)
    dust_mass = np.zeros(n)

    gas_accretion_rate = gas_inflow_rate(
        redshift, halo_mass, halo_mass_rate, UV_background
    )

    # delay_counter = None

    for j in range(1, n):
        t_span = [cosmic_time[j - 1], cosmic_time[j]]
        
# --- Delay index: find the timestep closest to t - t_d  ---
        t_delayed = cosmic_time[j] - t_d
        idx_delayed = int(np.searchsorted(cosmic_time, t_delayed))
        idx_delayed = max(0, min(idx_delayed, j - 1))

        if cosmic_time[j] <= tsn:
            feedback_type = "no"
            sfr_feedback = 0.0
            sfr_metals_feedback = sfr[j-1] 
        else:
            feedback_type = sn_type
            sfr_feedback = sfr[idx_delayed]
            sfr_metals_feedback = sfr[idx_delayed]
   
        # Delta stellar mass at the delayed index (used by dust RHS)
        past_stars_delta = stars_mass[idx_delayed] - stars_mass[max(0, idx_delayed - 1)]
            
        y0 = [
            gas_mass[j - 1],
            stars_mass[j - 1],
            gas_metals[j - 1],
            stars_metals[j - 1],
            dust_mass[j - 1],
        ]

        sol = solve_ivp(
            combined_rhs,
            t_span,                
            y0,
            method=method,
            args=(
                gas_accretion_rate[j - 1],
                halo_mass[j - 1],
                stellar_metallicity[j - 1],
                sfr_feedback,
                sfr_metals_feedback,         
                past_stars_delta,
                feedback_type,
            ),
        )

        (
            gas_mass[j],
            stars_mass[j],
            gas_metals[j],
            stars_metals[j],
            dust_mass[j],
        ) = sol.y[:, -1]

        # Star formation rate at current time
        sfr[j] = star_formation_rate(cosmic_time[j], gas_mass[j])

        # Enforce non-negativity
        gas_mass[j] = max(gas_mass[j], 0.0)
        gas_metals[j] = max(gas_metals[j], 0.0)
        stars_mass[j] = max(stars_mass[j], 0.0)
        stars_metals[j] = max(stars_metals[j], 0.0)
        dust_mass[j] = max(dust_mass[j], 0.0)

        # No gas metals or dust if no gas
        if gas_mass[j] <= 0:
            gas_metals[j] = 0.0
            dust_mass[j] = 0.0

        # Stellar metallicity
        if stars_mass[j] > 0:
            stellar_metallicity[j] = stars_metals[j] / stars_mass[j]
        else:
            stellar_metallicity[j] = 0.0
            stars_metals[j] = 0.0

    return {
        "gas_mass": gas_mass,
        "stars_mass": stars_mass,
        "gas_metals": gas_metals,
        "stars_metals": stars_metals,
        "dust_mass": dust_mass,
        "sfr": sfr,
        "cosmic_time": cosmic_time,
        "halo_mass": halo_mass,
        "halo_mass_rate": halo_mass_rate,
        "redshift": redshift,
    }


def run():
    print_config(PARAMS)

    halo_masses, halo_mass_rates, redshifts = utils.read_trees(
        file_path=PARAMS.io.tree_file, mass_bin=PARAMS.io.mass_bin
    )

    N_halos = np.shape(halo_masses)[0]
    with tqdm_joblib(tqdm(desc="Running halos", total=N_halos)):
        results = Parallel(n_jobs=-1)(
            delayed(run1)(halo_masses[i], halo_mass_rates[i], redshifts)
            for i in range(N_halos)
        )

    # Combine all properties across halos
    keys = results[0].keys()
    combined = {key: np.stack([res[key] for res in results], axis=0) for key in keys}

    # Output file
    output_file = PARAMS.io.dir_out + f"mass_bin_{PARAMS.io.mass_bin}_{sn_type}.hdf5"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with h5py.File(output_file, "w") as f:
        # Save common 1D arrays at the root
        f.create_dataset("cosmic_time", data=results[0]["cosmic_time"])
        f.create_dataset("redshift", data=results[0]["redshift"])

        # Group for halo properties
        grp = f.create_group(f"mass_bin_{PARAMS.io.mass_bin}")
        for key, val in combined.items():
            if key not in ["cosmic_time", "redshift"]:
                grp.create_dataset(key, data=val)

    tqdm.write(f"Saved outputs to {output_file}")
