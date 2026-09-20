from pathlib import Path
from dataclasses import dataclass, asdict
import yaml


@dataclass
class PaperReservoirParams:
    """
    Fiducial parameters for the paper-specific reservoir model
    (paper_reservoir.py), read from paper_reservoir_params.yaml at the
    project root. Centralises numbers that were previously hardcoded
    independently as module-level *_FIDUCIAL constants and
    function-signature defaults in paper_reservoir.py -- one such
    default (eta_acc) drifted out of sync with the paper's own stated
    fiducial for a period before being caught (see
    docs/2026_paper_session_code_catalogue.md); this file exists to
    prevent that class of bug recurring.

    Every run_reservoir_paper/critical_seed_paper caller can still
    override any of these per call via ordinary keyword arguments (e.g.
    for the sensitivity scans in the 2026 plots list) -- this file only
    fixes what "fiducial" means when no override is given.
    """
    # Cosmology / halo structure
    f_b: float = 0.156               # cosmic baryon fraction
    delta_vir: float = 200.0         # overdensity (relative to rho_crit) defining R_vir
    omega_m: float = 0.308
    omega_lambda: float = 0.692

    # Cold/hot-mode and angular-momentum selection (Eq. coldhot, low_spin_available_fraction)
    M_hot: float = 4.0e11            # Msun; cold/hot-mode transition mass
    phi_coldhot: float = 4.0         # dimensionless step sharpness
    lam_median: float = 0.035        # Bullock et al. (2001) median halo spin parameter
    sigma_lnj: float = 0.5           # log-normal spread of the spin-selected available fraction

    # Star formation
    epsilon_sf: float = 0.015        # nuclear star-formation efficiency (Menon, Balu & Power 2026 Table 1 fiducial)

    # Black hole growth (hobbs_slimdisk physics)
    epsilon: float = 0.1             # radiative efficiency, paper's own free-parameter fiducial
    chi_crit: float = 10.0           # graded super-Eddington cap, general trajectory integration
    strict_eddington_chi_crit: float = 1.0e12  # critical_seed_paper's boundary-finding default (no cap ever engages)
    epsilon_f: float = 5.0e-4        # King (2003) energy-driven AGN wind coupling
    eta_acc: float = 0.005           # nuclear accretion efficiency
    R_nuc_pc: float = 100.0          # fixed nuclear radius for the free-fall estimate, pc
    compaction_boost: float = 1.0    # fixed Phi_hat multiplier (>=1)

    # Critical-seed boundary target
    f_bh: float = 0.5                # overmassiveness target, M_BH(z_anchor) = f_bh*M_star(z_anchor)

    # Stellar-wind (supernova) outflow. `stellar_winds` is an option of the old run_reservoir_paper only and defaults to
    # False so that model behaves as at the last commit; the MVM always includes the stellar outflow.
    stellar_winds: bool = False      # include Ashvini's SN mass-loading outflow in run_reservoir_paper (metal-poor limit)
    eta_sn_scale: float = 1.0        # multiplier on supernovae_feedback.mass_loading_factor (used by the MVM)

    # Fields below are used by ashvini.reservoir_stock (the MVM) as follows: n_rd, c_nfw, f_mom, eta_sn_scale.
    # epsilon_sf_ext, wind_mode, wind_driver, R_nuc_rd, dm_nuclear, sf_ext_timescale, agn_wind_type, agn_wind_sink
    # exist only for the frozen pre-MVM reference (ashvini.reservoir_stock_premvm); the MVM has no such options.
    # Stock reservoir: galaxy-wide star formation and wind modes
    epsilon_sf_ext: float = 0.015    # galaxy-scale (outside R_nuc) star-formation efficiency
    wind_mode: str = "extranuclear"  # 'none' | 'extranuclear' | 'galaxy' (whole galaxy, in proportion to gas mass)
    wind_driver: str = "total"       # 'total' | 'nuclear': the SFR that drives the stellar wind
    R_nuc_rd: float | None = None    # if set, the nuclear scale is R_nuc_rd * R_d(t) (R_d = lam R_vir/sqrt2) instead of R_nuc_pc
    n_rd: float = 10.0               # galaxy scale in units of the disc scale length R_d = lam R_vir / sqrt(2); inf = no galaxy cut
    c_nfw: float = 4.0               # NFW concentration of the halo (assumed; not calibrated at z > 5)
    dm_nuclear: bool = False         # include the dark matter inside R_nuc in M_enc(R_nuc) (always included at the galaxy scale)
    sf_ext_timescale: str = "freefall"  # galaxy-scale star formation: 'freefall' (t_ff(R_gal)) | 'hubble' (0.141 t_H(z))
    agn_wind_type: str = "momentum"    # 'energy' (King 2003, epsilon_f) | 'momentum' (Mdot_out v_out = f_mom L/c, v_out = sigma)
    f_mom: float = 1.0               # momentum coupling of the momentum-driven wind (1 = L/c)
    agn_wind_sink: str = "galaxy"    # 'galaxy' | 'nuclear': gas the King AGN wind is drawn from (galaxy: in proportion to gas mass)

    # Reservoir initial conditions
    gas_mass0: float = 1.0e5
    stars_mass0: float = 1.0e3


def load_paper_reservoir_params(config_file=None) -> PaperReservoirParams:
    """
    config_file : str or Path, optional
        Defaults to the project root's paper_reservoir_params.yaml (the
        normal case, used by the module-level PAPER_PARAMS below).
        Overridable for testing config-parsing logic against a
        temporary file without touching the real one.
    """
    if config_file is None:
        root = Path(__file__).resolve().parents[1]
        config_file = root / "paper_reservoir_params.yaml"
    raw = yaml.safe_load(Path(config_file).read_text()) or {}
    return PaperReservoirParams(**raw)


def print_config(params: PaperReservoirParams):
    print("\nLoaded paper_reservoir.py fiducial parameters:\n")
    for k, v in asdict(params).items():
        print(f"  {k:<24} : {v}")


PAPER_PARAMS = load_paper_reservoir_params()
