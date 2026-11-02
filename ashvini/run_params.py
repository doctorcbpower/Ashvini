from pathlib import Path
from dataclasses import dataclass, asdict, field
import yaml


@dataclass
class PymctreesSourceParams:
    """
    basics.pymctrees in run_params.yaml -- only read when
    basics.tree_source is 'pymctrees'. Cosmology, dm_model,
    window_function_type etc. all live in the referenced pymctrees config
    file, not here; this only covers the tree-building call itself (see
    ashvini.pymctrees_adapter.build_forest_live, which takes exactly these
    as arguments).
    """
    config: str          # path to a pymctrees YAML config, e.g. config/planck2018_camb.yml
    n_halos: int = 100
    z0: float = 5.0
    z_max: float = 20.0
    dz: float = 0.1
    m_res: float = None  # Msun; None -> pymctrees_adapter's default (1e-3 * mass_bin)
    backend: str = "numpy"
    seed: int = None

    def __post_init__(self):
        # PyYAML's safe_load only parses bare-exponent numbers (1e10, or
        # even 1.0e10 without an explicit +/- on the exponent) as floats
        # if they also have a decimal point *and* an explicit sign on the
        # exponent (e.g. 1.0e+10) -- a YAML 1.1 spec quirk. Plain "1e10"
        # or "1e2" silently comes back as a str, which then blows up the
        # first time it's used arithmetically (e.g. 1e-3 * mass_bin) --
        # coerced explicitly here rather than trusting every future caller
        # to defensively float()/int() it themselves.
        self.n_halos = int(self.n_halos)
        self.z0 = float(self.z0)
        self.z_max = float(self.z_max)
        self.dz = float(self.dz)
        if self.m_res is not None:
            self.m_res = float(self.m_res)
        if self.seed is not None:
            self.seed = int(self.seed)


@dataclass
class IOParams:
    mass_bin: float
    dir_out: str
    tree_source: str = "file"  # "file" (default, read tree_file) or "pymctrees" (generate live)
    tree_file: str = None
    pymctrees: PymctreesSourceParams = None

    def __post_init__(self):
        # See PymctreesSourceParams.__post_init__ -- same YAML bare-exponent
        # gotcha applies to mass_bin (e.g. "mass_bin: 1e10" in
        # run_params.yaml parses as the *string* '1e10', not a float).
        self.mass_bin = float(self.mass_bin)


@dataclass
class StarFormationParams:
    efficiency: float


@dataclass
class SupernovaParams:
    type: str
    delay_time: float
    epsilon_p: float
    pi_fid: float
    # Supernova progenitor range and IMF: used only for the supernova rate in the dust destruction term
    # (the mass-loading factor above is a separate calibrated prescription).
    imf: str = "salpeter"
    m_sn_min: float = 8.0
    m_sn_max: float = 40.0
    # When the stars explode. "delta": every star explodes delay_time after it forms (the Ashvini model).
    # "imf": each progenitor explodes after its own main-sequence lifetime (stellar_lifetimes, grid metallicity
    # lifetime_Z as a mass fraction), with the progenitor masses spread over the IMF in n_progenitor_masses
    # equal-probability bins. Only used when type is "delayed".
    kernel: str = "delta"
    lifetime_Z: float = 0.002
    n_progenitor_masses: int = 400
    # How much gas the explosions remove. "ashvini": the calibrated mass loading epsilon_p, pi_fid above.
    # "momentum": each supernova removes f_p * p_SN / v_esc. "energy": each removes 2 eps_fb E_SN / v_esc^2.
    closure: str = "ashvini"
    f_p: float = 1.0
    p_SN_msun_kms: float = 3.0e5   # momentum per supernova, Msun km/s
    eps_fb: float = 0.1            # fraction of the supernova energy that unbinds gas
    E_SN_erg: float = 1.0e51


@dataclass
class ReionizationParams:
    UVB_enabled: bool
    z_reion: float
    gamma: float
    omega: float


@dataclass
class MetalsParams:
    Z_IGM_Zsun: float  # metallicity of inflowing gas, in units of Z_sun (constants.Z_SUN, a mass fraction)
    Z_yield: float     # mass of metals returned per unit mass of stars formed (mass ratio, not Z_sun)


@dataclass
class DustParams:
    m_swept: float
    dust_yield: float
    dust_alpha: float
    m_crit: float


@dataclass
class Pop3Seeding:
    enabled: bool
    z_min: float
    M_halo_min: float
    Z_gas_max_Zsun: float  # in units of Z_sun
    M_seed: float


@dataclass
class DirectCollapseSeeding:
    enabled: bool
    z_min: float
    M_halo_min: float
    Z_gas_max_Zsun: float  # in units of Z_sun
    M_seed: float


@dataclass
class HaloMassThresholdSeeding:
    enabled: bool
    M_halo_min: float
    M_seed: float


@dataclass
class SeedingParams:
    pop3: Pop3Seeding
    direct_collapse: DirectCollapseSeeding
    halo_mass_threshold: HaloMassThresholdSeeding


@dataclass
class SigmaFeedbackParams:
    """
    black_holes.sigma_feedback -- optional isothermal-sphere M-sigma
    self-regulation (King 2003, 2005; Power, Zubovas, Nayakshin & King 2011)
    replacing the constant eta_agn AGN wind coupling with one that switches
    on once M_BH crosses the self-regulation mass M_sigma(halo velocity
    dispersion). See ashvini.black_holes_growth.velocity_dispersion/
    m_sigma and ashvini.agn_feedback.coupling_switch, and MODELS.md's "AGN
    feedback" section for the equations. Off by default: enabled=False
    reproduces the prior constant-eta_agn behaviour exactly.
    """
    enabled: bool = False
    f_g: float = 0.16       # baryon fraction relative to dark matter (King 2003; PZNK11)
    kappa_es: float = None  # cm^2/g electron-scattering opacity; None -> sigma_thomson/m_p
    transition_width: float = 0.1  # dex width of the smooth M_BH/M_sigma switch

    # Hard growth cap at PZNK11's predicted M_BH ceiling (eq. 21-22):
    # M_sigma * (1 + 0.41*sigma_200/h(z)). The paper does not give an
    # explicit equation for *how* growth halts (their own words: "the most
    # unclear part of any theory") -- only this bound on the eventual
    # overshoot -- so a hard cutoff at the ceiling is a disclosed modelling
    # choice, not something read off the paper directly. See
    # black_holes_growth.growth_ceiling and MODELS.md's "AGN feedback"
    # section. Requires sigma_feedback.enabled=True (the cap without the
    # wind switch on is not a coherent self-regulation model). Off by
    # default: without it, sigma_feedback is a wind-strength switch only
    # (checked directly not to cap M_BH -- see
    # tests/test_sigma_feedback.py::test_sigma_feedback_wind_alone_does_not_cap_bh_growth_near_m_sigma).
    growth_cap_enabled: bool = False


@dataclass
class SlimDiskParams:
    """
    black_holes.slimdisk -- only read when black_holes.growth_model is
    "hobbs_slimdisk" (see ashvini.black_holes_growth_slimdisk and
    MODELS.md). Ported from the 2026 "Differential Growth" paper's
    standalone prototype; see
    docs/2026_paper_session_code_catalogue.md for provenance.
    """
    eta_acc: float = 0.005       # nuclear accretion efficiency (dimensionless)
    R_nuc_pc: float = 100.0      # fixed nuclear radius for the free-fall estimate, pc
    compaction_boost: float = 1.0  # fixed Phi_hat multiplier (>=1); see module docstring re: no time-varying/stochastic version here
    r_crit: float = 8.0          # supply/standard-Eddington ratio threshold for the slim-disc cap to engage
    epsilon_f: float = 5.0e-4    # King (2003) energy-driven AGN wind coupling
    a_star: float = 0.5          # BH spin, via ashvini.spin.epsilon_from_spin -> radiative efficiency epsilon


@dataclass
class BlackHoleParams:
    efficiency: float
    eta_agn: float
    seeding: SeedingParams
    eddington_multiplier: float = 1.0  # f_Edd: allows super-Eddington growth if > 1
    feedback_delay_time: float = 0.0   # Gyr; 0.0 = instantaneous AGN wind (default, matches prior behaviour)
    sigma_feedback: SigmaFeedbackParams = field(default_factory=SigmaFeedbackParams)
    # "pznk11_freefall" (default): the pre-existing model above (constant
    # free-fall time, hard eddington_multiplier cap, eta_agn/M-sigma
    # feedback). "hobbs_slimdisk": the 2026-paper alternative in
    # ashvini.black_holes_growth_slimdisk (enclosed-mass free-fall,
    # graded r_crit cap, King 2003 energy-driven feedback) -- see
    # MODELS.md. Selecting "hobbs_slimdisk" does not change any default
    # behaviour under "pznk11_freefall".
    growth_model: str = "pznk11_freefall"
    slimdisk: SlimDiskParams = field(default_factory=SlimDiskParams)


@dataclass
class Params:
    io: IOParams
    sf: StarFormationParams
    sn: SupernovaParams
    reion: ReionizationParams
    metals: MetalsParams
    dust: DustParams
    bh: BlackHoleParams


def _reject_legacy_keys(raw):
    """Keys whose meaning or unit changed: fail loudly rather than reinterpret an old file silently."""
    if "Z_IGM" in raw.get("metallicity", {}):
        raise ValueError("metallicity.Z_IGM is no longer accepted: it was documented in Z_sun units but used as a mass "
                         "fraction. Use metallicity.Z_IGM_Zsun (units of Z_sun; mass fraction = value * Z_SUN).")
    if "dust_gamma" in raw.get("dust", {}):
        raise ValueError("dust.dust_gamma is no longer accepted: the supernova rate in the dust destruction term is now "
                         "nu_SN * SFR with nu_SN computed from supernova.imf, m_sn_min and m_sn_max.")
    seeding = raw.get("black_holes", {}).get("seeding", {})
    for channel, d in seeding.items():
        if "Z_gas_max" in d:
            raise ValueError(f"black_holes.seeding.{channel}.Z_gas_max is no longer accepted: use Z_gas_max_Zsun "
                             "(units of Z_sun).")


def load_params(config_file=None) -> Params:
    """
    config_file : str or Path, optional
        Defaults to the project root's run_params.yaml (the normal case,
        used by the module-level PARAMS below). Overridable for testing
        config-parsing logic against a temporary file without touching the
        real one.
    """
    if config_file is None:
        root = Path(__file__).resolve().parents[1]
        config_file = root / "run_params.yaml"
    raw = yaml.safe_load(Path(config_file).read_text())

    _reject_legacy_keys(raw)

    basics = dict(raw["basics"])
    pymctrees_raw = basics.pop("pymctrees", None)

    params = Params(
        io=IOParams(
            **basics,
            pymctrees=PymctreesSourceParams(**pymctrees_raw) if pymctrees_raw else None,
        ),
        sf=StarFormationParams(**raw["star_formation"]),
        sn=SupernovaParams(**raw["supernova"]),
        reion=ReionizationParams(**raw["reionization"]),
        metals=MetalsParams(**raw["metallicity"]),
        dust=DustParams(**raw["dust"]),
        bh=BlackHoleParams(
            efficiency=raw["black_holes"]["efficiency"],
            eta_agn=raw["black_holes"]["eta_agn"],
            eddington_multiplier=raw["black_holes"].get("eddington_multiplier", 1.0),
            feedback_delay_time=raw["black_holes"].get("feedback_delay_time", 0.0),
            sigma_feedback=SigmaFeedbackParams(**raw["black_holes"].get("sigma_feedback", {})),
            growth_model=raw["black_holes"].get("growth_model", "pznk11_freefall"),
            slimdisk=SlimDiskParams(**raw["black_holes"].get("slimdisk", {})),
            seeding=SeedingParams(
                pop3=Pop3Seeding(**raw["black_holes"]["seeding"]["pop3"]),
                direct_collapse=DirectCollapseSeeding(
                    **raw["black_holes"]["seeding"]["direct_collapse"]
                ),
                halo_mass_threshold=HaloMassThresholdSeeding(
                    **raw["black_holes"]["seeding"]["halo_mass_threshold"]
                ),
            ),
        ),
    )

    return params


def print_config(params: Params):
    print("\n Loaded simulation parameters:\n")

    def section(title, d):
        print(f"\n[{title}]")
        for k, v in d.items():
            print(f"  {k:<16} : {v}")

    section("Basics", asdict(params.io))
    section("Star Formation", asdict(params.sf))
    section("Supernova Feedback", asdict(params.sn))
    section("Reionization", asdict(params.reion))
    section("Metals", asdict(params.metals))
    section("Dust", asdict(params.dust))
    section("Black Holes", asdict(params.bh))

    from .constants import Z_SUN
    from .imf import nu_SN
    print("\n[Resolved units]")
    print(f"  Z_sun (mass fraction)        : {Z_SUN}")
    print(f"  Z_IGM                        : {params.metals.Z_IGM_Zsun} Zsun = {params.metals.Z_IGM_Zsun * Z_SUN:.3e} (mass fraction)")
    s = params.bh.seeding
    print(f"  pop3 Z_gas_max               : {s.pop3.Z_gas_max_Zsun} Zsun = {s.pop3.Z_gas_max_Zsun * Z_SUN:.3e} (mass fraction)")
    print(f"  direct_collapse Z_gas_max    : {s.direct_collapse.Z_gas_max_Zsun} Zsun = {s.direct_collapse.Z_gas_max_Zsun * Z_SUN:.3e} (mass fraction)")
    nu = nu_SN(params.sn.imf, params.sn.m_sn_min, params.sn.m_sn_max)
    print(f"  nu_SN ({params.sn.imf}, {params.sn.m_sn_min}-{params.sn.m_sn_max} Msun) : {nu:.4e} per Msun formed")


PARAMS = load_params()
