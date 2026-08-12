from pathlib import Path
from dataclasses import dataclass, asdict
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


@dataclass
class ReionizationParams:
    UVB_enabled: bool
    z_reion: float
    gamma: float
    omega: float


@dataclass
class MetalsParams:
    Z_IGM: float
    Z_yield: float


@dataclass
class DustParams:
    m_swept: float
    dust_yield: float
    dust_gamma: float
    dust_alpha: float
    m_crit: float


@dataclass
class Pop3Seeding:
    enabled: bool
    z_min: float
    M_halo_min: float
    Z_gas_max: float
    M_seed: float


@dataclass
class DirectCollapseSeeding:
    enabled: bool
    z_min: float
    M_halo_min: float
    Z_gas_max: float
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
class BlackHoleParams:
    efficiency: float
    eta_agn: float
    seeding: SeedingParams


@dataclass
class Params:
    io: IOParams
    sf: StarFormationParams
    sn: SupernovaParams
    reion: ReionizationParams
    metals: MetalsParams
    dust: DustParams
    bh: BlackHoleParams


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


PARAMS = load_params()
