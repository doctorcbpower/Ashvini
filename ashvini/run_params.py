from pathlib import Path
from dataclasses import dataclass, asdict
import yaml


@dataclass
class IOParams:
    mass_bin: float
    tree_file: str
    dir_out: str


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


def load_params() -> Params:
    root = Path(__file__).resolve().parents[1]
    config_file = root / "run_params.yaml"
    raw = yaml.safe_load(config_file.read_text())

    params = Params(
        io=IOParams(**raw["basics"]),
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
