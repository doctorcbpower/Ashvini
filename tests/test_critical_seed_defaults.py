"""
Regression guards for the chi_crit/r_crit=1.0 "strict Eddington" bug
(2026-09-15; see docs/2026_paper_session_code_catalogue.md and the
ashvini_chi_crit_eddington_bug memory note). r_crit/chi_crit=1.0 does NOT
give strict Eddington-limited growth -- see
test_black_holes_growth_slimdisk.py::test_slimdisk_r_crit_equal_one_is_not_strict_eddington
for why -- so both modules that claim to compute a "strict Eddington
boundary" must default to a very large value instead. These tests exist
so a well-intentioned but wrong "simplification" back to 1.0 cannot land
silently in either module.
"""

import inspect

from ashvini import critical_seed
from ashvini import paper_reservoir


def test_critical_seed_default_r_crit_is_not_one():
    default = inspect.signature(critical_seed.critical_seed).parameters["r_crit"].default
    assert default == critical_seed.STRICT_EDDINGTON_R_CRIT
    assert default >= 1e6


def test_critical_seed_paper_default_chi_crit_is_not_one():
    default = inspect.signature(paper_reservoir.critical_seed_paper).parameters["chi_crit"].default
    assert default == paper_reservoir.STRICT_EDDINGTON_CHI_CRIT
    assert default >= 1e6
