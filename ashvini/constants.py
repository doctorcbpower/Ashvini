"""
Canonical physical and astronomical constants (CGS unless noted), shared
across every growth-model implementation in this package.

Single source of truth: black_holes_growth.py and
black_holes_growth_slimdisk.py previously each hard-coded their own
independent copy of these values. They happened to agree, but nothing
enforced that -- exactly the failure mode that produced a real (if small)
numerical mismatch when this session's paper-analysis prototype (a
separate, standalone reimplementation -- see
docs/2026_paper_session_code_catalogue.md -- frozen for reproducibility in
AshviniPapers/.../analysis/prototype_package/) was cross-checked against
this package: the prototype used its own independently-rounded constants
(m_p=1.673e-24 g, c=2.998e10 cm/s, sigma_T=6.652e-25 cm^2, 1 Gyr=3.156e16 s)
that differ from the values below at the third significant figure. That
mismatch was ~0.02-0.3% -- three orders of magnitude below any effect the
paper actually reports, so it does not change any conclusion -- but it is
exactly the kind of silent drift a shared module prevents from recurring
*within* this package going forward.

Values below are unchanged from what black_holes_growth.py and
black_holes_growth_slimdisk.py already used individually before this
consolidation -- this is a pure deduplication, not a precision upgrade, so
no existing output or test result changes as a result of this refactor.
"""

G = 6.674e-8              # gravitational constant, cm^3 g^-1 s^-2
c = 3.0e10                # speed of light, cm/s
m_p = 1.67e-24             # proton mass, g
sigma_thomson = 6.65e-25  # Thomson cross section, cm^2

Gyr_s = 3.15576e16              # seconds per Gyr
Msun_g = 1.98892e33             # grams per solar mass (IAU nominal)
pc_cm = 3.0856775814913673e18   # cm per parsec
