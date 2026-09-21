"""
Locate the foraois checkout used by the paper-figure and diagnostic scripts.

The scripts previously hard-coded a machine-specific path. They now read the environment variable
FORAOIS_ROOT, the directory of a foraois checkout that contains ``src/`` and ``config/``::

    export FORAOIS_ROOT=/path/to/foraois

Which foraois state produced the frozen production ensemble is recorded in docs/PRODUCTION_PROVENANCE.md;
this module only resolves paths and does not affect any calculation.
"""
import os
from pathlib import Path

_root = os.environ.get("FORAOIS_ROOT")
if not _root:
    raise RuntimeError(
        "Set FORAOIS_ROOT to the foraois checkout (the directory that contains src/ and config/). "
        "See docs/PRODUCTION_PROVENANCE.md for the state used for the production ensemble.")
FORAOIS_ROOT = Path(_root).expanduser().resolve()
FORAOIS_SRC = str(FORAOIS_ROOT / "src")
FORAOIS_CONFIG = str(FORAOIS_ROOT / "config" / "menon_power_2024.yml")
