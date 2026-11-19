"""The code blocks of docs/STUDENT_GUIDE.md run: documentation that does not run is wrong."""
import os
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = re.findall(r"```python\n(.*?)```", (ROOT / "docs" / "STUDENT_GUIDE.md").read_text(), flags=re.S)


@pytest.mark.parametrize("block", BLOCKS, ids=[f"block{i}" for i in range(len(BLOCKS))])
def test_guide_code_block_runs(block, monkeypatch):
    monkeypatch.chdir(ROOT)
    namespace = {}
    # blocks after the first use names the first one defines (halo_mass, rate, z)
    exec(BLOCKS[0], namespace)
    exec(block, namespace)
