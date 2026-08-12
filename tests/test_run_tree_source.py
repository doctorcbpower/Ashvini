"""
End-to-end test of main.py's run() with basics.tree_source='pymctrees' --
the actual feature this adds: generating a forest live via pymctrees
instead of reading a pre-built HDF5 file, then running the full baryonic
model on it and writing output, with no separate script step in between.

PARAMS is a module-level global built once at import time (see main.py's
module-level UV_background/t_d/sn_type/... derived from it), so this
monkeypatches PARAMS.io directly (a plain, unfrozen dataclass instance)
rather than trying to reload the whole module against a temporary
run_params.yaml -- run()'s only read of PARAMS.io happens before
run_forest() is called, and run_forest()'s own module-level globals are
all unrelated to io, so this is a safe, minimal patch surface for testing
the tree-loading branch dispatch specifically (the rest of the pipeline is
already covered by test_run1.py).
"""

import numpy as np
import pytest

pymctrees = pytest.importorskip("pymctrees")

import h5py

from ashvini import main
from ashvini.run_params import PymctreesSourceParams


@pytest.fixture
def pymctrees_config(tmp_path):
    path = tmp_path / "planck_like.yml"
    path.write_text(
        "Run:\n"
        "  mode: camb\n"
        "  pk_kmin: 1.0e-4\n"
        "  pk_kmax: 10.0\n"
        "  pk_npoints: 500\n"
        "Cosmology:\n"
        "  H0: 67.66\n"
        "  OmegaBar: 0.048\n"
        "  OmegaM: 0.3111\n"
        "  OmegaK: 0.0\n"
        "  As: 2.1e-9\n"
        "  ns: 0.9665\n"
        "  tau_reio: 0.0561\n"
        "  mnu: 0.0\n"
        "camb:\n"
    )
    return path


def test_run_generates_trees_live_via_pymctrees(pymctrees_config, tmp_path, monkeypatch):
    dir_out = str(tmp_path / "outputs") + "/"

    monkeypatch.setattr(main.PARAMS.io, "mass_bin", 1e10)
    monkeypatch.setattr(main.PARAMS.io, "dir_out", dir_out)
    monkeypatch.setattr(main.PARAMS.io, "tree_source", "pymctrees")
    monkeypatch.setattr(main.PARAMS.io, "tree_file", None)
    monkeypatch.setattr(main.PARAMS.io, "pymctrees", PymctreesSourceParams(
        config=str(pymctrees_config), n_halos=5, z0=0.0, z_max=3.0, dz=0.5,
        m_res=1e9, backend="numpy", seed=11,
    ))

    main.run()

    output_file = f"{dir_out}mass_bin_{main.PARAMS.io.mass_bin}_{main.sn_type}.hdf5"
    with h5py.File(output_file, "r") as f:
        grp = f[f"mass_bin_{main.PARAMS.io.mass_bin}"]
        assert grp["gas_mass"].shape[0] == 5  # n_halos
        for key in ["gas_mass", "stars_mass", "gas_metals", "stars_metals",
                    "dust_mass", "bh_mass", "sfr"]:
            arr = grp[key][:]
            assert np.all(np.isfinite(arr))
            assert np.all(arr >= 0)


def test_run_missing_pymctrees_block_raises_clear_error(monkeypatch):
    monkeypatch.setattr(main.PARAMS.io, "tree_source", "pymctrees")
    monkeypatch.setattr(main.PARAMS.io, "pymctrees", None)

    with pytest.raises(ValueError, match="pymctrees"):
        main.run()


def test_run_unknown_tree_source_raises_clear_error(monkeypatch):
    monkeypatch.setattr(main.PARAMS.io, "tree_source", "not_a_real_source")

    with pytest.raises(ValueError, match="not_a_real_source"):
        main.run()
