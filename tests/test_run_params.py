"""
Tests for run_params.py's config parsing -- specifically the
tree_source/pymctrees addition (basics.tree_source: "file" (default) or
"pymctrees", with basics.pymctrees holding the live-generation parameters).

load_params() accepts an optional config_file override (see its docstring)
specifically so these can point at a temporary file instead of the real
project run_params.yaml.
"""

import textwrap

import pytest

from ashvini.run_params import PymctreesSourceParams, load_params


def _write_config(tmp_path, basics_yaml):
    """A minimal but complete run_params.yaml, with `basics` swapped in."""
    config = textwrap.dedent(f"""\
        basics:
        {textwrap.indent(basics_yaml.strip(), '  ')}

        star_formation:
          efficiency: 0.015

        supernova:
          type: "delayed"
          delay_time: 0.015
          epsilon_p: 5
          pi_fid: 1

        reionization:
          UVB_enabled: True
          z_reion: 7
          gamma: 15
          omega: 2

        metallicity:
          Z_IGM: 1.0e-3
          Z_yield: 0.06

        dust:
          m_swept: 1.0e+3
          dust_yield: 0.004
          dust_gamma: 1.3e-4
          dust_alpha: 8
          m_crit: 1.0e+5

        black_holes:
          efficiency: 0.001
          eta_agn: 0.5
          seeding:
            pop3:
              enabled: True
              z_min: 15
              M_halo_min: 1.0e+6
              Z_gas_max: 1.0e-4
              M_seed: 1.0e+2
            direct_collapse:
              enabled: True
              z_min: 10
              M_halo_min: 1.0e+7
              Z_gas_max: 1.0e-5
              M_seed: 1.0e+5
            halo_mass_threshold:
              enabled: True
              M_halo_min: 1.0e+10
              M_seed: 1.0e+3
    """)
    path = tmp_path / "run_params_test.yaml"
    path.write_text(config)
    return path


def test_default_project_config_uses_file_tree_source():
    # The real run_params.yaml shipped in the repo -- must still parse to
    # tree_source='file' with pymctrees=None, the pre-existing behaviour,
    # confirming the new fields didn't break loading the actual project
    # config.
    params = load_params()
    assert params.io.tree_source == "file"
    assert params.io.pymctrees is None
    assert params.io.tree_file is not None


def test_missing_tree_source_defaults_to_file(tmp_path):
    # A config written before this feature existed (no tree_source key at
    # all) must still default cleanly to 'file' -- backward compatibility.
    config = _write_config(tmp_path, """
        mass_bin: 1e10
        tree_file: "./data/inputs/merger_trees.h5"
        dir_out: "./data/outputs/"
    """)
    params = load_params(config)
    assert params.io.tree_source == "file"
    assert params.io.pymctrees is None


def test_pymctrees_tree_source_parses_block(tmp_path):
    config = _write_config(tmp_path, """
        mass_bin: 1e9
        dir_out: "./data/outputs/"
        tree_source: "pymctrees"
        pymctrees:
          config: "/path/to/planck2018_camb.yml"
          n_halos: 50
          z0: 5.0
          z_max: 25.0
          dz: 0.01
          m_res: 100.0
          backend: "numba"
          seed: 42
    """)
    params = load_params(config)
    assert params.io.tree_source == "pymctrees"
    assert isinstance(params.io.pymctrees, PymctreesSourceParams)
    assert params.io.pymctrees.config == "/path/to/planck2018_camb.yml"
    assert params.io.pymctrees.n_halos == 50
    assert params.io.pymctrees.z_max == 25.0
    assert params.io.pymctrees.backend == "numba"
    assert params.io.pymctrees.seed == 42


def test_mass_bin_bare_exponent_yaml_coerced_to_float(tmp_path):
    # Regression test for a real bug: PyYAML's safe_load only parses
    # bare-exponent numbers as floats if they have *both* a decimal point
    # and an explicit +/- sign on the exponent (e.g. "1.0e+10") -- plain
    # "1e10" (as run_params.yaml has always written mass_bin) silently
    # comes back as the *string* '1e10', not a float. This never crashed
    # before because utils.read_trees() defensively did float(mass_bin)
    # itself; pymctrees_adapter.build_forest_live's `1e-3 * mass_bin` was
    # the first caller that didn't, and blew up with a TypeError.
    # IOParams.__post_init__ now coerces mass_bin explicitly so every
    # caller can rely on it being a real float.
    config = _write_config(tmp_path, """
        mass_bin: 1e10
        dir_out: "./data/outputs/"
    """)
    params = load_params(config)
    assert isinstance(params.io.mass_bin, float)
    assert params.io.mass_bin == pytest.approx(1e10)


def test_pymctrees_numeric_fields_bare_exponent_yaml_coerced(tmp_path):
    config = _write_config(tmp_path, """
        mass_bin: 1e9
        dir_out: "./data/outputs/"
        tree_source: "pymctrees"
        pymctrees:
          config: "/path/to/planck2018_camb.yml"
          m_res: 1e2
          z_max: 25
    """)
    params = load_params(config)
    assert isinstance(params.io.pymctrees.m_res, float)
    assert params.io.pymctrees.m_res == pytest.approx(1e2)
    assert isinstance(params.io.pymctrees.z_max, float)
    assert params.io.pymctrees.z_max == pytest.approx(25.0)


def test_pymctrees_block_optional_fields_default(tmp_path):
    # Only 'config' is required; everything else (n_halos, z0, z_max, dz,
    # m_res, backend, seed) should fall back to PymctreesSourceParams'
    # dataclass defaults if omitted.
    config = _write_config(tmp_path, """
        mass_bin: 1e9
        dir_out: "./data/outputs/"
        tree_source: "pymctrees"
        pymctrees:
          config: "/path/to/planck2018_camb.yml"
    """)
    params = load_params(config)
    assert params.io.pymctrees.n_halos == 100
    assert params.io.pymctrees.backend == "numpy"
    assert params.io.pymctrees.m_res is None
    assert params.io.pymctrees.seed is None
