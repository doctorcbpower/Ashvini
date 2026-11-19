"""
Direct tests of ashvini.reionization_physics -- the extracted, PARAMS-independent equilibrium
UV-suppression physics (Okamoto, Gnedin & Frenk 2008). Complements test_reionization.py's
regression tests of the wrapped, global-bound API in ashvini.reionization.
"""
import numpy as np
import pytest

from ashvini import reionization_physics as phys


def test_no_import_time_dependency_on_run_params():
    # Importing this module must not pull in ashvini.run_params (and therefore must not read
    # run_params.yaml) -- the entire point of the extraction. Checked in a fresh subprocess,
    # since this test file's own test session has likely already imported ashvini.run_params via
    # other fixtures/tests, which would make an in-process check meaningless.
    import pathlib
    import subprocess
    repo_root = pathlib.Path(__file__).resolve().parents[1]
    venv_python = repo_root / ".venv" / "bin" / "python"
    result = subprocess.run(
        [str(venv_python), "-c",
         "import sys; import ashvini.reionization_physics; "
         "print('ashvini.run_params' in sys.modules)"],
        capture_output=True, text=True, cwd=str(repo_root),
    )
    assert result.stdout.strip() == "False", result.stdout + result.stderr


def test_beta_matches_ashvini_fiducial():
    assert phys.beta(z_rei=7, gamma=15) == pytest.approx(6.498426358171733, rel=1e-12)


def test_c_omega_matches_ashvini_fiducial():
    assert phys.c_omega(omega=2) == pytest.approx(0.5874010519681994, rel=1e-12)


def test_M_c_matches_ashvini_regression_values():
    assert phys.M_c(0.0, z_rei=7, gamma=15) == pytest.approx(8450000000.0, rel=1e-10)
    assert phys.M_c(6.5, z_rei=7, gamma=15) == pytest.approx(75501002.8026284, rel=1e-8)


def test_s_at_mu_equals_one_is_one_half_for_any_omega():
    # Definition of M_c: baryon fraction suppressed by a factor of two at m_halo = M_c(z).
    for omega in (1.5, 2.0, 3.0, 10.0):
        assert phys.s(1.0, omega) == pytest.approx(0.5, rel=1e-12)


def test_s_bounded_in_unit_interval():
    for x in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
        for omega in (1.5, 2.0, 5.0):
            val = phys.s(x, omega)
            assert 0.0 <= val <= 1.0


def test_s_monotonic_increasing_in_x():
    omega = 2.0
    xs = np.array([0.1, 0.5, 1.0, 2.0, 10.0])
    vals = np.array([phys.s(x, omega) for x in xs])
    assert np.all(np.diff(vals) > 0)


def test_mu_c_consistent_with_M_c():
    z, m, z_rei, gamma = 6.0, 1e9, 7, 15
    assert phys.mu_c(z, m, z_rei, gamma) == pytest.approx(m / phys.M_c(z, z_rei, gamma), rel=1e-12)


def test_functions_are_pure_no_hidden_state():
    # Calling with different (z_rei, gamma, omega) side by side must not interfere -- there is no
    # module-level state to leak between calls (unlike ashvini.reionization's globals).
    a = phys.M_c(6.0, z_rei=7, gamma=15)
    b = phys.M_c(6.0, z_rei=5, gamma=15)
    c = phys.M_c(6.0, z_rei=7, gamma=15)
    assert a == c
    assert a != b
