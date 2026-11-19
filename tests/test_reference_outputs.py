"""Milestone M0: pin the outputs of the general chain and the frozen MVM (see tests/reference_cases.py)."""
import json
import os
import warnings

import numpy as np
import pytest

from tests import reference_cases as C

warnings.filterwarnings("ignore")


def _reference(name):
    with open(os.path.join(C.GOLDEN_DIR, name)) as f:
        return json.load(f)


GENERAL = _reference("general_chain.json")


def _compare(actual, expected, rtol=1e-9):
    for key, exp in expected.items():
        got = actual[key]
        got_final = [np.nan if v is None else v for v in got["final"]]
        exp_final = [np.nan if v is None else v for v in exp["final"]]
        np.testing.assert_allclose(got_final, exp_final, rtol=rtol, atol=1e-12, equal_nan=True, err_msg=f"{key} final")
        if exp["sum"] != 0.0:
            assert got["sum"] == pytest.approx(exp["sum"], rel=rtol), f"{key} sum"


@pytest.mark.parametrize("name", list(C.CASES))
def test_general_chain_matches_reference(name, monkeypatch):
    C.apply_case(monkeypatch, C.CASES[name])
    actual = C.run_case(C.CASES[name])
    for mass_bin, expected in GENERAL[name].items():
        _compare(actual[mass_bin], expected)


def test_scalar_reference_matches_reference():
    _compare(C.run_scalar_reference(), GENERAL["scalar-reference"], rtol=1e-8)


def test_mvm_source_is_the_frozen_file():
    # docs/PRODUCTION_PROVENANCE.md records this hash; the MVM must not change during the refactor.
    assert C.mvm_sha256() == C.MVM_SHA256


def test_mvm_outputs_match_reference():
    _compare(C.run_mvm(), _reference("mvm.json"))
