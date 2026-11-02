import pytest


def pytest_addoption(parser):
    parser.addoption("--stochocc", action="store_true", default=False,
                     help="also run the tests that rerun Stochastic_Occupation live (minutes)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--stochocc"):
        return
    skip = pytest.mark.skip(reason="rerunning Stochastic_Occupation is slow: use --stochocc")
    for item in items:
        if "stochocc" in item.keywords:
            item.add_marker(skip)
