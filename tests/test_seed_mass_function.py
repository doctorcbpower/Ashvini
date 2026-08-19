import numpy as np
import pytest

from ashvini.seed_mass_function import fraction_requiring_boost, boosted_fractions_by_channel, DEFAULT_CHANNELS


def test_fraction_requiring_boost_at_median_is_half():
    frac = fraction_requiring_boost(150.0, mu=np.log10(150.0), sigma=0.35)
    assert frac == pytest.approx(0.5, abs=1e-8)


def test_fraction_requiring_boost_increases_with_higher_boundary():
    mu, sigma = np.log10(150.0), 0.35
    low = fraction_requiring_boost(50.0, mu, sigma)
    high = fraction_requiring_boost(5000.0, mu, sigma)
    assert low < high


def test_boosted_fractions_by_channel_dict_input():
    boundary = {1e11: 500.0, 1e12: 5000.0, 1e13: 50000.0}
    out = boosted_fractions_by_channel(boundary)
    assert set(out.keys()) == set(DEFAULT_CHANNELS.keys())
    for channel_fracs in out.values():
        assert set(channel_fracs.keys()) == set(boundary.keys())
        for frac in channel_fracs.values():
            assert 0.0 <= frac <= 1.0
    # direct-collapse seeds (mu~2e5) should need boosting far less often
    # than Pop III (mu~150) for the same boundary
    assert out["direct_collapse"][1e12] < out["pop3"][1e12]
