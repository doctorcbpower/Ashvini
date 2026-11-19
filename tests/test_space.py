"""The parameter-space exploration: Saltelli design, Sobol indices (checked on the Ishigami function), and the parameter mapping."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import space_common as S  # noqa: E402


def ishigami(x, a=7.0, b=0.1):
    x = -np.pi + 2 * np.pi * x
    return np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2 + b * x[:, 2] ** 4 * np.sin(x[:, 0])


def test_sobol_indices_of_the_ishigami_function():
    """Known: S1 = 0.314, 0.442, 0; ST = 0.558, 0.442, 0.244."""
    from scipy.stats import qmc
    N, D = 2 ** 14, 3
    u = qmc.Sobol(d=2 * D, scramble=True, seed=1).random(N)
    A, B = u[:, :D], u[:, D:]
    blocks = [A, B]
    for i in range(D):
        AB = A.copy(); AB[:, i] = B[:, i]; blocks.append(AB)
    f = ishigami(np.vstack(blocks))
    S1, ST = S.sobol_indices(f, N, D)
    assert np.allclose(S1, [0.314, 0.442, 0.0], atol=0.03)
    assert np.allclose(ST, [0.558, 0.442, 0.244], atol=0.03)


def test_bootstrap_widths_shrink_with_more_samples():
    from scipy.stats import qmc
    widths = []
    for N in (2 ** 8, 2 ** 12):
        u = qmc.Sobol(d=6, scramble=True, seed=2).random(N)
        A, B = u[:, :3], u[:, 3:]
        blocks = [A, B] + [np.column_stack([B[:, i] if j == i else A[:, j] for j in range(3)]) for i in range(3)]
        widths.append(S.bootstrap_indices(ishigami(np.vstack(blocks)), N, 3, n_boot=60)[1].max())
    assert widths[1] < widths[0]


@pytest.mark.parametrize("subspace", list(S.SUBSPACES))
def test_design_shape_reproducibility_and_ranges(subspace):
    D = len(S.SUBSPACES[subspace][1])
    u = S.design(subspace)
    assert u.shape == (S.N_BASE * (D + 2), D) and np.all((u >= 0) & (u <= 1))
    assert np.array_equal(u, S.design(subspace))                                   # deterministic
    A, B = u[:S.N_BASE], u[S.N_BASE:2 * S.N_BASE]
    for i in range(D):                                                             # AB_i is A with column i from B
        AB = u[(2 + i) * S.N_BASE:(3 + i) * S.N_BASE]
        assert np.array_equal(np.delete(AB, i, axis=1), np.delete(A, i, axis=1)) and np.array_equal(AB[:, i], B[:, i])
    for row in (u[0], u[-1]):
        for name, v in S.physical(subspace, row).items():
            lo, hi, _, _ = S.PARAMETERS[name]
            assert lo <= v <= hi


def test_config_changes_are_valid_configs_and_carry_the_subspace():
    from ashvini.config import DEFAULT, override
    for subspace, (fixed, names) in S.SUBSPACES.items():
        row = np.full(len(names), 0.5)
        cfg = override(DEFAULT, **S.config_changes(subspace, S.physical(subspace, row)))
        assert cfg.reion.form == fixed["reion.form"]
        assert (cfg.sn.wind_model == "energy") == ("sn.wind_model" in fixed)
        assert cfg.disc.enabled == ("disc.enabled" in fixed)


def test_the_default_parameter_values_reproduce_the_default_configuration():
    from ashvini.config import DEFAULT, override
    defaults = {"eps": 0.015, "epsilon_p": 5.0, "delay": 0.015, "z_rei": 7.0, "gamma": 15.0, "omega": 2.0, "bh_eff": 0.001, "eta_agn": 0.5}
    cfg = override(DEFAULT, **S.config_changes("cur_ash", defaults))
    assert (cfg.sf.efficiency, cfg.sn.epsilon_p, cfg.sn.delay_time, cfg.reion.z_reion, cfg.reion.gamma, cfg.reion.omega,
            cfg.bh.efficiency, cfg.bh.eta_agn) == (DEFAULT.sf.efficiency, DEFAULT.sn.epsilon_p, DEFAULT.sn.delay_time, DEFAULT.reion.z_reion,
                                                    DEFAULT.reion.gamma, DEFAULT.reion.omega, DEFAULT.bh.efficiency, DEFAULT.bh.eta_agn)


def _fake_shards(tmp_path, sub="cur_ash"):
    """Write a complete set of shards for one sub-space whose stellar masses depend only on the first parameter (so its Sobol index is ~1)."""
    import shmr_z0_common as C
    n = S.n_samples(sub)
    u = S.design(sub)
    res = np.zeros((n, len(C.MASSES), len(S.FIELDS)))
    scale = 10 ** (-0.5 + 1.0 * u[:, 0])                                    # a factor 0.3 to 3 set by the first parameter
    for k, m in enumerate(C.MASSES):
        ref = np.sqrt(float(__import__("stellar_to_halo_mass_relation").behroozi2013_mstar(m)) * float(__import__("stellar_to_halo_mass_relation").moster2013_mstar(m)))
        res[:, k, 0] = ref * scale
    res[:, :, 8] = 0
    for shard in range(4):
        idx = np.arange(shard, n, 4)
        np.savez(tmp_path / f"{sub}_shard{shard}.npz", samples=idx, results=res[idx], done=np.ones(len(idx), dtype=bool))
    return res


def test_analysis_finds_the_driving_parameter_and_the_accepted_samples(tmp_path, monkeypatch):
    import space_analyse as A
    monkeypatch.setattr(A, "PLOT", str(tmp_path))
    res = _fake_shards(tmp_path)
    loaded = A.load_subspace("cur_ash")
    assert np.allclose(loaded, res)
    ok = A.accepted_mask(loaded)
    assert 0 < ok.sum() < len(ok)                                           # a factor 0.3-3 about the geometric mean: some inside the band, some outside
    D = len(S.SUBSPACES["cur_ash"][1])
    f = np.log10(loaded[:, 8, 0])
    S1, ST = S.sobol_indices(f, S.N_BASE, D)
    assert S1[0] > 0.95 and ST[0] > 0.95 and np.all(np.abs(S1[1:]) < 0.15)
    A.main()
    assert (tmp_path / "space_summary.txt").exists() and (tmp_path / "space_sobol_total.png").exists()
    assert "accepted:" in (tmp_path / "space_summary.txt").read_text()
