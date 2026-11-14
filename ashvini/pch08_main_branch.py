"""Main-progenitor histories from the Parkinson, Cole & Helly (2008) algorithm with PCH08's own adaptive step.

Why this module exists
----------------------
The Stochastic_Occupation paper uses PCH08 trees from foraois, with the step in redshift chosen branch by branch so that
the expected number of splits per step, N_upper, is at most 0.1 (PCH08 section 2.1).  The black hole paper needs only the
main progenitor of a halo of mass M0 at z = 5, followed back to the seeding redshift.  The PCH08 step is chosen per
branch, so secondary branches never influence the main branch and can be discarded without changing its statistics.  That
makes M0 up to 3e13 Msun with M_res = 1e4 Msun/h feasible (seconds to about half a minute per tree) where the full tree
would not be.

What is returned
----------------
A (n_trees, n_nodes) array of halo mass in Msun on a uniform redshift grid, ordered chronologically (highest z first), and
the grid itself.  This is the layout of `pymctrees_adapter.build_forest_for_bin`, so it can be passed straight to the
production script.  The mass between steps is interpolated linearly in ln M against z (the convention of
stochocc.trees); the adaptive steps are far finer than the output grid, so this is negligible.  Grid points earlier than
the formation of the main branch (the branch falls below M_res) are set to zero, as in the adapter.

Each tree is seeded with seed0 + tree index, so output is reproducible and independent of the number of threads.
"""
import numpy as np
import numba as nb


def make_tables(config_path, z_max, m_res_msun):
    """Foraois tables for the compiled PCH08 kernels (mirrors stochocc.trees.PCHTables.from_config).

    Power spectrum and sigma(M) are those at z = 0; the barrier delta_c(z) carries the redshift dependence, which is
    foraois' convention for trees anchored at any z0.  Returns (tables, h)."""
    from foraois import CosmoData, PCHMergerTree
    from foraois.utils import io
    params = io.get_params(str(config_path))
    h = float(params["Cosmology"]["h"])
    cd = CosmoData(params, redshift=[0.0])
    tg = PCHMergerTree(cd, params)
    tg._ensure_delta_col_covers(z_max + 5.0)
    cd.check_M_res(m_res_msun * h, lookup_factor=0.5, compiled=True)
    # alpha = |d ln sigma / d ln M| must be positive for the kernel; the stored derivative is signed (negative).
    tables = dict(dc_z=tg._dc_z_grid, dc_dc=tg._dc_dc_grid, logmass=cd._logmass, sigma=cd._sigma,
                  alpha=np.abs(cd._dlogsigma_dlogmass), j_u=tg._j_u_grid, j_val=tg._j_values,
                  G0=float(tg.G0), gamma1=float(tg.gamma1), gamma2=float(tg.gamma2))
    return tables, h


def _kernels():
    # Imported lazily so that importing this module does not require foraois.
    from foraois.pch_trees import _draw_progenitors_scalar, _interp_sorted, _pick_adaptive_step_scalar

    @nb.njit
    def walk(seed, M0, M_res, z0, zgrid, target_nupper, dz_max, min_dz, eps1,
             dc_z, dc_dc, logmass, sigma, alpha, j_u, j_val, G0, g1, g2, out):
        """Follow the main progenitor from (z0, M0) up the grid zgrid (ascending, zgrid[0] = z0).

        out[k] receives M at zgrid[k] (Msun/h), left at 0 once the branch has fallen below M_res.
        Returns (steps, splits, steps_at_min_dz)."""
        np.random.seed(seed)
        zmax = zgrid[zgrid.shape[0] - 1]
        z = z0
        m = M0
        d = _interp_sorted(z0, dc_z, dc_dc)
        out[0] = M0
        k = 1
        nstep = 0
        nsplit = 0
        nmin = 0
        while m >= M_res and z < zmax - 1e-12:
            ceiling = min(z + dz_max, zmax)
            z_next, d_next = _pick_adaptive_step_scalar(m, M_res, z, d, ceiling, target_nupper, min_dz, eps1,
                                                        dc_z, dc_dc, logmass, sigma, alpha, G0, g1, g2)
            n, a, c = _draw_progenitors_scalar(m, M_res, d, d_next - d, logmass, sigma, alpha, j_u, j_val, G0, g1, g2)
            nstep += 1
            if z_next - z <= 1.001 * min_dz and z_next < ceiling - 1e-12:
                nmin += 1
            if n == 0:
                break  # no resolved progenitor: the branch ends here
            if n >= 2:
                nsplit += 1
                m_new = max(a, c)
            else:
                m_new = a
            if m_new < M_res:
                break
            # write every grid redshift passed in this step, interpolating ln M linearly in z
            while k < zgrid.shape[0] and zgrid[k] <= z_next + 1e-12:
                f = (zgrid[k] - z) / (z_next - z)
                out[k] = np.exp((1.0 - f) * np.log(m) + f * np.log(m_new))
                k += 1
            m = m_new
            z = z_next
            d = d_next
        return nstep, nsplit, nmin

    @nb.njit(parallel=True)
    def forest(seed0, M0, n_trees, M_res, z0, zgrid, target_nupper, dz_max, min_dz, eps1,
               dc_z, dc_dc, logmass, sigma, alpha, j_u, j_val, G0, g1, g2):
        out = np.zeros((n_trees, zgrid.shape[0]))
        stats = np.zeros((n_trees, 3), dtype=np.int64)
        for i in nb.prange(n_trees):
            row = np.zeros(zgrid.shape[0])
            s = walk(seed0 + i, M0, M_res, z0, zgrid, target_nupper, dz_max, min_dz, eps1,
                     dc_z, dc_dc, logmass, sigma, alpha, j_u, j_val, G0, g1, g2, row)
            out[i, :] = row
            stats[i, 0] = s[0]
            stats[i, 1] = s[1]
            stats[i, 2] = s[2]
        return out, stats

    return forest


def build_main_branches(config_path, M0_msun, n_trees, z0, z_max, m_res_msun, dz_grid, seed0,
                        target_nupper=0.1, dz_max=0.05, min_dz=1e-6, eps1=0.1):
    """Return (halo_masses [Msun], redshifts, stats), chronological on a uniform grid of spacing dz_grid.

    stats is (n_trees, 3): steps taken, resolved splits, and steps that hit the floor min_dz (such steps may exceed the
    target N_upper, so that count should be zero)."""
    n_nodes = int(round((z_max - z0) / dz_grid)) + 1
    zgrid = z0 + dz_grid * np.arange(n_nodes)               # ascending, zgrid[0] = z0
    tab, h = make_tables(config_path, z_max, m_res_msun)
    forest = _kernels()
    out, stats = forest(int(seed0), float(M0_msun * h), int(n_trees), float(m_res_msun * h), float(z0), zgrid,
                        float(target_nupper), float(dz_max), float(min_dz), float(eps1),
                        tab["dc_z"], tab["dc_dc"], tab["logmass"], tab["sigma"], tab["alpha"], tab["j_u"], tab["j_val"],
                        tab["G0"], tab["gamma1"], tab["gamma2"])
    halo_masses = (out / h)[:, ::-1]                        # chronological, Msun
    redshifts = zgrid[::-1].copy()
    return halo_masses, redshifts, stats
