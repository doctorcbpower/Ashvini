"""
Adapter between pymctrees (Monte Carlo dark matter halo merger trees,
https://github.com/doctorcbpower/pymctrees) and Ashvini's own
(halo_masses, halo_growth_rates, redshifts) integration contract.

pymctrees is optional, not a hard Ashvini dependency (see _import_pymctrees
-- so CLASS/CAMB, which pymctrees needs for a real power spectrum, don't
become one transitively either). Two callers share the tree-generation
logic here:

* scripts/build_trees_from_pymctrees.py -- the "offline" path: builds a
  forest once, writes it to an HDF5 file in ashvini.utils.read_trees()'
  expected layout, for later runs to read like any other input catalogue.
* run() in main.py, when run_params.yaml's basics.tree_source is
  'pymctrees' instead of 'file' -- the "live" path: generates a forest for
  the configured mass bin on every run, no intermediate file.

Both call build_forest_for_bin()/compute_growth_rates() below, so a live
run and an offline-build-then-read run with the same parameters produce
identical trees -- there's exactly one implementation of the unit/ordering
conversions, not two copies that could drift apart.

Three conversions happen in build_forest_for_bin/compute_growth_rates that
are easy to get quietly wrong, so they're each called out explicitly there:

1. pymctrees' M0/mass_history are in Msun/h (it works in Mpc/h internally);
   Ashvini's are in plain Msun. Converted via the same 'h' pymctrees derived
   from the run's own H0.
2. pymctrees grows trees *backward* in time -- z_steps runs from z0 (low
   redshift, e.g. the present day) up to z_max (high redshift, the past) --
   whereas Ashvini's run1() integrates *forward* in cosmic time. The mass
   and redshift arrays are reversed so index 0 is the earliest (highest-z)
   point.
3. pymctrees only returns masses, not growth rates. halo_growth_rates is
   derived by differencing the (now chronological) mass array against
   Ashvini's own cosmic_time(z) mapping (astropy Planck18, via
   ashvini.utils.time_at_z) -- not against pymctrees' own H(z) -- since
   that's the expansion history run1() will actually use to interpret the
   rate. If the pymctrees run's Cosmology block doesn't match Planck18,
   this is a deliberate simplifying choice, not an oversight: it keeps the
   *output* dimensionally consistent with how Ashvini will integrate it,
   independent of whichever cosmology pymctrees used to grow the tree
   topology itself.

A fourth conversion, zeroing out each halo's pre-formation prefix (pymctrees
can't resolve a progenitor below M_res), is handled inside
build_forest_for_bin -- see _mask_unresolved_prefix's docstring.
"""

import numpy as np

from .utils import time_at_z


def _import_pymctrees():
    try:
        from pymctrees import cosmo_utils, PCHMergerTree
        from pymctrees.utils import io as pymctrees_io
    except ImportError as exc:
        raise ImportError(
            "pymctrees is required for tree_source='pymctrees' (or the "
            "build_trees_from_pymctrees.py script), and is not an Ashvini "
            "dependency (see the code audit, Section 7.4, for why CLASS/"
            "CAMB are kept optional). Install it with:\n"
            "    pip install -e /path/to/pymctrees[camb]   # or [class]"
        ) from exc
    return cosmo_utils, PCHMergerTree, pymctrees_io


def build_forest_for_bin(tree_generator, M0_msun, h, n_halos, z0, z_max, m_res_msun, dz, backend, rng_seed):
    """
    Grow n_halos trees for a single mass bin and return them in Ashvini's
    chronological, plain-Msun convention: (halo_masses, redshifts), both
    shape (n_halos, n_steps+1) / (n_steps+1,), index 0 = earliest (highest z).
    """
    if rng_seed is not None:
        np.random.seed(rng_seed)

    M0_array_hunits = np.full(n_halos, M0_msun * h)  # Msun -> Msun/h
    M_res_hunits = m_res_msun * h

    if backend == "numpy":
        mass_history, _split_events, z_steps = tree_generator.build_forest_numpy(
            M0_array=M0_array_hunits, z0=z0, z_max=z_max, M_res=M_res_hunits, dz=dz
        )
    elif backend == "numba":
        mass_history, z_steps = tree_generator.build_forest_numba(
            M0_array=M0_array_hunits, z0=z0, z_max=z_max, M_res=M_res_hunits, dz=dz
        )
    else:
        raise ValueError(f"Unknown backend '{backend}'")

    # (1) prepend M0 -- mass_history starts at the end of the first step,
    # not at z0 itself
    mass_hunits = np.concatenate([M0_array_hunits[:, None], mass_history], axis=1)

    # (2) reverse: pymctrees runs z0 (present) -> z_max (past); Ashvini
    # expects earliest (highest z) -> latest
    mass_hunits = mass_hunits[:, ::-1]
    redshifts = z_steps[::-1].copy()

    # convert Msun/h -> Msun
    halo_masses = mass_hunits / h

    # (4) zero out each halo's pre-formation prefix: pymctrees can't
    # resolve a progenitor below M_res, so a halo whose random walk
    # doesn't reach a resolved split until some z < z_max is held at a
    # frozen placeholder mass for every earlier (higher-z) step -- this is
    # a real tree-resolution limit (worse the smaller M_res is relative to
    # the halo's own formation history, e.g. z_max=25/M_res~100 for a
    # 1e7 Msun M0 can leave >75% of the grid frozen), not a discretization
    # artifact fixable by a finer dz (checked directly: it doesn't shrink
    # with dz). Left as a nonzero placeholder mass, this reads to Ashvini
    # as a real halo sitting inert for a large fraction of cosmic time,
    # and its finite-difference growth rate is then exactly 0 there --
    # which reionization.uv_suppression divides by, an Ashvini-side bug
    # fixed separately, but the deeper issue is that this stretch
    # shouldn't be presented as "a halo with mass" at all. Zeroing it
    # represents "this halo does not exist yet" instead, which is the
    # physically correct reading of an unresolved progenitor.
    is_frozen_prefix = _mask_unresolved_prefix(halo_masses)
    halo_masses = np.where(is_frozen_prefix, 0.0, halo_masses)

    return halo_masses, redshifts


def _mask_unresolved_prefix(halo_masses):
    """
    For each halo (row), find the run of steps from index 0 (earliest
    time) that are exactly equal to the row's own first value -- pymctrees'
    frozen-placeholder signature for "not yet resolved" -- and return a
    boolean mask of that prefix (True = pre-formation, to be zeroed).
    Stops at the first step that differs from halo_masses[:, 0].
    """
    N, n_steps = halo_masses.shape
    first_val = halo_masses[:, :1]
    same_as_first = halo_masses == first_val
    # cumulative-AND from the left: True only while every step so far has
    # matched the first value, i.e. the contiguous frozen run starting at
    # index 0 (a later step that happens to coincidentally re-equal the
    # first value, after the walk has already moved on, must not count).
    return np.minimum.accumulate(same_as_first, axis=1)


def compute_growth_rates(halo_masses, redshifts):
    """
    (3) halo_growth_rates via differencing against Ashvini's own
    cosmic_time(z) (astropy Planck18), not pymctrees' H(z). rate[:, -1] is
    never read by run1() (its loop only ever indexes rate[:, j-1] for
    j in [1, n-1]), so it's just filled by repeating the last real value.

    The single step where halo_masses jumps from 0 (the zeroed
    pre-formation prefix -- see _mask_unresolved_prefix) to its first
    resolved value would otherwise show a spuriously huge growth rate
    (a large mass change over one step's dt, an artifact of representing
    "formation" as instantaneous rather than a rate); that step's rate is
    zeroed too, so accretion switches on smoothly from the following step.
    """
    cosmic_time = time_at_z(redshifts)  # Gyr, same ordering as halo_masses
    dt = np.diff(cosmic_time)  # (n_steps,)
    dm = np.diff(halo_masses, axis=1)  # (N, n_steps)
    rates = dm / dt[None, :]

    was_zero = halo_masses[:, :-1] == 0.0
    now_nonzero = halo_masses[:, 1:] != 0.0
    formation_step = was_zero & now_nonzero
    rates = np.where(formation_step, 0.0, rates)

    rates = np.concatenate([rates, rates[:, -1:]], axis=1)  # pad to (N, n_steps+1)
    return rates


def build_forest_live(pymctrees_config_path, mass_bin, n_halos, z0, z_max, dz,
                       m_res=None, backend="numpy", seed=None):
    """
    Generate a single mass bin's forest live via pymctrees, returned in
    Ashvini's own (halo_masses, halo_growth_rates, redshifts) contract --
    the same shapes/ordering ashvini.utils.read_trees() returns, so this is
    a drop-in replacement for run()'s tree-loading call when
    basics.tree_source is 'pymctrees' instead of 'file'.

    Builds a fresh CosmoData/PCHMergerTree from pymctrees_config_path on
    every call -- no caching, since tree generation is fast (seconds, not
    minutes, for the halo counts run() typically uses) and caching would add
    invalidation logic (config changed? seed changed?) for little benefit.
    If you want to reuse a generated forest across multiple runs, use the
    offline scripts/build_trees_from_pymctrees.py + tree_source='file'
    instead.

    Parameters
    ----------
    pymctrees_config_path : str
        Path to a pymctrees YAML config (e.g. config/planck2018_camb.yml)
        -- cosmology, dm_model, window_function_type etc. all come from
        there, not from run_params.yaml.
    mass_bin : float
        Halo mass (Msun) at z0.
    n_halos : int
    z0, z_max, dz : float
    m_res : float or None
        Mass resolution (Msun); default 1e-3 * mass_bin.
    backend : {'numpy', 'numba'}
    seed : int or None

    Returns
    -------
    halo_masses : np.ndarray, shape (n_halos, n_steps+1), Msun
    halo_growth_rates : np.ndarray, shape (n_halos, n_steps+1), Msun/Gyr
    redshifts : np.ndarray, shape (n_steps+1,)
    """
    cosmo_utils, PCHMergerTree, pymctrees_io = _import_pymctrees()

    run_params = pymctrees_io.get_params(pymctrees_config_path)
    h = run_params["Cosmology"]["h"]
    m_res_msun = m_res if m_res is not None else 1e-3 * mass_bin

    cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[z0])
    tree_generator = PCHMergerTree(cosmo_data, run_params)

    halo_masses, redshifts = build_forest_for_bin(
        tree_generator, mass_bin, h, n_halos, z0, z_max, m_res_msun, dz, backend, seed,
    )
    halo_growth_rates = compute_growth_rates(halo_masses, redshifts)

    return halo_masses, halo_growth_rates, redshifts
