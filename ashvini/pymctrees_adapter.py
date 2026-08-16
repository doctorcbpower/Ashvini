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

Four conversions happen in build_forest_for_bin/compute_growth_rates that
are easy to get quietly wrong, so they're each called out explicitly there:

1. pymctrees' M0/mass_history are in Msun/h (it works in Mpc/h internally);
   Ashvini's are in plain Msun. Converted via the same 'h' pymctrees derived
   from the run's own H0.
2. pymctrees grows trees *backward* in time -- z_steps runs from z0 (low
   redshift, e.g. the present day) up to z_max (high redshift, the past) --
   whereas Ashvini's run1() integrates *forward* in cosmic time. The mass
   and redshift arrays are reversed so index 0 is the earliest (highest-z)
   point.
3. halo_growth_rates no longer comes from differencing the (now
   chronological) mass array. pymctrees' build_forest_numpy/numba (as of
   the smooth_accretion/merger_mass addition) already decompose each
   step's mass change analytically into a smooth (unresolved-accretion)
   channel and a discrete merger-event channel; differencing the total
   mass instead conflated the two into one derivative that looked "noisy"
   partly because it genuinely was (real, discrete merger jumps mixed in
   with continuous accretion). halo_growth_rates is (smooth_accretion +
   merger_mass)/dt -- both channels, not just the smooth one. This is a
   deliberate choice, not an oversight: Ashvini only tracks the main-
   progenitor branch, with no separate model of a merging companion's own
   gas/stars, so excluding merger-driven mass growth from the baryon-
   accretion term would just make that dark matter's cosmic baryon budget
   vanish from the model entirely, not represent it more correctly.
   merger_mass is *also* returned separately, sparse/event-like (0 except
   at an actual merger step), for any future merger-triggered physics
   (e.g. a starburst on top of the baseline accretion-driven rate) to
   build on -- not yet consumed that way by anything in Ashvini today.
4. Zeroing out each halo's pre-formation prefix (pymctrees can't resolve a
   progenitor below M_res) is handled inside build_forest_for_bin -- see
   _mask_unresolved_prefix's docstring. smooth_accretion/merger_mass need
   no separate masking for this: pymctrees' kernels already leave them at
   0 for every step in that same prefix (checked directly -- see
   pymctrees' own test_smooth_accretion_and_merger_mass_conserve_mass_*
   tests), so reversing them the same way as halo_masses lines up for free.

Only the *live* tree_source='pymctrees' path (build_forest_live, used by
main.py's run()) gets the smooth/merger decomposition; the offline
HDF5-file path (scripts/build_trees_from_pymctrees.py + tree_source='file')
still stores a single halo_growth_rates dataset with no separate merger
channel -- extending that file format is a separate, not-yet-made decision.
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
    chronological, plain-Msun convention:
    (halo_masses, redshifts, smooth_accretion, merger_mass) --
    halo_masses/redshifts shape (n_halos, n_steps+1)/(n_steps+1,), index 0 =
    earliest (highest z); smooth_accretion/merger_mass shape
    (n_halos, n_steps) (one entry per *gap* between consecutive
    halo_masses columns, not per column -- see pymctrees'
    build_forest_numpy/numba docstrings for what they mean).
    """
    if rng_seed is not None:
        np.random.seed(rng_seed)

    M0_array_hunits = np.full(n_halos, M0_msun * h)  # Msun -> Msun/h
    M_res_hunits = m_res_msun * h

    if backend == "numpy":
        mass_history, _split_events, z_steps, smooth_accretion_hunits, merger_mass_hunits = tree_generator.build_forest_numpy(
            M0_array=M0_array_hunits, z0=z0, z_max=z_max, M_res=M_res_hunits, dz=dz
        )
    elif backend == "numba":
        mass_history, z_steps, smooth_accretion_hunits, merger_mass_hunits = tree_generator.build_forest_numba(
            M0_array=M0_array_hunits, z0=z0, z_max=z_max, M_res=M_res_hunits, dz=dz
        )
    else:
        raise ValueError(f"Unknown backend '{backend}'")

    # (1) prepend M0 -- mass_history starts at the end of the first step,
    # not at z0 itself
    mass_hunits = np.concatenate([M0_array_hunits[:, None], mass_history], axis=1)

    # (2) reverse: pymctrees runs z0 (present) -> z_max (past); Ashvini
    # expects earliest (highest z) -> latest. smooth_accretion/merger_mass
    # need only a reversal, no prepend: pymctrees' raw index j is exactly
    # the gap between chronological (post-reversal) indices
    # (n_steps-1-j) and (n_steps-j), so flipping column order lines them
    # up with halo_masses' n_steps gaps directly (derived and verified
    # against pymctrees' own mass-conservation identity tests).
    mass_hunits = mass_hunits[:, ::-1]
    redshifts = z_steps[::-1].copy()
    smooth_accretion_hunits = smooth_accretion_hunits[:, ::-1]
    merger_mass_hunits = merger_mass_hunits[:, ::-1]

    # convert Msun/h -> Msun
    halo_masses = mass_hunits / h
    smooth_accretion = smooth_accretion_hunits / h
    merger_mass = merger_mass_hunits / h

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

    # smooth_accretion/merger_mass at gap k depend on halo_masses[:, k]
    # (the gap's earlier/chronologically-first endpoint): fully-interior
    # frozen gaps are already 0 from pymctrees itself (its kernels skip
    # the F/merger computation entirely below M_res -- verified by
    # pymctrees' own tests), but the single *formation* gap (frozen ->
    # first resolved value) is a real, pymctrees-computed nonzero value
    # that would otherwise read as a spurious "instant formation" event --
    # zeroed here for the same reason the old compute_growth_rates zeroed
    # its analogous formation-step rate: accretion should switch on
    # smoothly from the following step, not as one artificial spike.
    gap_starts_in_prefix = is_frozen_prefix[:, :-1]
    smooth_accretion = np.where(gap_starts_in_prefix, 0.0, smooth_accretion)
    merger_mass = np.where(gap_starts_in_prefix, 0.0, merger_mass)

    return halo_masses, redshifts, smooth_accretion, merger_mass


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


def compute_growth_rates(smooth_accretion, merger_mass, redshifts):
    """
    (3) halo_growth_rates as (smooth_accretion + merger_mass)/dt, against
    Ashvini's own cosmic_time(z) (astropy Planck18), not pymctrees' H(z) --
    consistent with run1()'s own expansion history. rate[:, -1] is never
    read by run1() (its loop only ever indexes rate[:, j-1] for
    j in [1, n-1]), so it's just filled by repeating the last real value.

    Formerly this differenced the *total* halo_masses array directly
    (dm/dt) -- numerically noisy, since a discrete merger jump divided by
    one small step's dt produces an artificially huge instantaneous rate,
    all concentrated into a single step rather than reflecting the merger
    event on any physically meaningful timescale. smooth_accretion and
    merger_mass are now computed analytically inside pymctrees' tree-
    building itself (build_forest_for_bin's smooth_accretion/merger_mass,
    already the formation-step-masked, chronologically-ordered channels --
    see that function's docstring), so this just sums and divides by dt --
    no differencing or formation-step special-casing needed here. Both
    channels are included deliberately (see the module docstring's point
    3): Ashvini's main-progenitor-only tracking has nowhere else to credit
    the baryon budget that should accompany merger-driven dark matter mass
    growth, so excluding it from this rate isn't more physically correct,
    just a smaller number.
    """
    cosmic_time = time_at_z(redshifts)  # Gyr, same ordering as halo_masses
    dt = np.diff(cosmic_time)  # (n_steps,)
    rates = (smooth_accretion + merger_mass) / dt[None, :]
    rates = np.concatenate([rates, rates[:, -1:]], axis=1)  # pad to (N, n_steps+1)
    return rates


def build_forest_live(pymctrees_config_path, mass_bin, n_halos, z0, z_max, dz,
                       m_res=None, backend="numpy", seed=None):
    """
    Generate a single mass bin's forest live via pymctrees, returned in
    Ashvini's own (halo_masses, halo_growth_rates, redshifts, merger_mass)
    contract -- a drop-in replacement for run()'s tree-loading call when
    basics.tree_source is 'pymctrees' instead of 'file' (the first three
    return values match the shapes/ordering ashvini.utils.read_trees()
    returns; merger_mass is new, not part of the file-based contract --
    see the module docstring).

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
        Total growth rate (smooth accretion + merger-driven mass growth) --
        see compute_growth_rates' docstring for why merger_mass is included
        rather than excluded.
    redshifts : np.ndarray, shape (n_steps+1,)
    merger_mass : np.ndarray, shape (n_halos, n_steps), Msun
        Per-halo, per-step mass gained via a resolved merger; 0 except at
        an actual merger step. Already folded into halo_growth_rates above
        -- returned separately too as the discrete-event record, for any
        future merger-triggered physics. One entry per *gap* between
        consecutive halo_masses columns (n_steps, not n_steps+1).
    """
    cosmo_utils, PCHMergerTree, pymctrees_io = _import_pymctrees()

    run_params = pymctrees_io.get_params(pymctrees_config_path)
    h = run_params["Cosmology"]["h"]
    m_res_msun = m_res if m_res is not None else 1e-3 * mass_bin

    cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[z0])
    tree_generator = PCHMergerTree(cosmo_data, run_params)

    halo_masses, redshifts, smooth_accretion, merger_mass = build_forest_for_bin(
        tree_generator, mass_bin, h, n_halos, z0, z_max, m_res_msun, dz, backend, seed,
    )
    halo_growth_rates = compute_growth_rates(smooth_accretion, merger_mass, redshifts)

    return halo_masses, halo_growth_rates, redshifts, merger_mass
