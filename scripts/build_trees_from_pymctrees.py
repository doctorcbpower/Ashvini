#!/usr/bin/env python3
"""
Build a pymctrees merger forest and adapt it into the HDF5 layout
ashvini.utils.read_trees expects:

    /redshifts                          (S,)      float64   -- shared redshift grid
    /<mass_bin_group>/halo_masses       (N, S)    float64   -- Msun
    /<mass_bin_group>/halo_growth_rates (N, S)    float64   -- Msun/Gyr
where mass_bin_group is e.g. "01e10" for a bin of 1e10 Msun.

This is the "adapter" approach from the code audit (Section 6): it requires
zero changes to Ashvini's main.py, and lets pymctrees-generated trees be
validated against the existing reference tree file before considering the
heavier "direct integration" path (calling pymctrees straight from run(),
which would make pymctrees -- and transitively CLASS/CAMB -- a hard
dependency of Ashvini).

Three conversions happen here that are easy to get quietly wrong, so they're
each called out explicitly in the code below:

1. pymctrees' M0/mass_history are in Msun/h (it works in Mpc/h internally);
   Ashvini's are in plain Msun. Converted via the same 'h' pymctrees derived
   from the run's own H0.
2. pymctrees grows trees *backward* in time -- z_steps runs from z0 (low
   redshift, e.g. the present day) up to z_max (high redshift, the past) --
   whereas Ashvini's run1() integrates *forward* in cosmic time. The mass
   and redshift arrays are reversed here so index 0 is the earliest
   (highest-z) point.
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

M0_array is built from --mass-bins interpreted as the halo mass (Msun) at
z0 -- matching the existing reference tree file's convention (100 haloes
per mass bin, linearly distributed within 1e6 <= Mh/Msun <= 1e11 at z=5,
per Ashvini's README).

Usage
-----
    python scripts/build_trees_from_pymctrees.py \\
        /path/to/pymctrees/config/planck2018.yml \\
        data/inputs/merger_trees_pymctrees.h5 \\
        --mass-bins 1e6 1e7 1e8 1e9 1e10 1e11 \\
        --n-halos 100 --z0 5.0 --z-max 20.0 --dz 0.1 \\
        --backend numpy

Requires pymctrees installed separately (not an Ashvini dependency):
    pip install -e /path/to/pymctrees[camb]   # or [class]
"""

import argparse
import sys

import h5py
import numpy as np

from ashvini.utils import time_at_z


def _mass_bin_group_name(mass_bin):
    """Matches ashvini.utils.read_trees' group-naming convention exactly."""
    mass_bin = float(mass_bin)
    exponent = int(np.log10(mass_bin))
    mantissa = int(mass_bin / 10**exponent)
    return f"{mantissa:02d}e{exponent:02d}"


def _import_pymctrees():
    try:
        from pymctrees import cosmo_utils, PCHMergerTree
        from pymctrees.utils import io as pymctrees_io
    except ImportError as exc:
        raise ImportError(
            "This script requires pymctrees, which is not an Ashvini "
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

    return halo_masses, redshifts


def compute_growth_rates(halo_masses, redshifts):
    """
    (3) halo_growth_rates via differencing against Ashvini's own
    cosmic_time(z) (astropy Planck18), not pymctrees' H(z). rate[:, -1] is
    never read by run1() (its loop only ever indexes rate[:, j-1] for
    j in [1, n-1]), so it's just filled by repeating the last real value.
    """
    cosmic_time = time_at_z(redshifts)  # Gyr, same ordering as halo_masses
    dt = np.diff(cosmic_time)  # (n_steps,)
    dm = np.diff(halo_masses, axis=1)  # (N, n_steps)
    rates = dm / dt[None, :]
    rates = np.concatenate([rates, rates[:, -1:]], axis=1)  # pad to (N, n_steps+1)
    return rates


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pymctrees_config", help="Path to a pymctrees YAML config (e.g. config/planck2018.yml)")
    parser.add_argument("output", help="Path to write the adapted merger-tree HDF5 file to")
    parser.add_argument("--mass-bins", nargs="+", type=float, required=True,
                         help="Halo masses (Msun) at z0 to build trees for, e.g. 1e6 1e7 1e10")
    parser.add_argument("--n-halos", type=int, default=100,
                         help="Haloes per mass bin (default: 100, matching the reference tree file)")
    parser.add_argument("--z0", type=float, default=5.0,
                         help="Starting (present-day, for tree-building purposes) redshift (default: 5.0)")
    parser.add_argument("--z-max", type=float, default=20.0,
                         help="Maximum (earliest) redshift to grow trees to (default: 20.0)")
    parser.add_argument("--dz", type=float, default=0.1,
                         help="pymctrees redshift step size (default: 0.1)")
    parser.add_argument("--m-res", type=float, default=None,
                         help="Mass resolution (Msun); default: 1e-3 x the smallest --mass-bins value")
    parser.add_argument("--backend", default="numpy", choices=["numpy", "numba"],
                         help="pymctrees tree-building backend (default: numpy, no extra deps)")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for reproducibility")
    parser.add_argument("--compression-level", type=int, default=4, help="gzip compression level 0-9")
    args = parser.parse_args()

    cosmo_utils, PCHMergerTree, pymctrees_io = _import_pymctrees()

    run_params = pymctrees_io.get_params(args.pymctrees_config)
    h = run_params["Cosmology"]["h"]

    m_res = args.m_res if args.m_res is not None else 1e-3 * min(args.mass_bins)

    cosmo_data = cosmo_utils.CosmoData(run_params, redshift=[args.z0])
    tree_generator = PCHMergerTree(cosmo_data, run_params)

    print(f"Cosmology: H0={run_params['Cosmology']['H0']}, h={h:.4f}")
    print(f"Building trees: backend={args.backend}, z0={args.z0}, z_max={args.z_max}, "
          f"dz={args.dz}, M_res={m_res:.3e} Msun, {args.n_halos} haloes/bin")

    redshifts_ref = None
    with h5py.File(args.output, "w") as fout:
        for mass_bin in args.mass_bins:
            group_name = _mass_bin_group_name(mass_bin)
            print(f"  {group_name} (M0={mass_bin:.3e} Msun) ...")

            halo_masses, redshifts = build_forest_for_bin(
                tree_generator, mass_bin, h, args.n_halos,
                args.z0, args.z_max, m_res, args.dz, args.backend, args.seed,
            )

            if redshifts_ref is None:
                redshifts_ref = redshifts
                fout.create_dataset(
                    "redshifts", data=redshifts_ref, compression="gzip",
                    compression_opts=args.compression_level, chunks=True,
                )
            elif not np.allclose(redshifts, redshifts_ref):
                raise RuntimeError(
                    f"Redshift grid for {group_name} doesn't match the first mass "
                    "bin's -- all bins in one output file must share a grid "
                    "(same z0/z_max/dz, matching ashvini.utils.read_trees' assumption "
                    "of a single top-level 'redshifts' dataset)."
                )

            halo_growth_rates = compute_growth_rates(halo_masses, redshifts)

            grp = fout.create_group(group_name)
            grp.create_dataset(
                "halo_masses", data=halo_masses, compression="gzip",
                compression_opts=args.compression_level, chunks=True,
            )
            grp.create_dataset(
                "halo_growth_rates", data=halo_growth_rates, compression="gzip",
                compression_opts=args.compression_level, chunks=True,
            )

            # index -1 is z0/M0 by construction (always "alive"); the
            # meaningful diagnostic is whether the progenitor was still
            # resolved (M > 0) all the way back to the earliest step, z_max
            n_resolved_to_zmax = int(np.sum(halo_masses[:, 0] > 0))
            print(f"    {n_resolved_to_zmax}/{args.n_halos} haloes have a resolved "
                  f"progenitor back to z={redshifts[0]:.2f}, {len(redshifts)} steps")

    print(f"\nSaved {args.output}")


if __name__ == "__main__":
    sys.exit(main())
