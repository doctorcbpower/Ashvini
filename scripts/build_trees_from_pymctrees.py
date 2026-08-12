#!/usr/bin/env python3
"""
Build a pymctrees merger forest and adapt it into the HDF5 layout
ashvini.utils.read_trees expects:

    /redshifts                          (S,)      float64   -- shared redshift grid
    /<mass_bin_group>/halo_masses       (N, S)    float64   -- Msun
    /<mass_bin_group>/halo_growth_rates (N, S)    float64   -- Msun/Gyr
where mass_bin_group is e.g. "01e10" for a bin of 1e10 Msun.

This is the "offline" path (write once, read like any other input
catalogue via tree_source: file in run_params.yaml): builds the forest,
writes it to a file, and never touches Ashvini's run() itself. For
generating a forest fresh on every run instead, see tree_source: pymctrees
in run_params.yaml's basics section (ashvini.pymctrees_adapter.
build_forest_live, used from main.py's run()) -- both paths share the same
tree-generation/unit-conversion logic in ashvini.pymctrees_adapter, so they
produce identical trees for the same parameters.

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

from ashvini.pymctrees_adapter import (
    _import_pymctrees,
    build_forest_for_bin,
    compute_growth_rates,
)
from ashvini.utils import mass_bin_group_name


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
            group_name = mass_bin_group_name(mass_bin)
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

            # index -1 is z0/M0 by construction (always "alive"); index 0
            # is always 0 by construction too (_mask_unresolved_prefix
            # zeroes every halo's earliest step unless it happened to
            # resolve on literally the first step). The meaningful
            # diagnostic is each halo's formation redshift -- the highest z
            # at which its mass is still nonzero -- and how much of the
            # grid that leaves as pre-formation.
            first_resolved_idx = np.argmax(halo_masses > 0, axis=1)  # per-halo
            formation_z = redshifts[first_resolved_idx]
            frac_preformation = first_resolved_idx / (len(redshifts) - 1)
            print(f"    formation z: median={np.median(formation_z):.2f} "
                  f"(range {formation_z.min():.2f}-{formation_z.max():.2f}); "
                  f"median pre-formation fraction of grid: {np.median(frac_preformation):.1%}")

    print(f"\nSaved {args.output}")


if __name__ == "__main__":
    sys.exit(main())
