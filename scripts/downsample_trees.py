#!/usr/bin/env python3
"""
Downsample a merger-tree HDF5 file (as produced/consumed by ashvini.utils.read_trees)
to a much smaller file for development, testing, or lighter-weight production runs.

The input format is:
    /redshifts                         (S,)      float64   -- shared redshift grid
    /<mass_bin_group>/halo_masses      (N, S)    float64
    /<mass_bin_group>/halo_growth_rates(N, S)    float64
where mass_bin_group is e.g. "01e10" for 1e10.

This script reduces file size along three independent axes, any of which can be
skipped by leaving the corresponding argument at its default:

  1. Time resolution  --step-factor K   : keep every K-th redshift step (default 1 = no change)
  2. Halo count        --n-halos N       : keep only the first N haloes per mass bin (default: all)
  3. Mass bins          --mass-bins ...  : keep only the listed bins (default: all)
  4. Precision          --dtype float32  : store as float32 instead of float64 (default: float32)

It also enables gzip+chunking on write, though note that for this kind of smooth,
high-dynamic-range mass-growth data, downsampling and float32 buy far more than
compression does (see audit notes) -- gzip is included mainly because it's free
and helps a little once N is small.

Usage
-----
    # Shrink the full production file for day-to-day dev runs (~10x fewer steps)
    python scripts/downsample_trees.py \\
        data/inputs/merger_trees.h5 data/inputs/merger_trees_dev.h5 \\
        --step-factor 10

    # Build a tiny fixture for the test suite (a couple of mass bins, a few haloes,
    # heavily decimated in time)
    python scripts/downsample_trees.py \\
        data/inputs/merger_trees.h5 tests/fixtures/merger_trees_fixture.h5 \\
        --step-factor 100 --n-halos 5 --mass-bins 01e08 01e10
"""

import argparse
import sys

import h5py
import numpy as np


def downsample(
    input_path,
    output_path,
    step_factor=1,
    n_halos=None,
    mass_bins=None,
    dtype="float32",
    compression_level=4,
):
    with h5py.File(input_path, "r") as fin:
        all_bins = [k for k in fin.keys() if k != "redshifts"]
        bins_to_use = mass_bins if mass_bins else all_bins
        missing = set(bins_to_use) - set(all_bins)
        if missing:
            raise KeyError(f"Requested mass bin(s) not found in input file: {sorted(missing)}")

        z_full = fin["redshifts"][:]
        idx = np.arange(0, len(z_full), step_factor)
        if idx[-1] != len(z_full) - 1:
            idx = np.append(idx, len(z_full) - 1)  # always keep the final step
        z_out = z_full[idx].astype(dtype)

        print(f"Redshift grid: {len(z_full)} -> {len(idx)} steps (factor {step_factor})")

        with h5py.File(output_path, "w") as fout:
            fout.create_dataset(
                "redshifts", data=z_out, compression="gzip",
                compression_opts=compression_level, chunks=True,
            )

            for bin_name in bins_to_use:
                grp_in = fin[bin_name]
                masses = grp_in["halo_masses"]
                rates = grp_in["halo_growth_rates"]

                n_available = masses.shape[0]
                n_keep = min(n_halos, n_available) if n_halos else n_available

                m_out = masses[:n_keep][:, idx].astype(dtype)
                r_out = rates[:n_keep][:, idx].astype(dtype)

                grp_out = fout.create_group(bin_name)
                grp_out.create_dataset(
                    "halo_masses", data=m_out, compression="gzip",
                    compression_opts=compression_level, chunks=True,
                )
                grp_out.create_dataset(
                    "halo_growth_rates", data=r_out, compression="gzip",
                    compression_opts=compression_level, chunks=True,
                )
                print(f"  {bin_name}: {n_available} -> {n_keep} haloes, "
                      f"{masses.shape[1]} -> {len(idx)} steps")

    import os
    in_size = os.path.getsize(input_path) / 1e6
    out_size = os.path.getsize(output_path) / 1e6
    print(f"\n{input_path}: {in_size:.1f} MB")
    print(f"{output_path}: {out_size:.1f} MB  ({in_size / out_size:.1f}x smaller)")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Path to the full-resolution merger_trees.h5")
    parser.add_argument("output", help="Path to write the downsampled file to")
    parser.add_argument("--step-factor", type=int, default=1,
                         help="Keep every K-th redshift step (default: 1, no time downsampling)")
    parser.add_argument("--n-halos", type=int, default=None,
                         help="Keep only the first N haloes per mass bin (default: all)")
    parser.add_argument("--mass-bins", nargs="+", default=None,
                         help="Only keep these mass-bin group names, e.g. 01e08 01e10 (default: all)")
    parser.add_argument("--dtype", default="float32", choices=["float32", "float64"],
                         help="Output precision (default: float32)")
    parser.add_argument("--compression-level", type=int, default=4,
                         help="gzip compression level 0-9 (default: 4)")
    args = parser.parse_args()

    downsample(
        args.input, args.output,
        step_factor=args.step_factor,
        n_halos=args.n_halos,
        mass_bins=args.mass_bins,
        dtype=args.dtype,
        compression_level=args.compression_level,
    )


if __name__ == "__main__":
    sys.exit(main())
