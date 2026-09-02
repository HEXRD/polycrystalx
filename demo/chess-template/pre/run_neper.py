#!/usr/bin/env python
"""Run neper to mesh a tessellation file at a given characteristic length.

To Do
-----
- Use the EXTENTS to make the "-domain" argument
"""
from set_path import DOTDOT
import sys
import argparse
import subprocess

import numpy as np

DATA_DIR = DOTDOT / "data"
SEEDS = DATA_DIR / "seeds.txt"
TESS =  DATA_DIR / "vor.tess"

# Neper numbers grains from 1.
seeds_a = np.loadtxt(SEEDS)
NUMGRAINS = len(seeds_a)


def tesselate(args):

    print("running tesselation")
    cmd = [
        "neper", "-T",
        "-domain", "cube(0.7, 0.5, 0.7):translate(-0.35, -0.25, -0.35)",
        "-morphooptiini",  f"coo:file({SEEDS})",
        "-morpho", "voronoi",
        "-o", TESS,
        "-n", f"{NUMGRAINS}"
    ]
    print("running command:\n", cmd)
    result = subprocess.run(cmd)
    return result.returncode


def mesh(args):

    print("running mesh generation")

    tess_path = TESS
    cl = args.characteristic_length
    clength = cl * 1e-3

    basename = tess_path.stem
    outname = tess_path.with_name(f"{basename}-{cl:03d}.msh")

    cmd = [
        "neper", "-M",
        "-n", "1",
        str(tess_path),
        "-cl", str(clength),
        "-order", "1",
        "-o", outname
    ]

    result = subprocess.run(cmd)
    return result.returncode


def parse_args():

    parser = argparse.ArgumentParser(
        description="Run neper"
    )
    parser.add_argument(
        "-T", "--tesselate", action=argparse.BooleanOptionalAction, default=False
    )

    parser.add_argument(
        "characteristic_length",
        type=int,
        help="Characteristic length input, 1-999 (actual length is this * 1e-3)",
    )
    args = parser.parse_args()

    if not 1 <= args.characteristic_length <= 999:
        parser.error("characteristic-length must be between 1 and 999 microns")

    return args


def main():
    args = parse_args()

    if args.tesselate:
        tesselate(args)
    result = mesh(args)

    sys.exit()


if __name__ == "__main__":
    main()
