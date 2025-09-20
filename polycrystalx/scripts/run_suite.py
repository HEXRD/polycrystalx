"""Run suite of Jobs"""
import sys
import argparse
import tempfile
import pathlib
import subprocess
import pickle

import numpy as np

from . import get_input_module


def main():
    """run a suite of jobs"""
    p = argparser(*sys.argv)
    args = p.parse_args()
    print("job keys: ", args.keys)

    user_module = get_input_module(args.input_module)
    if not hasattr(user_module, args.keys):
        emsg = f"job keys attribute '{args.keys}' was not found"
        raise AttributeError(emsg)

    for job in getattr(user_module, args.keys):
        # Pickle to a temp file.
        fp = tempfile.NamedTemporaryFile(delete=False)

        with open(fp.name, "wb") as f:
            pickle.dump(job, f)

        print("\n===== New Job Starting", flush=True)

        # Now run MPI job.
        cmd = [
            "mpirun", "-np", str(args.n),
            "pxx_mpijob", args.input_module, fp.name
        ]
        subprocess.run(cmd)
        pathlib.Path(fp.name).unlink()


def argparser(*args):

    p = argparse.ArgumentParser(
        description="run a suite of jobs"
    )
    p.add_argument(
        'input_module', type=str,
        help="module with inputs (name or path)"
    )
    p.add_argument(
        '-n', type=int,
        default=2,
        help="number of processes"
    )
    p.add_argument(
        '-k', '--keys', type=str,
        default="job_keys",
        help="name of attribute with job keys"
    )

    return p
