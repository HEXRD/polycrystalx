"""Run suite of Jobs"""
import sys
import argparse
import tempfile
import pathlib
import subprocess
import pickle

import numpy as np

from polycrystalx import utils


RUN_MPI = pathlib.Path(__file__).parent / "run_mpijob.py"


def main():
    """run a suite of jobs"""
    p = argparser(*sys.argv)
    args = p.parse_args()

    user_module = utils.get_input_module(args.input_module)
    if not hasattr(user_module, "job_keys"):
        raise AttributeError('module has no "job_keys" attribute')


    for job in user_module.job_keys:
        # Pickle to a temp file.
        fp = tempfile.NamedTemporaryFile(delete=False)
        with open(fp.name, "wb") as f:
            pickle.dump(job, f)

        print("\n===== New Job Starting", flush=True)
        # Now run MPI job.
        cmd = [
            "mpirun", "-np", str(args.n),
            "python", str(RUN_MPI), args.input_module, fp.name
        ]
        print(f"command is: {' '.join(cmd)}")
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

    return p


if __name__ == "__main__":
    #
    #  Run problem
    #
    main()
