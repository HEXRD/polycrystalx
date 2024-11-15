"""Run a job with MPI"""
import sys
import argparse
import pickle

from . import get_input_module
from polycrystalx import processes


def main():
    """Run a job in MPI"""
    p = argparser(*sys.argv)
    args = p.parse_args()

    user_module = get_input_module(args.input_module)
    if not hasattr(user_module, "get_job"):
        raise AttributeError('module has no "get_job" attribute')

    with open(args.key_file, "rb") as f:
        key = pickle.load(f)

    job = user_module.get_job(key)
    processes.run(job)


def argparser(*args):

    p = argparse.ArgumentParser(
        description="run a job"
    )
    p.add_argument(
        'input_module', type=str,
        help="name of batch input module"
    )

    p.add_argument(
        'key_file', type=str,
        help="name of file with serialized JobKey instance"
    )

    return p
