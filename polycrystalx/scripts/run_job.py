"""Run a job"""
import sys
import argparse
from pathlib import Path

from polycrystalx import processes
from . import get_input_module


def main():
    """run a single job"""
    p = argparser(*sys.argv)
    args = p.parse_args()

    user_module = get_input_module(args.input_module)
    if not hasattr(user_module, "job"):
        raise AttributeError('module has no "jobs" attribute')

    processes.run(user_module.job)


def argparser(*args):

    p = argparse.ArgumentParser(
        description="run a job"
    )
    p.add_argument(
        'input_module', type=str,
        help="module with inputs (name or path)"
    )

    return p
