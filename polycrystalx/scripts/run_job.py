"""Run a job"""
import sys
import argparse
from pathlib import Path

from dolfinx import log

from polycrystalx import utils, processes


def main():
    """run a single job"""
    p = argparser(*sys.argv)
    args = p.parse_args()

    log.set_log_level(log.LogLevel.INFO)
    user_module = utils.get_input_module(args.input_module)
    log.set_log_level(log.LogLevel.WARNING)
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


if __name__ == "__main__":
    #
    #  Run problem
    #
    main()
