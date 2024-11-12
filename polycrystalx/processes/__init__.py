"""This is the module for defining and executing model processes"""
import os

from .linear_elasticity import LinearElasticity
from .heat_transfer import HeatTransfer

from ..utils import setup_output


processes = (LinearElasticity, HeatTransfer)
process_dict = {}
for p in processes:
    process_dict[p.name] = p


def run(job):
    """Run a job

    PARAMETERS
    ----------
    job: inputs.job.Job
       the job to run
    """
    setup_output(job.output_directory)
    process = process_dict[job.process](job)
    process.run()
