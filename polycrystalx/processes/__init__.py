"""This is the module for defining and executing model processes"""
import os

from .linear_elasticity import LinearElasticity
from ..utils import setup_output


process_dict = {}
for p in (LinearElasticity,):
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
