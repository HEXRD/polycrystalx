"""Utilities for handling input"""
import sys
import os
from pathlib import Path
import time
from importlib import import_module

import numpy as np

from dolfinx import log
from dolfinx.fem import assemble_scalar

from mpi4py import MPI

from .xdmffile_ext import XDMFFile_Ext


myrank = MPI.COMM_WORLD.rank
commsize = MPI.COMM_WORLD.size

OFF = log.LogLevel.OFF
ERROR = log.LogLevel.ERROR
INFO = log.LogLevel.INFO
WARNING = log.LogLevel.WARNING
TAG_SYNC = 500


def mpi_sync(from_0=True):
    """Sync processes

    If we want processes to wait for process 0, we use from_0=True. If we
    want process 0 to wait for the others, we use from_0=False.

    Parameters
    ----------
    from_0: bool, default=True
       if true, node 0 sends to all other nodes, else other send to 0
    """
    # Send random number each time so we can check for matching values.
    loglevel = WARNING # debugging
    msg = np.random.random(1)
    for i in range(1, commsize):
        if from_0:
            if myrank == 0:
                MPI.COMM_WORLD.send(msg, i, tag=TAG_SYNC)
                log.log(loglevel,
                        f"0: sending sync from {myrank} to {i}: {msg}"
                        )
            if myrank == i:
                rmsg = MPI.COMM_WORLD.recv(source=0, tag=TAG_SYNC)
                log.log(loglevel, f"{myrank}: receiving sync from 0: ({rmsg})")
        else:
            if myrank == 0:
                rmsg = MPI.COMM_WORLD.recv(source=i, tag=TAG_SYNC)
                log.log(loglevel, f"0: receiving sync from {i}: ({rmsg})")
            if myrank == i:
                MPI.COMM_WORLD.send(msg, 0, tag=TAG_SYNC)
                log.log(loglevel, f"{myrank}: sending sync to 0: {msg}")


def get_input_module(imod):
    """Process input string to get module

    Parameters
    ----------
    imod: str or Path
       name of input module
    """
    # Expect the user to run from directory containing the input module.
    sys.path.insert(0, str(Path.cwd()))
    # Allow user to specify module by path
    if os.path.exists(imod):
        if imod.endswith("/"):
            imod = imod[:-1]
        path, ext = os.path.splitext(imod)
        imod = path.replace("/", ".")

    log.log(INFO, f"loading user input module: {imod}")
    user_inputs = import_module(imod)

    return user_inputs


def setup_output(outdir):
    """Make output directory if needed

    outdir: str or Path
        name of output directory
    """
    mpi_sync(from_0=False)
    if not os.path.exists(outdir):
        if myrank == 0:
            log.log(WARNING, f"creating output directory: {outdir}")
            os.makedirs(outdir)
            time.sleep(0.1)
        mpi_sync()
    log.log(WARNING, f"{myrank}: changing directory to {outdir}")

    try:
        os.chdir(outdir)
    except:
        log.log(WARNING, f"{myrank}: failed to find output directory")


def grain_volumes(comm, gv_form, indicator, grain_cells):
    """Compute grain volumes

    comm: MPI communicator
       main communicator for this problem
    gv_form: dolfinx Form
       form for computing grain volume
    indicator: dolfinx Function
        inidcator function for the gv_form
    grain_cells: dict
        dictionary giving array of cells for each grain

    Returns
    -------
    arary
       array of grain volumes
    """
    vols = np.zeros(ng := len(grain_cells))
    indicator.x.array[:] = 0.
    for g in range(ng):
        if g > 0:
            indicator.x.array[gcells] = 0
        gcells = grain_cells[g]
        indicator.x.array[gcells] = 1

        value = 0. if len(gcells) == 0 else assemble_scalar(gv_form)
        vols[g] = comm.allreduce(value, op=MPI.SUM)

    return vols


def grain_integrals(comm, gi_form, indicator, func, grain_cells, f):
    """Compute grain integrals of scalar function

    comm: MPI communicator
       main communicator for this problem
    gi_form: dolfinx Form
       form for computing grain volume
    indicator: dolfinx Function
        inidcator function for the gv_form
    func: dolfinx Function
        function to integrate over the grains
    grain_cells: dict
        dictionary giving list of cells for each grain

    Returns
    -------
    arary
       array of grain integrals
    """
    integrals = np.zeros(ng := len(grain_cells))
    func.x.array[:] = f.x.array
    indicator.x.array[:] = 0.
    for g in range(ng):
        if g > 0:
            indicator.x.array[gcells] = 0
        gcells = grain_cells[g]
        indicator.x.array[gcells] = 1

        value = 0. if len(gcells) == 0 else assemble_scalar(gi_form)
        integrals[g] = comm.allreduce(value, op=MPI.SUM)

    return integrals
