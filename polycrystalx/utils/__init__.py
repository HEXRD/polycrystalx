"""Utilities for handling input"""
import os
import time

import numpy as np

from dolfinx import log
from dolfinx.fem import assemble_scalar

from .xdmffile_ext import XDMFFile_Ext
from .mpi import MPI, mpi_sync, myrank


def setup_output(outdir):
    """Make output directory if needed

    outdir: str or Path
        name of output directory
    """
    mpi_sync(from_0=False)
    if not os.path.exists(outdir):
        if myrank == 0:
            log.log(log.LogLevel.INFO, f"creating output directory: {outdir}")
            os.makedirs(outdir)
            time.sleep(0.1)
        mpi_sync()
    log.log(log.LogLevel.INFO, f"{myrank}: changing directory to {outdir}")

    try:
        os.chdir(outdir)
    except:
        raise RuntimeError(f"{myrank}: failed to find output directory")


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
