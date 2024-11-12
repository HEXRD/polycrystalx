"""MPI Utilities"""
import numpy  as np

from dolfinx import log
from mpi4py import MPI


myrank = MPI.COMM_WORLD.rank
commsize = MPI.COMM_WORLD.size
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
    loglevel = log.LogLevel.INFO
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
