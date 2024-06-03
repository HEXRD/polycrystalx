"""Run time options"""
from collections import namedtuple


Options = namedtuple(
    "Options", ["name", "tolerance", "maxiter", "save_pvd", "save_hdf5", "outdir"]
)
Options.__doc__ = """Options
"""

_defaults = []
