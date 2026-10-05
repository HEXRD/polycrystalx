"""Run time options"""
from collections import namedtuple

_default_petsc={
    "ksp_type": "cg",
    "ksp_rtol": 1e-6,
    "ksp_atol": 1e-10,
    "ksp_max_it": 5000,
    "pc_type": "jacobi",
}

_fields = ["name", "solver_options"]
_defaults = [_default_petsc]

Options = namedtuple(
    "Options", _fields, defaults=_defaults
)
Options.__doc__ = """Options

Parameters
----------
name: str
    name of option group
solver_options: dict
    dictionary of PETSc options for Krylov solver
"""
