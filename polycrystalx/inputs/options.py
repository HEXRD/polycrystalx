"""Run time options

This is the first trial as an options input module. It only has solver options.
In the future we plan to have more an option for element degree and also options for
output selection.

The output selection is more complicated because it depends on the material process.
So the default set of options may in fact be a dictionary of default options for
each process.
"""
from collections import namedtuple

default_petsc={
    "ksp_type": "cg",
    "ksp_rtol": 1e-6,
    "ksp_atol": 1e-10,
    "ksp_max_it": 5000,
    "pc_type": "jacobi",
}

_fields = ["name", "solver_options"]
_defaults = ["dflt_opts", default_petsc]

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
