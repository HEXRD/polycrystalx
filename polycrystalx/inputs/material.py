"""Material input templates"""
from collections import namedtuple


LinearElasticity = namedtuple(
    "LinearElasticity", ["name", "materials"]
)
LinearElasticity.__doc__ = """Material input for Elasticity

Parameters
----------
name: str
    name of material
materials: list
    list of elastic single crystal instances
"""


Slip = namedtuple(
    "Slip", ["name", "materials"]
)
Slip.__doc__ = """Material input for Elasticity

Parameters
----------
name: str
    name of material
materials: list
    list of single crystal slip instances
"""
