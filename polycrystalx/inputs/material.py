"""Material input templates

Currently these are the same format for all processes: the name of the material
list and a list of instances of the appropriate material
"""
from collections import namedtuple


MaterialList = namedtuple(
    "MaterialLisst", ["name", "materials"]
)
MaterialList.__doc__ = """List of Materials

Parameters
----------
name: str
    name of material list
materials: list
    list of material instances
"""

LinearElasticity = MaterialList
HeatTransfer = MaterialList
