"""Polycrystal input templates"""
from collections import namedtuple


Polycrystal = namedtuple(
    "Polycrystal", ["name", "polycrystal", "use_meshtags"],
    defaults=[False]
)
Polycrystal.__doc__ = """Polycrystal input for Elasticity

The main input is the microstructure, which is an instance of Microstructure
(see `notes`_). However, the grain IDs may be already computed and can
be read in with the mesh, as in gmsh using neper. In that case the
`polycrystal` input has only the orientations.

Parameters
----------
name: str
    name of polycrystal input
polycrystal: (microstructure) Microstructure
    microstructure
use_meshtags: bool
    if True, grain IDs are read with the mesh as cell_tags (gmsh option)

.. _notes:

Notes
-----
Here, we use a companion package, *microstructure*, to define materials and
polycrystal configurations.
"""
