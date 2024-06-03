"""Mesh input class"""
from collections import namedtuple


BoundarySection = namedtuple("BoundarySection", ["name", "on_section"])
BoundarySection.__doc__ = """Section of boundary

Parameters
----------
name: str
    name of the boundary section
on_section: python function
    boolean valued function of position array x (of shape (d, n)),
    where `d` is the dimension and `n` is the number of points, giving
    True values for points on this boundary section
"""


Mesh = namedtuple(
    "Mesh", ["name", "source", "extents", "divisions", "celltype", "file",
             "boundary_sections"],
    defaults=4 * [None] + [[]]
)
Mesh.__doc__ = """Mesh inputs for elasticity

The mesh input is driven by the `source` keyword, which can be "box", "xdmf" or
"gmsh". For "box" source, you will need the `divisions`, `extents` and
`celltype` keywords. For "xdmf" or "gmsh", you will need the `file` keyword.

Parameters
----------
name: str
    name for this input specification
source: {"box", "xdmf", "gmsh"}
    source of mesh
extents: array(d, 2), where d is dimension (or None)
    for a box mesh, the min and max of each coordinate
divisions: {1-tuple | d-tuple of ints}
    for a box mesh, the subdivisions in each direction
celltype: str
    for a box mesh, cell type, one of "tetrahedron", "hexahedron",
    "quadrilateral", "triangle" or "interval"
file: str
    name of mesh XDMF file when `source` is "xdmf" or "gmsh"
boundary_sections: list
    list of `BoundarySection` instances
"""
