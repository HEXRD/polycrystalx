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


_Mesh = namedtuple(
    "_Mesh", ["name", "source", "extents", "divisions", "celltype", "file",
             "boundary_sections"],
    defaults=4 * [None] + [[]]
)


class Mesh(_Mesh):
    """Mesh inputs for elasticity

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

    def __init__(self, *args, **kwargs):
        super(__class__, self).__init__()
        self.check_inputs()

    def check_inputs(self):

        if self.source == "box":
            # Needs  "extents", "divisions", "celltype".
            bad = (
                (self.extents is None) or
                (self.divisions is None) or
                (self.celltype is None)
            )
            if bad:
                emsg = (
                    'For "box" mesh, you must also provide "extents", '
                    '"divisions", and "celltype".'
                )
                raise RuntimeError(emsg)

        elif self.source == "gmsh":
            if (self.file is None):
                emsg = 'For "gmsh" mesh, you must also provide "file".'
                raise RuntimeError(emsg)

        elif self.source == "xdmf":
            if (self.file is None):
                emsg = 'For "xdmf" mesh, you must also provide "file".'
                raise RuntimeError(emsg)

        else:
            emsg = (
                'no valid "source" type was found: '
                '"source" must be one of "box" | "gmsh" | "xdmf"'
            )
            raise RuntimeError(emsg)
