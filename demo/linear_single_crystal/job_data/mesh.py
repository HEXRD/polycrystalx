"""Mesh Input Module"""

import numpy as np

from polycrystalx import inputs


extents = np.array([
    [0, 1],
    [0, 2],
    [0, 3],
])


def get_mesh_input(key):
    """Return a named mesh input"""
    return MeshInput(key).mesh_input


class MeshInput:
    """Builds mesh input for meshx

    Parameters:
    ----------
    key: 3-tuple
       divisions in each coordinate direction
    """

    def __init__(self, key):
        self.divs = key

    @property
    def mesh_input(self):
        return inputs.mesh.Mesh(
            name=self.name,
            source="box",
            extents=extents,
            divisions=self.divs,
            celltype="tetrahedron",
            boundary_sections=[
                inputs.mesh.BoundarySection("not-zmax", self.not_on_zmax)
            ]
        )

    @property
    def name(self):
        return "x".join([str(d) for d in self.divs])

    def not_on_zmax(self, x):
        """Boundary not on zmax

        *** CHECk THIS
        is it simpler to say; logical_or(xmin, xmax, ymin, ymax, zmin)
        """
        xmin, xmax = extents[0]
        ymin, ymax = extents[1]
        zmax = extents[2, 1]

        on_zmax = np.isclose(x[2], zmax)
        near_x = np.logical_or(np.isclose(x[0], xmin), np.isclose(x[0], xmax))
        near_y = np.logical_or(np.isclose(x[1], ymin), np.isclose(x[1], ymax))
        near_zmax_edge = np.logical_and(np.logical_or(near_x, near_y), on_zmax)
        return np.logical_or(~on_zmax, near_zmax_edge)
