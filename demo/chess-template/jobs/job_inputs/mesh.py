"""Mesh Input Module"""

import numpy as np

from polycrystalx import inputs

from ..data import DATA_DIR


# We are using an abbreviated section of the data set. The data set is cylindrical, and
# we are using a box. The `REAL_EXTENTS` is that for the actual data set, and `EXTENTS`
# is what we are using in this demo.

REAL_EXTENTS = np.array([
    [-0.486, 0.486],
    [-0.250, 0.250],
    [-0.486, 0.486],
])

EXTENTS = np.loadtxt(DATA_DIR / "extents.txt")


def get_mesh_input(key):
    """Return a named mesh input"""
    return MeshInput(key).mesh_input


class MeshInput:
    """Builds mesh input for meshx

    Parameters:
    ----------
    key: str
         basename of the mesh file
    """
    def __init__(self, key):
        self.key = key

    @property
    def name(self):
        return self.key

    @property
    def file(self):
        return DATA_DIR / f"{self.name}.msh"

    @property
    def mesh_input(self):
        return inputs.mesh.Mesh(
            name=self.name,
            source="gmsh",
            file=self.file,
            extents=EXTENTS,
            boundary_sections=[]
        )
