"""Mesh and function spaces"""
import collections

import numpy as np
import dolfinx

from mpi4py import MPI

from . boundary import initialize_boundary_dict


class MeshLoader:
    """Load an input Mesh specification

    Parameters
    ----------
    userinput: inputs.mesh.Mesh
       input mesh specification
    """
    _CT = dolfinx.mesh.CellType
    celldict = {
        "tetrahedron": _CT.tetrahedron,
        "hexahedron": _CT.hexahedron,
        "pyramid": _CT.pyramid,
        "quadrilateral": _CT.quadrilateral,
        "triangle""triangle": _CT.triangle,
        "interval": _CT.interval,
        "prism": _CT.prism,
        "point": _CT.point,
    }

    def __init__(self, userinput):
        self.userinput = userinput

        self._extents = userinput.extents
        self._mesh = self._make_mesh()
        self.tdim = self._mesh.topology.dim
        self.bdim = self.tdim - 1
        self._bsecs = userinput.boundary_sections
        #
        self._boundary_dict = self._make_boundary_dict()

    def _make_mesh(self):
        self.cell_tags = None
        self.facet_tags = None

        if self.userinput.source == "box":
            return self._make_box_mesh()
        elif self.userinput.source == "xdmf":
            return self._read_xdmf()
        elif self.userinput.source == "gmsh":
            return self._read_gmsh()
        else:
            s = self.userinput.source
            raise ValueError(f"mesh source ({s}) not available")

    @property
    def extents(self):
        """Return extents of mesh"""
        if self._extents is not None:
            return np.array(self.userinput.extents)
        else:
            return None

    def _make_box_mesh(self):
        ext = self.extents
        udiv = self.userinput.divisions
        cell = self.userinput.celltype
        dim = len(ext)
        # `divsions` can an int, a tuple, a list or a numpy array.
        if isinstance(udiv, int):
            divs = dim * (x,)
        elif isinstance(udiv, collections.abc.Sequence) or\
             isinstance(udiv, np.ndarray):
            if len(udiv) == dim:
                divs = udiv
            elif len(udiv) == 1:
                divs = dim * udiv
            else:
                msg = "divisions has wrong length: "\
                    "should be 1 or d (mesh dimension)"
                raise ValueError(msg)
        else:
            raise ValueError(f"divisions spec not recognized: {udiv}")

        msh = dolfinx.mesh.create_box(
            comm=MPI.COMM_WORLD,
            points=(ext[:, 0], ext[:, 1]),
            n=divs,
            cell_type=self.celldict[cell],
        )
        return msh

    def _read_xdmf(self):
        fname = self.userinput.file
        with dolfinx.io.XDMFFile(MPI.COMM_WORLD, fname, "r") as f:
            m = f.read_mesh()
            try:
                self.cell_tags = f.read_meshtags(m, "grain-ids")
            except:
                pass

        return m

    def _read_gmsh(self):
        fname = self.userinput.file
        msh, self.cell_tags, self.facet_tags = dolfinx.io.gmshio.read_from_msh(
            fname, MPI.COMM_WORLD, 0
        )
        return msh


    def _make_boundary_dict(self):
        d = initialize_boundary_dict(self.mesh, self.extents)

        for bsec in self._bsecs:
            facets = dolfinx.mesh.locate_entities_boundary(
                self.mesh, self.bdim, bsec.on_section
            )
            d[bsec.name] = facets

        return d

    @property
    def mesh(self):
        """Return the dolfinx Mesh"""
        return self._mesh

    @property
    def boundary_dict(self):
        """Return the dicitonary of boundary sections"""
        return self._boundary_dict
