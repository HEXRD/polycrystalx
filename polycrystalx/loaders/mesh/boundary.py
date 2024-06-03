"""Tools for loading boundary sections"""
import numpy as np
import dolfinx


class XYZBoundary:
    """Boundary of simple box

    This is a convenience class for creating boundary definitions for
    rectangular geometries with normals in the X, Y, or Z direction. We can
    generalize this to any flat boundary by using a function of point and
    normal.

    Parameters
    ----------
    component: {self.X, self.Y, self.Z} or {0, 1, 2}
        direction of normal
    value: float
        value defining the surface
    sign: {self.MAX, self.MIN} or {1, -1}
        direction of normal
    EPS: float, default = 1e-14
        tolerance

    Returns
    -------
    function
        boolean function of position describing flat boundary
    """

    X, Y, Z = 0, 1, 2
    MAX, MIN = 1, -1

    def __init__(self, component, value, sign, EPS=1e-14):
        self.component = component
        self.value =  value
        self.sign = sign
        self.EPS = EPS

    def on_boundary(self, x):
        """Return True for  points lying on the flat boundary

        Parameters
        ----------
        x: array(d, n)
            array of `n` points in `d` dimensions

        Returns
        -------
        bool array (n)
            array True for points on the specified boundary
        """
        xi = x[self.component, :]

        return (xi - self.value) * self.sign > -self.EPS


def initialize_boundary_dict(msh, extents):
    """Initialize boundary dictionary

    This initializes the boundary dictionary to include the whole boundary
    and flat surfaces if extents are given. The whole boundary has key
    "boundary", and the six sides have keys "xmin", "ymin", "zmin", "xmax",
    "ymax" and "zmax".

    Parameters
    ----------
    extents: array (d, 2)
       ranges for each coordinate

    Returns
    -------
    dict
       dictionary with the boundary name as key and the array of facets as
       value
    """
    Bndry = XYZBoundary
    X, Y, Z = 0, 1, 2
    MIN, MAX = 0, 1

    bdim = msh.topology.dim - 1
    anycell = lambda x:  np.ones(x.shape[1], dtype=bool)
    allcells = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, anycell
    )
    d = {"boundary": allcells}

    if extents is None: return d

    # If extents are given, add one section for each flat surface.

    on_xmin = Bndry(Bndry.X, extents[X, MIN], Bndry.MIN).on_boundary
    on_ymin = Bndry(Bndry.Y, extents[Y, MIN], Bndry.MIN).on_boundary
    on_zmin = Bndry(Bndry.Z, extents[Z, MIN], Bndry.MIN).on_boundary
    #
    on_xmax = Bndry(Bndry.X, extents[X, MAX], Bndry.MAX).on_boundary
    on_ymax = Bndry(Bndry.Y, extents[Y, MAX], Bndry.MAX).on_boundary
    on_zmax = Bndry(Bndry.Z, extents[Z, MAX], Bndry.MAX).on_boundary

    facets_xmin = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_xmin
    )
    facets_ymin = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_ymin
    )
    facets_zmin = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_zmin
    )
    facets_xmax = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_xmax
    )
    facets_ymax = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_ymax
    )
    facets_zmax = dolfinx.mesh.locate_entities_boundary(
        msh, bdim, on_zmax
    )

    d["xmin"] = facets_xmin
    d["ymin"] = facets_ymin
    d["zmin"] = facets_zmin
    d["xmax"] = facets_xmax
    d["ymax"] = facets_ymax
    d["zmax"] = facets_zmax

    return d
