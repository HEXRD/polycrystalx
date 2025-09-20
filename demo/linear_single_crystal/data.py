"""Batch tools for linear single crystal"""
import numpy as np

from polycrystalx import inputs

from polycrystal.elasticity.single_crystal import SingleCrystal
from polycrystal.elasticity.moduli_tools import Cubic
from polycrystal.microstructure.single_crystal import (
    SingleCrystal as SingleCrystalMicro
)
from polycrystal.orientations.quaternions import from_exp, to_rmats

from . import defm_data


# Material Data

def get_material_input(name):
    """Return material input by name"""
    return inputs.material.LinearElasticity(
        name=name, materials=[matl_dict[name]]
    )


def cubic_from_KGdGs(K, Gd, Gs):
    c = Cubic.from_K_Gd_Gs(K, Gd, Gs)
    return SingleCrystal("cubic", (c.c11, c.c12, c.c44))


matl_dict = {
    "identity-iso": SingleCrystal.from_K_G(1/3, 1/2),
    "identity-cub": SingleCrystal("cubic", (1.0, 0.5, 0.25)),
    "iso-21": SingleCrystal.from_K_G(2, 1),
    "cubic-211": cubic_from_KGdGs(2, 1, 1),
    "cubic-121": cubic_from_KGdGs(1, 2, 1),
    "cubic-112": cubic_from_KGdGs(1, 1, 2),
}

# Polycrystal Data

def get_polycrystal_input(ax, ang_d):
    """Return polycrystal input for single crystal

    Parameters
    ----------
    ax: array(3)
       axis of rotation
    ang_d: float
       angle of rotation in degrees

    Returns
    -------
    inp: inputs.polycrystal.Polycrystal
       the polycrystal input
    """
    ms = get_micro(ax, ang_d)
    name = "-".join(("single-crystal", "xyz"[ax] + str(ang_d)))

    return inputs.polycrystal.Polycrystal(
        name=name,
        polycrystal=ms
    )


def rmatrix(angle, axis):
    a = np.atleast_2d(np.array(axis))
    n = a/np.linalg.norm(a,1)
    return to_rmats(from_exp(angle * n))


def get_micro(ax, ang):
    """Get single crystal microstructure with given angle/axis

    ax - 0, 1, 2 for x, y, z axis
    ang - angle in degrees
    """
    axis = np.zeros(3)
    axis[ax] = 1.
    ang_r = np.radians(ang)

    return SingleCrystalMicro(rmatrix(ang_r, axis))


# Mesh Data

extents = np.array([
    [0, 1],
    [0, 2],
    [0, 3],
])


def get_mesh_input(key):
    divs = key
    return  inputs.mesh.Mesh(
        name=name(divs),
        source="box",
        extents=extents,
        divisions=divs,
        celltype="tetrahedron",
        boundary_sections=[
            inputs.mesh.BoundarySection("not-zmax", not_on_zmax)
        ]
    )

def name(divs):
    """generate mesh name from divisions"""
    return "x".join([str(d) for d in divs])


def not_on_zmax(x):
    """Boundary not on zmax"""
    xmin, xmax = extents[0]
    ymin, ymax = extents[1]
    zmax = extents[2, 1]

    on_zmax = np.isclose(x[2], zmax)
    near_x = np.logical_or(np.isclose(x[0], xmin), np.isclose(x[0], xmax))
    near_y = np.logical_or(np.isclose(x[1], ymin), np.isclose(x[1], ymax))
    near_zmax_edge = np.logical_and(np.logical_or(near_x, near_y), on_zmax)
    return np.logical_or(~on_zmax, near_zmax_edge)


# Deformation Data

def get_deformation_input(key, matl_name):
    matl = get_material_input(matl_name).materials[0]
    return defm_data.get_input(key, matl)
