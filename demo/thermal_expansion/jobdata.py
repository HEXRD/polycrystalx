"""Batch tools for thermal expansion

References
----------
* Ti-64 data from matweb:
  Titanium Ti-6Al-4V (Grade 5), Annealed Bar
  https://www.matweb.com/search/DataSheet.aspx?MatGUID=10d463eb3d3d4ff48fc57e0ad1037434
"""
import numpy as np

from polycrystalx import inputs

from polycrystal.elasticity.single_crystal import SingleCrystal
from polycrystal.microstructure.single_crystal import (
    SingleCrystal as SingleCrystalMicro
)

from . import defm


# Material Data


def get_material_input(name):
    """Return material input by name"""
    return inputs.material.LinearElasticity(
        name=name, materials=[matl_dict[name]]
    )


""" CTE Data from matweb
20 - 100 C: 8.6 u/K
20 - 315 C: 9.2 u/K
20 - 650 C: 9.7 u/K
"""
E = 114 # GPa
nu = 0.33
cte_RT = 8.6e-6

matl_dict = {
    "identity-iso": SingleCrystal.from_K_G(1/3, 1/2, cte=cte_RT),
    "ti-64-bar-RT": SingleCrystal.from_E_nu(E, nu, cte=cte_RT),
}

# Polycrystal Data

def get_polycrystal_input(arg):
    """single crystal aligned with sample axes"""
    ms = SingleCrystalMicro(np.identity(3))
    return inputs.polycrystal.Polycrystal(
        name="sx",
        polycrystal=ms
    )


# Mesh Data

extents = np.array([
    [0, 1],
    [0, 1],
    [0, 1],
])


def get_mesh_input(key):
    divs = key
    return  inputs.mesh.Mesh(
        name=mesh_name(divs),
        source="box",
        extents=extents,
        divisions=divs,
        celltype="tetrahedron",
    )

def mesh_name(divs):
    """generate mesh name from divisions"""
    return "x".join([str(d) for d in divs])


# Deformation Data

def get_deformation_input(defm_key, matl_key):
    matl = matl_dict[matl_key]
    return defm.get_input(defm_key, matl)
