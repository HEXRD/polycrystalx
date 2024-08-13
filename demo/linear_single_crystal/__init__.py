"""Inputs Module for lienar single crystal"""
import numpy as np

from polycrystalx import inputs

from .data import (
    get_material_input, get_polycrystal_input, get_mesh_input,
    get_deformation_input, matl_dict
)
from .batch import suite, process


# ========== Material
#
# Materials are in matl_data.matl_dict.
matl_name = "cubic-211"
material_input = get_material_input(matl_name)

# ========== Microstructure

polycrystal_input = get_polycrystal_input(0, 0)

# ========== Mesh input

divs = (20, 20, 20)

mesh_input = get_mesh_input(divs)

# ========== Deformation Inputs
#
defm_key = ("full", 1, 1)
defm_key = ("zmax-traction", 1, 1)
defm_key = ("zmax-traction-z", 1, 1)
defm_key = ("zmax-traction-xy", 1, 1)

deformation_input = get_deformation_input(defm_key, matl_name)

# ========== Current Job

job = inputs.job.Job(
    suite = suite,
    process = process,
    mesh_input = mesh_input,
    material_input = material_input,
    polycrystal_input = polycrystal_input,
    deformation_input = deformation_input
)
