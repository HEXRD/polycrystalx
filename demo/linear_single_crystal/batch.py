"""Batch jobs"""
import itertools

from polycrystalx import inputs

from .data import (
    get_material_input, get_polycrystal_input, get_mesh_input,
    get_deformation_input, matl_dict
)

suite = "linear_single_crystal"
process = "linear-elasticity"


def get_job(key):
    matl, poly, mesh, defm = key
    matl_input = get_material_input(matl)
    poly_input = get_polycrystal_input(*poly)
    mesh_input = get_mesh_input(mesh)
    defm_input = get_deformation_input(defm, matl)

    return inputs.job.Job(
        suite = suite,
        process = process,
        mesh_input = mesh_input,
        material_input = matl_input,
        polycrystal_input = poly_input,
        deformation_input = defm_input
    )


# Below are some job suites.

matl_keys = list(matl_dict.keys())
poly_keys = [(0, 0)]
mesh_keys = [(10, 20, 30)]
defm_keys = [("full", 1, 1)]

# This suite runs all materials on the same mesh with displacements
# corresponding to an xx-strain field.
all_materials_xx = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)


matl_keys = ["cubic-211"]
poly_keys = [(0, 0)]
mesh_keys = [(10, 20, 30)]
defm_keys = [
    ("zmax-traction", 3, 3),
    ("zmax-traction-z", 3, 3),
    ("zmax-traction-xy", 3, 3)
]

# This suite runs various traction boundary conditions on the top surface with a
# fixed mesh and a cubic material.
cubic_traction_z = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)


# This is the default job suite.
job_keys = cubic_traction_z
