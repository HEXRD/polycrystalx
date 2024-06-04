"""Batch jobs"""
import itertools

from polycrystalx import inputs

from .data import (
    get_material_input, get_polycrystal_input, get_mesh_input,
    get_deformation_input, matl_dict
)

suite = "linear_single_crystal"
process = "linear-elasticity"

# The `get_job` function is required for running batch jobs.

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


# Next, we set up the iterator for job keys.

matl_keys = matl_dict.keys()
poly_keys = [(0, 0)]
mesh_keys = ((10, 20, 30),)
defm_keys = itertools.product(
    ("full", "zmax-traction", "zmax-traction-z", "zmax-traction-xy"),
    range(3), range(3)
)
job_keys = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)
