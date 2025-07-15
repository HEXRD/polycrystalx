"""Batch jobs"""

from polycrystalx import inputs

from .jobdata import (
    get_material_input,
    get_polycrystal_input,
    get_mesh_input,
    get_deformation_input,
    matl_dict
)

suite = "thermal_expansion"
process = "linear-elasticity"

# The `get_job` function is required for running batch jobs.

def get_job(key):
    matl, mesh, poly, defm = key
    matl_input = get_material_input(matl)
    mesh_input = get_mesh_input(mesh)
    poly_input = get_polycrystal_input(poly)
    defm_input = get_deformation_input(defm, matl)

    return inputs.job.Job(
        suite = suite,
        process = process,
        mesh_input = mesh_input,
        material_input = matl_input,
        polycrystal_input = poly_input,
        deformation_input = defm_input
    )


# Next, we set up the test suite.

ident = "identity-iso"
ti64 = "ti-64-bar-RT"
poly_key = None
mesh_key = (20, 20, 20)
defm_keys = ["zero", "match", "linear"]

job_keys = [
    (ident, mesh_key, poly_key, "zero"),
    (ident, mesh_key, poly_key, "match"),
    (ti64, mesh_key, poly_key, "linear"),
]
