"""Batch jobs"""
import itertools

from polycrystalx import inputs


from .job_inputs import (
    get_material_input,
    get_microstructure_input,
    get_mesh_input,
    get_deformation_input,
)


suite = "torsion"
process = "linear-elasticity"


def get_job(key):
    matl, poly, mesh, defm = key
    matl_input = get_material_input(matl)
    poly_input = get_microstructure_input(poly)
    mesh_input = get_mesh_input(mesh)
    defm_input = get_deformation_input(defm)

    return inputs.job.Job(
        suite = suite,
        process = process,
        mesh_input = mesh_input,
        material_input = matl_input,
        polycrystal_input = poly_input,
        deformation_input = defm_input
    )


# Define the suite of jobs.

matl_keys = ["lshr_660C"]
poly_keys = ["ms"]
mesh_keys = ["seeds-025"]
defm_keys = list(range(13))

job_keys = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)
