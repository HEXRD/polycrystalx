"""Batch jobs"""
import itertools

from polycrystalx import inputs

from .job_data import (
    get_material_input,
    get_polycrystal_input,
    get_mesh_input,
    get_deformation_input,
)


suite = "linear-single-crystal"
process = "linear-elasticity"


def get_job(key):
    matl, poly, mesh, defm = key
    matl_input = get_material_input(matl)
    poly_input = get_polycrystal_input(poly)
    mesh_input = get_mesh_input(mesh)
    defm_input = get_deformation_input(defm, matl_input, poly_input)

    return inputs.job.Job(
        suite = suite,
        process = process,
        mesh_input = mesh_input,
        material_input = matl_input,
        polycrystal_input = poly_input,
        deformation_input = defm_input
    )


# Define suites of jobs.


# All matrices A with a single cubic material, mesh and polcyrstal.

all_A = itertools.product(
    ["cubic-321"],
    [(2, 45)],
    [(30, 30, 30)],
    itertools.product(["full"], [0, 1, 2], [0, 1, 2])
)

# All materials for a single deformation.

all_materials =  itertools.product(
    ["identity-iso", "iso-21", "cubic-321", "cubic-211", "cubic-112", "cubic-321"],
    [(0, 0)],
    [(30, 30, 30)],
    [("full", 0, 0)],
)

# Vary cyrstal orientation for a single deformation with traction BCs.

vary_orientation = itertools.product(
    ["cubic-321"],
    [(0, 0), (0, 30), (1, 45), (2, 65)],
    [(30, 30, 30)],
    [("zmax-traction", 0, 0)],
)

# Vary boundary conditions.

vary_bcs = itertools.product(
    ["cubic-321"],
    [(1, 33)],
    [(30, 30, 30)],
    [
        ("full", 0, 0), ("zmax-traction", 0, 0), ("zmax-traction-z", 0, 0),
        ("zmax-traction-xy", 0, 0)
    ],
)
