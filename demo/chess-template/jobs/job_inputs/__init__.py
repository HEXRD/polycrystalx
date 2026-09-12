"""Batch tools for  ... """
from polycrystalx import inputs

from .material import get_material_input
from .microstructure import get_microstructure_input
from .mesh import get_mesh_input
from .deformation import get_deformation_input

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
