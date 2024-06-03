"""Tests for inputs"""
import numpy as np

from polycrystalx import inputs


def test_mesh_inputs():

    mesh_input = inputs.mesh.Mesh(
        name='test',
        source="box",
    )

    assert mesh_input.extents is None
    assert mesh_input.divisions is None
    assert mesh_input.celltype is None
    assert mesh_input.file is None
    assert mesh_input.boundary_sections == []


def test_deformation_inputs():
    defm_input = inputs.deformation.LinearElasticity(
        name=("test"),
    )
    assert defm_input.force_density is None
    assert defm_input.plastic_distortion is None
    assert defm_input.displacement_bcs == []
    assert defm_input.traction_bcs == []
