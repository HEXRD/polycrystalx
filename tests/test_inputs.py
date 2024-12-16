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


class TestDeformationInputs:

    def test_linear_elasticity(self):

        defm_input = inputs.deformation.LinearElasticity(
            name=("test"),
        )
        assert defm_input.name == "test"
        assert defm_input.force_density is None
        assert defm_input.plastic_distortion is None
        assert defm_input.displacement_bcs == []
        assert defm_input.traction_bcs == []

    def test_heat_transfer(self):

        defm_input = inputs.deformation.HeatTransfer(
            name=("test-heat-transfer"),
        )
        assert defm_input.name == "test-heat-transfer"
        assert defm_input.body_heat is None
        assert defm_input.temperature_bcs == []
        assert defm_input.flux_bcs == []
