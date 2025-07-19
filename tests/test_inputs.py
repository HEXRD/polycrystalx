"""Tests for inputs"""
import numpy as np
import pytest

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
        assert defm_input.thermal_expansion is None
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


class TestFunctionInputs:

    def test_check_inputs(self):

        with pytest.raises(RuntimeError, match="constant function"):
            inp = inputs.function.Function(source="constant")

        with pytest.raises(RuntimeError, match="interpolated function"):
            inp = inputs.function.Function(source="interpolation")

        with pytest.raises(RuntimeError, match="function from XDMF"):
            inp = inputs.function.Function(source="xdmf")

        with pytest.raises(RuntimeError, match='no valid'):
            inp = inputs.function.Function(source="something-else")
