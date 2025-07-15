"""Tests for loaders"""
import numpy as np
import pytest
from dolfinx import fem
import ufl

from polycrystalx import inputs
from polycrystalx.loaders.mesh import MeshLoader
from polycrystalx.loaders.function import FunctionLoader
from polycrystalx.loaders.deformation import (
    LinearElasticity, HeatTransfer
)


@pytest.fixture
def mesh_input():
    extents = np.array([
        [0, 1],
        [0, 2],
        [0, 3],
    ])
    extents = [[0, 1], [0, 2], [0, 3]]
    return inputs.mesh.Mesh(
        name="test-mesh",
        source="box",
        extents=extents,
        divisions=(2, 3, 5),
        celltype="tetrahedron",
    )


@pytest.fixture
def mesh_loader(mesh_input):
    return MeshLoader(mesh_input)


@pytest.fixture
def scalar_function_value():
    return 1.1


@pytest.fixture
def scalar_function(scalar_function_value):
    return inputs.function.Function(
        source="constant",
        value=scalar_function_value
    )


@pytest.fixture
def vector_function_value():
    return (1., 2, 3)


@pytest.fixture
def vector_function(vector_function_value):
    return inputs.function.Function(
        source="constant",
        value=vector_function_value
    )


@pytest.fixture
def tensor_function_value():
    return (1., 2, 3, 4, 5, 6, 7, 8, 9)


@pytest.fixture
def tensor_function(tensor_function_value):
    return inputs.function.Function(
        source="constant",
        value=tensor_function_value
    )


def test_function(
        mesh_loader,
        scalar_function, scalar_function_value,
        vector_function, vector_function_value
):
    """Test function loader"""
    msh = mesh_loader.mesh
    V = fem.functionspace(msh, ("P", 1))
    V3 = fem.functionspace(msh, ("P", 1, (3,)))

    # Check scalar constant
    f_c = FunctionLoader(scalar_function).load(V)
    assert np.all(f_c.x.array == scalar_function_value)

    # Check vector constant
    f_cr = FunctionLoader(vector_function).load(V3)
    arr = f_cr.x.array.reshape((len(f_cr.x.array) // 3, 3))
    assert np.all(arr == vector_function_value)


def test_mesh(mesh_loader):

    assert mesh_loader.tdim == 3
    assert mesh_loader.bdim == 2

    bd = mesh_loader.boundary_dict
    assert (
        "boundary" in bd
        and "xmin" in bd and "ymin" in bd and "zmin" in bd
        and "xmax" in bd and "ymax" in bd and "zmax" in bd
    )


class TestLinearElasticity:

    @pytest.fixture
    def defm_input(self, vector_function, tensor_function):
        defm_input = inputs.deformation.LinearElasticity(
            name="defm-input",
            force_density=vector_function,
            plastic_distortion=tensor_function,
            thermal_expansion=tensor_function,
        )
        return defm_input

    def test_force_density(self, mesh_loader, defm_input):

        ldr = LinearElasticity(defm_input)
        V3 = fem.functionspace(mesh_loader.mesh, ('P', 1, (3,)))

        assert isinstance(ldr.force_density(V3), fem.Function)

    def test_plastic_distortion(self, mesh_loader, defm_input):

        ldr = LinearElasticity(defm_input)
        T = fem.functionspace(mesh_loader.mesh, ('DG', 0, (3, 3)))

        assert isinstance(ldr.plastic_distortion(T), fem.Function)

    def test_thermal_expansion(self, mesh_loader, defm_input):

        ldr = LinearElasticity(defm_input)
        T = fem.functionspace(mesh_loader.mesh, ('DG', 0, (3, 3)))

        assert isinstance(ldr.thermal_expansion(T), fem.Function)


class TestHeatTransfer:

    @pytest.fixture
    def boundary_condition(self, scalar_function):
        return inputs.deformation.BoundaryCondition(
            section="xmin",
            value=lambda x: np.ones(x.shape[1]),
        )

    @pytest.fixture
    def defm_input(self, scalar_function, boundary_condition):
        return inputs.deformation.HeatTransfer(
            name="defm-input",
            body_heat=scalar_function,
            temperature_bcs=[boundary_condition],
            flux_bcs=[boundary_condition],
        )

    def test_body_heat(self, mesh_loader, defm_input):

        ldr = HeatTransfer(defm_input)
        V = fem.functionspace(mesh_loader.mesh, ('P', 1))

        assert isinstance(ldr.body_heat(V), fem.Function)

    def test_temperature_bcs(self, mesh_loader, defm_input):

        ldr = HeatTransfer(defm_input)
        V = fem.functionspace(mesh_loader.mesh, ('P', 1))
        bdict = mesh_loader.boundary_dict

        assert isinstance(ldr.temperature_bcs(V, bdict)[0], fem.DirichletBC)

    def test_flux_bcs(self, mesh_loader, defm_input):

        ldr = HeatTransfer(defm_input)
        V = fem.functionspace(mesh_loader.mesh, ('P', 1))
        bdict = mesh_loader.boundary_dict

        flux_bc0 = ldr.flux_bcs(V, bdict)[0]
        assert isinstance(flux_bc0.ds, ufl.Measure)
