"""Tests for loaders"""
import numpy as np
import pytest
from dolfinx import fem

from polycrystalx import inputs
from polycrystalx.loaders.mesh import MeshLoader
from polycrystalx.loaders.function import FunctionLoader
from polycrystalx.loaders.deformation import LinearElasticity, HeatTransfer


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

def test_function(mesh_loader, scalar_function, scalar_function_value):
    """Test function loader"""
    msh = mesh_loader.mesh
    V = fem.functionspace(msh, ("P", 1))
    V2 = fem.functionspace(msh, ("P", 1, (2,)))

    # Check scalar constant
    f_c = FunctionLoader(scalar_function).load(V)
    assert np.all(f_c.x.array == scalar_function_value)

    # Check vector constant
    inp_c2 = inputs.function.Function(
        source="constant",
        value=(value := (2.0, 3.0))
    )
    f_c2 = FunctionLoader(inp_c2).load(V2)
    result = np.vstack((
        2 * np.ones(len(f_c.x.array)),
        3 * np.ones(len(f_c.x.array)),
    )).T
    assert np.all(f_c2.x.array == result.flatten())


def test_mesh(mesh_loader):

    assert mesh_loader.tdim == 3
    assert mesh_loader.bdim == 2

    bd = mesh_loader.boundary_dict
    assert (
        "boundary" in bd
        and "xmin" in bd and "ymin" in bd and "zmin" in bd
        and "xmax" in bd and "ymax" in bd and "zmax" in bd
    )


# @pytest.mark.usefixtures(mesh_loader)
class TestDeformationLoader:

    def test_linear_elasticity(self, mesh_loader):

        pd = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
        pdfun = inputs.function.Function(
            source = "constant",
            value = pd,
        )
        defm_input = inputs.deformation.LinearElasticity(
            name="test-le",
            plastic_distortion=pdfun
        )

        ldr = LinearElasticity(defm_input)
        T = fem.functionspace(mesh_loader.mesh, ('DG', 0, (3, 3)))

        assert ldr.plastic_distortion(T) is not None


    def test_heat_transfer(self, mesh_loader):

        bh_fun = inputs.function.Function(
            source = "constant",
            value = 2.0,
        )
        defm_input = inputs.deformation.HeatTransfer(
            name="test-ht",
            body_heat=bh_fun
        )

        ldr = HeatTransfer(defm_input)
        V = fem.functionspace(mesh_loader.mesh, ('P', 1))

        assert ldr.body_heat(V) is not None
