"""Thermal Conductivity Process"""
import numpy as np

from ..loaders import mesh
from ..loaders import material
from ..loaders import polycrystal
from ..loaders import deformation


class ThermalConductivity:
    """Process of Thermal Conductivity

    Parameters
    ----------
    job: inputs.job.Job
       user inputs for this job

    """
    name = "thermal-conductivity"

    def __init__(self, job):
        self.loader = _Loader(job)
        self.mpirank = self.loader.mesh.comm.rank

    def run(self):
        pass


class _Loader:

    def __init__(self, job):

        self.job = job

        # Material Data
        self.material_data = material.LinearElasticity(input_mod.material_input)

        # Mesh Data and Function Spaces
        self.mesh_data = mesh.MeshLoader(job.mesh_input)

        self.problem = ThermalConductivityProblem(self.mesh)
        self.V = self.problem.V
        self.T = self.problem.T

        # Microstructure Data
        self.polycrystal_data = polycrystal.Polycrystal(
            input_mod.polycrystal_input
        )
        if self.polycrystal_data.use_meshtags:
            self.cell_tags = self.mesh_data.cell_tags
            print("using cell tags from gmsh input file")
        else:
            self.cell_tags = self.polycrystal_data.grain_cell_tags(
                self.mesh
            )
        self.grain_cells = self.polycrystal_data.grain_cell_dict(
            self.cell_tags
        )
        self.orientation_fld = self.polycrystal_data.orientation_field(
            self.T, self.grain_cells
        )
        self._stiffness_fld = self._make_stiffness_fld()

        # Deformation Data
        self.deformation_data = deformation.LinearElasticity(
            input_mod.deformation_input
        )
        self.force_density = self.deformation_data.force_density(self.V)

        self.plastic_distortion = self.deformation_data.plastic_distortion(
            self.T
        )

    @property
    def mesh(self):
        return self.mesh_data.mesh

    @property
    def stiffness_fld(self):
        return self._stiffness_fld

    def _make_stiffness_fld(self):
        stf_fld= fem.Function(self.T6)
        ms = self.polycrystal_data.polycrystal
        for gi in range(ms.num_grains):
            phase = int(ms.phase(np.array([gi])))
            matl = self.material_data.materials[phase]
            cells = self.grain_cells[gi]
            stf = matl.stiffness
            stf_fun = lambda x: self._stiffness(stf, x)
            stf_fld.interpolate(stf_fun, cells)
        return stf_fld

    def _stiffness(self, stf, x):
        return np.tile(stf.reshape(36,1), x.shape[1])

    @property
    def boundary_dict(self):
        return self.mesh_data.boundary_dict

    @property
    def displacement_bcs(self):
        return self.deformation_data.displacement_bcs(
            self.V, self.boundary_dict
        )

    @property
    def traction_bcs(self):
        return self.deformation_data.traction_bcs(self.V, self.boundary_dict)
