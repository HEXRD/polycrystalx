"""Heat Transfer"""
import numpy as np
from dolfinx import fem, log, io
from dolfinx.common import Timer
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

from ..loaders import mesh
from ..loaders import material
from ..loaders import polycrystal
from ..loaders import deformation

from ..forms.heat_transfer import HeatTransferProblem


class HeatTransfer:
    """Heat Transfer Process

    Parameters
    ----------
    job: inputs.job.Job
       user inputs for this job
    """
    name = "heat-transfer"

    def __init__(self, job):
        self.loader = _Loader(job)
        self.mpirank = self.loader.mesh.comm.rank

    def run(self):
        ldr = self.loader

        # Fill in the forms.

        coeffs = ldr.problem.coefficients
        coeffs.orientation.x.array[:] = ldr.orientation_fld.x.array
        coeffs.stiffness.x.array[:] = ldr.stiffness_fld.x.array
        coeffs.body_heat.x.array[:] = ldr.body_heat.x.array
        for fbc in ldr.flux_bcs:
            coeffs.tractions.append(fbc)
        a, L = ldr.problem.forms

        # Make temperature BCs.

        default_petsc_options={
            "ksp_type": "cg",
            "ksp_rtol": 1e-6,
            "ksp_atol": 1e-10,
            "ksp_max_it": 5000,
            "pc_type": "jacobi",
        }

        # Set up the linear problem and solve.

        mybcs = ldr.temperature_bcs
        linprob = LinearProblem(
            a, L, bcs=mybcs,
            petsc_options=default_petsc_options
        )

        with Timer() as t:
            print("starting linear solver", flush=True)
            uh = linprob.solve()
            print(f"linear solver time: {t.elapsed()}")

        solver = linprob.solver
        if solver.is_converged:
            print(f"solver converged: iterations = {solver.its}")
        else:
            msg = f"solver diverged: iterations = {solver.its}"

        print("postprocessing ...")
        with io.XDMFFile(ldr.mesh.comm, "output.xdmf", "w") as file:
            file.write_mesh(ldr.mesh)
            file.write_meshtags(ldr.cell_tags, ldr.mesh.geometry)
            file.write_function(uh)


class _Loader:

    def __init__(self, job):

        self.job = job

        # Material Data
        self.material_data = material.HeatTransfer(job.material_input)

        # Mesh Data and Function Spaces
        self.mesh_data = mesh.MeshLoader(job.mesh_input)

        self.problem = HeatTransferProblem(self.mesh)
        self.V = self.problem.V
        self.T = self.problem.T

        # Microstructure Data
        self.polycrystal_data = polycrystal.Polycrystal(
            job.polycrystal_input
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
        self.deformation_data = deformation.HeatTransfer(
            job.deformation_input
        )
        self.body_heat = self.deformation_data.body_heat(self.V)

    @property
    def mesh(self):
        return self.mesh_data.mesh

    @property
    def stiffness_fld(self):
        return self._stiffness_fld

    def _make_stiffness_fld(self):
        stf_fld= fem.Function(self.T)
        ms = self.polycrystal_data.polycrystal
        for gi in range(ms.num_grains):
            phase = int(ms.phase(np.array([gi])))
            matl = self.material_data.materials[phase]
            cells = self.grain_cells[gi]
            stf = matl.conductivity
            stf_fld.interpolate(
                lambda x: np.tile(stf.reshape(9,1), x.shape[1]), cells
            )
        return stf_fld

    @property
    def boundary_dict(self):
        return self.mesh_data.boundary_dict

    @property
    def temperature_bcs(self):
        return self.deformation_data.temperature_bcs(
            self.V, self.boundary_dict
        )

    @property
    def flux_bcs(self):
        return self.deformation_data.flux_bcs(self.V, self.boundary_dict)
