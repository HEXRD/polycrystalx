"""Heat Transfer"""
import pathlib
import xml.etree.ElementTree as ET

import numpy as np
from dolfinx import fem, log, io
from dolfinx.common import Timer
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl

from ..loaders import mesh
from ..loaders import material
from ..loaders import polycrystal
from ..loaders import deformation
from ..loaders import options

from ..forms.heat_transfer import HeatTransferProblem
from ..utils import grain_integrals


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

    def run(self, outdir):
        """Run the problem

        Parameters
        ----------
        outdir: pathlib.Path
            output directory for all result files
        """
        ldr = self.loader

        # Fill in the forms.

        coeffs = ldr.problem.coefficients
        coeffs.orientation.x.array[:] = ldr.orientation_fld.x.array
        coeffs.stiffness.x.array[:] = ldr.stiffness_fld.x.array
        coeffs.body_heat.x.array[:] = ldr.body_heat.x.array
        for fbc in ldr.flux_bcs:
            coeffs.fluxes.append(fbc)
        a, L = ldr.problem.forms

        # Make temperature BCs.

        default_petsc_options = {
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
            petsc_options=ldr.options.solver_options
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
        self.postprocess(uh, ldr, outdir)

    def postprocess(self, uh, ldr, outdir):
        """Write primary variables and compute grain averaged values

        Parameters
        ----------
        uh: dolfinx.fem.Function
            temperature solution
        ldr: _Loader
            loader holding mesh and field data
        outdir: pathlib.Path
            output directory for all result files
        """
        outdir = pathlib.Path(outdir)

        uh.name = "temperature"
        ldr.cell_tags.name = "grain-ids"

        # Compute flux field first.
        flux_form = ldr.problem.flux(uh)
        flux_expr = fem.Expression(
            flux_form, ldr.V3.element.interpolation_points()
        )
        flux_fun = fem.Function(ldr.V3, name="flux")
        flux_fun.interpolate(flux_expr)

        with io.XDMFFile(ldr.mesh.comm, str(outdir / "output.xdmf"), "w") as file:
            file.write_mesh(ldr.mesh)
            file.write_meshtags(ldr.cell_tags, ldr.mesh.geometry)
            file.write_function(uh)
            file.write_function(flux_fun)

        if self.mpirank == 0:
            print("Evaluating grain volumes and integrals")

        with Timer() as t:
            V = fem.functionspace(ldr.mesh, ("DG", 0))
            one = fem.Function(V)
            one.interpolate(lambda x: np.full_like(x[0], 1.0))

            g_volumes = grain_integrals(one, ldr.grain_cells)
            temp_ints = grain_integrals(uh, ldr.grain_cells)
            flux_ints = grain_integrals(flux_fun, ldr.grain_cells)

            elapsed = t.elapsed()

        if self.mpirank == 0:
            print(f"time for grain integrals calculation: {elapsed}")

            # Now find grain averages from integrals.
            nz = g_volumes > 0.
            nnz = np.count_nonzero(g_volumes > 0)
            gvnnz = g_volumes[nz].reshape(nnz)

            temp_avg = np.zeros_like(temp_ints)
            temp_avg[nz] = temp_ints[nz] / gvnnz

            flux_avg = np.zeros_like(flux_ints)
            flux_avg[nz] = flux_ints[nz] / gvnnz.reshape(nnz, 1)
            np.savez(outdir / "grain-averages.npz", volume=g_volumes,
                     temperature=temp_avg, flux=flux_avg)

            self.write_xdmf(outdir)

    def write_xdmf(self, outdir, output="output.xdmf", paraview="paraview.xdmf"):
        """This puts all the data into the same grid

        This writes two XDMF files--the usual output file written using the
        fenicsx writer and a second XDMF written specifically for viewing
        in paraview.

        Parameters
        ----------
        outdir: pathlib.Path
            directory containing output files
        output: str, default = "output.xdmf"
            name of output XDMF file (relative to outdir)
        paraview: str, default = "paraview.xdmf"
            name of XDMF file for paraview (relative to outdir)
        """
        outdir = pathlib.Path(outdir)

        ATTR = "Attribute"
        NAME = "Name"

        # Start with the mesh tags, which also includes the mesh.

        tree = ET.parse(outdir / output)
        root = tree.getroot()
        domain = root[0]
        meshgrid = domain[0]

        mtags = domain[1].find(ATTR)
        mtags.attrib[NAME] = "grain-ids"
        meshgrid.append(mtags)

        temperature = domain[2][0].find(ATTR)
        temperature.attrib[NAME] = "temperature"
        meshgrid.append(temperature)

        flux = domain[3][0].find(ATTR)
        flux.attrib[NAME] = "flux"
        meshgrid.append(flux)

        domain.remove(domain[3])
        domain.remove(domain[2])
        domain.remove(domain[1])

        # Write the modified tree.

        tree.write(outdir / paraview)


class _Loader:

    def __init__(self, job):

        self.job = job

        # Material Data
        self.material_data = material.HeatTransfer(job.material_input)

        # Mesh Data and Function Spaces
        self.mesh_data = mesh.MeshLoader(job.mesh_input)

        self.problem = HeatTransferProblem(self.mesh)
        self.V = self.problem.V
        self.V3 = self.problem.V3
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

        # Options data.
        self.options = options.OptionsLoader(job.options_input).options


    @property
    def mesh(self):
        return self.mesh_data.mesh

    @property
    def stiffness_fld(self):
        return self._stiffness_fld

    def _make_stiffness_fld(self):
        stf_fld = fem.Function(self.T)
        ms = self.polycrystal_data.polycrystal
        for gi in range(ms.num_grains):
            phase = int(ms.phase(np.array([gi]))[0])
            matl = self.material_data.materials[phase]
            cells = self.grain_cells[gi]
            stf = matl.conductivity
            stf_fld.interpolate(
                lambda x: np.tile(stf.reshape(9, 1), x.shape[1]), cells
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
