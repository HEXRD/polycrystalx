"""Elastic Process"""
import time

import xml.etree.ElementTree as ET

import numpy as np
from dolfinx import fem, log, io
from dolfinx.common import Timer
from dolfinx.fem import assemble_scalar, form
from dolfinx.mesh import locate_entities
from dolfinx.fem.petsc import LinearProblem
import ufl
from mpi4py import MPI

from ..loaders import mesh
from ..loaders import material
from ..loaders import polycrystal
from ..loaders import deformation
from ..forms.common import sigs_3x3, grain_volume, grain_integral
from ..forms.linear_elasticity import (
    LinearElasticity as LinearElasticityProblem
)
from ..utils import grain_volumes, grain_integrals


class LinearElasticity:
    """Linear elastic process

    Parameters
    ----------
    user_input: Job
       user inputs for this job

    """
    name = "linear-elasticity"

    def __init__(self, user_input):
        self.loader = _Loader(user_input)
        self.mpirank = self.loader.mesh.comm.rank
        print("My rank is ", self.mpirank)

    def run(self):
        """Run the problem"""
        print("running problem", flush=True)
        print("creating loader", flush=True)
        ldr = self.loader

        print("evaluating coefficients", flush=True)

        coeffs = ldr.problem.coefficients
        coeffs.orientation.x.array[:] = ldr.orientation_fld.x.array
        cstiff = ldr.stiffness_fld
        coeffs.stiffness.x.array[:] = cstiff.x.array
        coeffs.body_force.x.array[:] = ldr.force_density.x.array
        for tbc in ldr.traction_bcs:
            coeffs.tractions.append(tbc)
        a, L = ldr.problem.forms

        print("making displacement bcs", flush=True)
        mybcs = ldr.displacement_bcs
        print("setting up linear problem", flush=True)
        problem = LinearProblem(
            a, L, bcs=mybcs,
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
        )

        with Timer() as t:
            print("starting linear solver", flush=True)
            uh = problem.solve()
            print(f"linear solver time: {t.elapsed()}")
            print("postprocessing ...")
            self.postprocess(uh, ldr)

        return uh

    def postprocess(self, uh, ldr):
        """Compute strains and stresses and write output"""
        # Write the XDMF File.

        uh.name = "displacement"
        ldr.cell_tags.name = "grain-ids"

        strain_form = ufl.sym(ufl.grad(uh))
        strain_expr = fem.Expression(
            strain_form, ldr.T.element.interpolation_points()
        )
        strain = fem.Function(ldr.T, name="strain")
        strain.interpolate(strain_expr)

        stress_form = sigs_3x3(uh, ldr.stiffness_fld, ldr.orientation_fld)
        stress_expr = fem.Expression(
            stress_form, ldr.T.element.interpolation_points()
        )
        stress = fem.Function(ldr.T, name="stress")
        stress.interpolate(stress_expr)

        with io.XDMFFile(ldr.mesh.comm, "output.xdmf", "w") as file:
            file.write_mesh(ldr.mesh)
            file.write_meshtags(ldr.cell_tags, ldr.mesh.geometry)
            file.write_function(uh)
            file.write_function(strain)
            file.write_function(stress)

        # Compute grain volumes.

        print("finding grain volumes")
        with Timer() as t:
            gv_form, indic = grain_volume(ldr.mesh)
            g_volumes = grain_volumes(
                ldr.mesh.comm, gv_form, indic, ldr.grain_cells
            )
            elapsed = t.elapsed()
        if self.mpirank == 0:
            print(f"total volume: {np.sum(g_volumes)}", flush=True)
            print(f"time for volume calculation: {elapsed}")

        # Compute grain integrals.

        eps_int = np.zeros((num_grains := len(g_volumes), 6))

        V = fem.functionspace(ldr.mesh, ("DG", 0))
        gi_form, indicator, func = grain_integral(ldr.mesh, V)

        indmap = ldr.mesh.topology.index_map(ldr.mesh.topology.dim)
        allcells = np.arange(indmap.size_local).astype(np.int32)
        eps_fun = fem.Function(V)

        with Timer() as t:
            eps_expr = fem.Expression(
                strain[0, 0], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 0] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            eps_expr = fem.Expression(
                strain[1, 1], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 1] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            eps_expr = fem.Expression(
                strain[2, 2], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 2] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            eps_expr = fem.Expression(
                strain[1, 2], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 3] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            eps_expr = fem.Expression(
                strain[0, 2], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 4] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            eps_expr = fem.Expression(
                strain[0, 1], V.element.interpolation_points()
            )
            eps_fun.interpolate(eps_expr, allcells)
            eps_int[:, 5] = grain_integrals(
                ldr.mesh.comm, gi_form, indicator, func, ldr.grain_cells,
                eps_fun
            )
            elapsed = t.elapsed()
        if self.mpirank == 0:
            print(f"time for grain integrals calculation: {elapsed}")

        if self.mpirank == 0:
            eps_avg = np.zeros_like(eps_int)
            nz = g_volumes > 0.
            nnz = np.count_nonzero(g_volumes > 0)
            eps_avg[nz] = eps_int[nz]/g_volumes[nz].reshape(nnz, 1)
            np.savez("grain-averages.npz", volume=g_volumes, strain=eps_avg)
            #
            self.write_xdmf()

    def write_xdmf(self, output="output.xdmf", paraview="paraview.xdmf"):
        """This puts all the data into the same grid

        This writes two XDMF files--the usual output file written using the
        fenicsx writer and a second XDMF written specifically for viewing
        in paraview.

        Parameters
        ----------
        output: str or Path, default = "output.xdmf"
            name of output XDMF file
        paraview: str or Path, default = "paraview.xdmf"
            name of XMDF file for paraview
        """

        ATTR = "Attribute"
        NAME = "Name"

        # Start with displacement, which also includes meshtags.

        tree = ET.parse(output)
        root = tree.getroot()
        domain = root[0]
        meshgrid = domain[0]

        mtags = domain[1].find(ATTR)
        mtags.attrib[NAME] = "grain ID"
        meshgrid.append(mtags)

        disp = domain[2][0].find(ATTR)
        disp.attrib[NAME] = "displacement"
        meshgrid.append(disp)

        strn = domain[3][0].find(ATTR)
        strn.attrib[NAME] = "strain"
        meshgrid.append(strn)

        stress = domain[4][0].find(ATTR)
        stress.attrib[NAME] = "stress"
        meshgrid.append(stress)

        domain.remove(domain[4])
        domain.remove(domain[3])
        domain.remove(domain[2])
        domain.remove(domain[1])

        # Write the modified tree.

        tree.write(paraview)


class _Loader:

    def __init__(self, input_mod):

        self.input_module = input_mod

        # Material Data
        self.material_data = material.LinearElasticity(input_mod.material_input)

        # Mesh Data and Function Spaces
        self.mesh_data = mesh.MeshLoader(input_mod.mesh_input)

        self.problem = LinearElasticityProblem(self.mesh)
        self.V = self.problem.V
        self.T = self.problem.T
        self.T6 = self.problem.T6

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
