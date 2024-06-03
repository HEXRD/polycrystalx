"""Slip system hardness fields"""
import xml.etree.ElementTree as ET

import numpy as np

from dolfinx import fem, mesh, io

from ..loaders import mesh
from ..loaders import material
from ..loaders import polycrystal
from ..loaders import deformation
from ..loaders import function


class Slip:
    """Slip Modeling

    Parameters
    ----------
    user_input: Job
       user input for this job
    """
    name = "slip"

    def __init__(self, user_input):
        self.loader = _Loader(user_input)
        self.mpirank = self.loader.mesh.comm.rank
        print("My rank is ", self.mpirank)

    def run(self):
        """Run the integration"""
        nsteps = self.loader.nsteps
        dt = self.loader.dt
        stepscale = 1/nsteps

        s0_arr = self.loader.s0.x.array
        st_arr = s0_arr.copy()
        cstress_arr = self.loader.crystal_stress_array(self.loader.stress_0)
        dcstress = (
            self.loader.crystal_stress_array(self.loader.stress_t)
            - cstress_arr
        )
        for step in range(self.loader.nsteps):
            print(f"step: {step}", flush=True)
            sd_arr= self.state_derviative(st_arr, cstress_arr)
            st_arr += dt * sd_arr.flatten()
            cstress_arr += stepscale * dcstress

        sT = fem.Function(self.loader.Vstate, name="sT")
        sT.x.array[:] = st_arr
        self.postprocess(self.loader, sT)
        print("done")

    def state_derviative(self, s, cstress):
        """Evaluate state variable derivative

        PARAMETERS
        ----------
        s: array (n, nsv)
            state variable value
        cstress: array (n, 3, 3)
            crystal stresses

        RETURNS
        array (n, nsv)
            values of state variable derivatives
        """
        ntot = len(s)
        nsv = self.loader.num_statevar
        ms = self.loader.polycrystal_data.polycrystal
        if ms.num_phases == 1:
            matl = self.loader.material_data.materials[0]
            sdata = matl.get(
                cstress, s.reshape(ntot//nsv, nsv), state_derivative=True
            )
        else:
            raise NotImplementedError("multiphase not yet implemented for slip")

        # Review this:
        #for gi in range(ms.num_grains):
        #    phase = int(ms.phase(np.array([gi])))
        #    matl = self.loader.material_data.materials[phase]
        #    cells = self.loader.grain_cells[gi]
        #    sdata = matl.get(cstress[cells], s[cells], state_derivative=True)

        return sdata.state_derivative

    def postprocess(self, ldr, sdot):
        """Write output files"""

        with io.XDMFFile(ldr.mesh.comm, "output.xdmf", "w") as file:
            file.write_mesh(ldr.mesh)
            file.write_meshtags(ldr.cell_tags)
            file.write_function(self.loader.s0)
            file.write_function(sdot)

        if self.mpirank == 0:
            self.write_xdmf()

    def write_xdmf(self, output="output.xdmf", paraview="paraview.xdmf"):
        """This puts all the data into the same grid"""

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

        s0 = domain[2][0].find(ATTR)
        s0.attrib[NAME] = "state-0"
        meshgrid.append(s0)

        sT = domain[3][0].find(ATTR)
        sT.attrib[NAME] = "state-T"
        meshgrid.append(sT)

        domain.remove(domain[3])
        domain.remove(domain[2])
        domain.remove(domain[1])

        # Write the modified tree.

        tree.write(paraview)


class _Loader:

    def __init__(self, input_mod):

        self.input_module = input_mod

        # Material Data
        self.material_data = material.Slip(input_mod.material_input)

        # Mesh Data and Function Spaces
        self.mesh_data = mesh.MeshLoader(input_mod.mesh_input)

        self.Tstress = fem.TensorFunctionSpace(self.mesh, ("DG", 0))
        self.T = fem.TensorFunctionSpace(self.mesh, ("DG", 0))

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

        # Deformation Data
        self.deformation_data = input_mod.deformation_input

    @property
    def mesh(self):
        """The mesh"""
        return self.mesh_data.mesh

    @property
    def dim(self):
        """Mesh dimension"""
        return self.mesh.tdim

    @property
    def stress_0(self):
        """Stress at initial time"""
        stress_spec = self.deformation_data.stress_0
        return function.FunctionLoader(stress_spec).load(self.Tstress)

    @property
    def stress_t(self):
        """Stress at target time"""
        stress_spec = self.deformation_data.stress_t
        return function.FunctionLoader(stress_spec).load(self.Tstress)

    def crystal_stress_array(self, sample_stress):
        """Stress array in crystal frame"""
        n9 = len(sample_stress.x.array)
        n = n9 // 9
        shp = (n, 3, 3)
        oria = self.orientation_fld.x.array.reshape(shp)
        siga = sample_stress.x.array.reshape(shp)
        # This rotates the stress arrays into the crystal frame.
        csig_a = np.einsum(
            'lij,ljk->lik',
            np.einsum('lij,ljk->lik', oria.transpose(0, 2, 1), siga), oria
        )
        return csig_a

    @property
    def num_statevar(self):
        """Number of state variables"""
        m0 = self.material_data.materials[0]
        return m0.num_statevar

    @property
    def Vstate(self):
        """Fucntion space for state variable"""
        if not hasattr(self, "_Vstate"):
            nsv = self.num_statevar
            if nsv == 1:
                self._Vstate = fem.FunctionSpace(self.mesh, ("DG", 0))
            else:
                self._Vstate = fem.VectorFunctionSpace(
                    self.mesh, ("DG", 0), dim=nsv
                )
        return self._Vstate

    @property
    def s0(self):
        """Initial state (dolfinx Function)"""
        s0_spec = self.deformation_data.s0
        return function.FunctionLoader(s0_spec).load(self.Vstate, name="s0")

    @property
    def dt(self):
        """Time increment"""
        return self.deformation_data.dt

    @property
    def nsteps(self):
        """Number of time steps"""
        return self.deformation_data.nsteps
