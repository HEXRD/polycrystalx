"""Deformation Input Module"""

import numpy as np

from polycrystalx import inputs
from polycrystalx.inputs.tools import interpolate
from polycrystal.utils.tensor_data.mandel_system import MandelSystem


X_COMP, Y_COMP, Z_COMP = 0, 1, 2


def get_deformation_input(key, matl_inp, poly_inp):
    """Return a named deformation input"""
    return DeformationInput(key, matl_inp, poly_inp).deformation_input


class DeformationInput:
    """Builds deformation input for deformationx

    Parameters
    ----------
    key: 3-tuple
        tuple of (bcs, i, j), where bcs is a string describing type of
        boundary conditions, and i and j are integers between 1 and 3;
        bcs can be one of
            {"full", "zmax-traction", "zmax-traction-z", "zmax-traction-xy"}
    matl_inp: inputs.MaterialInput
        material input for this problem
    poly_inp: inputs.PolycrystalInput
        polycrystal input for this problem
    """

    def __init__(self, key, matl_inp, poly_inp):

        self.bcs, self.i, self.j = key
        self.matl_inp = matl_inp
        self.poly_inp = poly_inp
        self.set_boundary_conditions()

    @property
    def deformation_input(self):
        return inputs.deformation.LinearElasticity(
            name = self.name,
            force_density = self.body_force,
            displacement_bcs = self.displacement_bcs,
            traction_bcs = self.traction_bcs,
        )

    @property
    def name(self):
        return f"{self.bcs}-{self.i}-{self.j}"

    @property
    def body_force(self):
        return inputs.function.Function(
            source="constant",
            value=(0, 0, 0),
        )

    @property
    def A(self):
        """Deformation matrix A"""
        A = np.zeros((3, 3))
        A[self.i, self.j] = 1.0
        return A

    @property
    def orientation(self):
        """Crystal orientation

        Currently, this is simply the identity.
        """
        ### Fix this.
        return self.poly_inp.polycrystal.orientation_list[0]

    @property
    def crystal_stiffness(self):
        """Stiffness matrix in crystal frame"""
        return self.matl_inp.materials[0].stiffness

    @property
    def displacement_bcs(self):
        return self._displacement_bcs

    @property
    def traction_bcs(self):
        return self._traction_bcs

    def set_boundary_conditions(self):
        """Return displacement bcs and traction BCs"""
        A = self.A

        if self.bcs == "full":

            displacement_bcs = [
                inputs.deformation.DisplacementBC(
                    section = "boundary",
                    value = interpolate.linear_function(A),
                )
            ]
            traction_bcs = []

        elif self.bcs == "zmax-traction":

            displacement_bcs = [
                inputs.deformation.DisplacementBC(
                    section = "not-zmax",
                    value = interpolate.linear_function(A),
                )
            ]
            traction_bcs = [
                inputs.deformation.TractionBC(
                    section = "zmax",
                    value = interpolate.constant(self.traction((0, 0, 1))),
                )
            ]
        elif self.bcs == "zmax-traction-z":
            # Only z component of traction is applied.

            displacement_bcs = [
                inputs.deformation.DisplacementBC(
                    section = "not-zmax",
                    value = interpolate.linear_function(A),
                ),
                inputs.deformation.DisplacementBC(
                    section = "zmax",
                    value = interpolate.linear_function(A[X_COMP]),
                    component = X_COMP
                ),
                inputs.deformation.DisplacementBC(
                    section = "zmax",
                    value = interpolate.linear_function(A[Y_COMP]),
                    component = Y_COMP
                )
            ]
            trac_z = self.traction((0, 0, 1))[Z_COMP]
            traction_bcs = [
                inputs.deformation.TractionBC(
                    section = "zmax",
                    value = interpolate.constant((0, 0, trac_z)),
                )
            ]
        elif self.bcs == "zmax-traction-xy":
            # Tractions is x- and y- components on top surface.
            displacement_bcs = [
                inputs.deformation.DisplacementBC(
                    section = "not-zmax",
                    value = interpolate.linear_function(A),
                ),
                inputs.deformation.DisplacementBC(
                    section = "zmax",
                    value = interpolate.linear_function(A[Z_COMP]),
                    component = Z_COMP
                ),
            ]
            trac = self.traction((0, 0, 1))
            traction_bcs = [
                inputs.deformation.TractionBC(
                    section = "zmax",
                    value = interpolate.constant((trac[X_COMP], 0, 0)),
                ),
                inputs.deformation.TractionBC(
                    section = "zmax",
                    value = interpolate.constant((0, trac[Y_COMP], 0)),
                ),
            ]
        else:
            msg = f"BC key not found: {self.bcs}"
            raise RuntimeError(msg)

        self._displacement_bcs = displacement_bcs
        self._traction_bcs = traction_bcs

    def traction(self, n):
        """Traction on surface with normal n of def with grad u = A for matl

        n: 3-tuple or array(3)
           surface normal vector

        Returns
        -------
        t: array(3)
           traction vector
        """

        """The following section of code computes the stress in the sample frame:
        1. Convert sample strain (symm A) from sample to crystal
        2. Create the Mandel 6-vector for the crystal strain
        3. Apply the 6x6 stiffness to obtain the stress in the crystal frame
        4. Convert the stress 6-vector to a matrix, still in crystal frame
        5. Convert the stress matrix from crystal to sample
        6. Then you can apply the stress to the normal vector to get the traction.
        """
        C_c = self.crystal_stiffness
        A = self.A
        ori = self.orientation

        # import pdb; pdb.set_trace()
        eps_s = 0.5 * (A + A.T)
        eps_c = ori @ A @ ori.T
        eps_c6 = MandelSystem(eps_c).symm
        sig_c6 = C_c @ eps_c6.T
        sig_c = MandelSystem.from_parts(symm=sig_c6.T).matrices[0]
        sig = ori.T @ sig_c @ ori

        return sig @ n
