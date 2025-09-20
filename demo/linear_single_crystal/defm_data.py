"""Deformation inputs"""
import numpy as np

from polycrystalx import inputs
from polycrystalx.inputs.tools import interpolate
from polycrystal.utils.tensor_data.mandel_system import MandelSystem


X_COMP, Y_COMP, Z_COMP = 0, 1, 2
zero = inputs.function.Function(
    source="constant",
    value=(0, 0, 0),
)


def get_input(key, matl):
    """Return deformation input

    PARAMETERS
    ----------
    key: 3-tuple
        tuple of (bcs, i, j), where bcs is a string describing type of
        boundary conditions, and i and j are integers between 1 and 3;
        bcs can be one of {"full", "trac-z", "trac-zx", trac-zy", trac-zz"}
    matl_name: str
        name of material

    RETURNS
    -------
    inputs.deformation.Elastic
         input specification for this deformation
    """
    bcs, i, j = key
    defm_name = f"{bcs}-{i}-{j}"

    A = np.zeros((3, 3))
    A[i-1, j-1] = 1.0
    ori = np.identity(3)

    if bcs == "full":

        displacement_bcs = [
            inputs.deformation.DisplacementBC(
                section = "boundary",
                value = interpolate.linear_function(A),
            )
        ]
        traction_bcs = []

    elif bcs == "zmax-traction":

        displacement_bcs = [
            inputs.deformation.DisplacementBC(
                section = "not-zmax",
                value = interpolate.linear_function(A),
            )
        ]
        traction_bcs = [
            inputs.deformation.TractionBC(
                section = "zmax",
                value = interpolate.constant(traction(matl, ori, A, (0, 0, 1))),
            )
        ]
    elif bcs == "zmax-traction-z":
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
        trac_z = traction(matl, ori, A, (0, 0, 1))[Z_COMP]
        traction_bcs = [
            inputs.deformation.TractionBC(
                section = "zmax",
                value = interpolate.constant((0, 0, trac_z)),
            )
        ]
    elif bcs == "zmax-traction-xy":
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
        trac = traction(matl, ori, A, (0, 0, 1))
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
        msg = f"BC key not found: {bcs}"
        raise RuntimeError(msg)

    defm_input = inputs.deformation.LinearElasticity(
        name = defm_name,
        force_density = zero,
        displacement_bcs = displacement_bcs,
        traction_bcs = traction_bcs,
    )

    return defm_input


def traction(matl, ori, A, n):
    """Traction on surface with normal n of def with grad u = A for matl"""
    matl.output_system = "MANDEL"

    """The following section of code computes the stress in the sample frame:
    1. Convert strain from sample to crystal
    2. Create the Mandel 6-vector for the crystal strain
    3. Apply the 6x6 stiffness to obtain the stress in the crystal frame
    4. Convert the stress 6-vector to a matrix, still in crystal frame
    5. Convert the stress matrix from crystal to sample
    6. Then you can apply the stress to the normal vector to get the traction.
    """
    C_c = matl.stiffness # crystal frame
    eps_s = A
    eps_c = ori @ A @ ori.T
    eps_c6 = MandelSystem(eps_c).symm
    sig_c6 = C_c @ eps_c6.T
    sig_c = MandelSystem.from_parts(symm=sig_c6.T).matrices[0]
    sig = ori.T @ A @ ori

    return sig @ n
