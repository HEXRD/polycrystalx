"""Deformation inputs"""
import numpy as np

from polycrystalx import inputs
from polycrystalx.inputs.tools import interpolate
from polycrystal.utils.tensordata import TensorData


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
    C = matl.sample_stiffness(ori)
    eps = TensorData(A.reshape(1, 3, 3))
    sig6 = C @ eps.symm6[0]
    sig = TensorData.from_parts(symm6=sig6.reshape((1,6))).matrices[0]
    return sig @ n
