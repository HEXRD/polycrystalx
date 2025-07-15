"""Deformation input"""
import numpy as np

from polycrystalx import inputs
from polycrystalx.inputs.tools import interpolate


zero = inputs.function.Function(
    source="constant",
    value=(0, 0, 0),
)


def get_input(key, matl):
    """Return deformation inputs

    key: str
       key can be:
        "zero" for a zero displacement boundary conditions,
        "match" for boundary conditions that match thermal strains,
        "linear" for a linear temperature difference
    """
    expansion, displacement_bcs = _get_exp(key, matl)

    defm_input = inputs.deformation.LinearElasticity(
        name = key,
        force_density = zero,
        displacement_bcs = displacement_bcs,
        thermal_expansion = expansion,
        traction_bcs = [],
    )


    return defm_input


def _get_exp(defm_key, matl):

    linear_cte = matl.cte
    constant_diff = 300 * linear_cte

    def linear_diff(x):
        fac = 200 + x[0, :] * 50
        tmp = fac.reshape(-1, 1) * linear_cte.flatten()
        print("tmp: ", tmp.shape)
        return tmp.T


    constant_expansion = inputs.function.Function(
        source="constant",
        value=constant_diff.flatten(),
    )
    linear_expansion = inputs.function.Function(
        source="interpolation",
        function=linear_diff,
    )

    zero_bcs = [
        inputs.deformation.DisplacementBC(
            section = "boundary",
            value = interpolate.linear_function(np.zeros((3, 3))),
        )
    ]
    match_bcs = [
        inputs.deformation.DisplacementBC(
            section = "boundary",
            value = interpolate.linear_function(constant_diff),
        )
    ]

    if defm_key == "zero":

        texp = constant_expansion
        displacement_bcs = zero_bcs

    elif defm_key == "match":

        texp = constant_expansion
        displacement_bcs = match_bcs

    elif defm_key == "linear":

        texp = linear_expansion
        displacement_bcs = zero_bcs

    return texp, displacement_bcs
