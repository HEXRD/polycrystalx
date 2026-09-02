"""Deformation Input Module"""

import numpy as np

from scipy.spatial.transform import Rotation

from polycrystalx import inputs
from polycrystalx.inputs.tools import interpolate

from ..data import EXTENTS, twist_by_state


def get_deformation_input(key):
    """Return a named deformation input"""
    return DeformationInput(key).deformation_input


class DeformationInput:
    """Builds deformation input for deformationx

    Parameters:
    ----------
    key: int
         load state, between 0 and 12
    """
    def __init__(self, key):
        self.key = key
        self.state = key

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
        return f"twist-{self.key:0d}"

    @property
    def body_force(self):
        return inputs.function.Function(
            source="constant",
            value=(0, 0, 0),
        )

    @property
    def displacement_bcs(self):
        """Displacement boundary conditions"""
        rad_per_mm = twist_by_state(self.state)
        axis = np.array((0, 1.0, 0))
        ymin, ymax = EXTENTS[1]
        ang_bot, ang_top = ymin * rad_per_mm, ymax * rad_per_mm

        # Displacement operators on bottom and top.
        IDENTITY = np.identity(3)
        u_bot = Rotation.from_rotvec(ang_bot * axis).as_matrix() - IDENTITY
        u_top = Rotation.from_rotvec(ang_top * axis).as_matrix() - IDENTITY

        displacement_bcs = [
            inputs.deformation.DisplacementBC(
                section = "ymin",
                value = interpolate.linear_function(u_bot),
            ),
            inputs.deformation.DisplacementBC(
                section = "ymax",
                value = interpolate.linear_function(u_top),
            ),
        ]

        return displacement_bcs

    @property
    def traction_bcs(self):
        return []
