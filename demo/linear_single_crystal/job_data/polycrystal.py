"""Polycrystal Input Module"""

import numpy as np

from polycrystalx import inputs
from polycrystal.microstructure.single_crystal import (
    SingleCrystal as SingleCrystalMicro
)
from polycrystal.orientations.quaternions import from_exp, to_rmats


def get_polycrystal_input(key):
    """Return a named polycrystal input"""
    return PolycrystalInput(key).polycrystal_input


class PolycrystalInput:
    """Builds polycrystal input for polycrystalx

    Parameters:
    ----------
    key: (ax_num, ang_deg)
       coordinate axis number (0, 1, or 2) and angle (int) in degrees
    """

    def __init__(self, key):
        self.ax_num, self.ang_deg = key

    @property
    def polycrystal_input(self):
        return inputs.polycrystal.Polycrystal(
            name=self.name,
            polycrystal=self.microstructure
        )

    @property
    def name(self):
        return "-".join(("single-crystal", "xyz"[self.ax_num] + str(self.ang_deg)))

    @property
    def microstructure(self):
        axis = np.zeros(3)
        axis[self.ax_num] = 1.0
        ang_rad = np.radians(self.ang_deg)
        orientation = self.rmatrix(ang_rad, axis)

        return SingleCrystalMicro(orientation)


    @staticmethod
    def rmatrix(angle, axis):
        a = np.atleast_2d(np.array(axis))
        n = a / np.linalg.norm(a,1)
        return to_rmats(from_exp(angle * n))
