"""Polycrystal Input Module"""

import numpy as np

from polycrystalx import inputs
from polycrystal.microstructure.analytic import Analytic

from ..data import DATA_DIR


def get_microstructure_input(key):
    """Return a named polycrystal input"""
    return PolycrystalInput(key).polycrystal_input


class PolycrystalInput:
    """Builds polycrystal input for polycrystalx

    Parameters:
    ----------
    key: str
         name of the microstructure (doesn't actually do anything yet)
    """

    def __init__(self, key):
        self.key = key

    @property
    def polycrystal_input(self):
        return inputs.polycrystal.Polycrystal(
            name=self.name,
            polycrystal=self.microstructure,
            use_meshtags=True,
        )

    @property
    def name(self):
        return self.key

    @property
    def orientations(self):
        # Note that neper orientations start with 1, so we add a dummy orientation at 0.
        ori = np.load(DATA_DIR / "orientations.npy")
        return np.vstack((np.identity(3).reshape(1, 3, 3), ori))

    @property
    def microstructure(self):
        # This is just a dummy microstructure containing the orientations.
        return Analytic(None, self.orientations)
