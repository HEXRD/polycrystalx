"""Material Input Module"""

import numpy as np

from polycrystalx import inputs

from polycrystal.elasticity.single_crystal import SingleCrystal


matl_dict = {
    "identity-iso": SingleCrystal.from_K_G(1/3, 1/2),
    "iso-21": SingleCrystal.from_K_G(2, 1),
    "cubic-211": SingleCrystal.from_K_Gd_Gs(2.0, 1.0, 1.0, system="MANDEL"),
    "cubic-121": SingleCrystal.from_K_Gd_Gs(1.0, 2.0, 1.0, system="MANDEL"),
    "cubic-112": SingleCrystal.from_K_Gd_Gs(1.0, 1.0, 2.0, system="MANDEL"),
    "cubic-321": SingleCrystal.from_K_Gd_Gs(3.0, 2.0, 1.0, system="MANDEL"),
}


def get_material_input(key):
    """Return a named material input list"""
    return MaterialInput(key).material_input


class MaterialInput:
    """Builds material input for polycrystalx

    Parameters:
    ----------
    key: str (or list of strings)
       name of material in `material_dict`
    """

    OUTPUT_SYSTEM = "MANDEL"

    def __init__(self, key):
        self.name = key
        self.material = matl_dict[key]
        self.material.system = self.OUTPUT_SYSTEM

    @property
    def material_input(self):
        return inputs.material.LinearElasticity(
            name=self.name, materials=[self.material]
        )
