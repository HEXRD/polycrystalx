"""Material Input Module"""
from pathlib import Path

import numpy as np

from polycrystal.materials_database import MaterialsDataBase
from polycrystalx import inputs

from ..data import DATA_DIR

mdb = MaterialsDataBase(DATA_DIR / "materials.yaml")
PROCESS = "linear_elasticity"


def get_material_input(key):
    """Return a named material input list"""
    return MaterialInput(key).material_input


class MaterialInput:
    """Builds material input for polycrystalx

    Parameters:
    ----------
    key: str
       name of material
    """

    def __init__(self, key):
        self.key = key

    @property
    def material_input(self):
        return inputs.material.LinearElasticity(
            name=self.key,
            materials=self.materials
        )

    @property
    def materials(self):
        m = mdb.get_material(PROCESS, self.key)
        m.system = "MANDEL"
        return [m]
