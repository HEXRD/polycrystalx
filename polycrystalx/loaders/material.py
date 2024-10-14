"""Material input loader package

So far, all the material inputs are instances of certain material classes,
so the base class is just a list of instances.
"""
form polycrystal.thermal.single_crystal import (
    SingleCrystal as ThermalSingleCrystal
)


class MaterialList:

    def __init__(self, job):
        self.materials = job.materials
        self.check()

    def check(self):
        """Check that material instances are appropriate"""
        pass


class LinearElasticity:
    """Loader for LinearElasticity

    Parameters
    ----------
    userinput: module
       user input module
    """
    def __init__(self, userinput):
        self.materials = userinput.materials


class ThermalConductivity(MaterialList):

    def check(self):
        for m in self.materials:
            if not insinstance(m, ThermalSingleCrystal):
                raise ValueError("material is not a themal SingleCrystal")
