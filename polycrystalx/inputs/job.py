"""Input Job Module"""
from collections import namedtuple
import pathlib

_flds = [
    "suite", "process", "mesh_input", "material_input", "polycrystal_input",
    "deformation_input", "options_input"
]
_JobBase = namedtuple("_JobBase", _flds, defaults=[None])
del _flds


class Job(_JobBase):
    """Job Class

    PARAMETERS
    ----------
    suite: str
       name of the job suite
    process: name
       name of the process being modeled in this simulation
    mesh_input: inputs.mesh.Mesh
       input mesh specification
    material_input: inputs.material specification
        input specification for materials
    polycrystal_input: inputs.polcrystal.Polcrystal
       input specification for the polycrystal configuration
    deformation_input: inputs.deformation_input specification
       input deformation specification
    options_input: inputs.options.Option | None
       user options
    """

    @property
    def output_directory(self):
        """Name of output directory"""
        outbase = pathlib.Path("Outputs") / f"{self.suite}" / f"{self.process}"
        name = (
            f"{self.material_input.name}"
            f"-{self.mesh_input.name}"
            f"-{self.polycrystal_input.name}"
            f"-{self.deformation_input.name}"
        )
        if self.options_input is not None:
            name += f"-{self.options_input.name}"

        return outbase / name

    @property
    def log_file(self):
        """Name of log file"""
        return self.output_directory / "job.log"
