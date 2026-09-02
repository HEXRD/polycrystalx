# Chess Template

This is an example template to illustrate how to use the `polycrystalx` package on CHESS far field HEDM results. For the example, we use data from a torsion experiment at APS in 2019.  Measurements were taken at zero load and at 12 load steps. The full data set has over 1,200 grains, but we only use about 100 to keep the example simpler. The material used was LSHR, a nickel alloy. The example runs elastic simulations on any or all load states.

The `jobs` directory is a python package that contains the infrastructure to set up the inputs.  The main input is in the `jobs/__init__.py` file. You can set up suites of jobs in the `jobs/batch.py` file.   The core template (no just for CHESS data) provides the layout of the `jobs` package.

## Running

To run the simulation, use the `pxx_job` command with `mpirun`.

```
mpirun -n 2 pxx_job jobs
```
This uses processes and loads the `jobs` module.  Here is the `__init__.py` file.

```
"""Inputs Module"""
from . import batch

# Material key is an entry in the materials database (in jobs/data).
matl_key = "lshr_660C"

# For now, there is only one microstructure (from state 0).
poly_key = "ms"

# These are base names of neper meshes already generated.
mesh_key = "vor-025"
mesh_key = "vor-050"

# This is a load state, between 0 and 12 inclusive.
defm_key = 12

jobkey = (matl_key, poly_key, mesh_key, defm_key)
job = batch.get_job(jobkey)
```

Each job is defined by four inputs: the material, the microstructure, the mesh and the deformation.  Here the material is `lshr_660C`. The microstructure is a Voronoi tessellation based on the centroids and orientations in the unlaoded state.  There are two meshes: one with a resolution of 50 microns and one with 25. There is a command line script to generate meshes of any resolution.  Finally, the deformation is torsion with a twist angle depending on the load state.

When the job runs, it writes output files into the `Outputs` directory with a path and a name based on the inputs. It writes the finite element to an HDF5 file with a standard XDMF header file (`output.XDMF`) and another header for viewing with paraview (`paraview.XDMF`). It also writes a numpy grain data file (`grain-averages.npz`)  with grain volumes, grain averaged strain and stress tensors.

## Core Template

The core template provides basic constructors for the four polycrystalx inputs from the input keys. These are in the `job_inputs` subpackage. Each input has its own module. They are `material.py`, `mesh.py`, `microstructure.py`, and  `deformation.py`.	Let us consider each one.

### Material Input
The polycrystalx material input looks like:

```
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

    ...
```
It requires a name and a list of elastic `SingleCrystal` instances from `polycrystal`.  In this example, there is a YAML material database file in the `jobs/data` directory.  So the material name is the key, and the database entry is used to generate the instance.

### Mesh Input
For the `polycrystalx` mesh input, here is the beginning of the class that builds it:
```
class MeshInput:
    """Builds mesh input for meshx

    Parameters:
    ----------
    key: str
         basename of the mesh file
    """
    def __init__(self, key):
        self.key = key

    @property
    def mesh_input(self):
        return inputs.mesh.Mesh(
            name=self.name,
            source="gmsh",
            extents=EXTENTS,
            celltype="tetrahedron",
            boundary_sections=[]
        )
    ...
```
The `polycrystalx` mesh input object takes a name (as all do), the source (in this case a gmsh file), the extents (for determination of boundary sections), the cell type (tetrahedron) and boundary sections (for a box, these are determined from the extents).

For this demo example, the meshes are built using a script `run_neper.py` that uses `neper` and `gmsh` to make the mesh input file, which is placed in the `jobs/data` directory.  So the `key` here is the basename of the mesh.


### Microstructure Input
For the `polycrystalx` microstructure input, here is the code for the polycrystal (microstructure) input:

```
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
    ...
```

In this example, there is only one microstructure--from the unloaded state of the far measurements. Because this is a `neper` generated mesh, the grain IDs (cell tags) are read with the mesh, thus explaining the `use_meshtags` flag. The `polycrystal` component in this case is a dummy microstructure (polycrystal.microstructure.Microstructure) used only for the orientations.  The orientations are read from the far field data file, and the grain IDs are generated in `neper` using a Voronoi tessellation with far field centroids as seeds.

### Deformation Input

For the `polycrystalx` deformation input, here is the beginning of the class:

```
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
    ...

```

In this example, our deformation is torsion, and we have to provide the body force density and the boundary conditions, being some combination of specified displacements and tractions.  From the experimental data, there is a list of angular displacements for each state.  We construct the appropriate boundary conditions from that data.
