# Linear Single Crystal

This example uses simple test problems to illustrate how to set up a single job or a suite of jobs. The problems here are single crystal, meaning there is just one material and one orientation, so it's easy to set up problems where you know the answer.

## Jobs
To run a simulation, you will need to define a `inputs.Job` instance. A `Job` has six components.
```
job = inputs.job.Job(
    suite = suite,
    process = process,
    mesh_input = mesh_input,
    material_input = material_input,
    polycrystal_input = polycrystal_input,
    deformation_input = deformation_input
)
```
The `suite` and the `process` are simple strings. They give the name of the suite of simulations and the name of the material process being modeled. In this example, we have:
```
suite = "linear_single_crystal"
process = "linear-elasticity"
```
The other four components describe the mesh, the material, the polycrystal (microstructure) and the applied deformation. Each of these components will have `key` (hashable type)  that can be used to generate the input. The key may be a string, a number or a tuple, or some combination of these. Each input is then built from the key.

### Materials
In the `data` module, there is a dictionary of elastic materials. It defines two
isotropic materials and three cubic materials.  In this example, the materials are very simple and do not reflect actual materials. The material `key` is simply the name of the material. In this example, we set moduli that gives the identity stiffness matrix.
```
matl_key = "identity-iso"
```

### Polycrystal
In this example, we use a very simple polycrystal. In face, it is a single crystal. The only thing we need to assign is the crystal orientation. We chose to use a rotation through somae angle about one of the coordinate axes. The polycrystal key is a pair of numbers, the first being the axis of rotation: 0, 1, or 2.  The second is the angle of rotation. So for example, if the key were `(1, 45)`, that would represent a single crystal rotated 45 degrees about the y-axis. In this example, the material is isotropic so the orientation doesn't matter. This gives an orientation matrix of the identity.
```
poly_key = (0, 0)
```
### Meshes
These are simple box meshes, and the `mesh key` is the tuple of divsions is used to generate the mesh.  These are purposely small since we are running problems with an
exact answer. Here the mesh is a box with 40 subdivisions in each direction.
```
mesh_key = (40, 40, 40)
```
### Deformations
For the test problems, we set up the problem data to have exact solutions of
the form `u = Ax`.  We use nine choices for `A`; each one has zeros in all
positions except for 1 spot, which has value 1.0.  Body forces are zero. We apply
several types of boundary conditions:
* full displacement (Dirichlet) boundary conditions on the whole boundary
* traction boundary conditions on the top surface and displacement everywhere else
* traction boundary conditions on the top in the z-component only, and displacements everywhere else
* traction boundary conditions on the top in the `xy` directions only, and displacements everywhere else.

In this example, we apply traction on the top surface, and the strain will be zero except for the xx-component.
```
defm_key = ("zmax-traction", 1, 1)
```
## Running a Single Job
To run a single job, run the `pxx_job` script using `mpirun`.

```mpirun -n 2 pxx_job linear_single_crystal```

## Batch Jobs
To run in batch, use the `pxx_suite` script. In the `batch` module, iterators are set up to generate different combinations of the various inputs. Each iterator generates a sequence of job keys that are used to create each job. The `batch.get_job(key)` function is used to create the `Job` instance from the keys.

This example has two iterators. The first one is shown below. It's name is `all_materials_xx`. It creates jobs for every material in the material database, each job using the same meshe and microstructure and with displacement boundary conditions corresponding to an xx-strain field of 1.0. Note that the `itertools.product` is a python function that creates an iterator that takes all combinations of the four inputs. In this case, the material keys are the only keys that have more than item.
```
matl_keys = list(matl_dict.keys())
poly_keys = [(0, 0)]
mesh_keys = [(10, 20, 30)]
defm_keys = [("full", 1, 1)]

# This suite runs all materials on the same mesh with displacements
# corresponding to an xx-strain field.
all_materials_xx = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)
```
Run the suite like this. Note that `mpirun` is not used directly here, but it is called internally by the `pxx_suite` command.
```
pxx_suite -n 2 -k all_materials_xx linear_single_crystal.batch
```

The second batch suite is named `cubic_traction_z`.  It runs various types of traction boundary conditions for the same material, microstructure and mesh. The `zmax-traction` boundary condition applies traction on the top surface (`zmax`) and dipslacements eveverywhere else. The `zmax-traction-z` condition applies the z-component of traction on and x- and y-displacements on the top surface and displacments everywhere else. Finally, the `zmax-traction-xy` applies x- and y-tractions and a z-displacement on the top surface, again with dipslacements eveverywhere else.
```
matl_keys = ["cubic-211"]
poly_keys = [(0, 0)]
mesh_keys = [(10, 20, 30)]
defm_keys = [
    ("zmax-traction", 3, 3),
    ("zmax-traction-z", 3, 3),
    ("zmax-traction-xy", 3, 3)
]

# This suite runs various traction boundary conditions on the top surface
# with a fixed mesh and a cubic material.
cubic_traction_z = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)
```
Run the suite like this.
```
pxx_suite -n 2 -k cubic_traction_z linear_single_crystal.batch
```
