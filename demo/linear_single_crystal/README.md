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
```e
## Running a Single Job
To run a single job, use the `pxx_job` script.

```mpirun -n 2 pxx_job linear_single_crystal```

## Batch Job
To run in batch, use the `pxx_suite` script. In the `batch` module, an iterator
is set up to generate various combinations of the various inputs. The iterator
generates keys that are used to create the job, and a `get_job(key)` function
is used by the processing script to generate the actual job. Run like this:

```pxx_suite -n 2 linear_single_crystal.batch```
