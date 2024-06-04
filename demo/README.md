# Linear Single Crystal

This example uses simple test problems to illustrate how to set up a single job or a suite of jobs. The problems are single crystal, meaning there is just one
material and one orientation, so it's easy to set up problems where you know
the answer.

## Inputs
The first inputs are `suite` and `process`.  The suite is the name for the
whole group of simulations, `linear_single_crystal`. The process is the name
of the process model; here it is `linear-elasticity`.

### Materials
In the `data` module, there is a dictionary of materials. It defines two
isotropic materials and three cubic materials.  All are abstract and do not
reflect actual materials.

### Polycrystal
A simple polycrystal is set up. It has one crystal and so one orientation. For
the isotropic materials, the orientation doesn't matter at all (as long as
it's a rotation matrix), but it does matter for the cubic materials.  At this point, only the identity is used, but the mechanism to add more single orientaitons is there. The tuple `(axis, angle)` is the key to generate a
rotation through a given angle about one of the coordinate axes. As mentioned
above, only `(0, 0)`, representing the identity, is used currently.

### Meshes
These are simple box meshes, and the tuple of divsions is used to generate
the mesh.  These are purposely small since we are running problems with an
exact answer.

### Deformations
For the test problems, we set up the problem data to have exact solutions of
the for `u = Ax`.  We use nine choices for `A`; each one has zeros in all
positions except for 1 spot, which has value 1.  Body forces are zero. We apply
several types of boundary conditions:
* full displacement (Dirichlet) boundary conditions on the whole boundary
* traction boundary conditions on the top surface and displacement everywhere else
* traction boundary conditions on the top in the z-component only, and displacements everywhere else
* traction boundary conditions on the top in the `xy` directions only, and displacements everywhere else.

## Single Job

To set up a single job, you just create `Job` instance named `job`. It takes
the names of the suite and the process.  Then there are the four inputs
described above. To run the single job:

```mpirun -n 2 pxx_job linear_single_crystal```

## Batch Job
To run in batch, use the `pxx_suite` script. In the `batch` module, an iterator
is set up to generate various combinations of the various inputs. The iterator
generates keys that are used to create the job, and a `get_job(key)` function
is used by the processing script to generate the actual job. Run like this:

```pxx_suite -n 2 linear_single_crystal.batch```
