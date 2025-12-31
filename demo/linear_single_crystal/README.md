 ## Linear Single Crystal Demo

 The demo gives a suite of simple test problems on homogeneous deformations of a single crystal. The problems demonstrate how to specify materials, meshes and various boundary conditions. The materials are simple and unrealistic. The meshes are also simple --- rectangular meshes in three dimensions with regular subdivisions. Each deformation is defined by a linear displacement field, `u = Ax`.  The boundary conditions are determined from the deformation and can be specified as either displacements or a combination of displacements and tractions.

 There are four input specifications for each problem: material, mesh, polycrystal, and deformation.  Each input specification has a `key` that is used to full specify the inputs. We will discuss each input below.

### Material
 The demo has a dictionary of elastic materials. Each material is an elastic `SingleCrystal` instance, from the companion `polycrystal` package. The key for the material input is the name of the material, the key in the material dictionary.

 Choices are:     "identity-iso" "iso-21", "cubic-211", "cubic-121", "cubic-112" and "cubic-321". The first two are isotropic: one gives the identity matrix as the stiffness, and the other has a bulk modulus of 2 and a shear modulus of 1.  The last four are cubic materials and have variations on the bulk modulus and two shear moduli.

See the file `linear_single_crystal/job_data/material.py`.

### Polycrystal

The `polycrystal` specification essentially provides the material assignment. The body is divided in to a finite number of `crystals`, and each crystal is assigned a material (phase) and an orientation (rotation matrix). In this case, there is only one crystal, so all that needs to be specified is its orientation. The key for this input is a pair of numbers `(ax_num, ang_deg`), which gives a coordinate axis (0, 1 or 2) and an angle in degrees. The crystal orientation will be the rotation through that angle and about that axis.; *e.g.* (2, 30) is the key for a rotation of 30 degrees about the z-axis.

See the file `linear_single_crystal/job_data/polycrystal.py`.

### Mesh
The mesh is rectangular with extents of 1 in x-direction, 2 in the y, and 3 in the z. The key is the 3-tuple of divisions in each direction; *e.g.* `(10, 20, 30) is the key for the mesh with 10 divisions in x, 20 in y an d30 in z.

See the file `linear_single_crystal/job_data/mesh.py`.

### Deformation
For this suite, the deformation is the most detailed input.  As said above, the user specifies a matrix `A`, which will determine the solution.  Other inputs determine the boundary conditions to be used for the problem. They are determined from the the matrix `A`, the material stiffness and crystal orientation. The key for this input is a 3-tuple `(bcs, i, j)`.  The `bcs` item is a string describing the boundary condition type. It can be one of `"full"`, `"zmax-traction"`, `"zmax-traction-z"`, and `"zmax-traction-xy"`. For "full", the boundary conditions are completely displacements; for "zmax-traction", they are traction BCs on the top z-suraface and displacement BCs everywhere else; for "zmax-traction-z", the z-component of traction and x- and y-displacements are specified on the top and displacements everywhere else; for "zmax-traction-xy", the x- and y- traction components and the z-displacement are specified on the top, with displacements everywhere else. The matrix `A` is a matrix of all zeros, except the `i, j` component, which is 1.

See the file `linear_single_crystal/job_data/deformation.py`.

### Running

**Single Job**
To run a single job, set the job keys in the `linear_single_crystal/__init__.py` file and use the `pxx_job` command.  Here is an example `__init__.py` file. It uses a simple cubic material, aa single crystal rotated 45 degrees about the z-axis, a very coarse mesh (5 x 6 x 7), and yz-shear with purely displacement boundary conditions.
```"""Inputs Module"""
from . import batch

matl_key = "cubic-321"
poly_key = (2, 45)
mesh_key = (5, 6, 7)
defm_key = ("full", 1, 2)

jobkey = (matl_key, poly_key, mesh_key, defm_key)
job = batch.get_job(jobkey)
```
To run it with two processes, use the command:
```mpirun -n 2 pxx_job linear_single_crystal```

**Batch Job**
To run a batch job, use the `pxx_suite` command on a job iterator, which are intended to be set in the `batch.py` file. In the demo batch file, we have three jobs defined:

all_A
: This suite iterates on all nine possible deformation matrices `A` for a fixed cubic material, a fixed crystal orientation and a fixed mesh with 30 subdivisions in each direction, and full displacement boundary conditions.

all_materials
:  This suite iterates on all materials for a fixed mesh, fixed crystal orientation, and a fixed deformation with full displacement boundary conditions.

vary_orientation
: This varies the crystal orientation for a fixed mesh, fixed material and a fixed deformation with traction boundary conditions on the top z-surface and displacements everywhere else.

vary_bcs
: This varies the boundary conditions for a fixed deformation matrix `A` and a single material, a single orientation and a single mesh.

To run the `vary_bcs` suite, use:
```
pxx_suite -n 2 -k vary_bcs lsc.batch
```
