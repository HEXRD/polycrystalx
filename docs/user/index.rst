User Documentation
==================

.. toctree::
   :maxdepth: 2


Inputs
++++++

The inputs are python objects providing specifications for generating
the actual objects to be used in the FEniCSx simulation. There are
four required input types and one optional. Here, we shall walk through a
simple linear elasticity simulation.

Mesh Input
----------
The mesh input specified how to load the mesh. It can be loaded from an XDMF
file, from a gmsh file, or it can be a simple box to be generated inside
FEniCSx. Here is a simple example. The mesh goes from 0 to 1 in the x-dreiction,
from 0 to 2 in y, and from 0 to 3 in z.  It has 10 subdivisions in each
direction.

::

    from polycrystalx import inputs


    extents = np.array([
        [0, 1],
        [0, 2],
        [0, 3],
    ])

    mesh_input = inputs.mesh.Mesh(
       name="mesh-10x10x10",
       source="box",
       extents=extents,
       divisions=(10, 10, 10),
    )

Material Input
--------------
The material input is provided by the companion package, `polycrystal`. It is
a list of elastic materials. In this case, it is a very simple and abstract
material, which makes the stress and strain equal.

::

   from microstructure.elasticity.single_crystal import SingleCrystal
   from polycrystalx import inputs

   # ===== Stiffness is the identity, stress = strain.

   # crystal = SingleCrystal("cubic", (1.0, 0.5, 0.25))
   crystal = SingleCrystal.from_K_G(1/3, 1/2)

   matl_data = inputs.material.Elasticity(
        name = "identity",
        materials = [crystal]
   )


Polycrystal Input
-----------------
Like the material input, the polycrystal input comes from the companion material
package.

::

    from polycrystalx import inputs
    from microstructure.microstructure.single_crystal import (
        SingleCrystal as SingleCrystalMicro
    )

    orientation = np.idenitity(3)
    poly_data = inputs.polycrystal.Polycrystal(
        name = "single-crystal,
        polycrystal = SingleCrystalMicro(orientation)
    )



Deformation Input
-----------------

::

   scale = (2, 3, 4)
   grad_u = np.diag(scale)

   displacement_bcs = [
      inputs.deformation.DisplacementBC(
         section = "boundary",
	 value = interpolate.linear(grad_u),
      ),
   ]

   defm_data = inputs.deformation.Elasticity(
        name = "rescale-234",
        force_density = zero_vector_function,
        displacement_bcs = displacement_bcs,
        traction_bcs = [],
    )


Running Simulations
+++++++++++++++++++
To run a simulation, the key object is the Job. The Job class holds the
input specifications. In addition to the mesh, material, polycrystal and
deformation inputs, you also provide a name for the suite of simulations
this is part of and for the material process being modeled.


Running a Single Simulation
---------------------------
To run a single simulation, use the command line. For example:

::

    from polycrystalx import inputs
    from polycrystalx import processes

    suite = "linear_single_crystal"
    process = "linear-elasticity"

    job = inputs.job.Job(
        suite = suite,
        process = process,
        mesh_input = mesh_data,
        material_input = material_input,
        polycrystal_input = poly_data,
        deformation_input = deformation_data
    )
    processes.run(job)


There is also a command line script that is run using MPI and takes the
argument of `input_module`, which
is the name or path to a python module with a `job` attribute, as above.


Running a Suite of Simulations
------------------------------
To run multiple simulations, you need to set up a suite of jobs. The way to
do that is to set up a list or iterator of keys (input names for example.) that
can  be used to create jobs.  An example will be forthcoming. This works
with the comand line script `run_suite`.
