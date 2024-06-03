Developer Documentation
==============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


General Considerations
++++++++++++++++++++++

Try to adhere to `PEP 8 <https://peps.python.org/pep-0008/>`_.
For docstrings, we use the
`numpy style guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

Design
++++++

The package has four main components: ``inputs``, ``loaders``, ``forms`` and
``processes``. Their roles are shown in the diagram below. The user builds
**input specifications** using the ``inputs`` pacakge and the external material
library, ``polycrystal``.
When the simulation is run, the ``processes`` package determines which material
process is being used and generates blank forms using the ``forms`` library.
Then the process calls the appropriate loader, which then operates on the
input specifications to deliver the **form coefficients** for this problem.
Once the form coefficients have been set, the system is solved and the
**solution** is delivered.

.. image:: ../_static/polycrystalx-design.pdf

Inputs
------
Inputs are fairly simple specifications for the mesh and for coefficients
used in the equations being solved.

Loaders
-------
Loaders are python modules for delivering the dolfinx objects from the input
specifications.

Forms
-----
This sets up dolfinx forms for particular problems. On instantiation, the forms
are initialized to have generic coefficients, usually set to zero. After the
user input is read, the coefficients are set using the loaders.

Processes
---------
Processes refer to solving equations for particular models.
