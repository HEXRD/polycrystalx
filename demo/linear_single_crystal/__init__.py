"""Inputs Module for lienar single crystal"""

from . import batch

# Material has identity for stiffness matrix.
matl_key = "identity-iso"

# Single crystal with no rotation.
poly_key = (0, 0)

# Mesh is a box with 40 divisions in each directon.
mesh_key = (40, 40, 40)

# The applied deformation is a traction on the top and a displacement field
# u = Ax where
#     [1 0 0]
# A = [0 0 0]
#     [0 0 0]
defm_key = ("zmax-traction", 1, 1)

job_key = (matl_key, poly_key, mesh_key, defm_key)

job = batch.get_job(job_key)
