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
