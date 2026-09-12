"""Batch jobs"""
import itertools

from .job_inputs import get_job


# ==================== Single job.

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
job = get_job(jobkey)


# ==================== Suites of jobs


# Define the suite of jobs.

matl_keys = ["lshr_660C"]
poly_keys = ["ms"]
mesh_keys = ["vor-050"]
defm_keys = list(range(4))

# The `itertools.product` function generates all combinations of items, one from
# each list.  In this case, it will generate all the deformations for a single
# material, microstructure and mesh.
job_keys = itertools.product(matl_keys, poly_keys, mesh_keys, defm_keys)
