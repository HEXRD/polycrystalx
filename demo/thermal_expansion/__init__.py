"""Inputs module for thermal expansion"""
from . import batch


# matl_key = "ti-64-bar-RT"
matl_key = "identity-iso"
mesh_key = (15, 20, 25)
poly_key = None
defm_key = "linear"

jobkey = (matl_key, mesh_key, poly_key, defm_key)
job = batch.get_job(jobkey)
