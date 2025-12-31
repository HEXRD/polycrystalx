"""Inputs Module"""
from . import batch


matl_key = "cubic-321"
poly_key = (2, 45)
mesh_key = (5, 6, 7)
defm_key = ("full", 1, 2)

jobkey = (matl_key, poly_key, mesh_key, defm_key)
job = batch.get_job(jobkey)
