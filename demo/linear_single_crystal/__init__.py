"""Inputs Module for lienar single crystal"""
import numpy as np

from polycrystalx import inputs

from .batch import suite, process, get_job


matl = "identity-iso"
poly = (0, 0)
mesh = (20, 20, 20)
defm = ("full", 1, 1)

job = get_job((matl, poly, mesh, defm))
