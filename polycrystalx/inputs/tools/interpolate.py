"""Constant interpolation functions"""
import numpy as np


def constant(value):
    """Return python function for interpolation of a constant

    value: float | tuple
       value of the function

    RETURNS
    -------
    function
       function that returns value at each point
    """
    return lambda x: np.tile(value, (x.shape[1], 1)).T


constant_function = constant


def linear(A, b=None):
    """Return python function to interpolate linear or affine map

    Parameters
    ----------
    A: array (n, n)
       matrix
    b: array(n), optional
       vector

    Returns
    -------
    function
       function of `x` that returns `Ax`, or `Ax + b` if `b` is not None

    """
    if b is None:
        return lambda x: np.array(A) @ x
    else:
        return lambda x: np.array(A) @ x + np.array(b)


linear_function = linear
