"""Forms commonly used in polycrystal models"""
from collections import namedtuple

import numpy as np

from ufl import dot, sym, grad, as_vector, as_matrix, Measure
from dolfinx import fem


_s2 = np.sqrt(2)
_s2i = 1./_s2


# These are routines for representation of 3x3 symmetric matrices as 6-vectors.
def to6vector(w3x3):
    """Return 6-vector form of 3x3 matrix.

    Parameters
    ----------
    w3x3: Expression
       expression for 3X3 symmetric matrix

    Returns
    -------
    Expression:
       expression for 6-vector of unique components of `w3x3`
    """
    #
    return as_vector([w3x3[0, 0], w3x3[1, 1], w3x3[2, 2],
                      _s2*w3x3[1, 2], _s2*w3x3[2, 0], _s2*w3x3[0, 1]])


def totensor(w6):
    """Reconstruct symmetric matrix from 6-vector.

    Parameters
    ----------
    w6: Expression
       expression for a 6-vector

    Returns
    -------
    Expression:
       expression for symmetric matrix derived from `w6`
    """
    #
    return as_matrix([[     w6[0], _s2i*w6[5], _s2i*w6[4]],
                      [_s2i*w6[5],      w6[1], _s2i*w6[3]],
                      [_s2i*w6[4], _s2i*w6[3],      w6[2]]])


# These convert between cyrstal and sample reference frames.
def tocrystal(w3x3, orient):
    """Convert matrix from sample to cyrstal.

    Parameters
    ----------
    w3x3: Expression
       expression for matrix in sample coordinates
    orient: Expression
       3X3 rotation matrix for change of basis

    Returns
    -------
    Expression:
       expression for `w3x3` in crystal coordinates
    """
    return orient.T*w3x3*orient


def tosample(w3x3, orient):
    """Convert matrix from cyrstal to sample.

    Parameters
    ----------
    w3x3: Expression
       expression for matrix in crystal coordinates
    orient: Expression
       3X3 rotation matrix for change of basis

    Returns
    -------
    Expression:
       expression for `w3x3` in sample coordinates
    """
    return orient*w3x3*orient.T


# Following are stress forms for sample and crystal reference frames.
def sigc_6(w_s, stiff_c, orient):
    """Sample displacement field to crystal stress field as 6-vector

    Parameters
    ----------
    w_s: Expression
       displacement field in sample frame
    stiff_c: Expression
       stiffness matrix in crystal frame
    orient: Expression
       3X3 rotation matrix for change of basis

    Returns
    -------
    Expression:
       expression for stress in 6-vector form in crystal frame
    """
    return dot(stiff_c, to6vector(tocrystal(sym(grad(w_s)), orient)))


def sigc_3x3(w_s, stiff_c, orient):
    """Sample displacement field to crystal stress field as 3x3 matrix

    Parameters
    ----------
    w_s: Expression
       displacement field in sample frame
    stiff_c: Expression
       stiffness matrix in crystal frame
    orient: Expression
       3X3 rotation matrix for change of basis

    Returns
    -------
    Expression:
       expression for stress as 3x3 matrix in crystal frame
    """
    return totensor(sigc_6(w_s, stiff_c, orient))


def sigs_3x3(w_s, stiff_c, orient):
    """Sample displacement field to sample stress field as 3x3 matrix

    Parameters
    ----------
    w_s: Expression
       displacement field in sample frame
    stiff_c: Expression
       stiffness matrix in crystal frame
    orient: Expression
       3X3 rotation matrix for change of basis

    Returns
    -------
     Expression:
       expression for stress as 3x3 matrix in sample frame
   """
    return tosample(sigc_3x3(w_s, stiff_c, orient), orient)


# Here are two forms for integration over grains.

def grain_volume(msh):
    """Form for computing grain volumes

    Parameters
    ----------
    msh: (dolfinx) Mesh
        the mesh

    Returns
    -------
    Form
       generic form for integrating i(x) * dx, where `i` is an indicator
       function for a grain
    """
    V = fem.functionspace(msh, ("DG", 0))
    indicator = fem.Function(V)
    dx = Measure("dx", domain=msh)

    return fem.form(indicator * dx), indicator


def grain_integral(msh, V):
    """Form for computing grain integrals of f

    Parameters
    ----------
    msh: (dolfinx) Mesh
       the mesh
    V: (dolfinx) FunctionSpace
       the function space

    Returns
    -------
    form
       form for integrating i(x) * f(x) * dx, where `i` is an indicator
       function for a grain, and `f` is a function to be integrated over the
       grains
    """
    Vind = fem.functionspace(msh, ("DG", 0))
    indicator = fem.Function(Vind)
    func = fem.Function(V)
    dx = Measure("dx", domain=msh)

    return fem.form(indicator * func * dx), indicator, func
