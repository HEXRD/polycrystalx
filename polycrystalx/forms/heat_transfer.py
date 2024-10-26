"""Forms for Heat Transfer"""
from collections import namedtuple
from dolfinx import fem
from ufl import dot, inner, grad, sym, dx, TrialFunction, TestFunction

from .common import tosample


_coeffs = ["orientation", "stiffness", "body_heat", "fluxes"]
Coefficients = namedtuple("Coefficients", _coeffs)
Coefficients.__doc__ = """Coefficients for heat transfer forms

Parameters
----------
orientation: dolfinx Function
    orientation (rotation matrix) field
stiffness: dolfinx Function
    stiffness matrix field
body_heat: dolfinx Function
    body heat density field
fluxes: list of forms.heat_transfer.Flux
    list of applied tractions
"""

_flux_fields = ["value", "ds"]
Flux = namedtuple("Flux", _flux_fields)
Flux.__doc__ = """Flux specification

Parameters
----------
value: (dolfinx) Function
   flux value to apply
ds: Measure
   surface measure of boundary section
"""


class HeatTransferProblem:
    """Heat Transfer Forms

    This class provides the forms for anisotropic thermal conductivity. The
    forms are initialized with all coefficients being set to zero and fluxes
    set to an empty list.  When user input is loaded, the coefficient values
    and applied fluxes are reset to the desired values.

    Parameters
    ----------
    msh: Mesh instance
       the mesh
    opts:
       options (not yet active)
    """

    def __init__(self, msh, opts=None):

        self.msh = msh
        self.opts = opts

        self.V = fem.functionspace(self.msh, ("CG", 1))
        self.T = fem.functionspace(self.msh, ("DG", 0, (3, 3)))
        self.V3 = fem.functionspace(self.kmsh, ("DG", 0, (3,)))

        self._make_coefficients()

    def _make_coefficients(self):
        orient = fem.Function(self.T)
        stiff = fem.Function(self.T)
        bodyh = fem.Function(self.V)
        fluxes = []
        self._coeffs = Coefficients(orient, stiff, bodyh, fluxes)

    @property
    def coefficients(self):
        """Return a tuple of required coefficients"""
        return self._coeffs

    @property
    def forms(self):
        """Generate forms linear and bilinear form"""
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        c = self.coefficients

        a = inner(tosample(c.stiffness, c.orientation) * grad(u), grad(v)) * dx

        # Initialize the linear functional.

        if c.bodyh is None:
            L = dot(fem.Constant(self.msh, 0), v) * dx
        else:
            L = c.bodyh * v * dx

        # Add in the boundary fluxes.

        for flux in c.fluxes:
            L += flux.value * v * flux.ds

        return a, L
