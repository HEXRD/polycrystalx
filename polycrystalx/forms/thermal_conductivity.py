"""Forms for Thermal Conductivity"""
from dolfinx import fem
from ufl import dot, inner, grad, sym, dx, TrialFunction, TestFunction


_coeffs = ["orientation", "stiffness", "body_heat", "fluxes"]
Coefficients = namedtuple("Coefficients", _coeffs)
Coefficients.__doc__ = """Coefficients for thermal conductivity forms

Parameters
----------
orientation: dolfinx Function
    orientation (rotation matrix) field
stiffness: dolfinx Function
    stiffness matrix field
body_heat: dolfinx Function
    body heat density field
fluxes: list of forms.thermal_conductivity.Flux
    list of applied tractions
"""

_flux_fields = ["value", "ds"]
Flux = namedtuple("Flux", _trac_fields)
Flux.__doc__ = """Flux specification

Parameters
----------
value: (dolfinx) Function
   flux value to apply
ds: Measure
   surface measure of boundary section
"""


class ThermalConductivity:
    """Thermal Conductivity

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
        self.V3 = fem.functionspace(self.msh, ("DG", 0, (3,)))

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

        a = inner(sigs_3x3(u, c.stiffness, c.orientation), sym(grad(v))) * dx

        # Initialize the linear functional.
        if c.body_force is None:
            L = dot(fem.Constant(self.msh, (0, 0, 0)), v) * dx
        else:
            L = dot(c.body_force, v) * dx

        if c.plastic_distortion is not None:
            Cbeta_form = self._C_beta_s_form(
                c.plastic_distortion, c.stiffness, c.orientation
            )
            L += inner(Cbeta_form, sym(grad(v))) * dx

        # Add in the tractions.
        for trac in c.tractions:
            if trac.component is None:
                L += inner(trac.value, v) * trac.ds
            else:
                raise NotImplementedError(
                    "traction components not yet implemented--use full vector "
                    "with zero components as needed"
                )

        return a, L

    @staticmethod
    def _C_beta_s_form(beta_s_3x3, stiff_c, orient):
        """apply stiffness to a tensor field

        beta_s - the tensor field in sample frame (in sample frame)
        stiff_c - the stiffness matrix (in crystal frame)

        RETURNS

        Cbeta_s - form for the stiffness times the tensor field (in sample
                  frame)
        """
        # Note that the stiffness matrix operates only on symmetric matrices.
        beta_c_3x3 = tocrystal(sym(beta_s_3x3), orient)
        beta_c_6 = to6vector(beta_c_3x3)
        Cbeta_c_6 = dot(stiff_c, beta_c_6)
        Cbeta_c_3x3 = totensor(Cbeta_c_6)
        Cbeta_s_3x3 = tosample(Cbeta_c_3x3, orient)

        return Cbeta_s_3x3
