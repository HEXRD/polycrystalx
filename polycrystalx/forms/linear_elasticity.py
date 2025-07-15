"""Forms for Linear Elasticity"""
from .common import namedtuple, sigs_3x3
from .common import to6vector, totensor, tocrystal, tosample

from dolfinx import fem
from ufl import dot, inner, grad, sym, dx, TrialFunction, TestFunction


_coeffs = ["orientation", "stiffness", "body_force", "plastic_distortion",
           "thermal_expansion", "tractions"]
Coefficients = namedtuple("Coefficients", _coeffs)
Coefficients.__doc__ = """Coefficients for linear elasticty forms

Parameters
----------
orientation: dolfinx Function
    orientation (rotation matrix) field
stiffness: dolfinx Function
    stiffness tensor field
body_force: dolfinx Function
    body force density field
plastic_distortion: dolfinx Function
    plastic distortion tensor field
thermal_expansion: dolfinx Function
    function giving thermal expansion
tractions: list of forms.linear_elasticity.Traction instances
    list of applied tractions
"""

_trac_fields = ["value", "ds", "component"]
Traction = namedtuple("Traction", _trac_fields)
Traction.__doc__ = """Traction specification

Parameters
----------
value: (dolfinx) Function
   traction value to apply
ds: Measure
   surface measure of boundary section
component: int | None
   which component to apply traction (note, this is not yet implemented; just
   use a full vector function with zeros in the other components)
"""


class LinearElasticity:
    """Linear elasticity

    This class provides the forms for anisotropic linear elasticity. The
    forms are initialized with all coefficients being set to zero and tractions
    set to an empty list.  When user input is loaded, the coefficient values
    and applied tractions are reset to the desired values.

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

        self.V = fem.functionspace(self.msh, ("CG", 1, (3,)))
        self.T = fem.functionspace(self.msh, ("DG", 0, (3, 3)))
        self.T6 = fem.functionspace(self.msh, ("DG", 0, (6,6)))

        self._make_coefficients()

    def _make_coefficients(self):
        orient = fem.Function(self.T)
        stiff = fem.Function(self.T6)
        bodyf = fem.Function(self.V)
        pdist = fem.Function(self.T)
        texpand = fem.Function(self.T)
        tracs = []
        self._coeffs = Coefficients(orient, stiff, bodyf, pdist, texpand, tracs)

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

        if c.thermal_expansion is not None:
            # This is just like the plastic distortion, but the thermal
            # expansion tensor is in the crystal frame.

            te_form = self._expansion_form(
                c.thermal_expansion,  c.stiffness, c.orientation
            )
            L += inner(te_form, sym(grad(v))) * dx

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

        beta_s - the tensor field in sample frame
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

    @staticmethod
    def _expansion_form(expand_c_3x3, stiff_c, orient):
        """apply stiffness to a tensor field

        expand_c_3x3 - the expanion tensor field in crystal frame
        stiff_c - the stiffness matrix (in crystal frame)
        orient - the orientation field

        RETURNS

        expand_s_3x3 - form for the stiffness times the tensor field (in sample
                       frame)
        """
        # Note that the stiffness matrix operates only on symmetric matrices.
        expand_c_6 = to6vector(expand_c_3x3)
        stiff_expand_c_6 = dot(stiff_c, expand_c_6)
        stiff_expand_c_3x3 = totensor(stiff_expand_c_6)
        stiff_expand_s_3x3 = tosample(stiff_expand_c_3x3, orient)

        return stiff_expand_s_3x3
