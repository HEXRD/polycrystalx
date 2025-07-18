"""Function Input Specification"""

from collections import namedtuple


_Function = namedtuple(
    "_Function", ["source", "value", "function", "file", "name"], defaults=4 * [None]
)


class Function(_Function):
    """Input specification for a function

    This specifies inputs needed to create a dolfinx Function. The `source` input is
    required and has three possible values. If it is "constant", then the `value`
    parameter is needed.  It can be a single value or a tuple of values. If `source` is
    "interpolation", then the `function` argument is required, and its value is a
    python    function that can be used for interpolation. If `source` is "xdmf", both
    the `file` and `name` arguments are required. They give the name of the XDMF file
    and the name of the variable within it.

    Parameters
    ----------
    source: {"constant", "interpolation", "xdfm"}, required
       source of function data; "constant", "interpolation"  or "xdfm"
    value: float or tuple of floats (optional)
       the value for constant (if source == "constant")
    function: python function
       python function used to interpolate (if source == "interpolation")
    file: str
       name of file with saved function (if source == "xdmf")
    name: str
       name of function in the XDMF file (if source == "xdmf")
    """

    def __init__(self, *args, **kwargs):
        super(Function, self).__init__()
        self.check_inputs()

    def check_inputs(self):

        print("source: ", self.source)
        if self.source == "constant":
            if self.value is None:
                emsg = '"value" keyword not specified for constant function"'
                raise RuntimeError(emsg)
        elif self.source == "interpolation":
            if self.function is None:
                emsg = '"function" keyword not specified for interpolated function'
                raise RuntimeError(emsg)
        elif self.source == "xdfm":
            if (self.file is None) or (self.name is None):
                emsg = '"file" and "name" must both be specified for function from XDMF'
                raise RuntimeError(emsg)
        else:
            msg = (
                'no valid "source" type was found: '
                '"source" must one of "constant" | "interpolation" | "xdmf"'
            )
            raise RuntimeError(msg)
