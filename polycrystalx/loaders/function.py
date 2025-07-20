"""Load functions from input specifications"""
import numpy as np

from dolfinx import fem

from ..utils import XDMFFile_Ext


class FunctionLoader:
    """Loader for dolfinx Functions

    Parameters
    ----------
    userinput: inputs.function.Function
       input fucntion specification
    """

    def __init__(self, userinput):
        self.userinput = userinput
        self.source = userinput.source

    def load(self, V, name="f"):
        """Load function from input spec

        Parameters
        ----------
        V: FunctionSpace (dolfinx)
           the function space associated with this function
        name: str (defaults to "f"), required for xmdf input
           name to use in reading XDMF file
        """
        f = fem.Function(V, name=name)
        if self.source == "constant":
            self._load_constant(f)
        elif self.source == "interpolation":
            self._load_interpolate(f)
        elif self.source == "xdmf":
            f = self._load_xdmf(f, V.mesh)
            f.name = name
        else:
            msg = (
                'no valid "source" type was found: '
                '"source" must "constant" | "interpolation" | "xdmf"'
            )
            raise RuntimeError(msg)
        return f

    def _load_constant(self, f):
        if self.userinput.value is None:
            raise RuntimeError("value not specified for constant function")
        value = self.userinput.value
        f_interp = lambda x: np.tile(value, (x.shape[1], 1)).T
        f.interpolate(f_interp)

    def _load_interpolate(self, f):
        if self.userinput.function is None:
            msg = '"function" key not specified for  function'
            raise RuntimeError()
        f_interp = self.userinput.function
        f.interpolate(f_interp)

    def _load_xdmf(self, f, msh):
        filename = self.userinput.file
        funcname = self.userinput.name
        if (filename is None) or (funcname is None):
            msg = '"filename" and "funcname" must both be specified'
            raise RuntimeError(msg)

        with XDMFFile_Ext(msh.comm, filename, "r") as xfile:
            f = xfile.read_function(msh, funcname)

        return f
