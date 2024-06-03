"""Extended XDMFFile class to include read_function capability"""
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import h5py

from dolfinx import fem, io
from mpi4py import MPI


class XDMFFile_Ext(io.XDMFFile):
    """XDMFFile extension for reading functions"""

    fs_etype = {
        "Node": ("Lagrange", 1),
        "Cell": ("DG", 0),
    }

    def __init__(self, comm, filename, file_mode):
        self._comm = comm
        self._filename = filename
        self._file_mode = file_mode
        super().__init__(comm, filename, file_mode)

    def read_function(self, msh, name):
        """Read function"""
        # Read XML tree to find data element.
        tree = ET.parse(self._filename)
        root = tree.getroot()
        domain = root[0]
        att = self._find_att(domain, name)

        # Get data set info.
        att_type = att.get("AttributeType")
        att_cent = att.get("Center")

        data_item = att[0]
        data_dims = data_item.get("Dimensions")
        nv, nc = shp = (int(i) for i in data_dims.split())

        root_path = Path(self._filename).parent
        data_path = data_item.text
        f, p = data_path.split(":")

        # Now, load function data serially from HDF5.
        with h5py.File(root_path / f, "r") as hfile:
            dg = hfile[p]
            values = dg[:].flatten()

        # Create fenicsx Function.
        V = self._fspace(msh, att_type, att_cent)
        u = fem.Function(V)

        if att_cent == "Node":
            original_ind = np.array(msh.geometry.input_global_indices)
        elif att_cent == "Cell":
            original_ind = np.array(msh.topology.original_cell_index)
        else:
            raise ValueError("Attribute center must be 'Node' or 'Cell'")

        ind = np.tile(nc * original_ind, (nc, 1))
        if nc > 1:
            for i in range(1, nc):
                ind[i, :] += i

        u.x.array[:] = values[ind.T.flatten()]

        return u

    def _find_att(self, elem,  name):
        """find Attribute of elem with given name"""
        # First, find Grid by that name.
        grid = None
        for e in elem.findall("Grid"):
            if e.get("Name") == name:
                grid = e

        if grid is None:
            raise ValueError(f"no matching grid for {name}")

        att = grid[0].find("Attribute")
        return att

    def _fspace(self, msh, datatype, datacenter, shp):
        """Function space from mesh and data type and center"""
        etype = self.fs_etype[datacenter]
        if datatype == "Scalar":
            eshape = None
        elif datatype == "Vector":
            eshape = shp[1]
        elif datatype == "Tensor":
            eshape = (3, 3)

        if eshape is None:
            V = fem.functionspace(msh, etype)
        else:
            V = fem.functionspace(msh, etype, shape=eshape)
        return V
