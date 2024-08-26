"""Polycrystal loaders"""
import numpy as np
import dolfinx


class Polycrystal:

    def __init__(self, userinput, parent=None):
        """load polycrystal input

        Parameters
        ----------
        userinput: inputs.polycrystal.Polycrystal instance
           user input for polycrystal configuration
        """
        self.userinput = userinput
        self.use_meshtags = self.userinput.use_meshtags
        self._num_grains = self.polycrystal.num_grains

    # We could put a setter here and do checks on the input in the setter.
    @property
    def polycrystal(self):
        return self.userinput.polycrystal

    @property
    def orientation_list(self):
        return self.polycrystal.orientation_list

    @property
    def num_grains(self):
        """number of grains"""
        return self._num_grains

    def grain_cell_tags(self, msh):
        """Assign grain IDs to cells
        Parameters
        ----------
        msh: Mesh instance
           the mesh

        Returns
        -------
        cell_tags: meshtags instance
            grain IDs
        """
        indmap = msh.topology.index_map(msh.topology.dim)
        num_cells = indmap.size_local
        midpoints = dolfinx.mesh.compute_midpoints(
            msh, msh.topology.dim, np.arange(num_cells, dtype=np.int32)
        )
        gids = self.polycrystal.grain(midpoints).astype(np.int32)
        cell_tags = dolfinx.mesh.meshtags(
            msh, msh.topology.dim, np.arange(num_cells, dtype=np.int32), gids)

        return cell_tags

    def grain_cell_dict(self, cell_tags):
        # Partition cells by grain id
        # * Note: do not use [] as the value in dict.fromkeys because all will
        #   be the same list. Set each value to [] individually.
        gids = cell_tags.values
        num_cells = len(gids)

        gcell_d = dict.fromkeys(range(self.num_grains))
        for k in gcell_d:
            gcell_d[k] = []
        for i in range(num_cells):
            gcell_d[gids[i]].append(i)

        # For fenicsx 0.8.
        for k in gcell_d:
            gcell_d[k] = np.array(gcell_d[k], dtype=np.int32)

        return gcell_d

    def orientation_field(self, T, grain_cells):
        """Orientation Field

        PARAMETERS
        ----------
        T: tensor function space
           then function space for orientations
        grain_cells: dictionary
           map of cell ID arrays to grains

        RETURNS
        -------
        tensor valued function
           the orientation field on the mesh
        """
        ori_fld = dolfinx.fem.Function(T)
        for gi in range(self.polycrystal.num_grains):
            cells = grain_cells[gi]
            ori = self.orientation_list[gi]
            orif = lambda x: self._orif(ori, x)
            ori_fld.interpolate(orif, cells)

        return ori_fld

    def _orif(self, ori, x):
        return np.tile(ori.reshape(9,1), x.shape[1])
