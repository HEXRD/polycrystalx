"""Deformation input templates"""
import numpy as np
from dolfinx import fem, mesh
import ufl

from ..forms.linear_elasticity import Traction
from .function import FunctionLoader


class DefmLoader:
    """Base class for deformation loaders

    Parameters
    ----------
    job: inputs.Job instance
       input specification
    """

    def __init__(self, defm_input):
        self.defm_input = defm_input


class LinearElasticity(DefmLoader):

    def displacement_bcs(self, V, bdict):
        """Return list of Dirichlet BCs for this problem

        Parameters
        ----------
        V: dolfinx FunctionSpace
           the vector function space
        bdict: dict
           the boundary dictionary

        Returns
        -------
        list
           list of displacement (Dirichlet) boundary conditions
        """
        bdim = V.mesh.topology.dim - 1
        dbcs = []
        for dbc in self. userinput.displacement_bcs:
            facets = bdict[dbc.section]
            if dbc.component is None:
                dofs = fem.locate_dofs_topological(
                    V=V, entity_dim=bdim, entities=facets
                )
                ubc = fem.Function(V)
                ubc.interpolate(dbc.value)
                bc = fem.dirichletbc(value=ubc, dofs=dofs)
            else:
                Vbc = V.sub(dbc.component)
                Vc, _dofmap = Vbc.collapse()
                dofs_vector = fem.locate_dofs_topological(
                    V=V.sub(dbc.component), entity_dim=bdim, entities=facets
                )
                dofs_scalar = fem.locate_dofs_topological(
                    V=Vc, entity_dim=bdim, entities=facets
                )
                dofs = [dofs_vector, dofs_scalar]
                ubc = fem.Function(Vc)
                ubc.interpolate(dbc.value)
                bc = fem.dirichletbc(value=ubc, dofs=dofs, V=Vbc)
            dbcs.append(bc)

        return dbcs


    def traction_bcs(self, V, bdict):
        """Return list of traction BCs for this problem

        Parameters
        ----------
        V: dolfinx FunctionSpace
           the vector function space
        bdict: dict
           the boundary dictionary

        Returns
        -------
        list
           list of traction (natural) boundary conditions
        """
        if len(self.defm_input.traction_bcs) == 0:
            return []
        bdim = V.mesh.topology.dim - 1
        #
        # First, create the surface measure subdomain data using defined by
        # meshtags.
        #
        flist, vlist = [], []
        for i, tbc in enumerate(self.defm_input.traction_bcs):
            facets = bdict[tbc.section]
            flist.append(facets)
            vlist.append((i + 1) * np.ones(len(facets), dtype=np.int32))
        mtags = mesh.meshtags(
            V.mesh, bdim, np.hstack(flist), np.hstack(vlist)
        )
        ds = ufl.Measure("ds", subdomain_data=mtags)
        #
        # Next, create the array of traction forms.
        #
        tbcs = []
        for i, tbc in enumerate(self.defm_input.traction_bcs):
            Vbc = V if tbc.component is None else V.sub(tbc.component)
            ubc = fem.Function(Vbc)
            ubc.interpolate(tbc.value)
            t = Traction(ubc, ds(i + 1), tbc.component)
            tbcs.append(t)

        return tbcs


    def force_density(self, V):
        """Return force density function

        Parameters
        ----------
        V: dolfinx FunctionSpace
           vector function space for force density

        Returns
        -------
        dolfinx Function
           body force function as specified
        """
        if self.defm_input.force_density is not None:
            return FunctionLoader(self.defm_input.force_density).load(V)

    def plastic_distortion(self, T):
        """Return plastic distortion function

        Parameters
        ----------
        T: dolfinx FunctionSpace
           tensor function space for force density

        Returns
        -------
        dolfinx Function
           plastic distortion function as specified
        """
        if self.defm_input.plastic_distortion is not None:
            return FunctionLoader(self.defm_input.plastic_distortion).load(T)


class HeatTransfer(DefmLoader):
    """Loader for heat transfer inputs"""

    def body_heat(self, V):
        pass

    def temperature_bcs(self):
        pass

    def flux_bcs(self):
        pass
