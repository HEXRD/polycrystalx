"""This loads items from the data directory and sets up access tools as needed

The `GrainData` class is for accessing the far field data.

The `twist_by_state` function gives the twist angle (radians( per mm along the twist
axis.
"""
from pathlib import Path
from collections import namedtuple

import numpy as np
from scipy.spatial.transform import Rotation

DATA_DIR = Path(__file__).parents[1] / "data"

FF_FILE = DATA_DIR / "far-field-data.npy"
EXTENTS = np.loadtxt(DATA_DIR / "extents.txt")

# This is the GrainData  base class from hexrd, but we include it here so that
# we do not require the full hexrd package.

_flds = [
    "id", "completeness", "chisq", "expmap", "centroid", "inv_Vs", "ln_Vs"
]
GrainDataBase = namedtuple("GrainDataBase", _flds)


class GrainData(GrainDataBase):

    @property
    def orientations(self):
        """List of rotation matrices"""
        return Rotation.from_rotvec(self.expmap).as_matrix()


def load_ffdata(load):
    """Get far field data set for a load state

    Returns
    -------/
    GrainData
       far field data for 13 loads (0-12)
    """
    ffdata = np.load(FF_FILE)

    nstates, ngrains = ffdata.shape[:2]
    #
    # Initialize the grain data arrays.
    #
    gid = ffdata[load, :, 0]
    cmpl = ffdata[load, :, 1]
    chisq = ffdata[load, :, 2]
    expmap = ffdata[load, :, 3:6]
    cent = ffdata[load, :, 6:9]
    inv_V  = ffdata[load, :, 9:15]
    ln_Vs = ffdata[load, :, 15:21]

    return GrainData(gid, cmpl, chisq, expmap, cent, inv_V, ln_Vs)


def twist_by_state(state):
    """Get twist for given load state

    Parameters
    ----------
    state: int
           state number, between 0 and 12

    Returns
    -------
    float:
        twist per mm along axis (in radians / mm) for given state

    This sets the twist angle (deg) per mm by state number, taken from extensometry
    spreadsheet, "Radial Strains" tab column F (dTheat31/um).
    """
    TO_MM = 1000.0
    deg_twist_per_um = (
        0, 4.5676E-05, 1.0759E-04, 1.3128E-04, 1.7218E-04, 2.1938E-04, -1,
        3.2427E-04, 3.7651E-04, 4.1980E-04, 4.7603E-04, 5.2891E-04, 5.9447E-04
    )
    deg_twist_per_mm = TO_MM * np.array(deg_twist_per_um)

    return np.radians(deg_twist_per_mm[state])
