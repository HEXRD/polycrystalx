"""Set up microstructure from far field data

Notes
-----
* We are using a subset of the centroids at load 0 to seed the Voronoi. First, we
  are using a rectangular geometry that is inside of the cylindrical geometry
  of the actual experiment. Then we are only using a subset of these, just to keep
  the problem size small (around 100 grains) for the demo.
* We are also using the orientations from load 0.
"""
from set_path import DOTDOT
import numpy as np

from jobs.data import DATA_DIR, EXTENTS, load_ffdata


def in_box(x, a):
    """Return a boolean array indicating which rows of x lie inside the box a.

    Parameters
    ----------
    x : ndarray, shape (n, 3)
    a : ndarray, shape (3, 2)
        a[i] = (min, max) range for coordinate i

    Returns
    -------
    ndarray of bool, shape (n,)
    """
    return np.all((x >= a[:, 0]) & (x <= a[:, 1]), axis=1)


# We use the centroids at first load state, and we are just using a subset for this
# demo.

ffd0 = load_ffdata(0)
seeds = ffd0.centroid
s_in = in_box(seeds, EXTENTS)
seeds = seeds[s_in]
select = np.arange(0, len(seeds), 7) # gives 108 grains
np.savetxt(DATA_DIR / "seeds.txt", seeds[select])
print(f"wrote {len(select)} centroids to 'seeds.txt'")

# Next we write the orientations corresponding to the grain centroids above.
# Convert orientations to rotation matrices.
oris = ffd0.orientations[select]
np.save(DATA_DIR / "orientations", oris)
print(f"wrote orientations to 'orientations.npy'")
