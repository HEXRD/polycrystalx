# polycrystalx

This is framework for finite element modeling of polycrystalline materials using _fenicsx_. It provides wrappers for _fenicsx_ that simplify user inputs for materials,  meshes, functions, and boundary conditions. It also a provides a batch capability. The companion repository, _polycrystal_, implements the material models used in the  _polycrystalx_ package.
# Installation
##  fenicsx
First you need to install fenicsx.  The _main_ branch of _polycrystalx_ uses _fenicsx 0.9_, so these instructions are for that branch. A version of _polycrystalx_ for _fenicsx 0.10_, the newest version, is being verified and should be available soon.

For full details on installing `fenicsx`, see fenics project [downloads page](https://fenicsproject.org/download/). Here are the instructions using _conda_ for _fenicsx 0.9_.
```
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx=0.9 mpich pyvista
```
## polycrystalx and polycrystal
Next, install _polycrystalx_ and _polycrystal_ using _pip_. After cloning or downloading the packages install using _pip_.  The _polycrystal_ package is not directly used by the _polycrystalx_, but it will be needed to set up the actual problems. The commands below install both packages in place (so that any changes you make will take effect in this environment) and with any required packages. So go to the directory above the downloaded packages and run:

```
pip install -r polycrystal/requirements.txt
pip install -e polycrystal
pip install -r polycrystal/requirements.txt
pip install -e polycrystalx
```

# Citations

We do not yet have a reference for this work, but do not forget to cite the fenicsx project.  See [Citing FEniCS](https://fenicsproject.org/citing).

# See Also
[FEniCS Project](https://fenicsproject.org/)


# Acknowledgments

This work has been developed over years under contracts with the Air Force Research Laboratory, which has generously approved the work to be released as open source.
