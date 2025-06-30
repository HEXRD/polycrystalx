# Thermal Expansion Demo

This suite provides three simple examples illustrating how to use thermal expansion in the linear elasticity formulation for an isotropic material.

## Inputs
The mesh is a simple box mesh on the unit cube.  The microstructure is just a single crystal (material).  There are two materials available: an abstract material, named "identity-iso", that has a stiffness of the identity, _i.e._ the stress is the same as the strain, and a realistic material named "ti-64-bar-RT", which has properties for titanium (6 Al, 4V) bar at room temperature.
## Examples
There are three examples given.
1. This uses the "identity-iso" material with a temperature difference of 300C. It has displacement boundary conditions of zero displacement on the whole boundary. The result will be that the mechanical displacement (solution) will be zero. The stress will be the negative of the thermal expansion, because the stress is the stiffness applied to the mechanical strain less the thermal expansion.

2. This is similar to the first except that the boundary conditions are displacements designed to produce a mechanical strain that matches the constant thermal strain exactly. In this case, the stress field will be zero.

3. This example uses a thermal expansion that varies linearly in the x-coordinate.  The applied boundary condtions are again zero, so that the stress field will be the negative but scaled by the (three times) bulk modulus. This uses the "ti-64-bar-RT" material.
