# MHDG code, version 3.0
Transport code for plasma simulations based on a high-order HDG (hybrid discontinuous Galerkin) scheme.

# Models available
0- Laplace model: simple solver for anisotropic diffusion equation, 1 equation.

1- Isothemal model: model with 2 equations, unknowns density and parallel velocity.

2- Non-isothermal model: model with 4 equations, unknowns density, parallel velocity, ion temperature, electron temperature.

3- Isothermal + vorticity: model with 4 equations, unknowns density, parallel velocity, vorticity, electric potential (work in progress)

# Dimensions

1- 2D cartesian and axisymmetric simulations.

2- 3D cartesian (slab geometry) and axisymmetric (toroidal geometry) simulations.

# Parallelization

Parallelization is based on a hybrid OpenMP/MPI paradigm. Domain decomposition in the plane + toroidal divisions


