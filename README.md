# MHDG code, version 1
Transport code for plasma simulations based on a high-order HDG (hybrid discontinuous Galerkin) scheme.

# Models available
1- Isothemal model: model with 2 equations, unknowns density and parallel velocity.
2- Non-isothermal model: model with 4 equations, unknowns density, parallel velocity, ion temperature, electron temperature.

# Dimensions
1- 2D cartesian and axisymmetric simulations.
2- 3D cartesian (slab geometry) and axisymmetric (toroidal geometry) simulations.

# Parallelization
Parallelization is based on a hybrid OpenMP/MPI paradigm. 


