# SOLEDGE3X-HDG code

## Table of content
- [Introduction](#introduction)
    - [Models](#models-available)
    - [Dimensions](#dimensions)
    - [Parallelization](#parallelization)
- [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Use](#use)
    - [Examples](#examples)
- [Additional content](#additional-content)

## Introduction
MHDG or SOLEDGE3X-HDG is a 2D/3D electrostatic edge transport code for plasma simulations based on a high-order HDG (hybrid discontinuous Galerkin) scheme.
The code can run in realistic geometry and offer 3D magnetic perturbations under the vacuum approximation.
Electrostatic turbulence is a work in progress.

Related publications can be find [here](https://giorgianig.wixsite.com/giorgianig/publications).

### Models available
1. Laplace model: simple solver for anisotropic diffusion equation, 1 equation.

2. Isothermal model: model with 2 equations, unknowns density and parallel velocity.

3. Non-isothermal model: model with 4 equations, unknowns density, parallel velocity, ion temperature, electron temperature.

4. Isothermal + vorticity: model with 4 equations, unknowns density, parallel velocity, vorticity, electric potential
**Note:** Currently work in progress for thus last point.

### Dimensions

1- 2D cartesian and axisymmetric simulations.

2- 3D cartesian (slab geometry) and axisymmetric (toroidal geometry automatically generated from 2D mesh) simulations.

### Parallelization

Parallelization is based on a hybrid OpenMP/MPI paradigm. Domain decomposition in the plane + toroidal divisions

## Installation

### Prerequisites

An OpenMP/MPI installation is mandatory :p

MHDG uses the `hdf5` format. An installation of [hwloc](https://www.open-mpi.org/projects/hwloc/) (tested with HWLOC 2.0.4) is necessary

MHDG works with several libraries for the parallelization process:

1. [PaStiX](http://pastix.gforge.inria.fr/files/README-txt.html) (tested with PaStiX 5.2.3)

2. [PSBLAS](https://github.com/sfilippone/psblas3) (tested with last version)

3. [MLD2P-4](https://github.com/sfilippone/mld2p4-2) (tested with last version)

**Note:** not necessary (if not asked by the others libraries)

4. [SCOTCH](https://gitlab.inria.fr/scotch/scotch) (tested with SCOTCH 6.0.4)

**Note:** As a recommendation, you can install all the libraries in a `lib/libs/` directory. This could help for the compilation of the code in the following part. 

### Use

To compile MHDG, you go into `lib/Make.inc/` and check the `arch.make` file to set the correct path to each library.
It is possible to set the pathes in the `init_vars_lib_meso.sh` file, assuming each library is in the `lib/libs/` directory. Then you can do on a terminal (in the `lib/Make.inc/` directory:

```
    source init_vars_lib_meso.sh
```
Then don't forget to set the correct `arch.make` file in the `Makefile` in the `lib/` directory.

Finally you can compile the code with the command `make` in the `lib/` directory.

**Note:** the following option can be choosen in the `arch.make` file (uncomment):
- COMPTYPE: DEB for debbug, OPT for optimize, PRO for profiling. By default, the compiler is `ifort` (and related). The OPT/PRO option may be switched between `-O2` or `-O3` for the compilation.
- MDL: NGAMMA for isothermal model, NGAMMATITE for non-isothermal model, NGAMMAVORT for isothermal with vorticity model.
- DIM: 2D for the 2D version, 3D for 3D version of the code.

### Examples

In the `test/` directory you may try (after checking the `param.txt` file):

```
    ./MHDG-NGamma-serial-2D CircLimAlign_Quads_Nel208_P4
```
To start a OpenMP-only 2D N-Gamma simulation on a circular geometry with infinitely thin limiter.


```
    mpifort -n 4 ./MHDG-NGammaTiTe-serial-3D CircLimAlign_Quads_Nel208_P4
```
To start a OpenMP/MPI 3D N-Gamma-TiTe simulation on 4 processes, on a circular geometry with infinitely thin limiter.

**Note:** You can add a second argument after the mesh to do a restart from a former simulation (it is possible to do so from a 2D to a 3D simulation if the model of equations is the same).


## Additional content

The `matlab` files can be used for help:
- in `matlab/Preprocess`, meshes can be created with `GenerateHOMeshesCircGeomWithLim.m`.
- in `matlab/Preprocess`, meshes created with `GenerateHOMeshesCircGeomWithLim.m` or `[Gmsh](https://gmsh.info/)` can be separated for parallelization with `GenerateHDF5_meshes_Parallel.m`.
- in `matlab/MHD_NGammaTiTe/`, solution can be plotted with `plotFortranSolution.m` and similarly named files (depending on the model used).
