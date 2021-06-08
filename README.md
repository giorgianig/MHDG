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

1. [SCOTCH](https://gitlab.inria.fr/scotch/scotch) (tested with SCOTCH 6.0.4)

2. [PaStiX](http://pastix.gforge.inria.fr/files/README-txt.html) (tested with PaStiX 5.2.3)

3. [PSBLAS](https://github.com/sfilippone/psblas3) (tested with last version)

4. [MLD2P-4](https://github.com/sfilippone/mld2p4-2) (tested with last version)


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
The syntax to run the code is in serial (OpenMP only):

```
    ./executable* mesh* restart
```
**Note:** name with a `*` are mandatory.

- The executable name is always of the form `MHDG-[model]-[parallelization]-[dimension]`
- The mesh name is of the form `mesh_name.h5` in serial or `mesh_name_i_N.h5` in parallel (with `i = 1...N`). You must always put `mesh*` as `mesh_name`.
- The restart name must always be of the form `Sol[dimension]_mesh-name_DPe[value]_DPai[value]_DPae[value]` (following part, such as the extension must be ignore in the same manner as for the mesh).

In the `test/` directory you can run:

```
    ./MHDG-NGamma-serial-2D CircLimAlign_Quads_Nel208_P4
```
To start a OpenMP-only 2D N-Gamma simulation on a circular geometry with infinitely thin limiter.

For a parallel simulation (MPI/OpenMP), you must create a mesh subdivided in `N = number of task` part (see Additional content). 
As an example for a 4 part mesh, you can run in the `test/` directory:
```
    mpiifort -n 4 ./MHDG-NGammaTiTe-serial-3D CircLimAlign_Quads_Nel208_P4
```
To start an OpenMP/MPI 3D N-Gamma-TiTe simulation on 4 processes, on a circular geometry with infinitely thin limiter.

**Note:** This simulation uses the mesh files `CircLimAlign_Quads_Nel208_P4_i_4.h5` with `i = 1...4`.

In the `test/West/` directory, you may run a WEST simulation either in serial or on 16 tasks. Additional files are necessary (and present) such as the magnetic field `WEST_far_465.h5`or the position of the elements in a node (`positionFeketeNodesTri2D.h5`).

**Note:** In the param file, the WEST case is `testcase = 50`.

## Additional content

The `matlab` files can be used for help:
- in `matlab/Preprocess`, meshes can be created with `GenerateHOMeshesCircGeomWithLim.m`.
- in `matlab/Preprocess`, meshes created with `GenerateHOMeshesCircGeomWithLim.m` or `[Gmsh](https://gmsh.info/)` can be separated for parallelization with `GenerateHDF5_meshes_Parallel.m`.
- in `matlab/MHD_NGammaTiTe/`, solution can be plotted with `plotFortranSolution.m` and similarly named files (depending on the model used).
