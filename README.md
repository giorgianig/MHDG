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
    - [Mesh creation](#mesh-creation)
    - [Plot creation](#plot-creation)

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

An OpenMP/MPI installation is mandatory with either intel fortran compilers or GNU fortran compilers.
MHDG uses the `hdf5` format: the HDF5 library is necessary.
An installation of [hwloc](https://www.open-mpi.org/projects/hwloc/) (tested with HWLOC 2.0.4) is highly recommended (for PaStiX).

MHDG works with several libraries for inverting a sparse matrix in parallel:

1. [SCOTCH](https://gitlab.inria.fr/scotch/scotch) (tested with SCOTCH 6.0.4)

2. [PaStiX](http://pastix.gforge.inria.fr/files/README-txt.html) (tested with PaStiX 5.2.3)

3. [PSBLAS](https://github.com/sfilippone/psblas3) (tested with last version)

4. [MLD2P-4](https://github.com/sfilippone/mld2p4-2) (tested with last version)

**Note:** As a recommendation, you can install all the libraries in a `lib/libs/` directory. This could help for the compilation of the code in the following part.
**Note:** SCOTCH/PaStiX are for a direct sparse solver. PSBLAS/MLD2P-4 are for an iterative sparce solver. Only one of the two package is necessary.

### Use

To compile MHDG, you go into `lib/Make.inc/` and check the `arch.make` file to set the correct path to each library.
It is possible to set the pathes in the `init_vars_lib.sh` file, assuming each library is in the `lib/libs/` directory. Then you can do on a terminal (in the `lib/Make.inc/` directory):

```
    source init_vars_lib.sh
```
Then don't forget to set the correct `arch.make` file in the `Makefile` in the `lib/` directory. If you install the code on a local machine, 

Finally you can compile the code with the command `make` in the `lib/` directory.

**Note:** the following option can be choosen in the `arch.make` file (uncomment):
- COMPTYPE: DEB for debbug, OPT for optimize, PRO for profiling. By default, the compiler is `ifort` (and related). The OPT/PRO option may be switched between `-O2` or `-O3` for the compilation.
- MDL: NGAMMA for isothermal model, NGAMMATITE for non-isothermal model, NGAMMA(TITE)NEUTRAL for adding neutrals in (non-)isothermal model, NGAMMAVORT for isothermal with vorticity model.
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

**Note:** In the `param.txt` file, the WEST case is `testcase = 50`.

## Additional content

The `matlab` files can be used for help:
- in `matlab/Preprocess`, meshes can be created with either MATLAB only or Gmsh and MATLAB. (see related `README.md` in the directory).
- in `matlab/Postprocess/`, solutions can be plotted with `plotFortranSolution.m`.

**Very important note:** Once in the corresponding matlab directory, always execute first in the matlab console:
```
    >> setpath
```
**TODO**
- [] 1D plots (profiles)
- [] 3D plots
- [] Write the routine in python for better portability

### Mesh creation
Go into the `matlab/Preprocess` directory.
Meshes are saved in both `.mat` and `.h5` format but only the `.h5` is useful for SOLEDGE3X-HDG.

#### With MATLAB only
Use GenerateHOMeshesCircGeomWithLim.m to create a circular mesh with limiter

#### With Gmsh and MATLAB
Use [Gmsh](https://gmsh.info/) for any other mesh (example below)

1. Load a .geo file in Gmsh/ (you can create your own but be careful with the number of `IN`(1), `OUT(2)`, `DOM`(4), `LIM`(3) conditions)
2. On the left menu, unroll `Mesh` -> `Define` then click on `1D`, then `2D`
3. In the menu `Tools` -> `Options` go to the `Mesh` menu and `General` option.
    1. 2D algorithm: choose either `MeshAdapt` (recommended) or `Delaunay` or `Frontal-Delaunay`
    2. Element order: choose the polynomial order (generally `4` or `6`)
    3. If you want to refine the mesh, decrease `Element size factor` (esf) with 0 < esf < 1 OR use directly in the .geo the 4th value on each point
**Note:** `West_Mesh_YesHole_coarse.geo` can use the `lcX` variables to easily set the element size factor.
    4. If you modify the .geo file, reload it. Left menu: unroll `Geometry` -> `Physical groups` -> `Reload script`.
    5. Then do step 2.b again.
    6. Before saving, check the orientation of the mesh (normal toward you) in `Tools` -> `Options` -> `Geometry` -> `Visibility` and set Normals and tangents > 0.
        - If the normal is coming toward you, then the mesh should be correct for MHDG
        - Else, check the Curve Loop and change the signs.
    7. Save the file: `File` -> `Export`
        1. Save your file as a .msh
        2. Choose `ASCII 2` format
        3. Uncheck all the box
        4. Save
    8. Open with  matlab: `convertGmeshMeshToMat.m` (if triangle) or `convertGmeshMeshQuadsToMat.m` (if ONLY quads, not recommended for a first try)
    9. Set the filename and polynomial order (`P`) correctly and run (F5) to obtain a .mat and .h5 mesh file

#### Parallel mesh
1. Open DomainDecomposition.m
2. Set the meshName (.mat or .h5)
3. Choose the number of division div
4. Choose ONE of the `divType`: `'rot'` (recommended), `'x'`, `'y'`
5. Set `theta0` to shift the starting point (for `'rot'`)
6. Set `elemType` (1:triangle, 0:Quads, BE CAREFUL, it is the opposite in the mesh file...)
7. Set forbidden values with `theta_forb` (position of the limiter for example)
8. Run (F5)
**Note:** If the code loops infinitely, you can go to the line `Balancing divisions` and fine tune the if condition `(nel-nel_mean) >/< +/-var` by changing `var`.


### Plot creation
Go into the `matlab/Postprocess` directory and play with `plotFortranSolution.m`. The file is already commented.
