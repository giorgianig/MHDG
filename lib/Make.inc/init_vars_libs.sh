#!/bin/bash

# This file must be run with the command: 
#> source name_of_this_file
# before compiling libraries and MHDG
# See arch.make for local installation typical path (all commented path with $(HOME)...)

#OpenMP option necessary for parallel to avoid HWLOC warning
export OMP_NUM_THREADS=1

# HWLOC main directory
export MHDG_HWLOC_DIR=$HWLOC_DIR

# HDF5 main directory ($HDF5_DIR for mesocentre, $HDF5_HOME for marconi)
export MHDG_DF5_DIR=$HDF5_DIR

# Set where the libraries are (if common directory)
export MHDG_LIB_DIR=$(pwd)/../libs

# Set the libraries
export MHDG_SCOTCH_DIR=$MHDG_LIB_DIR/scotch_6.0.4
export MHDG_PASTIX_DIR=$MHDG_LIB_DIR/pastix_5.2.3
export MHDG_PSBLAS_DIR=$MHDG_LIB_DIR/psblas3
export MHDG_MLD2P4_DIR=$MHDG_LIB_DIR/mld2p4-2

echo "Libraries directory: $MHDG_LIB_DIR"
echo "HWLOC directory: $MHDG_HWLOC_DIR"
echo "HDF5 directory: $MHDG_HDF5_DIR"
echo "SCOTCH (https://gitlab.inria.fr/scotch/scotch) directory: $MHDG_SCOTCH_DIR"
echo "PASTIX (https://gitlab.inria.fr/solverstack/pastix) directory: $MHDG_PASTIX_DIR"
echo "PSBLAS (https://github.com/sfilippone/psblas3) directory: $MHDG_PSBLAS_DIR"
echo "MLD2P4 (https://github.com/sfilippone/mld2p4-2) directory: $MHDG_MLD2P4_DIR"
