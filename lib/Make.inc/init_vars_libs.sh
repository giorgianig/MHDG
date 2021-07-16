#!/bin/bash

# This file must be run with the command: 
#> source name_of_this_file
# before compiling libraries and MHDG
# See arch.make for local installation typical path (all commented path with $(HOME)...)

#Stuff to add for INTEL (I) lib or OpenMP (OMP) lib
#export I_MPI_DEBUG=5
#export I_MPI_PIN_DOMAIN=omp
#export I_MPI_PIN_ORDER=bunch
#export KMP_AFFINITY=verbose,compact
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# HWLOC main directory ($HWLOC_DIR for mesocentre, $HWLOC_HOME for marconi)
export MHDG_HWLOC_DIR=/usr/lib/x86_64-linux-gnu

# HDF5 main directory ($HDF5_DIR for mesocentre, $HDF5_HOME for marconi)
export MHDG_HDF5_DIR=/usr/lib/x86_64-linux-gnu

# Set where the libraries are (if common directory)
export MHDG_LIB_DIR=/usr/local #$(pwd)/../libs

# Set the libraries
export MHDG_SCOTCH_DIR=$MHDG_LIB_DIR #/scotch_6.0.4
export MHDG_PASTIX_DIR=$MHDG_LIB_DIR #/pastix_5.2.3
export MHDG_PSBLAS_DIR=$MHDG_LIB_DIR #/psblas3
export MHDG_MLD2P4_DIR=$MHDG_LIB_DIR #/mld2p4-2

echo "Libraries directory: $MHDG_LIB_DIR"
echo "HWLOC directory: $MHDG_HWLOC_DIR"
echo "HDF5 directory: $MHDG_HDF5_DIR"
echo "SCOTCH (https://gitlab.inria.fr/scotch/scotch) directory: $MHDG_SCOTCH_DIR"
echo "PASTIX (https://gitlab.inria.fr/solverstack/pastix) directory: $MHDG_PASTIX_DIR"
echo "PSBLAS (https://github.com/sfilippone/psblas3) directory: $MHDG_PSBLAS_DIR"
echo "MLD2P4 (https://github.com/sfilippone/mld2p4-2) directory: $MHDG_MLD2P4_DIR"
