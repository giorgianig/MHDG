#-------------------------------------------------------------------------------
# Source directory
#-------------------------------------------------------------------------------
SDIR=$(PWD)/../src/

#-------------------------------------------------------------------------------
# Compiling type: optimized or debugging
#-------------------------------------------------------------------------------
COMPTYPE_OPT = opt
COMPTYPE_DEB = deb
#COMPTYPE = $(COMPTYPE_DEB)
COMPTYPE = $(COMPTYPE_OPT)

#-------------------------------------------------------------------------------
# Mode: serial or parallel
#-------------------------------------------------------------------------------
MODE_SERIAL = serial
MODE_PARALL = parall
MODE = $(MODE_SERIAL)
MODE = $(MODE_PARALL)

#-------------------------------------------------------------------------------
# The compiler
#-------------------------------------------------------------------------------
FC = mpif90

#-------------------------------------------------------------------------------
# Model
#-------------------------------------------------------------------------------
# Available models
MDL_EULERISOTH=EulerIsoth
MDL_NGAMMA=NGamma

# Model chosen
#MDL=$(MDL_EULERISOTH)
MDL=$(MDL_NGAMMA)


ifeq ($(MDL),$(MDL_EULERISOTH))
 RMDL=EulerIsoth
 MACROS+= -DN3EQ
 MACROS+= -DEULISOTH
else ifeq ($(MDL),$(MDL_NGAMMA))
 RMDL=NGamma
 MACROS+= -DN2EQ
 MACROS+= -DNGAMMA
 ADDMOD+=hdg_LimitingTechniques.o
else
 abort Unsupported MDL==$(MDL)
 exit
endif

ifeq ($(MODE),$(MODE_SERIAL))
else ifeq ($(MODE),$(MODE_PARALL))
  MACROS+= -DPARALL
  ADDMOD+=Communications.o
endif
ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
 FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check 
 FCFLAGS += -Wall -Wextra -Wconversion -fcheck=all -Wuninitialized -Wtabs
else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
 FCFLAGS = -Ofast
endif



FCFLAGS += -cpp  -p -fopenmp
FCFLAGS += -fdefault-double-8 -fdefault-real-8 
FCFLAGS += -ffree-line-length-none -fimplicit-none -ffree-form

DEF = -DTHREAD_FUNNELED
# Includes
FCFLAGS += -I/usr/include
FCFLAGS += -I/home/giorgio/2.5.2/linux-x86_64/include/hdf5/include
FCFLAGS += -I/usr/include/X11
FCFLAGS += -I/home/giorgio/libs/pastix_5.2.3/install
FCFLAGS += -I/home/giorgio/libs/scotch_6.0.4/include

# Libraries needed for linking
LIB =  -L/usr/lib/libblas/
LIB += -L/usr/lib/lapack/
LIB += -L/usr/lib/x86_64-linux-gnu/
LIB += -lhdf5_fortran -lhdf5 -lz 
LIB += -lX11 -lXt 
LIB += -L/home/giorgio/libs/pastix_5.2.3/install -lpastix -lm -lrt
LIB += -L/home/giorgio/libs/scotch_6.0.4/lib/ -lscotch -lscotcherrexit  -lptscotchparmetis -lptscotch -lpthread
LIB += -lblas -llapack


# Add macros to the compiling flags
FCFLAGS += $(MACROS)
#-------------------------------------------------------------------------------
# Utilitaires
#-------------------------------------------------------------------------------

RM=/bin/rm -f
CP=/bin/cp -f
CMP=cmp
