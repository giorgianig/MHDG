#-------------------------------------------------------------------------------
# Source directory
#-------------------------------------------------------------------------------
SDIR=$(PWD)/../src/

#-------------------------------------------------------------------------------
# Compiling type: optimized - debugging - profiling
#-------------------------------------------------------------------------------
COMPTYPE_OPT = opt
COMPTYPE_DEB = deb
COMPTYPE_PRO = pro
COMPTYPE = $(COMPTYPE_DEB)
#COMPTYPE = $(COMPTYPE_OPT)
#COMPTYPE = $(COMPTYPE_PRO)

#-------------------------------------------------------------------------------
# Mode: serial or parallel
#-------------------------------------------------------------------------------
MODE_SERIAL = serial
MODE_PARALL = parall
MODE = $(MODE_SERIAL)
#MODE = $(MODE_PARALL)

#-------------------------------------------------------------------------------
# The compiler
#-------------------------------------------------------------------------------
FC = mpif90

#-------------------------------------------------------------------------------
# Model
#-------------------------------------------------------------------------------
MDL_LAPLACE=Laplace
MDL_NGAMMA=NGamma
MDL_NGAMMATITE=NGammaTiTe
MDL_NGAMMAVORT=NGammaVort
# Model chosen
MDL=$(MDL_NGAMMA)
#MDL=$(MDL_NGAMMATITE)
#MDL=$(MDL_NGAMMAVORT)
#MDL=$(MDL_LAPLACE)

#-------------------------------------------------------------------------------
# Dimensions 2D/3D
#-------------------------------------------------------------------------------
DIM_3D=3D
DIM_2D=2D
DIM=$(DIM_2D)

#-------------------------------------------------------------------------------
# Libraries for linear system solver
#-------------------------------------------------------------------------------
# Available libraries: put $(LIB_YES) to use the library, $(LIB_NO) to not use it
LIB_YES=yes
LIB_NO=no
PASTIX=$(LIB_YES)
#PASTIX=$(LIB_NO)
PSBLAS=$(LIB_YES)
#PSBLAS=$(LIB_NO)
PSBLMG=$(LIB_YES)
#PSBLMG=$(LIB_NO)


# MACROS FOR MODEL CHOICE
ifeq ($(MDL),$(MDL_NGAMMA))
 MACROS+= -DNGAMMA
 RMDL=NGamma
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITE))
 RMDL=NGammaTiTe
 MACROS+= -DNGAMMA 
 MACROS+= -DTEMPERATURE
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMAVORT))
 RMDL=NGammaVort
 MACROS+= -DNGAMMA 
 MACROS+= -DVORTICITY
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_LAPLACE))
 RMDL=Laplace
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMALAPLACE))
 RMDL=NGammaLaplace
 MACROS+= -DNGAMMA
 ADDMOD+=hdg_LimitingTechniques.o 
else
 abort Unsupported MDL==$(MDL)
 exit
endif


# MACROS FOR SERIAL/PARALLEL
ifeq ($(MODE),$(MODE_SERIAL))
else ifeq ($(MODE),$(MODE_PARALL))
  MACROS+= -DPARALL
  ADDMOD+=Communications.o  
endif

ifeq ($(DIM),$(DIM_3D))
  MACROS+= -DTOR3D
endif

# MACROS FOR LINEAR SYSTEM LIBRARIES
ifeq ($(PASTIX),$(LIB_YES))
 MACROS+= -DWITH_PASTIX
 ADDMOD+=solve_pastix.o
endif
ifeq ($(PSBLMG),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 MACROS+= -DWITH_MLD2P4
 ADDMOD+=solve_psblas.o
else ifeq ($(PSBLAS),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 ADDMOD+=solve_psblas.o
endif

# MACROS FOR COMPILING OPTIONS
ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
# FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid
 FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,invalid
 FCFLAGS += -Wall -Wextra -Wconversion -fcheck=all -Wuninitialized 
else ifeq ($(COMPTYPE),$(COMPTYPE_PRO))
 FCFLAGS = -Og -pg
else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
 FCFLAGS = -Ofast
endif



FCFLAGS += -cpp  -fopenmp
FCFLAGS += -fdefault-double-8 -fdefault-real-8 
FCFLAGS += -ffree-line-length-none -fimplicit-none -ffree-form -Wno-tabs

DEF = -DTHREAD_FUNNELED
# Includes
#FCFLAGS += -I/usr/include
#FCFLAGS += -I/home/giorgio/2.5.2/linux-x86_64/include/hdf5/include
FCFLAGS += -I/usr/include/hdf5/serial/
FCFLAGS += -I/usr/include/X11
ifeq ($(PASTIX),$(LIB_YES))
FCFLAGS += -I/home/giorgio/libs/pastix_5.2.3_rep/install
FCFLAGS += -I/home/giorgio/libs/scotch_6.0.4/include
endif
ifeq ($(PSBLAS),$(LIB_YES))
FCFLAGS += -I/home/giorgio/libs/psblas3/include/
FCFLAGS += -I//home/giorgio/libs/psblas3/modules/
endif
ifeq ($(PSBLMG),$(LIB_YES))
FCFLAGS += -I//home/giorgio/libs/mld2p4-2/modules/
FCFLAGS += -I/home/giorgio/libs/mld2p4-2/include/
endif
# Libraries needed for linking
LIB += -L/usr/lib/x86_64-linux-gnu/
LIB += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
LIB += -lhdf5_fortran -lhdf5 -lz 
LIB += -lX11 -lXt 
ifeq ($(PASTIX),$(LIB_YES))
LIB += -L/home/giorgio/libs/pastix_5.2.3_rep/install -lpastix -lm -lrt
LIB += -L/home/giorgio/libs/scotch_6.0.4/lib/ -lscotch -lscotcherrexit  -lptscotchparmetis -lptscotch -lpthread -lhwloc
endif
LIB += -lblas -llapack
ifeq  ($(PSBLMG),$(LIB_YES)) 
LIB += -L/home/giorgio/libs/psblas3/lib/ -L/home/giorgio/libs/mld2p4-2/lib/ 
LIB += -lpsb_krylov -lmld_prec -lpsb_prec -lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base
else ifeq ($(PSBLAS),$(LIB_YES)) 
LIB += -L/home/giorgio/libs/psblas3/lib/
LIB += -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
endif

ifeq ($(PSBLMG),$(LIB_YES))
include /home/giorgio/libs/mld2p4-2/include/Make.inc.mld2p4
else ifeq ($(PSBLAS),$(LIB_YES))
include /home/giorgio/libs/psblas3/include/Make.inc.psblas
endif
# Add macros to the compiling flags
FCFLAGS += $(MACROS)
#-------------------------------------------------------------------------------
# Utilitaires
#-------------------------------------------------------------------------------

RM=/bin/rm -f
CP=/bin/cp -f
CMP=cmp


## MACROS FOR COMPILING OPTIONS
#ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
# FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid
# FCFLAGS += -Wall -Wextra -Wconversion -fcheck=all -Wuninitialized -Wtabs
#else ifeq ($(COMPTYPE),$(COMPTYPE_PRO))
# FCFLAGS = -Og -pg
#else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
# FCFLAGS = -Ofast
#endif


#FCFLAGS += -cpp  -fopenmp
#FCFLAGS += -fdefault-double-8 -fdefault-real-8 
#FCFLAGS += -ffree-line-length-none -fimplicit-none -ffree-form -Wno-tabs

#DEF = -DTHREAD_FUNNELED
## Includes
#FCFLAGS += -I/usr/include
##FCFLAGS += -I/home/giorgio/2.5.2/linux-x86_64/include/hdf5/include
#FCFLAGS += -I/usr/include/hdf5/serial/
#FCFLAGS += -I/usr/include/X11
#FCFLAGS += -I/home/giorgio/libs/pastix_5.2.3_rep/install
#FCFLAGS += -I/home/giorgio/libs/scotch_6.0.4/include

## Libraries needed for linking
#LIB =  -L/usr/lib/libblas/
#LIB += -L/usr/lib/lapack/
#LIB += -L/usr/lib/x86_64-linux-gnu/
#LIB += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/
#LIB += -lhdf5_fortran -lhdf5 -lz 
#LIB += -lX11 -lXt 
#LIB += -L/home/giorgio/libs/pastix_5.2.3_rep/install -lpastix -lm -lrt
#LIB += -L/home/giorgio/libs/scotch_6.0.4/lib/ -lscotch -lscotcherrexit  -lptscotchparmetis -lptscotch -lpthread -lhwloc
#LIB += -lblas -llapack


## Add macros to the compiling flags
#FCFLAGS += $(MACROS)
##-------------------------------------------------------------------------------
## Utilitaires
##-------------------------------------------------------------------------------

#RM=/bin/rm -f
#CP=/bin/cp -f
#CMP=cmp
