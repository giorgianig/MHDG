#---------------------------------------
# project: MHDG
# Makefile
# date: 12/02/2016
#---------------------------------------
include Make.inc/arch.make

.PHONY : clean

all : allseq

allseq : MHDG

#allpar : checkmodepar plato geninitsol
#	@echo 'Makefile done'

default : MHDG

MHDG : MHDG-$(MDL)-$(MODE)-$(DIM)
	@echo ' '
	@echo '===================================================='
	@echo ' Executable done : $<'
	@echo '===================================================='
	@echo ' '

Convergence : Convergence-$(MDL)-$(DIM)
	@echo ' '
	@echo '===================================================='
	@echo ' Executable done : $<'
	@echo '===================================================='
	@echo ' '
	

clean: 
	@echo "Cleaning *.o *.mod *.F"
	$(RM) *.o *.mod *.kmo  *.F *.f90
	$(RM) MHDG-* Convergence *~

prec_const.F : $(SDIR)/Definitions/prec_const.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

prec_const.o : prec_const.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

types.F : $(SDIR)/Definitions/types.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

types.o : types.F prec_const.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

matrices_types.F : $(SDIR)/LinearAlgebra/matrices_types.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

matrices_types.o : matrices_types.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

globals.F : $(SDIR)/Definitions/globals.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

globals.o : globals.F types.o matrices_types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

MPI_OMP.F : $(SDIR)/MPI_OMP/MPI_OMP.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

MPI_OMP.o : MPI_OMP.F globals.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

magnetic_field.F : $(SDIR)/Models/magnetic_field.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

magnetic_field.o : magnetic_field.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

HDF5_io.F : $(SDIR)/InOut/HDF5_io.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

HDF5_io.o : HDF5_io.F prec_const.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

matrices_tools.F : $(SDIR)/LinearAlgebra/matrices_tools.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

matrices_tools.o : matrices_tools.F matrices_types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

PrintUtils.F : $(SDIR)/Utils/PrintUtils.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

PrintUtils.o : PrintUtils.F matrices_types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Communications.F : $(SDIR)/MPI_OMP/Communications.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Communications.o : Communications.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

interpolation.F : $(SDIR)/LinearAlgebra/interpolation.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

interpolation.o : interpolation.F PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adimensionalization.F : $(SDIR)/Models/adimensionalization.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adimensionalization.o : adimensionalization.F globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

LinearAlgebra.F : $(SDIR)/LinearAlgebra/LinearAlgebra.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

LinearAlgebra.o : LinearAlgebra.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
		
inout.F : $(SDIR)/InOut/inout.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

inout.o : inout.F HDF5_io.o globals.o PrintUtils.o LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
		
Debug.F : $(SDIR)/Utils/Debug.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Debug.o : Debug.F globals.o PrintUtils.o inout.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
		
physics.F : $(SDIR)/Models/$(RMDL)/physics.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

physics.o : physics.F globals.o interpolation.o HDF5_io.o magnetic_field.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

ifeq ($(PASTIX),$(LIB_YES))	
solve_pastix.F : $(SDIR)/LinearAlgebra/solve_pastix.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_pastix.o : solve_pastix.F matrices_types.o matrices_tools.o globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif

ifeq ($(PETSC),$(LIB_YES))	
solve_petsc.F : $(SDIR)/LinearAlgebra/solve_petsc.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@
	
solve_petsc.o : solve_petsc.F matrices_types.o matrices_tools.o globals.o types.o initialization.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif

ifeq ($(PSBLMG),$(LIB_YES))		
solve_psblas.F : $(SDIR)/LinearAlgebra/solve_psblas.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_psblas.o : solve_psblas.F matrices_types.o matrices_tools.o globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
else ifeq ($(PSBLAS),$(LIB_YES))
solve_psblas.F : $(SDIR)/LinearAlgebra/solve_psblas.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_psblas.o : solve_psblas.F matrices_types.o matrices_tools.o globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif
read_input.F : $(SDIR)/InOut/read_input.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

read_input.o : read_input.F prec_const.o globals.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
				
reference_element.F : $(SDIR)/HDG/reference_element.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

reference_element.o : reference_element.F globals.o LinearAlgebra.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
	
preprocess.F : $(SDIR)/HDG/preprocess.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

preprocess.o : preprocess.F globals.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
	
analytical.F : $(SDIR)/Models/$(RMDL)/analytical.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

analytical.o : analytical.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
	
initialization.F : $(SDIR)/HDG/initialization.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

initialization.o : initialization.F globals.o physics.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
	
hdg_PrecalculatedFirstEquation.F : $(SDIR)/HDG/hdg_PrecalculatedFirstEquation.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_PrecalculatedFirstEquation.o : hdg_PrecalculatedFirstEquation.F globals.o analytical.o LinearAlgebra.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

hdg_ComputeJacobian.F : $(SDIR)/HDG/hdg_ComputeJacobian.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_ComputeJacobian.o : hdg_ComputeJacobian.F globals.o analytical.o LinearAlgebra.o PrintUtils.o physics.o hdg_LimitingTechniques.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

hdg_Mapping.F : $(SDIR)/HDG/hdg_Mapping.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_Mapping.o : hdg_Mapping.F globals.o analytical.o LinearAlgebra.o PrintUtils.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	
	
hdg_BC.F : $(SDIR)/HDG/hdg_BC.f90 globals.o analytical.o LinearAlgebra.o PrintUtils.o
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_BC.o : hdg_BC.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

hdg_Assembly.F : $(SDIR)/HDG/hdg_Assembly.f90 globals.o analytical.o LinearAlgebra.o PrintUtils.o Debug.o
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_Assembly.o : hdg_Assembly.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

hdg_LimitingTechniques.F : $(SDIR)/HDG/hdg_LimitingTechniques.f90 globals.o  PrintUtils.o Debug.o
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_LimitingTechniques.o : hdg_LimitingTechniques.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	
	
solve_global_system.F : $(SDIR)/HDG/solve_global_system.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_global_system.o : solve_global_system.F globals.o analytical.o LinearAlgebra.o PrintUtils.o $(ADDMODSOLVER)
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

compute_element_solution.F : $(SDIR)/HDG/compute_element_solution.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

compute_element_solution.o : compute_element_solution.F globals.o analytical.o LinearAlgebra.o PrintUtils.o 
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@	

Postprocess.F : $(SDIR)/InOut/Postprocess.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Postprocess.o : Postprocess.F globals.o analytical.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

MHDG.F : $(SDIR)/MHDG.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@
	
MHDG.o : MHDG.F inout.o reference_element.o MPI_OMP.o PrintUtils.o Debug.o initialization.o Postprocess.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
	
Convergence.F : $(SDIR)/Convergence.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@
	
Convergence.o : Convergence.F inout.o reference_element.o MPI_OMP.o PrintUtils.o Debug.o initialization.o Postprocess.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@


MHDG-$(MDL)-$(MODE)-$(DIM) : MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o physics.o read_input.o HDF5_io.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o MHDG.o
	$(FC) $(FCFLAGS)   MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o physics.o read_input.o HDF5_io.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o MHDG.o -o MHDG-$(MDL)-$(MODE)-$(DIM) $(LIB)
	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test/
#	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_anal/
#	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_psblas/	
Convergence-$(MDL)-$(DIM) : MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o physics.o read_input.o HDF5_io.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o Convergence.o 	
	$(FC) $(FCFLAGS) $(OBJS) MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o physics.o read_input.o HDF5_io.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o Convergence.o -o Convergence-$(MDL)-$(MODE)-$(DIM) $(LIB)
	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test/
#	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test_anal/
