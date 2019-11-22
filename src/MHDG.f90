PROGRAM MHDG
 USE in_out
 USE reference_element
 USE preprocess
 USE MPI_OMP
 USE printutils
 USE debug
 USE initialization
#ifdef WITH_PASTIX 
 USE solve_pastix
#endif
#ifdef WITH_PSBLAS
 USE solve_psblas
#endif 
 USE Postprocess 
#ifdef PARALL
 USE Communications
#endif 
 USE HDG_LimitingTechniques 
 IMPLICIT NONE

 integer             :: Neq,Np,Nel,Ndim,Nfp,Nf,ntorloc
 integer             :: Nfl,Np1dPol,Np1dTor,Ng1dPol,Ng1dTor,Nthreads
 integer             :: it,ir,ifa,nts,nu,nut,nb_args,IERR,k,i,j,g
 integer,allocatable :: faceNodes1d(:)
 real*8              :: dt, errNR, errTime
 real*8              :: time_start, time_finish, tps, tpe
 real*8, allocatable :: uiter(:),u0(:),um2(:),L2err(:)
 character(LEN = 1024) :: mesh_name,namemat,save_name
 integer             :: cks,clock_rate,cke,clock_start,clock_end
 real*8              :: cputpre,cputmap,cputass,cputbcd,cputsol,cputcon,cputpdf,cputglb,cputtot
 real*8              :: runtpre,runtmap,runtass,runtbcd,runtsol,runtcon,runtpdf,runtglb,runttot
 integer             ::  OMP_GET_MAX_THREADS

write(6,* ) "STARTING"
 
  ! Check the number of input arguments
  nb_args = iargc()

  IF (nb_args .lt. 1) THEN
      PRINT *, " Error: the mesh name is needed"
      stop
  END IF
  
 CALL getarg(1, mesh_name)
 mesh_name = Adjustl(Trim(mesh_name))
 
 IF (nb_args .gt. 1) THEN
      CALL getarg(2, save_name)
      save_name = Adjustl(Trim(save_name))
      PRINT *, " Restart simulation with solution: ", save_name
 END IF     
 
 IF (nb_args .gt. 2) THEN
      PRINT *, " Too many arguments "
      stop
 END IF  
 
 ! Number of threads
 Nthreads = OMP_GET_MAX_THREADS()
 
! Start timing the code
 call cpu_time(time_start)
 call system_clock(clock_start,clock_rate)
 
 ! Initialize MPI
 CALL init_MPI_OMP()
 
 ! Read input file param.txt
 CALL read_input()

 ! Set parallelization division in toroidal and poloidal plane
#ifdef TOR3D
#ifdef PARALL
  MPIvar%ntor = numer%npartor
  IF (MPIvar%glob_size.lt.MPIvar%ntor) THEN
     WRITE(6,*) "Error: wrong number of MPI toroidal partition in input file or wrong number of MPI processes"
     WRITE(6,*) "Number of processes: ",MPIvar%glob_size
     WRITE(6,*) "Number of MPI toroidal partitions: ",MPIvar%ntor
     STOP
  ENDIF
  IF (mod(MPIvar%glob_size,MPIvar%ntor).ne.0) THEN
     WRITE(6,*) "Error: the number of MPI processes must be a multiple of the number of MPI toroidal divisions (set in input file)"
     WRITE(6,*) "Number of processes: ",MPIvar%glob_size
     WRITE(6,*) "Number of MPI toroidal partitions: ",MPIvar%ntor     
     STOP
  ENDIF

  MPIvar%npol = MPIvar%glob_size/MPIvar%ntor
  MPIvar%itor = 1+MPIvar%glob_id/MPIvar%npol
  MPIvar%ipol = 1+mod(MPIvar%glob_id,MPIvar%npol)
#endif
#endif

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!write(6,*) "Process: ", MPIvar%glob_id,"MPIvar%ntor",MPIvar%ntor
!write(6,*) "Process: ", MPIvar%glob_id,"MPIvar%npol",MPIvar%npol 
!write(6,*) "Process: ", MPIvar%glob_id,"MPIvar%itor",MPIvar%itor 
!write(6,*) "Process: ", MPIvar%glob_id,"MPIvar%ipol",MPIvar%ipol 
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!stop
 
 
 
 ! Load the mesh file
 CALL load_mesh(mesh_name)

 ! Linear solver: set the start to true
 matK%start = .true.

#ifndef MOVINGEQUILIBRIUM
 IF ( (switch%testcase .ge. 50) .and. (switch%testcase .le. 59) ) THEN
 ! Testcase is on the WEST mesh -> load magnetic field
	   CALL loadMagneticField()
 END IF 
#endif 

 ! Create the reference element based on the mesh type
 CALL create_reference_element(refElPol,2) 
#ifdef TOR3D
 !*************************************************
 !                REFERENCE ELEMENT 3D
 !*************************************************
 ! Create the reference element for the toroidal interpolation
 
 CALL create_reference_element(refElTor,1,numer%ptor)

  !***********************************************
  ! Generate 3D basis functions
  !*********************************************** 
  ALLOCATE(refElTor%N3D(refElPol%NGauss2D*refElTor%NGauss1D,refElPol%Nnodes2D*refElTor%Nnodes1D) ) 
  DO i = 1,refElTor%NGauss1D
     DO j = 1,refElPol%NGauss2D
        g = (i-1)*refElPol%NGauss2D+j
        refElTor%N3D(g,:) = col(TensorProduct(refElPol%N2D(j,:),refElTor%N1D(i,:)))
     END DO
  END DO


  !***********************************************
  ! Generate face nodes for toroidal faces
  !***********************************************  
  ALLOCATE(refElTor%faceNodes3(refElPol%Nfaces,size(refElPol%face_nodes,2)*size(refElTor%face_nodes,2)))
  ALLOCATE(faceNodes1d(refElTor%Nnodes1D))
		DO ifa = 1,refElPol%Nfaces
		    faceNodes1d =  (/ (i, i = 1, refElTor%Nnodes1D) /)
					 refElTor%faceNodes3(ifa,:) = colint(TensorSumInt(refElPol%face_nodes(ifa,:), (faceNodes1d-1)*refElPol%Nnodes2D))
		END DO
		DEALLOCATE(faceNodes1d)

  !***********************************************
  ! Generate shape functions for toroidal faces
  !***********************************************  
  Np1dPol = size(refElPol%N1d,2)
  Np1dTor = size(refElTor%N1d,2)
  Ng1dPol = size(refElPol%N1d,1)
  Ng1dTor = size(refElTor%N1d,1)
  refElTor%Nfl = Np1dPol*Np1dTor ! Number of nodes for lateral faces
  refElTor%Ngl = Ng1dPol*Ng1dTor ! Number of Gauss points for lateral faces
  refElTor%Nnodes3D = refElPol%Nnodes2D*refElTor%Nnodes1D ! Number of nodes in the 3D element
  refElTor%NGauss3D = refElPol%NGauss2D*refElTor%NGauss1D ! Number of Gauss points in the 3D element
  refElTor%Nft = refElTor%Nfl*refElPol%Nfaces+2*refElPol%Nnodes2d         ! Number of nodes in all the faces of the 3D element
  ALLOCATE(refElTor%sFTF(Ng1dPol*Ng1dTor,Np1dPol*Np1dTor))
  DO i = 1,Np1dTor
      DO j=1,Np1dPol
          k = (i-1)*Np1dPol+j
          refElTor%sFTF(:,k) = col(tensorProduct(refElPol%N1d(:,j),refElTor%N1d(:,i)))
      END DO
  END DO
#endif 
 ! Mesh preprocess: create the mesh related structures 
 ! used in the HDG scheme
 CALL mesh_preprocess()
 
  ! Initialization of the simulation parameters !TODO check if I need to pass dt to init_sim
 CALL init_sim(nts,dt)
 
 
 ! Variables for timing purposes 
 cputpre = 0.
 cputmap = 0.
 cputass = 0.
 cputbcd = 0.
 cputsol = 0.
 cputcon = 0.
 cputpdf = 0.
 cputglb = 0.
 cputtot = 0.

 runtpre = 0.
 runtmap = 0.
 runtass = 0.
 runtbcd = 0.
 runtsol = 0.
 runtcon = 0.
 runtpdf = 0. 
 runtglb = 0.
 runttot = 0. 
  
   
#ifdef TOR3D  
  Ndim        = 3                                               ! Number of dimensions
#ifdef PARALL
  ntorloc = numer%ntor/MPIvar%ntor
#else
  ntorloc = numer%ntor
#endif  
  Nel         = Mesh%Nelems*ntorloc                             ! Number of 3D elements
  Np          = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
  Nfl         = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
  Nfp         = refElPol%Nnodes2D*2+refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
  Nf          = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
#ifdef PARALL
  if (MPIvar%ntor.gt.1) then
     nut         = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
  else
     nut         = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
  endif
#else  
  nut         = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
#endif  
  Nfp         = refElPol%Nnodes2D*2+refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
  nu          = phys%Neq*Nel*Np
#else
  Ndim        = Mesh%ndim
  Nel         = Mesh%Nelems
  Np          = refElPol%Nnodes2D
  Nf          = refElPol%Nfaces 
  Nfp         = refElPol%Nfacenodes*Nf 
  nut         = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
  nu          = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif   

#ifdef PARALL
  CALL init_com()
#endif
  
 ! Allocation and initialization of the elemental matrices
 CALL init_elmat()
 
 ! Compute the elemental matrices that depend only on the 
 ! mesh, i.e. fixed in time
	if (utils%timing) then
		call cpu_time(tps)					
		call system_clock(cks,clock_rate)
	end if
 CALL hdg_precalculatedMatrices()
	if (utils%timing) then
	 call cpu_time(tpe)
	 call system_clock(cke,clock_rate)
	 runtpre = (cke-cks)/real(clock_rate)
	 cputpre = tpe-tps	 
!  write(6,*) "Precalculated matrices cpu-time: ", tpe-tps
!  write(6,*) "Precalculated matrices run-time: ", (cke-cks)/real(clock_rate)
	end if

 ! Add diffusion in corners
 if (switch%difcor.gt.0) then
!				if (utils%timing) then
!					call cpu_time(tps)					
!					call system_clock(cks,clock_rate)
!				end if
				CALL HDG_AddDiffusionCorner()! TODO : program for the 3D
!				if (utils%timing) then
!					call cpu_time(tpe)
!					call system_clock(cke,clock_rate)
!					write(6,*) "Add diffusion in corners cpu-time: ", tpe-tps
!					write(6,*) "Add diffusion in corners run-time: ", (cke-cks)/real(clock_rate)
!				end if  
	end if
	 ! Initialize limiting of rho
	 if (switch%limrho.gt.0) then
	    CALL initializeLimiting() ! TODO : program for the 3D
	 endif
	 ! Initialize shock capturing
	 if (switch%shockcp.gt.0) then
	    CALL initializeShockCapturing() ! TODO : program for the 3D
	 endif	 

 ! Initialize the solution
  IF (nb_args.eq.2) THEN
       ! restart simulation: load solution from file (the name is given in argument)
    	  CALL HDF5_load_solution(save_name)  
  ELSE
					  CALL init_sol()
  END IF


 ! Save solution
 CALL setSolName(save_name,mesh_name,0,.false.,.false.) 
 CALL HDF5_save_solution(save_name)

 ! Allocate and initialize uiter and u0
 ALLOCATE(uiter(nu))
 ALLOCATE(u0(nu))  
 if (time%tis==2) then
    ALLOCATE(um2(nu))
 else
    ALLOCATE(um2(1:2)) ! to avoid errors in debug mode..
 endif
 u0 = sol%u
 uiter = 0.
 
 !*******************************************************
 !                  TIME LOOP
 !*******************************************************
 DO it = 1,nts ! ************ TIME LOOP *********************
 
    ! Compute time step
    
    ! Actualization of time 
    time%t  = time%t + time%dt
		  time%it = time%it + 1 
		  time%ik = time%ik + 1
  		
  		IF (MPIvar%glob_id.eq.0) THEN
							WRITE(6, '(" *", 60("*"), "**")') 
						 WRITE(6, '(" *", 20X,    "Time iteration   = ", I5, 16X, " *")') time%it
							WRITE(6, '(" *", 60("*"), "**")') 
  		END IF
#ifdef MOVINGEQUILIBRIUM
    !*******************************************************************
    !                      MOVING EQUILIBRIUM SECTION
    !*******************************************************************
    ! Load a time variable magnetic field
				IF ( (switch%testcase .ge. 50) .and. (switch%testcase .le. 59) ) THEN
				! Testcase is on the WEST mesh -> load magnetic field
							CALL loadMagneticFieldTemporalEvolution()
				ELSE
				   WRITE(6,*) "Time variable magnetic field only available for West case"
				   STOP
				END IF   

				! Compute the elemental matrices that depends on the magnetic field
!				if (utils%timing) then
!					call cpu_time(tps)					
!					call system_clock(cks,clock_rate)
!				end if
				CALL hdg_MagneticDependingMatrices() ! TODO : program for the 3D
!				if (utils%timing) then
!				 call cpu_time(tpe)
!				 call system_clock(cke,clock_rate)
!     write(6,*) "Magnetic matrices cpu-time: ", tpe-tps
!     write(6,*) "Magnetic matrices run-time: ", (cke-cks)/real(clock_rate)
!				end if						 
#endif

				! Limiting of rho
				if (switch%limrho.gt.0) then
!							if (utils%timing) then
!								call cpu_time(tps)					
!								call system_clock(cks,clock_rate)
!							end if
							CALL HDG_LimitingRho() ! TODO : program for the 3D
!							if (utils%timing) then
!								call cpu_time(tpe)
!								call system_clock(cke,clock_rate)
!								write(6,*) "Limiting rho cpu-time: ", tpe-tps
!								write(6,*) "Limiting rho run-time: ", (cke-cks)/real(clock_rate)
!							end if  
				end if  		

			
				
!											    CALL setSolName(save_name,mesh_name,0,.false.) 
!															CALL HDF5_save_solution(save_name)
																			
									!*******************************************************
									!             Newton-Raphson iterations
									!*******************************************************
									uiter = u0
									DO ir = 1,numer%nrp ! ************ NEWTON-RAPHSON LOOP *********************
         
             IF (MPIvar%glob_id.eq.0) THEN
                WRITE(6,*) "***** NR iteration: ", ir, "*****"
									    ENDIF
									    
													! Shock capturing
													if (switch%shockcp.gt.0) then
!																if (utils%timing) then
!																	call cpu_time(tps)					
!																	call system_clock(cks,clock_rate)
!																end if
																CALL HDG_ShockCapturing() ! TODO : program for the 3D
!																if (utils%timing) then
!																	call cpu_time(tpe)
!																	call system_clock(cke,clock_rate)
!																	write(6,*) "Shock capturing cpu-time: ", tpe-tps
!																	write(6,*) "Shock capturing run-time: ", (cke-cks)/real(clock_rate)
!																end if  
													end if	
													    
   									! **** Compute the elemental matrices that depends on the N-R iteration ****
   									!***************************************************************************
   									! Convection matrices
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CALL hdg_ConvectionMatrices()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)
             runtcon = runtcon+(cke-cks)/real(clock_rate)
             cputcon = cputcon+tpe-tps   									 
!             write(6,*) "Convection cpu-time: ", tpe-tps
!             write(6,*) "Convection run-time: ", (cke-cks)/real(clock_rate)
   									end if	   		
   									
            
            
            
!if (ir==2) then
!call saveArray(elMat%Cv,'Cv')
!call saveArray(elMat%H,'H')
!call saveArray(elMat%Hf,'Hf')
!call saveArray(elMat%D,'D')
!call saveArray(elMat%E,'E')
!call saveArray(elMat%Ef,'Ef')
!call saveArray(elMat%Df,'Df')
!write(6,*) "Done saving"
!stop
!endif             

   									
!   									call saveMatrix(elMat%D(:,:,1),'D')
!   									call saveMatrix(elMat%E(:,:,1),'E')
!   									call saveMatrix(elMat%Df(:,:,1),'Df')
!   									call saveMatrix(elMat%Ef(:,:,1),'Ef')
!   									stop
 
#ifdef TEMPERATURE
   									! Parallel diffusion matrices
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CALL HDG_paralldiffmatrices()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)
!             write(6,*) "Parallel diffusion cpu-time: ", tpe-tps
!             write(6,*) "Parallel diffusion run-time: ", (cke-cks)/real(clock_rate)
             runtpdf = runtpdf+(cke-cks)/real(clock_rate)
             cputpdf = cputpdf+tpe-tps  
   									end if	   									
#endif
 
   									! Set boundary conditions
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if								
   									CALL hdg_BC()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)
             runtbcd = runtbcd+(cke-cks)/real(clock_rate)
             cputbcd = cputbcd+tpe-tps     									 
!             write(6,*) "Boundary conditions cpu-time: ", tpe-tps
!             write(6,*) "Boundary conditions run-time: ", (cke-cks)/real(clock_rate)
   									end if	 
   									  									  
   									  									  
   									  									   									
   									! Compute elemental mapping
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CAlL hdg_Mapping(u0,um2)
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)   									 
             runtmap = runtmap+(cke-cks)/real(clock_rate)
             cputmap = cputmap+tpe-tps   									 
!             write(6,*) "Mapping cpu-time: ", tpe-tps
!             write(6,*) "Mapping run-time: ", (cke-cks)/real(clock_rate)
   									end if   									



!WRITE(6,*) "LL: "
!CALL findNaNArr(elMat%LL)
!WRITE(6,*) "L0: "
!CALL findNaNMat(elMat%L0)
!WRITE(6,*) "UU: "
!CALL findNaNArr(elMat%UU)
!WRITE(6,*) "U0: "
!CALL findNaNMat(elMat%U0)

!call saveMatrix(elMat%Df(:,:,5637),'Df')
!call saveMatrix(elMat%Hf(:,:,5637),'Hf')
!call saveMatrix(elMat%Ef(:,:,5637),'Ef')
!call saveMatrix(elMat%Lf(:,:,5637),'Lf')
!call saveMatrix(elMat%Qf(:,:,5637),'Qf')
!call saveMatrix(elMat%TQhf(:,:,5637),'TQhf')
!stop

!call print_matrices_tex()
!stop

   									! Assembly the global matrix
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CALL hdg_Assembly()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)   									 
!		           write(6,*) "Assembly cpu-time: ", tpe-tps
!		           write(6,*) "Assembly run-time: ", (cke-cks)/real(clock_rate)
              runtass = runtass+(cke-cks)/real(clock_rate)
              cputass = cputass+tpe-tps
   									end if  
        
        
            
!call print_matrices_hdf5()
!call HDF5_save_CSR_matrix('matK')
!call HDF5_save_vector(rhs%vals,'rhs')
!stop        
        
!            call HDF5_save_matrix('HDG_mat')
!            stop

!if (ir==2) then
!write(namemat,'(A,I0,A,I0)') "MatK_", MPIvar%glob_id,"_",MPIvar%glob_size            
!call save_CSR_matrix_txt(matK,trim(adjustl(namemat)) )
!write(namemat,'(A,I0,A,I0)') "rhs_", MPIvar%glob_id,"_",MPIvar%glob_size 
!call save_CSR_vector_txt(rhs,trim(adjustl(namemat)))
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop
!end if


!write(namemat,'(A,I0,A,I0)') "utilde_", MPIvar%glob_id,"_",MPIvar%glob_size
!call saveVector(sol%u_tilde,trim(adjustl(namemat)))
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop

!if (ir==2) then



!#ifdef PARALL
!write(namemat,'(A,I0,A,I0)') "MatK_", MPIvar%glob_id,"_",MPIvar%glob_size            
!call save_CSR_matrix_txt(matK,trim(adjustl(namemat)) )
!write(namemat,'(A,I0,A,I0)') "rhs_", MPIvar%glob_id,"_",MPIvar%glob_size 
!call save_CSR_vector_txt(rhs,trim(adjustl(namemat)))
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop
!#else
!call save_CSR_matrix_txt(matK, 'MatK_seq' )
!call save_CSR_vector_txt(rhs,'rhs_seq')
!stop
!#endif

!#ifdef TOR3D
!call save_CSR_matrix_txt(matK, 'MatK_3D' )
!call save_CSR_vector_txt(rhs,'rhs_3D')
!stop
!#else
!call save_CSR_matrix_txt(matK, 'MatK_2D' )
!call save_CSR_vector_txt(rhs,'rhs_2D')
!stop
!#endif





!endif

!call saveArray(elMat%Lf,'Lf')
!call saveArray(elMat%TQhf, 'TQhf')
!stop

!write(6,*) "exting"
!call print_matrices_tex()
!call save_CSR_matrix_txt(matK,'MatK')
!call saveVector(rhs%vals,'f')
!stop


								! Solve linear system, compute face and element solution
!   									if (utils%timing) then
!    									call cpu_time(tps)					
!    									call system_clock(cks,clock_rate)
!   									end if
!   									CALL compute_solution()
!   									if (utils%timing) then
!   									 call cpu_time(tpe)
!   									 call system_clock(cke,clock_rate)   									 
!!		           write(6,*) "Linear system solve cpu-time: ", tpe-tps
!!		           write(6,*) "Linear system solve run-time: ", (cke-cks)/real(clock_rate)
!              runtsol = runtsol+(cke-cks)/real(clock_rate)
!              cputsol = cputsol+tpe-tps
!   									end if  	
   									
   									

            ! Solve linear system
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CALL solve_global_system()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)   									 
!		           write(6,*) "Linear system solve cpu-time: ", tpe-tps
!		           write(6,*) "Linear system solve run-time: ", (cke-cks)/real(clock_rate)
              runtsol = runtsol+(cke-cks)/real(clock_rate)
              cputsol = cputsol+tpe-tps
   									end if  




            ! Compute element-by-element solution
   									if (utils%timing) then
    									call cpu_time(tps)					
    									call system_clock(cks,clock_rate)
   									end if
   									CALL compute_element_solution()
   									if (utils%timing) then
   									 call cpu_time(tpe)
   									 call system_clock(cke,clock_rate)   									 
!		           write(6,*) "Linear system solve cpu-time: ", tpe-tps
!		           write(6,*) "Linear system solve run-time: ", (cke-cks)/real(clock_rate)
              runtglb = runtglb+(cke-cks)/real(clock_rate)
              cputglb = cputglb+tpe-tps
   									end if  

  									
!write(namemat,'(A,I0,A,I0)') "utilde_", MPIvar%glob_id,"_",MPIvar%glob_size
!call saveVector(sol%u_tilde,trim(adjustl(namemat)))
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop   									

!write(namemat,'(A,I0,A,I0)') "u_", MPIvar%glob_id,"_",MPIvar%glob_size
!call saveVector(sol%u,trim(adjustl(namemat)))
!call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!stop  


!call HDF5_save_vector(sol%u_tilde,'utilde')
!stop     									
								
!if (ir==15) then   									
!call print_matrices_singleElem_tex(1902)
!call saveVector(sol%q,'q')
!call saveVector(sol%u,'solu')
!call saveVector(sol%u_tilde,'solutilde')
!call saveVector(u0,'u0')
!stop
!endif
   									
												! Check for NaN (doesn't work with optimization flags)
												DO i=1,nu
														IF (sol%u(i) /= sol%u(i)) THEN
																		WRITE(6,*) "NaN detected"
																		STOP
															END IF
												END DO
  				   									
												! Apply threshold 
												CALL HDG_applyThreshold()	
!call saveVector(sol%u,'solu_afterThresh')
!call saveVector(sol%u_tilde,'solutilde_afterThresh')												
												! Apply threshold 
												CALL HDG_FilterSolution()													
!call saveVector(sol%u,'solu_afterFilter')
!call saveVector(sol%u_tilde,'solutilde_afterFilter')
!stop
!write(6,*) sol%u

!stop
											    ! Save solution
											    IF (switch%saveNR) THEN
   											    CALL setSolName(save_name,mesh_name,ir,.false.,.true.) 
			   												CALL HDF5_save_solution(save_name)
		             END IF

   									! Check convergence of Newton-Raphson
  									errNR = computeResidual(sol%u,uiter,1.)
  									errNR = errNR/numer%dumpnr
  									
  									 IF (MPIvar%glob_id.eq.0) THEN
   									   WRITE(6,*) "Error: ", errNR
   									END IF
   									IF (errNR.lt.numer%tNR) THEN
   									   EXIT
   									ELSEIF (errNR.gt.numer%div) THEN
   									   WRITE(6,*) 'Problem in the N-R procedure'
   									   STOP
   									ELSE
   									   uiter = sol%u
   									END IF
   									IF (MPIvar%glob_id.eq.0) THEN
               WRITE(6,*) "*********************************"			
               WRITE(6,*) " "
               WRITE(6,*) " "
            ENDIF
            
									END DO ! ************ END OF NEWTON-RAPHSON LOOP *********************
									
!							! Apply threshold 
!							CALL HDG_applyThreshold()	



									
	 				! Check convergence in time advancing and update
  				errTime = computeResidual(sol%u,u0,time%dt)
  				sol%tres(sol%Nt+1) = errTime
  				sol%time(sol%Nt+1) = time%t
  				sol%Nt       = sol%Nt + 1
  				
  				! Display results
  				IF (mod(time%it,utils%freqdisp).eq. 0) THEN
			   			CALL displayResults()
			   END IF
  				
  				! Check for NaN (doesn't work with optimization flags)
  				DO i=1,nu
									IF (sol%u(i) /= sol%u(i)) THEN
										  WRITE(6,*) "NaN detected"
										  STOP
									END IF
  				END DO
  				

			   
			   IF (.not.switch%steady) THEN
			     				! Save solution
				IF (.not.switch%steady) THEN
  				IF (mod(time%it,utils%freqsave).eq. 0) THEN
         CALL setSolName(save_name,mesh_name,time%it,.true.,.false.) 
			   			CALL HDF5_save_solution(save_name)
			   END IF
			   END IF
			   
									!****************************************
									! Check steady state and update or exit
									!****************************************
									IF (errTime.lt.numer%tTM) THEN
												IF (switch%psdtime) THEN
											 ! Pseudo-time simulation
											 
											    ! Save solution
											    CALL setSolName(save_name,mesh_name,it,.true.,.true.) 
															CALL HDF5_save_solution(save_name)
												
															! Update the diffusion, the elemental matrices and the solution
															IF (MPIvar%glob_id.eq.0) THEN
																		WRITE(6,*) "************************************************"
																		WRITE(6,*) "Reducing diffusion: ",phys%diff_n*switch%diffred
																		WRITE(6,*) "************************************************"
															END IF
															elMat%Lf = elMat%Lf*switch%diffred
												   elMat%Qf = elMat%Qf*switch%diffred
												   elMat%P  = elMat%P*switch%diffred
												   phys%diff_n = phys%diff_n*switch%diffred
															phys%diff_u = phys%diff_u*switch%diffred
#ifdef TEMPERATURE															
               phys%diff_e = phys%diff_e*switch%diffred						
               phys%diff_ee = phys%diff_ee*switch%diffred
#endif
															u0 = sol%u
															time%it = 0
															! compute dt
												ELSE
											 ! Time advancing simulation						   
															IF (MPIvar%glob_id.eq.0) THEN
																		WRITE(6,*) "**********************"
																		WRITE(6,*) "Time scheme converged!"
																		WRITE(6,*) "**********************"
															END IF
															EXIT ! Here I exit the time advancing scheme if I reach convergence
												END IF
									ELSEIF (errTime.gt.numer%div) THEN
												WRITE(6,*) 'Problem in the time advancing scheme'
												STOP
									ELSE
									   IF (time%tis.eq.2) THEN
              um2 = u0
            END IF
												u0 = sol%u
												! compute dt
!												CALL compute_dt(errTime)
									END IF			  
						END IF
 
 END DO ! ************ END OF THE TIME LOOP *********************

		! Save solution
		CALL setSolName(save_name,mesh_name,time%it,.true.,.true.) 
		CALL HDF5_save_solution(save_name)
		 
 call cpu_time(time_finish)
 call system_clock(clock_end,clock_rate)
 print '("Elapsed cpu-time = ",f10.3," seconds.")',time_finish-time_start
 print '("Elapsed run-time = ",f10.3," seconds.")',(clock_end-clock_start)/real(clock_rate)
 
 
 IF (MPIvar%glob_id.eq.0) THEN
	if (utils%timing) then
 cputtot = cputpre+cputcon+cputpdf+cputbcd+cputmap+cputass+cputglb+cputsol
 runttot = runtpre+runtcon+runtpdf+runtbcd+runtmap+runtass+runtglb+runtsol
#ifdef PARALL
 cputtot = cputtot + MPIvar%ctime1
 runttot = runttot + MPIvar%rtime1
#endif
	WRITE(6,*) " "
	WRITE(6,*) " "
	WRITE(6,*) " "
	WRITE(6, '(" *", 90("*"), "**")') 
 WRITE(6, '(" *", 36X, "CODE TIMING ( Nthreads = ",i2,")", 26X, " *")' ) Nthreads
 WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")' )
 WRITE(6, '(" *", 2X,  "Precal. matr     : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') cputpre,cputpre/cputtot*100,runtpre,runtpre/runttot*100,cputpre/runtpre/Nthreads
 WRITE(6, '(" *", 2X,  "Convection       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputcon,cputcon/cputtot*100,runtcon,runtcon/runttot*100,cputcon/runtcon/Nthreads
#ifdef TEMPERATURE
 WRITE(6, '(" *", 2X,  "Parall. diff.    : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputpdf,cputpdf/cputtot*100,runtpdf,runtpdf/runttot*100,cputpdf/runtpdf/Nthreads
#endif
 WRITE(6, '(" *", 2X,  "Mapping          : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputmap,cputmap/cputtot*100,runtmap,runtmap/runttot*100,cputmap/runtmap/Nthreads
 WRITE(6, '(" *", 2X,  "Boundary cond.   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputbcd,cputbcd/cputtot*100,runtbcd,runtbcd/runttot*100,cputbcd/runtbcd/Nthreads 
 WRITE(6, '(" *", 2X,  "Assembly         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputass,cputass/cputtot*100,runtass,runtass/runttot*100,cputass/runtass/Nthreads
 WRITE(6, '(" *", 2X,  "Solve glob. syst.: ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputglb,cputglb/cputtot*100,runtglb,runtglb/runttot*100,cputglb/runtglb/Nthreads 
 WRITE(6, '(" *", 2X,  "Element solution : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   cputsol,cputsol/cputtot*100,runtsol,runtsol/runttot*100,cputsol/runtsol/Nthreads 
#ifdef PARALL
 WRITE(6, '(" *", 2X,  "Communications   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') MPIvar%ctime1,MPIvar%ctime1/cputtot*100,MPIvar%rtime1,MPIvar%rtime1/runttot*100,MPIvar%ctime1/MPIvar%rtime1/Nthreads
#endif
  WRITE(6, '(" *", 2X, "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")')   cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads 
	WRITE(6, '(" *", 90("*"), "**")')
	WRITE(6,*) " "
	WRITE(6,*) " "
	WRITE(6,*) " "
	end if 		


	if (lssolver%timing) then
 cputtot = matK%ctime1+matK%ctime2+matK%ctime3+matK%ctime4+matK%ctime5+matK%ctime6
 runttot = matK%rtime1+matK%rtime2+matK%rtime3+matK%rtime4+matK%rtime5+matK%rtime6
	WRITE(6, '(" *", 90("*"), "**")')
if (lssolver%sollib.eq.1) then	 
 WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PASTIX ( Nthreads = ",i2,")", 14X, " *")' ) Nthreads
else if (lssolver%sollib.eq.2) then
 WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PSBLAS ( Nthreads = ",i2,")", 14X, " *")' ) Nthreads
endif 
 WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")' )
if (lssolver%sollib.eq.1) then
 WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') matK%ctime1,matK%ctime1/cputtot*100,matK%rtime1,matK%rtime1/runttot*100,matK%ctime1/matK%rtime1/Nthreads
 WRITE(6, '(" *", 2X,  "Check mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime2,matK%ctime2/cputtot*100,matK%rtime2,matK%rtime2/runttot*100,matK%ctime2/matK%rtime2/Nthreads
 WRITE(6, '(" *", 2X,  "Anal. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime3,matK%ctime3/cputtot*100,matK%rtime3,matK%rtime3/runttot*100,matK%ctime3/matK%rtime3/Nthreads
 WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime4,matK%ctime4/cputtot*100,matK%rtime4,matK%rtime4/runttot*100,matK%ctime4/matK%rtime4/Nthreads
 WRITE(6, '(" *", 2X,  "LU decomp.       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime5,matK%ctime5/cputtot*100,matK%rtime5,matK%rtime5/runttot*100,matK%ctime5/matK%rtime5/Nthreads
 WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime6,matK%ctime6/cputtot*100,matK%rtime6,matK%rtime6/runttot*100,matK%ctime6/matK%rtime6/Nthreads
  WRITE(6, '(" *", 2X,  "Total time      : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",7X,F5.1 , 10X, " *")') cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
elseif (lssolver%sollib.eq.2) then
 WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') matK%ctime1,matK%ctime1/cputtot*100,matK%rtime1,matK%rtime1/runttot*100,matK%ctime1/matK%rtime1/Nthreads
 WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime2,matK%ctime2/cputtot*100,matK%rtime2,matK%rtime2/runttot*100,matK%ctime2/matK%rtime2/Nthreads
 WRITE(6, '(" *", 2X,  "Build prec       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime3,matK%ctime3/cputtot*100,matK%rtime3,matK%rtime3/runttot*100,matK%ctime3/matK%rtime3/Nthreads
 WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime4,matK%ctime4/cputtot*100,matK%rtime4,matK%rtime4/runttot*100,matK%ctime4/matK%rtime4/Nthreads
 WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') matK%ctime5,matK%ctime5/cputtot*100,matK%rtime5,matK%rtime5/runttot*100,matK%ctime5/matK%rtime5/Nthreads
 WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
endif 
	WRITE(6, '(" *", 90("*"), "**")')
	WRITE(6,*) " "
	WRITE(6,*) " "
	WRITE(6,*) " "
	end if 	
 END IF
	 
 IF (switch%testcase<5) THEN
    ALLOCATE(L2err(phys%neq))
    CALL computeL2ErrorAnalyticSol(L2err)
    WRITE(6,*) " "
    DO i=1,phys%Neq
    WRITE(6,'(A,I1,A,ES16.5)') "L2 error in U(",i,") = ", L2err(i)
    END DO
    DEALLOCATE(L2err)
 END IF
 
 DEALLOCATE(uiter,u0)
 IF (time%tis.eq.2) THEN
    DEALLOCATE(um2)
 END IF
 
 IF (lssolver%sollib .eq. 1) THEN
#ifdef WITH_PASTIX 
    CALL terminate_mat_PASTIX()
    
    ! MPI finalization
    call MPI_finalize(IERR)
#endif
 ELSEIF (lssolver%sollib .eq. 2) THEN
#ifdef WITH_PSBLAS 
    CALL terminate_PSBLAS()
#endif
 ENDIF
 

 
 CONTAINS
 
    !************************************************
    ! Display results
    !************************************************ 
    SUBROUTINE displayResults()
    integer              :: ieq
    real*8,allocatable   :: uphy(:,:)
    real*8               :: Vmax(phys%npv),Vmin(phys%npv)
    
						 ALLOCATE(uphy(Mesh%Nnodesperelem*Mesh%Nelems,phys%npv))
						 
						 ! Compute physical variables
							CALL cons2phys(transpose(reshape(sol%u,(/ phys%Neq,Mesh%Nnodesperelem*Mesh%Nelems /))),uphy) 
		     DO ieq = 1, phys%npv
		        Vmax(ieq) = Maxval(uphy(:,ieq))
		        Vmin(ieq) = Minval(uphy(:,ieq))
		     END DO				
		     		
#ifdef PARALL
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vmax,phys%npv,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,Vmin,phys%npv,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
#endif		 
       IF (MPIvar%glob_id.eq.0) THEN
		     WRITE(6, '(" * ", 60("-"), "*")')						     
		     WRITE(6, '(" * Time (adimensional) = ", E12.5, 27X, " *")') time%t		     
		     WRITE(6, '(" * Dt (seconds)        = ", E12.5, 27X, " *")') time%dt*1.3736e-07
		     WRITE(6, '(" * ", 45("^"), 14X, " *")')
		     WRITE(6, '(" * ", 10("_"), "      Minimum   ", 4X, "     Maximum   ", 14X, " *")')
		     DO ieq = 1, phys%npv
		        WRITE(6, '(" * ", A7, " -->", ES16.8, 3X, ES16.8, 13X, " *")') &
		        & Trim(phys%phyVarNam(ieq)), Vmin(ieq), Vmax(ieq)
		     END DO
		     WRITE(6, '(" * ", 60("-"), "*")')
       WRITE(6, '(" * Time residual  = ", 1X, 2(E16.8, 2X), 13X, " *")') sol%tres(it)
       WRITE(6,*) '  '
       WRITE(6,*) '  '
       WRITE(6,*) '  '
			    END IF
#ifdef PARALL
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif			    
						  DEALLOCATE(uphy)         
    
    END SUBROUTINE displayResults
    
    !************************************************
    ! Compute the residual
    !************************************************
    FUNCTION computeResidual(u,uref,coeff) RESULT(res)
    real*8   :: u(nu),uref(nu)
    integer  :: nglo,ierr
    real     :: res,sum2,coeff
    
    sum2 = sum( (u-uref)**2)
    nglo = nu

#ifdef PARALL
    call mpi_allreduce(MPI_IN_PLACE,sum2,1,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)   
    call mpi_allreduce(nu,nglo,1,mpi_integer,mpi_sum,MPI_COMM_WORLD,ierr)    
#endif
    res = sqrt(sum2)/sqrt(dble(nglo))/coeff
    END FUNCTION computeResidual
    
    !************************************************
    ! Set the name of the solution
    !************************************************    
    SUBROUTINE setSolName(save_name,mesh_name,it,convNR,convT)
     character(LEN = 1024), intent(INOUT):: save_name
     character(LEN = 1024), intent(IN)   :: mesh_name
     integer,intent(IN)                  :: it
     logical,intent(IN)                  :: convNR,convT
     character(LEN = 20)                 :: Num
     integer                             :: l,i

     ! At the beginning, the save name is the mesh name..
     save_name = TRIM(ADJUSTL(mesh_name))
     
     ! Eliminate path info
     l = LEN(save_name)
     i = INDEX(save_name, '/' ,.true.)  
     save_name = save_name(i+1:l)

#ifdef TOR3D
     ! Add the number of toroidal elements
     WRITE(Num, "(i10)") numer%ntor
     Num = TRIM(ADJUSTL(Num))
     save_name = TRIM(ADJUSTL(save_name))//"_Ntor"//Num    
     
     ! Add the poloidal interpolation in the toroidal direction
     WRITE(Num, "(i10)") numer%ptor
     Num = TRIM(ADJUSTL(Num))
     save_name = TRIM(ADJUSTL(save_name))//"Ptor"//Num 
              
#endif  
     
     ! Diffusion
     WRITE(Num, "(E10.3)") phys%diff_n
     save_name = TRIM(ADJUSTL(save_name))//"_DPe"//TRIM(ADJUSTL(Num))
#ifdef TEMPERATURE
     WRITE(Num, "(E10.3)") phys%diff_pari
     save_name = TRIM(ADJUSTL(save_name))//"_DPai"//TRIM(ADJUSTL(Num))
     
     WRITE(Num, "(E10.3)") phys%diff_pare
     save_name = TRIM(ADJUSTL(save_name))//"_DPae"//TRIM(ADJUSTL(Num))     
#endif
     ! Complete the save name: if not converged the NR, I put NR + the iteration number
     IF (.not.convNR) THEN
     			WRITE(Num, "(i10)") it
 				   Num = TRIM(ADJUSTL(Num))
        k = INDEX(Num, " ") -1
        save_name = TRIM(ADJUSTL(save_name))//"_NR"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
     END IF

     ! Complete the save name: if not converged the time scheme, I put the iteration number
     IF (.not.convT) THEN
     			WRITE(Num, "(i10)") it/utils%freqsave
 				   Num = TRIM(ADJUSTL(Num))
        k = INDEX(Num, " ") -1
        save_name = TRIM(ADJUSTL(save_name))//"_"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
     END IF

     IF (switch%decoup) THEN
        save_name = trim(adjustl(save_name))//'_UNCP'
     ENDIF
     
     ! Add "Sol_"
#ifdef TOR3D     
     save_name = 'Sol3D_'//save_name
#else
     save_name = 'Sol2D_'//save_name
#endif 
    END SUBROUTINE setSolName 
    
    
    
    
    
    
    SUBROUTINE compute_dt(errtime)
    real*8, intent(in) :: errtime
    
    if ( (errtime*time%dt)<1e-3) then
       write(6,*) "******** Changing time step ***********"
       time%dt = time%dt*2.
    endif
    END SUBROUTINE compute_dt
    
END
