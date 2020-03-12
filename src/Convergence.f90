PROGRAM Convergence
 USE in_out
 USE reference_element
 USE preprocess
 USE MPI_OMP
 USE printutils
 USE debug
 USE initialization
 USE Postprocess
 USE solve_pastix
 IMPLICIT NONE
 
 integer             :: it, ir, nts,nu,nb_args, n_points, n_polyns
 integer             :: ieq
 integer             :: ipol,ipts
 real*8              :: dt, errNR, errTime
 real*8, allocatable :: uiter(:),u0(:)
 character(LEN = 1024) :: cname,mesh_name,n_points_c,n_polyns_c
 integer, allocatable:: Nelv(:,:),Nunk(:,:)
 real*8,allocatable  :: L2err(:,:,:),L2err_pt(:),elSize(:,:),slope(:,:,:)
 
  ! Check the number of input arguments
  nb_args = iargc()
  IF (nb_args .lt. 1) THEN
        PRINT *, " Error: please provide the case name, the number of points and the polynomial degrees"
      stop
  ELSEIF (nb_args .lt. 2) THEN
        PRINT *, " Error: please provide the number of points and the polynomial degrees"
      stop  
  ELSEIF (nb_args .lt. 3) THEN
      PRINT *, " Error: please provide the polynomial degrees"
      stop
  ELSEIF (nb_args .gt. 3) THEN
      PRINT *, " Error: too many input"
      stop
  END IF
  
  CALL getarg(1, cname)
  CALL getarg(2, n_points_c)
  CALL getarg(3, n_polyns_c)
  READ (n_points_c,*) n_points
  READ (n_polyns_c,*) n_polyns
  
	! Initialize MPI
	CALL init_MPI_OMP()

	! Read input file param.txt
	CALL read_input()
		
	! Initialization of the simulation parameters
	CALL init_sim(nts,dt)
			
	! Allocate result vectors
	ALLOCATE(L2err_pt(phys%Neq))
	ALLOCATE(L2err(n_points,n_polyns,phys%Neq))
	ALLOCATE(elSize(n_points,n_polyns))
	ALLOCATE(Nelv(n_points,n_polyns))
	ALLOCATE(Nunk(n_points,n_polyns))
	ALLOCATE(Slope(n_points,n_polyns,phys%Neq))
	
	! Initialization of the results
	L2err    = 0.
	elSize   = 0.
	Nelv     = 0
	Nunk     = 0
	slope    = 0.
	
	WRITE(6,*) "Computing a convergence on the mesh ", trim(cname), "with ", trim(n_points_c),  " points and ", trim(n_polyns_c), " polynomial degrees"
	
	
 DO ipol = 1, n_polyns   
    DO ipts = 1,n_points
    
       WRITE(mesh_name,'(A,A,I0,A,I0,A)') Trim(cname),'_',ipts,'_P',ipol
       WRITE(6,*) "Computing ", trim(mesh_name)

							! Load the mesh file
							CALL load_mesh(Trim(mesh_name))
		
							! Pastix: set the start to true
							matPASTIX%start = .true.
									
							! Create the reference element based on the mesh type
							CALL create_reference_element() 
		
							! Mesh preprocess: create the mesh related structures 
							! used in the HDG scheme
							CALL mesh_preprocess()
		
							! Allocation and initialization of the elemental matrices
							CALL init_elmat()
		
							! Compute the elemental matrices that depend only on the 
							! mesh, i.e. fixed in time
							CALL hdg_precalculatedMatrices()

							! Initialize the solution
							CALL init_sol()
				
							! Allocate and initialize uiter and u0
							nu = size(sol%u)
							ALLOCATE(uiter(1:nu))
							ALLOCATE(u0(1:nu))
							u0 = sol%u
							uiter = 0.
							!*******************************************************
							!                  TIME LOOP
							!*******************************************************
							DO it = 1,nts
		
										! Compute time step
								
										! Actualization of time 
										time%t  = time%t + dt
										time%it = time%it + 1 ! is it of any use??
										time%dt = dt 
										WRITE(6, '(" *", 60("*"), "**")') 
										WRITE(6, '(" *", 20X,    "Time iteration   = ", I5, 16X, " *")') time%it
										WRITE(6, '(" *", 60("*"), "**")') 
															!*******************************************************
															!             Newton-Raphson iterations
															!*******************************************************
															uiter = u0
															DO ir = 1,numer%nrp
									
																		! Compute the elemental matrices that depends on the N-R iteration
																		CALL hdg_ConvectionMatrices()
																		
#ifdef TEMPERATURE
                  ! Compute the elemental matrices relative to temperature
                  CALL HDG_paralldiffmatrices()
#endif																		
															
																		! Set boundary conditions
																		CALL hdg_BC()
																														   									
																		! Compute elemental mapping
																		CAlL hdg_Mapping(u0)

																		! Assembly the global matrix
																		CALL hdg_Assembly()
												      
																		! Solve linear system, compute face and element solution
																		CALL compute_solution()
															
																		! Check convergence of Newton-Raphson
																		errNR = norm2(sol%u-uiter)
															
																		WRITE(6,*) "NR iteration: ", ir, ' - error: ', errNR
																		IF (errNR.lt.numer%tNR) THEN
																					EXIT
																		ELSEIF (errNR.gt.numer%div) THEN
																					WRITE(6,*) 'Problem in the N-R procedure'
																					STOP
																		ELSE
																					uiter = sol%u
																		END IF
																		
															END DO ! End of NR loop
															WRITE(6,*) " "
									      WRITE(6,*) " "
									      
												! Check convergence in time advancing and update
												errTime = norm2(sol%u-u0)/dt
															
												IF (errTime.lt.numer%tTM) THEN
															WRITE(6,*) "**********************"
															WRITE(6,*) "Time scheme converged!"
															WRITE(6,*) "**********************"
															EXIT
												ELSEIF (errTime.gt.numer%div) THEN
															WRITE(6,*) 'Problem in the time advancing scheme'
															STOP
												ELSE
															u0 = sol%u
												END IF			
												WRITE(6,*) " "  
							END DO ! End time loop
							DEALLOCATE(uiter,u0)
							
							! Store results
							CALL computeL2ErrorAnalyticSol(L2err_pt)
				   L2err(ipts,ipol,:) = L2err_pt
				   elSize(ipts,ipol) = maxval(Mesh%elemSize)
				   Nelv(ipts,ipol) = Mesh%Nelems
				   Nunk(ipts,ipol) = MatK%n
				   IF (ipts.gt.1) THEN
				      slope(ipts,ipol,:) = (log10(L2err(ipts,ipol,:))-log10(L2err(ipts-1,ipol,:)))/ &
				                           (log10(elSize(ipts,ipol))-log10(elSize(ipts-1,ipol)))
!				      				      (log10(sqrt(dble(Nunk(ipts-1,ipol))))-log10(sqrt(dble(Nunk(ipts,ipol)))))
				   END IF
				   DO ieq = 1,phys%Neq
				      WRITE(6,* ) "Error for ", Trim(phys%conVarNam(ieq)), ": ",L2err_pt(ieq)
				   END DO
				   
				   ! Free allocated variables
				   CALL free_all()
				   time%it = 0
				   
				END DO ! Loop in number of convergence points
 END DO ! Loop in number of polynomial degrees
 
 !*******************************************
 ! Display results
 !*******************************************
 ! I reinitialize to allocate some deallocates 
 ! variables
 CALL init_sim(nts,dt)
 DO ieq = 1,phys%Neq
 
    WRITE(6, '(" * ", 55("-"), "*")')
    WRITE(6, '(" * Convergence for ", A7, 31X, " *")') Trim(phys%conVarNam(ieq))
    DO ipol = 1,n_polyns
       WRITE(6, '(" * ", 55("-"), "*")')
       WRITE(6, '(" * Interpolation: p", I0, 38X, "*")'),ipol
       WRITE(6, '(" * ",5X, " h ", 7X, " Nel ", 5X, " Nunk ", 2X, " L2err", 6X, "Slope", 5X, "*")')       
          DO ipts = 1,n_points
          WRITE(6, '(" * ",2X,ES10.2E3, 2X, I5, 2X, I8, 2X, ES10.2E3, 2X, ES10.2E3, 2X, "*")'), &
                    elSize(ipts,ipol),Nelv(ipts,ipol),Nunk(ipts,ipol),L2err(ipts,ipol,ieq),Slope(ipts,ipol,ieq)
          END DO    
    END DO
    WRITE(6, '(" * ", 55("-"), "*")')
    WRITE(6, *) " "
 END DO
 
 
 DEALLOCATE(L2err,Nunk,elSize,Slope,Nelv)

END PROGRAM
