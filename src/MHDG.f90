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

  integer             :: Neq, Np, Nel, Ndim, Nfp, Nf, ntorloc
  integer             :: Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor, Nthreads
  integer             :: it, ir, ifa, nts, nu, nut, nb_args, IERR, k, i, j, g, is
  integer, allocatable :: faceNodes1d(:)
  logical, allocatable :: mkelms(:)
  real*8              :: dt, errNR, errlstime
  real*8, allocatable :: uiter(:), L2err(:)
  character(LEN=1024) :: mesh_name, namemat, save_name
  real*8              :: cputtot, runttot
  integer             ::  OMP_GET_MAX_THREADS
  integer             :: cks, clock_rate, cke, clock_start, clock_end
  real*8              :: time_start, time_finish, tps, tpe
  INTEGER :: code, pr, aux
  INTEGER, PARAMETER :: etiquette = 1000
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut
  INTEGER             :: switch_save
  write (6, *) "STARTING"

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
  call system_clock(clock_start, clock_rate)

  ! Initialize MPI
  CALL init_MPI_OMP()

  ! Read input file param.txt
  CALL read_input()

  ! Set parallelization division in toroidal and poloidal plane
#ifdef TOR3D
#ifdef PARALL
  call set_divisions()
#endif
#endif

  ! Load the mesh file
  CALL load_mesh(mesh_name)


  ! Linear solver: set the start to true
  matK%start = .true.
  if (lssolver%timing) then
    call init_solve_timing
  endif
  ! Initialize marked elements for thresholds
  allocate (mkelms(Mesh%Nelems))
  mkelms = .false.


  ! Create the reference element based on the mesh type
  CALL create_reference_element(refElPol, 2)
#ifdef TOR3D
  !*************************************************
  !                REFERENCE ELEMENT 3D
  !*************************************************
  ! Create the reference element for the toroidal interpolation
  CALL create_reference_element(refElTor, 1, numer%ptor)
  CALL create_toroidal_structures(refElTor, refElPol)
#endif


  ! Mesh preprocess: create the mesh related structures
  ! used in the HDG scheme
  CALL mesh_preprocess()


#ifdef TOR3D
  ! Define toroidal discretization
  CALL define_toroidal_discretization
#endif



  ! Initialization of the simulation parameters !TODO check if I need to pass dt to init_sim
  CALL init_sim(nts, dt)

  ! Initialize magnetic field (the Mesh is needed)
  CALL initialize_magnetic_field()

  ! Load magnetic field
  CALL load_magnetic_field()

#ifdef TOR3D
  Ndim = 3                                               ! Number of dimensions
#ifdef PARALL
  ntorloc = numer%ntor/MPIvar%ntor
#else
  ntorloc = numer%ntor
#endif
  Nel = Mesh%Nelems*ntorloc                             ! Number of 3D elements
  Np = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
  Nfl = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
  Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
  Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
#ifdef PARALL
  if (MPIvar%ntor .gt. 1) then
    nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + &
      &phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
  else
    nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
  endif
#else
  nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
#endif
  Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
  nu = phys%Neq*Nel*Np
#else
  Ndim = Mesh%ndim
  Nel = Mesh%Nelems
  Np = refElPol%Nnodes2D
  Nf = refElPol%Nfaces
  Nfp = refElPol%Nfacenodes*Nf
  nut = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
  nu = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif

#ifdef PARALL
  CALL init_com()
#endif

  ! Allocation and initialization of the elemental matrices
  CALL init_elmat()

  ! Compute first equation (definition of the gradient)
  CALL HDG_precalculatedfirstequation()


  ! Initialize shock capturing
  if (switch%shockcp.gt.0) then
    CALL initializeShockCapturing()
  endif


  ! Initialize the solution
  IF (nb_args .eq. 2) THEN
    ! restart simulation: load solution from file (the name is given in argument)
    CALL HDF5_load_solution(save_name)
  ELSE
    CALL init_sol()
  END IF

  IF (switch%pertini .eq. 1) THEN
    call add_perturbation()
    write(6,*) "Adding perturbation to the initial solution"
  ELSE IF (switch%pertini .eq. 2) THEN
    call add_blob()
    write(6,*) "Adding density blob to initial solution"
  ENDIF

  ! Save solution
  CALL setSolName(save_name, mesh_name, 0, .true., .false.)
  CALL HDF5_save_solution(save_name)








  call mpi_barrier(mpi_comm_world,ierr)



  ! Allocate and initialize uiter and u0
  ALLOCATE (uiter(nu))
  ALLOCATE (sol%u0(nu, time%tis))
  sol%u0 = 0.
  sol%u0(:, 1) = sol%u
  uiter = 0.

  switch_save = 0
  !*******************************************************
  !                  TIME LOOP
  !*******************************************************
  DO it = 1, nts ! ************ TIME LOOP *********************

    ! Compute time step

    ! Actualization of time
    time%t = time%t + time%dt
    time%it = time%it + 1
    time%ik = time%ik + 1
    sol%Nt = sol%Nt + 1

    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, '(" *", 60("*"), "**")')
      WRITE (6, '(" *", 20X,    "Time iteration   = ", I5, 16X, " *")') time%it
      WRITE (6, '(" *", 60("*"), "**")')
    END IF

    !*******************************************************
    !             Newton-Raphson iterations
    !*******************************************************
    uiter = sol%u0(:, 1)
    DO ir = 1, numer%nrp ! ************ NEWTON-RAPHSON LOOP *********************

      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) "***** NR iteration: ", ir, "*****"
      ENDIF

      ! Compute Jacobian
      CALL HDG_computeJacobian()

      ! Set boundary conditions
      CALL hdg_BC()

      ! Compute elemental mapping
      CAlL hdg_Mapping()

      ! Assembly the global matrix
      CALL hdg_Assembly()
      !IF (switch_save.EQ.0) THEN
      !  WRITE (6, *) "Save matrix"
      !  call HDF5_save_CSR_matrix('Mat')
      !  call HDF5_save_CSR_vector('rhs')
      !  switch_save = 1
      !call displayMatrixInt(Mesh%F)
      !call displayMatrixInt(Mesh%extfaces)
      !call displayVectorInt(Mesh%periodic_faces)
      !stop
      !call print_matrices_hdf5
      !stop
      !ENDIF
      ! Solve linear system
      CALL solve_global_system()

      ! Compute element-by-element solution
      CALL compute_element_solution()

      ! Check for NaN (doesn't work with optimization flags)
      DO i = 1, nu
        IF (sol%u(i) /= sol%u(i)) THEN
          WRITE (6, *) "NaN detected"
          STOP
        END IF
      END DO

      ! Apply threshold
      CALL HDG_applyThreshold(mkelms)

      ! Apply filtering
      CALL HDG_FilterSolution()

      ! Save solution
      IF (switch%saveNR) THEN
        CALL setSolName(save_name, mesh_name, ir, .false., .true.)
        CALL HDF5_save_solution(save_name)
      END IF

      ! Check convergence of Newton-Raphson
      errNR = computeResidual(sol%u, uiter, 1.)
      errNR = errNR/numer%dumpnr

      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) "Error: ", errNR
      END IF
      IF (errNR .lt. numer%tNR) THEN
        EXIT
      ELSEIF (errNR .gt. numer%div) THEN
        WRITE (6, *) 'Problem in the N-R procedure'
        STOP
      ELSE
        uiter = sol%u
      END IF
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) "*********************************"
        WRITE (6, *) " "
        WRITE (6, *) " "
      ENDIF
    END DO ! ************ END OF NEWTON-RAPHSON LOOP *********************

    !  ! Apply threshold
    !  CALL HDG_applyThreshold()

    ! Check convergence in time advancing and update
    errlstime = computeResidual(sol%u, sol%u0(:, 1), time%dt)
    sol%tres(sol%Nt) = errlstime
    sol%time(sol%Nt) = time%t


    ! Display results
    IF (mod(time%it, utils%freqdisp) .eq. 0) THEN
      CALL displayResults()
    END IF

    ! Check for NaN (doesn't work with optimization flags)
    DO i = 1, nu
      IF (sol%u(i) /= sol%u(i)) THEN
        WRITE (6, *) "NaN detected"
        STOP
      END IF
    END DO

    IF (.not. switch%steady) THEN
      ! Save solution
      IF (.not. switch%steady) THEN
        IF (mod(time%it, utils%freqsave) .eq. 0) THEN
          CALL setSolName(save_name, mesh_name, time%it, .true., .false.)
          CALL HDF5_save_solution(save_name)
        END IF
      END IF

      !****************************************
      ! Check steady state and update or exit
      !****************************************
      IF (errlstime .lt. numer%tTM) THEN
        IF (switch%psdtime) THEN
          ! Pseudo-time simulation

          ! Save solution
          CALL setSolName(save_name, mesh_name, it, .true., .true.)
          CALL HDF5_save_solution(save_name)

          ! Update the diffusion, the elemental matrices and the solution
          IF (MPIvar%glob_id .eq. 0) THEN
            WRITE (6, *) "************************************************"
            WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion
            WRITE (6, *) "************************************************"
          END IF
          phys%diff_n = phys%diff_n*switch%diffred
          phys%diff_u = phys%diff_u*switch%diffred
#ifdef TEMPERATURE
          phys%diff_e = phys%diff_e*switch%diffred
          phys%diff_ee = phys%diff_ee*switch%diffred
#endif
#ifdef VORTICITY
          phys%diff_vort = phys%diff_vort*switch%diffred
          phys%diff_pot = phys%diff_pot*switch%diffred
#endif
#ifdef NEUTRAL
          phys%diff_nn = phys%diff_nn*switch%diffred
#endif

          !**********************************
          !           UPDATE SOLUTION
          !**********************************
          ! Update u0
          IF (time%tis .gt. 1 .and. it .lt. time%tis) THEN
            DO is = it, 1, -1
              sol%u0(:, is + 1) = sol%u0(:, is)
            END DO
          ELSEIF (time%tis > 1) THEN
            DO is = time%tis, 2, -1
              sol%u0(:, is) = sol%u0(:, is - 1)
            END DO
          END IF
          sol%u0(:, 1) = sol%u

          time%it = 0
          ! compute dt
        ELSE
          ! Time advancing simulation
          IF (MPIvar%glob_id .eq. 0) THEN
            WRITE (6, *) "**********************"
            WRITE (6, *) "Time scheme converged!"
            WRITE (6, *) "**********************"
          END IF
          EXIT ! Here I exit the time advancing scheme if I reach convergence
        END IF
      ELSEIF (errlstime .gt. numer%div) THEN
        WRITE (6, *) 'Problem in the time advancing scheme'
        STOP
      ELSE
        !**********************************
        !           UPDATE SOLUTION
        !**********************************
        ! Update u0
        IF (time%tis .gt. 1 .and. it .lt. time%tis) THEN
          DO is = it, 1, -1
            sol%u0(:, is + 1) = sol%u0(:, is)
          END DO
        ELSEIF (time%tis > 1) THEN
          DO is = time%tis, 2, -1
            sol%u0(:, is) = sol%u0(:, is - 1)
          END DO
        END IF
        sol%u0(:, 1) = sol%u

        ! compute dt
        ! CALL compute_dt(errlstime)
      END IF
    END IF

  END DO ! ************ END OF THE TIME LOOP *********************

  ! Save solution
  CALL setSolName(save_name, mesh_name, time%it, .true., .true.)
  CALL HDF5_save_solution(save_name)

  call cpu_time(time_finish)
  call system_clock(clock_end, clock_rate)
  print '("Elapsed cpu-time = ",f10.3," seconds.")', time_finish - time_start
  print '("Elapsed run-time = ",f10.3," seconds.")', (clock_end - clock_start)/real(clock_rate)

  ! Code timing
  IF (MPIvar%glob_id .eq. 0) THEN
    if (utils%timing) then
      cputtot = 1e-8
      runttot = 1e-8
      cputtot = timing%cputpre + timing%cputjac + timing%cputbcd + timing%cputmap + timing%cputass + timing%cputglb + timing%cputsol
      runttot = timing%runtpre + timing%runtjac + timing%runtbcd + timing%runtmap + timing%runtass + timing%runtglb + timing%runtsol
#ifdef PARALL
      cputtot = cputtot + timing%cputcom
      runttot = runttot + timing%runtcom
#endif
      WRITE(6, *) " "
      WRITE(6, *) " "
      WRITE(6, *) " "
      WRITE(6, '(" *", 90("*"), "**")')
      WRITE(6, '(" *", 36X, "CODE TIMING ( Nthreads = ",i2,")", 26X, " *")') Nthreads
      WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
      WRITE(6, '(" *", 2X,  "Precal. matr     : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
        &timing%cputpre,timing%cputpre/cputtot*100,timing%runtpre,timing%runtpre/runttot*100,timing%cputpre/timing%runtpre/Nthreads
      WRITE(6, '(" *", 2X,  "Jacobian         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputjac,timing%cputjac/cputtot*100,timing%runtjac,timing%runtjac/runttot*100,timing%cputjac/timing%runtjac/Nthreads
      WRITE(6, '(" *", 2X,  "Mapping          : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputmap,timing%cputmap/cputtot*100,timing%runtmap,timing%runtmap/runttot*100,timing%cputmap/timing%runtmap/Nthreads
      WRITE(6, '(" *", 2X,  "Boundary cond.   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputbcd,timing%cputbcd/cputtot*100,timing%runtbcd,timing%runtbcd/runttot*100,timing%cputbcd/timing%runtbcd/Nthreads
      WRITE(6, '(" *", 2X,  "Assembly         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputass,timing%cputass/cputtot*100,timing%runtass,timing%runtass/runttot*100,timing%cputass/timing%runtass/Nthreads
      WRITE(6, '(" *", 2X,  "Solve glob. syst.: ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputglb,timing%cputglb/cputtot*100,timing%runtglb,timing%runtglb/runttot*100,timing%cputglb/timing%runtglb/Nthreads
      WRITE(6, '(" *", 2X,  "Element solution : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputsol,timing%cputsol/cputtot*100,timing%runtsol,timing%runtsol/runttot*100,timing%cputsol/timing%runtsol/Nthreads
#ifdef PARALL
      WRITE(6, '(" *", 2X,  "Communications   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
        &timing%cputcom,timing%cputcom/cputtot*100,timing%runtcom,timing%runtcom/runttot*100,timing%cputcom/timing%runtcom/Nthreads
#endif
      WRITE(6, '(" *", 2X, "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")')   &
        cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
      WRITE(6, '(" *", 90("*"), "**")')
      WRITE(6, *) " "
      WRITE(6, *) " "
      WRITE(6, *) " "
    end if

    if (lssolver%timing) then
      cputtot = timing%clstime1 + timing%clstime2 + timing%clstime3 + timing%clstime4 + timing%clstime5 + timing%clstime6
      runttot = timing%rlstime1 + timing%rlstime2 + timing%rlstime3 + timing%rlstime4 + timing%rlstime5 + timing%rlstime6
      WRITE(6, '(" *", 90("*"), "**")')
      if (lssolver%sollib .eq. 1) then
        WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PASTIX ( Nthreads = ",i2,")", 14X, " *")') Nthreads
      else if (lssolver%sollib .eq. 2) then
        WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PSBLAS ( Nthreads = ",i2,")", 14X, " *")') Nthreads
      endif
      WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
      if (lssolver%sollib .eq. 1) then
        WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
          &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
        WRITE(6, '(" *", 2X,  "Check mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
        WRITE(6, '(" *", 2X,  "Anal. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
        WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
        WRITE(6, '(" *", 2X,  "LU decomp.       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
        WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime6,timing%clstime6/cputtot*100,timing%rlstime6,timing%rlstime6/runttot*100,timing%clstime6/timing%rlstime6/Nthreads
        WRITE(6, '(" *", 2X,  "Total time      : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",7X,F5.1 , 10X, " *")') &
          &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
      elseif (lssolver%sollib .eq. 2) then
        WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
          &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
        WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
        WRITE(6, '(" *", 2X,  "Build prec       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
        WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
        WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
          &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
        WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
          &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
      endif
      WRITE (6, '(" *", 90("*"), "**")')
      WRITE (6, *) " "
      WRITE (6, *) " "
      WRITE (6, *) " "
    end if
  END IF

  IF (switch%testcase < 5) THEN
    ALLOCATE (L2err(phys%neq))
    CALL computeL2ErrorAnalyticSol(L2err)
    WRITE (6, *) " "
    DO i = 1, phys%Neq
      WRITE (6, '(A,I1,A,ES16.5)') "L2 error in U(", i, ") = ", L2err(i)
    END DO
    DEALLOCATE (L2err)
  END IF

  DEALLOCATE (uiter)
  deallocate (mkelms)
  DEALLOCATE (sol%u0)

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
    real*8, allocatable   :: uphy(:, :)
    real*8               :: Vmax(phys%npv), Vmin(phys%npv)

    ALLOCATE (uphy(nu/phys%Neq, phys%npv))

    ! Compute physical variables
    CALL cons2phys(transpose(reshape(sol%u, (/phys%Neq, nu/phys%Neq/))), uphy)
    DO ieq = 1, phys%npv
      Vmax(ieq) = Maxval(uphy(:, ieq))
      Vmin(ieq) = Minval(uphy(:, ieq))
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmax, phys%npv, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmin, phys%npv, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, '(" * ", 60("-"), "*")')
      WRITE (6, '(" * Time (adimensional) = ", E12.5, 27X, " *")') time%t
      WRITE (6, '(" * Dt (adimensional)        = ", E12.5, 27X, " *")') time%dt
      WRITE (6, '(" * ", 45("^"), 14X, " *")')
      WRITE (6, '(" * ", 10("_"), "      Minimum   ", 4X, "     Maximum   ", 14X, " *")')
      DO ieq = 1, phys%npv
        WRITE (6, '(" * ", A7, " -->", ES16.8, 3X, ES16.8, 13X, " *")') &
          & Trim(phys%phyVarNam(ieq)), Vmin(ieq), Vmax(ieq)
      END DO
      WRITE (6, '(" * ", 60("-"), "*")')
      WRITE (6, '(" * Time residual  = ", 1X, 2(E16.8, 2X), 13X, " *")') sol%tres(it)
      WRITE (6, *) '  '
      WRITE (6, *) '  '
      WRITE (6, *) '  '
    END IF
#ifdef PARALL
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    DEALLOCATE (uphy)

  END SUBROUTINE displayResults

  !************************************************
  ! Compute the residual
  !************************************************
  FUNCTION computeResidual(u, uref, coeff) RESULT(res)
    real*8   :: u(nu), uref(nu)
    integer  :: nglo, ierr
    real     :: res, sum2, coeff

    sum2 = sum((u - uref)**2)
    nglo = nu

#ifdef PARALL
    call mpi_allreduce(MPI_IN_PLACE, sum2, 1, mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(nu, nglo, 1, mpi_integer, mpi_sum, MPI_COMM_WORLD, ierr)
#endif

    res = sqrt(sum2)/sqrt(dble(nglo))/coeff
  END FUNCTION computeResidual

  !************************************************
  ! Set the name of the solution
  !************************************************
  SUBROUTINE setSolName(save_name, mesh_name, it, convNR, convT)
    character(LEN=1024), intent(INOUT):: save_name
    character(LEN=1024), intent(IN)   :: mesh_name
    integer, intent(IN)                  :: it
    logical, intent(IN)                  :: convNR, convT
    character(LEN=20)                 :: Num
    integer                             :: l, i

    ! At the beginning, the save name is the mesh name..
    save_name = TRIM(ADJUSTL(mesh_name))

    ! Eliminate path info
    l = LEN(save_name)
    i = INDEX(save_name, '/', .true.)
    save_name = save_name(i + 1:l)

#ifdef TOR3D
    ! Add the number of toroidal elements
    WRITE (Num, "(i10)") numer%ntor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"_Ntor"//Num

    ! Add the poloidal interpolation in the toroidal direction
    WRITE (Num, "(i10)") numer%ptor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"Ptor"//Num

#endif

    ! Diffusion
    WRITE (Num, "(E10.3)") phys%diff_n*simpar%refval_diffusion
    save_name = TRIM(ADJUSTL(save_name))//"_DPe"//TRIM(ADJUSTL(Num))
#ifdef TEMPERATURE
    WRITE (Num, "(E10.3)") phys%diff_pari
    save_name = TRIM(ADJUSTL(save_name))//"_DPai"//TRIM(ADJUSTL(Num))

    WRITE (Num, "(E10.3)") phys%diff_pare
    save_name = TRIM(ADJUSTL(save_name))//"_DPae"//TRIM(ADJUSTL(Num))
#endif
    ! Complete the save name: if not converged the NR, I put NR + the iteration number
    IF (.not. convNR) THEN
      WRITE (Num, "(i10)") it
      Num = TRIM(ADJUSTL(Num))
      k = INDEX(Num, " ") - 1
      save_name = TRIM(ADJUSTL(save_name))//"_NR"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
    END IF

    ! Complete the save name: if not converged the time scheme, I put the iteration number
    IF (.not. convT) THEN
      WRITE (Num, "(i10)") it/utils%freqsave
      Num = TRIM(ADJUSTL(Num))
      k = INDEX(Num, " ") - 1
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

  SUBROUTINE compute_dt(errlstime)
    real*8, intent(in) :: errlstime

    if ((errlstime*time%dt) < 1e-3) then
      write (6, *) "******** Changing time step ***********"
      time%dt = time%dt*2.
    endif
  END SUBROUTINE compute_dt


#ifdef TOR3D
  !**********************************************
  ! Definition of the toroidal discretization
  !**********************************************
  SUBROUTINE define_toroidal_discretization
    integer :: i, ntorloc, itor, itorg, nnodes_toroidal
    integer :: ind(numer%ptor + 1)
    real*8  :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

    ! Toroidal discretization
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1, numer%ntor
      tdiv(i + 1) = i*htor
    END DO
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
      ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif

    nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
    Mesh%Nnodes_toroidal = nnodes_toroidal
    ALLOCATE (Mesh%toroidal(nnodes_toroidal))

    DO itor = 1, ntorloc
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
      ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
      Mesh%toroidal(ind) = tel
    END DO
  END SUBROUTINE define_toroidal_discretization
#endif
  END
