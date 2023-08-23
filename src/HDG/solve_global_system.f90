!*****************************************
! project: MHDG
! file: solve_global_system.f90
! date: 23/07/2019
! Solve the global system
!*****************************************

SUBROUTINE solve_global_system
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE in_out
#ifdef WITH_PASTIX
  USE solve_pastix
#endif
#ifdef WITH_PSBLAS
  USE solve_psblas
#endif
  USE MPI_OMP
#ifdef PARALL
  USE Communications
#endif
  IMPLICIT NONE

  integer*4             :: i, j
  real, allocatable      :: rhspert(:)
  real                  :: pertamp, errsol
#ifdef PARALL
  integer*4             :: ierr, Neq, Nfp
  integer*4             :: ct, indg(refElPol%Nfacenodes*phys%Neq), indl(refElPol%Nfacenodes*phys%Neq)
  real*8, allocatable    :: aux_sol(:)
#ifdef TOR3D
  integer*4             :: Nfl, N2d, Np2d, itor, Nfaces, Nfdir, Nghostf, Nghoste, dd, ddl, ntorass
  integer*4             :: indgp(refElPol%Nnodes2D*phys%Neq),indlp(refElPol%Nnodes2D*phys%Neq),indgt(refElTor%Nfl*phys%Neq),indlt(refElTor%Nfl*phys%Neq)
#endif
#endif
  !   integer             :: cks, clock_rate, cke, clock_start, clock_end
  !   real*8              :: time_start, time_finish, tps, tpe

  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1, timing%clock_rate1)
  end if

  !********************************
  ! Save rhs for accuracy check
  !********************************
  IF (switch%ckeramp .and. lssolver%sollib .eq. 1) THEN
    errsol = 0.
    pertamp = 1.e-6
    ALLOCATE (rhspert(matK%n))
    CALL RANDOM_NUMBER(rhspert)               ! rhspert is between 0-1
    rhspert = (1.+pertamp*2.*(0.5 - rhspert)) ! now it is between 1-eps and 1+eps
    rhspert = rhspert*rhs%vals
  END IF

  !********************************
  ! Compute the face solution
  !********************************
  IF (lssolver%sollib .eq. 1) THEN
#ifdef WITH_PASTIX
    ! Solver dependent part-->PASTIX
    IF (matK%start) THEN
      call displayMatrixInfo()
      call init_mat_PASTIX(matPASTIX)
      call check_mat_PASTIX(matPASTIX)
      call anal_mat_PASTIX(matPASTIX)
      matK%start = .false.
    ELSE
      call build_mat_PASTIX(matPASTIX)
      !#ifdef PARALL
      ! This needs to be redone in parallel (who knows why??)
      call check_mat_PASTIX(matPASTIX)
      !#endif
    END IF
    call LU_mat_pastix(matPASTIX)
    call solve_mat_PASTIX(matPASTIX)
#else
    WRITE (6, *) "Trying to use PASTIX but compiled without option -DWITH_PASTIX"
    STOP
#endif
  ELSE IF (lssolver%sollib .eq. 2) THEN
#ifdef WITH_PSBLAS
    ! Solver dependent part-->PSBLAS
    IF (matK%start) THEN
      call displayMatrixInfo()
      call init_mat_PSBLAS(matPSBLAS, matK)
      call build_mat_PSBLAS(matPSBLAS, matK)
      call build_prec_PSBLAS(matPSBLAS)
      call fill_vec_PSBLAS(matPSBLAS)
      call solve_mat_PSBLAS(matPSBLAS)
      matK%start = .false.
    ELSE
      call build_mat_PSBLAS(matPSBLAS, matK)
      call build_prec_PSBLAS(matPSBLAS)
      call fill_vec_PSBLAS(matPSBLAS)
      call solve_mat_PSBLAS(matPSBLAS)
    END IF
#else
    WRITE (6, *) "Trying to use PSBLAS but compiled without option -DWITH_PSBLAS"
    STOP
#endif
  ENDIF


  !********************************
  ! Store face solution
  !********************************
#ifdef PARALL
  ALLOCATE (aux_sol(matK%n))
  IF (lssolver%sollib .eq. 1) THEN
#ifdef WITH_PASTIX
    aux_sol = matPASTIX%rhs
#endif
  ELSE
#ifdef WITH_PSBLAS
    aux_sol = matPSBLAS%x%get_vect()
#endif
  ENDIF
#ifdef TOR3D
  ! ********* Parallel 3D ***********
  IF (MPIvar%ntor .gt. 1) THEN
    ntorass = numer%ntor/MPIvar%ntor
  ELSE
    ntorass = numer%ntor
  ENDIF

  Neq = phys%Neq
  Nfl = refElTor%Nfl
  N2d = Mesh%Nelems
  Np2d = refElPol%Nnodes2D
  Nfaces = Mesh%Nfaces
  Nfdir = Mesh%Ndir
  Nghostf = Mesh%nghostfaces
  Nghoste = Mesh%nghostelems
  IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .eq. MPIvar%ntor) THEN
    ct = 0
    dd = 1 + ntorass*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
    ddl = 1
    DO i = 1, N2D
      IF (Mesh%ghostelems(i) .eq. 1) CYCLE
      ct = ct + 1
      indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
    END DO
  END IF
  DO itor = 1, ntorass
    ct = 0
    dd = 1 + (itor - 1)*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
    ddl = 1 + (itor - 1)*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
    IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) THEN
      ddl = ddl - (N2D-Nghoste)*Np2d*Neq
    END IF
    DO i = 1, N2D
      IF (Mesh%ghostelems(i) .eq. 1) CYCLE
      IF (MPIvar%ntor .gt. 1 .and. itor == 1) CYCLE
      ct = ct + 1
      indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
    END DO
    ct = 0
    dd = dd + (N2D*Np2D)*Neq
    ddl = ddl + ((N2D-Nghoste)*Np2D)*Neq

    DO i = 1, Mesh%Nfaces
      IF (Mesh%ghostfaces(i) .eq. 1) CYCLE
      ct = ct + 1
      indgt = dd + (i - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
      indlt = ddl + (ct - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
      sol%u_tilde(indgt) = (1.-numer%dumpnr)*sol%u_tilde(indgt) + numer%dumpnr*aux_sol(indlt)
    END DO
  END DO
  IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) THEN
    ct = 0
    dd = 1 + ntorass*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
    ddl = 1 + ntorass*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
    ddl = ddl - (N2D-Nghoste)*Np2d*Neq
    DO i = 1, N2D
      IF (Mesh%ghostelems(i) .eq. 1) CYCLE
      ct = ct + 1
      indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
      sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
    END DO
  END IF

#else
  ! ********* Parallel 2D ***********
  ct = 0
  Neq = phys%Neq
  Nfp = Mesh%Nnodesperface
  DO i = 1, Mesh%Nfaces
    IF (Mesh%ghostfaces(i) .eq. 1) CYCLE
    ct = ct + 1
    indg = (i - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
    indl = (ct - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
    sol%u_tilde(indg) = (1.-numer%dumpnr)*sol%u_tilde(indg) + numer%dumpnr*aux_sol(indl)
  END DO
#endif

  DEALLOCATE (aux_sol)

#else
  ! ********* Sequential 2D and 3D ***********
  IF (lssolver%sollib .eq. 1) THEN
#ifdef WITH_PASTIX
    sol%u_tilde(1:matK%n) = (1.-numer%dumpnr)*sol%u_tilde(1:matK%n) + numer%dumpnr*matPASTIX%rhs
#else
    WRITE (6, *) "Trying to use PASTIX but compiled without option -DWITH_PASTIX"
    STOP
#endif
  ELSE IF (lssolver%sollib .eq. 2) THEN
#ifdef WITH_PSBLAS
    sol%u_tilde(1:matK%n) = (1.-numer%dumpnr)*sol%u_tilde(1:matK%n) + numer%dumpnr*matPSBLAS%x%get_vect()
#else
    WRITE (6, *) "Trying to use PSBLAS but compiled without option -DWITH_PSBLAS"
    STOP
#endif
  END IF
#endif


  !#ifdef PARALL
  !if (MPIvar%ntor>1) then
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_par_ptor')
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop
  !else
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_par')
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop
  !end if
  !#else
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_seq')
  !stop
  !#endif

  !********************************
  ! Check accuracy of the solution
  !********************************
  IF (switch%ckeramp .and. lssolver%sollib .eq. 1) THEN
    ! use the perturbed rhs
    matPASTIX%rhs = rhspert

    ! solve again..
    CALL solve_mat_PASTIX(matPASTIX)
#ifdef PARALL
    ct = 0
    DO i = 1, Mesh%Nfaces
      IF (Mesh%ghostfaces(i) .eq. 1) CYCLE
      ct = ct + 1
      indg = (i - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
      indl = (ct - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
      errsol = max(errsol, maxval(abs(sol%u_tilde(indg) - matPASTIX%rhs(indl))))
    END DO
#else
    DO i = 1, matPASTIX%n
      errsol = max(errsol, abs(sol%u_tilde(i) - matPASTIX%rhs(i)))
    ENDDO
#endif
    errsol = errsol/pertamp
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, errsol, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) "ERROR AMPLIFICATION: ", errsol
    ENDIF
    DEALLOCATE (rhspert)
  END IF

  !**********************************************
  ! Deallocate matK
  !**********************************************
  CALL free_mat()

  !**********************************************
  ! MPI communications
  !**********************************************


#ifdef PARALL
  if (utils%timing) then
    call cpu_time(timing%tps2)
    call system_clock(timing%cks2, timing%clock_rate2)
  end if
  CALL exchangeSol()

  if (utils%timing) then
    call cpu_time(timing%tpe2)
    call system_clock(timing%cke2, timing%clock_rate2)
    timing%runtcom = timing%runtcom + (timing%cke2 - timing%cks2)/real(timing%clock_rate2)
    timing%cputcom = timing%cputcom + timing%tpe2 - timing%tps2
  end if
#endif


  !#ifdef PARALL
  !if (MPIvar%ntor>1) then
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_par_ptor')
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop
  !else
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_par')
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop
  !end if
  !#else
  !call HDF5_save_vector(sol%u_tilde,'u_tilde_seq')
  !stop
  !#endif

  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1, timing%clock_rate1)
    timing%runtglb = timing%runtglb + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputglb = timing%cputglb + timing%tpe1 - timing%tps1
#ifdef PARALL
    ! Take out the communication time from the solve time
    timing%runtglb = timing%runtglb -(timing%cke2 - timing%cks2)/real(timing%clock_rate2)
    timing%cputglb = timing%cputglb - (timing%tpe2 - timing%tps2)
#endif
  end if

CONTAINS
  subroutine displayMatrixInfo
    WRITE (6, *) " "
    WRITE (6, '(" *", 41("*"), "**")')
    WRITE (6, '(" *", 2X,  "Linear system size : ", I12, 6X, " *")') MatK%n
    WRITE (6, '(" *", 2X,  "Number of nnz      : ", I12, 6X, " *")') MatK%nnz
    WRITE (6, '(" *", 41("*"), "**")')
    WRITE (6, *) " "
  end subroutine displayMatrixInfo

END SUBROUTINE solve_global_system

