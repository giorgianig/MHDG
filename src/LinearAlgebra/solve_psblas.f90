MODULE solve_psblas
  USE matrices_types
  USE types
  USE globals
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use mpi_omp
#ifdef WITH_MLD2P4
  use mld_prec_mod
#endif
  IMPLICIT NONE

  TYPE PSBLAS_STRUC ! A type to store structures for the PSBLAS library
    type(psb_dspmat_type) :: mat          ! sparse matrix for PSBLAS
#ifndef WITH_MLD2P4
    type(psb_dprec_type)  :: prec         ! preconditioner
#else
    type(mld_dprec_type)  :: prec
#endif
    type(psb_desc_type)   :: desc_a       ! descriptor
    type(psb_d_vect_type) :: x, b         ! dense vector
    integer(psb_ipk_)     :: ictxt, iam, np, info   ! parallel environment
  END TYPE

  TYPE(PSBLAS_STRUC) :: matPSBLAS

CONTAINS

  !                real*8  :: tps, tpe
  !                integer :: cks,clock_rate,cke
  !                if (lssolver%timing) then
  !                                call cpu_time(tps)
  !                                call system_clock(cks,clock_rate)
  !                endif

  !                if (lssolver%timing) then
  !                        call cpu_time(tpe)
  !                        call system_clock(cke,clock_rate)
  !                        timing%rlstime1 = timing%rlstime1+(cke-cks)/real(clock_rate)
  !                        timing%clstime1 = timing%clstime1+tpe-tps
  !                end if

  !***********************************************
  ! Initialization of matrix instance of the
  ! PSBLAS library
  !***********************************************
  SUBROUTINE init_mat_PSBLAS(matPSBLAS, matK)
    TYPE(MAT_CSR_TYP) :: matK
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_lpk_), allocatable   :: ind(:)
    integer(psb_lpk_), allocatable   :: rows(:), cols(:)
    integer(psb_ipk_)               :: info
    integer(psb_ipk_)               :: iam, np, ncoef
    integer                         :: irow, i
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    !******************************************
    !Initialize parallel environment for psblas
    !******************************************
    call psb_init(matPSBLAS%ictxt)

    ! get information about the PSBLAS parallel environment
    ! psb_info: in{icontxt} out{iam,np}
    ! iam is the current process
    ! np number of processes
    call psb_info(matPSBLAS%ictxt, iam, np)

    if (iam < 0) then
      ! This should not happen, but just in case
      call psb_exit(matPSBLAS%ictxt)
      stop
    endif
    !  if(psb_get_errstatus() /= 0) goto 9999
    !  name='mld_s_pde3d'
    call psb_set_errverbosity(itwo)

    !******************************************
    !  Allocate a communication descriptor
    !******************************************
    !   call psb_cdall(matPSBLAS%ictxt,matPSBLAS%desc_a,info,nl=matK%n)

    allocate (rows(matK%n))
    rows = matK%loc2glob
    call psb_cdall(matPSBLAS%ictxt, matPSBLAS%desc_a, info, vl=rows)
    deallocate (rows)

    !******************************************
    !  Create the communication descriptor
    !******************************************
    do irow = 1, matK%n
      ncoef = matK%rowptr(irow + 1) - matK%rowptr(irow)
      allocate (rows(ncoef), ind(ncoef), cols(ncoef))
      ind = (/(i, i=matK%rowptr(irow), matK%rowptr(irow + 1) - 1)/)
      rows = matK%loc2glob(irow)
      cols = matK%cols(ind)
      call psb_cdins(ncoef, rows, cols, matPSBLAS%desc_a, info)
      if (info .ne. psb_success_) write (6, *) "failed filling matrix!"
      deallocate (rows, ind, cols)
    end do

    !******************************************
    !  Assembly the communication descriptor
    !******************************************
    call psb_cdasb(matPSBLAS%desc_a, info)

    !******************************************
    !  Allocate matrix,unknowns and rhs
    !******************************************
    if (info == psb_success_) call psb_spall(matPSBLAS%mat, matPSBLAS%desc_a, info, nnz=matK%nnz)
    if (info == psb_success_) call psb_geall(matPSBLAS%x, matPSBLAS%desc_a, info)
    if (info == psb_success_) call psb_geall(matPSBLAS%b, matPSBLAS%desc_a, info)

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime1 = timing%rlstime1 + (cke - cks)/real(clock_rate)
      timing%clstime1 = timing%clstime1 + tpe - tps
    end if
  END SUBROUTINE init_mat_PSBLAS

  !***********************************************
  ! Here the matrix is passed to the PSBLAS
  ! library
  !***********************************************
  SUBROUTINE build_mat_PSBLAS(matPSBLAS, matK)
    TYPE(MAT_CSR_TYP) :: matK
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_lpk_), allocatable   :: ind(:)
    integer(psb_lpk_), allocatable   :: rows(:), cols(:)
    real(psb_dpk_), allocatable     :: vals(:)
    integer(psb_ipk_)               :: info, irow, i, ncoef
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    !******************************************
    !  Initialize matrix to zero
    !******************************************
    call psb_sprn(matPSBLAS%mat, matPSBLAS%desc_a, info)

    !******************************************
    !  Fill matrix
    !******************************************
    do irow = 1, matK%n
      ncoef = matK%rowptr(irow + 1) - matK%rowptr(irow)
      allocate (rows(ncoef), ind(ncoef), cols(ncoef), vals(ncoef))
      ind = (/(i, i=matK%rowptr(irow), matK%rowptr(irow + 1) - 1)/)
      rows = matK%loc2glob(irow)
      cols = matK%cols(ind)
      vals = matK%vals(ind)
      call psb_spins(ncoef, rows, cols, vals, matPSBLAS%mat, matPSBLAS%desc_a, info)
      if (info .ne. psb_success_) write (6, *) "failed filling matrix!"
      deallocate (rows, ind, cols, vals)
    end do

    !******************************************
    !  Assembly matrix
    !******************************************
    call psb_spasb(matPSBLAS%mat, matPSBLAS%desc_a, info, dupl=psb_dupl_err_)
    if (info .ne. psb_success_) then
      write (6, *) "failed assemblying the matrix"
      stop
    endif

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime2 = timing%rlstime2 + (cke - cks)/real(clock_rate)
      timing%clstime2 = timing%clstime2 + tpe - tps
    end if
  END SUBROUTINE build_mat_PSBLAS

  !***********************************************
  ! Build the preconditioner
  !
  !***********************************************
#ifndef WITH_MLD2P4
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_ipk_)  :: info
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    call matPSBLAS%prec%init(matPSBLAS%ictxt, lssolver%ptype, info)
    call matPSBLAS%prec%build(matPSBLAS%mat, matPSBLAS%desc_a, info)

    if (info .ne. psb_success_) then
      write (6, *) "Something wrong in the preconditioner construction"
      stop
    end if
    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime3 = timing%rlstime3 + (cke - cks)/real(clock_rate)
      timing%clstime3 = timing%clstime3 + tpe - tps
    end if
  END SUBROUTINE build_prec_PSBLAS
#else
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_ipk_)  :: info
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    call matPSBLAS%prec%init(matPSBLAS%ictxt, lssolver%ptype, info)
    select case (lssolver%ptype)
    case ('NONE')
      ! Do nothing, keep defaults

    case ('DIAG', 'JACOBI', 'GS', 'FBGS')
      ! 1-level sweeps from "outer_sweeps"
      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)

    case ('BJAC')
      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
      call matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
      call matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
      call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)

    case ('AS')
      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
      call matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
      call matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
      call matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
      call matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
      call matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
      call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)

    case ('ML')
      ! multilevel preconditioner

      call matPSBLAS%prec%set('ml_cycle', lssolver%mlcycle, info)
      call matPSBLAS%prec%set('outer_sweeps', lssolver%outer_sweeps, info)
      if (lssolver%csize > 0)&
        & call matPSBLAS%prec%set('min_coarse_size', lssolver%csize, info)
      if (lssolver%mncrratio > 1)&
        & call matPSBLAS%prec%set('min_cr_ratio', lssolver%mncrratio, info)
      if (lssolver%maxlevs > 0)&
        & call matPSBLAS%prec%set('max_levs', lssolver%maxlevs, info)
      if (lssolver%athres >= dzero) &
        & call matPSBLAS%prec%set('aggr_thresh', lssolver%athres, info)
      !    if (lssolver%thrvsz>0) then
      !      do k=1,min(lssolver%thrvsz,size(prec%precv)-1)
      !        call matPSBLAS%prec%set('aggr_thresh',     lssolver%athresv(k),  info,ilev=(k+1))
      !      end do
      !    end if

      call matPSBLAS%prec%set('aggr_prol', lssolver%aggr_prol, info)
      call matPSBLAS%prec%set('par_aggr_alg', lssolver%par_aggr_alg, info)
      call matPSBLAS%prec%set('aggr_ord', lssolver%aggr_ord, info)
      call matPSBLAS%prec%set('aggr_filter', lssolver%aggr_filter, info)

      if ((psb_toupper(lssolver%smther) /= 'ML').and.(psb_toupper(lssolver%smther) /= 'NONE')) then
        call matPSBLAS%prec%set('smoother_type', lssolver%smther, info)
      end if

      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)

      select case (psb_toupper(lssolver%smther))
      case ('NONE', 'DIAG', 'JACOBI', 'GS', 'FBGS', 'ML')
        ! do nothing
      case ('BJAC')
        call matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
      case ('AS')
        call matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
        call matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
        call matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
        call matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
      case default
        call matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
        call matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
        call matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
        call matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
      end select

      if ((psb_toupper(lssolver%smther2) /= 'ML').and.(psb_toupper(lssolver%smther2) /= 'NONE')) then
        call matPSBLAS%prec%set('smoother_type', lssolver%smther2, info, pos='post')
      end if

      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps2, info, pos='post')

      select case (psb_toupper(lssolver%smther2))
      case ('NONE', 'DIAG', 'JACOBI', 'GS', 'FBGS','ML')
        ! do nothing
      case ('BJAC')
        call matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
      case ('AS')
        call matPSBLAS%prec%set('sub_ovr', lssolver%novr2, info, pos='post')
        call matPSBLAS%prec%set('sub_restr', lssolver%restr2, info, pos='post')
        call matPSBLAS%prec%set('sub_prol', lssolver%prol2, info, pos='post')
        call matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
      case default
        call matPSBLAS%prec%set('sub_restr', lssolver%restr2, info, pos='post')
        call matPSBLAS%prec%set('sub_prol', lssolver%prol2, info, pos='post')
        call matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
        call matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
        call matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
      end select

      if (psb_toupper(lssolver%csolve) /= 'NONE') then
        call matPSBLAS%prec%set('coarse_solve', lssolver%csolve, info)
        if (psb_toupper(lssolver%csolve) == 'BJAC') &
          &  call matPSBLAS%prec%set('coarse_subsolve', lssolver%csbsolve, info)
        call matPSBLAS%prec%set('coarse_mat', lssolver%cmat, info)
        call matPSBLAS%prec%set('coarse_fillin', lssolver%cfill, info)
        call matPSBLAS%prec%set('coarse_iluthrs', lssolver%cthres, info)
        call matPSBLAS%prec%set('coarse_sweeps', lssolver%cjswp, info)
      endif

    end select

    ! build the preconditioner
    call matPSBLAS%prec%hierarchy_build(matPSBLAS%mat, matPSBLAS%desc_a, info)
    call matPSBLAS%prec%smoothers_build(matPSBLAS%mat, matPSBLAS%desc_a, info)

    if (info .ne. psb_success_) then
      write (6, *) "Something wrong in the preconditioner construction"
      stop
    end if
    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime3 = timing%rlstime3 + (cke - cks)/real(clock_rate)
      timing%clstime3 = timing%clstime3 + tpe - tps
    end if
  END SUBROUTINE build_prec_PSBLAS
#endif

  !***********************************************
  ! Fill RHS and initial guess
  !
  !***********************************************
  SUBROUTINE fill_vec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_lpk_), allocatable   :: rows(:)
    integer(psb_ipk_)  :: info, n, i
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    !     call matPSBLAS%b%set(rhs%vals)
    !     call matPSBLAS%x%set(sol%u_tilde)

    n = matPSBLAS%mat%get_nrows()
    call matPSBLAS%x%zero()
    call matPSBLAS%b%zero()
    allocate (rows(n))
    rows = matK%loc2glob
    call psb_geins(n, rows, rhs%vals, matPSBLAS%b, matPSBLAS%desc_a, info)
#ifdef PARALL
    call getFirstGuessInParallel(matPSBLAS, rows)
#else
    call psb_geins(n, rows, sol%u_tilde, matPSBLAS%x, matPSBLAS%desc_a, info)
#endif
    deallocate (rows)

    !     n = matPSBLAS%mat%get_nrows()
    !     call matPSBLAS%x%zero()
    !     call matPSBLAS%b%zero()
    !     allocate(ind(n))
    !     ind = (/(i,i=1,n)/)
    !     call psb_geins(n,ind,rhs%vals,matPSBLAS%b,matPSBLAS%desc_a,info)
    !     call psb_geins(n,ind,sol%u_tilde,matPSBLAS%x,matPSBLAS%desc_a,info)
    !     deallocate(ind)

    !******************************************
    !  Assembly rhs and initial guess
    !******************************************

    ! Init. guess
    call psb_geasb(matPSBLAS%x, matPSBLAS%desc_a, info)
    if (info .ne. psb_success_) then
      write (6, *) "Something wrong in the initial guess assembly"
      stop
    end if

    ! RHS
    call psb_geasb(matPSBLAS%b, matPSBLAS%desc_a, info)
    if (info .ne. psb_success_) then
      write (6, *) "Something wrong in the rhs assembly"
      stop
    end if

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime4 = timing%rlstime4 + (cke - cks)/real(clock_rate)
      timing%clstime4 = timing%clstime4 + tpe - tps
    end if
  END SUBROUTINE fill_vec_PSBLAS

  !***********************************************
  ! Solve call for PSBLAS, assign output to
  ! u_tilde
  !***********************************************
  SUBROUTINE solve_mat_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_ipk_)     :: info   ! parallel environment
    integer(psb_ipk_)     :: iter, itmax, itrace, istopc, irst
    real(psb_dpk_)        :: err
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke,ierr
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif


    !write(6,*) "saving matrix"
    !!call mm_mat_write(matPSBLAS%mat,'sparse_matrix' ,ierr,1,'psblas_mat')
    !!call mm_array_write(matPSBLAS%b%v%v,'rhs',ierr,1,'psblas_b')
    !call mm_array_write(matPSBLAS%x%v%v,'solution_pastix',ierr,1,'psblas_solpastix')
    !!call mm_array_write(rhs%vals,ierr,1,'psblas_b')
    !!call mm_array_write(sol%u_tilde,ierr,1,'psblas_x0')
    !write(6,*) "done saving matrix"
    !stop
    !******************************************
    !  Solve
    !******************************************
    call psb_krylov(lssolver%kmethd, matPSBLAS%mat, matPSBLAS%prec, matPSBLAS%b, matPSBLAS%x, lssolver%tol, matPSBLAS%desc_a, info,&
      & itmax=lssolver%itmax, iter=iter, err=err, itrace=lssolver%itrace, istop=lssolver%istop, irst=lssolver%rest)

    !******************************************
    !  Store solution
    !******************************************
    !   sol%u_tilde = (1.-numer%dumpnr)*sol%u_tilde + numer%dumpnr*matPSBLAS%x%get_vect()
    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime5 = timing%rlstime5 + (cke - cks)/real(clock_rate)
      timing%clstime5 = timing%clstime5 + tpe - tps
    end if
  END SUBROUTINE solve_mat_PSBLAS

  !******************************************
  !  Free memory
  !******************************************
  SUBROUTINE terminate_PSBLAS()
    integer(psb_ipk_)     :: info

    call psb_gefree(matPSBLAS%b, matPSBLAS%desc_a, info)
    call psb_gefree(matPSBLAS%x, matPSBLAS%desc_a, info)
    call psb_spfree(matPSBLAS%mat, matPSBLAS%desc_a, info)
    call matPSBLAS%prec%free(info)
    call psb_cdfree(matPSBLAS%desc_a, info)
    call psb_exit(matPSBLAS%ictxt)
  END SUBROUTINE terminate_PSBLAS

#ifdef PARALL
  !********************************************
  !  Get the first guess in parallel
  !********************************************
  SUBROUTINE getFirstGuessInParallel(matPSBLAS, rows)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    integer(psb_lpk_), intent(in)   :: rows(:)
    integer(psb_ipk_)     :: info

    integer*4             :: i, j, n
#ifdef PARALL
    integer*4             :: ierr, Neq, Nfp
    integer*4             :: ct, indg(refElPol%Nfacenodes*phys%Neq), indl(refElPol%Nfacenodes*phys%Neq)
    real*8, allocatable    :: aux_sol(:)
#ifdef TOR3D
    integer*4             :: Nfl, N2d, Np2d, itor, Nfaces, Nfdir, Nghostf, Nghoste, dd, ddl, ntorass
    integer*4             :: indgp(refElPol%Nnodes2D*phys%Neq),indlp(refElPol%Nnodes2D*phys%Neq),indgt(refElTor%Nfl*phys%Neq),indlt(refElTor%Nfl*phys%Neq)
#endif
#endif

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
      n = Np2D*Neq
      DO i = 1, N2D
        IF (Mesh%ghostelems(i) .eq. 1) CYCLE
        ct = ct + 1
        indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        call psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)

        !        aux_sol(indlp) = sol%u_tilde(indgp)
      END DO
    END IF
    DO itor = 1, ntorass
      ct = 0
      dd = 1 + (itor - 1)*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
      ddl = 1 + (itor - 1)*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
      IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) THEN
        ddl = ddl - (N2D-Nghoste)*Np2d*Neq
      END IF
      n = Np2D*Neq
      DO i = 1, N2D
        IF (Mesh%ghostelems(i) .eq. 1) CYCLE
        IF (MPIvar%ntor .gt. 1 .and. itor == 1) CYCLE
        ct = ct + 1
        indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        call psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)
        !        sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
      END DO
      ct = 0
      dd = dd + (N2D*Np2D)*Neq
      ddl = ddl + ((N2D-Nghoste)*Np2D)*Neq

      n = Neq*Nfl
      DO i = 1, Mesh%Nfaces
        IF (Mesh%ghostfaces(i) .eq. 1) CYCLE
        ct = ct + 1
        indgt = dd + (i - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
        indlt = ddl + (ct - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
        call psb_geins(n, rows(indlt), sol%u_tilde(indgt), matPSBLAS%x, matPSBLAS%desc_a, info)
        !        sol%u_tilde(indgt) = (1.-numer%dumpnr)*sol%u_tilde(indgt) + numer%dumpnr*aux_sol(indlt)
      END DO
    END DO
    IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) THEN
      ct = 0
      dd = 1 + ntorass*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
      ddl = 1 + ntorass*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
      ddl = ddl - (N2D-Nghoste)*Np2d*Neq
      n = Np2D*Neq
      DO i = 1, N2D
        IF (Mesh%ghostelems(i) .eq. 1) CYCLE
        ct = ct + 1
        indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
        call psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)
        !        sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
      END DO
    END IF

#else
    ! ********* Parallel 2D ***********
    ct = 0
    Neq = phys%Neq
    Nfp = Mesh%Nnodesperface
    n = Neq*Nfp
    DO i = 1, Mesh%Nfaces
      IF (Mesh%ghostfaces(i) .eq. 1) CYCLE
      ct = ct + 1
      indg = (i - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
      indl = (ct - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
      call psb_geins(n, rows(indl), sol%u_tilde(indg), matPSBLAS%x, matPSBLAS%desc_a, info)
      !     sol%u_tilde(indg) = (1.-numer%dumpnr)*sol%u_tilde(indg) + numer%dumpnr*aux_sol(indl)
    END DO
#endif

  END SUBROUTINE getFirstGuessInParallel
#endif

END MODULE solve_psblas
