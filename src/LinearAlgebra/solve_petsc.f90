MODULE solve_petsc
#include "petsc/finclude/petsc.h"

  USE matrices_types
  USE types
  USE globals
  USE MPI_OMP
  USE petsc
  use initialization
  IMPLICIT NONE


TYPE PETSC_STRUC ! A type to store structures for the PETSC library
  Mat                               :: matK          ! sparse matrix for PSBLAS
  KSP                               :: ksp           ! solver
  PC                                :: pc            ! preconditioner
  Vec                               :: solPETSC_vec      ! solution
  Vec                               :: rhs_vec           ! right hand side
  PetscInt                          :: n, nnz        ! matrix n and number of non-zeros
  PetscInt, dimension(:), pointer   :: rowptr => null()
  PetscInt, dimension(:), pointer   :: loc2glob => null()
  PetscInt, dimension(:), pointer   :: cols => null()
  PetscScalar, dimension(:), pointer :: vals_matK => null()
  PetscScalar, dimension(:), pointer :: vals_rhs => null()
  PetscScalar                       :: residue   ! UNPRECONDTIONED residue norm of the linear system
  PetscInt                          :: its      ! Previous number of iterations

END TYPE

TYPE(PETSC_STRUC) :: matPETSC


CONTAINS

  !***********************************************
  ! Initialization of matrix
  ! Part specific to PETSc
  !***********************************************
  SUBROUTINE init_mat_PETSC(matPETSC)
    use petsc
    use MPI_OMP

    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    PetscBool             :: initialized
    PetscScalar           :: alpha= 0
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke



    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Link the PETSC instance with the CSR matrix storage
    IF (associated(matPETSC%rowptr)) then
      deallocate (matPETSC%rowptr)
    endif
    IF (associated(matPETSC%loc2glob)) then
      deallocate (matPETSC%loc2glob)
    endif
    IF (associated(matPETSC%cols)) then
      deallocate (matPETSC%cols)
    endif
    IF (associated(matPETSC%vals_matK)) then
      deallocate (matPETSC%vals_matK)
    endif
    IF (associated(matPETSC%vals_rhs)) then
      deallocate (matPETSC%vals_rhs)
    endif

    matPETSC%n = matK%n
    matPETSC%nnz = matK%nnz
    ALLOCATE (matPETSC%rowptr(matPETSC%n + 1))
    ALLOCATE (matPETSC%loc2glob(matPETSC%n))
    ALLOCATE (matPETSC%cols(matPETSC%nnz))
    ALLOCATE (matPETSC%vals_matK(matPETSC%nnz))
    ALLOCATE (matPETSC%vals_rhs(matPETSC%n))
    matPETSC%rowptr    = matK%rowptr
    matPETSC%cols      = matK%cols
    matPETSC%vals_matK = matK%vals
    matPETSC%loc2glob  = matK%loc2glob
    matPETSC%vals_rhs  = rhs%vals


    IF (matK%start) THEN
      matPETSC%residue   = 0
      matPETSC%its       = 0
#ifdef PARALL
        ! create matrix and set it up
        call MatCreate(MPI_COMM_WORLD,matPETSC%matK , ierr)
        call MatSetType(matPETSC%matK, MATMPIAIJ, ierr)
        call MatSetSizes(matPETSC%matK, matPETSC%n,matPETSC%n, PETSC_DETERMINE,PETSC_DETERMINE,ierr);
        call MatSetUp(matPETSC%matK, ierr)

        ! create rhs and solution vector
        call MatCreateVecs(matPETSC%matK, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr)
        call VecAssemblyBegin(matPETSC%rhs_vec,ierr)
        call VecAssemblyEnd(matPETSC%rhs_vec,ierr)
        call VecAssemblyBegin(matPETSC%solPETSC_vec,ierr)
        call VecAssemblyEnd(matPETSC%solPETSC_vec,ierr)

#else
      ! create matrix and set it up
      call MatCreate(MPI_COMM_WORLD,matPETSC%matK , ierr)
      call MatSetType(matPETSC%matK, MATSEQAIJ, ierr)
      call MatSetSizes(matPETSC%matK, matPETSC%n,matPETSC%n, matPETSC%n,matPETSC%n,ierr);
      call MatSetUp(matPETSC%matK, ierr)

      ! Create solution vector / initial guess also
      call MatCreateVecs(matPETSC%matK, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr)
      call VecAssemblyBegin(matPETSC%rhs_vec,ierr)
      call VecAssemblyEnd(matPETSC%rhs_vec,ierr)
      call VecAssemblyBegin(matPETSC%solPETSC_vec,ierr)
      call VecAssemblyEnd(matPETSC%solPETSC_vec,ierr)
#endif
! create the ksp instance
    call KSPCreate(PETSC_COMM_WORLD, matPETSC%ksp, ierr);
    ENDIF

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime1 = timing%rlstime1 + (cke - cks)/real(clock_rate)
      timing%clstime1 = timing%clstime1 + tpe - tps
    end if
  END SUBROUTINE init_mat_PETSC


  SUBROUTINE preallocate_matrix_PETSC(matPETSC)
    use petsc
    TYPE(PETSC_STRUC)      :: matPETSC
    Integer, allocatable   :: counter_diag(:), counter_offdiag(:), istart_vec(:), iend_vec(:)
    PetscInt               :: istart, iend, r, rr, rrr,  M, N, I, J
    PetscScalar            :: K
    PetscErrorCode         :: ierr

    ! Get global sizes
    call MatGetSize(matPETSC%matK,M,N, ierr)

    ! Get local row ownership of the matrix
    call MatGetOwnershipRange(matPETSC%matK, istart, iend, ierr)

    ! allocate diagonal and offdiagonal non-zero counters and a vector containing the beginning and ends!
    ! of the row ownership, per process. Each process owns X rows of the matrix, a slice of matrix.
    allocate(counter_diag(N))
    allocate(counter_offdiag(N))
    allocate(istart_vec(MPIvar%glob_size))
    allocate(iend_vec(MPIvar%glob_size))

    counter_diag = 0
    counter_offdiag = 0
    istart_vec = 0
    iend_vec = 0
    istart_vec(MPIvar%glob_id+1) = istart
    iend_vec(MPIvar%glob_id+1) = iend

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,istart_vec,MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,iend_vec,MPIvar%glob_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    ! count non-zero elements
    DO r = 1, matPETSC%n
        DO rr = matPETSC%rowptr(r), (matPETSC%rowptr(r+1)-1)
          I = matPETSC%loc2glob(r)-1;
          J = matPETSC%cols(rr)-1;

          ! Find out on which process' slice I am (when found, exit the loop, rrr will be kept)
          DO rrr = 1, MPIvar%glob_size
            IF(((I+1) .ge. istart_vec(rrr)) .and. ((I+1) .le. iend_vec(rrr))) THEN
                EXIT
            ENDIF
          ENDDO

          ! put in or out diagonal (always on diagonal for serial)
          IF((J .ge. istart_vec(rrr)) .and. (J .le. iend_vec(rrr)-1)) THEN
              counter_diag(I+1) = counter_diag(I+1) + 1
            ELSE
              counter_offdiag(I+1) = counter_offdiag(I+1) + 1
            ENDIF
        END DO
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,counter_diag,N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,counter_offdiag,N, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MatMPIAIJSetPreallocation(matPETSC%matK,PETSC_DECIDE,counter_diag(istart+1:iend),PETSC_DECIDE,counter_offdiag(istart+1:iend), ierr); CHKERRA(ierr)
#else
    call MatSeqAIJSetPreallocation(matPETSC%matK, PETSC_DECIDE, counter_diag, ierr)
#endif

    DEALLOCATE(counter_diag)
    DEALLOCATE(counter_offdiag)
    DEALLOCATE(istart_vec)
    DEALLOCATE(iend_vec)
  END SUBROUTINE preallocate_matrix_PETSC

  !***********************************************
  ! From the generic CSR storage, fills the PETSC
  ! instance of the matrix
  !***********************************************
  SUBROUTINE build_mat_PETSC(matPETSC)
    TYPE(MAT_CSR_TYP) :: matK
    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    PetscInt              :: one = 1, r, rr
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke

    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! set the values in the matrix
    DO r = 1, matPETSC%n
        DO rr = matPETSC%rowptr(r), (matPETSC%rowptr(r+1)-1)
          call MatSetValues(matPETSC%matK,one, matPETSC%loc2glob(r)-1, one, matPETSC%cols(rr)-1, matPETSC%vals_matK(rr), INSERT_VALUES, ierr); CHKERRA(ierr)
        END DO
    END DO

    call MatAssemblyBegin(matPETSC%matK,MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(matPETSC%matK,MAT_FINAL_ASSEMBLY, ierr)

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime2 = timing%rlstime2 + (cke - cks)/real(clock_rate)
      timing%clstime2 = timing%clstime2 + tpe - tps
    end if
  END SUBROUTINE build_mat_PETSC


  SUBROUTINE fill_vec_PETSC(matPETSC)
    TYPE(PETSC_STRUC)     :: matPETSC
    PetscErrorCode        :: ierr
    PetscInt              :: one = 1, r, rr
    real*8                :: tps, tpe
    integer               :: cks, clock_rate, cke
    PetscInt              :: istart, iend
    PetscScalar           :: alpha= 0
    PetscViewer           :: viewer



    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! set values
    call VecSetValues(matPETSC%rhs_vec,matPETSC%n, matPETSC%loc2glob-1,rhs%vals, INSERT_VALUES, ierr);

    call VecAssemblyBegin(matPETSC%rhs_vec,ierr)
    call VecAssemblyEnd(matPETSC%rhs_vec,ierr)


    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime3 = timing%rlstime3 + (cke - cks)/real(clock_rate)
      timing%clstime3 = timing%clstime3 + tpe - tps
    end if
  END SUBROUTINE fill_vec_PETSC


  !***********************************************
  ! Solve problem with PETSC
  !***********************************************
  SUBROUTINE solve_mat_PETSC(matPETSC, ir)
    TYPE(PETSC_STRUC)     :: matPETSC
    Integer, INTENT(IN)   :: ir
    integer               :: Nall, irow, iproc
    real*8                :: tps, tpe
    integer               :: cks, clock_rate, cke
    PetscViewer           :: viewer
    KSPType               :: ksp_type
    PCType                :: pc_type
    PetscViewerAndFormat  :: vf
    PetscErrorCode        :: ierr

    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Build preconditioner based on linear system matrix
    call KSPSetOperators(matPETSC%ksp, matPETSC%matK, matPETSC%matK, ierr);

    ! set or not the initial guess to zero
    IF(.not. lssolver%igz)THEN
        call KSPSetInitialGuessNonzero(matPETSC%ksp, PETSC_TRUE, ierr)
    ENDIF

    ! set solver and preconditione type from param.txt
    call set_ksp_pc_type(ksp_type,pc_type)
    call KSPSetType(matPETSC%ksp,ksp_type, ierr);

    ! set the norm type
    IF(lssolver%kspnorm .eq. 1) THEN
      call KSPSetNormType(matPETSC%ksp,KSP_NORM_PRECONDITIONED, ierr)
    ELSEIF(lssolver%kspnorm .eq. 2) THEN
      call KSPSetNormType(matPETSC%ksp,KSP_NORM_UNPRECONDITIONED, ierr)
    ENDIF

    ! set tollerances
    call KSPSetTolerances(matPETSC%ksp, lssolver%rtol, lssolver%atol, PETSC_DEFAULT_REAL, lssolver%kspitmax, ierr)

    ! get the preconditioner and set it
    call KSPGetPC(matPETSC%ksp, matPETSC%pc, ierr);
    call PCSetType(matPETSC%pc, pc_type, ierr);
    ! in case of PCMG set the additional info
    select case (lssolver%pctype)
      case ('PCMG')
        call PCMGSetLevels(matPETSC%pc,lssolver%mglevels,MPI_COMM_WORLD, ierr)
        IF    (lssolver%mgtypeform .eq. 1) THEN
          call PCMGSetType(matPETSC%pc,PC_MG_MULTIPLICATIVE, ierr)
        ELSEIF(lssolver%mgtypeform .eq. 2) THEN
          call PCMGSetType(matPETSC%pc,PC_MG_ADDITIVE, ierr)
        ELSEIF(lssolver%mgtypeform .eq. 3) THEN
          call PCMGSetType(matPETSC%pc,PC_MG_FULL, ierr)
        ELSEIF(lssolver%mgtypeform .eq. 4) THEN
          call PCMGSetType(matPETSC%pc,PC_MG_KASKADE, ierr)
        ENDIF
    end select

    ! set preconditioner up
    call PCSetUp(matPETSC%pc, ierr);
    call KSPSetUp(matPETSC%ksp, ierr);

    ! print out the residue
    IF(lssolver%kspitrace .eq. 1) THEN
      call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,vf,ierr)
      call KSPMonitorSet(matPETSC%ksp,KSPMonitorTrueResidual,vf,PetscViewerAndFormatDestroy,ierr)
    ENDIF

    ! reuse or not preconditioner
    IF(lssolver%rprecond .eq. 1) THEN
      call KSPSetReusePreconditioner(matPETSC%ksp,PETSC_TRUE,ierr)
    ELSEIF(lssolver%rprecond .eq. 2) THEN
      IF(MOD(ir, lssolver%Nrprecond) .eq. 0) THEN
        call KSPSetReusePreconditioner(matPETSC%ksp,PETSC_FALSE,ierr)
        IF (MPIvar%glob_id .eq. 0) then
            WRITE(*,*) 'PETSC is re-computing preconditioner'
        ENDIF
      ENDIF
    ELSEIF(lssolver%rprecond .eq. 3) THEN
      IF(matPETSC%its .eq. lssolver%kspitmax) THEN
        call KSPSetReusePreconditioner(matPETSC%ksp,PETSC_FALSE,ierr)
        IF (MPIvar%glob_id .eq. 0) then
            WRITE(*,*) 'PETSC is re-computing preconditioner'
        ENDIF
      ENDIF
    ELSE
      IF (MPIvar%glob_id .eq. 0) then
          WRITE(*,*) 'Wrong choice, value must be between 1, 2 or 3, stopping. '
          STOP
      ENDIF
    ENDIF

    !call KSPView(matPETSC%ksp,PETSC_VIEWER_STDOUT_WORLD, ierr)
    call KSPSolve(matPETSC%ksp, matPETSC%rhs_vec, matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)

    call KSPGetIterationNumber(matPETSC%ksp,matPETSC%its,ierr)

    ! compute the unpreconditioned relative residue
    call compute_residue

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime4 = timing%rlstime4 + (cke - cks)/real(clock_rate)
      timing%clstime4 = timing%clstime4 + tpe - tps
    end if
  END SUBROUTINE solve_mat_PETSC


  SUBROUTINE PETSC_retrieve_array(vector, aux_sol, size)
    Vec, INTENT(IN)           :: vector
    Integer, INTENT(IN)       :: size
    PetscScalar, INTENT(INOUT)     :: aux_sol(size)
    PetscScalar, pointer      :: data_array(:)
    PetscErrorCode            :: ierr
    ! copy the vector values into the array
    call VecGetArrayF90(vector,data_array,ierr)
    aux_sol = data_array
    call VecRestoreArrayF90(vector,data_array,ierr)
  END SUBROUTINE PETSC_retrieve_array

  SUBROUTINE compute_residue
    Vec                       :: true_residue_vec
    PetscReal                 :: alpha=-1
    PetscScalar               :: true_residue, norm
    PetscErrorCode            :: ierr

    call MatCreateVecs(matPETSC%matK, true_residue_vec, true_residue_vec, ierr)
    call VecScale(matPETSC%rhs_vec, alpha, ierr)
    call MatMultAdd(matPETSC%matK, matPETSC%solPETSC_vec,matPETSC%rhs_vec, true_residue_vec, ierr)
    call VecScale(matPETSC%rhs_vec, alpha, ierr)
    call VecNorm(true_residue_vec, NORM_2, true_residue, ierr)
    call VecNorm(matPETSC%rhs_vec, NORM_2, norm, ierr)
    matPETSC%residue = true_residue/norm
    call VecDestroy(true_residue_vec, ierr)

  END SUBROUTINE compute_residue

  SUBROUTINE set_ksp_pc_type(ksp_type, pc_type)
    KSPType, INTENT(OUT)       :: ksp_type
    PCType, INTENT(OUT)       :: pc_type

    select case (lssolver%kspmethd)
      case ('KSPCG')
        ksp_type = KSPCG
      case ('KSPBCGS')
        ksp_type = KSPBCGS
      case ('KSPGMRES')
        ksp_type = KSPGMRES
      case ('KSPMINRES')
        ksp_type = KSPMINRES
      case ('KSPPREONLY')
        ksp_type = KSPPREONLY
      case ('KSPRICHARDSON')
        ksp_type = KSPRICHARDSON
      case default
        IF (MPIvar%glob_id .eq. 0) THEN
          WRITE(*,*) "KSP METHOD NOT IN LIST, add it and recompile"
          STOP
        ENDIF
      end select

      select case (lssolver%pctype)
        case ('PCNONE')
          pc_type = PCNONE
        case ('PCJACOBI')
          pc_type = PCJACOBI
        case ('PCSOR')
          pc_type = PCSOR
        case ('PCLU')
          pc_type = PCLU
        case ('PCILU')
          pc_type = PCILU
        case ('PCGAMG')
          pc_type = PCGAMG
        case ('PCMG')
          pc_type = PCMG
        case default
          IF (MPIvar%glob_id .eq. 0) THEN
            WRITE(*,*) "PC TYPE NOT IN LIST, add it and recompile"
            STOP
          ENDIF
        end select
  END SUBROUTINE set_ksp_pc_type



  SUBROUTINE terminate_PETSC()
    PetscErrorCode        :: ierr
    IF (associated(matPETSC%rowptr)) then
      deallocate (matPETSC%rowptr)
    endif
    IF (associated(matPETSC%loc2glob)) then
      deallocate (matPETSC%loc2glob)
    endif
    IF (associated(matPETSC%cols)) then
      deallocate (matPETSC%cols)
    endif
    IF (associated(matPETSC%vals_rhs)) then
      deallocate (matPETSC%vals_rhs)
    endif
    IF (associated(matPETSC%vals_matK)) then
      deallocate (matPETSC%vals_matK)
    endif

    ! destroy everything
    call VecDestroy(matPETSC%rhs_vec, ierr); CHKERRA(ierr)
    call VecDestroy(matPETSC%solPETSC_vec, ierr); CHKERRA(ierr)
    call MatDestroy(matPETSC%matK, ierr); CHKERRA(ierr)
    call KSPDestroy(matPETSC%ksp, ierr); CHKERRA(ierr)
  END SUBROUTINE terminate_PETSC
END MODULE solve_petsc
