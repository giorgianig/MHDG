MODULE solve_pastix
  USE matrices_types
  USE types
  USE globals
  IMPLICIT NONE

#include "pastix_fortran.h"

  TYPE PASTIX_STRUC ! A type to store matrices in CSR format
    ! Variables are defined in a PASTIX specific way
    pastix_int_t :: n    ! Number of (locally owned) columns in the matrix
    pastix_int_t :: nnz  ! Number of (locally owned) non-zeros
    pastix_int_t, dimension(:), pointer    :: rowptr => null()  ! Index of first element of each row in cols and vals
    pastix_int_t, dimension(:), pointer    :: cols => null()  ! Column of each element
    pastix_float_t, dimension(:), pointer  :: vals => null()  ! Value of each element
    pastix_int_t, dimension(:), pointer    :: loc2glob => null()  ! Global column number of the local columns
    ! From here below the variables are PASTIX specific
    pastix_data_ptr_t                     :: pastix_data ! Structure to keep information in PaStiX (0 for first CALL)
    integer                               :: pastix_comm ! MPI communicator used by pastix
    pastix_int_t, dimension(:), pointer :: perm => null() ! permutation tabular
    pastix_int_t, dimension(:), pointer :: invp => null() ! reverse permutation tabular
    pastix_float_t, dimension(:), pointer :: rhs => null() ! Right hand side
    pastix_int_t                          :: nrhs   ! right hand side number
    pastix_int_t, dimension(:), pointer :: iparm => null()! Integer parameters
    real(float), dimension(:), pointer :: dparm => null()! Floating point parameters
    Integer                               :: nbthread    ! Number of threads in PaStiX
    logical                               :: start  ! keeps track if it is the first solve
  END TYPE

  integer, parameter                   :: verbose_level_PASTIX = 0 ! Level of verbose (0 = Silent mode, no messages; 1 = Some messages;
  ! 2 = Many messages; 3 = Like a gossip; 4 = Really talking too much...)
  TYPE(PASTIX_STRUC) :: matPASTIX

CONTAINS

  !***********************************************
  ! Initialization of matrix
  ! Part specific to PASTIX
  !***********************************************
  SUBROUTINE init_mat_PASTIX(matPASTIX)
    use MPI_OMP

    !   TYPE(MAT_CSR_TYP) :: matCSR
    TYPE(PASTIX_STRUC) :: matPASTIX
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Link the PASTIX instance with the CSR matrix storage
    IF (associated(matPASTIX%rowptr)) then
      deallocate (matPASTIX%rowptr)
    endif
    IF (associated(matPASTIX%loc2glob)) then
      deallocate (matPASTIX%loc2glob)
    endif
    IF (associated(matPASTIX%cols)) then
      deallocate (matPASTIX%cols)
    endif
    IF (associated(matPASTIX%vals)) then
      deallocate (matPASTIX%vals)
    endif
    IF (associated(matPASTIX%rhs)) then
      deallocate (matPASTIX%rhs)
    endif
    matPASTIX%n = matK%n
    matPASTIX%nnz = matK%nnz
    ALLOCATE (matPASTIX%rowptr(matPASTIX%n + 1))
    ALLOCATE (matPASTIX%loc2glob(matPASTIX%n))
    ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
    ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
    ALLOCATE (matPASTIX%rhs(matPASTIX%n))
    matPASTIX%rowptr = matK%rowptr
    matPASTIX%cols = matK%cols
    matPASTIX%vals = matK%vals
    matPASTIX%loc2glob = matK%loc2glob
    matPASTIX%rhs = rhs%vals

    IF (associated(matPASTIX%iparm)) then
      deallocate (matPASTIX%iparm)
    endif
    IF (associated(matPASTIX%dparm)) then
      deallocate (matPASTIX%dparm)
    endif
    IF (associated(matPASTIX%perm)) then
      deallocate (matPASTIX%perm)
    endif
    IF (associated(matPASTIX%invp)) then
      deallocate (matPASTIX%invp)
    endif
    ! Allocation of arrays specific to PASTIX
    ALLOCATE (matPASTIX%iparm(IPARM_SIZE))
    ALLOCATE (matPASTIX%dparm(DPARM_SIZE))
    ALLOCATE (matPASTIX%perm(matPASTIX%n))
    ALLOCATE (matPASTIX%invp(matPASTIX%n))

    ! Initializes the pastix instance
    matPASTIX%pastix_comm = MPI_COMM_WORLD
    matPASTIX%pastix_data = 0
    matPASTIX%nrhs = 1
    matPASTIX%iparm(IPARM_MODIFY_PARAMETER) = API_NO
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_INIT
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_INIT
#ifdef THREAD_FUNNELED
    matPASTIX%iparm(IPARM_THREAD_COMM_MODE) = API_THREAD_FUNNELED
#endif

    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
      matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)
    matPASTIX%iparm(IPARM_THREAD_NBR) = OMPvar%Nthreads
    matPASTIX%iparm(IPARM_VERBOSE) = verbose_level_PASTIX
    matPASTIX%iparm(IPARM_SYM) = API_SYM_NO
    matPASTIX%iparm(IPARM_FACTORIZATION) = API_FACT_LU
#ifdef THREAD_FUNNELED
    matPASTIX%iparm(IPARM_THREAD_COMM_MODE) = API_THREAD_FUNNELED
#endif

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime1 = timing%rlstime1 + (cke - cks)/real(clock_rate)
      timing%clstime1 = timing%clstime1 + tpe - tps
    end if
  END SUBROUTINE init_mat_PASTIX

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
  ! From the generic CSR storage, fills the PASTIX
  ! instance of the matrix
  !***********************************************
  SUBROUTINE build_mat_PASTIX(matPASTIX)
    !   TYPE(MAT_CSR_TYP) :: matCSR
    TYPE(PASTIX_STRUC) :: matPASTIX
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Link the PASTIX instance with the CSR matrix storage
    IF (associated(matPASTIX%rowptr)) then
      deallocate (matPASTIX%rowptr)
    ENDif
    IF (associated(matPASTIX%loc2glob)) then
      deallocate (matPASTIX%loc2glob)
    ENDif
    IF (associated(matPASTIX%cols)) then
      deallocate (matPASTIX%cols)
    ENDif
    IF (associated(matPASTIX%vals)) then
      deallocate (matPASTIX%vals)
    ENDif
    IF (associated(matPASTIX%rhs)) then
      deallocate (matPASTIX%rhs)
    endif
    matPASTIX%n = matK%n
    matPASTIX%nnz = matK%nnz
    ALLOCATE (matPASTIX%rowptr(matPASTIX%n + 1))
    ALLOCATE (matPASTIX%loc2glob(matPASTIX%n))
    ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
    ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
    ALLOCATE (matPASTIX%rhs(matPASTIX%n))
    matPASTIX%rowptr = matK%rowptr
    matPASTIX%cols = matK%cols
    matPASTIX%vals = matK%vals
    matPASTIX%loc2glob = matK%loc2glob
    matPASTIX%rhs = rhs%vals

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime4 = timing%rlstime4 + (cke - cks)/real(clock_rate)
      timing%clstime4 = timing%clstime4 + tpe - tps
    end if
  END SUBROUTINE build_mat_PASTIX

  !***********************************
  ! PASTIX specific part: consists
  ! only in checking the matrix
  !***********************************
  SUBROUTINE check_mat_PASTIX(matPASTIX)

    TYPE(PASTIX_STRUC) :: matPASTIX

    integer :: nnz
    pastix_data_ptr_t :: data_check
    pastix_int_t :: flagcor
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Removes multiple entries and checks symmetry
    nnz = matPASTIX%rowptr(matPASTIX%n + 1) - 1
    flagcor = API_YES
    CALL pastix_fortran_checkMatrix(data_check, matPASTIX%pastix_comm, matPASTIX%iparm(IPARM_VERBOSE), &
      matPASTIX%iparm(IPARM_SYM), flagcor, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%iparm(IPARM_DOF_NBR))
    IF (nnz .NE. (matPASTIX%rowptr(matPASTIX%n + 1) - 1)) then
      deallocate (matPASTIX%cols, matPASTIX%vals)
      matPASTIX%nnz = matPASTIX%rowptr(matPASTIX%n + 1) - 1
      ALLOCATE (matPASTIX%cols(matPASTIX%nnz))
      ALLOCATE (matPASTIX%vals(matPASTIX%nnz))
      CALL pastix_fortran_checkMatrix_END(data_check, matPASTIX%iparm(IPARM_VERBOSE), &
        matPASTIX%cols, matPASTIX%vals, matPASTIX%iparm(IPARM_DOF_NBR))
    endif

    ! The matrix does not need to be checked anymore
    matPASTIX%iparm(IPARM_MATRIX_VERIFICATION) = API_NO

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime2 = timing%rlstime2 + (cke - cks)/real(clock_rate)
      timing%clstime2 = timing%clstime2 + tpe - tps
    end if
  END SUBROUTINE check_mat_PASTIX

  !***********************************************
  ! Analysis of matrix with PASTIX
  !***********************************************
  SUBROUTINE anal_mat_PASTIX(matPASTIX)
    use matrices_tools, only: dump_CSR

    TYPE(PASTIX_STRUC) :: matPASTIX
    integer :: iproc, ierr
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Analysis phase
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_ANALYSE
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
      matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime3 = timing%rlstime3 + (cke - cks)/real(clock_rate)
      timing%clstime3 = timing%clstime3 + tpe - tps
    end if
  END SUBROUTINE anal_mat_PASTIX

  !***********************************************
  ! LU factorization of matrix with PASTIX
  !***********************************************
  SUBROUTINE LU_mat_pastix(matPASTIX)

    TYPE(PASTIX_STRUC) :: matPASTIX

    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Factorization phase
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_NUMFACT
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
      matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime5 = timing%rlstime5 + (cke - cks)/real(clock_rate)
      timing%clstime5 = timing%clstime5 + tpe - tps
    end if

  END SUBROUTINE LU_mat_pastix

  !***********************************************
  ! Solve problem with PASTIX
  !***********************************************
  SUBROUTINE solve_mat_PASTIX(matPASTIX)
    TYPE(PASTIX_STRUC) :: matPASTIX

    integer :: Nall, irow, iproc, IERR
    real*8  :: tps, tpe
    integer :: cks, clock_rate, cke
    if (lssolver%timing) then
      call cpu_time(tps)
      call system_clock(cks, clock_rate)
    endif

    ! Solve phase
    matPASTIX%iparm(IPARM_TRANSPOSE_SOLVE) = API_YES
    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_SOLVE
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_SOLVE
    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
      matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    if (lssolver%timing) then
      call cpu_time(tpe)
      call system_clock(cke, clock_rate)
      timing%rlstime6 = timing%rlstime6 + (cke - cks)/real(clock_rate)
      timing%clstime6 = timing%clstime6 + tpe - tps
    end if
  END SUBROUTINE solve_mat_PASTIX

  !SUBROUTINE free_mat_PASTIX(matPASTIX)
  !   TYPE(PASTIX_STRUC) :: matPASTIX
  !
  !   IF (associated(matPASTIX%rowptr)) then
  !      deallocate(matPASTIX%rowptr)
  !   endif
  !   IF (associated(matPASTIX%loc2glob)) then
  !      deallocate(matPASTIX%loc2glob)
  !   endif
  !   IF (associated(matPASTIX%cols)) then
  !      deallocate(matPASTIX%cols)
  !   endif
  !   IF (associated(matPASTIX%vals)) then
  !      deallocate(matPASTIX%vals)
  !   endif
  !   IF (associated(matPASTIX%rhs)) then
  !      deallocate(matPASTIX%rhs)
  !   endif
  !END SUBROUTINE free_mat_PASTIX

  SUBROUTINE terminate_mat_PASTIX()

    matPASTIX%iparm(IPARM_START_TASK) = API_TASK_CLEAN
    matPASTIX%iparm(IPARM_END_TASK) = API_TASK_CLEAN

    CALL dpastix_fortran(matPASTIX%pastix_data, matPASTIX%pastix_comm, &
      matPASTIX%n, matPASTIX%rowptr, matPASTIX%cols, matPASTIX%vals, &
      matPASTIX%loc2glob, matPASTIX%perm, matPASTIX%invp, matPASTIX%rhs, &
      matPASTIX%nrhs, matPASTIX%iparm, matPASTIX%dparm)

    IF (associated(matPASTIX%rowptr)) then
      deallocate (matPASTIX%rowptr)
    endif
    IF (associated(matPASTIX%loc2glob)) then
      deallocate (matPASTIX%loc2glob)
    endif
    IF (associated(matPASTIX%cols)) then
      deallocate (matPASTIX%cols)
    endif
    IF (associated(matPASTIX%vals)) then
      deallocate (matPASTIX%vals)
    endif
    IF (associated(matPASTIX%rhs)) then
      deallocate (matPASTIX%rhs)
    endif
    IF (associated(matPASTIX%iparm)) then
      deallocate (matPASTIX%iparm)
    endif
    IF (associated(matPASTIX%dparm)) then
      deallocate (matPASTIX%dparm)
    endif
    IF (associated(matPASTIX%perm)) then
      deallocate (matPASTIX%perm)
    endif
    IF (associated(matPASTIX%invp)) then
      deallocate (matPASTIX%invp)
    endif
  END SUBROUTINE terminate_mat_PASTIX

END MODULE solve_pastix
