MODULE matrices_tools

  USE matrices_types

  implicit none

CONTAINS

  !********************************
  ! full_to_CSR - GG January 2016
  ! Convert a full matrix into a
  ! sparse matrix in CSR format
  !********************************
  SUBROUTINE full_to_CSR(M, S)
    REAL, DIMENSION(:, :), INTENT(IN) :: M ! full matrix
    TYPE(MAT_CSR_TYP), INTENT(OUT)  :: S ! sparse matrix in CSR format
    INTEGER  :: r, c, nnz, ic, i, j
    REAL, PARAMETER :: tol = 1e-10

    ! Matrix dimensions
    r = SIZE(M, 1)
    c = SIZE(M, 1)

    ! Find nnz
    nnz = 0
    DO i = 1, r
      DO j = 1, c
        IF (ABS(M(i, j)) > tol) nnz = nnz + 1
      END DO
    END DO

    S%n = c
    S%nnz = nnz
    ALLOCATE (S%vals(1:nnz), S%cols(1:nnz))
    ALLOCATE (S%rowptr(1:r + 1))

    ! Fill sparse matrix
    ic = 0
    S%rowptr(1) = 1
    DO i = 1, r
      DO j = 1, c
        IF (ABS(M(i, j)) > tol) THEN
          ic = ic + 1
          S%vals(ic) = M(i, j)
          S%cols(ic) = j
        END IF
      END DO
      S%rowptr(i + 1) = ic + 1
    END DO

  END SUBROUTINE full_to_CSR

  !********************************
  ! CSR_to_full - GG January 2016
  ! Convert a sparse matrix in CSR format
  ! into a full matrix
  !********************************
  SUBROUTINE CSR_to_full(S, M)
    TYPE(MAT_CSR_TYP), INTENT(IN)    :: S ! sparse matrix in CSR format
    REAL(8), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: M ! full matrix
    INTEGER :: r, c, i, j, nnz
    INTEGER, DIMENSION(:), POINTER :: cols, rowptr
    REAL(8), DIMENSION(:), POINTER :: vals

    nnz = S%nnz
    cols => S%cols
    rowptr => S%rowptr
    vals => S%vals
    r = size(rowptr) - 1
    c = S%n

    ALLOCATE (M(1:r, 1:c))
    M = 0.
    DO i = 1, r
      DO j = rowptr(i), rowptr(i + 1) - 1
        M(i, cols(j)) = vals(j)
      END DO
    END DO
  END SUBROUTINE CSR_to_full

  !********************************
  ! Converts a sparse matrix in CSR
  ! format into its equivalent in
  ! element IJV format
  !********************************
  subroutine CSR_to_IJV(mat_CSR, mat_IJV)

    TYPE(MAT_CSR_TYP) :: mat_CSR
    TYPE(MAT_IJV_TYP) :: mat_IJV

    integer :: k

    ! The size of the arrays are the same
    mat_IJV%n = mat_CSR%n
    mat_IJV%nnz = mat_CSR%nnz

    ! Re-allocate the IJV matrix if it does not have the right size
    if (associated(mat_IJV%rows)) then
      if (size(mat_IJV%rows) .NE. mat_IJV%nnz) then
        deallocate (mat_IJV%rows)
        allocate (mat_IJV%rows(mat_IJV%nnz))
      endif
    else
      allocate (mat_IJV%rows(mat_IJV%nnz))
    endif
    if (associated(mat_IJV%cols)) then
      if (size(mat_IJV%cols) .NE. mat_IJV%nnz) then
        deallocate (mat_IJV%cols)
        allocate (mat_IJV%cols(mat_IJV%nnz))
      endif
    else
      allocate (mat_IJV%cols(mat_IJV%nnz))
    endif
    if (associated(mat_IJV%vals)) then
      if (size(mat_IJV%vals) .NE. mat_IJV%nnz) then
        deallocate (mat_IJV%vals)
        allocate (mat_IJV%vals(mat_IJV%nnz))
      endif
    else
      allocate (mat_IJV%vals(mat_IJV%nnz))
    endif

    ! Duplicate the arrays that are in common between the two formats
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
    do k = 1, mat_IJV%nnz
      mat_IJV%cols(k) = mat_CSR%cols(k)
      mat_IJV%vals(k) = mat_CSR%vals(k)
    enddo
    !$OMP END PARALLEL DO

    ! Generates the rows array
    call CSR_to_IJV_partial(mat_IJV%n, mat_CSR%rowptr, mat_CSR%loc2glob, mat_IJV%rows)

  end subroutine CSR_to_IJV

  !********************************
  ! Converts a sparse matrix in CSR
  ! format into its equivalent in
  ! element IJV format. In fact, in
  ! this partial version, only
  ! the row arrays are given and
  ! modified
  !********************************
  subroutine CSR_to_IJV_partial(N, rowsptr, loc2glob, Irows)

    integer :: N
    integer, dimension(:) :: rowsptr
    integer, dimension(:) :: Irows
    integer, dimension(:) :: loc2glob

    integer :: irow

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(irow)
    do irow = 1, N
      Irows(rowsptr(irow):rowsptr(irow + 1) - 1) = loc2glob(irow)
    enddo
    !$OMP END PARALLEL DO

  end subroutine CSR_to_IJV_partial

  !**************************************************************
  ! AMUX multiply a sparse matrices in  CSR format by a vector
  ! Inspired by Youcef Saad - GG January 2016
  !**************************************************************
  SUBROUTINE amux(A, x, y)
    TYPE(MAT_CSR_TYP), INTENT(IN)                     :: A
    REAL(8), DIMENSION(:), INTENT(IN)                  :: x
    REAL(8), DIMENSION(size(A%rowptr) - 1), INTENT(OUT)  :: y
    INTEGER  :: i, k, r
    REAL(8)  :: t

    r = size(A%rowptr) - 1
    DO i = 1, r
      t = 0.0D+00
      DO k = A%rowptr(i), A%rowptr(i + 1) - 1
        t = t + A%vals(k)*x(A%cols(k))
      END DO
      y(i) = t
    END DO
  END SUBROUTINE amux

  subroutine bsort2(w, ind, n, ncut)

    !*****************************************************************************80
    !
    !! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
    !
    !  Discussion:
    !
    !    This routine carries out a simple bubble sort for getting the NCUT largest
    !    elements in modulus, in array W.  IND is sorted accordingly.
    !    (Ought to be replaced by a more efficient sort especially
    !    if NCUT is not that small).
    !
    !  Modified:
    !
    !    07 January 2004
    !
    !  Author:
    !
    !    Youcef Saad
    !
    !  Parameters:
    !
    implicit none

    integer(kind=4) n

    integer(kind=4) i
    integer(kind=4) ind(*)
    integer(kind=4) iswp
    integer(kind=4) j
    integer(kind=4) ncut
    logical test
    real(kind=8) w(n)
    real(kind=8) wswp

    i = 1

    do

      test = .false.

      do j = n - 1, i, -1

        if (abs(w(j)) < abs(w(j + 1))) then
          !
          !  Swap.
          !
          wswp = w(j)
          w(j) = w(j + 1)
          w(j + 1) = wswp
          !
          !  Reorder the original ind array accordingly.
          !
          iswp = ind(j)
          ind(j) = ind(j + 1)
          ind(j + 1) = iswp
          !
          !  Set indicator that sequence is still unsorted.
          !
          test = .true.

        end if

      end do

      i = i + 1

      if (.not. test .or. ncut < i) then
        exit
      end if

    end do

    return
  end subroutine bsort2

  !**************************************************************
  ! LSOL solves L*x = y ; L = lower triang in  CSR format
  ! by Youcef Saad - Modified by GG January 2016
  !**************************************************************
  FUNCTION lsol(A, b) RESULT(x)
    TYPE(MAT_CSR_TYP), INTENT(IN)    :: A ! sparse matrix in CSR format
    REAL(8), DIMENSION(:), INTENT(IN) :: b ! RHS
    REAL(8), DIMENSION(size(b))       :: x ! unknowns
    INTEGER                           :: n, j, k
    INTEGER, DIMENSION(:), POINTER      :: cols, rowptr
    REAL(8)                           :: t, d
    REAL(8), DIMENSION(:), POINTER      :: vals

    n = A%n
    vals => A%vals
    cols => A%cols
    rowptr => A%rowptr

    x = 0.

    x(1) = b(1)/vals(1)
    DO k = 2, n
      t = b(k)
      DO j = rowptr(k), rowptr(k + 1) - 1
        t = t - vals(j)*x(cols(j))
        IF (k == cols(j)) d = vals(j)
      END DO
      x(k) = t/d
    END DO
  END FUNCTION lsol

  !**************************************************************
  ! USOL solves   U x = y    U = upper triangular in CSR format
  ! by Youcef Saad - Modified by GG January 2016
  !**************************************************************
  FUNCTION usol(A, b) RESULT(x)
    TYPE(MAT_CSR_TYP), INTENT(IN)    :: A ! sparse matrix in CSR format
    REAL(8), DIMENSION(:), INTENT(IN) :: b ! RHS
    REAL(8), DIMENSION(size(b))       :: x ! unknowns
    INTEGER                           :: n, j, k, nnz
    INTEGER, DIMENSION(:), POINTER      :: cols, rowptr
    REAL(8)                           :: t, d
    REAL(8), DIMENSION(:), POINTER      :: vals

    n = A%n
    nnz = A%nnz
    vals => A%vals
    cols => A%cols
    rowptr => A%rowptr

    x = 0.
    x(n) = b(n)/vals(nnz)
    DO k = n - 1, 1, -1
      t = b(k)
      DO j = rowptr(k), rowptr(k + 1) - 1
        t = t - vals(j)*x(cols(j))
        IF (k == cols(j)) d = vals(j)
      END DO
      x(k) = t/d
    END DO
  END FUNCTION usol

  !******************************************************
  ! Transpose of a sparse matrix in CSR format
  ! Inspired by Youcef Saad - GG January 2016
  !  A is n x m
  !  T is m x n
  !******************************************************
  SUBROUTINE trspspr(A, T)
    TYPE(MAT_CSR_TYP), INTENT(IN)    :: A
    TYPE(MAT_CSR_TYP), INTENT(OUT)   :: T

    INTEGER  :: i, n, k, j, next, m

    m = A%n               ! Number of columns of A
    n = size(A%rowptr) - 1  ! Number of rows of A
    T%n = n
    T%nnz = A%nnz
    ALLOCATE (T%rowptr(1:m + 1))
    ALLOCATE (T%cols(1:T%nnz))
    ALLOCATE (T%vals(1:T%nnz))

    !  Compute lengths of rows of A'.
    T%rowptr = 0
    DO i = 1, n
      DO k = A%rowptr(i), A%rowptr(i + 1) - 1
        j = A%cols(k) + 1
        T%rowptr(j) = T%rowptr(j) + 1
      END DO
    END DO

    !  Compute pointers from lengths.
    T%rowptr(1) = 1
    DO i = 1, m
      T%rowptr(i + 1) = T%rowptr(i) + T%rowptr(i + 1)
    END DO

    !  Do the actual copying.
    DO i = 1, n
      DO k = A%rowptr(i), A%rowptr(i + 1) - 1
        j = A%cols(k)
        next = T%rowptr(j)
        T%vals(next) = A%vals(k)
        T%cols(next) = i
        T%rowptr(j) = next + 1
      END DO
    END DO

    !  Reshift IAO and leave.
    DO i = m, 1, -1
      T%rowptr(i + 1) = T%rowptr(i)
    END DO
    T%rowptr(1) = 1

  END SUBROUTINE trspspr

  !**************************************************************
  ! AMUB multiply two sparse matrices in  CSR format
  ! Inspired by Youcef Saad - GG January 2016
  !**************************************************************
  ! A: n x m
  ! B: m x p
  ! C: n x p
  ! TODO: the output matrix is not correctly ordered, -->need to fix this!
  ! cols are not in crescent order: cannot use for ILU0!
  ! workaround: transpose and transpose back with trspspr
  SUBROUTINE amub(A, B, C)
    TYPE(MAT_CSR_TYP), INTENT(IN)  :: A, B ! INPUT: sparse matrix in CSR format
    TYPE(MAT_CSR_TYP), INTENT(OUT) :: C   ! OUTPUT:sparse matrix in CSR format
    INTEGER                          :: ncol, nrow, nnz, ii, jj, len, jpos, jcol, k, ka, kb
    INTEGER, DIMENSION(:), ALLOCATABLE :: iw
    REAL(8)                          :: scal

    ncol = B%n
    nrow = size(A%rowptr) - 1
    ALLOCATE (iw(1:ncol))

    ! Compute the nnz of the resulting matrix
    nnz = 0
    iw(1:ncol) = 0
    DO ii = 1, nrow
      DO ka = A%rowptr(ii), A%rowptr(ii + 1) - 1
        jj = A%cols(ka)
        DO kb = B%rowptr(jj), B%rowptr(jj + 1) - 1
          jcol = B%cols(kb)
          jpos = iw(jcol)
          IF (jpos == 0) THEN
            nnz = nnz + 1
            iw(jcol) = 1
          END IF
        END DO
      END DO
      iw = 0
    END DO

    ! Compute the actual result matrix
    C%n = B%n
    ALLOCATE (C%vals(1:nnz), C%cols(1:nnz))
    ALLOCATE (C%rowptr(1:nrow + 1))
    len = 0
    iw(1:ncol) = 0
    C%rowptr(1) = 1
    DO ii = 1, nrow
      DO ka = A%rowptr(ii), A%rowptr(ii + 1) - 1
        scal = A%vals(ka)
        jj = A%cols(ka)
        DO kb = B%rowptr(jj), B%rowptr(jj + 1) - 1
          jcol = B%cols(kb)
          jpos = iw(jcol)
          IF (jpos == 0) THEN
            len = len + 1
            C%cols(len) = jcol
            iw(jcol) = len
            C%vals(len) = scal*B%vals(kb)
          ELSE
            C%vals(jpos) = C%vals(jpos) + scal*B%vals(kb)
          END IF
        END DO
      END DO
      DO k = C%rowptr(ii), len
        iw(C%cols(k)) = 0
      END DO
      C%rowptr(ii + 1) = len + 1
    END DO
    C%nnz = nnz
  END SUBROUTINE amub

  !********************************
  ! Dump a CSR matrix to a text file
  ! in I,J,val format
  !********************************
  subroutine dump_CSR(prefix, rowptr, loc2glob, cols, vals)
    use MPI_OMP

    character(len=*) :: prefix
    integer, dimension(:) :: rowptr, loc2glob, cols
    real*8, dimension(:) :: vals

    character(len=100) :: tempstr, fname_complete
    integer :: Nrows, krow, kNZ

    Nrows = ubound(rowptr, 1) - 1

    write (tempstr, *) MPIvar%glob_id
    fname_complete = trim(adjustl(prefix))//'_p'//trim(adjustl(tempstr))//'.txt'
    open (UNIT=12, FILE=fname_complete, ACTION="write", STATUS="replace")
    !write(12,*) '%%MatrixMarket matrix coordinate real general'
    !write(12,*) maxval(loc2glob), maxval(loc2glob), rowptr(size(rowptr))-1
    do krow = 1, Nrows
      do kNZ = rowptr(krow), rowptr(krow + 1) - 1
        write (12, *) loc2glob(krow), cols(kNZ), vals(kNZ)
      enddo
    enddo
    close (12)

  end subroutine dump_CSR

  !********************************
  ! Dump a I,J,V matrix to a text file
  ! in I,J,val format
  !********************************
  subroutine dump_IJV(prefix, Irows, Jcols, Avals)
    use MPI_OMP

    character(len=*) :: prefix
    integer, dimension(:) :: Irows, Jcols
    real*8, dimension(:) :: Avals

    character(len=100) :: tempstr, fname_complete
    integer :: kNZ, NNZ

    NNZ = ubound(Avals, 1)

    write (tempstr, *) MPIvar%glob_id
    fname_complete = trim(adjustl(prefix))//'_p'//trim(adjustl(tempstr))//'.txt'
    open (UNIT=12, FILE=fname_complete, ACTION="write", STATUS="replace")
    do kNZ = 1, NNZ
      write (12, *) Irows(kNZ), Jcols(kNZ), Avals(kNZ)
    enddo
    close (12)

  end subroutine dump_IJV

  !********************************
  ! Dump a vector to a text file
  !********************************
  subroutine dump_vec(prefix, vec)
    use MPI_OMP

    character(len=*) :: prefix
    real*8, dimension(:) :: vec

    character(len=100) :: tempstr, fname_complete
    integer :: Nrows, krow

    Nrows = size(vec)

    write (tempstr, *) MPIvar%glob_id
    fname_complete = trim(adjustl(prefix))//'_p'//trim(adjustl(tempstr))//'.txt'
    open (UNIT=12, FILE=fname_complete, ACTION="write", STATUS="replace")
    do krow = 1, Nrows
      write (12, *) vec(krow)
    enddo
    close (12)

  end subroutine dump_vec

END MODULE matrices_tools
