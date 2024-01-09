!*****************************************
! project: MHDG
! file: LinearAlgebra.f90
! date: 04/01/2016
! Module collecting different routines for
! linear algebra operations
!*****************************************

MODULE LinearAlgebra

CONTAINS

  SUBROUTINE compute_lup_decomp(A, n, L, U, P)
    integer*4, intent(in) :: n
    real*8, intent(in) :: A(:, :)
    real*8, intent(inout) :: L(:, :), U(:, :)
    integer*4, intent(inout) :: P(:, :)
    integer*4 :: i, j, itemp(1:n)
    real*8 :: temp(1:n, 1:n)

    temp = A
    CALL dgetrf(n, n, temp, n, itemp, i)
    L = temp
    U = temp
    DO i = 1, n
      L(i, i) = 1.d0
      DO j = 1, n
        IF (i .gt. j) THEN
          U(i, j) = 0.d0
        ELSEIF (i .lt. j) THEN
          L(i, j) = 0.d0
        END IF
      END DO
    END DO
    CALL find_permutation(n, P, itemp)
  END SUBROUTINE compute_lup_decomp

  SUBROUTINE find_permutation(n, P, per)
    integer*4, intent(in) :: n
    integer*4, intent(inout) :: per(1:n)
    integer*4, intent(out) :: P(1:n, 1:n)
    integer*4 :: i, temp, auxper(1:n)

    P = 0
    auxper = (/(i, i=1, n)/)

    DO i = 1, n
      temp = per(i)
      per(i) = auxper(per(i))
      auxper(temp) = auxper(i)
    END DO

    DO i = 1, n
      P(i, per(i)) = 1
    END DO
    !  CALL check_permutation(P,A,L,U,n)
  END SUBROUTINE find_permutation

  SUBROUTINE check_permutation(P, A, L, U, n)
    integer*4, intent(in) :: n, P(1:n, 1:n)
    real*8, intent(in) :: A(1:n, 1:n), L(1:n, 1:n), U(1:n, 1:n)
    real*8 :: temp1(1:n, 1:n), temp2(1:n, 1:n)
    integer*4 :: i, j

    temp1 = matmul(dble(P), A)
    temp2 = matmul(L, U)

    DO i = 1, n
      DO j = 1, n
        IF (temp1(i, j) - temp2(i, j) .gt. 1.d-10) THEN
          write (*, *) 'temp1-temp2 = ', temp1(i, j) - temp2(i, j), ' at i,j = ', i, j
          write (*, *) 'P'
          write (*, *) P(1, :)
          write (*, *) P(2, :)
          write (*, *) P(3, :)
          write (*, *) 'P*A'
          write (*, *) temp1(1, :)
          write (*, *) temp1(2, :)
          write (*, *) temp1(3, :)
          write (*, *) 'L*U'
          write (*, *) temp2(1, :)
          write (*, *) temp2(2, :)
          write (*, *) temp2(3, :)
          error stop "P*A != L*U"
        END IF
      END DO
    END DO

  END SUBROUTINE check_permutation

  !***************************
  ! Invert a square matrix
  !***************************
  SUBROUTINE invert_matrix(M, I)
    real*8, intent(in)   :: M(:, :)
    real*8, intent(out)  :: I(:, :)
    integer*4            :: n
    real*8               :: work(1:size(M, 1)) ! work array for LAPACK
    integer*4            :: ipiv(1:size(M, 1)) ! pivot indices
    integer*4            :: info

    external DGETRF
    external DGETRI

    n = size(M, 1)

    ! store the matrix in inverse to prevent it from being overwritten by LAPACK
    I = M

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL DGETRF(n, n, I, n, ipiv, info)
    IF (info .ne. 0) THEN
      stop 'Matrix is numerically singular!'
    END IF

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF
    CALL DGETRI(n, I, n, ipiv, work, n, info)
    IF (info .ne. 0) THEN
      stop 'Matrix inversion failed!'
    END IF

  END SUBROUTINE invert_matrix

  !***************************
  ! Solve a linear system with
  ! multiple rhs
  !***************************
  SUBROUTINE solve_linear_system(A, b, x)
    real*8, intent(in) :: A(:, :)
    real*8, intent(in) :: b(:, :)
    real*8, intent(out):: x(:, :)
    integer :: n, m, info
    integer :: p(1:size(b, 1))
    real*8,allocatable :: temp(:,:)

    external DGESV

    allocate(temp(size(A, 1), size(A, 2)))
    n = size(b, 1) ! size of the matrix (n x n)
    m = size(b, 2) ! number of RHS
    temp = A
    x = b
    CALL DGESV(n, m, temp, n, p, x, n, info)

    IF (INFO .GT. 0) THEN
      WRITE (*, *) 'The diagonal element of the triangular factor of A,'
      WRITE (*, *) 'U(', INFO, ',', INFO, ') is zero, so that'
      WRITE (*, *) 'A is singular; the solution could not be computed.'
      STOP
    END IF
    deallocate(temp)
  END SUBROUTINE solve_linear_system

  !***************************
  ! Solve a linear system with
  ! single rhs
  !***************************
  SUBROUTINE solve_linear_system_sing(A, b, x)
    real*8, intent(in) :: A(:, :)
    real*8, intent(in) :: b(:)
    real*8, intent(out):: x(:)
    integer :: n, m, info
    integer :: p(1:size(b, 1))
    real*8 :: temp(size(A, 1), size(A, 2))
    external DGESV

    n = size(b, 1) ! size of the matrix (n x n)
    m = 1 ! number of RHS
    temp = A
    x = b
    CALL DGESV(n, m, temp, n, p, x, n, info)

    IF (INFO .GT. 0) THEN
      WRITE (*, *) 'The diagonal element of the triangular factor of A,'
      WRITE (*, *) 'U(', INFO, ',', INFO, ') is zero, so that'
      WRITE (*, *) 'A is singular; the solution could not be computed.'
      STOP
    END IF
  END SUBROUTINE solve_linear_system_sing



  !***************************
  ! Matrix multiplication using
  ! dgemm (avoid stack overflow)
  !***************************
  SUBROUTINE mymatmul(A,B,C)
    real*8, intent(in) :: A(:, :)
    real*8, intent(in) :: B(:,:)
    real*8, intent(out):: C(:,:)
    integer :: a1,a2,b1,b2,c1,c2
    external DGEMM

    a1 = size(A,1)
    a2 = size(A,2)
    b1 = size(B,1)
    b2 = size(B,2)
    c1 = size(C,1)
    c2 = size(C,2)

    ! Little check
    if (a2.ne.b1) then
      write(6,*) "Error: number of colums of A different from number of rows of B"
      stop
    endif
    if (a1.ne.c1) then
      write(6,*) "Error: number of rows of A different from number of rows of C"
      stop
    endif
    if (b2.ne.c2) then
      write(6,*) "Error: number of colums of B different from number of columns of C"
      stop
    endif

    CALL DGEMM('N','N',a1,b2,a2,1.,A,a1,B,a2,0.,C,a1)


  END SUBROUTINE mymatmul



  !***************************
  ! Compute eigenvalues and
  ! right eigenvectors
  !***************************
  SUBROUTINE eig(A, V, D)
    real*8, intent(in)  :: A(:, :)
    real*8, intent(out) :: V(:, :)
    real*8, intent(out) :: D(:, :)
    integer            :: m, lwork, i, info
    Integer, Parameter :: nb = 64
    real*8,allocatable    :: work(:)
    real*8, allocatable :: Vmat(:, :), wr(:), wi(:)
    real*8             :: dummy(size(A,1), size(A,1))

    external dgeev

    m = size(A, 1)

    allocate (wr(1:m), wi(1:m))
    allocate (Vmat(1:m, 1:m))
    D = 0.d0
    Vmat = A

    ! Query the optimal workspace
    lwork = -1
    call dgeev('N', 'V', m, Vmat, m, wr, wi, dummy, 1, V, m, dummy, lwork, info)

    ! Compute the eigenvalues and right eigenvectors of A
    lwork = max((nb+2)*m, nint(dummy(1,1)))
    Allocate (work(lwork))
    call dgeev('N', 'V', m, Vmat, m, wr, wi, dummy, 1, V, m, work, lwork, info)

    ! Only real eigenvalues are supposed
    do i = 1, m
      D(i, i) = wr(i)
    end do
    deallocate (wr, wi, Vmat,work)
  END SUBROUTINE eig

  !*******************************************
  ! Tensor product
  !*******************************************
  FUNCTION TensorProduct(Vl, Vr) RESULT(Tens)
    REAL, DIMENSION(:) :: Vl, Vr
    REAL, DIMENSION(SIZE(Vl), SIZE(Vr)) :: Tens

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
      DO i = 1, Nl
        Tens(i, j) = Vl(i)*Vr(j)
      END DO
    END DO
  END FUNCTION TensorProduct

  !*******************************************
  ! Tensor product cumulative
  !*******************************************
  SUBROUTINE TensorProductCumul(M, Vl, Vr)
    REAL, DIMENSION(:, :) :: M
    REAL, DIMENSION(:)   :: Vl, Vr

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
      DO i = 1, Nl
        M(i, j) = M(i, j) + Vl(i)*Vr(j)
      END DO
    END DO
  END SUBROUTINE TensorProductCumul

  !*******************************************
  ! Tensor sum for integers
  !*******************************************
  FUNCTION TensorSumInt(Vl, Vr) RESULT(Tens)
    INTEGER, DIMENSION(:) :: Vl, Vr
    INTEGER, DIMENSION(SIZE(Vl), SIZE(Vr)) :: Tens

    INTEGER :: i, j, Nl, Nr

    Nl = SIZE(Vl); Nr = SIZE(Vr)

    DO j = 1, Nr
      DO i = 1, Nl
        Tens(i, j) = Vl(i) + Vr(j)
      END DO
    END DO
  END FUNCTION TensorSumInt

  !****************************************************
  ! Transform matrices in vectors for integers: colint
  !****************************************************
  FUNCTION colint(M) RESULT(V)
    INTEGER, DIMENSION(:, :) :: M
    INTEGER, DIMENSION(SIZE(M, 1)*SIZE(M, 2)) :: V

    INTEGER :: i, j

    DO j = 1, size(M, 2)
      DO i = 1, size(M, 1)
        V((j - 1)*size(M, 1) + i) = M(i, j)
      END DO
    END DO
  END FUNCTION colint

  !****************************************************
  ! Transform matrices in vectors for reals: col
  !****************************************************
  FUNCTION col(M) RESULT(V)
    REAL, DIMENSION(:, :) :: M
    REAL, DIMENSION(SIZE(M, 1)*SIZE(M, 2)) :: V

    INTEGER :: i, j

    DO j = 1, size(M, 2)
      DO i = 1, size(M, 1)
        V((j - 1)*size(M, 1) + i) = M(i, j)
      END DO
    END DO
  END FUNCTION col
END MODULE LinearAlgebra

SUBROUTINE cross_product(a, b, c)
  real*8, intent(in) :: a(3), b(3)
  real*8, intent(out):: c(3)

  c = 0.
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

END SUBROUTINE cross_product

SUBROUTINE ijk_cross_product(i, j, k)
  integer, intent(in) :: i
  integer, intent(out):: j, k

  SELECT CASE (i)
  CASE (1)
    j = 3
    k = 2
  CASE (2)
    j = 1
    k = 3
  CASE (3)
    j = 2
    k = 1
  CASE DEFAULT
    WRITE (*, *) "Error, wrong case in cross-product"
  END SELECT
END SUBROUTINE ijk_cross_product

