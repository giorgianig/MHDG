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
      real*8 :: temp(size(A, 1), size(A, 2))

      external DGESV

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
! Compute eigenvalues and
! right eigenvectors
!***************************
   SUBROUTINE eig(A, V, D)
      real*8, intent(in)  :: A(:, :)
      real*8, intent(out) :: V(:, :)
      real*8, intent(out) :: D(:, :)
      integer            :: m, lwork, i, info
      integer, parameter  :: lwmax = 1000
      real               :: work(lwmax)
      real*8, allocatable :: Vmat(:, :), wr(:), wi(:)
      real*8             :: dummy(1, 1)

      external dgeev

      m = size(A, 1)
      lwork = -1

      ! Compute the eigenvalues (stored in the diagonal of D) and eigenvectors (stored in the columns of V) of A
      allocate (wr(1:m), wi(1:m))
      allocate (Vmat(1:m, 1:m))
      D = 0.d0
      Vmat = A

      ! Query the optimal workspace
      call dgeev('N', 'V', m, Vmat, m, wr, wi, V, m, dummy, 1, info)

      ! Solve eigenproblem.
      lwork = min(lwmax, int(work(1)))
      call dgeev('N', 'V', m, Vmat, m, wr, wi, V, m, dummy, 1, info)

      ! Only real eigenvalues are supposed
      do i = 1, m
         D(i, i) = wr(i)
      end do
      deallocate (wr, wi, Vmat)
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

