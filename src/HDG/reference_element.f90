!************************************************************
! project: MHDG
! file: reference_element.f90
! date: 16/12/2016
! Definition of the reference element
!************************************************************

MODULE reference_element
  USE globals
  USE LinearAlgebra
  USE printutils
  USE MPI_OMP

  IMPLICIT NONE

CONTAINS

  SUBROUTINE create_reference_element(refEl, ndim, degree, ngauss_1d)
    TYPE(Reference_element_type), intent(OUT) :: refEl
    integer, intent(in) :: ndim
    integer, intent(in), optional :: degree, ngauss_1d
    integer :: i

    !*********************************
    ! Initialize the reference element
    !*********************************
    refEl%ndim = ndim

    ! If 1d, I generate a refEl of type quad
    IF (ndim > 1) THEN
      refEl%elemType = Mesh%elemType
    ELSE
      refEl%elemType = 1
    ENDIF
    ! If the degree is passed as argument, I generate
    ! a refEl of that degree, otherwise the degree
    ! matches the one of the grid
    IF (PRESENT(degree)) THEN
      IF (refEl%elemType == 0) THEN
        refEl%Nnodes2D = 0.5*(degree + 1)*(degree + 2)
      ELSEIF (refEl%elemType == 1) THEN
        refEl%Nnodes2D = (degree + 1)**2
      ELSEIF (refEl%elemType == 2) THEN
        WRITE (6, *) "Thetrahedra not implemented yet"
        STOP
      ELSEIF (refEl%elemType == 3) THEN
        WRITE (6, *) "Hexahedra not implemented yet"
        STOP
      END IF
    ELSE
      IF (Mesh%ndim == 3) THEN
        refEl%Nnodes3D = Mesh%Nnodesperelem
      ELSEIF (Mesh%ndim == 2) THEN
        refEl%Nnodes2D = Mesh%Nnodesperelem
      ELSE
        WRITE (6, *) "Error! Wrong number of dimensions"
        STOP
      END IF
    END IF

    IF (PRESENT(ngauss_1d)) THEN
      CALL initialize_reference_element(refEl, ngauss_1d)
    ELSE
      CALL initialize_reference_element(refEl)
    END IF

    !***********************************************
    ! Compute the shape functions and the quadrature
    !***********************************************
    DO i = 1, ndim
      CALL compute_shape_functions(i, refEl)
    END DO

    !***********************************************
    ! Printing out stuff...
    !***********************************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*                REFERENCE ELEMENT              *'
        WRITE (6, *) '*************************************************'
        WRITE (6, '(A,I10)') " Number of dimensions              :    ", refEl%ndim
        WRITE (6, '(A,I10)') " Polynomial interpolation of degree:    ", refEl%nDeg
        IF (refEl%ndim == 3) THEN
          WRITE (6, '(A,I10)') " Number of shape functions 3D:          ", refEl%Nnodes3D
          WRITE (6, '(A,I10)') " Number of shape functions 2D:          ", refEl%Nnodes2D
          WRITE (6, '(A,I10)') " Number of shape functions 1D:          ", refEl%Nnodes1D
        ELSE IF (refEl%ndim == 2) THEN
          WRITE (6, '(A,I10)') " Number of shape functions 2D:          ", refEl%Nnodes2D
          WRITE (6, '(A,I10)') " Number of shape functions 1D:          ", refEl%Nnodes1D
        ELSE
          WRITE (6, '(A,I10)') " Number of shape functions 1D:          ", refEl%Nnodes1D
        END IF

        IF (refEl%ndim == 3) THEN
          WRITE (6, '(A,I10)') " Number of Gauss points 3D:             ", refEl%NGauss3D
          WRITE (6, '(A,I10)') " Number of Gauss points 2D:             ", refEl%NGauss2D
          WRITE (6, '(A,I10)') " Number of Gauss points 1D:             ", refEl%NGauss1D
        ELSE IF (refEl%ndim == 2) THEN
          WRITE (6, '(A,I10)') " Number of Gauss points 2D:             ", refEl%NGauss2D
          WRITE (6, '(A,I10)') " Number of Gauss points 1D:             ", refEl%NGauss1D
        ELSE
          WRITE (6, '(A,I10)') " Number of Gauss points 1D:             ", refEl%NGauss1D
        END IF

        print *, '        '
      END IF

      IF (utils%printint > 1) THEN
        IF (Mesh%ndim == 3) THEN
          WRITE (6, *) "Nodes coordinates 3D:"
          CALL displayMatrix(refEl%coord3d)
        END IF
        WRITE (6, *) "Nodes coordinates 2D:"
        CALL displayMatrix(refEl%coord2d)
        WRITE (6, *) "Nodes coordinates 1D:"
        CALL displayVector(refEl%coord1d)
        IF (Mesh%ndim == 3) THEN
          WRITE (6, *) "Gauss points 3D coordinates:"
          CALL displayMatrix(refEl%gauss_points3D)
        END IF
        WRITE (6, *) "Gauss points 2D coordinates:"
        CALL displayMatrix(refEl%gauss_points2D)
        WRITE (6, *) "Gauss points 1D coordinates:"
        CALL displayVector(refEl%gauss_points1D)
        WRITE (6, *) "Face nodes:"
        CALL displayMatrixInt(refEl%face_nodes)
        IF (Mesh%ndim == 3) THEN
          WRITE (6, *) "Shape functions 3D at Gauss points:"
          CALL displayMatrix(refEl%N3D)
        END IF
        WRITE (6, *) "Shape functions 2D at Gauss points:"
        CALL displayMatrix(refEl%N2D)
        WRITE (6, *) "Shape functions 1D at Gauss points:"
        CALL displayMatrix(refEl%N1D)
      END IF
    ENDIF
  END SUBROUTINE create_reference_element

  !****************************************************************************************
  !
  !                        INITIALIZATION OF THE REFERENCE ELEMENT
  !
  !****************************************************************************************
  SUBROUTINE initialize_reference_element(refEl, ngauss_1d)

    TYPE(Reference_element_type), intent(INOUT) :: refEl
    integer*4, intent(in), optional :: ngauss_1d
    integer*4 :: i, j
    integer, allocatable :: aux_coord2d(:, :), perm(:), ind(:)

    IF (refEl%elemType .eq. 0) THEN
      !*************************************
      ! Triangular elements
      !*************************************
      refEl%Nvertices = 3
      refEl%Nfacenodeslin = 2
      refEl%Nfaces = 3
      ALLOCATE (refEl%coord2d(1:refEl%Nnodes2D, 1:2))
      SELECT CASE (refEl%Nnodes2D)
      CASE (3)
        refEl%Ndeg = 1
        refEl%coord2d(1, :) = (/-1.d0, -1.d0/)
        refEl%coord2d(2, :) = (/1.d0, -1.d0/)
        refEl%coord2d(3, :) = (/-1.d0, 1.d0/)
      CASE (6)
        refEl%Ndeg = 2
        refEl%coord2d = reshape((/-1.d0, 1.d0, -1.d0, 0.d0, 0.d0, -1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 0.d0, 0.d0/), &
          (/refEl%Nnodes2D, 2/))
      CASE (10)
        refEl%Ndeg = 3
      CASE (15)
        refEl%Ndeg = 4
      CASE (21)
        refEl%Ndeg = 5
      CASE (28)
        refEl%Ndeg = 6
      CASE (36)
        refEl%Ndeg = 7
      CASE (45)
        refEl%Ndeg = 8
      CASE (55)
        refEl%Ndeg = 9
      CASE (66)
        refEl%Ndeg = 10
      CASE (78)
        refEl%Ndeg = 11
      CASE (91)
        refEl%Ndeg = 12
      CASE (105)
        refEl%Ndeg = 13
      CASE (120)
        refEl%Ndeg = 14
      CASE (136)
        refEl%Ndeg = 15
      CASE (153)
        refEl%Ndeg = 16
      CASE (171)
        refEl%Ndeg = 17
      CASE (190)
        refEl%Ndeg = 18
      CASE (210)
        refEl%Ndeg = 19
      CASE (231)
        refEl%Ndeg = 20
      END SELECT
      refEl%Nextnodes = refEl%Nvertices + refEl%Nfaces*(refEl%Ndeg - 1)
      refEl%Nfacenodes = refEl%Ndeg + 1
      refEl%Nnodes1D = refEl%Ndeg + 1
      refEl%Ninnernodes = refEl%Nnodes2D-refEl%Nextnodes

      IF (refEl%Ndeg .gt. 2) THEN ! EZ4U rules imply fekete nodes
        CALL read_fekete_nodes_2D(refEl)
      END IF
      CALL feketeNodes1D(refEl)

    ELSEIF (refEl%elemType .eq. 1) THEN

      !*************************************
      ! Quadrilateral elements
      !*************************************
      refEl%Nvertices = 4
      refEl%Nfaces = 4
      refEl%Nfacenodeslin = 2

      ALLOCATE (refEl%coord2d(1:refEl%Nnodes2D, 1:2))
      ALLOCATE (perm(1:refEl%Nnodes2D), aux_coord2d(1:refEl%Nnodes2D, 1:2))
      SELECT CASE (refEl%Nnodes2D)
      CASE (4)
        refEl%Ndeg = 1
        perm = (/1, 2, 4, 3/)
      CASE (9)
        refEl%Ndeg = 2
        perm = (/1, 5, 2, 8, 9, 6, 4, 7, 3/)
      CASE (16)
        refEl%Ndeg = 3
        perm = (/1, 5, 6, 2, 12, 13, 14, 7, 11, 15, 16, 8, 4, 10, 9, 3/)
      CASE (25)
        refEl%Ndeg = 4
        perm = (/1, 5, 6, 7, 2, 16, 17, 18, 19, 8, 15, 20, 21, 22, 9, 14,&
          &23, 24, 25, 10, 4, 13, 12, 11, 3/)
      CASE (36)
        refEl%Ndeg = 5
        perm = (/1, 5, 6, 7, 8, 2, 20, 21, 22, 23, 24, 9, 19, 25, 26, 27, 28,&
          &10, 18, 29, 30, 31, 32, 11, 17, 33, 34, 35, 36, 12, 4, 16, 15, 14, 13, 3/)
      CASE (49)
        refEl%Ndeg = 6
        perm = (/1, 5, 6, 7, 8, 9, 2, 24, 25, 26, 27, 28, 29, 10, 23, 30, 31,&
          &32, 33, 34, 11, 22, 35, 36, 37, 38, 39, 12, 21, 40, 41, 42,&
          &43, 44, 13, 20, 45, 46, 47, 48, 49, 14, 4, 19, 18, 17, 16, 15, 3/)
      CASE (64)
        refEl%Ndeg = 7
        perm = (/1, 5, 6, 7, 8, 9, 10, 2, 28, 29, 30, 31, 32, 33, 34, 11, 27, 35, 36,&
          &37, 38, 39, 40, 12, 26, 41, 42, 43, 44, 45, 46, 13, 25, 47, 48,&
          &49, 50, 51, 52, 14, 24, 53, 54, 55, 56, 57, 58, 15, 23, 59, 60,&
          &61, 62, 63, 64, 16, 4, 22, 21, 20, 19, 18, 17, 3/)
      CASE (81)
        refEl%Ndeg = 8
        perm = (/1, 5, 6, 7, 8, 9, 10, 11, 2, 32, 33, 34, 35, 36, 37, 38, 39, 12,&
          &31, 40, 41, 42, 43, 44, 45, 46, 13, 30, 47, 48, 49,&
          &50, 51, 52, 53, 14, 29, 54, 55, 56, 57, 58, 59, 60, 15,&
          &28, 61, 62, 63, 64, 65, 66, 67, 16, 27, 68, 69, 70, 71,&
          &72, 73, 74, 17, 26, 75, 76, 77, 78, 79, 80, 81, 18, 4,&
          &25, 24, 23, 22, 21, 20, 19, 3/)
      CASE DEFAULT
        WRITE (6, *) "Error, quadrilateral element not available"
        stop
      END SELECT

      refEl%Nextnodes = refEl%Nvertices + refEl%Nfaces*(refEl%Ndeg - 1)
      refEl%Nfacenodes = refEl%Ndeg + 1
      refEl%Nnodes1D = refEl%Ndeg + 1
      refEl%Ninnernodes = refEl%Nnodes2D-refEl%Nextnodes

      CALL feketeNodes1D(refEl) ! fill refEl%coord1d

      ALLOCATE (ind(refEl%Ndeg + 1))
      DO i = 1, refEl%Ndeg + 1
        ind = (i - 1)*(refEl%Ndeg + 1) + (/(j, j=1, refEl%Ndeg + 1)/)
        aux_coord2d(ind, 1) = refEl%coord1d
        aux_coord2d(ind, 2) = refEl%coord1d(i)
      END DO
      DEALLOCATE (ind)

      DO i = 1, refEl%Nnodes2D
        refEl%coord2d(perm(i), :) = aux_coord2d(i, :)
      END DO
      DEALLOCATE (aux_coord2d, perm)
    ELSEIF (refEl%elemType .eq. 2) THEN
      refEl%Nvertices = 4
      refEl%Nfaces = 4
      refEl%Nfacenodeslin = 3
      !*************************************
      ! Thetrahedra elements
      !*************************************
    ELSEIF (refEl%elemType .eq. 3) THEN
      refEl%Nvertices = 8
      refEl%Nfaces = 6
      refEl%Nfacenodeslin = 4
      !*************************************
      ! Hexahedra elements
      !*************************************
    ELSE
      WRITE (6, *) 'Error! Element type not valid'
      STOP
    END IF

    !****************************************************************
    ! Create face_nodes
    ALLOCATE (refEl%Face_nodes(1:refEl%Nfaces, 1:refEl%Nfacenodes))
    refEl%Face_nodes = 0
    IF (refEl%Ninnernodes .gt. 0) THEN
      ALLOCATE (refEl%inner_nodes(1:refEl%Ninnernodes))
      DO i = 1, refEl%Ninnernodes
        refEl%inner_nodes(i) = refEl%Nextnodes + i
      END DO
    END IF

    DO i = 1, refEl%Nfaces
      refEl%Face_nodes(i, 1) = i
      refEl%Face_nodes(i, refEl%Nfacenodes) = i + 1
      IF (i .eq. refEl%Nfaces) THEN
        refEl%Face_nodes(i, refEl%Nfacenodes) = 1
      END IF
      IF (refEl%Nfacenodes .gt. 2) THEN
        IF (i .eq. 1) THEN
          DO j = 2, refEl%Nfacenodes - 1
            !              refEl%Face_nodes(i,j) = refEl%Face_nodes(refEl%Nfaces,1)+j-1
            refEl%Face_nodes(i, j) = refEl%Nfaces + j - 1
          END DO
        ELSE
          DO j = 2, refEl%Nfacenodes - 1
            refEl%Face_nodes(i, j) = refEl%Face_nodes(i - 1, j) + refEl%Nfacenodes - 2
          END DO
        END IF
      END IF
    END DO
    !****************************************************************

    !    refEl%NGauss1D = refEl%Ndeg + 2

    IF (PRESENT(ngauss_1d)) THEN
      refEl%NGauss1D = ngauss_1d
    ELSE
      refEl%NGauss1D = refEl%Ndeg + 2
    END IF

  END SUBROUTINE initialize_reference_element

  !****************************************************************************************
  !
  !                        DEFINITION OF THE SHAPE FUNCTIONS
  !
  !****************************************************************************************

  SUBROUTINE compute_shape_functions(dimensions, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    integer*4, intent(in) :: dimensions
    integer*4 :: eltype
    eltype = refEl%elemType

    SELECT CASE (dimensions)

    CASE (3)
      ! 3D shape functions
      IF (eltype == 3) THEN
        ! Hexahedra
        !TODO code the following routine
        !           CALL compute_shape_functions_hexa()
      ELSEIF (eltype == 2) THEN
        ! Thetrahedra
        CALL compute_shape_functions_thetra(refEl)
      ELSE
        WRITE (6, *) "Error! Wrong element type"
        STOP
      END IF

    CASE (2)
      ! 2D shape functions
      IF (eltype == 1) THEN
        ! Quadrilaterals
        CALL compute_shape_functions_quads(refEl)
      ELSEIF (eltype == 0) THEN
        ! Triangles
        CALL compute_shape_functions_triangles(refEl)
      ELSE
        WRITE (6, *) "Error! Wrong element type"
        STOP
      END IF

    CASE (1)
      ! 1D shape functions
      CALL compute_shape_functions_1d(refEl)

    CASE DEFAULT
      WRITE (6, *) "Error! Wrong number of dimensions in shape functions creation"
      STOP
    END SELECT

  END SUBROUTINE compute_shape_functions

  !**********************
  ! Shape functions in 1D
  !**********************
  SUBROUTINE compute_shape_functions_1d(refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    real*8, allocatable :: V(:, :), L(:, :), U(:, :), p(:), dp(:)
    real*8, allocatable :: Linv(:, :), Uinv(:, :), A(:, :), temp(:, :)
    integer*4, allocatable :: Pmatrix(:, :)
    integer*4 :: i, j

    ALLOCATE (refEl%gauss_points1D(1:refEl%NGauss1D))
    ALLOCATE (refEl%gauss_weights1D(1:refEl%NGauss1D))
    ALLOCATE (refEl%N1D(1:refEl%NGauss1D, 1:refEl%Nfacenodes))
    ALLOCATE (refEl%Nxi1D(1:refEl%NGauss1D, 1:refEl%Nfacenodes))

    CALL gauss_legendre(refEl%NGauss1D, refEl%gauss_points1D, refEl%gauss_weights1D)
    !!!!!!!!!!!!!!!!!!!!!! TEMPORARY FIX
    refEl%gauss_points1D = -refEl%gauss_points1D
    !  CALL gauss_legendre_new((refEl%NGauss1D,-1.d0,1.d0,refEl%gauss_points1D,refEl%gauss_weights1D)

    ALLOCATE (V(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))

    CALL vandermonde_1d(V, refEl)

    ALLOCATE (L(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))
    ALLOCATE (U(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))
    ALLOCATE (Pmatrix(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))
    L = 0.d0
    U = 0.d0
    Pmatrix = 0
    CALL compute_lup_decomp(transpose(V), refEl%Nfacenodes, L, U, Pmatrix)

    ALLOCATE (Linv(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))
    CALL invert_matrix(L, Linv)
    ALLOCATE (Uinv(1:refEl%Nfacenodes, 1:refEl%Nfacenodes))
    CALL invert_matrix(U, Uinv)

    ALLOCATE (p(1:refEl%Nfacenodes), dp(1:refEl%Nfacenodes))
    ALLOCATE (A(1:refEl%Nfacenodes, 1:2))
    ALLOCATE (temp(1:refEl%Nfacenodes, 1:2))

    DO i = 1, refEl%NGauss1D
      p = 0.d0
      dp = 0.d0
      A = 0.d0
      temp = 0.d0
      CALL orthopoly1D_deriv(refEl%gauss_points1D(i), refEl%Ndeg, p, dp)
      DO j = 1, refEl%Nfacenodes
        A(j, 1) = p(j)
        A(j, 2) = dp(j)
      END DO
      temp = matmul(Uinv, matmul(Linv, matmul(Pmatrix, A)))
      refEl%N1D(i, :) = temp(:, 1)
      refEl%Nxi1D(i, :) = temp(:, 2)
    END DO

    DEALLOCATE (V, L, U, p, dp)
    DEALLOCATE (Linv, Uinv, A, temp)
    DEALLOCATE (Pmatrix)

  END SUBROUTINE compute_shape_functions_1d

  !*********************************
  ! Shape functions 2D for triangles
  !*********************************
  SUBROUTINE compute_shape_functions_triangles(refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    integer*4 :: order_cub, n, i, j, k
    real*8, allocatable :: V(:, :), L(:, :), U(:, :), Linv(:, :), Uinv(:, :)
    real*8, allocatable :: p(:), dpxi(:), dpeta(:), A(:, :), temp(:, :)
    real*8, allocatable :: x(:), w(:)
    integer*4, allocatable :: Pmatrix(:, :)

    IF ((refEl%Ndeg .lt. 12) .and. (refEl%NGauss1D-2 .lt. 11)) THEN
      IF (refEl%Ndeg + 2 .eq. refEl%NGauss1D) THEN
        SELECT CASE (refEl%Ndeg)
        CASE (1)
          n = 5
        CASE (2, 3)
          n = 10
        CASE (4, 5, 6, 7)
          n = 15
        CASE (8, 9, 10, 11)
          n = 25
        END SELECT
      ELSE
        SELECT CASE (refEl%NGauss1D-2)
        CASE (1)
          n = 5
        CASE (2, 3)
          n = 10
        CASE (4, 5, 6, 7)
          n = 15
        CASE (8, 9, 10, 11)
          n = 25
        END SELECT
      END IF

      !     ALLOCATE(refEl%gauss_points2D(1:n,1:2),refEl%gauss_weights2D(1:n))
      !     CALL gauss_legendre_cubature2d(n,refEl%gauss_points2D,refEl%gauss_weights2D)
      CALL gauss_legendre_cubature2d_new(n, refEl)

      refEl%NGauss2D = size(refEl%gauss_points2D, 1)

      DO i = 1, refEl%NGauss2D
        refEl%gauss_points2D(i, :) = 2.d0*refEl%gauss_points2D(i, :) - 1.d0
        refEl%gauss_weights2D(i) = 2.d0*refEl%gauss_weights2D(i)
      END DO

      ALLOCATE (V(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      CALL vandermonde_2d(V, refEl)

      ALLOCATE (L(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      ALLOCATE (U(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      ALLOCATE (Pmatrix(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      L = 0.d0
      U = 0.d0
      Pmatrix = 0
      CALL compute_lup_decomp(transpose(V), refEl%Nnodes2D, L, U, Pmatrix)

      ALLOCATE (Linv(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      Linv = 0.d0
      CALL invert_matrix(L, Linv)
      ALLOCATE (Uinv(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      CALL invert_matrix(U, Uinv)

      ALLOCATE (refEl%N2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
      ALLOCATE (refEl%Nxi2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
      ALLOCATE (refEl%Neta2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))

      ALLOCATE (p(1:refEl%Nnodes2D), dpxi(1:refEl%Nnodes2D), dpeta(1:refEl%Nnodes2D))
      ALLOCATE (A(1:refEl%Nnodes2D, 1:3))
      ALLOCATE (temp(1:refEl%Nnodes2D, 1:3))
      DO i = 1, refEl%NGauss2D
        p = 0.d0
        dpxi = 0.d0
        dpeta = 0.d0
        A = 0.d0
        temp = 0.d0
        CALL orthopoly2d_deriv(refEl%gauss_points2D(i, 1), refEl%gauss_points2D(i, 2), &
          refEl%Ndeg, refEl%Nnodes2D, p, dpxi, dpeta)
        A(:, 1) = p(:)
        A(:, 2) = dpxi(:)
        A(:, 3) = dpeta(:)
        temp = matmul(Uinv, matmul(Linv, matmul(Pmatrix, A)))
        refEl%N2D(i, :) = temp(:, 1)
        refEl%Nxi2D(i, :) = temp(:, 2)
        refEl%Neta2D(i, :) = temp(:, 3)
      END DO
    ELSE
      ALLOCATE (x(1:refEl%NGauss1D), w(1:refEl%NGauss1D))
      CALL gauss_legendre(refEl%NGauss1D, x, w)
      refEl%NGauss2D = refEl%NGauss1D**2

      ALLOCATE (V(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      CALL vandermonde_2d(V, refEl)

      ALLOCATE (L(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      ALLOCATE (U(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      ALLOCATE (Pmatrix(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      L = 0.d0
      U = 0.d0
      Pmatrix = 0
      CALL compute_lup_decomp(transpose(V), refEl%Nnodes2D, L, U, Pmatrix)

      ALLOCATE (Linv(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      CALL invert_matrix(L, Linv)
      ALLOCATE (Uinv(1:refEl%Nnodes2D, 1:refEl%Nnodes2D))
      CALL invert_matrix(U, Uinv)

      ALLOCATE (refEl%gauss_points2D(1:refEl%NGauss2D, 1:2))
      ALLOCATE (refEl%gauss_weights2D(1:refEl%NGauss2D))
      ALLOCATE (refEl%N2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
      ALLOCATE (refEl%Nxi2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
      ALLOCATE (refEl%Neta2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))

      ALLOCATE (p(1:refEl%Nnodes2D), dpxi(1:refEl%Nnodes2D), dpeta(1:refEl%Nnodes2D))
      ALLOCATE (A(1:refEl%Nnodes2D, 1:3))
      ALLOCATE (temp(1:refEl%Nnodes2D, 1:3))
      k = 1
      DO i = 1, refEl%NGauss1D
        DO j = 1, refEl%NGauss1D
          refEl%gauss_points2D(k, 1) = (1.d0 + x(i))*(1.d0 - x(j))/2.d0 - 1.d0
          refEl%gauss_points2D(k, 2) = x(j)
          refEl%gauss_weights2D(k) = w(i)*w(j)*(1.d0 - x(j))/2.d0
          p = 0.d0
          dpxi = 0.d0
          dpeta = 0.d0
          A = 0.d0
          temp = 0.d0
          CALL orthopoly2d_deriv(x(i), x(j), refEl%Ndeg, refEl%Nnodes2D, p, dpxi, dpeta)
          A(:, 1) = p(:)
          A(:, 2) = dpxi(:)
          A(:, 3) = dpeta(:)
          temp = matmul(Uinv, matmul(Linv, matmul(Pmatrix, A)))
          refEl%N2D(k, :) = temp(:, 1)
          refEl%Nxi2D(k, :) = temp(:, 2)
          refEl%Neta2D(k, :) = temp(:, 3)
          k = k + 1
        END DO
      END DO
      DEALLOCATE (x, w)
    END IF

    DEALLOCATE (V, L, U, Linv, Uinv)
    DEALLOCATE (p, dpxi, dpeta, A, temp)
    DEALLOCATE (Pmatrix)

  END SUBROUTINE compute_shape_functions_triangles

  !*********************************
  ! Shape functions 2D for quads
  !*********************************
  SUBROUTINE compute_shape_functions_quads(refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    real*8 :: x(1:refEl%NGauss1D), w(1:refEl%NGauss1D)
    integer*4 :: i, j, k, n
    integer, allocatable :: perm(:)

    CALL gauss_legendre(refEl%NGauss1D, x, w)
    refEl%NGauss2D = refEl%NGauss1D**2
    ALLOCATE (refEl%gauss_points2D(1:refEl%NGauss2D, 1:2), refEl%gauss_weights2D(1:refEl%NGauss2D))
    k = 1
    DO i = 1, refEl%NGauss1D
      DO j = 1, refEl%NGauss1D
        refEl%gauss_points2D(k, 1) = x(i)
        refEl%gauss_points2D(k, 2) = x(j)
        refEl%gauss_weights2D(k) = w(i)*w(j)
        k = k + 1
      END DO
    END DO

    ALLOCATE (refEl%N2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
    ALLOCATE (refEl%Nxi2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))
    ALLOCATE (refEl%Neta2D(1:refEl%NGauss2D, 1:refEl%Nnodes2D))

    ALLOCATE (perm(1:refEl%Nnodes2D))
    SELECT CASE (refEl%Nnodes2D)
    CASE (4)
      perm = (/1, 2, 4, 3/)
    CASE (9)
      perm = (/1, 5, 2, 8, 9, 6, 4, 7, 3/)
    CASE (16)
      perm = (/1, 5, 6, 2, 12, 13, 14, 7, 11, 15, 16, 8, 4, 10, 9, 3/)
    CASE (25)
      perm = (/1, 5, 6, 7, 2, 16, 17, 18, 19, 8, 15, 20, 21, 22, 9, 14, 23, 24,&
        &25, 10, 4, 13, 12, 11, 3/)
    CASE (36)
      perm = (/1, 5, 6, 7, 8, 2, 20, 21, 22, 23, 24, 9, 19, 25, 26, 27, 28, 10, 18,&
        &29, 30, 31, 32, 11, 17, 33, 34, 35, 36, 12, 4, 16, 15, 14, 13, 3/)
    CASE (49)
      perm = (/1, 5, 6, 7, 8, 9, 2, 24, 25, 26, 27, 28, 29, 10, 23, 30, 31, 32, 33,&
        &34, 11, 22, 35, 36, 37, 38, 39, 12, 21, 40, 41, 42, 43, 44, 13,&
        &20, 45, 46, 47, 48, 49, 14, 4, 19, 18, 17, 16, 15, 3/)
    CASE (64)
      perm = (/1, 5, 6, 7, 8, 9, 10, 2, 28, 29, 30, 31, 32, 33, 34, 11, 27, 35, 36,&
        &37, 38, 39, 40, 12, 26, 41, 42, 43, 44, 45, 46, 13, 25, 47, 48,&
        &49, 50, 51, 52, 14, 24, 53, 54, 55, 56, 57, 58, 15, 23, 59, 60,&
        &61, 62, 63, 64, 16, 4, 22, 21, 20, 19, 18, 17, 3/)
    CASE (81)
      perm = (/1, 5, 6, 7, 8, 9, 10, 11, 2, 32, 33, 34, 35, 36, 37, 38, 39, 12, 31, 40, 41, 42,&
        &43, 44, 45, 46, 13, 30, 47, 48, 49, 50, 51, 52, 53, 14, 29,&
        &54, 55, 56, 57, 58, 59, 60, 15, 28, 61, 62, 63, 64, 65, 66,&
        &67, 16, 27, 68, 69, 70, 71, 72, 73, 74, 17, 26, 75, 76, 77,&
        &78, 79, 80, 81, 18, 4, 25, 24, 23, 22, 21, 20, 19, 3/)
    CASE DEFAULT
      WRITE (6, *) "Error, quadrilateral element not available"
      stop
    END SELECT
    k = 1
    DO i = 1, refEl%Nnodes1D
      DO j = 1, refEl%Nnodes1D
        refEl%N2D(:, perm(k)) = reshape(tensorProduct(refEl%N1d(:, j), refEl%N1d(:, i)), (/refEl%NGauss2D/))
        refEl%Nxi2D(:, perm(k)) = reshape(tensorProduct(refEl%Nxi1d(:, j), refEl%N1d(:, i)), (/refEl%NGauss2D/))
        refEl%Neta2D(:, perm(k)) = reshape(tensorProduct(refEl%N1d(:, j), refEl%Nxi1d(:, i)), (/refEl%NGauss2D/))
        k = k + 1
      END DO
    END DO
    DEALLOCATE (perm)

  END SUBROUTINE compute_shape_functions_quads

  FUNCTION permutation_quads(refEl) RESULT(perm)
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    integer :: perm(refEl%Nnodes2D)
    SELECT CASE (refEl%Nnodes2D)
    CASE (4)
      perm = (/1, 2, 4, 3/)
    CASE (9)
      perm = (/1, 5, 2, 8, 9, 6, 4, 7, 3/)
    CASE (16)
      perm = (/1, 5, 6, 2, 12, 13, 14, 7, 11, 15, 16, 8, 4, 10, 9, 3/)
    CASE (25)
      perm = (/1, 5, 6, 7, 2, 16, 17, 18, 19, 8, 15, 20, 21, 22, 9, 14, 23, 24, 25, 10, 4, 13, 12, 11, 3/)
    CASE (36)
      perm = (/1, 5, 6, 7, 8, 2, 20, 21, 22, 23, 24, 9, 19, 25, 26, 27, 28, 10, 18, 29, 30, 31, 32, 11, 17,&
        &33, 34, 35, 36, 12, 4, 16, 15, 14, 13, 3/)
    CASE (49)
      perm = (/1, 5, 6, 7, 8, 9, 2, 24, 25, 26, 27, 28, 29, 10, 23, 30, 31, 32, 33, 34, 11, 22, 35, 36, 37, 38,&
        &39, 12, 21, 40, 41, 42, 43, 44, 13, 20, 45, 46, 47, 48, 49, 14, 4, 19, 18, 17, 16, 15, 3/)
    CASE (64)
      perm = (/1, 5, 6, 7, 8, 9, 10, 2, 28, 29, 30, 31, 32, 33, 34, 11, 27, 35, 36, 37, 38, 39, 40, 12, 26, 41,&
        &42, 43, 44, 45, 46, 13, 25, 47, 48, 49, 50, 51, 52, 14, 24, 53, 54, 55, 56, 57, 58, 15, 23,&
        &59, 60, 61, 62, 63, 64, 16, 4, 22, 21, 20, 19, 18, 17, 3/)
    CASE (81)
      perm = (/1, 5, 6, 7, 8, 9, 10, 11, 2, 32, 33, 34, 35, 36, 37, 38, 39, 12, 31, 40, 41, 42, 43, 44, 45, 46,&
        &13, 30, 47, 48, 49, 50, 51, 52, 53, 14, 29, 54, 55, 56, 57, 58, 59, 60, 15, 28, 61,&
        &62, 63, 64, 65, 66, 67, 16, 27, 68, 69, 70, 71, 72, 73, 74, 17, 26, 75, 76, 77, 78,&
        &79, 80, 81, 18, 4, 25, 24, 23, 22, 21, 20, 19, 3/)
    CASE DEFAULT
      WRITE (6, *) "Error, quadrilateral element not available"
      stop
    END SELECT
  END FUNCTION permutation_quads

  !**********************************
  ! Shape functions 3D for tetrahedra
  !**********************************
  SUBROUTINE compute_shape_functions_thetra(refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    real*8, allocatable :: V(:, :), L(:, :), U(:, :), Linv(:, :), Uinv(:, :)
    real*8, allocatable :: p(:), dpxi(:), dpeta(:), dpzeta(:)
    real*8, allocatable :: A(:, :), temp(:, :)
    integer*4, allocatable :: Pmatrix(:, :)
    real*8 :: x(1:refEl%NGauss1D), w(1:refEl%NGauss1D)
    real*8 :: r, s, t, xi, eta
    integer*4 :: i, j, k, il, n

    CALL gauss_legendre(refEl%NGauss1D, x, w)
    ALLOCATE (V(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    CALL vandermonde_3d(V, refEl)
    ALLOCATE (L(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    ALLOCATE (U(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    ALLOCATE (Pmatrix(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    L = 0.d0
    U = 0.d0
    Pmatrix = 0
    CALL compute_lup_decomp(transpose(V), refEl%Nnodes3D, L, U, Pmatrix)
    ALLOCATE (Linv(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    CALL invert_matrix(L, Linv)
    ALLOCATE (Uinv(1:refEl%Nnodes3D, 1:refEl%Nnodes3D))
    CALL invert_matrix(U, Uinv)

    refEl%NGauss3D = refEl%NGauss1D**3
    ALLOCATE (refEl%N3D(1:refEl%NGauss3D, 1:refEl%Nnodes3D))
    ALLOCATE (refEl%Nxi3D(1:refEl%NGauss3D, 1:refEl%Nnodes3D))
    ALLOCATE (refEl%Neta3D(1:refEl%NGauss3D, 1:refEl%Nnodes3D))
    ALLOCATE (refEl%Nzeta3D(1:refEl%NGauss3D, 1:refEl%Nnodes3D))

    ALLOCATE (refEl%gauss_points3D(1:refEl%NGauss3D, 1:3), refEl%gauss_weights3D(1:refEl%NGauss3D))

    ALLOCATE (p(1:refEl%Nnodes3D), dpxi(1:refEl%Nnodes3D))
    ALLOCATE (dpeta(1:refEl%Nnodes3D), dpzeta(1:refEl%Nnodes3D))
    ALLOCATE (A(1:refEl%Nnodes3D, 1:4), temp(1:refEl%Nnodes3D, 1:4))

    il = 1
    DO i = 1, refEl%NGauss1D
      DO j = 1, refEl%NGauss1D
        DO k = 1, refEl%NGauss1D
          p = 0.d0
          dpxi = 0.d0
          dpeta = 0.d0
          dpzeta = 0.d0
          A = 0.d0
          temp = 0.d0
          CALL orthopoly3d_deriv(x(i), x(j), x(k), refEl%Ndeg, refEl%Nnodes3D, &
            p, dpxi, dpeta, dpzeta)
          A(:, 1) = p(:)
          A(:, 2) = dpxi(:)
          A(:, 3) = dpeta(:)
          A(:, 4) = dpzeta(:)
          temp = matmul(Uinv, matmul(Linv, matmul(Pmatrix, A)))
          refEl%N3D(il, :) = temp(:, 1)
          refEl%Nxi3D(il, :) = temp(:, 2)
          refEl%Neta3D(il, :) = temp(:, 3)
          refEl%Nzeta3D(il, :) = temp(:, 4)
          refEl%gauss_weights3D(il) = w(i)*w(j)*w(k)*((1.d0 - x(j))/2.d0)* &
            ((1.d0 - x(k))/2.d0)**2.d0
          r = x(i)
          s = x(j)
          t = x(k)
          eta = 1.d0/2.d0*(s - s*t - 1.d0 - t)
          xi = -1.d0/2.d0*(r + 1.d0)*(eta + t) - 1.d0
          refEl%gauss_points3D(il, 1) = xi
          refEl%gauss_points3D(il, 2) = eta
          refEl%gauss_points3D(il, 3) = t
          il = il + 1
        END DO
      END DO
    END DO

    DEALLOCATE (V, L, U, Linv, Uinv)
    DEALLOCATE (p, dpxi, dpeta, dpzeta)
    DEALLOCATE (A, temp)
    DEALLOCATE (Pmatrix)

  END SUBROUTINE compute_shape_functions_thetra

  !****************************************************************************************
  !
  !                                     FEKETE NODES
  !
  !****************************************************************************************

  ! Read 2D fekete nodes from an external hd5 file
  SUBROUTINE read_fekete_nodes_2D(refEl)
    USE HDF5
    USE HDF5_io_module
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    integer(HID_T) :: file_id
    integer :: IERR
    character(50) :: datasetname

    write (datasetname, "(I2)") refEl%Ndeg
    datasetname = "P"//adjustl(trim(datasetname))

    CALL HDF5_open('positionFeketeNodesTri2D.h5', file_id, IERR)
    CALL HDF5_array2D_reading(file_id, refEl%coord2d, datasetname)
  END SUBROUTINE read_fekete_nodes_2D

  ! Hard coded positions for the fekete nodes in 1D
  SUBROUTINE feketeNodes1D(refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    ALLOCATE (refEl%coord1d(1:refEl%Nfacenodes))
    refEl%coord1d(1) = -1.d0
    refEl%coord1d(refEl%Nfacenodes) = 1.d0
    SELECT CASE (refEl%Ndeg)
    CASE (2)
      refEl%coord1d(2) = 0.d0
    CASE (3)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-4.472135955000d-1, 4.472135955000d-1/)
    CASE (4)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-6.546536707080d-1, 1.096484628801d-16, &
        6.546536707080d-1/)
    CASE (5)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-7.650553239295d-1, -2.852315164806d-1, &
        2.852315164806d-1, 7.650553239295d-1/)
    CASE (6)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-8.302238962786d-1, -4.688487934707d-1, &
        1.045388153545d-16, 4.688487934707d-1, 8.302238962786d-1/)
    CASE (7)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-8.717401485096d-1, -5.917001814331d-1, &
        -2.092992179025d-1, 2.092992179025d-1, 5.917001814331d-1, 8.717401485096d-1/)
    CASE (8)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-8.997579954115d-1, -6.771862795107d-1, &
        -3.631174638262d-1, 1.178993670744d-16, 3.631174638262d-1, 6.771862795107d-1, &
        8.997579954115d-1/)
    CASE (9)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.195339081665d-1, -7.387738651055d-1, &
        -4.779249498104d-1, -1.652789576664d-1, 1.652789576664d-1, 4.779249498104d-1, &
        7.387738651055d-1, 9.195339081665d-1/)
    CASE (10)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.340014304081d-1, -7.844834736631d-1, &
        -5.652353269962d-1, -2.957581355869d-1, 1.147670565004d-16, 2.957581355869d-1, &
        5.652353269962d-1, 7.844834736631d-1, 9.340014304081d-1/)
    CASE (11)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.448992722229d-1, -8.192793216440d-1, &
        -6.328761530319d-1, -3.995309409653d-1, -1.365529328549d-1, 1.365529328549d-1, &
        3.995309409653d-1, 6.328761530319d-1, 8.192793216440d-1, 9.448992722229d-1/)
    CASE (12)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.533098466422d-1, -8.463475646519d-1, &
        -6.861884690818d-1, -4.829098210913d-1, -2.492869301062d-1, 1.203298976554d-16, &
        2.492869301062d-1, 4.829098210913d-1, 6.861884690818d-1, 8.463475646519d-1, &
        9.533098466422d-1/)
    CASE (13)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.599350452673d-1, -8.678010538303d-1, &
        -7.288685990913d-1, -5.506394029286d-1, -3.427240133427d-1, -1.163318688837d-1, &
        1.163318688837d-1, 3.427240133427d-1, 5.506394029286d-1, 7.288685990913d-1, &
        8.678010538303d-1, 9.599350452673d-1/)
    CASE (14)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-9.652459265038d-1, -8.850820442230d-1, &
        -7.635196899518d-1, -6.062532054698d-1, -4.206380547137d-1, -2.153539553638d-1, &
        1.105861828404d-16, 2.153539553638d-1, 4.206380547137d-1, 6.062532054698d-1, &
        7.635196899518d-1, 8.850820442230d-1, 9.652459265038d-1/)
    CASE (15)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.969568046270195d0, -0.899200533093470d0, &
        -0.792008291861815d0, -0.652388702882493d0, -0.486059421887138d0, -0.299830468900763d0, &
        -0.101326273521949d0, 0.101326273521949d0, 0.299830468900763d0, 0.486059421887138d0, &
        0.652388702882493d0, 0.792008291861815d0, 0.899200533093470d0, 0.969568046270195d0/)
    CASE (16)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.973132176631391d0, -0.910879995915571d0, &
        -0.815696251221770d0, -0.691028980627685d0, -0.541385399330101d0, -0.372174433565477d0, &
        -0.189511973518318d0, 0.d0, 0.189511973518318d0, 0.372174433565477d0, 0.541385399330101d0, &
        0.691028980627685d0, 0.815696251221770d0, 0.910879995915571d0, 0.973132176631391d0/)
    CASE (17)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.976105557412166d0, -0.920649185347529d0, &
        -0.835593535218090d0, -0.723679329283243d0, -0.588504834318662d0, -0.434415036912124d0, &
        -0.266362652878281d0, -0.089749093484652d0, 0.089749093484652d0, 0.266362652878281d0, &
        0.434415036912124d0, 0.588504834318662d0, 0.723679329283243d0, 0.835593535218090d0, &
        0.920649185347529d0, 0.976105557412166d0/)
    CASE (18)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.978611766222042d0, -0.928901528152579d0, &
        -0.852460577796646d0, -0.751494202552613d0, -0.628908137265221d0, -0.488229285680714d0, &
        -0.333504847824499d0, -0.169186023409281d0, 0.d0, 0.169186023409281d0, 0.333504847824499d0, &
        0.488229285680714d0, 0.628908137265221d0, 0.751494202552613d0, 0.852460577796646d0, &
        0.928901528152579d0, 0.978611766222042d0/)
    CASE (19)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.980743704893872d0, -0.935934498812655d0, &
        -0.866877978089950d0, -0.775368260952056d0, -0.663776402290311d0, -0.534992864031886d0, &
        -0.392353183713909d0, -0.239551705922986d0, -0.080545937238822d0, 0.080545937238822d0, &
        0.239551705922987d0, 0.392353183713909d0, 0.534992864031886d0, 0.663776402290311d0, &
        0.775368260952056d0, 0.866877978089950d0, 0.935934498812655d0, 0.980743704893872d0/)
    CASE (20)
      refEl%coord1d(2:refEl%Nfacenodes - 1) = (/-0.982572296604501d0, -0.941976296959732d0, &
        -0.879294755323590d0, -0.796001926077712d0, -0.694051026062223d0, -0.575831960261831d0, &
        -0.444115783279002d0, -0.301989856508765d0, -0.152785515802186d0, 0.d0, 0.152785515802186d0, &
        0.301989856508765d0, 0.444115783279002d0, 0.575831960261831d0, 0.694051026062223d0, &
        0.796001926077712d0, 0.879294755323590d0, 0.941976296959732d0, 0.982572296604501d0/)
    END SELECT

  END SUBROUTINE feketeNodes1D

  !****************************************************************************************
  !
  !                                 GAUSS QUADRATURES
  !
  !****************************************************************************************

  !******************************
  ! Gauss Legendre quadratures 1D
  !******************************
  SUBROUTINE gauss_legendre(N, x, w)

    integer*4, intent(in) :: N
    real*8, intent(out) :: x(1:N), w(1:N)
    integer*4 :: m, i, j
    real*8 :: z, z1, p1, p2, p3, pp
    real*8, parameter :: eps = 3.d-15

    ! Roots are symmetrical in the interval (-1,1), so only need to find half of them
    m = (N + 1)/2
    DO i = 1, m
      z = cos(pi*(dble(i) - 0.25d0)/(dble(N) + 0.5d0)) ! initial guess of ith root
      ! refine guess by Newton's method
      100      p1 = 1.d0
      p2 = 0.d0
      DO j = 1, N
        p3 = p2
        p2 = p1
        p1 = ((2.d0*dble(j) - 1.d0)*z*p2 - (dble(j) - 1.d0)*p3)/dble(j)
      END DO
      ! p1 is now the Legendre polynomial
      ! Calculate now pp (its derivative) using p2, the polynomial of one lower order
      pp = dble(N)*(z*p1 - p2)/(z*z - 1.d0)
      z1 = z
      z = z1 - p1/pp
      IF (dabs(z - z1) .gt. eps) goto 100

      x(i) = -z
      x(N + 1 - i) = z
      w(i) = 2.d0/((1.d0 - z*z)*pp*pp)
      w(N + 1 - i) = w(i)
    END DO

  END SUBROUTINE gauss_legendre

  !********************************
  ! Gauss cubature 2D for triangles
  !********************************
  SUBROUTINE gauss_legendre_cubature2d_new(order, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    integer*4, intent(in) :: order
    real*8 :: a, b, c, F(1:2, 1:2), u
    real*8, allocatable :: pospg(:, :), pespg(:)
    integer*4 :: gpts, i

    ! order is the Cubature order
    ! gpts is the number of Gauss points
    SELECT CASE (order)
    CASE (1)
      gpts = 1
    CASE (2)
      gpts = 3
    CASE (3)
      gpts = 4
    CASE (4)
      gpts = 6
    CASE (5)
      gpts = 7
    CASE (6)
      gpts = 12
    CASE (7)
      gpts = 13
    CASE (8)
      gpts = 16
    CASE (10)
      gpts = 25
    CASE (15)
      gpts = 54
    CASE (20)
      gpts = 85
    CASE (25)
      gpts = 126
    CASE (30)
      gpts = 175
    END SELECT

    ALLOCATE (refEl%gauss_points2D(1:gpts, 1:2))
    ALLOCATE (refEl%gauss_weights2D(1:gpts))
    refEl%gauss_points2D = 0.d0
    refEl%gauss_weights2D = 0.d0

    SELECT CASE (gpts)
    CASE (1)
      a = 1.d0/3.d0
      refEl%gauss_points2D(1, 1:2) = a
      refEl%gauss_weights2D = 1.d0
    CASE (3)
      a = 1.d0/6.d0
      refEl%gauss_points2D(:, :) = a
      a = 2.d0/3.d0
      refEl%gauss_points2D(1, 1) = a
      refEl%gauss_points2D(2, 2) = a
      a = 1.d0/3.d0
      refEl%gauss_weights2D(:) = a
    CASE (4)
      a = 1.d0/3.d0
      refEl%gauss_points2D(:, :) = 0.2d0
      refEl%gauss_points2D(2, 1) = 0.6d0
      refEl%gauss_points2D(3, 2) = 0.6d0
      refEl%gauss_points2D(1, :) = a
      a = 25.d0/48.d0
      refEl%gauss_weights2D(:) = a
      a = -27.d0/48.d0
      refEl%gauss_weights2D(1) = a
    CASE (6)
      a = 0.659027622374092d0
      b = 0.231933368553031d0
      c = 0.109039009072877d0
      refEl%gauss_points2D(1, :) = (/a, b/)
      refEl%gauss_points2D(2, :) = (/b, a/)
      refEl%gauss_points2D(3, :) = (/a, c/)
      refEl%gauss_points2D(4, :) = (/c, a/)
      refEl%gauss_points2D(5, :) = (/b, c/)
      refEl%gauss_points2D(6, :) = (/c, b/)
      a = 1.d0/6.d0
      refEl%gauss_weights2D(:) = a
    CASE (12)
      a = 0.873821971016996d0
      b = 0.063089014491502d0
      refEl%gauss_points2D(1, :) = (/a, b/)
      refEl%gauss_points2D(2, :) = (/b, b/)
      refEl%gauss_points2D(3, :) = (/b, a/)
      a = 0.501426509658179d0
      b = 0.249286745170910d0
      refEl%gauss_points2D(4, :) = (/a, b/)
      refEl%gauss_points2D(5, :) = (/b, b/)
      refEl%gauss_points2D(6, :) = (/b, a/)
      a = 0.636502499121399d0
      b = 0.310352451033785d0
      c = 0.053145049844816d0
      refEl%gauss_points2D(7, :) = (/a, b/)
      refEl%gauss_points2D(8, :) = (/b, a/)
      refEl%gauss_points2D(9, :) = (/a, c/)
      refEl%gauss_points2D(10, :) = (/c, a/)
      refEl%gauss_points2D(11, :) = (/b, c/)
      refEl%gauss_points2D(12, :) = (/c, b/)
      a = 0.050844906370207d0
      b = 0.116786275726379d0
      c = 0.082851075618374d0
      refEl%gauss_weights2D(1:3) = a
      refEl%gauss_weights2D(4:6) = b
      refEl%gauss_weights2D(7:12) = c
    CASE (7, 25, 54, 85, 126, 175)
      ALLOCATE (pospg(1:gpts, 1:2), pespg(1:gpts))
      CALL quadrature_tri_new(order, gpts, pospg, pespg)
      u = 1.d0/3.d0
      F = 0.d0
      F(1, 1) = 2.d0*u
      F(2, 1) = -u
      F(2, 2) = sqrt(u)
      pospg = matmul(pospg, transpose(F))
      DO i = 1, gpts
        pospg(i, :) = pospg(i, :) + u
      END DO
      pespg(:) = pespg(:)/sum(pespg)
      refEl%gauss_points2D(:, :) = pospg(:, :)
      refEl%gauss_weights2D(:) = pespg(:)
      DEALLOCATE (pospg, pespg)
    END SELECT

  END SUBROUTINE gauss_legendre_cubature2d_new

  SUBROUTINE quadrature_tri_new(order, n1, pospg, pespg)

    integer*4, intent(in) :: order, n1
    real*8, intent(out) :: pospg(1:n1, 1:2), pespg(1:n1)
    integer*4 :: i, j, k, l, nn, C(1:2, 1:2), E(1:2, 1:2)
    real*8 :: fi, A(1:2, 1:2), B(1:2, 1:2), D(1:2, 1:2), F(1:2, 1:2), coord(1:2)
    real*8, allocatable :: X(:), Y(:), w(:)
    integer*4, allocatable :: M(:)

    fi = 2.d0*pi/3.d0

    A(1, 1) = cos(fi)
    A(1, 2) = sin(fi)
    A(2, 1) = -A(1, 2)
    A(2, 2) = A(1, 1)

    B(1, 1) = cos(2.d0*fi)
    B(1, 2) = sin(2.d0*fi)
    B(2, 1) = -B(1, 2)
    B(2, 2) = B(1, 1)

    C(1, 1) = 1
    C(1, 2) = 0
    C(2, 1) = 0
    C(2, 2) = 1

    D = A
    D(1, 2) = -D(1, 2)
    D(2, 2) = -D(2, 2)

    E = C
    E(2, 2) = -E(2, 2)

    F = B
    F(1, 2) = -F(1, 2)
    F(2, 2) = -F(2, 2)

    SELECT CASE (order)
    CASE (5)
      nn = 3
    CASE (10)
      nn = 7
    CASE (15)
      nn = 12
    CASE (20)
      nn = 19
    CASE (25)
      nn = 26
    CASE (30)
      nn = 36
    END SELECT

    ALLOCATE (M(1:nn), X(1:nn), Y(1:nn), w(1:nn))
    M = 0
    X = 0.d0
    Y = 0.d0
    w = 0.d0

    SELECT CASE (order)
    CASE (5)
      M = (/1, 3, 3/)
      X = (/0.d0, -0.4104261923153453d0, 0.6961404780296310d0/)
      Y = (/0.d0, 0.d0, 0.d0/)
      w = (/0.225, 0.13239415278850623d0, 0.12593918054482713d0/)
    CASE (10)
      M = (/1, 3, 3, 3, 3, 6, 6/)
      X = (/0.d0, -0.49359629886342453d0, -0.28403734918716863d0, 0.44573076177032633d0, &
        0.9385563442849673d0, -0.4474955151540920d0, -0.4436763946123360d0/)
      Y = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -0.5991595522781586d0, -0.2571781329392130d0/)
      w = (/0.08352339980519638d0, 0.007229850592056743d0, 0.07449217792098051d0, &
        0.07864647340310853d0, 0.006928323087107504d0, 0.0295183203347794d0, 0.03957936719606124d0/)
    CASE (15)
      M = (/3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6/)
      X = (/-0.3748423891073751d0, -0.2108313937373917d0, 0.12040849626092393d0, 0.5605966391716812d0, &
        0.8309113970031897d0, 0.95027461942488903d0, -0.4851316950361628d0, -0.4762943440546580d0, &
        -0.4922845867745440d0, -0.4266165113705168d0, -0.3968468770512212d0, -0.2473933728129512d0/)
      Y = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -0.4425551659467111d0, -0.1510682717598242d0, &
        -0.6970224211436132d0, -0.5642774363966393d0, -0.3095105740458471d0, -0.2320292030461791d0/)
      w = (/0.03266181884880529d0, 0.027412818031364363d0, 0.02651003659870330d0, 0.02921596213648611d0, &
        0.01058460806624399d0, 0.003614643064092035d0, 0.008527748101709436d0, 0.01391617651669193d0, &
        0.004291932940734835d0, 0.01623532928177489d0, 0.02560734092126239d0, 0.033088195531645673d0/)
    CASE (20)
      M = (/1, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/)
      X = (/0.d0, -0.4977490260133565d0, -0.35879037209157373d0, -0.19329181386571043d0, &
        0.20649939240163803d0, 0.3669431077237697d0, 0.6767931784861860d0, 0.8827927364865920d0, &
        0.9664768608120111d0, -0.4919755727189941d0, -0.4880677744007016d0, -0.4843664025781043d0, &
        -0.4835533778058150d0, -0.4421499318718065d0, -0.44662923827417273d0, -0.4254937754558538d0, &
        -0.4122204123735024d0, -0.31775331949340863d0, -0.2889337325840919d0/)
      Y = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -0.7513212483763635d0, -0.5870191642967427d0, &
        -0.17172709841143283d0, -0.3833898305784408d0, -0.6563281974461070d0, -0.06157647932662624d0, &
        -0.47831240826600273d0, -0.2537089901614676d0, -0.3996183176834929d0, -0.18441839672339823d0/)
      w = (/0.02761042699769952d0, 0.00177902954732674d0, 0.02011239811396117d0, 0.02681784725933157d0, &
        0.02452313380150201d0, 0.01639457841069539d0, 0.01479590739864960d0, 0.004579282277704251d0, &
        0.001651826515576217d0, 0.002349170908575584d0, 0.004465925754181793d0, 0.006099566807907972d0, &
        0.006891081327188203d0, 0.007997475072478163d0, 0.0073861342853360243d0, 0.01279933187864826d0, &
        0.01725807117569655d0, 0.01867294590293547d0, 0.02281822405839526d0/)
    CASE (25)
      M = (/3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/)
      X = (/-0.4580802753902387, -0.3032320980085228, -0.1696674057318916, 0.1046702979405866, 0.2978674829878846, &
        0.5455949961729473, 0.6617983193620190, 0.7668529237254211, 0.8953207191571090, 0.9782254461372029, &
        -0.4980614709433367, -0.4919004480918257, -0.4904239954490375, -0.4924576827470104, -0.4897598620673272, &
        -0.4849757005401057, -0.4613632802399150, -0.4546581528201263, -0.4542425148392569, -0.4310651789561460, &
        -0.3988357991895837, -0.3949323628761341, -0.3741327130398251, -0.3194366964842710, -0.2778996512639500, -0.2123422011990124/)
      Y = (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.4713592181681879, -0.1078887424748246, -0.3057041948876942, -0.7027546250883238, &
        -0.7942765584469995, -0.5846826436376921, -0.4282174042835178, -0.2129434060653430, -0.6948910659636692, -0.5691146659505208, &
        -0.3161666335733065, -0.1005941839340892, -0.4571406037889341, -0.2003599744104858, -0.3406754571040736, -0.1359589640107579/)
      w =(/ 0.008005581880020417,  0.01594707683239050,   0.01310914123079553,  0.019583000965635623,  0.01647088544153727,  0.008547279074092100, &
        0.008161885857226492,  0.0061211465399837793, 0.002908498264936665, 0.0006922752456619963, 0.001248289199277397, 0.003404752908803022, &
        0.003359654326064051,  0.001716156539496754,  0.001480856316715606, 0.003511312610728685,  0.007393550149706484, 0.007983087477376558, &
        0.0043559626131580413, 0.007365056701417832, 0.01096357284641955, 0.01174996174354112, 0.01001560071379857, 0.01330964078762868, &
        0.01415444650522614, 0.01488137956116801/)
    CASE (30)
      M = (/1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6/)
      X = (/0., -0.489004825350852, -0.375506486295553, -0.273528565811884, -0.146141210161750, 0.157036462611772, &
        0.317953072437897, 0.476322665473811, 0.630224718395690, 0.759747313323409, 0.856676597776304, &
        0.934838455959575, 0.985705967153689, -0.498611943209980, -0.497921111216654, -0.494476376816134, &
        -0.494145164863761, -0.495150127767484, -0.490298851831645, -0.495128786763001, -0.486987363789869, &
        -0.476604460299029, -0.473034918119472, -0.474313631969166, -0.465674891980127, -0.450893604068350, &
        -0.449268481486489, -0.446678578309977, -0.424190314539700, -0.414477927626402, -0.403770790368195, &
        -0.379248277568562, -0.343449397798204, -0.329232658356873, -0.281954768426714, -0.215081520767032/)
      Y = (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.145911499458133, -0.758841124126978, -0.577206108525577, &
        -0.819283113385993, -0.333106124712369, -0.674968075724015, -0.464914848460198, -0.274747981868076, &
        -0.755078734433048, -0.153390877058151, -0.429173048901523, -0.559744628102069, -0.665677920960733, &
        -0.302435402004506, -0.037339333379264, -0.443257445349149, -0.159839002260082, -0.562852040975635, &
        -0.304872368029416, -0.434881627890658, -0.151014758677329, -0.290117766854826, -0.143940337075373/)
      w = (/0.015579960202899, 0.003177233700534, 0.010483426635731, 0.013209459577744, 0.014975006966272, 0.014987904443384, &
        0.013338864741022, 0.010889171113902, 0.008189440660893, 0.005575387588608, 0.003191216473412, 0.001296715144327, &
        0.000298262826135, 0.000998905685079, 0.000462850849173, 0.001234451336382, 0.000570719852243, 0.001126946125878, &
        0.001747866949407, 0.001182818815032, 0.001990839294675, 0.001900412795036, 0.004498365808817, 0.003478719460275, &
        0.004102399036724, 0.004021761549744, 0.006033164660795, 0.003946290302130, 0.006644044537680, 0.008254305856078, &
        0.006496056633406, 0.009252778144147, 0.009164920726294, 0.011569524628098, 0.011761116467609, 0.013824702182165/)
    CASE DEFAULT
      error stop "Quadrature not yet implemented"
    END SELECT

    k = 1
    l = 1
    DO i = 1, nn
      coord(1) = X(i)
      coord(2) = Y(i)
      pospg(k, :) = matmul(dble(C), coord)
      k = k + 1
      IF (M(i) .gt. 1) THEN
        pospg(k, :) = matmul(A, coord)
        k = k + 1
        pospg(k, :) = matmul(B, coord)
        k = k + 1
      END IF
      IF (M(i) .gt. 3) THEN
        pospg(k, :) = matmul(D, coord)
        k = k + 1
        pospg(k, :) = matmul(dble(E), coord)
        k = k + 1
        pospg(k, :) = matmul(F, coord)
        k = k + 1
      END IF
      DO j = 1, M(i)
        pespg(l) = w(i)*sqrt(27.d0)/4.d0
        l = l + 1
      END DO
    END DO

    DEALLOCATE (M, X, Y, w)

  END SUBROUTINE quadrature_tri_new

  !****************************************************************************************
  !
  !                                 VANDERMONDE MATRIX
  !
  !****************************************************************************************
  SUBROUTINE vandermonde_1d(V, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    real*8, intent(out) :: V(1:refEl%Nfacenodes, 1:refEl%Nfacenodes)
    real*8 :: p(1:refEl%Nfacenodes)
    integer*4 :: i

    V = 0.d0
    DO i = 1, refEl%Nfacenodes
      CALL orthopoly1d(refEl%coord1d(i), refEl%Ndeg, p)
      V(i, :) = p(:)
    END DO

  END SUBROUTINE vandermonde_1d

  SUBROUTINE vandermonde_2d(V, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    real*8, intent(out) :: V(1:refEl%Nnodes2D, 1:refEl%Nnodes2D)
    real*8 :: r, s, p(1:refEl%Nnodes2D)
    integer*4 :: i

    V = 0.d0
    DO i = 1, refEl%Nnodes2D
      r = refEl%coord2d(i, 1)
      s = refEl%coord2d(i, 2)
      CALL orthopoly2d(r, s, refEl%Nnodes2D, refEl%Ndeg, p)
      V(i, :) = p(:)
    END DO

  END SUBROUTINE vandermonde_2d

  SUBROUTINE vandermonde_3d(V, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl

    real*8, intent(out) :: V(1:refEl%Nnodes3D, 1:refEl%Nnodes3D)
    real*8 :: r, s, t, p(1:refEl%Nnodes3D)
    integer*4 :: i

    V = 0.d0
    DO i = 1, refEl%Nnodes3D
      r = refEl%coord3d(i, 1)
      s = refEl%coord3d(i, 2)
      t = refEl%coord3d(i, 3)
      CALL orthopoly3d(r, s, t, refEl%Nnodes3D, refEl%Ndeg, p)
      V(i, :) = p(:)
    END DO

  END SUBROUTINE vandermonde_3d

  ! Vandermonde matrix for quadrilateral elements.
  ! Added 12/09/2018
  SUBROUTINE vandermonde_qua(V, refEl)
    TYPE(Reference_element_type), intent(INOUT) :: refEl
    real*8, intent(out) :: V(refEl%Nnodes2D, refEl%Nnodes2D)
    real*8 :: px(refEl%Nnodes1D), py(refEl%Nnodes1D)
    real*8 :: aux(refEl%Nnodes2D, refEl%Nnodes2D)
    integer*4 :: i, j, ind, perm(refEl%Nnodes2D)

    aux = 0.d0
    DO j = 1, refEl%Nnodes1D
      CALL orthopoly1d(refEl%coord1d(j), refEl%Ndeg, py)
      DO i = 1, refEl%Nnodes1D
        CALL orthopoly1d(refEl%coord1d(i), refEl%Ndeg, px)
        ind = (j - 1)*refEl%Nnodes1d+i
        aux(ind, :) = reshape(tensorProduct(px, py), (/refEl%Nnodes1D**2/))
      END DO
    END DO

    perm = permutation_quads(refEl)

    DO i = 1, refEl%Nnodes2D
      V(perm(i), :) = aux(i, :)
    END DO
  END SUBROUTINE vandermonde_qua

  !****************************************************************************************
  !
  !                                 ORTHOGONAL POLYNOMIALS
  !
  !****************************************************************************************
  !******************************
  ! Jacobi polynomials
  !******************************
  SUBROUTINE Jacobi_P(x, a, b, N, p)

    integer*4, intent(in) :: N
    real*8, intent(in) :: x, a, b
    real*8, intent(out) :: p
    integer*4 :: m
    real*8 :: c, d, pm1, pm2, temp1, temp2

    IF (N .eq. 0) THEN
      p = 1
    ELSEIF (N .eq. 1) THEN
      p = 1.d0/2.d0*(a - b + (2.d0 + a + b)*x)
    ELSE
      c = a + b
      d = (a + b)*(a - b)
      pm2 = 1.d0
      pm1 = 1.d0/2.d0*(a - b + (2.d0 + a + b)*x)
      DO m = 2, N
        temp1 = 2.d0*dble(m) + c
        temp2 = dble(m)*(dble(m) + c)*(temp1 - 2.d0)
        p = ((temp1 - 1.d0)*(d + x*(temp1 - 2.d0)*temp1)/(2.d0*temp2))*pm1 - &
          ((dble(m) + a - 1.d0)*(dble(m) + b - 1.d0)*temp1/temp2)*pm2
        pm2 = pm1
        pm1 = p
      END DO
    END IF

  END SUBROUTINE Jacobi_P

  SUBROUTINE orthopoly1d(x, N, p)

    integer*4, intent(in):: N
    real*8, intent(in) :: x
    real*8, intent(out) :: p(1:N + 1)
    integer*4 :: i

    p = 0.d0
    DO i = 0, N
      CALL Jacobi_P(x, 0.d0, 0.d0, i, p(i + 1))
      p(i + 1) = p(i + 1)*sqrt((2.d0*dble(i) + 1.d0)/2.d0)
    END DO

  END SUBROUTINE orthopoly1d

  SUBROUTINE orthopoly2d(xi, eta, N, m, p)

    real*8, intent(in) :: xi, eta
    integer*4, intent(in) :: N, m
    real*8, intent(out) :: p(1:N)
    real*8 :: r, s
    real*8, parameter :: tol = 1.e-12

    IF (abs(eta - 1.d0) .le. tol) THEN
      r = -1.d0
      s = 1.d0
    ELSE
      r = 2.d0*(1.d0 + xi)/(1.d0 - eta) - 1.d0
      s = eta
    END IF
    CALL orthopoly2d_rst(r, s, N, m, p)

  END SUBROUTINE orthopoly2d

  SUBROUTINE orthopoly2d_rst(r, s, N, m, p)

    integer*4, intent(in) :: N, m
    real*8, intent(in) :: r, s
    real*8, intent(out) :: p(1:N)
    integer*4 :: ncount, i, j, k
    real*8 :: pi, qi, pj, factor

    p = 0.d0
    ncount = 0
    DO k = 0, m
      DO i = 0, k
        IF (i .eq. 0) THEN
          pi = 1.d0
          qi = 1.d0
        ELSE
          CALL Jacobi_P(r, 0.d0, 0.d0, i, pi)
          qi = qi*(1.d0 - s)/2.d0
        END IF
        j = k - i
        IF (j .eq. 0) THEN
          pj = 1.d0
        ELSE
          CALL Jacobi_P(s, 2.d0*dble(i) + 1.d0, 0.d0, j, pj)
        END IF
        ncount = ncount + 1
        factor = sqrt((2.d0*dble(i) + 1.d0)*(dble(i) + dble(j) + 1.d0)/2.d0)
        p(ncount) = pi*qi*pj*factor
      END DO
    END DO

  END SUBROUTINE orthopoly2d_rst

  SUBROUTINE orthopoly3d(xi, eta, zeta, N, m, p)

    integer*4, intent(in) :: N, m
    real*8, intent(in) :: xi, eta, zeta
    real*8, intent(out) :: p(1:N)
    real*8 :: r, s
    real*8, parameter :: tol = 1e-12

    IF (abs(eta + zeta) .le. tol) THEN
      r = -1.d0
      s = 1.d0
    ELSEIF (abs(zeta - 1.0) .le. tol) THEN
      r = -1.d0
      s = 1.d0
    ELSE
      r = -2.d0*(1.d0 + xi)/(eta + zeta) - 1.d0
      s = 2.d0*(1.d0 + eta)/(1.d0 - zeta) - 1.d0
    END IF

    CALL orthopoly3d_rst(r, s, zeta, N, m, p)

  END SUBROUTINE orthopoly3d

  SUBROUTINE orthopoly3d_rst(r, s, t, N, nn, p)

    integer*4, intent(in) :: N, nn
    real*8, intent(in) :: r, s, t
    real*8, intent(out) :: p(1:N)
    integer*4 :: i, j, k, l, m
    real*8 :: pi, qi, pj, qj, pk, factor

    l = 0
    DO m = 0, nn
      DO i = 0, m
        IF (i .eq. 0) THEN
          pi = 1.d0
          qi = 1.d0
        ELSE
          CALL Jacobi_P(r, 0.d0, 0.d0, i, pi)
          qi = qi*(1.d0 - s)/2.d0
        END IF
        DO j = 0, m - i
          IF (j .eq. 0) THEN
            pj = 1.d0
            qj = ((1.d0 - t)/2.d0)**dble(i)
          ELSE
            CALL Jacobi_P(s, 2.d0*dble(i) + 1.d0, 0.d0, j, pj)
            qj = qj*(1.d0 - t)/2.d0
          END IF
          k = m - (i + j)
          IF (k .eq. 0) THEN
            pk = 1.d0
          ELSE
            CALL Jacobi_P(t, 2.d0*dble(i + j) + 2.d0, 0.d0, k, pk)
          END IF
          l = l + 1
          factor = sqrt((2.d0*dble(i) + 1.d0)*(dble(i) + dble(j) + 1.d0)*(2.d0*dble(i + j + k) + 3.d0)/4.d0)
          p(l) = pi*qi*pj*qj*pk*factor
        END DO
      END DO
    END DO

  END SUBROUTINE orthopoly3d_rst

  SUBROUTINE orthopoly1d_deriv(x, n, a, b)

    real*8, intent(in) :: x
    integer*4, intent(in) :: n
    real*8, intent(out) :: a(1:n + 1), b(1:n + 1)
    real*8 :: factor
    integer*4 :: i

    a = 0.d0
    b = 0.d0
    a(1) = 1.d0/sqrt(2.d0)

    DO i = 1, n
      factor = sqrt((2.d0*dble(i) + 1.d0)/2.d0)
      CALL Jacobi_P(x, 0.d0, 0.d0, i, a(i + 1))
      a(i + 1) = a(i + 1)*factor
      CALL Jacobi_P(x, 1.d0, 1.d0, i - 1, b(i + 1))
      b(i + 1) = b(i + 1)*((dble(i) + 1.d0)/2.d0)*factor
    END DO

  END SUBROUTINE orthopoly1d_deriv

  SUBROUTINE orthopoly2d_deriv(xi, eta, n, m, p, dp_dxi, dp_deta)

    real*8, intent(in) :: xi, eta
    integer*4, intent(in) :: n, m
    real*8 :: p(1:m), dp_dxi(1:m), dp_deta(1:m)
    real*8 :: r, s
    real*8, parameter :: tol = 1e-12

    IF (abs(eta - 1.0) .le. tol) THEN
      r = -1.d0
      s = 1.d0
    ELSE
      r = 2.d0*(1.d0 + xi)/(1.d0 - eta) - 1.d0
      s = eta
    END IF

    CALL orthopoly2d_deriv_rst(r, s, n, m, p, dp_dxi, dp_deta)

  END SUBROUTINE orthopoly2d_deriv

  SUBROUTINE orthopoly2d_deriv_rst(r, s, n, m, p, dp_dxi, dp_deta)

    real*8, intent(in) :: r, s
    integer*4, intent(in) :: n, m
    real*8, intent(out) :: p(1:m), dp_dxi(1:m), dp_deta(1:m)
    integer*4 :: ncount, i, j, k
    real*8 :: xi, eta, dr_dxi, dr_deta, pi, pj, dpi, dpj, qi, dqi, factor, dp_dr, dp_ds

    xi = (1.d0 + r)*(1.d0 - s)/2.d0 - 1.d0
    eta = s

    dr_dxi = 2.d0/(1.d0 - eta)
    dr_deta = 2.d0*(1.d0 + xi)/(1.d0 - eta)**2.d0

    ncount = 0
    DO k = 0, n
      DO i = 0, k
        IF (i .eq. 0) THEN
          pi = 1.d0
          qi = 1.d0
          dpi = 0.d0
          dqi = 0.d0
        ELSE
          CALL Jacobi_P(r, 0.d0, 0.d0, i, pi)
          CALL Jacobi_P(r, 1.d0, 1.d0, i - 1, dpi)
          dpi = dpi*(dble(i) + 1.d0)/2.d0
          qi = qi*(1 - s)/2.d0
          dqi = qi*(-dble(i))/(1.d0 - s)
        END IF
        j = k - i
        IF (j .eq. 0) THEN
          pj = 1.d0
          dpj = 0.d0
        ELSE
          CALL Jacobi_P(s, 2.d0*dble(i) + 1.d0, 0.d0, j, pj)
          CALL Jacobi_P(s, 2.d0*dble(i) + 2.d0, 1.d0, j - 1, dpj)
          dpj = dpj*(dble(j) + 2.d0*dble(i) + 2.d0)/2.d0
        END IF
        ncount = ncount + 1
        factor = sqrt((2.d0*dble(i) + 1.d0)*(dble(i) + dble(j) + 1.d0)/2.d0)
        p(ncount) = pi*qi*pj*factor
        dp_dr = dpi*qi*pj*factor
        dp_ds = pi*(dqi*pj + qi*dpj)*factor
        dp_dxi(ncount) = dp_dr*dr_dxi
        dp_deta(ncount) = dp_dr*dr_deta + dp_ds
      END DO
    END DO

  END SUBROUTINE orthopoly2d_deriv_rst

  SUBROUTINE orthopoly3d_deriv(r, s, t, n, m, p, dp_dxi, dp_deta, dp_dzeta)

    integer*4, intent(in) :: n, m
    real*8, intent(in) :: r, s, t
    real*8, intent(out) :: p(1:m), dp_dxi(1:m), dp_deta(1:m), dp_dzeta(1:m)
    integer*4 :: i, j, k, l, z
    real*8 :: xi, eta, zeta, dr_dxi, dr_deta, dr_dzeta, ds_deta, ds_dzeta
    real*8 :: pi, qi, dpi, dqi, pj, qj, dpj, dqj, pk, dpk, factor, dp_dr, dp_ds, dp_dt

    eta = 1.d0/2.d0*(s - s*t - 1.d0 - t)
    xi = -1.d0/2.d0*(r + 1.d0)*(eta + t) - 1.d0
    zeta = t

    dr_dxi = -2.d0/(eta + zeta)
    dr_deta = 2.d0*(1.d0 + xi)/(eta + zeta)**2.d0
    dr_dzeta = dr_deta

    ds_deta = 2.d0/(1.d0 - zeta)
    ds_dzeta = 2.d0*(1.d0 + eta)/(1.d0 - zeta)**2.d0

    z = 0
    DO l = 0, n
      DO i = 0, l
        IF (i .eq. 0) THEN
          pi = 1.d0
          qi = 1.d0
          dpi = 0.d0
          dqi = 0.d0
        ELSE
          CALL Jacobi_P(r, 0.d0, 0.d0, i, pi)
          CALL Jacobi_P(r, 1.d0, 1.d0, i - 1, dpi)
          dpi = dpi*(dble(i) + 1.d0)/2.d0
          qi = qi*(1.d0 - s)/2.d0
          dqi = -qi*dble(i)/(1.d0 - s)
        END IF
        DO j = 0, l - i
          IF (j .eq. 0) THEN
            pj = 1.d0
            qj = ((1.d0 - t)/2.d0)**i
            dpj = 0.d0
            dqj = -qj*(dble(i) + dble(j))/(1.d0 - t)
          ELSE
            CALL Jacobi_P(s, 2.d0*dble(i) + 1.d0, 0.d0, j, pj)
            CALL Jacobi_P(s, 2.d0*dble(i) + 2.d0, 1.d0, j - 1, dpj)
            dpj = dpj*(dble(j) + 2.d0*dble(i) + 2.d0)/2.d0
            qj = qj*(1.d0 - t)/2.d0
            dqj = -qj*(dble(i) + dble(j))/(1.d0 - t)
          END IF
          k = l - (i + j)
          IF (k .eq. 0) THEN
            pk = 1.d0
            dpk = 0.d0
          ELSE
            CALL Jacobi_P(t, 2.d0*dble(i + j) + 2.d0, 0.d0, k, pk)
            CALL Jacobi_P(t, 2.d0*dble(i + j) + 3.d0, 1.d0, k - 1, dpk)
            dpk = dpk*(k + 2.d0*dble(i) + 2.d0*dble(j) + 3.d0)/2.d0
          END IF
          z = z + 1
          factor = sqrt((2.d0*dble(i) + 1.d0)*(dble(i) + dble(j) + 1.d0)*(2.d0*dble(i + j + k) + 3.d0)/4.d0)
          p(z) = pi*qi*pj*qj*pk*factor
          dp_dr = dpi*qi*pj*qj*pk*factor
          dp_ds = pi*(dqi*pj + qi*dpj)*qj*pk*factor
          dp_dt = pi*qi*pj*(dqj*pk + qj*dpk)*factor
          dp_dxi(z) = dp_dr*dr_dxi
          dp_deta(z) = dp_dr*dr_deta + dp_ds*ds_deta
          dp_dzeta(z) = dp_dr*dr_dzeta + dp_ds*ds_dzeta + dp_dt
        END DO
      END DO
    END DO

  END SUBROUTINE orthopoly3d_deriv

  SUBROUTINE create_toroidal_structures(refElTor, refElPol)
    TYPE(Reference_element_type), intent(INOUT) :: refElTor
    TYPE(Reference_element_type), intent(IN)   :: refElPol
    integer :: i, j, k, g, ifa, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor
    integer, allocatable :: faceNodes1d(:)
    !***********************************************
    ! Generate 3D basis functions
    !***********************************************
    ALLOCATE (refElTor%N3D(refElPol%NGauss2D*refElTor%NGauss1D, refElPol%Nnodes2D*refElTor%Nnodes1D))
    DO i = 1, refElTor%NGauss1D
      DO j = 1, refElPol%NGauss2D
        g = (i - 1)*refElPol%NGauss2D+j
        refElTor%N3D(g, :) = col(TensorProduct(refElPol%N2D(j, :), refElTor%N1D(i, :)))
      END DO
    END DO

    !***********************************************
    ! Generate face nodes for toroidal faces
    !***********************************************
    ALLOCATE (refElTor%faceNodes3(refElPol%Nfaces, size(refElPol%face_nodes, 2)*size(refElTor%face_nodes, 2)))
    ALLOCATE (faceNodes1d(refElTor%Nnodes1D))
    DO ifa = 1, refElPol%Nfaces
      faceNodes1d = (/(i, i=1, refElTor%Nnodes1D)/)
      refElTor%faceNodes3(ifa, :) = colint(TensorSumInt(refElPol%face_nodes(ifa, :), (faceNodes1d-1)*refElPol%Nnodes2D))
    END DO
    DEALLOCATE (faceNodes1d)

    !***********************************************
    ! Generate shape functions for toroidal faces
    !***********************************************
    Np1dPol = size(refElPol%N1d, 2)
    Np1dTor = size(refElTor%N1d, 2)
    Ng1dPol = size(refElPol%N1d, 1)
    Ng1dTor = size(refElTor%N1d, 1)
    refElTor%Nfl = Np1dPol*Np1dTor ! Number of nodes for lateral faces
    refElTor%Ngl = Ng1dPol*Ng1dTor ! Number of Gauss points for lateral faces
    refElTor%Nnodes3D = refElPol%Nnodes2D*refElTor%Nnodes1D ! Number of nodes in the 3D element
    refElTor%NGauss3D = refElPol%NGauss2D*refElTor%NGauss1D ! Number of Gauss points in the 3D element
    refElTor%Nft = refElTor%Nfl*refElPol%Nfaces + 2*refElPol%Nnodes2d         ! Number of nodes in all the faces of the 3D element
    ALLOCATE (refElTor%sFTF(Ng1dPol*Ng1dTor, Np1dPol*Np1dTor))
    DO i = 1, Np1dTor
      DO j = 1, Np1dPol
        k = (i - 1)*Np1dPol + j
        refElTor%sFTF(:, k) = col(tensorProduct(refElPol%N1d(:, j), refElTor%N1d(:, i)))
      END DO
    END DO
  END SUBROUTINE create_toroidal_structures

END MODULE reference_element

