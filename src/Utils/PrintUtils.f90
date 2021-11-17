!************************************************************
! project: MHDG
! file: PrintUtils.f90
! date: 19/09/2016
! Some printing routines
!************************************************************

MODULE Printutils

CONTAINS

  !***************************
  ! Print a vector of integers
  !***************************
  SUBROUTINE displayVectorInt(Vec)
    INTEGER, DIMENSION(:), INTENT(IN) :: Vec
    INTEGER :: di, i

    di = SIZE(Vec)

    DO i = 1, di
      WRITE (6, '(I7)') Vec(i)
    END DO   ! end i
    WRITE (6, *)
    WRITE (6, *)
  END SUBROUTINE displayVectorInt

  !***************************
  ! Print a vector of reals
  !***************************
  SUBROUTINE displayVector(Vec)
    REAL, DIMENSION(:), INTENT(IN) :: Vec
    INTEGER :: di, i

    di = SIZE(Vec)

    DO i = 1, di
      WRITE (6, '(ES25.15)') Vec(i)
    END DO   ! end i
    WRITE (6, *)
    WRITE (6, *)
  END SUBROUTINE displayVector

  !***************************
  ! Print a matrix of integers
  !***************************
  SUBROUTINE displayMatrixInt(Mat)
    INTEGER, DIMENSION(:, :), INTENT(IN) :: Mat
    INTEGER :: di, dj, i, j

    di = SIZE(Mat, dim=1)
    dj = SIZE(Mat, dim=2)

    DO i = 1, di
      DO j = 1, dj
        WRITE (6, '(I5)', advance="no") Mat(i, j)
      END DO  ! end j
      WRITE (6, *)
    END DO   ! end i
    WRITE (6, *)
    WRITE (6, *)
  END SUBROUTINE displayMatrixInt

  !***************************
  ! Print a matrix of logicals
  !***************************
  SUBROUTINE displayMatrixLog(Mat)
    LOGICAL, DIMENSION(:, :), INTENT(IN) :: Mat
    INTEGER :: di, dj, i, j

    di = SIZE(Mat, dim=1)
    dj = SIZE(Mat, dim=2)

    DO i = 1, di
      DO j = 1, dj
        IF (Mat(i, j)) THEN
          WRITE (6, '(A)', advance="no") "1 "
        ELSE
          WRITE (6, '(A)', advance="no") "0 "
        END IF
      END DO  ! end j
      WRITE (6, *)
    END DO   ! end i
    WRITE (6, *)
    WRITE (6, *)
  END SUBROUTINE displayMatrixLog

  !***************************
  ! Print a matrix of reals
  !***************************
  SUBROUTINE displayMatrix(Mat)
    REAL, DIMENSION(:, :), INTENT(IN) :: Mat
    INTEGER :: di, dj, i, j

    di = SIZE(Mat, dim=1)
    dj = SIZE(Mat, dim=2)

    DO i = 1, di
      DO j = 1, dj
        !           WRITE(6,'(ES25.15)',advance="no") Mat(i,j)
        !            WRITE(6,'(ES25.15)',advance="no") Mat(i,j)
        WRITE (6, '(ES15.5)', advance="no") Mat(i, j)
      END DO  ! end j
      WRITE (6, *)
    END DO   ! end i
    WRITE (6, *)
    WRITE (6, *)
  END SUBROUTINE displayMatrix

  !*******************************
  ! Print an array 3D of integers
  !*******************************
  SUBROUTINE displayArrayInt(Array)
    INTEGER, DIMENSION(:, :, :), INTENT(IN) :: Array
    INTEGER :: di, dj, dk, i, j, k

    di = SIZE(Array, 1)
    dj = SIZE(Array, 2)
    dk = SIZE(Array, 3)

    DO k = 1, dk
      DO i = 1, di
        DO j = 1, dj
          print *, Array(i, j, k)
        END DO  ! end j
        WRITE (6, *)
      END DO   ! end i
      WRITE (6, *) " "
    END DO     ! end k
  END SUBROUTINE displayArrayInt

  !***************************
  ! Print an array 3D of reals
  !***************************
  SUBROUTINE displayArray(Array)
    REAL, DIMENSION(:, :, :), INTENT(IN) :: Array
    INTEGER :: di, dj, dk, i, j, k

    di = SIZE(Array, 1)
    dj = SIZE(Array, 2)
    dk = SIZE(Array, 3)

    DO k = 1, dk
      DO i = 1, di
        DO j = 1, dj
          IF (abs(Array(i, j, k)) .ge. 1e-10) THEN
            WRITE (6, '(ES16.8,$)') Array(i, j, k)
          ELSE
            WRITE (6, '(ES16.1,$)') 0.0
          END IF
        END DO  ! end j
        WRITE (6, *)
      END DO   ! end i
      WRITE (6, *) " "
    END DO     ! end k

  END SUBROUTINE displayArray

  !***************************
  ! Save an array of reals to
  ! a .txt file
  !***************************
  SUBROUTINE saveArray(Arr, name)
    REAL, DIMENSION(:, :, :), INTENT(IN) :: Arr
    INTEGER :: di, dj, dk, i, j, k
    CHARACTER(LEN=*) :: name

    di = SIZE(Arr, 1)
    dj = SIZE(Arr, 2)
    dk = SIZE(Arr, 3)

    OPEN (33, file=trim(name)//".txt")
    WRITE (33, *) di
    WRITE (33, *) dj
    WRITE (33, *) dk
    DO k = 1, dk
      DO j = 1, dj
        DO i = 1, di
          IF (abs(Arr(i, j, k)) .ge. 1e-10) THEN
            WRITE (33, '(ES25.15E3)') Arr(i, j, k)
          ELSE
            WRITE (33, '(ES16.1)') 0.0
          END IF
        END DO  ! end j
      END DO   ! end i
    END DO     ! end k

    CLOSE (33)
  END SUBROUTINE saveArray

  !***************************
  ! Save a matrix of reals to
  ! a .txt file
  !***************************
  SUBROUTINE saveMatrix(Mat, name)
    REAL, DIMENSION(:, :), INTENT(IN) :: Mat
    INTEGER :: di, dj, dk, i, j, k
    CHARACTER(LEN=*) :: name

    di = SIZE(Mat, 1)
    dj = SIZE(Mat, 2)
    dk = 1

    OPEN (33, file=trim(name)//".txt")
    WRITE (33, *) di
    WRITE (33, *) dj
    WRITE (33, *) dk
    DO j = 1, dj
      DO i = 1, di
        IF (abs(Mat(i, j)) .ge. 1e-10) THEN
          WRITE (33, '(ES25.15E3)') Mat(i, j)
        ELSE
          WRITE (33, '(ES16.1)') 0.0
        END IF
      END DO  ! end j
    END DO   ! end i

    CLOSE (33)
  END SUBROUTINE saveMatrix

  !***************************
  ! Save a vector of reals to
  ! a .txt file
  !***************************
  SUBROUTINE saveVector(Vec, name)
    REAL, DIMENSION(:), INTENT(IN) :: Vec
    INTEGER :: di, dj, dk, i, j, k
    CHARACTER(LEN=*) :: name

    di = SIZE(Vec)
    dj = 1
    dk = 1

    OPEN (33, file=trim(name)//".txt")
    WRITE (33, *) di
    WRITE (33, *) dj
    WRITE (33, *) dk
    DO i = 1, di
      IF (abs(Vec(i)) .ge. 1e-10) THEN
        WRITE (33, '(ES25.15E3)') Vec(i)
      ELSE
        WRITE (33, '(ES16.1)') 0.0
      END IF
    END DO  ! end j

    CLOSE (33)
  END SUBROUTINE saveVector

  SUBROUTINE save_CSR_matrix_txt(Mat, name)
    USE matrices_types
    TYPE(MAT_CSR_TYP), intent(in)   :: Mat
    INTEGER :: i, n, nnz
    CHARACTER(LEN=*) :: name
    n = Mat%n
    nnz = Mat%nnz
    OPEN (33, file=trim(name)//".txt")
    WRITE (33, *) Mat%n
    WRITE (33, *) Mat%nnz
    WRITE (33, *) "vals"
    DO i = 1, nnz
      IF (abs(Mat%vals(i)) .ge. 1e-10) THEN
        WRITE (33, '(ES25.15E3)') Mat%vals(i)
      ELSE
        WRITE (33, '(ES16.1)') 0.0
      END IF
    END DO  ! end j
    WRITE (33, *) "cols"
    DO i = 1, nnz
      WRITE (33, '(I10)') Mat%cols(i)
    END DO
    WRITE (33, *) "rowptr"
    DO i = 1, n + 1
      WRITE (33, '(I10)') Mat%rowptr(i)
    END DO
    WRITE (33, *) "loc2glob"
    DO i = 1, n
      WRITE (33, '(I10)') Mat%loc2glob(i)
    END DO
    CLOSE (33)
  END SUBROUTINE save_CSR_matrix_txt

  SUBROUTINE save_CSR_vector_txt(vec, name)
    USE matrices_types
    TYPE(RHS_TYP), intent(in)   :: vec
    INTEGER :: i, n, nnz
    CHARACTER(LEN=*) :: name
    n = vec%n
    OPEN (33, file=trim(name)//".txt")
    WRITE (33, *) vec%n
    WRITE (33, *) "vals"
    DO i = 1, n
      IF (abs(vec%vals(i)) .ge. 1e-10) THEN
        WRITE (33, '(ES25.15E3)') vec%vals(i)
      ELSE
        WRITE (33, '(ES16.1)') 0.0
      END IF
    END DO
    WRITE (33, *) "loc2glob"
    DO i = 1, n
      WRITE (33, '(I10)') vec%loc2glob(i)
    END DO
    CLOSE (33)
  END SUBROUTINE save_CSR_vector_txt

#ifdef PARALL
  SUBROUTINE syncroprint_vector_int(Vec)
    USE MPI_OMP
    INTEGER, DIMENSION(:), INTENT(IN) :: Vec
    INTEGER :: i, code, pr, aux
    INTEGER, PARAMETER :: etiquette = 1000
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    aux = 1

    IF (MPIvar%glob_id == 0) THEN
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVectorInt(Vec)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    ELSE IF (MPIvar%glob_id == MPIvar%glob_size - 1) THEN
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVectorInt(Vec)
    ELSE
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVectorInt(Vec)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    END IF
    WRITE (6, *) "DONE WRITING FOR Process ", MPIvar%glob_id
    WRITE (6, *) " "
    CALL mpi_barrier(MPI_COMM_WORLD, aux)
  END SUBROUTINE syncroprint_vector_int

  SUBROUTINE syncroprint_vector(Vec)
    USE MPI_OMP
    REAL, DIMENSION(:), INTENT(IN) :: Vec
    INTEGER :: i, code, pr, aux
    INTEGER, PARAMETER :: etiquette = 1000
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    aux = 1

    IF (MPIvar%glob_id == 0) THEN
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVector(Vec)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    ELSE IF (MPIvar%glob_id == MPIvar%glob_size - 1) THEN
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVector(Vec)
    ELSE
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayVector(Vec)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    END IF
    WRITE (6, *) "DONE WRITING FOR Process ", MPIvar%glob_id

    CALL mpi_barrier(MPI_COMM_WORLD, aux)
  END SUBROUTINE syncroprint_vector

  SUBROUTINE syncroprint_matrix(Mat)
    USE MPI_OMP
    REAL, DIMENSION(:, :), INTENT(IN) :: Mat
    INTEGER :: i, code, pr, aux
    INTEGER, PARAMETER :: etiquette = 1000
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    aux = 1

    IF (MPIvar%glob_id == 0) THEN
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayMatrix(Mat)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    ELSE IF (MPIvar%glob_id == MPIvar%glob_size - 1) THEN
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayMatrix(Mat)
    ELSE
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayMatrix(Mat)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    END IF
    WRITE (6, *) "DONE WRITING FOR Process ", MPIvar%glob_id

    CALL mpi_barrier(MPI_COMM_WORLD, aux)

  END SUBROUTINE syncroprint_matrix

  SUBROUTINE syncroprint_array(Mat)
    USE MPI_OMP
    REAL, DIMENSION(:, :, :), INTENT(IN) :: Mat
    INTEGER :: i, code, pr, aux
    INTEGER, PARAMETER :: etiquette = 1000
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: statut

    aux = 1

    IF (MPIvar%glob_id == 0) THEN
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayArray(Mat)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    ELSE IF (MPIvar%glob_id == MPIvar%glob_size - 1) THEN
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayArray(Mat)
    ELSE
      CALL MPI_RECV(aux, 1, MPI_INTEGER, MPIvar%glob_id - 1, etiquette, MPI_COMM_WORLD, statut, code)
      WRITE (6, *) "Process ", MPIvar%glob_id
      CALL displayArray(Mat)
      CALL MPI_SEND(aux, 1, MPI_INTEGER, MPIvar%glob_id + 1, etiquette, MPI_COMM_WORLD, code)
    END IF
    WRITE (6, *) "DONE WRITING FOR Process ", MPIvar%glob_id

    CALL mpi_barrier(MPI_COMM_WORLD, aux)

  END SUBROUTINE syncroprint_array
#endif
END MODULE Printutils
