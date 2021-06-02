!*****************************************
! project: MHDG
! file: debug.f90
! date: 01/02/2017
! Collection of debugging routines
!*****************************************
MODULE debug

  USE globals
  USE Printutils
  USE in_out

  IMPLICIT NONE

CONTAINS

  SUBROUTINE print_matrices_tex()

    write (6, *) "Saving elemental matrices"
    call saveArray(elMat%iAqq, 'iAqq')
    call saveArray(elMat%Aqu, 'Aqu')
    call saveArray(elMat%Aql, 'Aql')
    call saveArray(elMat%Auq, 'Auq')
    call saveArray(elMat%Auu, 'Auu')
    call saveArray(elMat%Aul, 'Aul')
    call saveArray(elMat%Alq, 'Alq')
    call saveArray(elMat%Alu, 'Alu')
    call saveArray(elMat%All, 'All')
    IF (Mesh%ndir .gt. 0) THEN
      call saveMatrix(elMat%Aql_dir, 'Aql_dir')
      call saveMatrix(elMat%Aul_dir, 'Aul_dir')
    END IF

    call saveArray(elMat%UU, 'UU')
    call saveMatrix(elMat%U0, 'U0')
    call saveArray(elMat%LL, 'LL')
    call saveMatrix(elMat%L0, 'L0')
  END SUBROUTINE print_matrices_tex

  SUBROUTINE print_matrices_hdf5()

    write (6, *) "Saving elemental matrices in hdf5"
    call HDF5_save_array(elMat%iAqq, 'iAqq')
    call HDF5_save_array(elMat%Aqu, 'Aqu')
    call HDF5_save_array(elMat%Aql, 'Aql')
    call HDF5_save_array(elMat%Auq, 'Auq')
    call HDF5_save_array(elMat%Auu, 'Auu')
    call HDF5_save_array(elMat%Aul, 'Aul')
    call HDF5_save_array(elMat%Alq, 'Alq')
    call HDF5_save_array(elMat%Alu, 'Alu')
    call HDF5_save_array(elMat%All, 'All')
    call HDF5_save_matrix(elMat%S, 'S')
    call HDF5_save_matrix(elMat%fh, 'fh')
    IF (Mesh%ndir .gt. 0) THEN
      call HDF5_save_matrix(elMat%Aql_dir, 'Aql_dir')
      call HDF5_save_matrix(elMat%Aul_dir, 'Aul_dir')
    END IF
    call HDF5_save_array(elMat%UU, 'UU')
    call HDF5_save_matrix(elMat%U0, 'U0')
    call HDF5_save_array(elMat%LL, 'LL')
    call HDF5_save_matrix(elMat%L0, 'L0')
  END SUBROUTINE print_matrices_hdf5

  !  SUBROUTINE print_matrices_hdf5_parall()
  !  USE MPI_OMP
  !  CHARACTER(LEN = 60)  :: num
  !
  !  write(6,*) "Saving elemental matrices in hdf5"
  !     call HDF5_save_array(elMat%M, 'M' )
  !     call HDF5_save_array(elMat%Cv,'Cv')
  !     call HDF5_save_array(elMat%H, 'H')
  !     IF (Mesh%ndir.gt.0) THEN
  !        call HDF5_save_matrix(elMat%Hdir, 'Hdir')
  !        call HDF5_save_matrix(elMat%Edir,'Edir')
  !     END IF
  !     call HDF5_save_array(elMat%D, 'D')
  !     call HDF5_save_array(elMat%E,'E')
  !     call HDF5_save_matrix(elMat%S,'S')
  !     call HDF5_save_array(elMat%UU,'UU')
  !     call HDF5_save_matrix(elMat%U0,'U0')
  !     call HDF5_save_array(elMat%Hf,'Hf')
  !     call HDF5_save_array(elMat%Df,'Df')
  !     call HDF5_save_array(elMat%Ef,'Ef')
  !     call HDF5_save_matrix(elMat%fH,'fH')
  !     call HDF5_save_array(elMat%B, 'B')
  !     call HDF5_save_array(elMat%C, 'C')
  !     IF (Mesh%ndir.gt.0) THEN
  !        call HDF5_save_matrix(elMat%Cdir, 'Cdir')
  !     END IF
  !     call HDF5_save_array(elMat%P, 'P')
  !     call HDF5_save_array(elMat%G, 'G')
  !     call HDF5_save_array(elMat%iL, 'iL')
  !     call HDF5_save_array(elMat%Lf, 'Lf')
  !     call HDF5_save_array(elMat%Qf, 'Qf')
  !     call HDF5_save_array(elMat%LL,'LL')
  !     call HDF5_save_matrix(elMat%L0,'L0')
  !#ifdef TEMPERATURE
  !     call HDF5_save_array(elMat%TQ, 'TQ')
  !     call HDF5_save_array(elMat%TQhf, 'TQhf')
  !     call HDF5_save_matrix(elMat%Tfhf, 'Tfhf')
  !     call HDF5_save_matrix(elMat%Tf, 'Tf')
  !     IF (Mesh%ndir.gt.0) THEN
  !        call HDF5_save_matrix(elMat%Thdir, 'Thdir')
  !     END IF
  !#endif
  !  END SUBROUTINE print_matrices_hdf5_parall

  SUBROUTINE print_matrices_singleElem_tex(iel)
    integer, intent(in) :: iel

    write (6, *) "Saving elemental matrices"
    call saveMatrix(elMat%M(:, :, iel), 'M')
    call saveMatrix(elMat%Cv(:, :, iel), 'Cv')
    call saveMatrix(elMat%H(:, :, iel), 'H')
    IF (Mesh%ndir .gt. 0) THEN
      call saveVector(elMat%Hdir(:, iel), 'Hdir')
      call saveVector(elMat%Edir(:, iel), 'Edir')
    END IF
    call saveMatrix(elMat%D(:, :, iel), 'D')
    call saveMatrix(elMat%E(:, :, iel), 'E')
    call saveVector(elMat%S(:, iel), 'S')
    call saveMatrix(elMat%UU(:, :, iel), 'UU')
    call saveVector(elMat%U0(:, iel), 'U0')
    call saveMatrix(elMat%Hf(:, :, iel), 'Hf')
    call saveMatrix(elMat%Df(:, :, iel), 'Df')
    call saveMatrix(elMat%Ef(:, :, iel), 'Ef')
    call saveVector(elMat%fH(:, iel), 'fH')
    call saveMatrix(elMat%B(:, :, iel), 'B')
    call saveMatrix(elMat%C(:, :, iel), 'C')
    IF (Mesh%ndir .gt. 0) THEN
      call saveVector(elMat%Cdir(:, iel), 'Cdir')
    END IF
    call saveMatrix(elMat%P(:, :, iel), 'P')
    call saveMatrix(elMat%G(:, :, iel), 'G')
    call saveMatrix(elMat%iL(:, :, iel), 'iL')
    call saveMatrix(elMat%Lf(:, :, iel), 'Lf')
    call saveMatrix(elMat%Qf(:, :, iel), 'Qf')
    call saveMatrix(elMat%LL(:, :, iel), 'LL')
    call saveVector(elMat%L0(:, iel), 'L0')
#ifdef TEMPERATURE
    call saveMatrix(elMat%TQ(:, :, iel), 'TQ')
    call saveMatrix(elMat%TQhf(:, :, iel), 'TQhf')
    call saveVector(elMat%Tfhf(:, iel), 'Tfhf')
    call saveVector(elMat%Tf(:, iel), 'Tf')
    IF (Mesh%ndir .gt. 0) THEN
      call saveVector(elMat%Thdir(:, iel), 'Thdir')
    END IF
#endif
  END SUBROUTINE print_matrices_singleElem_tex

  SUBROUTINE print_elmat_tex(K, f, iel)
    REAL, INTENT(IN)     :: K(:, :), f(:)
    INTEGER, INTENT(IN)   :: iel
    CHARACTER(LEN=10)  :: num

    WRITE (num, "(I10)") iel
    CALL saveMatrix(K, "K_"//ADJUSTL(num))
    CALL saveVector(f, "f_"//ADJUSTL(num))

  END SUBROUTINE print_elmat_tex

  SUBROUTINE print_elmat_tex_par(K, f, iel)
    USE MPI_OMP
    REAL, INTENT(IN)     :: K(:, :), f(:)
    INTEGER, INTENT(IN)   :: iel
    CHARACTER(LEN=60)  :: num

    WRITE (num, "(A,I0,A,I0,A,I0)") "iel", iel, "_ipr", MPIvar%glob_id, "_npr", MPIvar%glob_size
    CALL saveMatrix(K, "K_"//ADJUSTL(num))
    CALL saveVector(f, "f_"//ADJUSTL(num))

  END SUBROUTINE print_elmat_tex_par

  FUNCTION maxmat(V) RESULT(m)
    REAL     :: V(:, :)
    REAL     :: m
    INTEGER  :: i, j

    m = 0.
    DO i = 1, size(V, 1)
      DO j = 1, size(V, 2)
        m = max(V(i, j), m)
      END DO
    END DO

  END FUNCTION maxmat

  SUBROUTINE findNaNArr(A)
    real*8, intent(in) :: A(:, :, :)
    integer :: i, j, k

    DO k = 1, size(A, 3)
      DO j = 1, size(A, 2)
        DO i = 1, size(A, 1)
          IF (A(i, j, k) .ne. A(i, j, k)) THEN
            WRITE (6, *) "NaN detected in array at pos: i=", i, " j=", j, " k=", k
            stop
          ENDIF
        END DO
      END DO
    END DO
  END SUBROUTINE findNaNArr

  SUBROUTINE findNaNMat(M)
    real*8, intent(in) :: M(:, :)
    integer :: i, j

    DO j = 1, size(M, 2)
      DO i = 1, size(M, 1)
        IF (M(i, j) .ne. M(i, j)) THEN
          WRITE (6, *) "NaN detected in matrix at pos: i=", i, " j=", j
          stop
        ENDIF
      END DO
    END DO
  END SUBROUTINE findNaNMat

  SUBROUTINE findNaNVec(V)
    real*8, intent(in) :: V(:)
    integer :: i

    DO i = 1, size(V, 1)
      IF (V(i) .ne. V(i)) THEN
        WRITE (6, *) "NaN detected in vector at pos: i=", i
        stop
      ENDIF
    END DO
  END SUBROUTINE findNaNVec

END MODULE debug
