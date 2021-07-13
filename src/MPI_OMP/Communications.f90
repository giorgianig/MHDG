!*****************************************
! project: MHDG
! file: Communications.f90
! date: 10/03/2017
! Exchange solution between ghost faces
!*****************************************
MODULE Communications
  USE MPI_OMP
  USE globals
  USE PrintUtils

CONTAINS

  SUBROUTINE init_Com()
    integer, parameter                   :: etq = 100
    integer, dimension(MPI_STATUS_SIZE)  :: stat
    integer                             :: code, req, ierr
    integer                             :: npro, rbuf, lpro, i, nf2sd, ct, ncon, ne2sd
    integer                             :: psd, prv, ifa, fcount, ecount
    integer, dimension(MPIvar%glob_size) :: proext, proint, proc
    integer, allocatable                :: faces_ext(:), faces_int(:)!,connpro(:)
    integer, allocatable                :: elems_ext(:), elems_int(:)

    ! Number of processes
    npro = MPIvar%glob_size

    ! Local process
    lpro = MPIvar%glob_id

    ! faces and processes to receive
    ALLOCATE (Mesh%pr2rv(Mesh%nghostfaces))
    ALLOCATE (Mesh%fc2rv(Mesh%nghostfaces))
    Mesh%pr2rv = Mesh%ghostpro
    Mesh%fc2rv = 0
    ct = 1
    DO i = 1, Mesh%Nfaces
      IF (Mesh%ghostfaces(i) .eq. 1) THEN
        Mesh%fc2rv(ct) = i
        ct = ct + 1
      END IF
    END DO
    ! ***********************************************************
    ! proext contains the number of ghost faces that each process
    ! compute for the local process, proint contains the number of
    ! ghost faces that the local process compute for the other
    ! processes
    ! ***********************************************************
    proext = 0
    proint = 0
    DO i = 1, Mesh%nghostfaces
      IF (Mesh%ghostpro(i).eq.0) CYCLE ! Add by Benjamin
      proext(Mesh%ghostpro(i)) = proext(Mesh%ghostpro(i)) + 1
    END DO
    DO i = 1, npro
      rbuf = 0
      CALL MPI_Scatter(proext, 1, MPI_INTEGER, rbuf, 1, MPI_INTEGER, i - 1, MPI_COMM_WORLD, code)
      proint(i) = rbuf
    END DO

    ! Number of faces to send
    nf2sd = 0
    DO i = 1, npro
      nf2sd = nf2sd + proint(i)
    END DO
    ALLOCATE (Mesh%fc2sd(nf2sd))
    ALLOCATE (Mesh%pr2sd(nf2sd))

    ! Connected processes
    ! ncon = 0
    ! DO i=1,npro
    !    IF (proext(i).eq.0 .and. proint(i).eq.0) CYCLE
    !    ncon=ncon+1
    ! END DO
    ! ALLOCATE(connpro(ncon))
    ! connpro=0
    ! proc = 0
    ! ct = 0
    ! DO i=1,npro
    !   IF (proext(i).eq.0 .and. proint(i).eq.0) CYCLE
    !   IF (proc(i).eq.0) THEN
    !      ct = ct+1
    !      connpro(ct)=i
    !      proc(i) = 1
    !   END IF
    ! END DO

    ! call syncroprint_vector_int(Mesh%connpro)
    ! stop

    ! Communicate to each process the local number of the faces he needs to send,
    ! and to whom
    fcount = 1
    DO i = 1, npro

      IF (MPIvar%glob_id .eq. (i - 1)) CYCLE
      ! Process to send
      IF (proext(i) .eq. 0) THEN
        psd = MPI_PROC_NULL
      ELSE
        psd = i - 1
      END IF
      ! Process to receive
      IF (proint(i) .eq. 0) THEN
        prv = MPI_PROC_NULL
      ELSE
        prv = i - 1
      END IF
      ALLOCATE (faces_ext(proext(i)))
      ALLOCATE (faces_int(proint(i)))
      faces_ext = 0
      faces_int = 0

      ! prepare the vector with the index of the faces the local
      ! process need to receive from a connected process
      ct = 1
      DO ifa = 1, Mesh%nghostfaces
        IF (Mesh%ghostpro(ifa) .eq. i) THEN
          faces_ext(ct) = Mesh%ghostLoc(ifa)
          ct = ct + 1
        END IF
      END DO

      CALL MPI_SEND(faces_ext, proext(i), MPI_INTEGER, psd, etq, MPI_COMM_WORLD, code)
      CALL MPI_RECV(faces_int, proint(i), MPI_INTEGER, prv, etq, MPI_COMM_WORLD, stat, code)

      Mesh%fc2sd(fcount:fcount + proint(i) - 1) = faces_int
      Mesh%pr2sd(fcount:fcount + proint(i) - 1) = i
      fcount = fcount + proint(i)

      DEALLOCATE (faces_ext, faces_int)
    END DO

#ifdef TOR3D
    ! elements and processes to receive
    IF (Mesh%nghostelems .gt. 0) THEN
      ALLOCATE (Mesh%pe2rv(Mesh%nghostelems))
      ALLOCATE (Mesh%el2rv(Mesh%nghostelems))
      Mesh%pe2rv = Mesh%ghelspro
      Mesh%el2rv = 0
      ct = 1
      DO i = 1, Mesh%Nelems
        IF (Mesh%ghostelems(i) .eq. 1) THEN
          Mesh%el2rv(ct) = i
          ct = ct + 1
        END IF
      END DO
    ENDIF
    ! ***********************************************************
    ! proext contains the number of ghost elements that each process
    ! compute for the local process, proint contains the number of
    ! ghost elements that the local process compute for the other
    ! processes
    ! ***********************************************************
    proext = 0
    proint = 0
    DO i = 1, Mesh%nghostelems
      proext(Mesh%ghelspro(i)) = proext(Mesh%ghelspro(i)) + 1
    END DO

    DO i = 1, npro
      rbuf = 0
      CALL MPI_Scatter(proext, 1, MPI_INTEGER, rbuf, 1, MPI_INTEGER, i - 1, MPI_COMM_WORLD, code)
      proint(i) = rbuf
    END DO

    ! Number of elements to send
    ne2sd = 0
    DO i = 1, npro
      ne2sd = ne2sd + proint(i)
    END DO
    ALLOCATE (Mesh%el2sd(ne2sd))
    ALLOCATE (Mesh%pe2sd(ne2sd))

    ! Communicate to each process the local number of the faces he needs to send,
    ! and to whom
    ecount = 1
    DO i = 1, npro

      IF (MPIvar%glob_id .eq. (i - 1)) CYCLE
      ! Process to send
      IF (proext(i) .eq. 0) THEN
        psd = MPI_PROC_NULL
      ELSE
        psd = i - 1
      END IF
      ! Process to receive
      IF (proint(i) .eq. 0) THEN
        prv = MPI_PROC_NULL
      ELSE
        prv = i - 1
      END IF
      ALLOCATE (elems_ext(proext(i)))
      ALLOCATE (elems_int(proint(i)))
      elems_ext = 0
      elems_int = 0

      ! prepare the vector with the index of the faces the local
      ! process need to receive from a connected process
      ct = 1
      DO ifa = 1, Mesh%nghostelems
        IF (Mesh%ghelspro(ifa) .eq. i) THEN
          elems_ext(ct) = Mesh%ghelsLoc(ifa)
          ct = ct + 1
        END IF
      END DO

      CALL MPI_SEND(elems_ext, proext(i), MPI_INTEGER, psd, etq, MPI_COMM_WORLD, code)
      CALL MPI_RECV(elems_int, proint(i), MPI_INTEGER, prv, etq, MPI_COMM_WORLD, stat, code)

      Mesh%el2sd(ecount:ecount + proint(i) - 1) = elems_int
      Mesh%pe2sd(ecount:ecount + proint(i) - 1) = i
      ecount = ecount + proint(i)

      DEALLOCATE (elems_ext, elems_int)
    END DO
#endif
    ! call syncroprint_vector_int(Mesh%fc2sd)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! call syncroprint_vector_int(Mesh%pr2sd)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! call syncroprint_vector_int(Mesh%fc2rv)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! stop
  END SUBROUTINE init_com

#ifdef TOR3D
  SUBROUTINE exchangeSol()
    integer, parameter  :: etq = 100
    integer            :: Neq, Nfl, Np2d, i, j, Fi, itor, N2d, Nfdir
    integer            :: indt(refElTor%Nfl*phys%Neq)
    integer            :: indp(refElPol%Nnodes2D*phys%Neq)
    integer            :: nf2sd, ne2sd
    real*8, allocatable :: buffrv(:, :), buffsd(:, :)
    real*8             :: fbufsd(refElTor%Nfl*phys%Neq)
    real*8             :: fbufrv(refElTor%Nfl*phys%Neq)
    integer            :: code, ierr
    integer, allocatable:: req(:), stat(:, :)
    integer            :: dd, delta
    integer            :: prtorsd, prtorrv, iel, ifa
    integer :: Ntorloc,auxel(refElPol%Nnodes2D*phys%neq),auxfl(phys%Neq*refElTor%Nfl)
    ! integer            :: perm(1:refElPol%Nfacenodes*phys%Neq)

#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
      !       ntorass = ntorloc-1
    ELSE
      ntorloc = numer%ntor
      !       ntorass = ntorloc
    ENDIF
#else
    ntorloc = numer%ntor
    ntorass = ntorloc
#endif
    Neq = phys%Neq
    Nfl = refElTor%Nfl
    Np2D = refElPol%Nnodes2D
    N2d = Mesh%Nelems
    Nfdir = Mesh%Ndir
    nf2sd = size(Mesh%fc2sd)
    ne2sd = size(Mesh%el2sd)

    !************************************************************************************************
    !
    !       COMMUNICATIONS FOR POLOIDAL DISTRIBUTION OF THE MESH
    !
    !************************************************************************************************
    !*************************************************
    !  Communication for poloidal faces
    !*************************************************
    ALLOCATE (req(ne2sd + Mesh%nghostelems))
    ALLOCATE (stat(MPI_STATUS_SIZE, ne2sd + Mesh%nghostelems))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Np2D*phys%Neq, Mesh%nghostelems))
    ALLOCATE (buffsd(Np2D*phys%Neq, ne2sd))

    auxel = (/(j, j=0, Np2D*Neq - 1)/)
    DO itor = 1, ntorloc
      buffrv = 0.
      buffsd = 0.
      dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq

      ! Filling the send buffer
      DO i = 1, ne2sd
        Fi = Mesh%el2sd(i)
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        buffsd(:, i) = sol%u_tilde(indp)
      END DO

      ! Receiving
      DO i = 1, Mesh%nghostelems
        CALL MPI_IRECV(buffrv(:, i), Neq*Np2D, MPI_DOUBLE_PRECISION, Mesh%pe2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
      END DO

      ! Sending
      DO i = 1, ne2sd
        CALL MPI_ISEND(buffsd(:,i),Neq*Np2D,MPI_DOUBLE_PRECISION,Mesh%pe2sd(i)-1,etq,&
          &MPI_COMM_WORLD,req(Mesh%nghostelems+i),code)
      END DO

      CALL MPI_WAITALL(size(req), req, stat, code)

      ! Storing at the right place
      DO i = 1, Mesh%nghostelems
        Fi = Mesh%el2rv(i)!<--here Fi is the poloidal element
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        sol%u_tilde(indp) = buffrv(:, i)
      END DO

    END DO
    DEALLOCATE (buffrv, buffsd, req, stat)

    !*************************************************
    !  Communication for toroidal faces
    !*************************************************
    ALLOCATE (req(nf2sd + Mesh%nghostfaces))
    ALLOCATE (stat(MPI_STATUS_SIZE, nf2sd + Mesh%nghostfaces))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Nfl*phys%Neq, Mesh%nghostfaces))
    ALLOCATE (buffsd(Nfl*phys%Neq, nf2sd))

    auxfl = (/(j, j=0, Nfl*Neq - 1)/)
    DO itor = 1, ntorloc
      buffrv = 0.
      buffsd = 0.
      dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq

      ! Filling the send buffer
      DO i = 1, nf2sd
        Fi = Mesh%fc2sd(i)
        delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
        indt = delta + auxfl
        buffsd(:, i) = sol%u_tilde(indt)
      END DO
      ! Receiving
      DO i = 1, Mesh%nghostfaces
        CALL MPI_IRECV(buffrv(:, i), Neq*Nfl, MPI_DOUBLE_PRECISION, Mesh%pr2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
      END DO

      ! Sending
      DO i = 1, nf2sd
        CALL MPI_ISEND(buffsd(:, i), Neq*Nfl, MPI_DOUBLE_PRECISION, Mesh%pr2sd(i) - 1, &
          &etq, MPI_COMM_WORLD, req(Mesh%nghostfaces + i), code)
      END DO

      CALL MPI_WAITALL(size(req), req, stat, code)
      ! Storing at the right place
      DO i = 1, Mesh%nghostfaces
        Fi = Mesh%fc2rv(i)
        delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
        indt = delta + auxfl
        sol%u_tilde(indt) = buffrv(:, i)
      END DO

    END DO
    DEALLOCATE (buffrv, buffsd, req, stat)

    !************************************************************************************************
    !
    !       COMMUNICATIONS FOR TOROIDAL DISTRIBUTION OF THE MESH
    !
    !************************************************************************************************
    IF (MPIvar%ntor .gt. 1) THEN

      !******************************************
      !
      !************ Sending backward ************
      !
      !******************************************
      ! Process to send to
      IF (MPIvar%itor == 1) THEN
        prtorsd = (MPIvar%ntor - 1)*MPIvar%npol + MPIvar%ipol - 1
      ELSE
        prtorsd = (MPIvar%itor - 2)*MPIvar%npol + MPIvar%ipol - 1
      END IF

      ! Process to receive from
      IF (MPIvar%itor == MPIvar%ntor) THEN
        prtorrv = MPIvar%ipol - 1
      ELSE
        prtorrv = MPIvar%itor*MPIvar%npol + MPIvar%ipol - 1
      END IF

      !******************************************
      ! Sending backward toroidal faces
      !******************************************
      ALLOCATE (req(2*(Mesh%Nfaces - Nfdir)))
      ALLOCATE (stat(MPI_STATUS_SIZE, 2*(Mesh%Nfaces - Nfdir)))
      req = mpi_request_null
      ! Allocate buffers
      ALLOCATE (buffrv(Nfl*phys%Neq, (Mesh%Nfaces - Nfdir)))
      ALLOCATE (buffsd(Nfl*phys%Neq, (Mesh%Nfaces - Nfdir)))
      buffrv = 0.
      buffsd = 0.

      dd = 1
      auxfl = (/(j, j=0, Nfl*Neq - 1)/)
      ! Filling the send buffer
      DO Fi = 1, Mesh%Nfaces
        IF (Fi .gt. Mesh%Nintfaces) THEN
          iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
          ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
          IF (Mesh%Fdir(iel, ifa)) CYCLE
        END IF
        delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
        indt = delta + auxfl
        buffsd(:, Fi) = sol%u_tilde(indt)
      END DO

      ! Receiving
      DO Fi = 1, Mesh%Nfaces
        IF (Fi .gt. Mesh%Nintfaces) THEN
          iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
          ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
          IF (Mesh%Fdir(iel, ifa)) CYCLE
        END IF
        CALL MPI_IRECV(buffrv(:, Fi), Neq*Nfl, MPI_DOUBLE_PRECISION, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
      END DO

      ! Sending
      DO Fi = 1, Mesh%Nfaces
        IF (Fi .gt. Mesh%Nintfaces) THEN
          iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
          ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
          IF (Mesh%Fdir(iel, ifa)) CYCLE
        END IF
        CALL MPI_ISEND(buffsd(:, Fi), Neq*Nfl, MPI_DOUBLE_PRECISION, prtorsd, &
          &etq, MPI_COMM_WORLD, req(Mesh%Nfaces - Nfdir + i), code)
      END DO

      CALL MPI_WAITALL(size(req), req, stat, code)

      ! Storing at the right place
      dd = 1 + (ntorloc - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
      DO Fi = 1, Mesh%Nfaces
        IF (Fi .gt. Mesh%Nintfaces) THEN
          iel = Mesh%extfaces(Fi - Mesh%Nintfaces, 1)
          ifa = Mesh%extfaces(Fi - Mesh%Nintfaces, 2)
          IF (Mesh%Fdir(iel, ifa)) CYCLE
        END IF
        delta = dd + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
        indt = delta + auxfl
        sol%u_tilde(indt) = buffrv(:, Fi)
      END DO

      DEALLOCATE (buffrv, buffsd, req, stat)

      !******************************************
      ! Sending backward poloidal faces
      !******************************************
      ALLOCATE (req(2*Mesh%Nelems))
      ALLOCATE (stat(MPI_STATUS_SIZE, 2*Mesh%Nelems))
      req = mpi_request_null

      ! Allocate buffers
      ALLOCATE (buffrv(Np2D*phys%Neq, Mesh%Nelems))
      ALLOCATE (buffsd(Np2D*phys%Neq, Mesh%Nelems))
      buffsd = 0.
      buffrv = 0.

      ! Filling the send buffer
      auxel =  (/(j, j=0, Np2D*Neq - 1)/)
      dd = 1 + (N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
      DO Fi = 1, Mesh%Nelems
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        buffsd(:, Fi) = sol%u_tilde(indp)
      END DO

      ! Receiving
      DO Fi = 1, Mesh%Nelems
        CALL MPI_IRECV(buffrv(:, Fi), Neq*Np2D, MPI_DOUBLE_PRECISION, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
      END DO

      ! Sending
      DO Fi = 1, Mesh%Nelems
        CALL MPI_ISEND(buffsd(:, Fi), Neq*Np2D, MPI_DOUBLE_PRECISION, prtorsd, etq, MPI_COMM_WORLD, req(Mesh%Nelems + Fi), code)
      END DO

      CALL MPI_WAITALL(size(req), req, stat, code)

      ! Storing at the right place
      dd = 1 + ntorloc*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
      DO Fi = 1, Mesh%Nelems
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        sol%u_tilde(indp) = buffrv(:, Fi)
      END DO

      !******************************************
      !
      !************ Sending forward ************
      !
      !******************************************

      ! Process to receive from
      IF (MPIvar%itor == 1) THEN
        prtorrv = (MPIvar%ntor - 1)*MPIvar%npol + MPIvar%ipol - 1
      ELSE
        prtorrv = (MPIvar%itor - 2)*MPIvar%npol + MPIvar%ipol - 1
      END IF

      ! Process to send to
      IF (MPIvar%itor == MPIvar%ntor) THEN
        prtorsd = MPIvar%ipol - 1
      ELSE
        prtorsd = MPIvar%itor*MPIvar%npol + MPIvar%ipol - 1
      END IF
      !******************************************
      ! Sending forward poloidal faces
      !******************************************
      buffsd = 0.
      buffrv = 0.

      ! Filling the send buffer
      dd = 1 + (ntorloc - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)*Neq
      DO Fi = 1, Mesh%Nelems
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        buffsd(:, Fi) = sol%u_tilde(indp)
      END DO

      ! Receiving
      DO Fi = 1, Mesh%Nelems
        CALL MPI_IRECV(buffrv(:, Fi), Neq*Np2D, MPI_DOUBLE_PRECISION, prtorrv, etq, MPI_COMM_WORLD, req(Fi), code)
      END DO

      ! Sending
      DO Fi = 1, Mesh%Nelems
        CALL MPI_ISEND(buffsd(:, Fi), Neq*Np2D, MPI_DOUBLE_PRECISION, prtorsd, etq, MPI_COMM_WORLD, req(Mesh%Nelems + Fi), code)
      END DO

      CALL MPI_WAITALL(size(req), req, stat, code)

      ! Storing at the right place
      dd = 1
      DO Fi = 1, Mesh%Nelems
        delta = dd + (Fi - 1)*Np2D*Neq !<--here Fi is the poloidal element
        indp = delta + auxel
        sol%u_tilde(indp) = buffrv(:, Fi)
      END DO

      DEALLOCATE (buffrv, buffsd, req, stat)
    END IF

  END SUBROUTINE exchangeSol
#else
  SUBROUTINE exchangeSol()
    integer, parameter  :: etq = 100
    integer            :: Neq, Nfp, i, j, Fi
    integer            :: ind(Mesh%Nnodesperface*phys%Neq)
    integer            :: nf2sd
    real*8, allocatable :: buffrv(:, :), buffsd(:, :)
    real*8             :: fbufsd(Mesh%Nnodesperface*phys%Neq)
    real*8             :: fbufrv(Mesh%Nnodesperface*phys%Neq)
    integer            :: code, ierr
    integer, allocatable:: req(:), stat(:, :)
    integer            :: perm(refElPol%Nfacenodes*phys%Neq),aux(phys%Neq*Mesh%Nnodesperface)

    Neq = phys%Neq
    Nfp = Mesh%Nnodesperface
    nf2sd = size(Mesh%fc2sd)

    ! Set permutations for ghostfaces that need to be flipped
    CALL set_permutations(Neq*Nfp, Neq, perm)

    ALLOCATE (req(nf2sd + Mesh%nghostfaces))
    ALLOCATE (stat(MPI_STATUS_SIZE, nf2sd + Mesh%nghostfaces))
    req = mpi_request_null

    ! Allocate buffers
    ALLOCATE (buffrv(Mesh%Nnodesperface*phys%Neq, Mesh%nghostfaces))
    ALLOCATE (buffsd(Mesh%Nnodesperface*phys%Neq, nf2sd))
    buffrv = 0.
    buffsd = 0.
    ! Filling the send buffer
    aux =  (/(j, j=1, Neq*Nfp)/)
    DO i = 1, nf2sd
      Fi = Mesh%fc2sd(i)
      ind = (Fi - 1)*Neq*Nfp + aux
      buffsd(:, i) = sol%u_tilde(ind)
    END DO

    ! Receiving
    DO i = 1, Mesh%nghostfaces
      CALL MPI_IRECV(buffrv(:, i), Neq*Nfp, MPI_DOUBLE_PRECISION, Mesh%pr2rv(i) - 1, etq, MPI_COMM_WORLD, req(i), code)
    END DO

    ! Sending
    DO i = 1, nf2sd
      CALL MPI_ISEND(buffsd(:, i), Neq*Nfp, MPI_DOUBLE_PRECISION, Mesh%pr2sd(i) - 1, &
        &etq, MPI_COMM_WORLD, req(Mesh%nghostfaces + i), code)
    END DO

    CALL MPI_WAITALL(size(req), req, stat, code)





    ! call syncroprint_matrix(buffsd)
    ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ! call syncroprint_matrix(buffrv)
    ! stop

    ! Storing at the right place
    DO i = 1, Mesh%nghostfaces
      Fi = Mesh%fc2rv(i)
      ind = (Fi - 1)*Neq*Nfp + aux
      sol%u_tilde(ind) = buffrv(:, i)
    END DO

    DEALLOCATE (buffrv, buffsd, req, stat)

  CONTAINS
    !*****************************************
    ! Set permutations for flipping faces
    !****************************************
    SUBROUTINE set_permutations(n, m, perm)
      integer, intent(IN)  :: n, m
      integer, intent(OUT) :: perm(:)
      integer              :: i
      integer              :: temp(m, n/m), templr(m, n/m)

      IF (mod(n, m) .ne. 0) then
        WRITE (6, *) 'Error! n must be a multiple of m'
        STOP
      END IF

      templr = 0
      temp = reshape((/(i, i=1, n)/), (/m, n/m/))
      DO i = 1, n/m
        templr(:, i) = temp(:, n/m - i + 1)
      END DO
      perm = reshape(templr, (/n/))
    END SUBROUTINE set_permutations

  END SUBROUTINE exchangeSol
#endif

END MODULE

