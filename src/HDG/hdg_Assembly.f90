!*****************************************
! project: MHDG
! file: hdg_Assembly.f90
! date: 20/12/2016
! Assembly the global matrix
!*****************************************

#ifdef TOR3D

!***********************************************************************
!
!                            VERSION 3D TOROIDAL
!
!***********************************************************************
SUBROUTINE HDG_assembly()
  USE globals
  USE matrices_types
  USE printUtils
  USE debug
  USE LinearAlgebra, ONLY: tensorSumInt
  USE MPI_OMP

  integer                        :: i, j, jj, k, iel, ifa, ifl, nnz, nn, cont, ieloc, ct_sc, iel3, itor, itorg
  integer                        :: Neq, Nf, Nfl, Nfg, Nftot, Nfaces, Nfp, Nel, Nintf, Nextf, Ndirf, Nfunk, blk, blkp, blkt, Ndim, Np
  integer                        :: Np1Dpol, Np1Dtor, Np2d, Npfl, auxnn
  integer                        :: nColsPosPerFaceInt, nColsPosPerFaceExt, nColsPosPerElem
  integer                        :: Fe(refElPol%Nfaces), Fi, Fig
  integer, allocatable            :: indglo(:), indpos(:)
  integer*4                      :: ind_loc(refElPol%Nfaces, refElPol%Nfacenodes*phys%Neq)
  integer                        :: ind_dim(refElPol%Nfaces + 2), ind_sta(refElPol%Nfaces + 2)
  integer                        :: delta(refElPol%Nfaces + 2)
  integer, allocatable            :: linew(:), shift(:)
  integer                        :: pos_start(refElPol%Nfaces + 2)
  integer, pointer, dimension(:)   :: cols, rows, loc2glob
  real*8, pointer, dimension(:)    :: vals, rhsvec
  logical                        :: Fd(refElPol%Nfaces)
  real*8, allocatable             :: Kel(:,:,:), fel(:,:)
  real*8, pointer                 :: Df(:,:), Hf(:,:), Ef(:,:)
  real*8, pointer                 :: UU(:,:), U0(:)
  real*8, allocatable             :: fh(:)
  real*8, allocatable             :: Lf(:,:)
  real*8, pointer                 :: Qf(:,:)
  real*8, pointer                 :: LL(:,:), L0(:)
#ifdef PARALL
  integer                        ::        Nghostf, Nghoste, Fasind(Mesh%Nfaces), Easind(Mesh%Nelems)
  integer                        :: shiftGlob3d(refElPol%Nfaces + 2), rept
#endif
  integer :: ierr
  integer :: Ntorloc, Ntorass

  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1, timing%clock_rate1)
  end if

#ifdef PARALL
  IF (MPIvar%ntor .gt. 1) THEN
    ntorloc = numer%ntor/MPIvar%ntor + 1
    ntorass = ntorloc-1
  ELSE
    ntorloc = numer%ntor
    ntorass = ntorloc
  ENDIF
#else
  ntorloc = numer%ntor
  ntorass = ntorloc
#endif

  Neq = phys%Neq
  Np2d = refElPol%Nnodes2D
  Np = refElPol%Nnodes2D*refElTor%Nnodes1D       ! Number of nodes for each 3D element
  Npfl = refElPol%Nnodes1D*refElTor%Nnodes1D       ! Number of nodes for lateral faces
  Ndim = Mesh%Ndim
  Nf = refElPol%Nfaces
  Nel = Mesh%Nelems
  Nfaces = Mesh%Nfaces
  Nintf = Mesh%Nintfaces                            ! Number of interior faces in the mesh
  Nextf = Mesh%Nextfaces                            ! Number of exterior faces in the mesh
  Ndirf = Mesh%Ndir                                 ! Number of Dirichlet faces in the mesh
  Nfunk = Mesh%ukf                                  ! Number of unknown faces in the mesh (Nfaces-Ndirf)
  Nftot = (Nel + Nfaces)*Ntorass                      ! Number of total faces
  Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Nfl = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfp = Np2D*2 + refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
#ifdef PARALL
  Nghostf = Mesh%nghostfaces
  Nghoste = Mesh%nghostelems
#endif

  ALLOCATE (linew(Nftot))
  ALLOCATE (shift(Nftot))
  ! Compute the dimension and the number of nonzero of the global matrix
  blkp = Np2d*Neq
  blkt = Npfl*Neq

  nn = blkt*Nfunk*Ntorass + blkp*Nel*Ntorass   ! K = nn x nn

  ! Filling ind_dim and ind_sta
  ind_sta(1) = 1
  DO i = 1, Nf + 2
    IF (i == 1 .or. i == Nf + 2) THEN
      blk = blkp
    ELSE
      blk = blkt
    END IF
    ind_dim(i) = blk
    IF (i < Nf + 2) THEN
      ind_sta(i + 1) = ind_sta(i) + ind_dim(i)
    END IF
  END DO

#ifdef PARALL
  nn = nn-blkt*Nghostf*Ntorass
  nn = nn-blkp*Nghoste*Ntorass
  Fasind = 0
  rept = 0
  DO i = 1, Nfaces
    IF (Mesh%ghostfaces(i) .eq. 1) THEN
      rept = rept + 1
      CYCLE
    END IF
    Fasind(i) = i-rept
  END DO
  Easind = 0
  rept = 0
  DO i = 1, Nel
    IF (Mesh%ghostelems(i) .eq. 1) THEN
      rept = rept + 1
      CYCLE
    END IF
    Easind(i) = i-rept
  END DO
#endif


  ! Compute nnz and line width
  CALL computennz()

  ! Deallocate and reallcate pointers of  MatK struct and rhs struct initialize to 0
  CALL init_mat()

  cols => MatK%cols
  rows => MatK%rowptr
  vals => MatK%vals
  loc2glob => MatK%loc2glob
  rhsvec => rhs%vals

  ALLOCATE (indpos(2*blkp + Nf*blkt))
  ALLOCATE (indglo(2*blkp + Nf*blkt))


  ALLOCATE (Kel(Nfp*Neq, Nfp*Neq,Nel*ntorloc))
  ALLOCATE (fel(Nfp*Neq,Nel*ntorloc))
  Kel = 0.
  fel = 0.
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel3)
  !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO itor=1,Ntorloc
    DO iel=1,Nel
      iel3 = (itor-1)*Nel + iel
      CALL computeElementalMatrix(iel3)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  ! Counter
  cont = 0
  !****************************
  ! Loop in toroidal divisions
  !****************************
  DO itor = 1, Ntorloc
#ifdef PARALL
    itorg = itor + (MPIvar%itor-1)*numer%ntor/MPIvar%ntor
    if (itorg == numer%ntor + 1) itorg = 1
#else
    itorg = itor
#endif
    !**************************
    ! Loop in the mesh elements
    !**************************
    DO iel = 1, Nel

      iel3 = (itor-1)*Nel + iel
      Fe = Mesh%F(iel,:)
      Fd = Mesh%Fdir(iel,:)

      !         CALL computeElementalMatrix()

      !#ifdef PARALL
      !    call print_elmat_tex_par(Kel,fel,iel3)
      !#else
      !    call print_elmat_tex(Kel,fel,iel3)
      !#endif

#ifdef PARALL
      delta = 1 + (itor-1)*((Nel-Nghoste)*blkp + (Nfaces-Ndirf-Nghostf)*blkt) +&
        &(/ (Easind(iel)-1)*blkp,(Nel-Nghoste)*blkp+(Fasind(Fe)-1)*blkt,(Nel-Nghoste)*blkp+(Nfaces-Ndirf-Nghostf)*blkt+(Easind(iel)-1)*blkp /)
      shiftGlob3d = 1 + (itorg-1)*(Mesh%Nel_glob*blkp + (Mesh%Nfa_glob-Mesh%Ndir_glob)*blkt) + (/(Mesh%loc2glob_el(iel)-1)*blkp, &
        & Mesh%Nel_glob*blkp+(Mesh%loc2glob_fa(Fe)-1)*blkt, Mesh%Nel_glob*blkp+(Mesh%Nfa_glob-Mesh%Ndir_glob )*blkt+(Mesh%loc2glob_el(iel)-1)*blkp/)
      IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) delta = delta-(Nel-Nghoste)*blkp
#else
      delta = 1+(itor-1)*(Nel*blkp+(Nfaces-Ndirf)*blkt)+(/ (iel-1)*blkp,Nel*blkp+(Fe-1)*blkt,Nel*blkp+(Nfaces-Ndirf)*blkt+(iel-1)*blkp /)
#endif

#ifdef PARALL
      IF (itorg == numer%ntor) THEN
        delta(Nf + 2) = 1 + (Easind(iel)-1)*blkp
        shiftGlob3d(Nf + 2) = 1 + (Mesh%loc2glob_el(iel)-1)*blkp
      ELSE IF (itor == ntorass + 1 .and. MPIvar%itor == MPIvar%ntor) THEN
        !            I enter here only with toroidal parallelization
        delta(1) = 1 + (Easind(iel)-1)*blkp
      ENDIF
#else
      IF (itor == numer%ntor) THEN
        delta(Nf + 2) = 1 + (iel-1)*blkp
      ENDIF
#endif

      DO ifa = 1, Nf + 2
#ifdef PARALL
        indglo(ind_sta(ifa):ind_sta(ifa) + ind_dim(ifa)-1) = shiftGlob3d(ifa) + (/(i, i=0, ind_dim(ifa)-1)/)
#else
        indglo(ind_sta(ifa):ind_sta(ifa) + ind_dim(ifa)-1) = delta(ifa) + (/(i, i=0, ind_dim(ifa)-1)/)
#endif
      END DO
      !***********************
      ! Loop in local faces
      !***********************
      DO ifa = 1, Nf + 2

        IF (ifa == 1) THEN
          !******************************************************
          ! ******** First poloidal face of the element *********
          !******************************************************
#ifdef PARALL
          IF (Mesh%ghostElems(iel) .eq. 1) CYCLE
          IF ((MPIvar%ntor .gt. 1) .and. (itor == 1)) CYCLE ! toroidal ghost face
#endif

          Fig = (itor-1)*(Nel + Nfaces) + iel ! Global numbering of the face
#ifdef PARALL
          IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) Fig = Fig-Nel
#endif
#ifdef PARALL
          ! Case of extra toroidal element for parallel purposes (only toroidal parallelization)
          IF (MPIvar%itor == MPIvar%ntor .and. itor == ntorass + 1) Fig = iel
#endif

          blk = blkp                        ! Number of nodes of the face
          CALL getindpospol(indpos)

        ELSEIF (ifa == Nf + 2) THEN
          !******************************************************
          ! ********* Second poloidal face of the element *******
          !******************************************************
#ifdef PARALL
          IF (Mesh%ghostElems(iel) .eq. 1) CYCLE
          IF (itor == ntorass + 1) CYCLE
#endif

          Fig = itor*(Nel + Nfaces) + iel ! Global numbering of the face
#ifdef PARALL
          IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) Fig = Fig-Nel
#endif
#ifdef PARALL
          ! TODO: Fig è sempre iel????
          IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor == MPIvar%ntor .and. itor == ntorass) THEN
            ! Toroidal parallelization
            Fig = iel
          ELSE IF (MPIvar%ntor .eq. 1 .and. itor == ntorass) THEN
            ! No toroidal parallelization
            Fig = iel
          END IF
#else
          IF (itor == numer%ntor) THEN
            Fig = iel
          END IF

#endif
          blk = blkp                       ! Number of nodes of the face
          CALL getindpospol(indpos)

        ELSE
#ifdef PARALL
          IF (itor == ntorass + 1) CYCLE
#endif
          !******************************************************
          ! ****************** Toroidal face ********************
          !******************************************************
          Fi = Fe(ifa-1)                      ! Face number in the 2D mesh
          Fig = (itor-1)*(Nel + Nfaces) + Nel + Fi ! Global numbering of the face
#ifdef PARALL
          IF (MPIvar%ntor .gt. 1 .and. MPIvar%itor .ne. MPIvar%ntor) Fig = Fig-Nel
#endif
          blk = blkt                          ! Number of nodes of the face
#ifdef PARALL
          IF (Mesh%ghostFaces(Fi) .eq. 1) CYCLE
#endif
          IF (Fi .le. Nintf) THEN
            !**************************
            ! Interior face
            !**************************
            CALL getindposint(indpos)

          ELSE
            !**************************
            ! Exterior face
            !**************************
            IF (Mesh%Fdir(iel, ifa-1)) CYCLE

            CALL getindposext(indpos)

          END IF
        END IF

        CALL fill_cols_vals_rowptr_loc2glob()
      END DO

    END DO ! Loop in elements
  END DO ! Loop in toroidal divisions

  rows(nn + 1) = nnz + 1
  rhs%loc2glob = loc2glob
  DEALLOCATE (Kel, fel)
  DEALLOCATE (linew, shift, indpos, indglo)
  NULLIFY (cols, rows, vals, loc2glob, rhsvec)


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1, timing%clock_rate1)
    timing%runtass = timing%runtass + (timing%cke1-timing%cks1)/real(timing%clock_rate1)
    timing%cputass = timing%cputass + timing%tpe1-timing%tps1
  end if


CONTAINS

  !**************************************************
  ! Compute the nnz of the matrix
  !**************************************************
  SUBROUTINE computennz()
    integer::shift_inc

    nnz = 0
    linew = 0
    shift = 0
    Fig = 0
    shift_inc = 0
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      IF (MPIvar%itor .eq. MPIvar%ntor) THEN
        DO iel = 1, Nel
          ! Increment face count
          Fig = Fig + 1
          IF (Mesh%ghostelems(iel) .eq. 1) THEN
            IF (Fig .gt. 1) THEN
              shift(Fig) = shift(Fig-1)
              linew(Fig) = linew(Fig-1)
            ENDIF
            CYCLE
          END IF
          IF (numer%ntor == 2) THEN
            nnz = nnz + 2*blkp**2
            linew(Fig) = linew(Fig) + 2*blkp
          ELSE
            nnz = nnz + 3*blkp**2
            linew(Fig) = linew(Fig) + 3*blkp
          ENDIF
          ! add nnz for connections with lateral faces
          Fe = Mesh%F(iel,:)
          Fd = Mesh%Fdir(iel,:)

          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + 2*blkp*blkt
            linew(Fig) = linew(Fig) + 2*blkt
          END DO
          IF (Fig .gt. 1) THEN
            shift(Fig) = shift(Fig-1) + shift_inc
          ENDIF
          shift_inc = linew(Fig)*blkp
        END DO
      END IF
    END IF
#endif
    ! ***************************
    ! Loop in toroidal divisions
    !****************************
    DO itor = 1, Ntorass
      ! ************************************************
      ! Loop in 2D elements to count the poloidal faces
      ! ************************************************
      DO iel = 1, Nel
#ifdef PARALL
        IF ((MPIvar%ntor .gt. 1) .and. (itor == 1)) CYCLE
#endif
        ! Increment face count
        Fig = Fig + 1
#ifdef PARALL
        IF (Mesh%ghostelems(iel) .eq. 1) THEN
          IF (Fig .gt. 1) THEN
            shift(Fig) = shift(Fig-1)
            linew(Fig) = linew(Fig-1)
          END IF
          CYCLE
        END IF
#endif
        ! count nnz for this poloidal face
        nnz = nnz + blkp**2
        linew(Fig) = linew(Fig) + blkp

        ! add nnz for other poloidal faces connected
        IF (Ntorloc > 1) THEN                  ! if only 1 toroidal division, nothing to add
          IF (Ntorloc == 2) THEN
            nnz = nnz + blkp**2         !  only 1 connection (with 2nd set)
            linew(Fig) = linew(Fig) + blkp
          ELSE
            nnz = nnz + 2*blkp**2       ! general case, 2 connections
            linew(Fig) = linew(Fig) + 2*blkp
          END IF
        END IF

        ! add nnz for connections with lateral faces
        Fe = Mesh%F(iel,:)
        Fd = Mesh%Fdir(iel,:)
        IF (Ntorloc > 1) THEN
          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + 2*blkp*blkt
            linew(Fig) = linew(Fig) + 2*blkt
          END DO
        ELSE
          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + blkp*blkt
            linew(Fig) = linew(Fig) + blkt
          END DO
        END IF

        IF (Fig .gt. 1) THEN
          shift(Fig) = shift(Fig-1) + shift_inc
        END IF
        shift_inc = linew(Fig)*blkp
      END DO
      ! ************************
      ! Loop in lateral faces
      ! ***********************
      DO Fi = 1, Nfaces
        Fig = Fig + 1
#ifdef PARALL
        IF (Mesh%ghostfaces(Fi) .eq. 1) THEN
          IF (Fig .gt. 1) THEN
            shift(Fig) = shift(Fig-1)
            linew(Fig) = linew(Fig-1)
          END IF
          CYCLE
        END IF
#endif
        IF (Fi .le. Nintf) THEN
          !**************************
          ! Interior face
          !**************************
          ! Connection with poloidal faces
          nnz = nnz + 2*blkp*blkt
          linew(Fig) = linew(Fig) + 2*blkp

          ! add nnz for other poloidal faces connected
          IF (Ntorloc > 1) THEN                  ! if only 1 toroidal division, nothing to add
            nnz = nnz + 2*blkp*blkt
            linew(Fig) = linew(Fig) + 2*blkp
          END IF

          ! add nnz for lateral faces connected
          DO ieloc = 1, 2
            iel = Mesh%intfaces(Fi, (ieloc-1)*2 + 1)
            Fe = Mesh%F(iel,:)
            Fd = Mesh%Fdir(iel,:)
            DO ifl = 1, Nf
              IF (Fd(ifl)) CYCLE
              IF (Fe(ifl) == Fi .and. ieloc == 2) CYCLE
              nnz = nnz + blkt**2
              linew(Fig) = linew(Fig) + blkt
            END DO
          END DO
        ELSE
          !**************************
          ! Exterior face
          !**************************
          iel = Mesh%extfaces(Fi-Nintf, 1)
          ifa = Mesh%extfaces(Fi-Nintf, 2)
          Fd = Mesh%Fdir(iel,:)

          ! Skip Dirichlet faces
          IF (Fd(ifa)) THEN
            IF (Fig .gt. 1) THEN
              shift(Fig) = shift(Fig-1)
              linew(Fig) = linew(Fig-1)
            END IF
            CYCLE
          END IF

          ! Connection with poloidal faces
          nnz = nnz + blkp*blkt
          linew(Fig) = linew(Fig) + blkp

          ! add nnz for other poloidal faces connected
          IF (Ntorloc > 1) THEN                  ! if only 1 toroidal division, nothing to add
            nnz = nnz + blkp*blkt
            linew(Fig) = linew(Fig) + blkp
          END IF

          ! Add nnz for lateral faces connected
          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + blkt**2
            linew(Fig) = linew(Fig) + blkt
          END DO

        END IF

        IF (Fig .gt. 1) THEN
          shift(Fig) = shift(Fig-1) + shift_inc
        END IF
        shift_inc = linew(Fig)*blkt
      END DO
    END DO
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      IF (MPIvar%itor .ne. MPIvar%ntor) THEN
        DO iel = 1, Nel
          ! Increment face count
          Fig = Fig + 1
          IF (Mesh%ghostelems(iel) .eq. 1) THEN
            shift(Fig) = shift(Fig-1)
            linew(Fig) = linew(Fig-1)
            CYCLE
          END IF
          IF (numer%ntor == 2) THEN
            nnz = nnz + 2*blkp**2
            linew(Fig) = linew(Fig) + 2*blkp
          ELSE
            nnz = nnz + 3*blkp**2
            linew(Fig) = linew(Fig) + 3*blkp
          ENDIF
          ! add nnz for connections with lateral faces
          Fe = Mesh%F(iel,:)
          Fd = Mesh%Fdir(iel,:)

          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + 2*blkp*blkt
            linew(Fig) = linew(Fig) + 2*blkt
          END DO
          shift(Fig) = shift(Fig-1) + shift_inc
          shift_inc = linew(Fig)*blkp
        END DO
      END IF
    END IF
#endif
  END SUBROUTINE computennz

  !**************************************************
  ! Generate the elemental matrix
  !**************************************************
  SUBROUTINE computeElementalMatrix(iel3)
    USE LinearAlgebra, ONLY: mymatmul
    integer :: iel3
    real*8, pointer                 :: A_lq(:,:), A_lu(:,:), A_ll(:,:)
    real*8, pointer                 :: LL(:,:), L0(:), UU(:,:), U0(:), f(:)

    real*8, allocatable :: auxAlqLL(:,:),auxAluUU(:,:)


    LL => elMat%LL(:,:, iel3)
    L0 => elMat%L0(:, iel3)
    UU => elMat%UU(:,:, iel3)
    U0 => elMat%U0(:, iel3)

    A_lq => elMat%Alq(:,:, iel3)
    A_lu => elMat%Alu(:,:, iel3)
    A_ll => elMat%All(:,:, iel3)
    f => elMat%fh(:, iel3)

    ! Elemental matrix and elemental RHS

    allocate(auxAlqLL(size(A_lq,1),size(LL,2)))
    allocate(auxAluUU(size(A_lu,1),size(UU,2)))
    call mymatmul(A_lq,LL,auxAlqLL)
    call mymatmul(A_lu,UU,auxAluUU)
    Kel(:,:,iel3) = auxAlqLL+auxAluUU+A_ll
    deallocate(auxAlqLL,auxAluUU)
    !      Kel(:,:,iel3) = matmul(A_lq, LL) + matmul(A_lu, UU) + A_ll




    !write(6,*) "TEST 1"
    !write(6,*) "size(A_lq)",size(A_lq,1),size(A_lq,2)
    !write(6,*) "size(LL)",size(LL,1),size(LL,2)
    !write(6,*) "size(Kel)",size(Kel,1),size(Kel,2)
    !Kel =  matmul(A_lq, LL)
    !write(6,*) "TEST 2"
    !Kel = Kel + matmul(A_lu, UU)
    !write(6,*) "TEST 3"
    !Kel = Kel + A_ll
    !write(6,*) "TEST 4"
    fel(:,iel3) = -(matmul(A_lu, U0) + matmul(A_lq, L0)) + f

    !call HDF5_save_vector(fel,"f1")
    !stop
    NULLIFY (LL, L0, UU, U0)
    NULLIFY (A_lq, A_lu, A_ll)
  END SUBROUTINE computeElementalMatrix
  !                                        SUBROUTINE computeElementalMatrix()
  !                                        Df => elMat%Df(:,:,iel3)
  !     UU => elMat%UU(:,:,iel3)
  !     Hf => elMat%Hf(:,:,iel3)
  !     Ef => elMat%Ef(:,:,iel3)
  !     U0 => elMat%U0(:,iel3)
  !     ALLOCATE(Lf(Neq*Nfp,Neq*Ndim*Np))
  !     ALLOCATE(fh(Neq*Nfp))
  !     LL => elMat%LL(:,:,iel3)
  !     Lf = elMat%Lf(:,:,iel3)
  !     fh = elMat%fh(:,iel3)
  !#ifdef TEMPERATURE
  !     Lf = Lf+elMat%TQhf(:,:,iel3)
  !     fh = fh-elMat%Tfhf(:,iel3)
  !#endif
  !     Qf => elMat%Qf(:,:,iel3)
  !     L0 => elMat%L0(:,iel3)
  !     IF (switch%shockcp.gt.0) THEN
  !                             IF (Mesh%flag_elems_sc(iel3).ne.0) THEN
  !                                Lf = Lf+elMat%Lf_sc(:,:,Mesh%flag_elems_sc(iel3))
  !                             END IF
  !     END IF
  !     ! Faces numbering and Dirichlet booleans
  !     Fe = Mesh%F(iel,:)
  !     Fd = Mesh%Fdir(iel,:)
  !
  !     ! Elemental matrix and elemental RHS
  !     Kel =  matmul(Df,UU) + Hf-Ef
  !     fel = -matmul(Df,U0) + fh
  !     Kel = Kel-matmul((Lf-Qf),LL)
  !     fel = fel + matmul((Lf-Qf),L0)
  !
  !     DEALLOCATE(Lf,fh)
  !     NULLIFY(LL,Qf,L0)
  !     NULLIFY(Df,UU,Hf,Ef,U0)
  !                                        END SUBROUTINE computeElementalMatrix

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for poloidal faces
  !**************************************************
  SUBROUTINE getindpospol(indpos)
    integer :: i, j, indpos(:)
    integer :: pos(1:Nf)
    integer ::  qq, blkl
    logical :: skip_first

    ! Filling pos: needed to determine the relative position of the
    ! toroidal faces in the assembling line of the matrix for the
    ! face Fig
    pos = 1
    DO i = 1, Nf
      IF (Fd(i)) CYCLE
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
    END DO

    ! qq gives here the number of lateral faces assebled for this element
    ! (take out Dirichlet faces)
    qq = maxval(pos)

    ! Filling pos_start: needed to determine the starting position of assembling
    ! the coefficients in the cols and vals vectors
    pos_start = 0
    pos_start(2:Nf + 1) = (pos-1)*blkt
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      !***************************************************************************************
      !
      !                  CASE TOROIDAL +  POLOIDAL PARALLELLIZATION
      !
      !***************************************************************************************
      IF (numer%ntor == 2) THEN
        ! ONLY 2 toroidal divisions. MPIvar%ntor must be 2
        IF (MPIvar%itor == 1) THEN

          IF (itor == 1) THEN
            ! first set of poloidal faces
            ! GLOB. EQ.    0 x x x 0     x x x
            !              O x x x O  [O x x x O]
            pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
            pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
          ELSE IF (itor == 2) THEN
            ! case last set with only 2 toroidal divisions
            ! GLOB. EQ.    0 x x x 0     x x x
            !             [O x x x O]  O x x x O
            pos_start(1) = blkp + qq*blkt
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
            pos_start(Nf + 2) = 0
          END IF

        ELSE IF (MPIvar%itor == 2) THEN
          IF (itor == 1) THEN
            ! GLOB. EQ.    0 x x x 0     x x x
            !             [O x x x O]  O x x x O
            pos_start(1) = blkp + qq*blkt
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
            pos_start(Nf + 2) = 0
          ELSE IF (itor == 2) THEN
            ! GLOB. EQ.    0 x x x 0     x x x
            !              O x x x O  [O x x x O]
            pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
            pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
          END IF
        END IF

      ELSE

        IF (itorg == 1) THEN
          ! first set of poloidal faces
          ! GLOB. EQ.    0 x x x O   O x x x
          !              O x x x O  [O x x x O]
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt

        ELSE IF (itorg == numer%ntor) THEN
          ! case last set with more then 2 toroidal divisions
          ! GLOB. EQ.    0 x x x O   O x x x
          !             [O x x x O]  O x x x O
          pos_start = pos_start + 2*blkp + qq*blkt
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = 0
        ELSE IF ((ifa == Nf + 2) .and. itorg == numer%ntor-1) THEN
          ! case second-last set
          ! GLOB. EQ.    O O x x x 0    x x x
          !                O x x x O [O x x x O]
          pos_start = pos_start + blkp
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt

        ELSE
          ! general case
          IF (ifa == Nf + 2) THEN
            ! GLOB. EQ.    O x x x 0    x x x O
            !              O x x x O [O x x x O]
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + blkp
            pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt + blkp
          ELSE
            ! GLOB. EQ.    O x x x    0 x x x O
            !             [O x x x O] O x x x O
            pos_start(1) = blkp + qq*blkt
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
            pos_start(Nf + 2) = 2*blkp + 2*qq*blkt
          END IF
        END IF
      END IF

    ELSE

      !***************************************************************************************
      !
      !                  CASE ONLY POLOIDAL PARALLELLIZATION
      !
      !***************************************************************************************

      IF (Ntorloc == 1) THEN
        ! ****************************
        ! Case only 1 toroidal division
        ! ****************************
        pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
        pos_start(Nf + 2) = 0
      ELSE
        ! ****************************************
        ! Case only more than 1 toroidal divisions
        ! ****************************************
        IF (itor == 1) THEN
          ! first set of poloidal faces
          ! GLOB. EQ.    0 x x x O   O x x x
          !              O x x x O  [O x x x O]
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
        ELSE IF (itor == Ntorloc .and. Ntorloc == 2) THEN
          ! case last set with only 2 toroidal divisions
          ! GLOB. EQ.    0 x x x 0     x x x
          !             [O x x x O]  O x x x O
          pos_start(1) = blkp + qq*blkt
          pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
          pos_start(Nf + 2) = 0
        ELSE IF (itor == Ntorloc) THEN
          ! case last set with more then 2 toroidal divisions
          ! GLOB. EQ.    0 x x x O   O x x x
          !             [O x x x O]  O x x x O
          pos_start = pos_start + 2*blkp + qq*blkt
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = 0
        ELSE IF ((ifa == Nf + 2) .and. itor == Ntorloc-1) THEN
          ! case second-last set
          ! GLOB. EQ.    O O x x x 0    x x x
          !                O x x x O [O x x x O]
          pos_start = pos_start + blkp
          pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
          pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
        ELSE
          ! general case
          IF (ifa == Nf + 2) THEN
            ! GLOB. EQ.    O x x x 0    x x x O
            !              O x x x O [O x x x O]
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + blkp
            pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt + blkp
          ELSE
            ! GLOB. EQ.    O x x x    0 x x x O
            !             [O x x x O] O x x x O
            pos_start(1) = blkp + qq*blkt
            pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
            pos_start(Nf + 2) = 2*blkp + 2*qq*blkt
          END IF
        END IF
      END IF

    END IF
#else
    !***************************************************************************************
    !
    !                  CASE NO PARALLELLIZATION
    !
    !***************************************************************************************
    IF (Ntorloc == 1) THEN
      ! ****************************
      ! Case only 1 toroidal division
      ! ****************************
      pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
      pos_start(Nf + 2) = 0
    ELSE
      ! ****************************************
      ! Case only more than 1 toroidal divisions
      ! ****************************************
      IF (itor == 1) THEN
        ! first set of poloidal faces
        ! GLOB. EQ.    0 x x x O   O x x x
        !              O x x x O  [O x x x O]
        pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
        pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
      ELSE IF (itor == Ntorloc .and. Ntorloc == 2) THEN
        ! case last set with only 2 toroidal divisions
        ! GLOB. EQ.    0 x x x 0     x x x
        !             [O x x x O]  O x x x O
        pos_start(1) = blkp + qq*blkt
        pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
        pos_start(Nf + 2) = 0
      ELSE IF (itor == Ntorloc) THEN
        ! case last set with more then 2 toroidal divisions
        ! GLOB. EQ.    0 x x x O   O x x x
        !             [O x x x O]  O x x x O
        pos_start = pos_start + 2*blkp + qq*blkt
        pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
        pos_start(Nf + 2) = 0
      ELSE IF ((ifa == Nf + 2) .and. itor == Ntorloc-1) THEN
        ! case second-last set
        ! GLOB. EQ.    O O x x x 0    x x x
        !                O x x x O [O x x x O]
        pos_start = pos_start + blkp
        pos_start(2:Nf + 2) = pos_start(2:Nf + 2) + blkp
        pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt
      ELSE
        ! general case
        IF (ifa == Nf + 2) THEN
          ! GLOB. EQ.    O x x x 0    x x x O
          !              O x x x O [O x x x O]
          pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + blkp
          pos_start(Nf + 2) = pos_start(Nf + 2) + qq*blkt + blkp
        ELSE
          ! GLOB. EQ.    O x x x    0 x x x O
          !             [O x x x O] O x x x O
          pos_start(1) = blkp + qq*blkt
          pos_start(2:Nf + 1) = pos_start(2:Nf + 1) + 2*blkp + qq*blkt
          pos_start(Nf + 2) = 2*blkp + 2*qq*blkt
        END IF
      END IF
    END IF

#endif

    ! Filling ind_pos: determines the assembly position in cols and vals
    qq = 1
    DO i = 1, Nf + 2
      IF (i == 1 .or. i == Nf + 2) THEN
        blkl = blkp
      ELSE
        blkl = blkt
      END IF
      indpos(qq:qq + blkl-1) = (/(i, i=1, blkl)/) + pos_start(i)
      qq = qq + blkl
    END DO

  END SUBROUTINE getindpospol

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for toroidal interior faces
  !**************************************************
  SUBROUTINE getindposint(indpos)
    integer :: i, j, indpos(:), ieln, blkl
    integer :: pos(1:Nf)
    integer :: Fen(1:Nf), intface(5), qq
    logical :: Fdn(1:Nf)

    ! Extract infos on the element faces and the next element
    ! faces (the element connected to the actual one by the
    ! face Fig)
    intface = Mesh%intfaces(Fi,:)
    IF (intface(1) == iel) THEN
      ieln = intface(3)
      Fen = Mesh%F(ieln,:)
      Fdn = Mesh%Fdir(ieln,:)
    ELSE
      ieln = intface(1)
      Fen = Mesh%F(ieln,:)
      Fdn = Mesh%Fdir(ieln,:)
    END IF

    ! Filling pos: needed to determine the relative position of the
    ! toroidal faces in the assembling line of the matrix for the
    ! face Fig
    pos = 1
    DO i = 1, Nf
      IF (Fd(i)) CYCLE
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
      DO j = 1, Nf
        IF (Fen(j) == Fi) CYCLE
        IF (Fdn(j)) CYCLE
        IF (Fe(i) > Fen(j)) THEN
          pos(i) = pos(i) + 1
        ENDIF
      END DO
    END DO

    ! qq gives here the number of lateral faces assebled for this element
    ! (take out Dirichlet faces)
    qq = 0
    DO i = 1, Nf
      IF (Fd(i)) CYCLE
      qq = qq + 1
    END DO
    DO i = 1, Nf
      IF (Fdn(i)) CYCLE
      IF (Fen(i) == Fi) CYCLE
      qq = qq + 1
    END DO

    ! Filling pos_start: needed to determine the starting position of assembling
    ! the coefficients in the cols and vals vectors
    pos_start = 0
    pos_start(2:Nf + 1) = (pos-1)*blkt + 2*blkp
    pos_start(Nf + 2) = blkt*qq + 2*blkp
    IF (iel > ieln) THEN
      pos_start(1) = blkp
      pos_start(Nf + 2) = pos_start(Nf + 2) + blkp
    END IF

    IF (itor == numer%ntor .and. numer%ntor > 1) THEN
      pos_start = pos_start + 2*blkp
      pos_start(Nf + 2) = 0
      IF (iel > ieln) THEN
        pos_start(Nf + 2) = pos_start(Nf + 2) + blkp
      END IF
    ELSE IF (numer%ntor == 1) THEN
      IF (iel > ieln) THEN
        pos_start(Nf + 2) = blkp
      ELSE
        pos_start(Nf + 2) = 0
      END IF
    END IF

    ! Filling ind_pos: determines the assembly position in cols and vals
    qq = 1
    DO i = 1, Nf + 2
      IF (i == 1 .or. i == Nf + 2) THEN
        blkl = blkp
      ELSE
        blkl = blkt
      END IF
      indpos(qq:qq + blkl-1) = (/(i, i=1, blkl)/) + pos_start(i)
      qq = qq + blkl
    END DO
  END SUBROUTINE getindposint

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for toroidal exterior faces
  !**************************************************
  SUBROUTINE getindposext(indpos)
    integer :: i, j, indpos(:)
    integer :: pos(1:Nf), qq, blkl

    ! Filling pos: needed to determine the relative position of the
    ! toroidal faces in the assembling line of the matrix for the
    ! face Fig
    pos = 1
    DO i = 1, Nf
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
    END DO
    ! qq gives here the number of lateral faces assembled for this element
    ! (take out Dirichlet faces)
    qq = maxval(pos)

    ! Filling pos_start: needed to determine the starting position of assembling
    ! the coefficients in the cols and vals vectors
    pos_start = 0
    pos_start(2:Nf + 1) = (pos-1)*blkt + blkp
    pos_start(Nf + 2) = blkp + blkt*qq
    IF (itor == numer%ntor .and. numer%ntor > 1) THEN
      pos_start = pos_start + blkp
      pos_start(Nf + 2) = 0
    ELSE IF (numer%ntor == 1) THEN
      pos_start(Nf + 2) = 0
    END IF

    ! Filling ind_pos: determines the assembly position in cols and vals
    qq = 1
    DO i = 1, Nf + 2
      IF (i == 1 .or. i == Nf + 2) THEN
        blkl = blkp
      ELSE
        blkl = blkt
      END IF
      indpos(qq:qq + blkl-1) = (/(i, i=1, blkl)/) + pos_start(i)
      qq = qq + blkl
    END DO

  END SUBROUTINE getindposext

  SUBROUTINE fill_cols_vals_rowptr_loc2glob()
    integer :: sl, blkl
    CHARACTER(LEN=10)  :: num

    sl = delta(ifa)-1
    DO i = 1, blk
      rows(sl + i) = (i-1)*linew(Fig) + shift(Fig) + 1
      rhsvec(sl + i) = rhsvec(sl + i) + fel(i + ind_sta(ifa)-1,iel3)
#ifdef PARALL
      loc2glob(sl + i) = shiftGlob3d(ifa)-1 + i
#else
      loc2glob(sl + i) = sl + i
#endif

      DO ifl = 1, Nf + 2
        IF (ifl > 1 .and. ifl < Nf + 2) THEN
          IF (Fd(ifl-1)) CYCLE
          blkl = blkt
        ELSE
          blkl = blkp
        END IF
        DO jj = 1, blkl
          j = jj + ind_sta(ifl)-1
          k = (i-1)*linew(Fig) ! local shift
          cols(shift(Fig) + k + indpos(j)) = indglo(j)
          vals(shift(Fig) + k + indpos(j)) = vals(shift(Fig) + k + indpos(j)) + Kel(i + ind_sta(ifa)-1, j,iel3)
        END DO
      END DO
    END DO

  END SUBROUTINE fill_cols_vals_rowptr_loc2glob

  !*********************************************************
  ! Deallocate and reallocate MatK and rhs. Initialize to 0
  !*********************************************************
  SUBROUTINE init_mat()
    IF (ASSOCIATED(MatK%cols)) THEN
      DEALLOCATE (MatK%cols)
    END IF
    IF (ASSOCIATED(MatK%rowptr)) THEN
      DEALLOCATE (MatK%rowptr)
    END IF
    IF (ASSOCIATED(MatK%vals)) THEN
      DEALLOCATE (MatK%vals)
    END IF
    IF (ASSOCIATED(MatK%loc2glob)) THEN
      DEALLOCATE (MatK%loc2glob)
    END IF
    IF (ASSOCIATED(rhs%loc2glob)) THEN
      DEALLOCATE (rhs%loc2glob)
    END IF
    IF (ASSOCIATED(rhs%vals)) THEN
      DEALLOCATE (rhs%vals)
    END IF

    MatK%n = nn
    MatK%nnz = nnz
    rhs%n = nn
    ALLOCATE (MatK%cols(nnz))
    ALLOCATE (MatK%vals(nnz))
    ALLOCATE (MatK%rowptr(nn + 1))
    ALLOCATE (MatK%loc2glob(nn))
    ALLOCATE (rhs%vals(nn))
    ALLOCATE (rhs%loc2glob(nn))
    MatK%cols = 0
    MatK%rowptr = 0
    MatK%vals = 0.
    MatK%loc2glob = 0
    rhs%vals = 0.
    rhs%loc2glob = 0

  END SUBROUTINE init_mat

END SUBROUTINE HDG_assembly
#else

!***********************************************************************
!
!                            VERSION 2D
!
!***********************************************************************
SUBROUTINE HDG_assembly()
  USE globals
  USE matrices_types
  USE printUtils
  USE debug
  USE LinearAlgebra, ONLY: tensorSumInt
  USE MPI_OMP

  integer                        :: i, j, jj, k, iel, ifa, ifl, nnz, nn, cont, ieloc, ct_sc
  integer                        :: ieln,ifan,Fi_per,Fe_aux(refElPol%Nfaces)
  integer                        :: Neq, Nf, Nfaces, Nfp, Nel, Nintf, Nextf, Ndirf, Nfunk, blk, Ndim, Np
  integer                        :: nColsPosPerFaceInt, nColsPosPerFaceExt, nColsPosPerElem
  integer                        :: Fe(refElPol%Nfaces), Fi
  integer                        :: indglo(1:refElPol%Nfacenodes*refElPol%Nfaces*phys%Neq)
  integer                        :: indpos(1:refElPol%Nfacenodes*refElPol%Nfaces*phys%Neq)
  integer*4                      :: ind_loc(1:refElPol%Nfaces, 1:refElPol%Nfacenodes*phys%Neq)
  integer                        :: n_periodic_faces
  integer, dimension(Mesh%Nfaces) :: linew, shift
  integer, pointer, dimension(:)   :: cols, rows, loc2glob
  real*8, pointer, dimension(:)    :: vals, rhsvec
  logical                        :: Fd(refElPol%Nfaces)
  real*8, allocatable             :: Kel(:,:), fel(:)
  real*8, pointer                 :: Df(:,:), Hf(:,:), Ef(:,:)
  real*8, pointer                 :: UU(:,:), U0(:)
  real*8, allocatable             :: fh(:)
  real*8, allocatable             :: Lf(:,:)
  real*8, pointer                 :: Qf(:,:)
  real*8, pointer                 :: LL(:,:), L0(:)
#ifdef PARALL
  integer                        ::        Fasind(Mesh%Nfaces)
  integer                        :: rept
#endif
  integer :: ierr


  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1, timing%clock_rate1)
  end if

  ! Number of periodic faces
  n_periodic_faces=0
  do i=1,Mesh%nextfaces
    if (Mesh%periodic_faces(i).ne.0) then
      n_periodic_faces = n_periodic_faces+1
    endif
  enddo

  Neq = phys%Neq
  Np = refElPol%Nnodes2D
  Ndim = Mesh%Ndim
  Nf = refElPol%Nfaces
  Nfp = refElPol%Nfacenodes
  Nel = Mesh%Nelems
  Nfaces = Mesh%Nfaces
  Nintf = Mesh%Nintfaces ! Number of interior faces in the mesh
  Nextf = Mesh%Nextfaces ! Number of exterior faces in the mesh
  Ndirf = Mesh%Ndir      ! Number of Dirichlet faces in the mesh
  Nfunk = Mesh%ukf       ! Number of unknown faces in the mesh (Nfaces-Ndirf)

  ! Compute the dimension and the number of nonzero of the global matrix
  blk = Nfp*Neq
  nn = blk*Nfunk   ! K = nn x nn

#ifdef PARALL
  nn = nn-blk*Mesh%nghostfaces
  Fasind = 0
  rept = 0
  DO i = 1, Nfaces
    IF (Mesh%ghostfaces(i) .eq. 1) THEN
      rept = rept + 1
      CYCLE
    END IF
    Fasind(i) = i-rept
  END DO
#endif

  ! Compute nnz and line width
  CALL computennz()
  !call displayVectorint(shift)
  !call displayVectorint(linew)

  ! Deallocate and reallcate pointers of  MatK struct and rhs struct initialize to 0
  CALL init_mat()
  cols => MatK%cols
  rows => MatK%rowptr
  vals => MatK%vals
  loc2glob => MatK%loc2glob
  rhsvec => rhs%vals

  ! Number of cols position filled by each face (int & ext)
  !   nColsPosPerFaceInt = Nfp**2*Neq**2*(2*Nf-1)
  !   nColsPosPerFaceExt = Nfp**2*Neq**2*Nf
  !   nColsPosPerElem = Nfp**2*Neq**2*Nf

  ind_loc = 0
  DO i = 1, Nf
    DO j = 1, Neq*Nfp
      ind_loc(i, j) = Neq*Nfp*(i-1) + j
    END DO
  END DO

  ALLOCATE (Kel(1:Nf*Nfp*Neq, 1:Nf*Nfp*Neq))
  ALLOCATE (fel(1:Nf*Nfp*Neq))

  ! Counter
  cont = 0
  !**************************
  ! Loop in the mesh elements
  !**************************
  DO iel = 1, Nel

    Fe = Mesh%F(iel,:)
    Fd = Mesh%Fdir(iel,:)

    CALL computeElementalMatrix()

    ! call print_elmat_tex(Kel,fel,iel)

    !***********************
    ! Loop in local faces
    !***********************
    DO ifa = 1, Nf

      Fi = Fe(ifa)
#ifdef PARALL
      IF (Mesh%ghostFaces(Fi) .eq. 1) CYCLE
#endif
      IF (Fi .le. Nintf) THEN
        !**************************
        ! Interior face
        !**************************
        CALL getindposint(indpos)

      ELSE
        !**************************
        ! Exterior face
        !**************************
        IF (Mesh%Fdir(iel, ifa)) CYCLE
        if (Mesh%periodic_faces(Fi-Nintf).eq.0) then
          CALL getindposext(indpos)
        else
          ieln = Mesh%extfaces(Mesh%periodic_faces(Fi-Nintf),1)
          ifan = Mesh%extfaces(Mesh%periodic_faces(Fi-Nintf),2)
          Fi_per = Mesh%F(ieln,ifan)
          !write(6,*) "iel:",iel
          !write(6,*) "ieln:",ieln
          !write(6,*) "Fi:",Fi
          !write(6,*) "Fi_per:",Fi_per
          CALL getindposextperiodic(indpos)
          CALL fill_cols_vals_rowptr_loc2glob()
          Fe_aux = Fe
          call exchange_Fi ! Fi--->Fi_per ; Fi_per--->Fi
          CALL getindposextperiodic(indpos)
          CALL fill_cols_vals_rowptr_loc2glob()
          Fe = Fe_aux
          !stop
          CYCLE
        endif

      END IF

      CALL fill_cols_vals_rowptr_loc2glob()
    END DO

  END DO ! Loop in elements

  rows(nn + 1) = nnz + 1
  rhs%loc2glob = loc2glob

  DEALLOCATE (Kel, fel)
  NULLIFY (cols, rows, vals, loc2glob, rhsvec)
  !call displayVectorInt(loc2glob)
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1, timing%clock_rate1)
    timing%runtass = timing%runtass + (timing%cke1-timing%cks1)/real(timing%clock_rate1)
    timing%cputass = timing%cputass + timing%tpe1-timing%tps1
  end if


CONTAINS

  subroutine exchange_Fi
    integer :: aux,i
    do i=1,Nf
      if (Fe(i)==Fi) then
        Fe(i)=Fi_per
      endif
    end do
    aux = Fi
    Fi = Fi_per
    Fi_per = aux
  end subroutine exchange_Fi

  !**************************************************
  ! Compute the nnz of the matrix
  !**************************************************
  SUBROUTINE computennz()
    nnz = 0
    linew = 0
    shift = 0
    ! Loop in faces
    DO Fi = 1, Nfaces
#ifdef PARALL
      IF (Mesh%ghostfaces(Fi) .eq. 1) THEN
        IF (Fi .gt. 1) THEN
          shift(Fi) = shift(Fi-1)
          linew(Fi) = linew(Fi-1)
        END IF
        CYCLE
      END IF
#endif
      IF (Fi .le. Nintf) THEN
        !**************************
        ! Interior face
        !**************************
        DO ieloc = 1, 2
          iel = Mesh%intfaces(Fi, (ieloc-1)*2 + 1)
          Fe = Mesh%F(iel,:)
          Fd = Mesh%Fdir(iel,:)
          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            IF (Fe(ifl) == Fi .and. ieloc == 2) CYCLE
            nnz = nnz + blk**2
            linew(Fi) = linew(Fi) + blk
          END DO
        END DO
      ELSE
        !**************************
        ! Exterior face
        !**************************
        if (Mesh%periodic_faces(Fi-Nintf).eq.0) then
          iel = Mesh%extfaces(Fi-Nintf, 1)
          ifa = Mesh%extfaces(Fi-Nintf, 2)
          Fd = Mesh%Fdir(iel,:)

          ! Skip Dirichlet faces
          IF (Fd(ifa)) THEN
            IF (Fi .gt. 1) THEN
              shift(Fi) = shift(Fi-1)
              linew(Fi) = linew(Fi-1)
            END IF
            CYCLE
          END IF

          DO ifl = 1, Nf
            IF (Fd(ifl)) CYCLE
            nnz = nnz + blk**2
            linew(Fi) = linew(Fi) + blk
          END DO
        else
          ! periodic face
          ieln = Mesh%extfaces(Mesh%periodic_faces(Fi-Nintf),1)
          ifan = Mesh%extfaces(Mesh%periodic_faces(Fi-Nintf),2)
          Fi_per = Mesh%F(ieln,ifan)
          DO ieloc = 1, 2
            if (ieloc==1) then
              iel=Mesh%extfaces(Fi-Nintf, 1)
            else
              iel=Mesh%extfaces(Mesh%periodic_faces(Fi-Nintf),1)
            endif
            Fe = Mesh%F(iel,:)
            Fd = Mesh%Fdir(iel,:)
            DO ifl = 1, Nf
              IF (Fd(ifl)) CYCLE
              !						               IF ( Fe(ifl) == Fi  .and. ieloc == 2) CYCLE
              IF ( Fe(ifl) == Fi_per) CYCLE
              nnz = nnz + blk**2
              linew(Fi) = linew(Fi) + blk
            END DO
          END DO
        endif

      END IF

      IF (Fi .gt. 1) THEN
        shift(Fi) = shift(Fi-1) + linew(Fi-1)*blk
      END IF
    END DO
  END SUBROUTINE computennz

  !**************************************************
  ! Generate the elemental matrix
  !**************************************************
  SUBROUTINE computeElementalMatrix_old()
    Df => elMat%Df(:,:, iel)
    UU => elMat%UU(:,:, iel)
    Hf => elMat%Hf(:,:, iel)
    Ef => elMat%Ef(:,:, iel)
    U0 => elMat%U0(:, iel)
    ALLOCATE (Lf(Neq*Nf*Nfp, Neq*Ndim*Np))
    ALLOCATE (fh(Neq*Nf*Nfp))
    LL => elMat%LL(:,:, iel)
    Lf = elMat%Lf(:,:, iel)
    fh = elMat%fh(:, iel)
#ifdef TEMPERATURE
    Lf = Lf + elMat%TQhf(:,:, iel)
    fh = fh-elMat%Tfhf(:, iel)
#endif
    Qf => elMat%Qf(:,:, iel)
    L0 => elMat%L0(:, iel)
    IF (switch%shockcp .gt. 0) THEN
      IF (Mesh%flag_elems_sc(iel) .ne. 0) THEN
        Lf = Lf + elMat%Lf_sc(:,:, Mesh%flag_elems_sc(iel))
      END IF
    END IF
    ! Faces numbering and Dirichlet booleans
    Fe = Mesh%F(iel,:)
    Fd = Mesh%Fdir(iel,:)

    ! Elemental matrix and elemental RHS
    Kel = matmul(Df, UU) + Hf-Ef
    fel = -matmul(Df, U0) + fh
    Kel = Kel-matmul((Lf-Qf), LL)
    fel = fel + matmul((Lf-Qf), L0)
    !if (iel==2) THEN
    !call saveMatrix(Kel,"Kel")
    !call saveVector(fel,"fel")
    !call saveMatrix(Hf,"Hf")
    !call saveMatrix(Lf,"Lf")
    !call saveMatrix(Qf,"Qf")
    !call saveMatrix(Ef,"Ef")
    !call saveMatrix(Df,"Df")
    !stop
    !endif
    DEALLOCATE (Lf, fh)
    NULLIFY (LL, Qf, L0)
    NULLIFY (Df, UU, Hf, Ef, U0)
  END SUBROUTINE computeElementalMatrix_old

  SUBROUTINE computeElementalMatrix()
    real*8, pointer                 :: A_lq(:,:), A_lu(:,:), A_ll(:,:)
    real*8, pointer                 :: LL(:,:), L0(:), UU(:,:), U0(:), f(:)
    LL => elMat%LL(:,:, iel)
    L0 => elMat%L0(:, iel)
    UU => elMat%UU(:,:, iel)
    U0 => elMat%U0(:, iel)

    A_lq => elMat%Alq(:,:, iel)
    A_lu => elMat%Alu(:,:, iel)
    A_ll => elMat%All(:,:, iel)
    f => elMat%fh(:, iel)

    ! Elemental matrix and elemental RHS
    Kel = matmul(A_lq, LL) + matmul(A_lu, UU) + A_ll
    fel = -(matmul(A_lu, U0) + matmul(A_lq, L0)) + f

    !call HDF5_save_vector(fel,"f1")
    !stop
    NULLIFY (LL, L0, UU, U0)
    NULLIFY (A_lq, A_lu, A_ll)
  END SUBROUTINE computeElementalMatrix

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for interior faces
  !**************************************************
  SUBROUTINE getindposint(indpos)
    integer :: i, j, indpos(1:neq*Nfp*Nf)
    integer :: pos(1:Nf)
    integer :: Fen(1:Nf), intface(5)
    logical :: Fdn(1:Nf)

    intface = Mesh%intfaces(Fi,:)

    IF (intface(1) == iel) THEN
      Fen = Mesh%F(intface(3),:)
      Fdn = Mesh%Fdir(intface(3),:)
    ELSE
      Fen = Mesh%F(intface(1),:)
      Fdn = Mesh%Fdir(intface(1),:)
    END IF

    pos = 1
    DO i = 1, Nf
      IF (Fd(i)) CYCLE
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
      DO j = 1, Nf
        IF (Fen(j) == Fi) CYCLE
        IF (Fdn(j)) CYCLE
        IF (Fe(i) > Fen(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
    END DO
    indpos = reshape(TensorSumInt((/(i, i=1, neq*Nfp)/), (pos-1)*neq*Nfp), (/neq*Nfp*Nf/))

  END SUBROUTINE getindposint

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for exterior faces
  !**************************************************
  SUBROUTINE getindposext(indpos)
    integer :: i, j, indpos(1:neq*Nfp*Nf)
    integer :: pos(1:Nf)

    pos = 1
    DO i = 1, Nf
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
    END DO
    indpos = reshape(TensorSumInt((/(i, i=1, neq*Nfp)/), (pos-1)*neq*Nfp), (/neq*Nfp*Nf/))
  END SUBROUTINE getindposext

  !**************************************************
  ! This routine computes the indices to
  ! determine the assembly positions of cols and vals
  ! for exterior periodic faces
  !**************************************************
  SUBROUTINE getindposextperiodic(indpos)
    integer :: i, j, indpos(1:neq*Nfp*Nf)
    integer :: pos(1:Nf)
    integer :: Fen(1:Nf)
    logical :: Fdn(1:Nf)

    !      intface = Mesh%intfaces(Fi,:)

    !      IF (intface(1) == iel) THEN
    !         Fen = Mesh%F(intface(3),:)
    !         Fdn = Mesh%Fdir(intface(3),:)
    !      ELSE
    !         Fen = Mesh%F(intface(1),:)
    !         Fdn = Mesh%Fdir(intface(1),:)
    !      END IF

    Fen = Mesh%F(ieln,:)
    Fdn = Mesh%Fdir(ieln,:)

    do i=1,Nf
      if (Fen(i)==Fi_per) then
        Fen(i) = Fi
      endif
    end do

    !write(6,*) "Fi:",Fi
    !write(6,*) "Fe:",Fe
    !write(6,*) "Fen:",Fen
    pos = 1
    DO i = 1, Nf
      IF (Fd(i)) CYCLE
      DO j = 1, Nf
        IF (i == j) CYCLE
        IF (Fd(j)) CYCLE
        IF (Fe(i) > Fe(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
      DO j = 1, Nf
        IF (Fen(j) == Fi) CYCLE
        IF (Fdn(j)) CYCLE
        IF (Fe(i) > Fen(j)) THEN
          pos(i) = pos(i) + 1
        END IF
      END DO
    END DO
    indpos = reshape(TensorSumInt((/(i, i=1, neq*Nfp)/), (pos-1)*neq*Nfp), (/neq*Nfp*Nf/))

    !write(6,*) "pos:",pos
    !call displayVectorInt(indpos)
  END SUBROUTINE getindposextperiodic

  !**************************************************
  ! Fill cols, vals, rowptr, loc2glob
  !**************************************************
  SUBROUTINE fill_cols_vals_rowptr_loc2glob()
    integer :: sl
    ! Find global indices to fill cols and vals
#ifdef PARALL
    indglo = reshape(TensorSumInt((/(i, i=1, neq*Nfp)/), (Mesh%loc2glob_fa(Fe)-1)*neq*Nfp), (/neq*Nfp*Nf/))
#else
    indglo = reshape(TensorSumInt((/(i, i=1, neq*Nfp)/), (Fe-1)*neq*Nfp), (/neq*Nfp*Nf/))
#endif

#ifdef PARALL
    sl = (Fasind(Fi)-1)*blk ! shift linear
#else
    sl = (Fi-1)*blk ! shift linear
#endif
    DO i = 1, blk
      rows(sl + i) = (i-1)*linew(Fi) + shift(Fi) + 1
      rhsvec(sl + i) = rhsvec(sl + i) + fel(i + (ifa-1)*blk)
#ifdef PARALL
      loc2glob(sl + i) = (Mesh%loc2glob_fa(Fi)-1)*blk + i
#else
      loc2glob(sl + i) = sl + i
#endif

      DO ifl = 1, Nf
        IF (Fd(ifl)) CYCLE
        DO jj = 1, blk
          j = jj + blk*(ifl-1)
          k = (i-1)*linew(Fi) ! local shift
          cols(shift(Fi) + k + indpos(j)) = indglo(j)
          vals(shift(Fi) + k + indpos(j)) = vals(shift(Fi) + k + indpos(j)) + Kel(i + (ifa-1)*blk, j)
        END DO
      END DO
    END DO
  END SUBROUTINE fill_cols_vals_rowptr_loc2glob

  !*********************************************************
  ! Deallocate and reallocate MatK and rhs. Initialize to 0
  !*********************************************************
  SUBROUTINE init_mat()
    IF (ASSOCIATED(MatK%cols)) THEN
      DEALLOCATE (MatK%cols)
    END IF
    IF (ASSOCIATED(MatK%rowptr)) THEN
      DEALLOCATE (MatK%rowptr)
    END IF
    IF (ASSOCIATED(MatK%vals)) THEN
      DEALLOCATE (MatK%vals)
    END IF
    IF (ASSOCIATED(MatK%loc2glob)) THEN
      DEALLOCATE (MatK%loc2glob)
    END IF
    IF (ASSOCIATED(rhs%loc2glob)) THEN
      DEALLOCATE (rhs%loc2glob)
    END IF
    IF (ASSOCIATED(rhs%vals)) THEN
      DEALLOCATE (rhs%vals)
    END IF

    MatK%n = nn
    MatK%nnz = nnz
    rhs%n = nn
    ALLOCATE (MatK%cols(nnz))
    ALLOCATE (MatK%vals(nnz))
    ALLOCATE (MatK%rowptr(nn + 1))
    ALLOCATE (MatK%loc2glob(nn))
    ALLOCATE (rhs%vals(nn))
    ALLOCATE (rhs%loc2glob(nn))
    MatK%cols = 0
    MatK%rowptr = 0
    MatK%vals = 0.
    MatK%loc2glob = 0
    rhs%vals = 0.
    rhs%loc2glob = 0

  END SUBROUTINE init_mat

END SUBROUTINE HDG_assembly
#endif
