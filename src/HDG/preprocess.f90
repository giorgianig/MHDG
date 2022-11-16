!************************************************************
! project: MHDG
! file: preprocess.f90
! date: 16/12/2016
! Preprocess of the mesh: creates Tlin, F, N, flipFaces,
! innerfaces, outerfaces
!************************************************************
MODULE preprocess
  USE types
  USE globals
  USE printutils
  USE MPI_OMP

  IMPLICIT NONE
CONTAINS

  SUBROUTINE mesh_preprocess

    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        IF (MPIvar%glob_id .eq. 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*                MESH PREPROCESS                *'
          WRITE (6, *) '*************************************************'
        ENDIF
      END IF
    ENDIF
    !**************************************
    ! Create the linear connectivity matrix
    !**************************************
    ALLOCATE (Mesh%Tlin(Mesh%Nelems, 1:refElPol%Nvertices))
    Mesh%Tlin = Mesh%T(:, 1:refElPol%Nvertices)

    !**********************************
    ! Create nodal connectivity matrix
    !**********************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Creating nodal connectivity "
      END IF
    ENDIF
    CALL createNodalConnectivity()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
      WRITE (6, *) "Nodal connectivity matrix N: "
      CALL displayMatrixInt(Mesh%N)
    END IF
    !**********************************
    ! Create outerfaces and innerfaces
    !**********************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Creating faces structures "
      END IF
    END IF
    CALL GetFaces()
    ! Print out stuff...
#ifdef PARALL
    IF (utils%printint > 0) THEN
      WRITE (6, '(A,I13)') "Process: ", MPIvar%glob_id, " Number of interior faces:           ", Mesh%Nintfaces
      WRITE (6, '(A,I13)') "Process: ", MPIvar%glob_id, " Number of faces:                    ", Mesh%Nfaces
    END IF
#else
    IF (utils%printint > 0) THEN
      WRITE (6, '(A,I13)') " Number of interior faces:           ", Mesh%Nintfaces
      WRITE (6, '(A,I13)') " Number of faces:                    ", Mesh%Nfaces
    END IF
#endif
    IF (utils%printint > 1) THEN
      WRITE (6, *) "Exterior faces: "
      CALL displayMatrixInt(Mesh%extfaces)
      WRITE (6, *) "Interior faces: "
      CALL displayMatrixInt(Mesh%intfaces)
    END IF

    !*******************************
    ! Boundary condition preprocess
    !*******************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Boundary conditions preprocess "
      END IF
    ENDIF
    CALL Bc_preprocess()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
      WRITE (6, *) "Exterior faces after preprocess: "
      CALL displayMatrixInt(Mesh%extfaces)
      WRITE (6, *) "Boundary connectivity matrix Tb after preprocess: "
      CALL displayMatrixInt(Mesh%Tb)
    END IF

    !*****************************************
    ! Create face connectivity F and flipFaces
    !*****************************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Creating face connectivity "
      END IF
    ENDIF
    CALL CreateFaceConnectivity()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
      WRITE (6, *) "Face connectivity matrix F: "
      CALL displayMatrixInt(Mesh%F)
      WRITE (6, *) "F-dir: "
      CALL displayMatrixLog(Mesh%Fdir)
    END IF

    !*****************************************
    ! Compute max and min element size
    !*****************************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Computing element size "
      END IF
    ENDIF
    CALL computeElementSize()
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Done! "
      END IF
    ENDIF

    !*****************************************
    ! Compute puff area
    !*****************************************
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Computing puff area"
      END IF
    ENDIF
! It does not work in parallel (TODO)
#ifndef PARALL
    CALL computePuffArea()
#endif
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Puff area:  ", Mesh%puff_area*phys%lscale*phys%lscale, " m^2"
      END IF
    ENDIF
! It does not work in parallel (TODO)
#ifndef PARALL
    CALL computeCoreArea()
#endif
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Core area:  ", Mesh%core_area*phys%lscale*phys%lscale, " m^2"
        WRITE (6, *) "Done! "
      END IF
    ENDIF    
  END SUBROUTINE mesh_preprocess

  !********************
  ! Nodal Connectivity
  !********************
  SUBROUTINE createNodalConnectivity()

    integer :: nc, k, iel, nn_Te(1:Mesh%Nnodesperelem)
    integer :: nn(1:Mesh%Nnodes), Te(1:Mesh%Nnodesperelem)
    ! Find the number of colums of the N matrix
    nn = 1
    DO iel = 1, Mesh%Nelems
      Te = Mesh%T(iel, :)
      nn(Te) = nn(Te) + 1
    END DO
    nc = MAXVAL(nn)

    ALLOCATE (Mesh%N(Mesh%Nnodes, nc))
    Mesh%N = 0
    nn = 1
    nn_Te = 0
    DO iel = 1, Mesh%Nelems
      Te = Mesh%T(iel, :)
      nn_Te = nn(Te)
      DO k = 1, Mesh%Nnodesperelem
        Mesh%N(Te(k), nn_Te(k)) = iel
      END DO
      nn(Te) = nn(Te) + 1
    END DO
  END SUBROUTINE createNodalConnectivity

  !********************************************
  ! Create extFaces and intFaces: the exterior
  ! faces have the same numbering than in Tb
  !********************************************
  SUBROUTINE GetFaces()

    integer :: iel, jel, ifa, jfa, ie, ii, node1
    integer :: Efaces(refElPol%Nfaces, refElPol%Nfacenodeslin), nf(1:refElPol%Nfacenodeslin)
    integer, allocatable :: temp_intFaces(:, :)
    logical, allocatable :: markE(:, :), markF(:)

    !Definition of the faces in the reference element
    SELECT CASE (refElPol%elemType)
    CASE (0) !triangle
      Efaces = reshape((/1, 2, 3, 2, 3, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (1) !quadrilaterals
      Efaces = reshape((/1, 2, 3, 4, 2, 3, 4, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (2) !thetrahedra
      !            Efaces = computeNodesFacesTetra(1);
    CASE (3) !hexahedra
      !        Efaces = computeNodesFacesHexa(1);
    END SELECT

    ALLOCATE (markE(1:Mesh%Nelems, 1:refElPol%Nfaces))
    ALLOCATE (markF(1:Mesh%Nextfaces))
    ALLOCATE (Mesh%extFaces(1:Mesh%Nextfaces, 2))
    ALLOCATE (temp_intFaces(2*Mesh%Nelems, 5))

    ie = 0
    ii = 0
    markE = .false.
    markF = .false.
    DO iel = 1, Mesh%Nelems
      DO ifa = 1, refElPol%Nfaces
        IF (.not. markE(iel, ifa)) THEN
          markE(iel, ifa) = .true.
          nf = Mesh%Tlin(iel, Efaces(ifa, :))
          CALL FindElem(nf, iel, jel)
          IF (jel .ne. 0) THEN
            ! Inner element
            ii = ii + 1
            CALL FindFace(nf, Mesh%Tlin(jel, :), Efaces, jfa, node1);
            temp_intFaces(ii, 1) = iel
            temp_intFaces(ii, 2) = ifa
            temp_intFaces(ii, 3) = jel
            temp_intFaces(ii, 4) = jfa
            temp_intFaces(ii, 5) = node1
            markE(jel, jfa) = .true.
          ELSE
            ! Boundary element
            CALL FindCorrespondingFace(iel, ifa, markF, ie)
            !                                                                                           ie = ie + 1
            Mesh%extFaces(ie, 1) = iel
            Mesh%extFaces(ie, 2) = ifa
          END IF
        END IF
      END DO
    END DO

    Mesh%Nintfaces = ii
    ALLOCATE (Mesh%intFaces(ii, 5))
    Mesh%intFaces = temp_intFaces(1:ii, :)
    Mesh%Nfaces = Mesh%Nextfaces + Mesh%Nintfaces
    DEALLOCATE (markE, markF, temp_intFaces)

  CONTAINS

    !*************************************************
    ! Find element: find the neighboring element to the
    ! element i connected by the nodes defined in nf
    !*************************************************
    SUBROUTINE FindElem(nf, iel, jel)
      integer, intent(IN) :: nf(1:refElPol%Nfacenodeslin)
      integer, intent(IN) :: iel
      integer, intent(OUT):: jel
      integer :: i, z, kel, check, k, zel
      integer, dimension(1:size(Mesh%N, 2))::elems1, elemsi

      jel = 0
      elems1 = Mesh%N(nf(1), :)

      ! Loop in the element list associated to the 1st node
      DO k = 1, size(Mesh%N, 2)
        kel = elems1(k)

        IF (kel .eq. 0) EXIT
        IF (kel .eq. iel) CYCLE ! I skip the element i

        ! Loop in the connecting node list
        DO i = 2, refElPol%Nfacenodeslin
          elemsi = Mesh%N(nf(i), :)
          check = 0

          ! Loop in the element list of the i-th node
          DO z = 1, size(Mesh%N, 2)
            zel = elemsi(z)

            IF (zel .eq. 0) EXIT
            IF (zel .eq. iel) CYCLE ! I skip the element i
            IF (zel == kel) THEN
              jel = kel
              check = 1
              EXIT
            END IF
          END DO

          IF (check == 0) jel = 0 ! if check=0 means I didn't find any equal element in this node
        END DO
        IF (jel .ne. 0) EXIT ! I found the element
      END DO
    END SUBROUTINE FindElem

    !*************************************************
    ! Find face: find the local face number in the
    ! element and the matching node for the first node
    !*************************************************
    SUBROUTINE FindFace(nfi, nodesE, Efaces, jfa, node1)
      integer, intent(IN) :: nfi(:)
      integer, intent(IN) :: nodesE(:)
      integer, intent(IN) :: Efaces(:, :)
      integer, intent(OUT):: jfa
      integer, intent(OUT):: node1
      integer :: check, fn, i, j
      integer :: nfj(1:refElPol%Nfacenodeslin)

      ! Find the corresponding face
      DO jfa = 1, refElPol%Nfaces
        nfj = nodesE(Efaces(jfa, :))
        ! check if the set of nodes nfj is equal to the set of nodes nfi
        fn = 0
        DO i = 1, refElPol%Nfacenodeslin
          check = 0
          DO j = 1, refElPol%Nfacenodeslin
            IF (nfi(i) .eq. nfj(j)) THEN
              check = 1
              fn = fn + 1
              EXIT
            END IF
          END DO
          IF (check .eq. 0) THEN ! it is not the right face
            EXIT
          END IF
        END DO
        IF (fn .eq. refElPol%Nfacenodeslin) THEN
          EXIT
        END IF
      END DO

      IF (fn .ne. refElPol%Nfacenodeslin) THEN
        write (6, *) "Error! Corresponding face not found"
        STOP
      END IF

      ! Find the corresponding node
      DO i = 1, refElPol%Nfacenodeslin
        IF (nfj(i) .eq. nfi(1)) EXIT
      END DO
      IF (i .eq. refElPol%Nfacenodeslin + 1) THEN
        write (6, *) "Error! Corresponding node not found"
        STOP
      END IF
      node1 = i

    END SUBROUTINE FindFace

    !*************************************************
    ! Find where to place the current exterior face in
    ! order to match the order of Tb
    !*************************************************
    SUBROUTINE FindCorrespondingFace(iel, ifa, markF, iextf)
      integer, intent(IN)    :: iel, ifa
      integer, intent(OUT)   :: iextf
      logical, intent(INOUT) :: markF(:)
      integer               :: i
      integer               :: nodes(1:Mesh%Nnodesperface), nodes_check(1:Mesh%Nnodesperface)
      logical               :: isok

      nodes = Mesh%T(iel, refElPol%face_nodes(ifa, :))
      DO iextf = 1, Mesh%Nextfaces

        IF (markF(iextf)) CYCLE
        nodes_check = Mesh%Tb(iextf, :)! <--- Here I use Tb!!
        isok = .true.

        ! Check if "nodes" and "nodes_check" have the same nodes
        DO i = 1, Mesh%Nnodesperface
          IF (nodes(i) .ne. nodes_check(i)) THEN
            isok = .false.
            EXIT
          END IF
        END DO
        IF (isok) THEN
          markF(iextf) = .true.
          EXIT
        END IF
      END DO

      IF (iextf == Mesh%Nextfaces + 1) THEN
        WRITE (6, *) 'Error! Corresponding face in Tb not found'
        STOP
      END IF

    END SUBROUTINE FindCorrespondingFace

  END SUBROUTINE GetFaces

  !********************************************
  ! Identify the number of boundaries and place
  ! Dirichlet boundaries at the bottom of the
  ! list
  !********************************************
  SUBROUTINE Bc_preprocess()
    integer :: fl, nb, ndir, ifa, bt, idir
    integer :: df(max_num_diff_bc)
    integer, allocatable :: aux_Tb(:, :), aux_extFaces(:, :), aux_boundaryflag(:)

    nb = 0
    df = 0
    ndir = 0

    ! Find how many different boundary types are present (different
    ! boundary conditions) and how many Dirichlet faces are used
    DO ifa = 1, Mesh%Nextfaces
      fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face
  
      
#ifdef PARALL
      IF (fl .eq. 0) CYCLE
#endif
      IF (df(phys%bcflags(fl)) == 0) THEN
        df(phys%bcflags(fl)) = 1
        nb = nb + 1
      END IF
      IF (phys%bcflags(fl) == bc_dirichlet) THEN
        ndir = ndir + 1
      END IF
    END DO


    ! Here I place Dirichlet boundary faces at the bottom of Tb: I change
    ! the order of the faces only in case that at least 1 Dirichlet
    ! face exists and there are more then one type of boundary conditions.
    ! This guarantees compatibility with the Matlab version of the code
    IF (nb > 1 .and. ndir > 0) THEN
      ALLOCATE (aux_Tb(size(Mesh%Tb, 1), size(Mesh%Tb, 2)))
      ALLOCATE (aux_extFaces(size(Mesh%extFaces, 1), size(Mesh%extFaces, 2)))
      ALLOCATE (aux_boundaryflag(Mesh%NextFaces))
      aux_boundaryflag = 0
      aux_Tb = 0
      aux_extFaces = 0
      idir = 0
      DO ifa = 1, Mesh%Nextfaces
        fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face
#ifdef PARALL
        IF (fl .eq. 0) CYCLE ! Add by Benjamin
#endif
        IF (phys%bcflags(fl) == bc_dirichlet) THEN
          aux_Tb(ifa + Mesh%Nextfaces - ndir, :) = Mesh%Tb(ifa, :)
          aux_extFaces(ifa + Mesh%Nextfaces - ndir, :) = Mesh%extFaces(ifa, :)
          aux_boundaryflag(ifa + Mesh%Nextfaces - ndir) = Mesh%boundaryFlag(ifa)
          idir = idir + 1
        ELSE
          aux_Tb(ifa - idir, :) = Mesh%Tb(ifa, :)
          aux_extFaces(ifa - idir, :) = Mesh%extFaces(ifa, :)
          aux_boundaryflag(ifa - idir) = Mesh%boundaryFlag(ifa)
        END IF
      END DO
      Mesh%Tb = aux_Tb
      Mesh%extFaces = aux_extFaces
      Mesh%boundaryFlag = aux_boundaryflag
      DEALLOCATE (aux_Tb, aux_extFaces, aux_boundaryflag)
    END IF

    CALL set_boundary_flag_names()

    ! Periodic faces preprocess
    allocate(Mesh%periodic_faces(Mesh%Nextfaces))
    Mesh%periodic_faces = 0
    do ifa = 1, Mesh%Nextfaces
      fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
      IF (fl .eq. 0) CYCLE ! Add by Benjamin
#endif
      if (phys%bcflags(fl) == bc_periodic) then
        call periodic_faces_preprocess()
        exit
      end if
    end do

    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) '*************************************************'
      WRITE (6, *) '*          BOUNDARY CONDITIONS                  *'
      WRITE (6, *) '*************************************************'
    ENDIF
    bt = 0
    DO ifa = 1, Mesh%Nextfaces
      IF (Mesh%boundaryFlag(ifa) .eq. bt) CYCLE
      fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
      IF (fl .eq. 0) CYCLE
#endif
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) 'Boundary type: ', Adjustl(Trim(bc_flag_name(phys%bcflags(fl)))), &
          ' on boundary flagged as: ', Adjustl(Trim(bc_flag_type(fl)))
      ENDIF
      bt = fl
    END DO
    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) " "
    ENDIF

    Mesh%ndir = ndir
    Mesh%ukf = Mesh%Nfaces - ndir

    IF (utils%printint > 0) THEN
      WRITE (6, '(A,I13)') ' Number of different boundaries:     ', nb
      WRITE (6, '(A,I13)') ' Number of Dirichlet faces:          ', ndir
#ifdef PARALL
      WRITE (6, '(A,I13)') ' Number of Ghost faces:              ', Mesh%nghostfaces
#endif
      WRITE (6, *) ' '
    END IF

  contains


    subroutine periodic_faces_preprocess()
      integer*4 :: i,j,ifa,ifan,fl,iel,ieln,ni,ne,nin,nen
      real*8    :: xi,xe,yi,ye,xin,xen,yin,yen
      do i = 1, Mesh%Nextfaces
        fl = Mesh%boundaryFlag(i)
        if (phys%bcflags(fl) == bc_periodic) then
          ! Find corresponding face
          iel = Mesh%extFaces(i,1)
          ifa = Mesh%extFaces(i,2)
          ni = refElPol%face_nodes(ifa,1)
          ne = refElPol%face_nodes(ifa,refElPol%Nfacenodes)
          xi=Mesh%X(Mesh%T(iel,ni),1)
          xe=Mesh%X(Mesh%T(iel,ne),1)
          yi=Mesh%X(Mesh%T(iel,ni),2)
          ye=Mesh%X(Mesh%T(iel,ne),2)
          do j = 1, Mesh%Nextfaces
            if (i==j) cycle
            ieln = Mesh%extFaces(j,1)
            ifan = Mesh%extFaces(j,2)
            nin = refElPol%face_nodes(ifan,1)
            nen = refElPol%face_nodes(ifan,refElPol%Nfacenodes)
            xin=Mesh%X(Mesh%T(ieln,nin),1)
            xen=Mesh%X(Mesh%T(ieln,nen),1)
            yin=Mesh%X(Mesh%T(ieln,nin),2)
            yen=Mesh%X(Mesh%T(ieln,nen),2)
            if ( (xi.ne.xe .and. xi==xen .and. xe==xin) .or. (yi.ne.ye .and. yi==yen .and. ye==yin) ) then
              Mesh%periodic_faces(i)=j
              exit
            endif
          end do
        endif
      end do

    end subroutine periodic_faces_preprocess



  END SUBROUTINE Bc_preprocess

  SUBROUTINE CreateFaceConnectivity()
    USE MPI_OMP
    integer :: ifa, igh
    integer :: infoFace(5), infoFace_ex(2)
    logical :: isdir
    integer :: ifan,ieln,Fi,Fi_per

    igh = 0
    ALLOCATE (Mesh%F(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%Fdir(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%flipFace(Mesh%Nelems, refElPol%Nfaces))
    Mesh%F = 0
    Mesh%Fdir = .false.
    Mesh%flipFace = .false.
    isdir = .false.
    DO ifa = 1, Mesh%Nintfaces
      infoFace = Mesh%intFaces(ifa, :)
      Mesh%F(infoFace(1), infoFace(2)) = ifa
      Mesh%F(infoFace(3), infoFace(4)) = ifa
      IF (infoFace(1) < infoFace(3)) THEN
        Mesh%flipFace(infoFace(3), infoFace(4)) = .true.
      ELSE
        Mesh%flipFace(infoFace(1), infoFace(2)) = .true.
      END IF
#ifdef PARALL
      IF (Mesh%ghostFaces(ifa) .eq. 1) THEN
        igh = igh + 1
      END IF
#endif
    END DO

    DO ifa = 1, Mesh%Nextfaces
      infoFace_ex = Mesh%extFaces(ifa, :)
      IF (Mesh%boundaryFlag(ifa) == 0) THEN
        isdir = .false.
      ELSE
        isdir = (phys%bcflags(Mesh%boundaryFlag(ifa)) == bc_dirichlet)
      ENDIF
      Mesh%F(infoFace_ex(1), infoFace_ex(2)) = ifa + Mesh%Nintfaces
      Mesh%Fdir(infoFace_ex(1), infoFace_ex(2)) = isdir
#ifdef PARALL
      IF (Mesh%ghostFaces(ifa + Mesh%Nintfaces) .eq. 1) THEN
        igh = igh + 1
        IF (Mesh%ghostFlp(igh) .eq. 1) THEN
          Mesh%flipFace(infoFace_ex(1), infoFace_ex(2)) = .true.
        END IF
      END IF
#endif
    END DO
    ! Modify flipface for periodic faces
    DO ifa = 1, Mesh%Nextfaces
      if (Mesh%periodic_faces(ifa).ne.0) then
        infoFace_ex = Mesh%extFaces(ifa, :)
        Fi=Mesh%F(infoFace_ex(1), infoFace_ex(2))
        ieln = Mesh%extfaces(Mesh%periodic_faces(ifa),1)
        ifan = Mesh%extfaces(Mesh%periodic_faces(ifa),2)
        Fi_per = Mesh%F(ieln,ifan)
        if (Fi>Fi_per) then
          Mesh%flipFace(infoFace_ex(1), infoFace_ex(2))=.true.
        endif

      endif
    END DO

  END SUBROUTINE CreateFaceConnectivity

  !********************************************
  ! Compute the element size in the
  ! mesh (the length of the edge 2D-3D is
  ! considered)
  !********************************************
  SUBROUTINE computeElementSize()
    integer                         :: i, ierr
    real*8                          :: h1, h2, h3, h4
    real*8, dimension(Mesh%ndim)     :: p1, p2, p3, p4

    ALLOCATE (Mesh%elemSize(Mesh%Nelems))
    Mesh%elemSize = 0.

    ! Loop in elements
    IF (refElPol%elemType .eq. 0) THEN
      !$OMP PARALLEL PRIVATE(i,p1,p2,p3,h1,h2,h3)
      !$OMP DO
      DO i = 1, Mesh%Nelems
        p1 = Mesh%X(Mesh%Tlin(i, 1), :)
        p2 = Mesh%X(Mesh%Tlin(i, 2), :)
        p3 = Mesh%X(Mesh%Tlin(i, 3), :)
        h1 = norm2(p1 - p2)
        h2 = norm2(p1 - p3)
        h3 = norm2(p3 - p2)
        Mesh%elemSize(i) = min(h1, h2, h3)
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    ELSEIF (refElPol%elemType .eq. 1) THEN
      !$OMP PARALLEL PRIVATE(i,p1,p2,p3,p4,h1,h2,h3,h4)
      !$OMP DO
      DO i = 1, Mesh%Nelems
        p1 = Mesh%X(Mesh%Tlin(i, 1), :)
        p2 = Mesh%X(Mesh%Tlin(i, 2), :)
        p3 = Mesh%X(Mesh%Tlin(i, 3), :)
        p4 = Mesh%X(Mesh%Tlin(i, 4), :)
        h1 = norm2(p1 - p2)
        h2 = norm2(p1 - p3)
        h3 = norm2(p3 - p2)
        h4 = norm2(p4 - p3)
        Mesh%elemSize(i) = min(h1, h2, h3, h4)
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
    END IF

  END SUBROUTINE computeElementSize
  
  
  SUBROUTINE computePuffArea()
  real*8   :: Xf(refElPol%Nfacenodes,2),xyg(refElPol%NGauss1D,2),xyg_d(refElPol%NGauss1D,2),dline
  integer  :: i,el,fa,fl,g
  real*8   :: xyDerNorm_g
  
  
  Mesh%puff_area = 0.
  
		DO i = 1, Mesh%Nextfaces
		
			  fl = Mesh%boundaryFlag(i) 
			  
					IF (phys%bcflags(fl) .ne. bc_BohmPuff) THEN
					  CYCLE
					END IF
					
					el = Mesh%extfaces(i,1)
					fa = Mesh%extfaces(i,2)
					Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
					xyg = matmul(refElPol%N1D,Xf)
					xyg_d = matmul(refElPol%Nxi1D,Xf)
					DO g = 1, refElPol%NGauss1D
					
						  xyDerNorm_g = norm2(xyg_d(g,:))
						  dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
						  dline = dline*xyg(g,1)
						  Mesh%puff_area = Mesh%puff_area + 2*pi*dline	
					END DO
		END DO

  END SUBROUTINE computePuffArea
  
  
  SUBROUTINE allocate_load_JTOR_PUFF
    IF (switch%testcase .eq. 54) THEN
      CALL loadJtorMap()
#ifdef NEUTRAL
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE(6,*) 'Puff is analytical'
      ENDIF
#else
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE(6,*) 'This test case should be ran with neutrals'
        STOP
      ENDIF
#endif
    ENDIF
    IF (switch%testcase .eq. 59) THEN
      CALL loadJtorMap()
#ifdef NEUTRAL
      IF(switch%time_init) THEN
        IF (MPIvar%glob_id .eq. 0) THEN
          WRITE(6,*) 'Puff is analytical'
        ENDIF
      ELSE
        IF (MPIvar%glob_id .eq. 0) THEN
          WRITE(6,*) 'Puff is experimental'
        ENDIF
          CALL SetPuff()
      ENDIF
#else
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE(6,*) 'This test case should be ran with neutrals'
      ENDIF
      STOP
#endif
    ENDIF

  END SUBROUTINE allocate_load_JTOR_PUFF

  ! Below are routines from Manuel MHDG v2.1. Copy as it without any check: TODO adapt it to global magnetic field

  SUBROUTINE loadJtorMap()
    USE interpolation
    USE HDF5
    USE HDF5_io_module
    !USE MPI_OMP
    integer        :: i,ierr,ip,jp,ind,j,k
    integer(HID_T) :: file_id

    character(LEN=1000)    :: fname = 'Jtor_exp/WEST_54487_Jtor'
    character(70)        :: npr,nid,nit

    real*8,pointer,dimension(:,:) :: r2D,z2D,Jtor
    real*8,allocatable,dimension(:)   :: xvec,yvec
    real*8                            :: x,y
    real*8                            :: Jtor_temp

    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE(6,*) "******* Loading Toroidal Current *******"
    ENDIF
#ifdef TOR3D
    ! Allocate storing space in phys
    ALLOCATE(phys%Jtor(Mesh%Nnodes*Mesh%Nnodes_toroidal))
#else
    ALLOCATE(phys%Jtor(Mesh%Nnodes))
#endif
    phys%Jtor = 0.

    ! Dimensions of the file storing the magnetic field for West
    !ip = 200
    !jp = 200
    ip = 457;
    jp = 457;
    ALLOCATE(r2D(ip,jp))
    ALLOCATE(z2D(ip,jp))
    ALLOCATE(Jtor(ip,jp))

    ! Read file
    IF(switch%testcase .eq. 54) THEN
        fname = trim(adjustl(fname)) //'_0099.h5'
    ELSE
      ! File name
        write(nit, "(i10)") time%it
        nit = trim(adjustl(nit))
        k = INDEX(nit, " ") -1
        fname = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
    ENDIF
    IF (MPIvar%glob_id .eq. 0) THEN
      write(6,*) 'Toroidal current loaded from file: ', trim(adjustl(fname))
    ENDIF

    CALL HDF5_open(fname,file_id,IERR)

    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
    CALL HDF5_array2D_reading(file_id,Jtor,'Jtor')
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Interpolate
    ALLOCATE(xvec(jp))
    ALLOCATE(yvec(ip))
    xvec = r2D(1,:)
    yvec = z2D(:,1)

    DO i = 1,Mesh%Nnodes
      x = Mesh%X(i,1)
      y = Mesh%X(i,2)
      ind = i
      Jtor_temp = interpolate(ip, yvec,jp, xvec,Jtor, y,x, 1e-12)
#ifdef TOR3D
      DO j = 1, Mesh%Nnodes_toroidal
        ind = (j - 1)*Mesh%Nnodes + i
#endif
        phys%Jtor(ind) = Jtor_temp
#ifdef TOR3D
      END DO
#endif
    END DO

    ! Free memory
    DEALLOCATE(r2D,z2D,Jtor,xvec,yvec)

  END SUBROUTINE loadJtorMap

  SUBROUTINE SetPuff()
    USE HDF5
    USE HDF5_io_module
    integer        :: ierr
    character(LEN=25) :: fname = 'Puff_54487.h5'
    integer(HID_T) :: file_id

    ! Allocate storing space in phys (puff for WEST, 403 entries)
    ALLOCATE(phys%puff_exp(403))
    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) "******* Loading puff *******"
    ENDIF

    ! Read file
    CALL HDF5_open(fname,file_id,IERR)
    CALL HDF5_array1D_reading(file_id,phys%puff_exp,'puff')
    IF (MPIvar%glob_id .eq. 0) THEN
      write(6,*) 'Puff loaded from file: ', trim(adjustl(fname))
    ENDIF
    CALL HDF5_close(file_id)

  END SUBROUTINE SetPuff

!   SUBROUTINE computePuffDensityWest()
!     USE globals
!     integer               :: ierr
!     integer*4             :: g, iel, N2D, counter
!     real*8,allocatable    :: Xel(:,:)
!     real*8                :: xy(refElPol%Ngauss2d,2)
!
!     allocate(Xel(Mesh%Nnodesperelem,2))
!     N2D = Mesh%Nelems
!     counter = 0
!
!     DO iel = 1,N2D
!       ! Coordinates of the nodes of the element
!       Xel = Mesh%X(Mesh%T(iel,:),:)
!       ! Gauss points position
!       xy = matmul(refElPol%N2D,Xel)
!         ! Compute number of Gauss points per element in the box (don't consider ghost faces)
! #ifdef PARALL
!       IF (Mesh%ghostElems(iel) .eq. 0) THEN
! #endif
!         DO g = 1,refElPol%Ngauss2d
!           IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!             counter = counter + 1
!           ENDIF
!         ENDDO
! #ifdef PARALL
!       ENDIF
! #endif
!     ENDDO
!
! #ifdef PARALL
!         CALL MPI_ALLREDUCE(MPI_IN_PLACE, counter, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
! #endif
!   IF (counter .eq. 0) THEN
!     stop "PROBLEM: no elements are present in the buffer. STOPPING."
!   ELSE
!     phys%puff_density = phys%puff/counter
!
!     IF((switch%testcase .eq. 59) .and. (.not. switch%time_init))  THEN
!       allocate (phys%puff_exp_density(size(phys%puff_exp)))
!       phys%puff_exp_density(time%it+1) = phys%puff_exp(time%it+1)/counter
!     ENDIF
!   ENDIF
! END SUBROUTINE computePuffDensityWest


  
  SUBROUTINE computeCoreArea()
  real*8   :: Xf(refElPol%Nfacenodes,2),xyg(refElPol%NGauss1D,2),xyg_d(refElPol%NGauss1D,2),dline
  integer  :: i,el,fa,fl,g
  real*8   :: xyDerNorm_g
  
  
  Mesh%core_area = 0.
  
		DO i = 1, Mesh%Nextfaces
		
			  fl = Mesh%boundaryFlag(i) 
			  
					IF (phys%bcflags(fl) > 10) THEN
					  CYCLE
					END IF
					
					el = Mesh%extfaces(i,1)
					fa = Mesh%extfaces(i,2)
					Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
					xyg = matmul(refElPol%N1D,Xf)
					xyg_d = matmul(refElPol%Nxi1D,Xf)
					DO g = 1, refElPol%NGauss1D
					
						  xyDerNorm_g = norm2(xyg_d(g,:))
						  dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
						  dline = dline*xyg(g,1)
						  Mesh%core_area = Mesh%core_area + 2*pi*dline
					
					END DO
					
		END DO

  END SUBROUTINE computeCoreArea
END MODULE preprocess
