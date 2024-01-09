!*****************************************
! project: MHDG
! file: HDG_precalculatedfirstequation.f90
! date: 02/01/2020
! Generate the matrices that remain constant
! during the whole computation
!*****************************************

SUBROUTINE HDG_precalculatedfirstequation()
  USE globals
  USE analytical
  USE physics
  USE LinearAlgebra
  USE printUtils
  USE MPI_OMP
  USE Debug

  IMPLICIT NONE

  !*******************************************************************************
  !
  !    PRECALCULATED MATRICES FOR THE FIRST EQUATION (DEFINITION OF THE GRADIENT)
  !
  !*******************************************************************************
  integer*4             :: Ndim,Neq,N2D,Npel,Npfl,Ngfl,Ngvo
  integer*4             :: iel,ifa,i,j,iface
  integer*4             :: els(2),fas(2)
  real*8,allocatable    :: Xel(:,:),Xfl(:,:)
#ifdef TOR3D
  integer*4             :: itor,itorg,iel3
  integer*4             :: ntorloc
  integer*4             :: Np2D,Np1dpol,Np1dtor,Npfp,Ng2D,Ng1dpol,Ng1dtor,Ngfp
  real*8                :: tdiv(numer%ntor + 1)
  real*8                :: htor,tel(refElTor%Nnodes1d)
#endif

  integer :: OMP_GET_THREAD_NUM

  IF (MPIvar%glob_id .eq. 0) THEN
    IF (utils%printint > 0) THEN
      WRITE (6,*) '*************************************************'
      WRITE (6,*) '*           PRECALCULATED MATRICES              *'
      WRITE (6,*) '*************************************************'
    END IF
  END IF
  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1,timing%clock_rate1)
  end if

#ifdef TOR3D
  !*************************************************************
  !               3D stuff
  !*************************************************************
  Ndim = 3                         ! Number of dimensions
  Neq = phys%Neq                  ! Number of equations
  N2D = Mesh%Nelems               ! Number of 2D elements
  Np2D = refElPol%Nnodes2D         ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D         ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D         ! Number of nodes in the 1D toroidal segments
  Npel = Np2D*refElTor%Nnodes1D    ! Number of nodes of each element
  Npfl = Np1Dpol*Np1Dtor           ! Number of nodes of each lateral face
  Npfp = Np2D                      ! Number of nodes of each poloidal face
  Ng2D = refElPol%Ngauss2d         ! Number of Gauss points in the 2D element
  Ng1Dpol = refElPol%Ngauss1d         ! Number of Gauss points in the 1D poloidal segments
  Ng1Dtor = refEltor%Ngauss1d         ! Number of Gauss points in the 1D toroidal segments
  Ngvo = Ng2D*Ng1dtor              ! Number of Gauss points for volume computations
  Ngfl = Ng1dtor*Ng1dpol           ! Number of Gauss points for toroidal faces computations
  Ngfp = Ng2D                      ! Number of Gauss points for poloidal faces computations

  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i = 1,numer%ntor
    tdiv(i + 1) = i*htor
  END DO
#ifdef PARALL
  IF (MPIvar%ntor .gt. 1) THEN
    ntorloc = numer%ntor/MPIvar%ntor + 1
  ELSE
    ntorloc = numer%ntor
  ENDIF
#else
  ntorloc = numer%ntor
#endif


  !*****************
  ! Loop in elements
  !*****************
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,iel3,ifa,iface,itor,itorg,tel,Xel,Xfl)
  allocate(Xel(Mesh%Nnodesperelem,2))
  allocate(Xfl(refElPol%Nfacenodes,2))
  !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO itor = 1,ntorloc
    DO iel = 1,N2D
      ! I made a perfectly nested loop for enabling omp parallelization
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))

      ! Index of 3D element
      iel3 = (itor - 1)*N2d+iel

      ! Coordinates of the nodes of the element
      Xel = Mesh%X(Mesh%T(iel,:),:)

      ! Compute the matrices for the element
      CALL elemental_matrices_volume(iel3,Xel,tel)

      ! First poloidal face
      ifa = 1
      CALL elemental_matrices_pol_faces(iel3,ifa,Xel)

      ! Toroidal faces
      DO ifa=1,refElPol%Nfaces
        IF (Mesh%Fdir(iel,ifa)) CYCLE
        iface = Mesh%F(iel,ifa)
        Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)
        if (iface.le.Mesh%Nintfaces) then
          call elemental_matrices_int_faces(iel3,ifa+1,Xfl,tel)
        else
          call elemental_matrices_ext_faces(iel3,ifa+1,Xfl,tel)
        endif
      END DO

      ! Second poloidal face
      ifa = refElPol%Nfaces + 2
      CALL elemental_matrices_pol_faces(iel3,ifa,Xel)

    END DO
  END DO
  !$OMP END DO
  deallocate(Xel,Xfl)
  !$OMP END PARALLEL



  !****************************
  ! Loop in poloidal faces
  !***************** ***********
  !allocate(Xfl(Mesh%Nnodesperelem,2))
  !   DO itor = 1,ntorloc
  !   DO iface = 1,N2D


  !      els(1) = (itor - 1)*N2D+iface
  !      fas(1) = 1
  !      els(2) = (itor - 2)*N2D+iface
  !      if (itor == 1 .and. MPIvar%ntor == 1) then
  !         els(2) = N2D*(ntorloc - 1) + iface
  !      endif
  !      fas(2) = refElPol%Nfaces + 2

  !      ! Coordinates of the nodes of the face
  !      Xfl = Mesh%X(Mesh%T(iface,:),:)

  !      ! Compute the matrices for the element
  !      CALL elemental_matrices_pol_faces(els,fas,Xfl)

  !   END DO
  !   END DO
  !deallocate(Xfl)








  !**********************************
  ! Loop in toroidal interior faces
  !**********************************
  !allocate(Xfl(refElPol%Nfacenodes,2))
  !   DO itor = 1,ntorloc
  !#ifdef PARALL
  !      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
  !      if (itorg == numer%ntor + 1) itorg = 1
  !#else
  !      itorg = itor
  !#endif
  !      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))

  !      DO iface = 1,Mesh%Nintfaces

  !         els = Mesh%intfaces(iFace,(/1,3/))
  !         fas = Mesh%intfaces(iFace,(/2,4/))

  !         ! Coordinates of the nodes of the face
  !         Xfl = Mesh%X(Mesh%T(els(1),refElPol%face_nodes(fas(1),:)),:)

  !         ! Compute the matrices for the element
  !         CALL elemental_matrices_int_faces((itor - 1)*N2d+els,fas + 1,Xfl,tel)

  !      END DO
  !   END DO
  !deallocate(Xfl)


  !*********************************
  ! Loop in toroidal exterior faces
  !*********************************
  !allocate(Xfl(refElPol%Nfacenodes,2))
  !   DO itor = 1,ntorloc
  !#ifdef PARALL
  !      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
  !      if (itorg == numer%ntor + 1) itorg = 1
  !#else
  !      itorg = itor
  !#endif
  !      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))

  !      DO iface = 1,Mesh%Nextfaces

  !         iel = Mesh%extfaces(iFace,1)
  !         ifa = Mesh%extfaces(iFace,2)

  !         IF (Mesh%Fdir(iel,ifa)) CYCLE

  !         ! Index of 3D element
  !         iel3 = (itor - 1)*N2d+iel

  !         ! Coordinates of the nodes of the element
  !         Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

  !         ! Compute the matrices for the element
  !         CALL elemental_matrices_ext_faces(iel3,ifa + 1,Xfl,tel)
  !      END DO
  !   END DO
  !deallocate(Xfl)


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1,timing%clock_rate1)
    timing%runtpre = timing%runtpre + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputpre = timing%cputpre + timing%tpe1 - timing%tps1
  end if


CONTAINS

  !*****************************************
  ! Volume computations in 3D
  !*****************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel,tel)
    integer*4,intent(IN)  :: iel
    real*8,intent(IN)     :: Xel(:,:),tel(:)
    integer*4             :: g,NGaussPol,NGaussTor,igtor,igpol,i
    real*8                :: dvolu,dvolu1d,htor
    real*8                :: xy(Ng2D,2),teg(Ng1Dtor)
    real*8                :: J11(Ng2D),J12(Ng2D)
    real*8                :: J21(Ng2D),J22(Ng2D)
    real*8                :: iJ11(Ng2D),iJ12(Ng2D)
    real*8                :: iJ21(Ng2D),iJ22(Ng2D)
    real*8                :: detJ(Ng2D)
    real*8                :: Nxg(Np2D),Nyg(Np2D),Nx_ax(Np2D)
    real*8,pointer        :: N1g(:),N2g(:)
    real*8                :: N1xg(Np1Dtor),N1xg_cart(Np1Dtor)
    real*8                :: t_g(2),n_g(2)
    real*8                :: NN(Npel,Npel)
    real*8                :: NxNy_ax(Npel,Npel,Ndim)
    real*8                :: Ni(Npel),Nidvolu(Npel)
    real*8,allocatable    :: Aqq(:,:),Aqu(:,:,:)
    integer*4             :: ind_ass(Npel),ind_asq(Npel)

    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

    allocate(Aqq(Npel,Npel))
    allocate(Aqu(Npel,Npel,Ndim))
    Aqq = 0.
    Aqu = 0.

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Gauss points position
    xy = matmul(refElPol%N2D,Xel)

    ! Gauss points position in the toroidal direction
    teg = matmul(refElTor%N1d,tel)

    ! Loop in 2D Gauss points
    J11 = matmul(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
    J12 = matmul(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
    J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
    J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ

    NgaussPol = Ng2D
    NgaussTor = Ng1Dtor
    DO igtor = 1,NGaussTor
      N1g => refElTor%N1D(igtor,:)                   ! Toroidal shape function
      N1xg_cart = refElTor%Nxi1D(igtor,:)*2/htor       ! Toroidal shape function derivative
      dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume

      DO igpol = 1,NgaussPol
        g = (igtor - 1)*NGaussPol + igpol

        ! Poloidal shape functions and derivatives
        N2g => refElPol%N2D(igpol,:)
        Nxg = iJ11(igpol)*refElPol%Nxi2D(igpol,:) + iJ12(igpol)*refElPol%Neta2D(igpol,:)
        Nyg = iJ21(igpol)*refElPol%Nxi2D(igpol,:) + iJ22(igpol)*refElPol%Neta2D(igpol,:)

        IF (switch%axisym) THEN
          Nx_ax = Nxg + 1./xy(igpol,1)*N2g
          N1xg = N1xg_cart/xy(igpol,1)
        ELSE
          Nx_ax = Nxg
          N1xg = N1xg_cart
        END IF

        ! 3D integration weight
        dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(igpol,1)
        END IF

        ! 3D shape functions
        Ni = col(TensorProduct(N2g,N1g))      ! 3D shape function
        Nidvolu = Ni*dvolu
        NN = tensorProduct(Ni,Nidvolu)
        NxNy_ax(:,:,1) = TensorProduct(col(TensorProduct(Nx_ax,N1g)),Nidvolu)
        NxNy_ax(:,:,2) = TensorProduct(col(TensorProduct(Nyg,N1g)),Nidvolu)
        NxNy_ax(:,:,3) = TensorProduct(col(TensorProduct(N2g,N1xg)),Nidvolu)

        ! Local assembly
        CALL assemblyVolumeContribution(NxNy_ax,NN,Aqq,Aqu)
      END DO
    END DO

    CALL do_assembly(Aqq,Aqu,ind_ass,ind_asq,iel)
    DEALLOCATE (Aqq,Aqu)
    NULLIFY (N1g,N2g)

  END SUBROUTINE elemental_matrices_volume

  !*****************************************
  ! Poloidal faces computations
  !*****************************************
  SUBROUTINE elemental_matrices_pol_faces(iel,ifa,Xfp)
    integer*4,intent(IN)  :: iel,ifa
    real*8,intent(IN)     :: Xfp(:,:)
    integer*4             :: g,NGauss,i,j
    real*8                :: dsurf(Ng2d)
    real*8                :: xy(Ng2d,2)
    real*8                :: J11(Ng2d),J12(Ng2d)
    real*8                :: J21(Ng2d),J22(Ng2d)
    real*8                :: detJ(Ng2d)
    real*8                :: iJ11(Ng2d),iJ12(Ng2d)
    real*8                :: iJ21(Ng2d),iJ22(Ng2d)
    real*8                :: NiNi(Np2D,Np2D)
    real*8                :: n_g(Ng2d,3)
    integer*4             :: ind_asf(Np2D),ind_ash(Np2D)
    integer*4             :: ind_ff(neq*Np2D),ind_fG(neq*Np2D*Ndim)
    real*8,parameter     :: tol = 1e-12

    ind_asf = (/(i,i=0,Neq*(Np2D-1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Np2D-1)*Ndim,Neq*Ndim)/)

    ! Gauss points position
    xy = matmul(refElPol%N2D,Xfp)

    ! Loop in 2D Gauss points
    Ngauss = Ng2d

    J11 = matmul(refElPol%Nxi2D,Xfp(:,1))                           ! ng x 1
    J12 = matmul(refElPol%Nxi2D,Xfp(:,2))                           ! ng x 1
    J21 = matmul(refElPol%Neta2D,Xfp(:,1))                          ! ng x 1
    J22 = matmul(refElPol%Neta2D,Xfp(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ
    dsurf = refElPol%gauss_weights2D*detJ

    !      DO lel = 1,2

    !         iel = els(lel)
    !         if (iel < 1) CYCLE
    !         ifa = fas(lel)
    IF (ifa == 1) THEN
      ind_ff = (/(i,i=1,Np2D*Neq)/)
      ind_fg = (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = -1
    ELSE
      ind_ff = Np2D*Neq + refElPol%Nfaces*Npfl*Neq + (/(i,i=1,Np2D*Neq)/)
      ind_fg = Np2d*(Np1dTor - 1)*Ndim*Neq + (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = 1
    ENDIF
    DO g = 1,NGauss

      ! Shape functions product
      NiNi = TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dsurf(g)

      ! Local assembly
      CALL assemblyFacesContribution(iel,NiNi,n_g(g,:),ind_asf,ind_ash,ind_fG,ind_ff)
    END DO
    !      END DO

  END SUBROUTINE elemental_matrices_pol_faces

  !*****************************************
  !  Toroidal interior faces computations
  !*****************************************
  SUBROUTINE elemental_matrices_int_faces(iel,ifa,Xfl,tel)
    integer*4,intent(IN)  :: iel,ifa
    real*8,intent(IN)     :: Xfl(:,:),tel(:)
    integer*4             :: g,igtor,igpol,i,j
    real*8                :: xyf(Ng1Dpol,2)
    real*8                :: xyDer(Ng1Dpol,2)
    integer*4             :: ind_ff(neq*Npfl),ind_fG(neq*Npfl*Ndim)
    real*8                :: t_g(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    real*8                :: n_g(Ngfl,3),dsurf(Ngfl),dsurfg
    real*8                :: NiNi(Npfl,Npfl)
    real*8,pointer        :: Nfg(:)
    integer*4             :: ind(Ng1dPol),ind_asf(Npfl),ind_ash(Npfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Local assembly indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Coordinates and derivatives at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)
    ! Shape function derivatives at Gauss points
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    ! Compute dsurf
    xydNorm_g = sqrt(xyDer(:,1)**2 + xyDer(:,2)**2)
    dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))

    ! Compute exterior normal
    t_g(:,1) = xyDer(:,1)/xydNorm_g
    t_g(:,2) = xyDer(:,2)/xydNorm_g
    n_g = 0.
    DO i = 1,Ng1dTor
      ind = (i - 1)*Ng1dPol + (/(j,j=1,Ng1dPol)/)
      n_g(ind,1) = t_g(:,2)
      n_g(ind,2) = -t_g(:,1)
    END DO

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,Ng1dTor
      DO igpol = 1,Ng1dPol

        g = (igtor - 1)*Ng1dPol + igpol

        ! Face shape functions
        Nfg => refElTor%sFTF(g,:)

        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

        ! Shape functions product
        NiNi = tensorProduct(Nfg,Nfg)*dsurfg

        ! Local assembly
        CALL assemblyFacesContribution(iel,NiNi,n_g(g,:),ind_asf,ind_ash,ind_fG,ind_ff)
      END DO ! Gauss points
    END DO

  END SUBROUTINE elemental_matrices_int_faces

  !*****************************************
  ! Toroidal exterior faces computations
  !*****************************************
  SUBROUTINE elemental_matrices_ext_faces(iel,ifa,Xfl,tel)
    integer*4,intent(IN)  :: iel,ifa
    real*8,intent(IN)     :: Xfl(:,:),tel(:)
    integer*4             :: g,igtor,igpol,i,j
    real*8                :: xyf(Ng1Dpol,2)
    real*8                :: xyDer(Ng1Dpol,2)
    integer*4             :: ind_ff(neq*Npfl),ind_fG(neq*Npfl*Ndim)
    real*8                :: t_g(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    real*8                :: n_g(Ngfl,3),dsurf(Ngfl),dsurfg
    real*8                :: NiNi(Npfl,Npfl)
    real*8,pointer        :: Nfg(:)
    integer*4             :: ind(Ng1dPol),ind_asf(Npfl),ind_ash(Npfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Local assembly indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Coordinates and derivatives at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    ! Compute dsurf
    xydNorm_g = sqrt(xyDer(:,1)**2 + xyDer(:,2)**2)
    dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))

    ! Compute exterior normal
    t_g(:,1) = xyDer(:,1)/xydNorm_g
    t_g(:,2) = xyDer(:,2)/xydNorm_g
    n_g = 0.
    DO i = 1,Ng1dTor
      ind = (i - 1)*Ng1dPol + (/(j,j=1,Ng1dPol)/)
      n_g(ind,1) = t_g(:,2)
      n_g(ind,2) = -t_g(:,1)
    END DO

    !*****************************
    ! Loop in face Gauss points
    !*****************************
    DO igtor = 1,Ng1dTor
      DO igpol = 1,Ng1dPol

        g = (igtor - 1)*Ng1dPol + igpol

        ! Face shape functions
        Nfg => refElTor%sFTF(g,:)

        IF (switch%axisym) THEN
          dsurfg = dsurf(g)*xyf(igpol,1)
        ELSE
          dsurfg = dsurf(g)
        END IF

        ! Shape functions product
        NiNi = tensorProduct(Nfg,Nfg)*dsurfg

        ! Local assembly
        CALL assemblyFacesContribution(iel,NiNi,n_g(g,:),ind_asf,ind_ash,ind_fG,ind_ff)
      END DO ! Gauss points
    END DO

  END SUBROUTINE elemental_matrices_ext_faces

#else

  !*************************************************************
  !               2D stuff
  !*************************************************************
  Ndim = 2                         ! Number of dimensions
  Neq = phys%Neq                  ! Number of equations
  N2D = Mesh%Nelems               ! Number of 2D elements
  Npel = refElPol%Nnodes2D         ! Number of nodes of each element
  Npfl = refElPol%Nfacenodes       ! Number of nodes of each lateral face
  Ngvo = refElPol%Ngauss2d         ! Number of Gauss points for volume computations
  Ngfl = refElPol%Ngauss1d         ! Number of Gauss points for faces computations

  !*****************
  ! Loop in elements
  !*****************
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,ifa,iface,Xel,Xfl)
  allocate(Xel(Mesh%Nnodesperelem,2))
  allocate(Xfl(refElPol%Nfacenodes,2))
  !$OMP DO SCHEDULE(STATIC)
  DO iel = 1,N2D

    ! Coordinates of the nodes of the element
    Xel = Mesh%X(Mesh%T(iel,:),:)

    ! Compute the matrices for the element
    CALL elemental_matrices_volume(iel,Xel)

    ! Loop in local faces
    DO ifa=1,refElPol%Nfaces
      IF (Mesh%Fdir(iel,ifa)) CYCLE
      iface = Mesh%F(iel,ifa)
      Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)
      if (iface.le.Mesh%Nintfaces) then
        call elemental_matrices_int_faces(iel,ifa,Xfl)
      else
        if (Mesh%periodic_faces(iface-Mesh%Nintfaces).eq.0) then
          call elemental_matrices_ext_faces(iel,ifa,Xfl)
        else
          ! periodic face
          call elemental_matrices_int_faces(iel,ifa,Xfl)
        endif
      endif
    END DO

  END DO
  !$OMP END DO
  deallocate(Xel,Xfl)
  !$OMP END PARALLEL

  IF (MPIvar%glob_id .eq. 0) THEN
    IF (utils%printint > 0) THEN
      WRITE (6,*) "Done!"
    END IF
  END IF


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1,timing%clock_rate1)
    timing%runtpre = timing%runtpre + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputpre = timing%cputpre + timing%tpe1 - timing%tps1
  end if

CONTAINS

  !*****************************************
  ! Volume computations in 2D
  !*****************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel)
    integer*4,intent(IN)  :: iel
    real*8,intent(IN)     :: Xel(:,:)
    integer*4             :: g,NGauss,i,j,k,idm
    real*8                :: dvolu
    real*8                :: xy(Ngvo,2)
    real*8                :: J11(Ngvo),J12(Ngvo)
    real*8                :: J21(Ngvo),J22(Ngvo)
    real*8                :: detJ(Ngvo)
    real*8                :: iJ11(Ngvo),iJ12(Ngvo)
    real*8                :: iJ21(Ngvo),iJ22(Ngvo)
    real*8                :: Nxg(Npel),Nyg(Npel),Nx_ax(Npel)
    real*8                :: t_g(2),n_g(2)
    real*8                :: NN(Npel,Npel)
    real*8                :: NxNy_ax(Npel,Npel,Ndim)
    real*8,allocatable    :: Aqq(:,:),Aqu(:,:,:)
    integer*4             :: ind_ass(Npel),ind_asq(Npel)
    real*8,parameter     :: tol = 1e-12

    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

    ! Gauss points position
    xy = matmul(refElPol%N2D,Xel)

    allocate(Aqq(Npel,Npel))
    allocate(Aqu(Npel,Npel,Ndim))
    Aqq = 0.
    Aqu = 0.

    ! Loop in 2D Gauss points
    Ngauss = Ngvo
    J11 = matmul(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
    J12 = matmul(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
    J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
    J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ
    DO g = 1,NGauss

      IF (detJ(g) < tol) THEN
        write(6,*) "iel: ", iel, "detJ(g): ", detJ(g)
        error stop "Negative jacobian in element: "
      END if

      ! x and y derivatives of the shape functions
      Nxg = iJ11(g)*refElPol%Nxi2D(g,:) + iJ12(g)*refElPol%Neta2D(g,:)
      Nyg = iJ21(g)*refElPol%Nxi2D(g,:) + iJ22(g)*refElPol%Neta2D(g,:)
      IF (switch%axisym) THEN
        Nx_ax = Nxg + 1./xy(g,1)*refElPol%N2D(g,:)
      ELSE
        Nx_ax = Nxg
      END IF
      ! Integration weight
      dvolu = refElPol%gauss_weights2D(g)*detJ(g)
      IF (switch%axisym) THEN
        dvolu = dvolu*xy(g,1)
      END IF

      NN = TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu
      NxNy_ax(:,:,1) = TensorProduct(Nx_ax,refElPol%N2D(g,:))*dvolu
      NxNy_ax(:,:,2) = TensorProduct(Nyg,refElPol%N2D(g,:))*dvolu

      ! Local assembly
      CALL assemblyVolumeContribution(NxNy_ax,NN,Aqq,Aqu)
    END DO
    CALL do_assembly(Aqq,Aqu,ind_ass,ind_asq,iel)

    DEALLOCATE (Aqq,Aqu)

  END SUBROUTINE elemental_matrices_volume

  !*****************************************
  ! Interior faces computations
  !*****************************************
  SUBROUTINE elemental_matrices_int_faces(iel,ifa,Xfl)
    integer*4,intent(IN)  :: iel,ifa
    real*8,intent(IN)     :: Xfl(:,:)
    integer*4             :: g,NGauss,i,j,k,idm
    real*8                :: dline,xyDerNorm_g
    real*8                :: xyf(Ngfl,2)
    real*8                :: xyDer(Ngfl,2)
    integer*4             :: ind_ff(neq*Npfl),ind_fG(neq*Npfl*Ndim)
    real*8                :: t_g(2),n_g(2)
    real*8                :: NiNi(Npfl,Npfl)
    integer*4             :: ind_asf(Npfl),ind_ash(Npfl)
    !***********************************
    !    Volume computation
    !***********************************

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    ! Faces computations
    !***********************************
    NGauss = Ngfl

    ind_fG = reshape(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*ndim*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)

    ! Coordinates and derivatives at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)

    ! Shape function derivatives at Gauss points
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    ! Loop in 1D Gauss points
    DO g = 1,NGauss

      ! Calculate the integration weight
      xyDerNorm_g = norm2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

      IF (switch%axisym) THEN
        dline = dline*xyf(g,1)
      END IF
      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]

      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      CALL assemblyFacesContribution(iel,NiNi,n_g,ind_asf,ind_ash,ind_fG,ind_ff)
    END DO ! Gauss points

  END SUBROUTINE elemental_matrices_int_faces

  !*****************************************
  ! Exterior faces computations
  !*****************************************
  SUBROUTINE elemental_matrices_ext_faces(iel,ifa,Xfl)
    integer*4,intent(IN)  :: iel,ifa
    real*8,intent(IN)     :: Xfl(:,:)
    integer*4             :: g,NGauss,i,j,k,idm
    real*8                :: dline,xyDerNorm_g
    real*8                :: xyf(Ngfl,2)
    real*8                :: xyDer(Ngfl,2)
    integer*4             :: ind_ff(neq*Npfl),ind_fG(neq*Ndim*Npfl)
    real*8                :: t_g(2),n_g(2)
    real*8                :: NiNi(Npfl,Npfl)
    integer*4             :: ind_asf(Npfl),ind_ash(Npfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    NGauss = Ngfl

    ! Assembly indices
    ind_fG = reshape(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*ndim*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)

    ! Coordinates and derivatives at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    ! Loop in 1D Gauss points
    DO g = 1,NGauss

      ! Calculate the integration weight
      xyDerNorm_g = norm2(xyDer(g,:))
      dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

      IF (switch%axisym) THEN
        dline = dline*xyf(g,1)
      END IF
      ! Unit normal to the boundary
      t_g = xyDer(g,:)/xyDerNorm_g
      n_g = [t_g(2),-t_g(1)]

      NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      CALL assemblyFacesContribution(iel,NiNi,n_g,ind_asf,ind_ash,ind_fG,ind_ff)
    END DO ! Gauss points

  END SUBROUTINE elemental_matrices_ext_faces

#endif

  SUBROUTINE assemblyVolumeContribution(NxNy_ax,NN,Aqq,Aqu)
    real*8                  :: NN(:,:)
    real*8                  :: NxNy_ax(:,:,:)
    real*8,intent(inout)    :: Aqq(:,:),Aqu(:,:,:)
    integer*4   :: i,k,idm
    integer*4   :: ind(Npel),ind_1(Npel),ind_2(Npel)
    Aqq = Aqq + NN
    DO k = 1,Ndim
      Aqu(:,:,k) = Aqu(:,:,k) +  NxNy_ax(:,:,k)
    END DO
  END SUBROUTINE assemblyVolumeContribution

  SUBROUTINE do_assembly(Aqq,Aqu,ind_ass,ind_asq,iel)
    real*8,intent(in)    :: Aqq(:,:),Aqu(:,:,:)
    integer*4,intent(in) :: iel,ind_ass(:),ind_asq(:)
    integer :: i,j,k,z
    integer*4,dimension(Npel) :: ind_i,ind_j,ind_k
    real*8,allocatable :: iAqq(:,:),exAqq(:,:)
    allocate(iAqq(Neq*Ndim*Npel,Neq*Ndim*Npel))
    allocate(exAqq(Neq*Ndim*Npel,Neq*Ndim*Npel))
    iAqq = 0.
    exAqq = 0.
    DO j = 1,Neq
      ind_j = j + ind_ass
      DO k = 1,Ndim
        ind_i = ind_asq + k + (j - 1)*Ndim
        elMat%Aqu(ind_i,ind_j,iel)=elMat%Aqu(ind_i,ind_j,iel)+Aqu(:,:,k)
        exAqq(ind_i,ind_i) = exAqq(ind_i,ind_i) + Aqq
      END DO
    END DO
    CALL invert_matrix(exAqq,iAqq)
    elMat%iAqq(:,:,iel)=elMat%iAqq(:,:,iel)+iAqq
    deallocate(iAqq,exAqq)
  END SUBROUTINE do_assembly

  SUBROUTINE assemblyFacesContribution(iel,NiNi,n_g,ind_asf,ind_ash,ind_fG,ind_ff)
    integer*4,intent(in)   :: iel,ind_asf(:),ind_ash(:),ind_fG(:),ind_ff(:)
    real*8,intent(in)      :: n_g(:)
    real*8,intent(in)      :: NiNi(:,:)
    integer*4              :: indf(size(ind_asf)),ind_1f(size(ind_asf)),ind_2f(size(ind_asf))
    integer*4              :: k,idm
    DO k = 1,Neq
      DO idm = 1,Ndim
        ind_1f = ind_ash + idm + (k - 1)*Ndim
        ind_2f = ind_asf + k
        elMat%Aql(ind_fG(ind_1f),ind_ff(ind_2f),iel) = elMat%Aql(ind_fG(ind_1f),ind_ff(ind_2f),iel) - NiNi*n_g(idm)
      END DO
    END DO
  END SUBROUTINE assemblyFacesContribution





END SUBROUTINE HDG_precalculatedfirstequation

