!*****************************************
! project: MHDG
! file: hdg_ConvectionMatrices.f90
! date: 20/12/2016
! Generate the matrices that change
! during the NR iterative process
!*****************************************

SUBROUTINE HDG_computeJacobian()
  USE globals
  USE LinearAlgebra
  USE analytical
  USE physics
  USE printUtils
  USE MPI_OMP
  USE Debug
  USE hdg_limitingtechniques, ONLY:HDG_ShockCapturing

  IMPLICIT NONE

  !***********************************************************************
  !
  !              COMPUTATION OF THE JACOBIANS
  !
  !***********************************************************************
  integer*4             :: Ndim,Neq,N2D,Npel,Npfl,Nfre,Ng1d,Ng2d,ierr
  integer*4             :: iel,ifa,iface,i,j,els(2),fas(2)
  integer*4             :: sizeu,sizel
  real*8,allocatable    :: ures(:,:),lres(:,:),u0res(:,:,:)
  real*8,allocatable    :: Xel(:,:),Xfl(:,:)
  real*8,allocatable    :: tau_save(:,:),DnnVol(:)
  real*8,allocatable    :: xy_g_save(:,:)
  logical               :: isdir
#ifdef TEMPERATURE
  real*8                :: coefi,coefe
#endif
#ifdef TOR3D
  ! Definitions in 3D
  integer*4             :: itor,iel3,ntorloc,itorg,Np1Dpol,Np1Dtor,Np2D,Npfp,Ng1dpol,Ng1dtor,Ngfp,Ngfl,Nfdir,Ngvo
  integer*4             :: inde(refElTor%Nnodes3D)
  integer*4             :: indft(refElTor%Nfl),indfp(Mesh%Nnodesperelem),indl(Mesh%Nnodesperelem)
  integer*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  integer               :: dd
  integer               :: ind_dim(refElPol%Nfaces + 2),ind_sta(refElPol%Nfaces + 2),aux
  real*8                :: tdiv(numer%ntor + 1),tel(refElTor%Nnodes1d),tg(1),htor
  real*8,allocatable    :: ue(:,:),u0e(:,:,:)
  real*8,allocatable    :: ufp(:,:),uft(:,:)
  real*8,allocatable    :: uefp(:,:),ueft(:,:)
  real*8,allocatable    :: qe(:,:)
  real*8,allocatable    :: qefp(:,:),qeft(:,:)
  real*8,allocatable    :: qres(:,:)
  real*8,allocatable    :: Bel(:,:),fluxel(:),Bfl(:,:),Bfp(:,:)
  integer               :: indbe(refElTor%Nnodes3d),indbp(Mesh%Nnodesperelem),indbt(refElTor%Nfl)
  real*8                :: Jtorel(refElTor%Nnodes3d)
#else
  ! Definitions in 2D
  logical               :: save_tau
  integer*4             :: inde(Mesh%Nnodesperelem)
  integer*4             :: indf(refElPol%Nfacenodes)
  integer*4             :: ind_loc(refElPol%Nfaces,refElPol%Nfacenodes*phys%Neq),perm(refElPol%Nfacenodes*phys%Neq)
  real*8                :: ue(Mesh%Nnodesperelem,phys%Neq),u0e(Mesh%Nnodesperelem,phys%Neq,time%tis)
  real*8                :: uf(refElPol%Nfacenodes,phys%Neq),uef(refElPol%Nfacenodes,phys%Neq)
  integer               :: indtausave(refElPol%Nfaces*refElPol%Ngauss1d)
  integer               :: inddiff_nn_Vol(refElPol%NGauss2D)
  real*8                :: tau_save_el(refElPol%Nfaces*refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Nfaces*refElPol%Ngauss1d,2)
  real*8                :: qe(Mesh%Nnodesperelem,phys%Neq*2),qef(refElPol%Nfacenodes,phys%Neq*2)
  real*8,allocatable    :: qres(:,:)
  real*8                :: Bel(refElPol%Nnodes2d,3),fluxel(refElPol%Nnodes2d),psiel(refElPol%Nnodes2d),Bfl(refElPol%Nfacenodes,3),psifl(refElPol%Nfacenodes)
  real*8                :: Jtorel(refElPol%Nnodes2d)
  real*8                :: n,El_n,nn,El_nn,totaln
  real*8                :: diff_nn_Vol_el(refElPol%NGauss2D),v_nn_Vol_el(refElPol%NGauss2D,Mesh%Ndim),Xg_el(refElPol%NGauss2D,Mesh%Ndim)
  real*8                :: diff_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D),v_nn_Fac_el(refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim)
#endif

  IF (utils%printint > 1) THEN
    WRITE (6,*) '*************************************************'
    WRITE (6,*) '*          COMPUTING JACOBIAN                   *'
    WRITE (6,*) '*************************************************'
  END IF

  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1,timing%clock_rate1)
  end if

  ! Reset matrices
  elMat%Auq = 0.
  elMat%Auu = 0.
  elMat%Aul = 0.
  elMat%Alq = 0.
  elMat%Alu = 0.
  elMat%All = 0.
  elMat%S = 0.
  elMat%fh = 0.

#ifdef TEMPERATURE
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)
#endif

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
  Nfre = refElPol%Nfaces           ! Number of faces in the reference element
  Nfdir = Mesh%Ndir
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

  ! Local indices
  ind_loc = 0
  DO i = 1,Nfre
    DO j = 1,Npfl*Neq
      ind_loc(i,j) = refElPol%Nnodes2D*Neq + Neq*Npfl*(i - 1) + j
    END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol,Np1Dtor,Neq,perm)

#else
  !*************************************************************
  !               2D stuff
  !*************************************************************
  Ndim = 2
  Neq = phys%Neq
  N2D = Mesh%Nelems
  Npel = refElPol%Nnodes2D
  Npfl = refElPol%Nfacenodes
  Nfre = refElPol%Nfaces
  Ng1d = refElPol%Ngauss1d
  Ng2d = refElPol%Ngauss2d

  ind_loc = 0
  DO i = 1,Nfre
    DO j = 1,Neq*Npfl
      ind_loc(i,j) = Neq*Npfl*(i - 1) + j
    END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Npfl,Neq,perm)

  save_tau = switch%saveTau
  if (save_tau) then
    allocate (tau_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,phys%neq))
    allocate (xy_g_save(refElPol%Nfaces*Mesh%Nelems*refElPol%Ngauss1d,2))    
    tau_save = 0.
    xy_g_save = 0.
  endif

  ! Compute shock capturing diffusion
  if (switch%shockcp.gt.0) then
    call HDG_ShockCapturing()
  end if

#endif
  ! reshape
  sizeu = size(sol%u)
  sizel = size(sol%u_tilde)
  allocate (ures(sizeu/Neq,Neq))
  allocate (lres(sizel/Neq,Neq))
  allocate (u0res(sizeu/Neq,Neq,time%tis))
  ures = transpose(reshape(sol%u,[Neq,sizeu/Neq]))
  lres = transpose(reshape(sol%u_tilde,[Neq,sizel/Neq]))
  do i = 1,time%tis
    u0res(:,:,i) = transpose(reshape(sol%u0(:,i),[Neq,sizeu/Neq]))
  end do
  allocate (qres(sizeu/Neq,Neq*Ndim))
  qres = transpose(reshape(sol%q,[Neq*Ndim,sizeu/Neq]))

#ifdef TOR3D
  !********************************************
  !
  !                 3D routines
  !
  !********************************************
  !************************************
  !   Loop in elements
  !************************************
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(itor,iel,itorg,iel3,tel,Xel,indbe,Bel,fluxel,inde,qe,ue,u0e,ifa,indbp,Bfp,dd,indfp,ufp,indl)&
  !$OMP PRIVATE(uefp,qefp,iface,Xfl,indbt,Bfl,isdir,indft,uft,ueft,qeft,i,Jtorel)
  allocate(Xel(Mesh%Nnodesperelem,2))
  allocate(Xfl(refElPol%Nfacenodes,2))
  allocate(Bel(refElTor%Nnodes3d,3),fluxel(refElTor%Nnodes3d),Bfl(refElTor%Nfl,3),Bfp(Mesh%Nnodesperelem,3))
  allocate(ue(refElTor%Nnodes3D,phys%Neq),u0e(refElTor%Nnodes3D,phys%Neq,time%tis))
  allocate(ufp(Mesh%Nnodesperelem,phys%Neq),uft(refElTor%Nfl,phys%Neq))
  allocate(uefp(Mesh%Nnodesperelem,phys%Neq),ueft(refElTor%Nfl,phys%Neq))
  allocate(qe(refElTor%Nnodes3D,phys%Neq*3))
  allocate(qefp(Mesh%Nnodesperelem,phys%Neq*3),qeft(refElTor%Nfl,phys%Neq*3))
  !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
  DO itor = 1,ntorloc
    DO iel = 1,N2D
      ! I made a perfectly nested loop to enable omp parallelization
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))


      !call mpi_barrier(mpi_comm_world,ierr)
      !if (mpivar%glob_id.eq.0) then
      !write(6,*) "TEST 3"
      !endif
      !flush(6)
      !call mpi_barrier(mpi_comm_world,ierr)

      ! Index of 3D element
      iel3 = (itor - 1)*N2d+iel

      ! Coordinates of the nodes of the element
      Xel = Mesh%X(Mesh%T(iel,:),:)

      ! Magnetic field of the nodes of the element
      indbe = colint(tensorSumInt(Mesh%T(iel,:),(itor - 1)*(Np1Dtor - 1)*Mesh%Nnodes + &
        &Mesh%Nnodes*((/(i,i=1,Np1Dtor)/) - 1)))

      Bel = phys%B(indbe,:)
      fluxel = phys%magnetic_flux(indbe)

      ! Ohmic heating (toroidal current)
      IF (switch%ohmicsrc) THEN
        Jtorel = phys%Jtor(indbe)
      ELSE
        Jtorel = 0.
      END IF

      ! Indices to extract the elemental solution
      inde = (iel3 - 1)*Npel + (/(i,i=1,Npel)/)

      qe = qres(inde,:)
      ue = ures(inde,:)
      u0e = u0res(inde,:,:)

      ! Compute the matrices for the element
      CALL elemental_matrices_volume(iel3,Xel,tel,Bel,fluxel,qe,ue,u0e,Jtorel)

      !------------- First poloidal face-----------------
      ifa = 1

      ! Magnetic field of the nodes of the face
      indbp = (itor - 1)*Mesh%Nnodes*(Np1Dtor - 1) + Mesh%T(iel,:)
      Bfp = phys%B(indbp,:)

      ! Face solution
      dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)
      indfp = dd + (iel - 1)*Np2D+(/(i,i=1,Np2d)/)
      ufp = lres(indfp,:)

      ! Elements solution
      indl = (/(i,i=1,Np2d)/)
      uefp = ures(inde(indl),:)
      qefp = qres(inde(indl),:)


      ! Compute the matrices for the element
      CALL elemental_matrices_faces_pol(iel3,ifa,iel,Xel,tdiv(itorg),Bfp,qefp,uefp,ufp)

      !-------------- Toroidal faces ---------------------
      DO ifa=1,refElPol%Nfaces
        iface = Mesh%F(iel,ifa)
        Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

        ! Magnetic field of the nodes of the face
        indbt = colint(tensorSumInt(Mesh%T(iel,refElPol%face_nodes(ifa,:)),(itor-1)*&
          &(Np1Dtor-1)*Mesh%Nnodes+Mesh%Nnodes*((/(i,i=1,Np1Dtor) /)-1)))
        Bfl = phys%B(indbt,:)

        ! Face solution
        isdir = Mesh%Fdir(iel,ifa)
        IF (isdir) THEN
          uft = 0.
        ELSE
          dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl) + N2D*Np2D
          indft = dd + (iface - 1)*Npfl + (/(i,i=1,Npfl)/)
          uft = lres(indft,:)
        ENDIF

        ! Elements solution
        ueft = ures(inde(refElTor%faceNodes3(ifa,:)),:)
        qeft = qres(inde(refElTor%faceNodes3(ifa,:)),:)

        if (iface.le.Mesh%Nintfaces) then
          call elemental_matrices_faces_int(iel3,ifa+1,iel,Xfl,tel,Bfl,qeft,ueft,uft)
        else
          call elemental_matrices_faces_ext(iel3,ifa+1,iel,Xfl,tel,Bfl,qeft,ueft,uft,isdir)
        endif
        ! Flip faces
        IF (Mesh%flipface(iel,ifa)) then
          elMat%Alq(ind_loc(ifa,:),:,iel3) = elMat%Alq(ind_loc(ifa,perm),:,iel3)
          elMat%Alu(ind_loc(ifa,:),:,iel3) = elMat%Alu(ind_loc(ifa,perm),:,iel3)
          elMat%All(ind_loc(ifa,:),:,iel3) = elMat%All(ind_loc(ifa,perm),:,iel3)
          elMat%All(:,ind_loc(ifa,:),iel3) = elMat%All(:,ind_loc(ifa,perm),iel3)
          elMat%fh(ind_loc(ifa,:),iel3) = elMat%fh(ind_loc(ifa,perm),iel3)
        END if
      END DO

      !------------- Second poloidal face-----------------
      ifa = refElPol%Nfaces + 2

      ! Magnetic field of the nodes of the face
      indbp = itor*Mesh%Nnodes*(Np1Dtor - 1) + Mesh%T(iel,:)
      Bfp = phys%B(indbp,:)

      ! Face solution

      if (itor==numer%ntor) then
        dd = 0
      else
        dd = itor*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)
      endif

      indfp = dd + (iel - 1)*Np2D+(/(i,i=1,Np2d)/)
      ufp = lres(indfp,:)

      ! Elements solution
      indl = (/(i,i=Npel - Np2d+1,Npel)/)
      uefp = ures(inde(indl),:)
      qefp = qres(inde(indl),:)

      ! Compute the matrices for the element
      CALL elemental_matrices_faces_pol(iel3,ifa,iel,Xel,tdiv(itorg+1),Bfp,qefp,uefp,ufp)
    END DO
  END DO
  !$OMP END DO
  deallocate(Xel,Xfl,Bel,fluxel,Bfl,Bfp)
  deallocate(ue,u0e,ufp,uft,uefp,ueft,qe,qefp,qeft)
  !$OMP END PARALLEL

  deallocate (ures,lres,u0res)
  deallocate (qres)

  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1,timing%clock_rate1)
    timing%runtjac = timing%runtjac + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputjac = timing%cputjac + timing%tpe1 - timing%tps1
  end if


CONTAINS

  !***************************************************
  ! Volume computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel,tel,Bel,fluxel,qe,ue,u0e,Jtorel)
    integer,intent(IN)          :: iel
    real*8,intent(IN)           :: Xel(:,:),tel(:)
    real*8,intent(IN)           :: Bel(:,:),fluxel(:),Jtorel(:)
    real*8,intent(IN)           :: qe(:,:)
    real*8,intent(IN)           :: ue(:,:)
    real*8,intent(IN)           :: u0e(:,:,:)
    integer*4                   :: g,NGaussPol,NGaussTor,igtor,igpol,i,j,k,iord
    real*8                      :: dvolu,dvolu1d,htor
    real*8                      :: J11(Ng2d),J12(Ng2d)
    real*8                      :: J21(Ng2d),J22(Ng2d)
    real*8                      :: detJ(Ng2d)
    real*8                      :: iJ11(Ng2d),iJ12(Ng2d)
    real*8                      :: iJ21(Ng2d),iJ22(Ng2d)
    real*8                      :: fluxg(Ng2d)
    real*8                      :: xy(Ng2d,2),teg(Ng1dtor)
    real*8                      :: ueg(Ngvo,neq),upg(Ngvo,phys%npv),u0eg(Ngvo,neq,time%tis)
    real*8                      :: qeg(Ngvo,neq*Ndim)
    real*8                      :: force(Ngvo,Neq)
    integer*4,dimension(Npel)   :: ind_ass,ind_asq
    real*8                      :: ktis(time%tis + 1)
    real*8                      :: Nxg(Np2D),Nyg(Np2D),Nx_ax(Np2D)
    real*8,pointer              :: N1g(:),N2g(:),N3g(:)
    real*8                      :: N1xg(Np1Dtor),N1xg_cart(Np1Dtor)
    real*8,dimension(Npel)      :: Nrr,Nr,Nz,Nt,Ni,NNbb
    real*8,dimension(Npel,Npel) :: NNi,NxNi,NyNi,NtNi,NNxy
    real*8                      :: NxyzNi(Npel,Npel,Ndim),Nxyzg(Npel,Ndim)

    real*8                      :: kmult(Npfl,Npfl)
    real*8,parameter            :: tol = 1e-12
    integer*4                   :: ind2(Ng2d)
    real*8                      :: Bmod_nod(Npel),b_nod(Npel,3),b(Ngvo,3),Bmod(Ngvo),divbg,driftg(3),gradbmod(3)
    real*8                      :: bg(3),Jtor(Ngvo)
    real*8                      :: diff_iso_vol(Neq,Neq,Ngvo),diff_ani_vol(Neq,Neq,Ngvo)
    real*8                      :: sigma,sigmax,sigmay,x0,y0,A

    real*8,allocatable  :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)

    !index i from 2nd to 3rd term with a 4th term as step
    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    !    Volume computation
    !***********************************

    ! Gauss points position
    xy = matmul(refElPol%N2D,Xel)

    ! Gauss points position in the toroidal direction
    teg = matmul(refElTor%N1d,tel)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bel(:,1)**2 + Bel(:,2)**2 + Bel(:,3)**2)
    b_nod(:,1) = Bel(:,1)/Bmod_nod
    b_nod(:,2) = Bel(:,2)/Bmod_nod
    b_nod(:,3) = Bel(:,3)/Bmod_nod

    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElTor%N3D,Bmod_nod)
    b = matmul(refElTor%N3D,b_nod)

    ! toroidal current at Gauss points
    IF (switch%ohmicsrc) THEN
      Jtor = matmul(refElPol%N3D,Jtorel)
    ELSE
      Jtor = 0.
    END IF

    ! Solution at Gauss points
    ueg = matmul(refElTor%N3D,ue)
    qeg = matmul(refElTor%N3D,qe)

    ! Compute diffusion at Gauss points
    CALL setLocalDiff(xy,ueg,diff_iso_vol,diff_ani_vol,Bmod)

    ! Solution at previous time steps,at Gauss points
    do i = 1,time%tis
      u0eg(:,:,i) = matmul(refElTor%N3D,u0e(:,:,i))
    end do

    ! Physical variables at Gauss points
    CALL cons2phys(ueg,upg)

    ! Constant sources
    ! Body force at the integration points
    CALL body_force(xy(:,1),xy(:,2),teg,force)

#ifdef NEUTRAL
!    ! Some neutral for WEST
!    IF (switch%testcase .ge. 50 .and. switch%testcase .le. 59) THEN
!      DO g = 1,Ng2d
!        DO igtor =1,Ng1Dtor
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!            ! Case moving equilibrium
!            ! force(g,5) = phys%puff_exp(time%it+1)
!            i = (igtor-1)*Ng2d+g
!            force(i,5) = phys%puff
!          ENDIF
!          !x0 = 2.213
!          !y0 = -0.6968
!          !sigmax = 0.02
!          !sigmay = 0.01
!          !A = phys%lscale**2/(pi*sigmax*sigmay)
!          !force(g,5) = phys%puff*A*exp(-((xy(g,1)*phys%lscale - x0)**2)/(2*sigmax**2) - ((xy(g,2)*phys%lscale - y0)**2)/(2*sigmay**2))
!        ENDDO
!      ENDDO
!    ENDIF
#endif
    ! Some sources for West cases
    IF (switch%testcase .ge. 51 .and. switch%testcase .le. 55) THEN
      fluxg = matmul(refElPol%N2D,fluxel)
      DO g = 1,Ngvo
        IF (switch%testcase == 51) THEN
          IF (fluxg(g) .le. -0.88 .and. fluxg(g) .ge. -0.90) THEN
            !force(g,1) = 3.20119388718018e-05
            force(g,1) = 4.782676673609557e-05
          END IF
        ELSE IF (switch%testcase == 52) THEN
          sigma = 0.1
          x0 = 0.
          A = (phys%lscale**2)/(2*PI*sigma**2)
#ifndef NEUTRAL
          force(g,1) = 0.4*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
#ifdef TEMPERATURE
          force(g,3) = 0.
          force(g,4) = 30.*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
          !IF (fluxg(g) .le. -0.90 .and. fluxg(g) .ge. -1.) THEN
          !  force(g,1) = 9.45155008295538e-06
          !END IF
        ELSE IF (switch%testcase == 53) THEN
          IF (fluxg(g) .le. -0.90) THEN
            force(g,1) = 7.24032211339971e-06
          END IF
#ifdef TEMPERATURE
          force(g,3) = 18*force(g,1)
          force(g,4) = force(g,3)
#endif
        ELSE IF (switch%testcase == 54) THEN
#ifndef NEUTRAL
          sigma = 0.22
          x0 = 0.
          A = (phys%lscale**2)/(2*PI*sigma**2)
          IF (fluxg(g) .le. 0.35) THEN
            force(g,1) = 1.*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
          ENDIF
          !IF (fluxg(g) .le. -1.03) THEN
          !  force(g,1) = 0.000115575293741846
          !END IF
#endif
        ELSE IF (switch%testcase == 55) THEN
          IF (fluxg(g) .le. -0.88 .and. fluxg(g) .ge. -0.90) THEN
            force(g,1) = 10
          END IF
#ifdef TEMPERATURE
          force(g,3) = 18*force(g,1)
          force(g,4) = force(g,3)
#endif
        END IF
      END DO
    END IF
    !! end sources

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

    ! Coefficient time integration scheme
    call setTimeIntegrationCoefficients(ktis)

    NgaussPol = refElPol%NGauss2D
    NgaussTor = refElTor%NGauss1D


    ! Allocate temporary matrices
    allocate(Auq(Npel,Npel, neq*neq*ndim  ))
    allocate(Auu(Npel,Npel, neq*neq  ))
    allocate(rhs(Npel,Neq))
    Auq = 0.
    Auu = 0.
    rhs = 0.
    DO igtor = 1,NGaussTor
      N1g => refElTor%N1D(igtor,:)         ! Toroidal shape function
      N1xg_cart = refElTor%Nxi1D(igtor,:)*2/htor       ! Toroidal shape function derivative
      dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume

      DO igpol = 1,NGaussPol
        g = (igtor - 1)*NGaussPol + igpol

        ! Poloidal shape functions and derivatives
        N2g => refElPol%N2D(igpol,:)
        Nxg = iJ11(igpol)*refElPol%Nxi2D(igpol,:) + iJ12(igpol)*refElPol%Neta2D(igpol,:)
        Nyg = iJ21(igpol)*refElPol%Nxi2D(igpol,:) + iJ22(igpol)*refElPol%Neta2D(igpol,:)

        ! 3D integration weight
        dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d

        IF (switch%axisym) THEN
          dvolu = dvolu*xy(igpol,1)
          N1xg = N1xg_cart/xy(igpol,1)
          Nx_ax = Nxg + 1./xy(igpol,1)*N2g
        ELSE
          N1xg = N1xg_cart
          Nx_ax = Nxg
        END IF

        ! 3D shape functions
        N3g => refElTor%N3D(g,:)               ! 3D shape function
        Nrr = col(TensorProduct(Nx_ax,N1g))    ! 3D shape function,derivative in r for computing the divergence
        Nr = col(TensorProduct(Nxg,N1g))      ! 3D shape function,derivative in r
        Nz = col(TensorProduct(Nyg,N1g))      ! 3D shape function,derivative in z
        Nt = col(TensorProduct(N2g,N1xg))     ! 3D shape function,derivative in t

        ! Shape functions products
        Ni = N3g*dvolu                                                ! Npel x 1
        NNi = tensorProduct(N3g,Ni)                                    ! Npel x Npel
        NxNi = tensorProduct(Nr,Ni)                                     ! Npel x Npel
        NyNi = tensorProduct(Nz,Ni)                                     ! Npel x Npel
        NtNi = tensorProduct(Nt,Ni)                                     ! Npel x Npel
        NNxy = b(g,1)*NxNi + b(g,2)*NyNi + b(g,3)*NtNi                      ! Npel x Npel
        NxyzNi(:,:,1) = NxNi
        NxyzNi(:,:,2) = NyNi
        NxyzNi(:,:,3) = NtNi                                           ! Npel x Npel x 3
        NNbb = (Nr*b(g,1) + Nz*b(g,2) + Nt*b(g,3))*dvolu                   ! Npel x 1
        Nxyzg(:,1) = Nr*dvolu
        Nxyzg(:,2) = Nz*dvolu
        Nxyzg(:,3) = Nt*dvolu

        divbg = dot_product(Nrr,b_nod(:,1)) + dot_product(Nz,b_nod(:,2)) + dot_product(Nt,b_nod(:,3))

        ! Diamagnetic drift !TODO: verify drift intensity in isothermal and non-isothermal cases
        driftg = 0.
        gradbmod = 0.
        gradbmod(1) = dot_product(Nr,Bmod_nod)
        gradbmod(2) = dot_product(Nz,Bmod_nod)
        gradbmod(3) = dot_product(Nt,Bmod_nod)
        bg = b(g,:)
        call cross_product(bg,gradbmod,driftg)
        driftg = phys%dfcoef*driftg/Bmod(g)

        CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),divbg,driftg,Bmod(g),force(g,:),&
          &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
          &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g))
      END DO ! END loop in volume Gauss points
    END DO
    call do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
    deallocate(Auq,Auu,rhs)

  END SUBROUTINE elemental_matrices_volume

  !*****************************************
  ! Poloidal faces computations in 3D
  !*****************************************
  SUBROUTINE elemental_matrices_faces_pol(iel,ifa,iel2,Xfp,tg,Bfp,qef,uef,uf)
    integer*4,intent(IN)  :: iel,ifa,iel2
    real*8,intent(IN)     :: Xfp(:,:),tg(1)
    real*8,intent(IN)     :: Bfp(:,:)
    real*8,intent(IN)     :: qef(:,:)
    real*8,intent(IN)     :: uef(:,:),uf(:,:)
    real*8                :: ufg(Ng2d,Neq),uefg(Ng2d,Neq),upgf(Ng2d,phys%npv)
    real*8                :: qfg(Ng2d,Neq*Ndim)
    integer*4             :: ind_asf(Np2D),ind_ash(Np2D)
    integer*4             :: g,NGauss,i,j,lel
    real*8                :: dsurf(Ng2d),bn
    real*8                :: xyf(Ng2d,2)
    real*8                :: J11(Ng2d),J12(Ng2d)
    real*8                :: J21(Ng2d),J22(Ng2d)
    real*8                :: detJ(Ng2d)
    real*8                :: iJ11(Ng2d),iJ12(Ng2d)
    real*8                :: iJ21(Ng2d),iJ22(Ng2d)
    real*8,pointer        :: Nfg(:)
    real*8                :: NNif(Np2D,Np2D),Nif(Np2D),Nfbn(Np2D)
    real*8                :: n_g(Ng2d,3)
    integer*4             :: ind_ff(Neq*Np2D),ind_fe(Neq*Np2D),ind_fg(Neq*Ndim*Np2D)
    real*8,parameter     :: tol = 1e-12
    real*8                :: Bmod_nod(Np2D),b_nod(Np2D,3),b(Ng2d,3),Bmod(Ng2d)
    real*8                :: tau(Neq,Neq)
    real*8                :: diff_iso_fac(Neq,Neq,Ng2D),diff_ani_fac(Neq,Neq,Ng2D)

    ind_asf = (/(i,i=0,Neq*(Np2D-1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Np2D-1)*Ndim,Neq*Ndim)/)

    ! Gauss points position
    xyf = matmul(refElPol%N2D,Xfp)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bfp(:,1)**2 + Bfp(:,2)**2 + Bfp(:,3)**2)
    b_nod(:,1) = Bfp(:,1)/Bmod_nod
    b_nod(:,2) = Bfp(:,2)/Bmod_nod
    b_nod(:,3) = Bfp(:,3)/Bmod_nod
    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElPol%N2D,Bmod_nod)
    b = matmul(refElPol%N2D,b_nod)

    ! Trace solution at face Gauss points
    ufg = matmul(refElPol%N2D,uf)

    ! Compute diffusion at Gauss points
    CALL setLocalDiff(xyf,ufg,diff_iso_fac,diff_ani_fac,Bmod)

    ! Loop in 2D Gauss points
    Ngauss = Ng2d

    ! Physical variables related to the trace solution
    CALL cons2phys(ufg,upgf)

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
      ind_fe = (/(i,i=1,Np2D*Neq)/)
      ind_fg = (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = -1
    ELSE
      ind_ff = Np2D*Neq + refElPol%Nfaces*Npfl*Neq + (/(i,i=1,Np2D*Neq)/)
      ind_fe = Np2d*(Np1dTor - 1)*Neq + (/(i,i=1,Np2D*Neq)/)
      ind_fg = Np2d*(Np1dTor - 1)*Ndim*Neq + (/(i,i=1,Np2D*Ndim*Neq)/)
      ! Exterior normal
      n_g = 0.; n_g(:,3) = 1
    ENDIF

    ! Element solution at face Gauss points
    uefg = matmul(refElPol%N2D,uef)
    ! Gradient solution at face gauss points
    qfg = matmul(refElPol%N2D,qef)

    DO g = 1,NGauss

      ! Shape functions product
      Nfg => refElPol%N2D(g,:)
      bn = dot_product(b(g,:),n_g(g,:))
      NNif = tensorProduct(Nfg,Nfg)*dsurf(g)
      Nif = Nfg*dsurf(g)
      Nfbn = bn*Nfg*dsurf(g)

      ! Compute the stabilization term
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,ifa,0.,xyf(g,:),tau)
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
        ENDIF
      END IF

      ! Assembly local contributions
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
        n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
    END DO
    !      END DO

  END SUBROUTINE elemental_matrices_faces_pol

  !***************************************************
  ! Interior faces computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_int(iel,ifa,iel2,Xfl,tel,Bfl,qef,uef,uf)
    integer,intent(IN)        :: iel,ifa,iel2
    real*8,intent(IN)         :: Xfl(:,:),tel(:)
    real*8,intent(IN)         :: Bfl(:,:)
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(IN)         :: uef(:,:),uf(:,:)

    integer*4                 :: g,i,j,k,igtor,igpol
    real*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
    real*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    real*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
    real*8                    :: dsurf(Ngfl),dsurfg
    real*8                    :: qfg(Ngfl,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl) :: ind_asf,ind_ash
    integer*4,dimension(Npfl) :: indf,ind_if,ind_jf,ind_kf
    integer*4                 :: ind(Ng1Dpol)
    integer*4                 :: permsing(Npfl)
    real*8                    :: t_g(Ng1dpol,2),n_g(Ngfl,3),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8,pointer            :: Nfg(:)
    real*8                    :: tau(Neq,Neq)
    real*8,parameter         :: tol = 1e-12
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ngfl,3),Bmod(Ngfl)
    real*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Gauss points position in the toroidal direction
    teg = matmul(refElTor%N1d,tel)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Loop in the faces of the element
    !      DO lel = 1,2

    !         iel = els(lel)
    !         ifa = fas(lel)
    !         if (lel == 1) then
    !            ieln = els(2)
    !         else
    !            ieln = els(1)
    !         endif

    ! Indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*Neq))
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Coordinates,derivatives and trace solution at face Gauss points
    IF (Mesh%flipFace(iel2,ifa - 1)) THEN
      call set_permutations(Np1Dpol,Np1Dtor,1,permsing)
      ufg = matmul(refElTor%sFTF,uf(permsing,:))
    ELSE
      ufg = matmul(refElTor%sFTF,uf)
    END IF

    xyf = matmul(refElPol%N1D,Xfl)
    ! Shape function derivatives at Gauss points
    xyDer = matmul(refElPol%Nxi1D,Xfl)
    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElTor%sFTF,Bmod_nod)
    b = matmul(refElTor%sFTF,b_nod)

    ! Element solution at face Gauss points
    uefg = matmul(refElTor%sFTF,uef)
    ! Gradient solution at face gauss points
    qfg = matmul(refElTor%sFTF,qef)
    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac,Bmod)

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

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
        bn = dot_product(b(g,:),n_g(g,:))
        NNif = tensorProduct(Nfg,Nfg)*dsurfg
        Nif = Nfg*dsurfg
        Nfbn = bn*Nfg*dsurfg

        ! Compute the stabilization term
        tau = 0.
        IF (numer%stab == 1) THEN
          ! Constant stabilization
          DO i = 1,Neq
            tau(i,i) = numer%tau(i)
          END DO
        ELSE
          ! Non constant stabilization
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,ifa,0.,xyf(g,:),tau)
          ELSE
            CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
          ENDIF
        END IF

        ! Assembly local contributions
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)

      END DO ! Gauss points
    END DO
    !      END DO ! 2 elements

  END SUBROUTINE elemental_matrices_faces_int

  !***************************************************
  ! Exterior faces computation in 3D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,iel2,Xfl,tel,Bfl,qef,uef,uf,isdir)
    integer,intent(IN)                             :: iel,ifa,iel2
    real*8,intent(IN)         :: Xfl(:,:),tel(:)
    real*8,intent(IN)         :: Bfl(:,:)
    logical,intent(IN)        :: isdir
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(INOUT)      :: uef(:,:),uf(:,:)

    integer*4                 :: i,j,k,g,igtor,igpol
    real*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
    real*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)
    real*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
    real*8                    :: dsurf(Ngfl),dsurfg
    real*8                    :: qfg(Ngfl,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl) :: ind_asf,ind_ash
    integer*4,dimension(Npfl) :: indf,ind_if,ind_jf,ind_kf
    integer*4                 :: ind(Ng1Dpol)
    integer*4                 :: permsing(Npfl)
    real                      :: isext
    real*8                    :: t_g(Ng1Dpol,2),n_g(Ngfl,3),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8,pointer            :: Nfg(:)
    real*8                    :: tau(Neq,Neq)

    real*8,parameter         :: tol = 1e-12
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ngfl,3),Bmod(Ngfl)
    real*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    ! toroidal element size
    htor = tel(Np1Dtor) - tel(1)

    ! Gauss points position in the toroidal direction
    teg = matmul(refElTor%N1d,tel)

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElTor%sFTF,Bmod_nod)
    b = matmul(refElTor%sFTF,b_nod)

    ! Indices
    ind_ff = Np2d*Neq + (ifa - 2)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
    ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*Neq))
    ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa - 1,:) - 1)*3*Neq))

    ! Trace solution at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    IF (isdir) THEN
      CALL analytical_solution(iel,xyf(:,1),xyf(:,2),teg,ufg)
    ELSE
#ifdef PARALL
      IF (Mesh%flipface(iel2,ifa - 1)) THEN
        call set_permutations(Np1Dpol,Np1Dtor,1,permsing)
        uf = uf(permsing,:)
      END IF
      ! TODO: VERIFY IF I NEED TO FLIP ALSO xyf,b and Bmod in this case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      ufg = matmul(refElTor%sFTF,uf)
    END IF

    ! Element solution at face Gauss points
    uefg = matmul(refElTor%sFTF,uef)
    ! Gradient solution at face gauss points
    qfg = matmul(refElTor%sFTF,qef)

    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,diff_iso_fac,diff_ani_fac,Bmod)

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

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
        bn = dot_product(b(g,:),n_g(g,:))
        NNif = tensorProduct(Nfg,Nfg)*dsurfg
        Nif = Nfg*dsurfg
        Nfbn = bn*Nfg*dsurfg

        ! Compute the stabilization term
        isext = 1.
#ifdef PARALL
        IF (Mesh%boundaryFlag(Mesh%F(iel2,ifa - 1) - Mesh%Nintfaces) .eq. 0) THEN
          isext = 0.
        END IF
#endif
        tau = 0.
        IF (numer%stab == 1) THEN
          ! Constant stabilization
          DO i = 1,Neq
            tau(i,i) = numer%tau(i)
          END DO
        ELSE
          ! Non constant stabilization
          ! Compute tau in the Gauss points
          IF (numer%stab < 6) THEN
            CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b(g,:),Bmod(g),n_g(g,:),iel2,ifa,isext,xyf(g,:),tau)
          ELSE
            CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),isext,iel2,tau)
          ENDIF
        END IF

        ! Assembly local contributions
#ifdef PARALL
        IF (Mesh%boundaryFlag(Mesh%F(iel2,ifa - 1) - Mesh%Nintfaces) .eq. 0) THEN
          ! Ghost face: assembly it as interior
          CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
        ELSE
          CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
            n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)
        ENDIF
#else
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),&
          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,ifa)

#endif
      END DO ! Gauss points
    END DO

  END SUBROUTINE elemental_matrices_faces_ext

  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(Np1Dpol,Np1Dtor,Neq,perm)
    integer,intent(IN)  :: Np1Dpol,Np1Dtor,Neq
    integer,intent(OUT) :: perm(:)
    integer              :: i,j,k
    integer              :: temp_pol(Np1Dpol),temp_tor(Np1Dtor),aux(Np1Dpol,Np1Dtor)

    temp_pol = (/(i,i=1,Np1Dpol)/)
    temp_tor = Np1Dpol*((/(i,i=1,Np1Dtor)/) - 1)
    aux = TensorSumInt(temp_pol,temp_tor)
    DO j = 1,Np1Dtor
      DO i = 1,Np1Dpol
        DO k = 1,Neq
          perm((j - 1)*Np1Dpol*Neq + (i - 1)*Neq + k) = (aux(Np1Dpol - i + 1,j) - 1)*Neq + k
        END DO
      END DO
    END DO
  END SUBROUTINE set_permutations

#else
!TOR3D

  !********************************************
  !
  !                 2D routines
  !
  !********************************************

  !************************************
  !   Loop in elements in 2D
  !************************************
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(iel,ifa,iface,inde,indf,Xel,Xfl,i,qe,qef,ue,uef,uf,u0e,Bel,Bfl,fluxel,psiel,psifl,isdir,Jtorel,El_n,El_nn) &
  !$OMP PRIVATE(Xg_el,diff_nn_Vol_el,diff_nn_Fac_el,v_nn_Vol_el,v_nn_Fac_el,xy_g_save,xy_g_save_el,tau_save,tau_save_el)  
  allocate(Xel(Mesh%Nnodesperelem,2))
  allocate(Xfl(refElPol%Nfacenodes,2))

  n = 0.
  nn = 0.
  !$OMP DO SCHEDULE(STATIC) REDUCTION(+:n,nn)
  DO iel = 1,N2D

    ! Coordinates of the nodes of the element
    Xel = Mesh%X(Mesh%T(iel,:),:)

    ! Magnetic field of the nodes of the element
    Bel = phys%B(Mesh%T(iel,:),:)
    fluxel = phys%magnetic_flux(Mesh%T(iel,:))
    
    ! Normalized magnetic flux of the nodes of the element: PSI el
    psiel = phys%magnetic_psi(Mesh%T(iel,:))
    
    ! Ohmic heating (toroidal current)
    IF (switch%ohmicsrc) THEN
      Jtorel = phys%Jtor(Mesh%T(iel,:))
    ELSE
      Jtorel = 0.
    END IF

    ! Indices to extract the elemental and face solution
    inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)

    qe = qres(inde,:)
    ue = ures(inde,:)
    u0e = u0res(inde,:,:)

    ! Compute the matrices for the element
    CALL elemental_matrices_volume(iel,Xel,Bel,fluxel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
    if (save_tau) then
       inddiff_nn_Vol = (iel - 1)*refElPol%NGauss2D+(/(i,i=1,refElPol%NGauss2D)/)
       phys%diff_nn_Vol(inddiff_nn_Vol) = diff_nn_Vol_el
       phys%v_nn_Vol(inddiff_nn_Vol,:) = v_nn_Vol_el
       Mesh%Xg(inddiff_nn_Vol,:) = Xg_el
    endif

    ! Compute total plasma and neutral density (don't add contribution of ghost elements)
#ifdef PARALL
    IF (Mesh%ghostElems(iel) .eq. 0) THEN
#endif
      n  = n + El_n
      nn = nn + El_nn
#ifdef PARALL
    ENDIF
#endif
    ! Loop in local faces
    if (save_tau) then
       diff_nn_Fac_el = 0.
       v_nn_Fac_el = 0.
       tau_save_el = 0.
       xy_g_save_el = 0;
    endif
    
    DO ifa=1,refElPol%Nfaces
      iface = Mesh%F(iel,ifa)
      isdir = Mesh%Fdir(iel,ifa)

      ! Coordinates of the nodes of the face
      Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

      ! Magnetic field of the nodes of the face
      Bfl = phys%B(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)
   
      ! Normalized magnetic flux of the nodes of the face: PSI fl
      psifl = phys%magnetic_psi(Mesh%T(iel,refElPol%face_nodes(ifa,:)))
            
      ! Face solution
      indf = (iface-1)*Npfl + (/(i,i=1,Npfl)/)
      uf = lres(indf,:)

      ! Elements solution
      inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)
      uef = ures(inde(refElPol%face_nodes(ifa,:)),:)
      qef = qres(inde(refElPol%face_nodes(ifa,:)),:)

      if (iface.le.Mesh%Nintfaces) then
        CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)		 
      else
        if (Mesh%periodic_faces(iface-Mesh%Nintfaces).eq.0) then
          CALL elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        else
          ! periodic face
          CALL elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
        endif
      endif
	 
      ! Flip faces
      if (Mesh%flipface(iel,ifa)) then
        elMat%Alq(ind_loc(ifa,:),:,iel) = elMat%Alq(ind_loc(ifa,perm),:,iel)
        elMat%Alu(ind_loc(ifa,:),:,iel) = elMat%Alu(ind_loc(ifa,perm),:,iel)
        elMat%All(ind_loc(ifa,:),:,iel) = elMat%All(ind_loc(ifa,perm),:,iel)
        elMat%All(:,ind_loc(ifa,:),iel) = elMat%All(:,ind_loc(ifa,perm),iel)
        elMat%fh(ind_loc(ifa,:),iel) = elMat%fh(ind_loc(ifa,perm),iel)
      end if
    END DO
	
    if (save_tau) then
       indtausave = (iel - 1)*refElPol%Nfaces*refElPol%Ngauss1d+(/(i,i=1,refElPol%Nfaces*refElPol%Ngauss1d)/)
	   phys%diff_nn_Fac(indtausave) = diff_nn_Fac_el
	   phys%v_nn_Fac(indtausave,:) = v_nn_Fac_el
	   tau_save(indtausave,:) = tau_save_el
	   Mesh%Xgf(indtausave,:) = xy_g_save_el
       xy_g_save(indtausave,:) = xy_g_save_el
    end if


  END DO
  !$OMP END DO
  deallocate(Xel,Xfl)
  !$OMP END PARALLEL

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, nn, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

   IF (MPIvar%glob_id.eq.0) THEN
     if(switch%ME .eq. .TRUE. .AND. switch%testcase .ge. 80) then
        WRITE(6,*) 'D_n = ', phys%ME_diff_n*simpar%refval_length**2/simpar%refval_time
     endif 
     totaln = n + nn
     WRITE(6,*) 'n = ',n
     WRITE(6,*) 'nn = ',nn
     WRITE(6,*) 'total n = ',totaln
   ENDIF

  deallocate (ures,lres,u0res)
  deallocate (qres)

  if (save_tau) then
    write (6,*) "Saving Dnn in 2D gauss points"
    write (6,*) "Saving tau in the faces"
    call saveMatrix(tau_save,'tau_save')
    call saveMatrix(xy_g_save,'xy_g_save')
    deallocate (tau_save,xy_g_save)
    write (6,*) "Done saving tau!"
  endif

  IF (utils%printint > 1) THEN
    WRITE (6,*) "Done!"
  END IF


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1,timing%clock_rate1)
    timing%runtjac = timing%runtjac + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputjac = timing%cputjac + timing%tpe1 - timing%tps1
  end if


CONTAINS

  !***************************************************
  ! Volume computation in 2D
  !***************************************************
  SUBROUTINE elemental_matrices_volume(iel,Xel,Bel,fluxel,psiel,qe,ue,u0e,Jtorel,El_n,El_nn,diff_nn_Vol_el,v_nn_Vol_el,Xg_el)
    integer,intent(IN)            :: iel
    real*8,intent(IN)             :: Xel(:,:)
    real*8,intent(IN)             :: Bel(:,:),fluxel(:),psiel(:),Jtorel(:)
    real*8,intent(IN)             :: qe(:,:)
    real*8,intent(IN)             :: ue(:,:),u0e(:,:,:)
    real*8,intent(OUT)            :: El_n,El_nn
    real*8,intent(OUT)            :: diff_nn_Vol_el(Ng2D),v_nn_Vol_el(Ng2D,ndim),Xg_el(Ng2D,ndim)
    integer*4                     :: g,NGauss,i,j,k
    real*8                        :: dvolu
    real*8                        :: xy(Ng2d,ndim),ueg(Ng2d,neq),u0eg(Ng2d,neq,time%tis)
    real*8                        :: force(Ng2d,Neq)
    real*8                        :: qeg(Ng2d,neq*Ndim)
    real*8                        :: J11(Ng2d),J12(Ng2d)
    real*8                        :: J21(Ng2d),J22(Ng2d)
    real*8                        :: detJ(Ng2d)
    real*8                        :: iJ11(Ng2d),iJ12(Ng2d)
    real*8                        :: iJ21(Ng2d),iJ22(Ng2d)
    real*8                        :: fluxg(Ng2d),max_flux2D,min_flux2D, Psig(Ng2d)
    integer*4,dimension(Npel)     :: ind_ass,ind_asq
    real*8                        :: ktis(time%tis + 1)
    real*8,dimension(Npel)        :: Ni,Nxg,Nyg,NNbb,Nx_ax
    real*8,dimension(Npel,Npel)   :: NxNi,NyNi,NNxy,NNi
    real*8                        :: NxyzNi(Npel,Npel,3),Nxyzg(Npel,3)
    real*8                        :: upg(Ng2d,phys%npv)
    real*8                        :: Bmod_nod(Npel),b_nod(Npel,3),b(Ng2d,3),Bmod(Ng2d),divbg,driftg(3),gradbmod(3)
    real*8                        :: bg(3), Jtor(Ng2d)
    real*8                        :: diff_iso_vol(Neq,Neq,Ng2d),diff_ani_vol(Neq,Neq,Ng2d)
    real*8,allocatable            :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
    real*8                        :: auxdiffsc(Ng2d)
    real*8                        :: Pi,sigma,sigmax,sigmay,x0,y0,A,r
    real*8                        :: th_n = 1.e-14
    real*8                        :: Vnng(Ndim)
    
    if (save_tau) then
       Xg_el = 0.
       diff_nn_Vol_el = 0.
       v_nn_Vol_el = 0.
    endif

    ind_ass = (/(i,i=0,Neq*(Npel - 1),Neq)/)
    ind_asq = (/(i,i=0,Neq*(Npel - 1)*Ndim,Neq*Ndim)/)

    force = 0.
    El_n  = 0.
    El_nn  = 0.
    Pi = 3.1415926535
    !***********************************
    !    Volume computation
    !***********************************

    ! Gauss points position
    xy = matmul(refElPol%N2D,Xel)
    if (save_tau) then
       Xg_el = xy;
    endif

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bel(:,1)**2 + Bel(:,2)**2 + Bel(:,3)**2)
    b_nod(:,1) = Bel(:,1)/Bmod_nod
    b_nod(:,2) = Bel(:,2)/Bmod_nod
    b_nod(:,3) = Bel(:,3)/Bmod_nod

    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElPol%N2D,Bmod_nod)
    b = matmul(refElPol%N2D,b_nod)
   
    ! Normalized magnetic flux at Gauss points: PSI
    Psig = matmul(refElPol%N2D,psiel) 
  
    ! toroidal current at Gauss points
    IF (switch%ohmicsrc) THEN
      Jtor = matmul(refElPol%N2D,Jtorel)
    ELSE
      Jtor = 0.
    END IF

    ! Solution at Gauss points
    ueg = matmul(refElPol%N2D,ue)
    qeg = matmul(refElPol%N2D,qe)

    ! Compute diffusion at Gauss points
    CALL setLocalDiff(xy,ueg,qeg,diff_iso_vol,diff_ani_vol)
    if (save_tau) then
       diff_nn_Vol_el = diff_iso_vol(5,5,:)
    endif
    
    if (switch%shockcp.gt.0) then
      auxdiffsc = matmul(refElPol%N2D,Mesh%scdiff_nodes(iel,:))
      do i=1,Neq
        diff_iso_vol(i,i,:) = diff_iso_vol(i,i,:)+auxdiffsc
      end do
    endif

    ! Solution at previous time steps,at Gauss points
    do i = 1,time%tis
      u0eg(:,:,i) = matmul(refElPol%N2D,u0e(:,:,i))
    end do

    ! Physical variables at Gauss points
    CALL cons2phys(ueg,upg)

    ! Constant sources
    ! Body force at the integration points
    CALL body_force(xy(:,1),xy(:,2),force)
    
    ! Some sources to limit low density and temeprauture values
    !DO g=1, Ng2d
    !   IF (ueg(g,1) .lt. 1.e-7) force(g,1) = 1.e-7 - ueg(g,1) !Th at n = 1.00E+12 [m^(-3)]
    !   IF (upg(g,7) .lt. 6.e-4) force(g,3) = 3./2.*ueg(g,1)*(1.e-3 - upg(g,7)) !Th at Ti 0.03 eV 
    !   IF (upg(g,8) .lt. 6.e-4) force(g,4) = 3./2.*ueg(g,1)*(1.e-3 - upg(g,8)) !Th at Te 0.03 eV
    !   IF (ueg(g,5) .gt. 1.e+0) force(g,5) = 1.e+0 - ueg(g,5) !Th at nEe
    !    IF (ueg(g,3) .lt. 2.e-4) force(g,3) = 2.e-5
    !    IF (ueg(g,4) .lt. 2.e-4) force(g,4) = 2.e-5
    !END DO
   
    ! Some sources for West cases
    IF (switch%testcase .ge. 50 .and. switch%testcase .le. 59) THEN
      ! Compute flux surfaces and normalise them
      fluxg = matmul(refElPol%N2D,fluxel)
      max_flux2D = maxval(phys%magnetic_flux)
      min_flux2D = minval(phys%magnetic_flux)
      fluxg = (fluxg - min_flux2D)/(max_flux2D - min_flux2D)
      DO g = 1,Ng2d
        IF (ueg(g,1) .lt. th_n) THEN
          force(g,1) = th_n - ueg(g,1)
        ENDIF
        ! WEST CASE with analytical Gaussian sources on density and energies, no puff.
        IF (switch%testcase == 52) THEN
          sigma = phys%sigma_source
          x0 = 0.
          A = (phys%lscale**2)/sqrt((2*Pi*sigma**2))
          IF (fluxg(g) .le. phys%fluxg_trunc) THEN
            force(g,1) = phys%density_source*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
#ifdef TEMPERATURE
            force(g,3) = phys%ener_source_e*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
            force(g,4) = phys%ener_source_ee*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
          ENDIF
!#ifdef NEUTRAL
!        ! WEST CASE only with analytical puff. Energy source is given by JTOR.
!        ELSE IF (switch%testcase == 54) THEN
!          ! location of the buffer differs for the mesh 3895_P4 and for the mesh 26830_P4
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!            force(g,5) = phys%puff
!          ENDIF
!        ELSE IF(switch%testcase .eq. 59) THEN
!          ! location of the buffer differs for the mesh 3895_P4 and for the mesh 26830_P4
!          !IF (xy(g,1)*phys%lscale .gt. 2.36 .and. xy(g,2)*phys%lscale .lt. -0.69 ) THEN
!          IF (xy(g,1)*phys%lscale .gt. 2.446 .and. xy(g,1)*phys%lscale .lt. 2.59 .and. xy(g,2)*phys%lscale .gt. -0.7964 .and. xy(g,2)*phys%lscale .lt. -0.7304 ) THEN
!           ! If it is a time initialization simulation, the puff is analytical (from param.txt) otherwise experimental
!            IF(switch%time_init) THEN
!              force(g,5) = phys%puff
!            ELSE
!              force(g,5) = phys%puff_exp(time%it+1)
!            ENDIF
!          ENDIF
!#endif
        ENDIF
      END DO
    END IF
    
    ! Some sources for ITER cases    
    IF (switch%testcase .ge. 80) THEN
       ! Compute flux surfaces and normalise them
       fluxg = matmul(refElPol%N2D,fluxel)
       max_flux2D = maxval(phys%magnetic_flux)
       min_flux2D = minval(phys%magnetic_flux)
       fluxg = (fluxg - min_flux2D)/(max_flux2D - min_flux2D)
       IF (switch%testcase == 81) THEN
          sigma = phys%sigma_source
          x0 = 0.
          A = (phys%lscale**2)/sqrt((2*Pi*sigma**2))
          ! Only energy sources: density from neutral model
          IF (fluxg(g) .le. phys%fluxg_trunc) THEN
#ifdef NEUTRAL   
#ifdef TEMPERATURE
          force(g,3) = phys%ener_source_e*A*exp(-((fluxg(g)-x0)**2)/(2*sigma**2))
          force(g,4) = phys%ener_source_ee*A*exp(-((fluxg(g)-x0)**2)/(2*sigma**2))
#endif
#endif
          ENDIF
       ELSE IF (switch%testcase == 82) THEN
          sigma = phys%sigma_source
          A = (phys%lscale**2)/sqrt((2*Pi*sigma**2))
          ! Only energy sources: density from neutral model
#ifdef NEUTRAL
#ifdef TEMPERATURE
          force(g,3) = phys%ener_source_e*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
          force(g,4) = phys%ener_source_ee*A*exp(-((fluxg(g) - x0)**2)/(2*sigma**2))
#endif
#endif
       ENDIF
    ENDIF 
    ! end sources

    ! Some sources for Circular cases
    IF (switch%testcase .ge. 60) THEN
      DO g=1,Ng2D
			IF (switch%testcase==61) THEN
			     r =   sqrt ( (xy(g,1)*phys%lscale-geom%R0)**2 +(xy(g,2)*phys%lscale)**2 )
					 IF (r .le. 0.4) THEN
					        force(g,1) = 0. !2.838272668283863e-05
#ifdef TEMPERATURE
                  force(g,3) = 8*2.838272668283863e-05 !force(g,1)
                  force(g,4) = 8*2.838272668283863e-05
#endif
					 END IF
      ELSE IF (switch%testcase==62) THEN
#ifdef TEMPERATURE
            Pi = 3.1415926535
            sigma = 0.3
            x0 = geom%R0
            A = (phys%lscale**2)/(2*Pi*sigma**2)
            force(g,1) = 0.
            force(g,3) = 18.*A*exp(-((xy(g,1)*phys%lscale - x0)**2 + (xy(g,2)*phys%lscale)**2)/(2*sigma**2))
            force(g,4) = force(g,3)
#endif
			END IF
			END DO
		END IF
		!! end sources


    ! Loop in 2D Gauss points
    Ngauss = Ng2d
    J11 = matmul(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
    J12 = matmul(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
    J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
    J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
    detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
    iJ11 = J22/detJ
    iJ12 = -J12/detJ
    iJ21 = -J21/detJ
    iJ22 = J11/detJ

    ! Coefficient time integration scheme
    call setTimeIntegrationCoefficients(ktis)
    ! Allocate temporary matrices
    allocate(Auq(Npel,Npel, neq*neq*ndim  ))
    allocate(Auu(Npel,Npel, neq*neq  ))
    allocate(rhs(Npel,Neq))
    Auq = 0.
    Auu = 0.
    rhs = 0.
    ! Loop in 2D Gauss points
    DO g = 1,NGauss

      ! Integration weight
      dvolu = refElPol%gauss_weights2D(g)*detJ(g)
      IF (switch%axisym) THEN
      	dvolu = dvolu*xy(g,1)
      END IF

      ! Check if total density is costant
      El_n  = El_n  + ueg(g,1)*2*3.1416*dvolu*phys%lscale**3
      El_nn = El_nn + ueg(g,5)*2*3.1416*dvolu*phys%lscale**3

      ! x and y derivatives of the shape functions
      Nxg = iJ11(g)*refElPol%Nxi2D(g,:) + iJ12(g)*refElPol%Neta2D(g,:)
      Nyg = iJ21(g)*refElPol%Nxi2D(g,:) + iJ22(g)*refElPol%Neta2D(g,:)

      ! Shape functions products
      Ni = refElPol%N2D(g,:)*dvolu
      NNi = tensorProduct(Ni,refElPol%N2D(g,:))                        ! Npel x Npel
      NxNi = tensorProduct(Nxg,Ni)                                     ! Npel x Npel
      NyNi = tensorProduct(Nyg,Ni)                                     ! Npel x Npel
      NNxy = b(g,1)*NxNi + b(g,2)*NyNi                                                        ! Npel x Npel
      NxyzNi = 0.
      NxyzNi(:,:,1) = NxNi
      NxyzNi(:,:,2) = NyNi                                            ! Npel x Npel x 2
      NNbb = (Nxg*b(g,1) + Nyg*b(g,2))*dvolu                             ! Npel x 1
      Nxyzg = 0.
      Nxyzg(:,1) = Nxg*dvolu
      Nxyzg(:,2) = Nyg*dvolu

      ! Divergence of b at the Gauss points
      IF (switch%axisym) THEN
        Nx_ax = Nxg + 1./xy(g,1)*refElPol%N2D(g,:)
      ELSE
        Nx_ax = Nxg
      END IF
      divbg = dot_product(Nx_ax,b_nod(:,1)) + dot_product(Nyg,b_nod(:,2))

      ! Diamagnetic drift !TODO: verify drift intensity in isothermal and non-isothermal cases
      driftg = 0.
      gradbmod = 0.
      gradbmod(1) = dot_product(Nxg,Bmod_nod)
      gradbmod(2) = dot_product(Nyg,Bmod_nod)
      bg = b(g,:)
      call cross_product(bg,gradbmod,driftg)
      driftg = phys%dfcoef*driftg/Bmod(g)

      CALL assemblyVolumeContribution(Auq,Auu,rhs,b(g,:),Psig(g),divbg,driftg,Bmod(g),force(g,:),&
        &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upg(g,:),&
        &ueg(g,:),qeg(g,:),u0eg(g,:,:),xy(g,:),Jtor(g),Vnng)
        
      if (save_tau) then
         v_nn_Vol_el(g,:) = Vnng
      endif
      
    END DO ! END loop in volume Gauss points
    call do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
    deallocate(Auq,Auu,rhs)

  END SUBROUTINE elemental_matrices_volume

  !***************************************************
  ! Interior faces computation in 2D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_int(iel,ifa,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(IN)         :: uef(:,:),uf(:,:)
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,j,k,ieln,indsave(Ng1d)
    real*8                    :: dline,xyDerNorm_g
    real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
    real*8                    :: xyf(Ng1d,ndim)
    real*8                    :: xyDer(Ng1d,ndim)
    real*8                    :: qfg(Ng1d,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl)  :: ind_asf,ind_ash
    real*8                    :: t_g(ndim),n_g(ndim),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8                    :: upgf(Ng1d,phys%npv)
    real*8                    :: tau(Neq,Neq),Vnng(Ndim)
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d),Psig(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)

    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    ! Faces computations
    !***********************************
    NGauss = Ng1d

    ! Indices
    ind_fe = reshape(tensorSumInt((/(i,i=1,neq)/),neq*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)
    ind_fg = reshape(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl*ndim/))

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod

    ! Trace solution at face Gauss points
    IF (Mesh%flipFace(iel,ifa)) THEN
      ufg = matmul(refElPol%N1D,uf((/(i,i=Npfl,1,-1)/),:))
    ELSE
      ufg = matmul(refElPol%N1D,uf)
    END IF


    ! Gauss points position and derivatives
    xyf = matmul(refElPol%N1D,Xfl)
    xyDer = matmul(refElPol%Nxi1D,Xfl)

    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElPol%N1D,Bmod_nod)
    b = matmul(refElPol%N1D,b_nod)
   
    ! Normalaized magnetic flux at Gauss points: PSI
    Psig = matmul(refElPol%N1d,psifl) 

    ! Element solution at face Gauss points
    uefg = matmul(refElPol%N1D,uef)
    ! Gradient solution at face gauss points
    qfg = matmul(refElPol%N1D,qef)

    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac)
    if (save_tau) then
       indsave = (ifa - 1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(5,5,:)
    end if

    if (switch%shockcp.gt.0) then
      auxdiffsc = matmul(refElPol%N1D,Mesh%scdiff_nodes(iel,refElPol%face_nodes(ifa,:)))
      do i=1,Neq
        diff_iso_fac(i,i,:) = diff_iso_fac(i,i,:)+auxdiffsc
      end do
    endif

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

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

      ! Shape functions products
      bn = dot_product(b(g,1:2),n_g)
      NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Nif = refElPol%N1D(g,:)*dline
      Nfbn = bn*refElPol%N1D(g,:)*dline

      ! Compute the stabilization term
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,0.,xyf(g,:),tau)
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g,xyf(g,:),0.,iel,tau)
        ENDIF
      END IF

      ! Assembly local contributions
      CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)

      if (save_tau) then
        DO i = 1,Neq
          tau_save_el((ifa - 1)*Ngauss + g,i) = tau(i,i)
        END DO
        v_nn_Fac_el((ifa -1)*Ngauss + g,:) = Vnng
        xy_g_save_el((ifa - 1)*Ngauss + g,:) = xyf(g,:)
      endif
      

    END DO ! Gauss points
!stop


  END SUBROUTINE elemental_matrices_faces_int

  !***************************************************
  ! Exterior faces computation in 2D
  !***************************************************
  SUBROUTINE elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,Bfl,psifl,qef,uef,uf,diff_nn_Fac_el,v_nn_Fac_el,tau_save_el,xy_g_save_el)
    integer,intent(IN)        :: iel,ifa
    real*8,intent(IN)         :: Xfl(:,:)
    real*8,intent(IN)         :: Bfl(:,:), psifl(:)
    logical,intent(IN)        :: isdir
    real*8,intent(IN)         :: qef(:,:)
    real*8,intent(INOUT)      :: uef(:,:),uf(:,:)
    real*8,intent(out)        :: diff_nn_Fac_el(:),v_nn_Fac_el(:,:),tau_save_el(:,:),xy_g_save_el(:,:)
    integer*4                 :: g,NGauss,i,j,k,indsave(Ng1d)
    real*8                    :: dline,xyDerNorm_g
    real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
    real*8                    :: xyf(Ng1d,ndim)
    real*8                    :: xyDer(Ng1d,ndim)
    real*8                    :: qfg(Ng1d,neq*Ndim)
    integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
    integer*4,dimension(Npfl)  :: ind_asf,ind_ash
    real                      :: isext
    real*8                    :: t_g(ndim),n_g(ndim),bn
    real*8                    :: NNif(Npfl,Npfl),Nif(Npfl),Nfbn(Npfl)
    real*8                    :: tau(Neq,Neq)
    real*8                    :: upgf(Ng1d,phys%npv)
    real*8                    :: Bmod_nod(Npfl),b_nod(Npfl,3),b(Ng1d,3),Bmod(Ng1d), Psig(Ng1d)
    real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)
    real*8                    :: auxdiffsc(Ng1d)
    real*8                    :: Vnng(Ndim)
	
    ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
    ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

    !***********************************
    ! Faces computations
    !***********************************
    NGauss = Ng1d

    ! Indices
    ind_fe = reshape(tensorSumInt((/(i,i=1,neq)/),neq*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl/))
    ind_ff = (ifa - 1)*neq*Npfl + (/(i,i=1,neq*Npfl)/)
    ind_fg = reshape(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifa,:) - 1)),(/neq*Npfl*ndim/))

    !****************************************************
    !                      Magnetic field
    !****************************************************
    ! Magnetic field norm and direction at element nodes
    Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
    b_nod(:,1) = Bfl(:,1)/Bmod_nod
    b_nod(:,2) = Bfl(:,2)/Bmod_nod
    b_nod(:,3) = Bfl(:,3)/Bmod_nod
    ! Magnetic field norm and direction at Gauss points
    Bmod = matmul(refElPol%N1D,Bmod_nod)
    b = matmul(refElPol%N1D,b_nod)

    ! Normalaized magnetic flux at Gauss points: PSI
    Psig = matmul(refElPol%N1D,psifl)

    ! Trace solution at face Gauss points
    xyf = matmul(refElPol%N1D,Xfl)
    IF (isdir) THEN
      CALL analytical_solution(iel,xyf(:,1),xyf(:,2),ufg)
    ELSE
#ifdef PARALL
      IF (Mesh%flipFace(iel,ifa)) THEN
        uf = uf((/(i,i=Npfl,1,-1)/),:)
      ENDIF
      ! TODO: VERIFY IF I NEED TO FLIP ALSO xyf,b and Bmod in this case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif
      ufg = matmul(refElPol%N1D,uf)
    END IF

    ! Element solution at face Gauss points
    uefg = matmul(refElPol%N1D,uef)
    ! Gradient solution at face gauss points
    qfg = matmul(refElPol%N1D,qef)

    ! Compute diffusion at faces Gauss points
    CALL setLocalDiff(xyf,uefg,qfg,diff_iso_fac,diff_ani_fac)
    if (save_tau) then
       indsave = (ifa -1)*Ngauss + (/(i,i=1,Ngauss)/)
       diff_nn_Fac_el(indsave) = diff_iso_fac(5,5,:)
    end if
	
    if (switch%shockcp.gt.0) then
      auxdiffsc = matmul(refElPol%N1D,Mesh%scdiff_nodes(iel,refElPol%face_nodes(ifa,:)))
      do i=1,Neq
        diff_iso_fac(i,i,:) = diff_iso_fac(i,i,:)+auxdiffsc
      end do
    endif

    ! Physical variables at face Gauss points
    CALL cons2phys(ufg,upgf)

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

      ! Compute the stabilization term
      isext = 1.
#ifdef PARALL
      IF (Mesh%boundaryFlag(Mesh%F(iel,ifa) - Mesh%Nintfaces) .eq. 0) THEN
        isext = 0.
      END IF
#endif
      tau = 0.
      IF (numer%stab == 1) THEN
        ! Constant stabilization
        DO i = 1,Neq
          tau(i,i) = numer%tau(i)
        END DO
      ELSE
        ! Non constant stabilization
        ! Compute tau in the Gauss points
        IF (numer%stab < 6) THEN
          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),qfg(g,:),b(g,:),n_g,iel,ifa,isext,xyf(g,:),tau)
        ELSE
          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b(g,:),n_g,xyf(g,:),isext,iel,tau)
        ENDIF
      END IF

      ! Shape functions products
      bn = dot_product(b(g,1:2),n_g)
      NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
      Nfbn = bn*refElPol%N1D(g,:)*dline
      Nif = refElPol%N1D(g,:)*dline

      !call displayMatrix(NNif)
      !call displayVector(Nif)
      !call displayMatrix(uf)
      !call displayMatrix(ufg)
      !stop

      ! Assembly local contributions
#ifdef PARALL
      IF (Mesh%boundaryFlag(Mesh%F(iel,ifa) - Mesh%Nintfaces) .eq. 0) THEN
        ! Ghost face: assembly it as interior
        CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)
      ELSE
        CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau)
      ENDIF
#else
      CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b(g,:),Bmod(g),Psig(g),&
        n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nif,Nfbn,ufg(g,:),upgf(g,:),qfg(g,:),tau,Vnng)
#endif
      if (save_tau) then
        DO i = 1,Neq
          tau_save_el((ifa - 1)*Ngauss + g,i) = tau(i,i)
        END DO
        v_nn_Fac_el((ifa - 1)*Ngauss + g,:) = Vnng
        xy_g_save_el((ifa - 1)*Ngauss + g,:) = xyf(g,:)
      endif

    END DO ! Gauss points

  END SUBROUTINE elemental_matrices_faces_ext

  !*******************************************
  !           AUXILIARY ROUTINES
  !*******************************************

  !*****************************************
  ! Set permutations for flipping faces
  !****************************************
  SUBROUTINE set_permutations(n,m,perm)
    integer,intent(IN)  :: n,m
    integer,intent(OUT) :: perm(:)
    integer              :: i
    integer              :: temp(m,n/m),templr(m,n/m)

    IF (mod(n,m) .ne. 0) then
      WRITE (6,*) 'Error! n must be a multiple of m'
      STOP
    END IF

    templr = 0
    temp = reshape((/(i,i=1,n)/),(/m,n/m/))
    DO i = 1,n/m
      templr(:,i) = temp(:,n/m - i + 1)
    END DO
    perm = reshape(templr,(/n/))
  END SUBROUTINE set_permutations

#endif
!TOR3D

  !********************************************************************************************
  !
  !
  !                             ROUTINES FOR 2D AND 3D COMPUTATIONS
  !
  !
  !********************************************************************************************
  SUBROUTINE setTimeIntegrationCoefficients(ktis)
    real*8,intent(out) :: ktis(:)
    integer :: it

    ktis = 0.

    if (time%ik .lt. time%tis) then
      it = time%ik
    else
      it = time%tis
    end if

    SELECT CASE (it)
    case (1)
      ktis(1) = 1.
      ktis(2) = 1.
    case (2)
      ktis(1) = 1.5
      ktis(2) = 2
      ktis(3) = -0.5
    case (3)
      ktis(1) = 11./6.
      ktis(2) = 3.
      ktis(3) = -1.5
      ktis(4) = 1./3.
    case (4)
      ktis(1) = 25./12.
      ktis(2) = 4.
      ktis(3) = -3.
      ktis(4) = 4./3.
      ktis(4) = -0.25
    case (5)
      ktis(1) = 137./60.
      ktis(2) = 5.
      ktis(3) = -5.
      ktis(4) = 10./3.
      ktis(5) = -1.25
      ktis(6) = 0.2
    case (6)
      ktis(1) = 147./60.
      ktis(2) = 6.
      ktis(3) = -7.5
      ktis(4) = 20./3.
      ktis(5) = -3.75
      ktis(6) = 1.2
      ktis(7) = -1.6
    case default
      write (6,*) 'Formula not available'
      stop
    END SELECT
  END SUBROUTINE setTimeIntegrationCoefficients

  !********************************************************************
  !
  !         ASSEMBLY VOLUME CONTRIBUTION
  !
  !********************************************************************
  SUBROUTINE assemblyVolumeContribution(Auq,Auu,rhs,b3,psi,divb,drift,Bmod,f,&
      &ktis,diffiso,diffani,Ni,NNi,Nxyzg,NNxy,NxyzNi,NNbb,upe,ue,qe,u0e,xy,Jtor,Vnng)
    real*8,intent(inout)      :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
    real*8,intent(IN)         :: b3(:),psi,divb,drift(:),f(:),ktis(:),Bmod
    real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
    real*8,intent(IN)         :: Ni(:),NNi(:,:),Nxyzg(:,:),NNxy(:,:),NxyzNi(:,:,:),NNbb(:)
    real*8,intent(IN)         :: upe(:),ue(:),xy(:),Jtor
    real*8,intent(INOUT)      :: u0e(:,:),Vnng(:)
    real*8,intent(IN)         :: qe(:)
    integer*4                 :: i,j,k,iord,ii,alpha,beta,z
    integer*4,dimension(Npel) :: ind_i,ind_j,ind_k
    real*8,dimension(neq,neq) :: A
    real*8,dimension(neq,Ndim):: APinch
    real*8                    :: kcoeff
    real*8                    :: Qpr(Ndim,Neq),exb(3),bb(3)
    real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
    real*8                    :: qq(3,Neq),b(Ndim)
    real*8                    :: grad_n(3),gradpar_n
    real*8                    :: auxvec(Neq)
#ifdef TEMPERATURE
    real*8,dimension(neq,neq) :: GG
    real*8                    :: Telect
    real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
    real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
    real*8                    :: W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
    real*8                    :: Sohmic,dSohmic_dU(Neq) ! Ohmic heating
    real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
    real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#endif
#ifdef NEUTRAL
    real*8                    :: Ax(Neq,Neq),Ay(Neq,Neq)
    real*8                    :: niz,nrec,fGammacx,fGammarec
    real*8                    :: dniz_dU(Neq),dnrec_dU(Neq),dfGammacx_dU(Neq),dfGammarec_dU(Neq)
#ifdef TEMPERATURE
    real*8                    :: sigmaviz,sigmavrec,sigmavcx,Tloss,Tlossrec,fEiiz,fEirec,fEicx
    real*8                    :: dsigmaviz_dU(Neq),dsigmavrec_dU(Neq),dsigmavcx_dU(Neq),dTloss_dU(Neq),dTlossrec_dU(Neq)
    real*8                    :: dfEiiz_dU(Neq),dfEirec_dU(Neq),dfEicx_dU(Neq)
#ifdef NEUTRALP
    real*8                    :: Dnn,Dpn,Alphanp,Betanp,GammaLim,Gammaredpn,Tmin
    real*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),dDpn_dU(Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq)
#endif
#endif
    real*8                    :: Sn(Neq,Neq),Sn0(Neq)
#endif



    real*8 :: kmult(size(Auq,1),size(Auq,2))



    b = b3(1:Ndim)

    bb = 0.
    bb = b3

    ! Jacobian for convection term
    CALL jacobianMatrices(ue,A)
    
    ! Jacobian for pinch term
    CALL computePinch(b,psi,APinch)

    ! Compute Q^T^(k-1)
    Qpr = reshape(qe,(/Ndim,Neq/))

    ! Split diffusion matrices/vectors for the momentum equation
    CALL compute_W2(ue,W2)
    CALL compute_dW2_dU(ue,dW2_dU)
    QdW2 = matmul(Qpr,dW2_dU)

    qq = 0.
    qq(1:Ndim,:) = Qpr

    ! Perpendicular gradient of density
    grad_n=qq(:,1)
    gradpar_n = dot_product(grad_n,b3)

#ifdef TEMPERATURE
    ! Jacobian for the curvature term
    CALL GimpMatrix(ue,divb,GG)

    ! Compute V(U^(k-1))
    call computeVi(ue,Vveci)
    call computeVe(ue,Vvece)

    ! Compute dV_dU (k-1)
    call compute_dV_dUi(ue,dV_dUi)
    call compute_dV_dUe(ue,dV_dUe)

    ! Compute Alpha(U^(k-1))
    Alphai = computeAlphai(ue)
    Alphae = computeAlphae(ue)

    ! Compute dAlpha/dU^(k-1)
    call compute_dAlpha_dUi(ue,dAlpha_dUi)
    call compute_dAlpha_dUe(ue,dAlpha_dUe)

    gmi = dot_product(matmul(Qpr,Vveci),b)    ! scalar
    gme = dot_product(matmul(Qpr,Vvece),b)    ! scalar
    Taui = matmul(Qpr,dV_dUi)                 ! Ndim x Neq
    Taue = matmul(Qpr,dV_dUe)                 ! Ndim x Neq

    ! Parallel current term
    ! Compute W(U^(k-1))
    call compute_W(ue,W)
    ! Compute dW_dU(U^(k-1))
    call compute_dW_dU(ue,dW_dU)

    ! Split diffusion matrices/vectors for the energies equations
    CALL compute_W3(ue,W3)
    CALL compute_dW3_dU(ue,dW3_dU)
    QdW3 = matmul(Qpr,dW3_dU)

    CALL compute_W4(ue,W4)
    CALL compute_dW4_dU(ue,dW4_dU)
    QdW4 = matmul(Qpr,dW4_dU)

    ! Temperature exchange terms
    ! s(U^(k-1))
    call compute_S(ue,s)
    ! ds_du(U^(k-1))
    call compute_dS_dU(ue,ds_dU)

    !Ohmic Heating
    IF (switch%ohmicsrc) THEN
      !Compute Sohmic(U^(k-1))
      call compute_Sohmic(ue,Sohmic)
      !Compute dSohmic_dU(U^(k-1))
      call compute_dSohmic_dU(ue,dSohmic_dU)
    ENDIF

    Zet = matmul(Qpr,dW_dU)       ! Ndim x Neq
    
#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(ue,Vpn)
	   ! Compute dVpn_dU(U^(k-1))
	   CALL compute_dVpn_dU(ue,dVpn_dU)
	   gmpn = matmul(Qpr,Vpn)                      ! Ndim x 1
	   Taupn = matmul(Qpr,dVpn_dU)                 ! Ndim x Neq
	   ! Compute Dpn(U^(k-1))
    CALL computeDpn(ue,Qpr,Vpn,Dpn)
    ! Compute dDpn_dU(U^(k-1))
    CALL compute_dDpn_dU(ue,Qpr,Vpn,dDpn_dU) 
    ! Reduce Grad Pn for low collision regime 
    ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
    Gammaredpn = 1.
    Tmin = 0.2/simpar%refval_temperature
    IF (Tmin/upe(7) .le. 1.) Gammaredpn = Gammaredpn*Tmin/upe(7)
    Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*Dpn 
    ! Compute Gammaredpn(U^(k-1))
    !CALL computeGammared(ue,Gammaredpn) 
    !gmipn = matmul(Qpr,Vveci) 
    !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
    ! Set Grad Ti = 0. for low collision regime 
    ! (back to diffusion equation for neutral density)
    !CALL computeAlphaCoeff(ue,Qpr,Vpn,Alphanp)  
    !CALL computeBetaCoeff(ue,Qpr,Vpn,Betanp)  
    !Dpn = Alphanp*Dpn
    !dDpn_dU = Alphanp*dDpn_dU
    !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*Dpn 
    !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = Alphanp*Dpn !0.
    !   dDpn_dU = Alphanp*dDpn_dU !0.
    !   Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upe(7)*phys%diff_nn
    !   END IF
    ! Set Gamma Convective = cs_n*n_n for low collision regime 
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = 0.
    !   dDpn_dU = 0.
    !   CALL jacobianMatricesNP(ue,Anp)
    !ELSE
    !   Anp = 0.
    !END IF
#endif
#endif



#ifdef NEUTRAL
    !Neutral Source Terms needed in the plasma and neutral density equations
    call compute_niz(ue,niz)
    call compute_nrec(ue,nrec)
    call compute_dniz_dU(ue,dniz_dU)
    call compute_dnrec_dU(ue,dnrec_dU)
#ifdef TEMPERATURE
    call compute_sigmaviz(ue,sigmaviz)
    call compute_sigmavrec(ue,sigmavrec)
    call compute_dsigmaviz_dU(ue,dsigmaviz_dU)
    call compute_dsigmavrec_dU(ue,dsigmavrec_dU)
#endif
    !Neutral Source Terms needed in the plasma momentum equation
    call compute_fGammacx(ue,fGammacx)
    call compute_dfGammacx_dU(ue,dfGammacx_dU)
    call compute_fGammarec(ue,fGammarec)
    call compute_dfGammarec_dU(ue,dfGammarec_dU)
#ifdef TEMPERATURE
    call compute_sigmavcx(ue,sigmavcx)
    call compute_dsigmavcx_dU(ue,dsigmavcx_dU)
    !Neutral Source Terms needed in the ion energy equation
    call compute_fEiiz(ue,fEiiz)
    call compute_dfEiiz_dU(ue,dfEiiz_dU)
    call compute_fEirec(ue,fEirec)
    call compute_dfEirec_dU(ue,dfEirec_dU)
    call compute_fEicx(ue,fEicx)
    call compute_dfEicx_dU(ue,dfEicx_dU)
    !Neutral Source Terms needed in the electron energy equation
    call compute_Tloss(ue,Tloss)
    call compute_dTloss_dU(ue,dTloss_dU)
    call compute_Tlossrec(ue,Tlossrec)
    call compute_dTlossrec_dU(ue,dTlossrec_dU)
#endif

    !Assembly the matrix for neutral sources
#ifdef TEMPERATURE
    call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Tloss,dTloss_dU,Tlossrec,dTlossrec_dU,Sn,Sn0)
#else
    call assemblyNeutral(ue,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
#endif
!NEUTRAL

    ! Assembly local matrix
    ! Loop in equations
    DO i = 1,Neq
      ! ind_i = i + ind_ass
      IF (.not. switch%steady) THEN
        ! Time derivative contribution
#ifdef VORTICITY
        ! In the vorticity model I don't assemble the mass matrix for the potential equation
        IF (i .ne. 4) THEN
#endif
          z = i+(i-1)*Neq
          Auu(:,:,z)= Auu(:,:,z)+ ktis(1)*NNi/time%dt
#ifdef VORTICITY
        END IF
#endif
      END IF
#ifndef TEMPERATURE
      IF (i == 2) THEN
        ! Curvature contribution (isothermal)
        z = i+(i-2)*Neq
        if (switch%logrho) then
          Auu(:,:,z) = Auu(:,:,z) - phys%a*divb*upe(1)*NNi
          rhs(:,i)=rhs(:,i)+Ni*phys%a*divb*upe(1)*(1-ue(1))
        else
          Auu(:,:,z) = Auu(:,:,z) - phys%a*divb*NNi
        endif

        ! split diffusion momentum equation (LU) (isothermal)
        DO j = 1,Neq
           z = i+(j-1)*Neq
           DO k = 1,Ndim
              Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*QdW2(k,j))
           END DO
              Auu(:,:,z) = Auu(:,:,z) - (dot_product(QdW2(:,j),b))*NNxy
        END DO
      END IF
#ifdef VORTICITY
      IF (switch%driftdia .and. i.ne.4) THEN
#else
        IF (switch%driftdia) THEN
#endif
          ! B x GradB drift (isothermal)
          DO k = 1,Ndim
            z = i+(i-1)*Neq
            Auu(:,:,z)= Auu(:,:,z) +transpose(NxyzNi(:,:,k))*drift(k)
          END DO
        END IF
#else
        IF (i == 2) THEN
          DO j = 1,Neq
            z = i+(j-1)*Neq
            ! Curvature contribution (non-isothermal)
            Auu(:,:,z) = Auu(:,:,z) - GG(i,j)*NNi

            ! split diffusion momentum equation (LU) (non-isothermal)
            DO k = 1,Ndim
               Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*QdW2(k,j))
            END DO
            Auu(:,:,z) = Auu(:,:,z) - (dot_product(QdW2(:,j),b))*NNxy
          END DO
        END IF

        IF (switch%driftdia) THEN
          Telect = upe(8)
          ! B x GradB drift (non-isothermal)
          DO k = 1,Ndim
            z = i+(i-1)*Neq
            Auu(:,:,z)= Auu(:,:,z) + Telect*transpose(NxyzNi(:,:,k))*drift(k)
          END DO
        END IF

        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
            Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq)+coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),b)))*NNxy + &
              &(dot_product(Zet(:,j),b) + ds_dU(j))*NNi
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z) = Auq(:,:,z)+ coefi*Alphai*Vveci(j)*b(k)*NNxy
              IF (j == 4) THEN
                Auq(:,:,z) = Auq(:,:,z)+W*NNi*b(k)
              END IF
            END DO
          END DO
          rhs(:,i) = rhs(:,i) + coefi*Alphai*(dot_product(matmul(transpose(Taui),b),ue))*NNbb + s*Ni
        ELSEIF (i == 4) THEN
          DO j = 1,4
            z = i+(j-1)*Neq
            Auu(:,:,z)=Auu(:,:,z)+coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),b)))*NNxy - (dot_product(Zet(:,j),b) + ds_dU(j))*NNi
            IF (switch%ohmicsrc) THEN
              Auu(:,:,z) = Auu(:,:,z) - dSohmic_dU(j)*(Jtor**2)*NNi
            ENDIF
            DO k = 1,Ndim
              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
              Auq(:,:,z)=Auq(:,:,z)+coefe*Alphae*Vvece(j)*b(k)*NNxy
              IF (j == 4) THEN
                Auq(:,:,z)=Auq(:,:,z)- W*NNi*b(k)
              END IF
            END DO
            ! split diffusion electron energy equation (LU)
            z = i+(j-1)*Neq
            DO k = 1,Ndim
              Auu(:,:,z) = Auu(:,:,z) + (NxyzNi(:,:,k)*QdW4(k,j))
            END DO
            Auu(:,:,z) = Auu(:,:,z) - (dot_product(QdW4(:,j),b))*NNxy
          END DO
          rhs(:,i) = rhs(:,i)+coefe*Alphae*(dot_product(matmul(transpose(Taue),b),ue))*NNbb - s*Ni
          IF (switch%ohmicsrc) THEN
            rhs(:,i) = rhs(:,i) + Sohmic*(Jtor**2)*Ni
          ENDIF
#ifdef NEUTRALP		  
		      ELSEIF (i == 5) THEN
		         DO j = 1,5
			           DO k = 1,Ndim
			              z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                    Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq) + (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*NxyzNi(:,:,k) 
                    !Auu(:,:,i+(j-1)*Neq) = Auu(:,:,i+(j-1)*Neq) - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*NxyzNi(:,:,k) 
                    Auq(:,:,z) = Auq(:,:,z) + Dpn*Vpn(j)*NxyzNi(:,:,k)
                    !Auq(:,:,z) = Auq(:,:,z) - Gammaredpn*Dpn*Vveci(j)*NxyzNi(:,:,k)
                    IF (j == 5) THEN
                       Auq(:,:,z) = Auq(:,:,z) + Dnn*NxyzNi(:,:,k)
                    END IF
				         END DO
           END DO
           DO k = 1,Ndim
              rhs(:,i) = rhs(:,i) + dot_product(dDpn_dU,ue)*gmpn(k)*Nxyzg(:,k)
              !rhs(:,i) = rhs(:,i) - Gammaredpn*(Dpn*dot_product(Taui(k,:),ue) + dot_product(dDpn_dU,ue)*gmipn(k))*Nxyzg(:,k)
           END DO
#endif
		END IF
#endif
       
	! Convection contribution
        DO j = 1,Neq
           z = i+(j-1)*Neq
           Auu(:,:,z)= Auu(:,:,z) - A(i,j)*NNxy
#ifdef NEUTRAL
!#ifdef NEUTRALP
!          IF (i == 5) Auu(:,:,z) = Auu(:,:,z) - Anp(j)*NNxy
!#endif
          !Sources
          Auu(:,:,z) = Auu(:,:,z) + Sn(i,j)*NNi
#endif
        END DO

        ! Pinch contribution
        z = i+(i-1)*Neq
        Auu(:,:,z) =  Auu(:,:,z) - (APinch(i,1)*NxyzNi(:,:,1) + APinch(i,2)*NxyzNi(:,:,2))

#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
        if (switch%logrho) then
          call logrhojacobianVector(ue,upe,auxvec)
          rhs(:,i)=rhs(:,i)+NNbb*auxvec(i)

          do j = 1,Neq
            if (i==1 .and. j==1) then
              z = i+(j-1)*Neq
              Auu(:,:,z)=Auu(:,:,z)+upe(2)*transpose(NNxy) ! TODO: this is a first order linearization!!!
            endif
          end do
        endif
#endif

        DO k = 1,Ndim
        ! split diffusion contributions (LQ)
	        IF (i==2) THEN
            DO j = 1,Neq
                z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auq(:,:,z) = Auq(:,:,z) + W2(j)*(NxyzNi(:,:,k) -NNxy*b(k))
            END DO
#ifdef TEMPERATURE
          ELSEIF(i==3) THEN
            DO j = 1,Neq
	         z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                 Auq(:,:,z) = Auq(:,:,z) + W3(j)*(NxyzNi(:,:,k) -NNxy*b(k))
           END DO
          ELSEIF(i==4) THEN
            DO j = 1,Neq
                z = i+(k-1)*Neq+(j-1)*Neq*Ndim
                Auq(:,:,z) = Auq(:,:,z) + W4(j)*(NxyzNi(:,:,k) -NNxy*b(k))
            END DO
#endif
           ENDIF

          ! Diagonal terms for perpendicular diffusion
          z = i+(k-1)*Neq+(i-1)*Neq*Ndim
          kmult = diffiso(i,i)*NxyzNi(:,:,k) - diffani(i,i)*NNxy*b(k)

#ifdef VORTICITY
          if (i==4) then
            kmult = kmult/ue(1)
          endif
#endif

          Auq(:,:,z)=Auq(:,:,z)+kmult

#ifndef TEMPERATURE
          ! Added term for n=exp(x) change of variable \Grad \chi **2
          if (switch%logrho) then
            if (i==1) then
              Auq(:,:,z)=Auq(:,:,z)-2*(diffiso(i,i)*grad_n(k) - diffani(i,i)*b(k)*gradpar_n)*NNi
              rhs(:,i)=rhs(:,i)-Ni*(diffiso(i,i)*grad_n(k)*grad_n(k) - diffani(i,i)*gradpar_n**2/Ndim )
            endif
          endif
#endif

#ifdef VORTICITY
          if (switch%bxgradb) then
            ! B x GradB current in the vorticity equation
            IF (i==3) THEN
              ii =1
              z = i+(ii-1)*Neq
              if (switch%logrho) then
                Auu(:,:,z)= Auu(:,:,z) +2*upe(1)*transpose(NxyzNi(:,:,k))*drift(k) ! TODO: first order linearization
              else
                Auu(:,:,z)= Auu(:,:,z) +2*transpose(NxyzNi(:,:,k))*drift(k)
              endif
            ENDIF
          endif
          IF (switch%driftexb .and. i .ne. 4) THEN
            ! ExB terms
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
            call ijk_cross_product(k,alpha,beta)
            ii = 4
            z = i+(k-1)*Neq+(ii-1)*Neq*Ndim
            Auq(:,:,z)=Auq(:,:,z)+kcoeff*(NxyzNi(:,:,alpha)*b3(beta) - NxyzNi(:,:,beta)*b3(alpha))*ue(i)
            call cross_product(qq(:,ii),bb,exb)
            z = i+(i-1)*Neq
            Auu(:,:,z)=Auu(:,:,z)+kcoeff*exb(k)*NxyzNi(:,:,k)
            rhs(:,i)=rhs(:,i)+kcoeff*exb(k)*ue(i)*Nxyzg(:,k)
          ENDIF

          ! Non-diagonal terms for perpendicular diffusion
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal already assembled
            IF (abs(diffiso(i,ii)) < 1e-12 .and. abs(diffani(i,ii)) < 1e-12) CYCLE
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)

            !IF ((i == 3 .or. i == 4) .and. ii == 1) then
            IF ((i == 3) .and. ii == 1) then
              z = i+(ii-1)*Neq
              kcoeff = 1./ue(1)
              Auu(:,:,z)=Auu(:,:,z)-kcoeff**2*(diffiso(i,ii)*Qpr(k,ii)*NxyzNi(:,:,k) - diffani(i,ii)*Qpr(k,ii)*b(k)*NNxy)
              rhs(:,i)=rhs(:,i)-kcoeff*(diffiso(i,ii)*Qpr(k,ii)*Nxyzg(:,k) - diffani(i,ii)*Qpr(k,ii)*b(k)*NNbb)
            ENDIF
            z=i+(k-1)*Neq+(ii-1)*Ndim*Neq
            Auq(:,:,z)=Auq(:,:,z)+kcoeff*diffiso(i,ii)*NxyzNi(:,:,k) - kcoeff*diffani(i,ii)*NNxy*b(k)

            !write(6,*) "kcoeff",kcoeff
            !write(6,*) "i:",i,"ii:",ii, "diffiso(i,ii)",diffiso(i,ii)
            !write(6,*) "i:",i,"ii:",ii, "diffani(i,ii)",diffani(i,ii)
            !call HDF5_save_matrix(kcoeff*diffiso(i,ii)*NxyzNi(:,:,k) - kcoeff*diffani(i,ii)*NNxy*b(k),'fava')
            !stop

          END DO
#endif
        END DO ! loop in k: 1-Ndim

#ifdef VORTICITY
        ! The vorticity is the source term in the potential equation
        IF (i == 4) THEN
          j=3
          z = i+(j-1)*Neq
          Auu(:,:,z)=Auu(:,:,z) +NNi
          !rhs(:,i)=rhs(:,i)-ue(j)*Ni ! doens't work very well like this
        ENDIF
#endif

#ifdef VORTICITY
        !if (switch%testcase.eq.7 .and. switch%logrho .and. i.eq.1 .and. upe(1).gt.1) then
        !  rhs(:,i)=rhs(:,i)-100*(upe(1)-1.)*Ni
        !endif
        if (switch%testcase.eq.7) then
          if ( (xy(1)-geom%R0)/phys%lscale .gt. 0.4 ) then
            ! Implicit sources to take into account parallel losses
            if (switch%logrho) then
              if (i==1) then
                rhs(:,i)=rhs(:,i)-phys%diagsource(i)*Ni
              else if (i==3) then
                j=4
                z = i+(j-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
              endif
            else
              if (i==1) then
                z = i+(i-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
              else if (i==3) then
                j=4
                z = i+(j-1)*Neq
                Auu(:,:,z)=Auu(:,:,z) +phys%diagsource(i)*NNi
              endif
            endif
          endif
        endif
#endif
      END DO ! Loop in equations

      ! Assembly RHS
      IF (.not. switch%steady) THEN
#ifdef VORTICITY
        u0e(4,:) = 0.
#endif
        DO iord = 1,time%tis
          ! Time derivative contribution
          rhs=rhs+ktis(iord + 1)*tensorProduct(Ni,u0e(:,iord))/time%dt
        END DO
      END IF

      ! Linear body force contribution
      rhs = rhs+tensorProduct(Ni,f)
#ifdef NEUTRAL
      rhs = rhs-tensorProduct(Ni,Sn0)
#endif
!#ifdef NEUTRALP
!      rhs = rhs+tensorProduct(Ni,fth)
!#endif
    END SUBROUTINE assemblyVolumeContribution

    !********************************************************************
    !
    !         ASSEMBLY INTERIOR FACES CONTRIBUTION
    !
    !********************************************************************
    SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,Bmod,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8,intent(IN)         :: b3(:),n(:),Bmod, psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:),upf(:)
      real*8,intent(IN)         :: qf(:)
      real*8,optional,intent(INOUT) :: tau(:,:),Vnng(:)
      real*8                     :: kcoeff
      real*8                     :: b(Ndim)
      integer*4,optional         :: ifa
      integer*4                  :: i,j,k,ii,alpha,beta
      integer*4,dimension(size(ind_asf))  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
      real*8                    :: auxvec(Neq)
      real*8                    :: nn(3),qq(3,Neq),bb(3)
      real*8                    :: bn,kmult(size(ind_asf),size(ind_asf)),kmultf(size(ind_asf))
      real*8                    :: Qpr(Ndim,Neq),exb(3)
      real*8                    :: W2(Neq),dW2_dU(Neq,Neq),QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: W3(Neq),dW3_dU(Neq,Neq),QdW3(Ndim,Neq)
      real*8                    :: W4(Neq),dW4_dU(Neq,Neq),QdW4(Ndim,Neq)
#ifdef NEUTRALP
      real*8                    :: Dnn,Dpn,GammaLim,Alphanp,Betanp,Gammaredpn,Tmin
   	  real*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq),dDpn_dU(Neq)
#endif
#endif

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
      bn = dot_product(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian for pinch term
      CALL computePinch(b,psi,APinch)  
   
      ! Compute Q^T^(k-1)
      Qpr = reshape(qf,(/Ndim,Neq/))

      ! Split diffusion vector/matrix for momentum equation
      CALL compute_W2(uf,W2)
      CALL compute_dW2_dU(uf,dW2_dU)
      QdW2 = matmul(Qpr,dW2_dU)

      nn = 0.
      qq = 0.
      nn(1:Ndim) = n
      qq(1:Ndim,:) = Qpr
      
#ifdef TEMPERATURE
      ! Compute V(U^(k-1))
      call computeVi(uf,Vveci)
      call computeVe(uf,Vvece)

      ! Compute dV_dU (k-1)
      call compute_dV_dUi(uf,dV_dUi)
      call compute_dV_dUe(uf,dV_dUe)

      ! Split diffusion vector/matrix for the energies equations
      CALL compute_W3(uf,W3)
      CALL compute_dW3_dU(uf,dW3_dU)
      QdW3 = matmul(Qpr,dW3_dU)

      CALL compute_W4(uf,W4)
      CALL compute_dW4_dU(uf,dW4_dU)
      QdW4 = matmul(Qpr,dW4_dU)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(uf)
      Alphae = computeAlphae(uf)

      ! Compute dAlpha/dU^(k-1)
      call compute_dAlpha_dUi(uf,dAlpha_dUi)
      call compute_dAlpha_dUe(uf,dAlpha_dUe)

      gmi = dot_product(matmul(Qpr,Vveci),b)  ! scalar
      gme = dot_product(matmul(Qpr,Vvece),b)
      Taui = matmul(Qpr,dV_dUi)      ! 2x3
      Taue = matmul(Qpr,dV_dUe)      ! 2x3
#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(uf,Vpn)
    gmpn = matmul(Qpr,Vpn)                     ! Ndim x 1
	  ! Compute dVpn-dU(U^(k-1))
	  CALL compute_dVpn_dU(uf,dVpn_dU)
	  Taupn = matmul(Qpr,dVpn_dU)                 ! Ndim x Neq
	  ! Compute Dpn(U^(k-1))
    CALL computeDpn(uf,Qpr,Vpn,Dpn)
    ! Compute dDpn_dU(U^(k-1))
    CALL compute_dDpn_dU(uf,Qpr,Vpn,dDpn_dU)
    ! Reduce Grad Pn for low collision regime 
    ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
    Gammaredpn = 1.
    Tmin = 0.2/simpar%refval_temperature
    IF (Tmin/upf(7) .le. 1.) Gammaredpn = Gammaredpn*Tmin/upf(7)
    Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
    ! Comput Gammaredpn(U^(k-1))
    !CALL computeGammared(uf,Gammaredpn)
    !gmipn = matmul(Qpr,Vveci)
    !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
    ! Set Grad Ti = 0. for low collision regime 
    ! (back to diffusion equation for neutral density)
    !CALL computeAlphaCoeff(uf,Qpr,Vpn,Alphanp)  
    !CALL computeBetaCoeff(uf,Qpr,Vpn,Betanp)  
    !Dpn = Alphanp*Dpn
    !dDpn_dU = Alphanp*dDpn_dU
    !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
    !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = Alphanp*Dpn !0.
    !   dDpn_dU = Alphanp*dDpn_dU !0.
    !   Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*phys%diff_nn
    !   END IF
    ! Set Gamma Convective = cs_n*n_n for low collision regime 
    !IF (Dpn .gt. phys%diff_nn) THEN
    !   Dpn = 0.
    !   dDpn_dU = 0.
    !   CALL jacobianMatricesNP(uf,Anp)
    !ELSE
    !   Anp = 0.
    !END IF
#endif
#endif

      ! Assembly local matrix
      DO i = 1,Neq
        ind_if = ind_asf + i
        DO k = 1,Ndim
          ind_kf = ind_ash + k + (i - 1)*Ndim
          kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))

#ifdef VORTICITY
          if (i==4) then
            kmult=kmult/uf(1)
          endif
#endif
          ! Diagonal terms for perpendicular diffusion
          elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
          elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult

          ! Split diffusion contribution
          IF(i == 2) THEN
            ! Assembly LQ momentum equatuation
            j = 1
            ind_kf = ind_ash + k + (j - 1)*Ndim
             kmult = NNif*W2(j)*(n(k) - bn*b(k))
             elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
             elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
           !Assembly LU momentum equation
             DO j=1,Neq
               ind_jf = ind_asf+j
               kmult = QdW2(k,j)*NNif*(n(k)-bn*b(k))
               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
               elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
             END DO
#ifdef TEMPERATURE
         ELSEIF(i==3) THEN
            ! Assembly LQ ion energy
            DO j=1,Neq
              ind_kf = ind_ash + k + (j - 1)*Ndim
               kmult = NNif*W3(j)*(n(k) - bn*b(k))
               elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            END DO
           ! Assembly LU ion energy
             DO j=1,Neq
               ind_jf = ind_asf+j
               kmult = QdW3(k,j)*NNif*(n(k)-bn*b(k))
               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
               elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
             END DO
          ELSEIF(i == 4) THEN
             ! Assembly LQ electron energy
             j = 1
             ind_kf = ind_ash + k + (j - 1)*Ndim
              kmult = NNif*W4(j)*(n(k) - bn*b(k))
              elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
              elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
             !Assembly LU electron energy
              DO j=1,Neq
                ind_jf = ind_asf+j
                kmult = QdW4(k,j)*NNif*(n(k)-bn*b(k))
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
              END DO
#endif
          ENDIF
#ifdef VORTICITY

          IF (switch%driftexb .and. i .ne. 4) THEN
            ! ExB terms
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
            ii = 4
            call ijk_cross_product(k,alpha,beta)
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kmult = kcoeff*NNif*(nn(alpha)*b3(beta) - nn(beta)*b3(alpha))*uf(i)
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
            call cross_product(qq(:,ii),bb,exb)
            kmult = kcoeff*exb(k)*NNif*b(k)
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
            elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) - kmult
            kmultf = kcoeff*exb(k)*uf(i)*Nif*b(k)
            elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf

            !tau(i,i) = tau(i,i) + abs(kcoeff*exb(k)*b(k))
            tau(i,i) = 100
          ENDIF
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal alredy assembled
            IF (abs(diffiso(i,ii)) < 1e-12 .and. abs(diffani(i,ii)) < 1e-12) CYCLE
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)

            !IF ((i == 3 .or. i == 4) .and. ii == 1) then
            IF ((i == 3 ) .and. ii == 1) then
              ! Non-linear term in the vorticity equation (\Grad// n/n b)
              ind_jf = ind_asf + ii
              kcoeff = 1./uf(1)
              kmult = kcoeff**2*(diffiso(i,ii)*Qpr(k,1)*n(k)*NNif - diffani(i,ii)*Qpr(k,1)*b(k)*NNif*bn)
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
              elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult
              kmultf = kcoeff*(diffiso(i,ii)*Qpr(k,1)*n(k)*Nif - diffani(i,ii)*Qpr(k,1)*b(k)*Nfbn)
              elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) + kmultf
              elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) + kmultf
            ENDIF
            kmult = NNif*kcoeff*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii))
            elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          END DO
#endif
        END DO ! k-loop

        ! Convection contribution
        DO j = 1,Neq
          ind_jf = ind_asf + j
          kmult = bn*A(i,j)*NNif
!#ifdef NEUTRALP
!          IF (i == 5) kmult = kmult + bn*Anp(j)*NNif
!#endif
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
          elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult
!#ifdef TEMPERATURE
!#ifdef NEUTRAL
!          !X component neutral convective velocity
!          kmult = n(1)*Ax(i,j)*NNif
!          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!          elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
!          !Y component neutral convective velocity
!          kmult = n(2)*Ay(i,j)*NNif
!          elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!          elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
!#endif
!#endif
        END DO ! j-loop
        ! Pinch contribution
        elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif
        elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif

#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
        if (switch%logrho) then
          call logrhojacobianVector(uf,upf,auxvec)
          kmultf = Nfbn*auxvec(i)
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel)-kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) -kmultf
        endif
#endif

#ifdef TEMPERATURE
        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
            kmult = coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),b)))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
          kmultf = coefi*Alphai*(dot_product(matmul(transpose(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
            kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),b)))*NNif*bn
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
              elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
          kmultf = coefe*Alphae*(dot_product(matmul(transpose(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
#ifdef NEUTRALP
       ELSEIF (i == 5) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*n(k)*NNif 
                !kmult = kmult - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*n(k)*NNif 
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
                kmult = Dpn*Vpn(j)*n(k)*NNif
                !kmult = kmult - Gammaredpn*Dpn*Vveci(j)*n(k)*NNif
                IF (j == 5) THEN
                   kmult = kmult + Dnn*n(k)*NNif
                END IF 
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
                elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
          kmultf = dot_product(dDpn_dU,uf)*(gmpn(1)*n(1) + gmpn(2)*n(2))*Nif              
          !kmultf = kmultf - Gammaredpn*(Dpn*(dot_product(Taui(1,:),uf)*n(1) + dot_product(Taui(2,:),uf)*n(2)) + dot_product(dDpn_dU,uf)*(gmipn(1)*n(1) + gmipn(2)*n(2)))*Nif             
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
          elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf 
#endif                                                                           
       END IF
#endif
      END DO  ! i-Loop

      ! Assembly stabilization terms
      IF (numer%stab < 6) THEN
        DO i = 1,Neq
          ind_if = i + ind_asf
          kmult = tau(i,i)*NNif
          elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
          elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) - kmult
          elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) + kmult
        END DO
      ELSE
        DO i = 1,Neq
          ind_if = i + ind_asf
          DO j = 1,Neq
            ind_jf = j + ind_asf
            kmult = tau(i,j)*NNif
            elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
            elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) + kmult
          END DO
        END DO
      ENDIF
      !************* End stabilization terms************************

    END SUBROUTINE assemblyIntFacesContribution

    !********************************************************************
    !
    !         ASSEMBLY EXTERIOR FACES CONTRIBUTION
    !
    !********************************************************************
    SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,&
        &ind_fg,b3,Bmod,psi,n,diffiso,diffani,NNif,Nif,Nfbn,uf,upf,qf,tau,Vnng,ifa)
      integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      logical                   :: isdir
      real*8,intent(IN)         :: b3(:),n(:),Bmod, psi
      real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
      real*8,intent(IN)         :: NNif(:,:),Nif(:),Nfbn(:)
      real*8,intent(IN)         :: uf(:),upf(:)
      real*8,intent(IN)         :: qf(:)
      real*8,optional,intent(INOUT) :: tau(:,:),Vnng(:)
      integer*4,optional         :: ifa
      real*8                    :: kcoeff
      integer*4                 :: i,j,k,ii,alpha,beta
      integer*4,dimension(Npfl)  :: ind_if,ind_jf,ind_kf
      real*8,dimension(neq,neq) :: A
      real*8,dimension(neq,Ndim):: APinch
      real*8                    :: auxvec(neq)
      real*8                    :: bn,kmult(Npfl,Npfl),kmultf(Npfl)
      real*8                    :: Qpr(Ndim,Neq),exb(3)
      real*8                    :: nn(3),qq(3,Neq),b(Ndim),bb(3)
      real*8                    :: W2(Neq), dW2_dU(Neq,Neq), QdW2(Ndim,Neq)
#ifdef TEMPERATURE
      real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
      real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
      real*8                    :: W3(Neq), dW3_dU(Neq,Neq), QdW3(Ndim,Neq)
      real*8                    :: W4(Neq), dW4_dU(Neq,Neq), QdW4(Ndim,Neq)
#ifdef NEUTRALP
      real*8                    :: Dnn,Dpn,GammaLim,Alphanp,Betanp,Gammaredpn,Tmin
      real*8                    :: Anp(Neq),Vpn(Neq),dVpn_dU(Neq,Neq),gmpn(Ndim),gmipn(Ndim),Taupn(Ndim,Neq),dDpn_dU(Neq)
#endif
#endif

      b = b3(1:Ndim)
      bb = b3
      ! Jacobian matrices
      bn = dot_product(b,n)
      CALL jacobianMatrices(uf,A)

      ! Jacobian matrices Pinch
      CALL computePinch(b,psi,APinch)

      ! Compute Q^T^(k-1)
      Qpr = reshape(qf,(/Ndim,Neq/))

      ! Split diffusion vector/matrix for the momentum equation
      CALL compute_W2(uf,W2)
      CALL compute_dW2_dU(uf,dW2_dU)
      QdW2 = matmul(Qpr,dW2_dU)

      nn = 0.
      qq = 0.
      nn(1:Ndim) = n
      qq(1:Ndim,:) = Qpr

#ifdef TEMPERATURE

      ! Compute V(U^(k-1))
      call computeVi(uf,Vveci)
      call computeVe(uf,Vvece)

      ! Compute dV_dU (k-1)
      call compute_dV_dUi(uf,dV_dUi)
      call compute_dV_dUe(uf,dV_dUe)

      ! Split diffusion vector/matrix for the energies equations
      CALL compute_W3(uf,W3)
      CALL compute_dW3_dU(uf,dW3_dU)
      QdW3 = matmul(Qpr,dW3_dU)

      CALL compute_W4(uf,W4)
      CALL compute_dW4_dU(uf,dW4_dU)
      QdW4 = matmul(Qpr,dW4_dU)

      ! Compute Alpha(U^(k-1))
      Alphai = computeAlphai(uf)
      Alphae = computeAlphae(uf)

      ! Compute dAlpha/dU^(k-1)
      call compute_dAlpha_dUi(uf,dAlpha_dUi)
      call compute_dAlpha_dUe(uf,dAlpha_dUe)

      gmi = dot_product(matmul(Qpr,Vveci),b)  ! scalar
      gme = dot_product(matmul(Qpr,Vvece),b)
      Taui = matmul(Qpr,dV_dUi)               ! 2x3
      Taue = matmul(Qpr,dV_dUe)               ! 2x3
#ifdef NEUTRALP
      ! Compute Vpn(U^(k-1))
      CALL computeVpn(uf,Vpn)
	    gmpn = matmul(Qpr,Vpn)                       ! Ndim x 1
      ! Compute dVpn-dU(U^(k-1))
	    CALL compute_dVpn_dU(uf,dVpn_dU)
	    Taupn = matmul(Qpr,dVpn_dU)                 ! Ndim x Neq
	    ! Compute Dpn(U^(k-1))
      CALL computeDpn(uf,Qpr,Vpn,Dpn)
      ! Compute dDpn_dU(U^(k-1))
      CALL compute_dDpn_dU(uf,Qpr,Vpn,dDpn_dU)
      ! Reduce Grad Pn for low collision regime 
      ! Threshold set at 0.5xGradPn for Ti = 0.2 eV 
      Gammaredpn = 1.
      Tmin = 0.2/simpar%refval_temperature
      IF (Tmin/upf(7) .le. 1.) Gammaredpn = Gammaredpn*Tmin/upf(7)
      Dnn = Gammaredpn*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn 
      ! Comput Gammaredpn(U^(k-1))
      !CALL computeGammared(uf,Gammaredpn)
      !gmipn = matmul(Qpr,Vveci)
      !CALL computeGammaLim(ue,Qpr,Vpn,GammaLim)
      ! Set Grad Ti = 0. for low collision regime 
      ! (back to diffusion equation for neutral density)
      !CALL computeAlphaCoeff(uf,Qpr,Vpn,Alphanp)  
      !CALL computeBetaCoeff(uf,Qpr,Vpn,Betanp)  
      !Dpn = Alphanp*Dpn
      !dDpn_dU = Alphanp*dDpn_dU
      !Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*Dpn      
      !IF (Dnn .gt. phys%diff_nn) Dnn = phys%diff_nn    
      !IF (Dpn .gt. phys%diff_nn) THEN
      ! Dpn = Alphanp*Dpn !0.
      ! dDpn_dU = Alphanp*dDpn_dU !0.
      ! Dnn = Betanp*(simpar%refval_time**2/simpar%refval_length**2*simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass)*upf(7)*phys%diff_nn
      ! END IF  
      ! Set Gamma Convective = cs_n*n_n for low collision regime 
      !IF (Dpn .gt. phys%diff_nn) THEN
      !   Dpn = 0.
      !   dDpn_dU = 0.
      !   CALL jacobianMatricesNP(uf,Anp)
      !ELSE
      !   Anp = 0.
      !END IF   
#endif
#endif

      ! Assembly local matrix
      DO i = 1,Neq
        ind_if = ind_asf + i
        DO k = 1,Ndim
          ind_kf = ind_ash + k + (i - 1)*Ndim
          kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))
#ifdef VORTICITY
          if (i==4) then
            kmult=kmult/uf(1)
          endif
#endif
          ! Diagonal terms for perpendicular diffusion
          elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          IF(i == 2) THEN
            ! Assembly LQ
            j = 1 ! other terms are 0 anyway in vector W2
            ind_kf = ind_ash + k + (j - 1)*Ndim
             kmult = NNif*W2(j)*(n(k) - bn*b(k))
             elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          ! Assembly LU
             DO j=1,Neq
               ind_jf = ind_asf+j
               IF (.not. isdir) THEN
                kmult = QdW2(k,j)*NNif*(n(k)-bn*b(k))
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
               ENDIF
             END DO
#ifdef TEMPERATURE
         ELSEIF (i ==3) THEN
              ! Assembly LQ
              DO j=1,Neq ! here there are 2 non-zero elements in vector W3
               ind_kf = ind_ash + k + (j - 1)*Ndim
               kmult = NNif*W3(j)*(n(k) - bn*b(k))
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
             END DO
            ! Assembly LU
               DO j=1,Neq
                 ind_jf = ind_asf+j
                 IF (.not. isdir) THEN
                   kmult = QdW3(k,j)*NNif*(n(k)-bn*b(k))
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                 ENDIF
              END DO
          ELSEIF (i ==4) THEN
              ! Assembly LQ
              j = 1 !other terms are 0 anyway in vector W4
              ind_kf = ind_ash + k + (j - 1)*Ndim
              kmult = NNif*W4(j)*(n(k) - bn*b(k))
               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            ! Assembly LU
               DO j=1,Neq
                 ind_jf = ind_asf+j
                 IF (.not. isdir) THEN
                   kmult = QdW4(k,j)*NNif*(n(k)-bn*b(k))
                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                  ENDIF
                 !ind_jf = ind_asf + j
                 !IF (.not. isdir) THEN
                !    kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),b)))*NNif*bn
                !    elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                ! END IF
               END DO

#endif
END IF



#ifdef VORTICITY
          IF (switch%driftexb .and. i .ne. 4) THEN
            ! ExB terms
            ii = 4
            kcoeff = phys%dfcoef*numer%exbdump/Bmod
            call ijk_cross_product(k,alpha,beta)
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kmult = kcoeff*NNif*(nn(alpha)*b3(beta) - nn(beta)*b3(alpha))*uf(i)
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
            IF (.not. isdir) THEN
              call cross_product(qq(:,ii),bb,exb)
              kmult = kcoeff*exb(k)*NNif*b(k)
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
              kmultf = kcoeff*exb(k)*uf(i)*Nif*b(k)
              elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
            ENDIF
            tau(i,i) = tau(i,i) + abs(kcoeff*exb(k)*b(k))
            tau(i,i) = 100
          ENDIF
          DO ii = 1,Neq
            IF (ii == i) CYCLE ! diagonal alredy assembled
            IF (abs(diffiso(i,ii)) < 1e-12 .and. abs(diffani(i,ii)) < 1e-12) CYCLE
            ind_kf = ind_ash + k + (ii - 1)*Ndim
            kcoeff = 1.
            ! Non-linear correction for non-linear diffusive terms.
            ! TODO: find a smarter way to include it,avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)



            !               IF ((i == 3 .or. i == 4) .and. ii == 1) then
            IF ((i == 3 ) .and. ii == 1) then


              kcoeff = 1./uf(1)
              IF (.not. isdir) THEN
                ind_jf = ind_asf + ii
                kmult = kcoeff**2*(diffiso(i,ii)*Qpr(k,1)*n(k)*NNif - diffani(i,ii)*Qpr(k,1)*b(k)*NNif*bn)
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
                elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) + kcoeff*(diffiso(i,ii)*Qpr(k,1)*n(k)*Nif - &
                  &diffani(i,ii)*Qpr(k,1)*b(k)*Nfbn)
              END IF
            ENDIF
            kmult = NNif*kcoeff*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii))
            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
          END DO
#endif
        END DO ! k-loop

        ! Convection contribution
        IF (.not. isdir) THEN
          DO j = 1,Neq
            ind_jf = ind_asf + j
            kmult = bn*A(i,j)*NNif
!#ifdef NEUTRALP
!            IF (i == 5) kmult = kmult + bn*Anp(j)*NNif
!#endif
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
!#ifdef TEMPERATURE
!#ifdef NEUTRAL
!            !X component neutral convective velocity
!            kmult = n(1)*Ax(i,j)*NNif
!            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!            !Y component neutral convective velocity
!            kmult = n(2)*Ay(i,j)*NNif
!            elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
!#endif
!#endif
          END DO ! j-loop
          ! Pinch Contribution
          elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) + (APinch(i,1)*n(1) + APinch(i,2)*n(2))*NNif
        ENDIF


#ifndef TEMPERATURE
        ! Added term for n=exp(x) change of variable
        if (.not. isdir) then
          if (switch%logrho) then
            call logrhojacobianVector(uf,upf,auxvec)
            kmultf = Nfbn*auxvec(i)
            elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel)-kmultf
          endif
        endif
#endif



#ifdef TEMPERATURE
        ! Parallel diffusion for the temperature
        IF (i == 3) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
            IF (.not. isdir) THEN
              kmult = coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),b)))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
          kmultf = coefi*Alphai*(dot_product(matmul(transpose(Taui),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
        ELSEIF (i == 4) THEN
          DO j = 1,4
            ind_jf = ind_asf + j
            IF (.not. isdir) THEN
              kmult = coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),b)))*NNif*bn
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            END IF
            DO k = 1,Ndim
              ind_kf = k + (j - 1)*Ndim + ind_ash
              kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
              elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
            END DO
          END DO
          kmultf = coefe*Alphae*(dot_product(matmul(transpose(Taue),b),uf))*Nfbn
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#ifdef NEUTRALP
       ELSEIF (i == 5) THEN
          DO j = 1,Neq
             ind_jf = ind_asf + j
             DO k = 1,Ndim
                ind_kf = k + (j - 1)*Ndim + ind_ash
                kmult = (Dpn*Taupn(k,j) + dDpn_dU(j)*gmpn(k))*n(k)*NNif 
                !kmult = kmult - Gammaredpn*(Dpn*Taui(k,j) + dDpn_dU(j)*gmipn(k))*n(k)*NNif 
                elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                kmult = Dpn*Vpn(j)*n(k)*NNif
                !kmult = kmult - Gammaredpn*Dpn*Vveci(j)*n(k)*NNif
                IF (j == 5) THEN
                   kmult = kmult + Dnn*n(k)*NNif
                END IF
                elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) - kmult
             END DO
          END DO
          kmultf = dot_product(dDpn_dU,uf)*(gmpn(1)*n(1) + gmpn(2)*n(2))*Nif        
          !kmultf = kmultf - Gammaredpn*(Dpn*(dot_product(Taui(1,:),uf)*n(1) + dot_product(Taui(2,:),uf)*n(2)) + dot_product(dDpn_dU,uf)*(gmipn(1)*n(1) + gmipn(2)*n(2)))*Nif                             
          elMat%S(ind_fe(ind_if),iel) = elMat%S(ind_fe(ind_if),iel) - kmultf
#endif
      END IF
#endif
      END DO  ! i-Loop

      ! Assembly stabilization terms
      IF (numer%stab < 6) THEN
        DO i = 1,Neq
          ind_if = i + ind_asf
          kmult = tau(i,i)*NNif
          elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
          IF (.not. isdir) THEN
            elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
          ENDIF
        END DO
      ELSE
        DO i = 1,Neq
          ind_if = i + ind_asf
          DO j = 1,Neq
            ind_jf = j + ind_asf
            kmult = tau(i,j)*NNif
            elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
            IF (.not. isdir) THEN
              elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
            ENDIF
          END DO
        END DO
      ENDIF
      !************* End stabilization terms************************
    END SUBROUTINE assemblyExtFacesContribution

    SUBROUTINE do_assembly(Auq,Auu,rhs,ind_ass,ind_asq,iel)
      real*8,intent(in)    :: Auq(:,:,:),Auu(:,:,:),rhs(:,:)
      integer*4,intent(in) :: iel,ind_ass(:),ind_asq(:)
      integer :: i,j,k,z
      integer*4,dimension(Npel) :: ind_i,ind_j,ind_k

      DO i = 1,Neq
        ind_i = i + ind_ass
        elMat%S(ind_i,iel)=elMat%S(ind_i,iel)+rhs(:,i)
        DO j = 1,Neq
          ind_j = j + ind_ass
          z = i+(j-1)*Neq
          elMat%Auu(ind_i,ind_j,iel)=elMat%Auu(ind_i,ind_j,iel)+Auu(:,:,z)
          DO k = 1,Ndim
            z = i+(k-1)*Neq+(j-1)*Ndim*Neq
            ind_k = ind_asq + k + (j - 1)*Ndim
            elMat%Auq(ind_i,ind_k,iel)=elMat%Auq(ind_i,ind_k,iel)+Auq(:,:,z)
          END DO
        END DO
      END DO

    END SUBROUTINE do_assembly

#ifdef NEUTRAL
  !********************************************************************
  !
  !         ASSEMBLY NEUTRAL SOURCE MATRIX
  !
  !********************************************************************
#ifdef TEMPERATURE
  SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,sigmaviz,dsigmaviz_dU,sigmavrec,dsigmavrec_dU,&
      &fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,sigmavcx,dsigmavcx_dU,fEiiz,&
      &dfEiiz_dU,fEirec,dfEirec_dU,fEicx,dfEicx_dU,Tloss,dTloss_dU,Tlossrec,dTlossrec_dU,Sn,Sn0)
#else
    SUBROUTINE assemblyNeutral(U,niz,dniz_dU,nrec,dnrec_dU,fGammacx,dfGammacx_dU,fGammarec,dfGammarec_dU,Sn,Sn0)
#endif
      real*8, intent(IN) :: niz,nrec,fGammacx,fGammarec
      real*8, intent(IN) :: U(:),dniz_dU(:),dnrec_dU(:),dfGammacx_dU(:),dfGammarec_dU(:)
#ifndef TEMPERATURE
      real*8             :: sigmaviz,sigmavrec,sigmavcx
#else
      real*8, intent(IN) :: sigmaviz,sigmavrec,sigmavcx,fEiiz,fEirec,fEicx,Tloss,Tlossrec
      real*8, intent(IN) :: dsigmaviz_dU(:),dsigmavrec_dU(:),dsigmavcx_dU(:),dTloss_dU(:),dTlossrec_dU(:)
      real*8, intent(IN) :: dfEiiz_dU(:),dfEirec_dU(:),dfEicx_dU(:)
#endif
      real*8             :: ad,ad4,RE,Sn(:,:),Sn0(:),Ti,Te

      Sn   = 0.
      Sn0  = 0.
      !RE   = 0.2
      RE   = 1.
      ad   = 1.3737e12
      ad4  = (ad*1.6e-19)/((1.3839e4**2)*3.35e-27)

#ifndef TEMPERATURE
      sigmaviz   = 3.01e-14
      sigmavrec  = 1.3638e-20
      sigmavcx   = 4.0808e-15
#endif


      !Assembly Source Terms in plasma density equation
      Sn(1,1)   = ad*(-dniz_dU(1)*sigmaviz + dnrec_dU(1)*sigmavrec)
      Sn(1,Neq) = ad*(-dniz_dU(Neq)*sigmaviz)
#ifdef TEMPERATURE
      Sn(1,1)   = Sn(1,1) + ad*(-niz*dsigmaviz_dU(1) + nrec*dsigmavrec_dU(1))
      Sn(1,4)   = Sn(1,4) + ad*(-niz*dsigmaviz_dU(4) + nrec*dsigmavrec_dU(4))
#endif
      !Assembly Source Terms in plasma momentum equation
      Sn(2,1)   = ad*(dfGammarec_dU(1)*sigmavrec)
      Sn(2,2)   = ad*(dfGammacx_dU(2)*sigmavcx + dfGammarec_dU(2)*sigmavrec)
      Sn(2,Neq) = ad*(dfGammacx_dU(Neq)*sigmavcx)
#ifdef TEMPERATURE
      Sn(2,1)   = Sn(2,1) + ad*( fGammacx*dsigmavcx_dU(1) + fGammarec*dsigmavrec_dU(1))
      Sn(2,4)   = Sn(2,4) + ad*( fGammacx*dsigmavcx_dU(4) + fGammarec*dsigmavrec_dU(4))
      !Assembly Source Terms in ion energy equation
      Sn(3,1)   = ad*( - RE*fEiiz*dsigmaviz_dU(1) + dfEirec_dU(1)*sigmavrec + fEirec*dsigmavrec_dU(1) +&
        &dfEicx_dU(1)*sigmavcx + fEicx*dsigmavcx_dU(1))
      Sn(3,2)   = ad*(dfEicx_dU(2)*sigmavcx )
      Sn(3,3)   = ad*(-RE*dfEiiz_dU(3)*sigmaviz + dfEirec_dU(3)*sigmavrec)
      Sn(3,4)   = ad*(-RE*fEiiz*dsigmaviz_dU(4) + fEirec*dsigmavrec_dU(4) + fEicx*dsigmavcx_dU(4))
      Sn(3,Neq) = ad*(-RE*dfEiiz_dU(Neq)*sigmaviz + dfEicx_dU(Neq)*sigmavcx)
      !Assembly Source Terms in electron energy equation
      Sn(4,1)   = ad4*(dniz_dU(1)*sigmaviz*Tloss + niz*dsigmaviz_dU(1)*Tloss + niz*sigmaviz*dTloss_dU(1) +&
        &dnrec_dU(1)*sigmavrec*Tlossrec + nrec*dsigmavrec_dU(1)*Tlossrec + nrec*sigmavrec*dTlossrec_dU(1))
      Sn(4,4)   = ad4*(niz*dsigmaviz_dU(4)*Tloss + niz*sigmaviz*dTloss_dU(4) +&
        &nrec*dsigmavrec_dU(4)*Tlossrec + nrec*sigmavrec*dTlossrec_dU(4))
      Sn(4,Neq) = ad4*(dniz_dU(Neq)*sigmaviz*Tloss )
#endif
      !Assembly Source Terms in neutral density equation
      Sn(Neq,:) = -Sn(1,:)

      !Assembly RHS Neutral Source Terms
      Sn0(1)    = ad*(niz*sigmaviz - nrec*sigmavrec)
      Sn0(2)    = ad*(-fGammacx*sigmavcx - fGammarec*sigmavrec)
#ifdef TEMPERATURE
      Sn0(3)    = ad*(RE*fEiiz*sigmaviz - fEirec*sigmavrec - fEicx*sigmavcx)
      Sn0(4)    = ad4*(-niz*sigmaviz*Tloss - nrec*sigmavrec*Tlossrec)
#endif
      Sn0(Neq)  = -Sn0(1)

      !Thresholds: 
#ifdef TEMPERATURE
      Ti = 2./(3.*phys%Mref)*(U(3)/U(1) - 1./2.*(U(2)/U(1))**2)
      Te = 2./(3.*phys%Mref)*U(4)/U(1)
      if (Ti .le. 6.e-4 .or. Te .le. 6.e-4) then 
!         Sn(1,:) = - abs(Sn(1,:))
!         Sn0(1) = - abs(Sn0(1))
!         Sn(2,:) = - abs(Sn(2,:))
!         Sn0(2) = - abs(Sn0(2))
         Sn(3,1) = ad*( - RE*fEiiz*dsigmaviz_dU(1))
         Sn(3,2) = 0.
         Sn(3,3) = ad*(-RE*dfEiiz_dU(3)*sigmaviz)
         Sn(3,4) = ad*(-RE*fEiiz*dsigmaviz_dU(4))
         Sn(3,5) = ad*(-RE*dfEiiz_dU(Neq)*sigmaviz)
         Sn0(3) = ad*(RE*fEiiz*sigmaviz)
         Sn(4,1) = 0. !3./2.*6.e-10
         Sn(4,2) = 0.
         Sn(4,3) = 0.
         Sn(4,4) = 0. !-1. 
         Sn(4,5) = 0.
         Sn0(4) = 0.
!         Sn(5,:) = - abs(Sn(5,:))
!         Sn0(5) = - abs(Sn0(5))
!         Sn(:,:) = 0.
!         Sn0(:) = 0.
      endif
#endif

    END SUBROUTINE assemblyNeutral
#endif


  END SUBROUTINE HDG_computeJacobian
