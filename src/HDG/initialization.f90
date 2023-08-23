!************************************************************
! project: MHDG
! file: initialization.f90
! date: 03/01/2017
! Initialization of the physics depending on the model,of
! the elemental matrices,of the solution
!************************************************************
MODULE initialization
  USE globals
  USE printutils
  USE analytical
  USE physics
  USE MPI_OMP
  USE LinearAlgebra, ONLY: col,tensorProduct,solve_linear_system


  IMPLICIT NONE
CONTAINS

  !*********************************************
  ! Initialization of the simulation parameters
  !*********************************************
  SUBROUTINE init_sim(nts,dt)

    integer,intent(out) :: nts
    real,intent(out) :: dt

    CALL initPhys()

    ! the time is set to zero here but in case of a restart this value is
    ! overwritten by the loaded solution final time
    time%t = 0.

    ! time step is initialized to dt0.
    dt = time%dt0
    time%dt = dt

    ! time count
    time%it = 0
    time%ik = 0

    ! number of time step of the simulation: if it is a steady state simulation
    ! ndt is set to 1,otherwise the input value is used
    IF (switch%steady) THEN
      nts = 1
    ELSE
    nts = time%nts
    END IF

    ! Allocate and initialize time residual
    ALLOCATE (sol%tres(nts))
    ALLOCATE (sol%time(nts))
    sol%tres = 0.
    sol%time = 0.
    sol%Nt = 0
  END SUBROUTINE init_sim

  !*******************************************
  ! Initialization of the elemental matrices
  !*******************************************
  SUBROUTINE init_elmat
    integer :: Neq,Ndim,Ntorloc,N2d,Nel,Np1d,Np2d,Np,Nfl,Nfp,Nfg,Nf

    Neq = phys%Neq                      ! N. of equations
#ifdef TOR3D
    Ndim = 3                             ! N. of dimensions
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
      ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                   ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
#endif

    ALLOCATE (elmat%iAqq(Neq*Ndim*Np,Neq*Ndim*Np,Nel))
    ALLOCATE (elmat%Aqu(Neq*Ndim*Np,Neq*Np,Nel))
    ALLOCATE (elmat%Aql(Neq*Ndim*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%Auq(Neq*Np,Ndim*Neq*Np,Nel))
    ALLOCATE (elmat%Auu(Neq*Np,Neq*Np,Nel))
    ALLOCATE (elmat%Aul(Neq*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%Alq(Neq*Nfg,Neq*Ndim*Np,Nel))
    ALLOCATE (elmat%Alu(Neq*Nfg,Neq*Np,Nel))
    ALLOCATE (elmat%All(Neq*Nfg,Neq*Nfg,Nel))
    ALLOCATE (elmat%Aql_dir(Neq*Np*Ndim,Nel))
    ALLOCATE (elmat%Aul_dir(Neq*Np,Nel))
    ALLOCATE (elmat%S(Neq*Np,Nel))
    ALLOCATE (elmat%fH(Neq*Nfg,Nel))
    IF (switch%ME) THEN
       ALLOCATE(phys%puff_exp(time%nts))
    END IF
    IF (switch%saveTau) THEN
       ALLOCATE(phys%diff_nn_Vol(Mesh%Nelems*refElPol%NGauss2D))
       ALLOCATE(phys%diff_nn_Fac(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D))
       ALLOCATE(phys%diff_nn_Bou(Mesh%Nextfaces*refElPol%NGauss1D))
       ALLOCATE(phys%v_nn_Vol(Mesh%Nelems*refElPol%NGauss2D,Mesh%Ndim))
       ALLOCATE(phys%v_nn_Fac(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(phys%v_nn_Bou(Mesh%Nextfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(Mesh%Xg(Mesh%Nelems*refElPol%NGauss2D,Mesh%Ndim))
       ALLOCATE(Mesh%Xgf(Mesh%Nelems*refElPol%Nfaces*refElPol%NGauss1D,Mesh%Ndim))
       ALLOCATE(Mesh%Xgb(Mesh%Nextfaces*refElPol%NGauss1D,Mesh%Ndim))
    END IF

    elmat%iAqq = 0.
    elmat%Aqu = 0.
    elmat%Aql = 0.
    elmat%Auq = 0.
    elmat%Auu = 0.
    elmat%Aul = 0.
    elmat%Alq = 0.
    elmat%Alu = 0.
    elmat%All = 0.
    elmat%Aql_dir = 0.
    elmat%Aul_dir = 0.
    elmat%S = 0.
    elmat%fh = 0.
    ALLOCATE (elmat%UU(Neq*Np,Neq*Nfg,Nel))
    ALLOCATE (elmat%U0(Neq*Np,Nel))
    ALLOCATE (elmat%LL(Neq*Np*Ndim,Neq*Nfg,Nel))
    ALLOCATE (elmat%L0(Neq*Np*Ndim,Nel))
    elmat%UU = 0.d0
    elmat%U0 = 0.d0
    elmat%LL = 0.d0
    elmat%L0 = 0.d0
    IF (switch%ME) THEN
       phys%puff_exp = 0.
    END IF
    IF (switch%saveTau) THEN
       phys%diff_nn_Vol = 0.
       phys%diff_nn_Fac = 0.
       phys%diff_nn_Bou = 0.
       phys%v_nn_Vol = 0.
       phys%v_nn_Fac = 0.
       phys%v_nn_Bou = 0.
       Mesh%Xg = 0.
       Mesh%Xgf = 0.
       Mesh%Xgb = 0.
    END IF
	
  END SUBROUTINE init_elmat

  subroutine init_solve_timing
    use globals
    timing%cputpre=1.e-8
    timing%cputmap=1.e-8
    timing%cputass=1.e-8
    timing%cputbcd=1.e-8
    timing%cputsol=1.e-8
    timing%cputjac=1.e-8
    timing%cputglb=1.e-8
    timing%cputcom=1.e-8

    timing%runtpre=1.e-8
    timing%runtmap=1.e-8
    timing%runtass=1.e-8
    timing%runtbcd=1.e-8
    timing%runtsol=1.e-8
    timing%runtjac=1.e-8
    timing%runtglb=1.e-8
    timing%runtcom=1.e-8

    timing%clstime1=1.e-8
    timing%clstime2=1.e-8
    timing%clstime3=1.e-8
    timing%clstime4=1.e-8
    timing%clstime5=1.e-8
    timing%clstime6=1.e-8

    timing%rlstime1=1.e-8
    timing%rlstime2=1.e-8
    timing%rlstime3=1.e-8
    timing%rlstime4=1.e-8
    timing%rlstime5=1.e-8
    timing%rlstime6=1.e-8
  end subroutine

  !********************************
  ! Initialization of the solution
  !********************************
  SUBROUTINE init_sol
    integer :: Neq,Ndim,Ntorloc,N2d,Nel,Np1d,Np2d,Np,Nfl,Nfp,Nfg,Nf,sizeutilde,sizeu
    integer :: Ngvo
    Neq = phys%Neq
#ifdef TOR3D
    Ndim = 3                             ! N. of dimensions
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
      ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                   ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
    sizeu = Neq*Nel*Np                    ! Size of u
    Ngvo = refElPol%Ngauss2d*refEltor%Ngauss1d
#ifdef PARALL
    if (MPIvar%ntor .gt. 1) then
      sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d) + Neq*Np2d*N2d! Size of utilde
    else
      sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
    endif
#else
    sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
#endif
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
    sizeu = Neq*Nel*Np
    sizeutilde = Neq*Mesh%Nfaces*Mesh%Nnodesperface
    Ngvo = refElPol%Ngauss2d
#endif

    ! Allocation of the solution vector
    ALLOCATE (sol%u(sizeu))
    ALLOCATE (sol%u_tilde(sizeutilde))
    ALLOCATE (sol%q(sizeu*Ndim))
    sol%u = 0.
    sol%u_tilde = 0.
    sol%q = 0.
    ! Initialize the solution
    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        write (6,*) "*** Initializing the solution"
      END IF
    ENDIF
    if (switch%init.eq.1) then
      ! The solution is intialized in each node to the analytical solution
      IF (MPIvar%glob_id .eq. 0) THEN
        IF (utils%printint > 0) THEN
          write (6,*) "******* Initializing the solution to the analytic solution"
        END IF
      ENDIF
      CALL init_sol_analytic()
    elseif (switch%init.eq.2) then
      ! The solution is intialized in each node to the analytical solution
      IF (MPIvar%glob_id .eq. 0) THEN
        IF (utils%printint > 0) THEN
          write (6,*) "******* Initializing the solution with L2 projection"
        END IF
      ENDIF
      CALL init_sol_l2proj()
    else
      write(6,*) "Wrong initialization type"
      stop
    endif
    ! Extract the face solution from the elemental one
    CALL extractFaceSolution()

    IF (MPIvar%glob_id .eq. 0) THEN
      IF (utils%printint > 0) THEN
        write (6,*) "Done!"
      END IF
    ENDIF
  CONTAINS
    !***********************************************************
    ! Initialization of the solution using the analytic solution
    !***********************************************************
    SUBROUTINE init_sol_analytic
      integer             :: itor,itorg,iel,iel3,i
      integer             :: ind(Np)
      real*8              :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
      real*8              :: ue(Np,Neq)
      real*8,allocatable :: u(:,:)
      real*8              :: tdiv(numer%ntor + 1)
#ifdef TOR3D
      real*8              :: htor,tel(refElTor%Nnodes1d)
#endif
      real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)
      real*8              :: uex(Np,Neq),uey(Np,Neq)
#ifdef TOR3D
      real*8,allocatable :: qt(:,:)
      real*8              :: uet(Np,Neq)
#endif
      ALLOCATE (u(Nel*Np,phys%Neq))
      u = 0.
      ALLOCATE (qx(Nel*Np,phys%Neq))
      ALLOCATE (qy(Nel*Np,phys%Neq))
      ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
      qx = 0.; qy = 0.
#ifdef TOR3D
      ALLOCATE (qt(Nel*Np,phys%Neq))
      qt = 0.
#endif

#ifdef TOR3D

      !****************************************
      !          3D
      !****************************************
      htor = numer%tmax/numer%ntor
      tdiv = 0.
      DO i = 1,numer%ntor
        tdiv(i + 1) = i*htor
      END DO

      DO itor = 1,ntorloc
#ifdef PARALL
        itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
        if (itorg == numer%ntor + 1) itorg = 1
#else
        itorg = itor
#endif
        tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
        DO iel = 1,Mesh%Nelems
          iel3 = (itor - 1)*N2d+iel
          ind = (iel3 - 1)*Np + (/(i,i=1,Np)/)
          Xe = Mesh%X(Mesh%T(iel,:),:)
          CALL analytical_solution(iel,Xe(:,1),Xe(:,2),tel,ue)
          CALL analytical_gradient(Xe(:,1),Xe(:,2),tel,ue,uex,uey,uet)
          qx(ind,:) = uex
          qy(ind,:) = uey
          qt(ind,:) = uet
          u(ind,:) = ue
        END DO
      END DO
#else
      !****************************************
      !          2D
      !****************************************
      DO iel = 1,Mesh%Nelems

        ind = (iel - 1)*Np + (/(i,i=1,Np)/)
        Xe = Mesh%X(Mesh%T(iel,:),:)
        CALL analytical_solution(iel,Xe(:,1),Xe(:,2),ue)
        CALL analytical_gradient(Xe(:,1),Xe(:,2),ue,uex,uey)
        qx(ind,:) = uex
        qy(ind,:) = uey
        u(ind,:) = ue
      END DO
#endif

      !****************************************
      !          common
      !****************************************
      sol%u = reshape(transpose(u),(/Nel*Np*phys%Neq/))

      auxq(:,1) = reshape(transpose(qx),(/Nel*Np*phys%Neq/))
      auxq(:,2) = reshape(transpose(qy),(/Nel*Np*phys%Neq/))
#ifdef TOR3D
      auxq(:,3) = reshape(transpose(qt),(/Nel*Np*phys%Neq/))
#endif
      sol%q = reshape(transpose(auxq),(/Nel*Np*phys%Neq*Ndim/))
      DEALLOCATE (qx,qy,auxq)
#ifdef TOR3D
      DEALLOCATE (qt)
#endif
      DEALLOCATE (u)
    END SUBROUTINE init_sol_analytic







    !***********************************************************
    ! Initialization of the solution using an L2 projection
    !***********************************************************
    SUBROUTINE init_sol_l2proj
      integer             :: itor,itorg,iel,iel3,i,g
      integer             :: ind(Np)
      real*8              :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
      real*8              :: ue(Np,Neq)
      real*8              :: uex(Np,Neq),uey(Np,Neq)
      real*8,allocatable  :: u(:,:)
      real*8              :: tdiv(numer%ntor + 1)
#ifdef TOR3D
      real*8              :: htor,tel(refElTor%Nnodes1d)
#endif
      real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)

#ifdef TOR3D
      real*8,allocatable :: qt(:,:)
      real*8              :: uet(Np,Neq)
#endif
      real*8                :: J11(Ngvo),J12(Ngvo)
      real*8                :: J21(Ngvo),J22(Ngvo)
      real*8                :: detJ(Ngvo)
      real*8                :: iJ11(Ngvo),iJ12(Ngvo)
      real*8                :: iJ21(Ngvo),iJ22(Ngvo)
      real*8                :: dvolu,M(Np,Np),rhs_u(Np,Neq),rhs_ux(Np,Neq),rhs_uy(Np,Neq)
      real*8                :: xyg(Ngvo,2)
      real*8                :: ug(Ngvo,Neq)
      real*8                :: ugx(Ngvo,Neq),ugy(Ngvo,Neq)
      real*8,parameter     :: tol = 1e-12

      ALLOCATE (u(Nel*Np,phys%Neq))
      u = 0.
      ALLOCATE (qx(Nel*Np,phys%Neq))
      ALLOCATE (qy(Nel*Np,phys%Neq))
      ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
      qx = 0.; qy = 0.
#ifdef TOR3D
      ALLOCATE (qt(Nel*Np,phys%Neq))
      qt = 0.
#endif

#ifdef TOR3D
      write(6,*) "Not coded yet"
      stop
      !         !****************************************
      !         !          3D
      !         !****************************************
      !         htor = numer%tmax/numer%ntor
      !         tdiv = 0.
      !         DO i = 1,numer%ntor
      !            tdiv(i + 1) = i*htor
      !         END DO

      !         DO itor = 1,ntorloc
      !#ifdef PARALL
      !            itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      !            if (itorg == numer%ntor + 1) itorg = 1
      !#else
      !            itorg = itor
      !#endif
      !            tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
      !            DO iel = 1,Mesh%Nelems
      !               iel3 = (itor - 1)*N2d+iel
      !               ind = (iel3 - 1)*Np + (/(i,i=1,Np)/)
      !               Xe = Mesh%X(Mesh%T(iel,:),:)
      !               CALL analytical_solution(Xe(:,1),Xe(:,2),tel,ue)
      !               CALL analytical_gradient(Xe(:,1),Xe(:,2),tel,ue,uex,uey,uet)
      !               qx(ind,:) = uex
      !               qy(ind,:) = uey
      !               qt(ind,:) = uet
      !               u(ind,:) = ue
      !            END DO
      !         END DO
#else
      !****************************************
      !          2D
      !****************************************
      DO iel = 1,Mesh%Nelems

        ind = (iel - 1)*Np + (/(i,i=1,Np)/)
        Xe = Mesh%X(Mesh%T(iel,:),:)
        J11 = matmul(refElPol%Nxi2D,Xe(:,1))                           ! ng x 1
        J12 = matmul(refElPol%Nxi2D,Xe(:,2))                           ! ng x 1
        J21 = matmul(refElPol%Neta2D,Xe(:,1))                          ! ng x 1
        J22 = matmul(refElPol%Neta2D,Xe(:,2))                          ! ng x 1
        detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
        iJ11 = J22/detJ
        iJ12 = -J12/detJ
        iJ21 = -J21/detJ
        iJ22 = J11/detJ

        ! Solution at Gauss points
        xyg = matmul(refElPol%N2D,Xe)
        CALL analytical_solution(iel,xyg(:,1),xyg(:,2),ug)
        CALL analytical_gradient(xyg(:,1),xyg(:,2),ug,ugx,ugy)


        ! Initialize mass matrix and rhs
        M=0.
        rhs_u=0.
        rhs_ux=0.
        rhs_uy=0.
        DO g=1,Ngvo


          IF (detJ(g) < tol) THEN
            error stop "Negative jacobian"
          END if

          ! Integration weight
          dvolu = refElPol%gauss_weights2D(g)*detJ(g)
          IF (switch%axisym) THEN
            dvolu = dvolu*xyg(g,1)
          END IF


          M = M + TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu
          rhs_u = rhs_u+tensorProduct(refElPol%N2D(g,:),ug(g,:))*dvolu
          rhs_ux = rhs_ux+tensorProduct(refElPol%N2D(g,:),ugx(g,:))*dvolu
          rhs_uy = rhs_uy+tensorProduct(refElPol%N2D(g,:),ugy(g,:))*dvolu
        END DO
        call solve_linear_system(M,rhs_u,ue)
        call solve_linear_system(M,rhs_ux,uex)
        call solve_linear_system(M,rhs_uy,uey)
        qx(ind,:) = uex
        qy(ind,:) = uey
        u(ind,:) = ue
      END DO




#endif

      !****************************************
      !          common
      !****************************************
      sol%u = reshape(transpose(u),(/Nel*Np*phys%Neq/))

      auxq(:,1) = reshape(transpose(qx),(/Nel*Np*phys%Neq/))
      auxq(:,2) = reshape(transpose(qy),(/Nel*Np*phys%Neq/))
#ifdef TOR3D
      auxq(:,3) = reshape(transpose(qt),(/Nel*Np*phys%Neq/))
#endif
      sol%q = reshape(transpose(auxq),(/Nel*Np*phys%Neq*Ndim/))
      DEALLOCATE (qx,qy,auxq)
#ifdef TOR3D
      DEALLOCATE (qt)
#endif
      DEALLOCATE (u)
    END SUBROUTINE init_sol_l2proj


  END SUBROUTINE init_sol



    !     subroutine reset_variables()
    !         integer      :: iel,i,iface,ifa,ieq
    !         integer      :: Np,Neq,Nf,Nel,Nut
    !         integer      :: ind(Mesh%Nnodesperelem)
    !         real*8       :: Xe(Mesh%Nnodesperelem,Mesh%Ndim),Xf(Mesh%Nnodesperface,Mesh%Ndim)
    !         real*8       :: ue(refEl%Nnodes2D,phys%Neq)
    !         integer      :: ind_uf(Mesh%Nnodesperface),faceNodes(Mesh%Nnodesperface)
    !         real*8,allocatable :: u(:,:)
    !         real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)
    !         real*8,allocatable :: u_tilde(:,:)




    !!         integer :: neq,Ne,Nf,Nfe,unkF,Np,Nfp,iElem,ifa,iFace,i
    !!         integer :: ind_uf(1:Mesh%Nnodesperface),faceNodes(1:Mesh%Nnodesperface)
    !!         integer :: ind_ue(1:Mesh%Nnodesperelem)
    !!         logical :: alreadydone(1:Mesh%ukf)
    !
    !!         real*8,allocatable :: u(:,:)
    !!         real*8              :: tdiv(numer%ntor + 1)
    !!#ifdef TOR3D
    !!         real*8              :: htor,tel(refElTor%Nnodes1d)
    !!#endif
    !!         real*8,allocatable :: qx(:,:),qy(:,:),auxq(:,:)
    !!         real*8              :: uex(Np,Neq),uey(Np,Neq)
    !!#ifdef TOR3D
    !!         real*8,allocatable :: qt(:,:)
    !!         real*8              :: uet(Np,Neq)
    !!#endif
    !!         ALLOCATE (u(Nel*Np,phys%Neq))
    !!         u = 0.
    !!         ALLOCATE (qx(Nel*Np,phys%Neq))
    !!         ALLOCATE (qy(Nel*Np,phys%Neq))
    !!         ALLOCATE (auxq(Nel*Np*phys%Neq,Ndim))
    !!         qx = 0.; qy = 0.
    !!#ifdef TOR3D
    !!         ALLOCATE (qt(Nel*Np,phys%Neq))
    !!         qt = 0.
    !!#endif
    !         neq = phys%Neq
    !         Nel = Mesh%Nelems
    !         Np = Mesh%Nnodesperelem
    !         Nfp = Mesh%Nnodesperface
    !         Nut = size(sol%u_tilde)
    !
    !         ALLOCATE (u(Nel*Np,Neq))
    !         ALLOCATE (qx(Nel*Np,Neq))
    !         ALLOCATE (qy(Nel*Np,Neq))
    !         ALLOCATE (auxq(Nel*Np*Neq,Ndim))
    !         ALLOCATE (u_tilde(Nut,neq))
    !
    !         u = transpose(reshape(sol%u,[Neq,Nel*Np]))
    !         u_tilde = transpose(reshape(sol%u_tilde,[Neq,Nut]))
    !         q = transpose(reshape(sol%u,[Neq*Ndim,Nel*Np]))
    !
    !
    !
    !#ifdef TOR3D
    ! write(6,*) "Reset variables error in 3D: Not coded yet"
    ! stop
    !#else
    !         !****************************************
    !         !          2D
    !         !****************************************
    !         DO iel = 1,Mesh%Nelems

    !            ind = (iel - 1)*Np + (/(i,i=1,Np)/)
    !            Xe = Mesh%X(Mesh%T(iel,:),:)
    !            CALL analytical_solution(Xe(:,1),Xe(:,2),ue)
    !            CALL analytical_gradient(Xe(:,1),Xe(:,2),ue,uex,uey)
    !            do ieq = 1,phys%neq
    !               if (switch%reset_eqs(ieq).ne.0) then
    !                  qx(ind,ieq) = uex
    !                  qy(ind,ieq) = uey
    !                  u(ind,ieq) = ue
    !
    !               endif
    !            end do
    !         END DO
    !
    !
    !         DO iFace = 1,Mesh%Nintfaces
    !            iel = Mesh%intfaces(iFace,1)
    !            ifa = Mesh%intfaces(iFace,2)
    !            faceNodes = refElPol%Face_nodes(ifa,:)
    !            Xf = Mesh%X(Mesh%T(iel,faceNodes),:)
    !            CALL analytical_solution(Xf(:,1),Xf(:,2),uf)
    !            ind_uf = (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
    !            do ieq = 1,phys%neq
    !               if (switch%reset_eqs(ieq).ne.0) then
    !                  u_tilde(ind_uf,ieq) = uf
    !               endif
    !            end do
    !         END DO

    !         DO iFace = 1,Mesh%Nextfaces
    !            iel = Mesh%extfaces(iFace,1)
    !            ifa = Mesh%extfaces(iFace,2)
    !            faceNodes = refElPol%Face_nodes(ifa,:)
    !            Xf = Mesh%X(Mesh%T(iel,faceNodes),:)
    !            CALL analytical_solution(Xf(:,1),Xf(:,2),uf)
    !            IF (.not. Mesh%Fdir(iel,ifa)) THEN
    !               ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
    !               do ieq = 1,phys%neq
    !                  if (switch%reset_eqs(ieq).ne.0) then
    !                     u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
    !                  endif
    !               end do
    !            END IF
    !         END DO
    !
    !#endif

    !!         !****************************************
    !!         !          common
    !!         !****************************************
    !!         sol%u = reshape(transpose(u),(/Nel*Np*phys%Neq/))

    !!         auxq(:,1) = reshape(transpose(qx),(/Nel*Np*phys%Neq/))
    !!         auxq(:,2) = reshape(transpose(qy),(/Nel*Np*phys%Neq/))
    !!#ifdef TOR3D
    !!         auxq(:,3) = reshape(transpose(qt),(/Nel*Np*phys%Neq/))
    !!#endif
    !!         sol%q = reshape(transpose(auxq),(/Nel*Np*phys%Neq*Ndim/))
    !!         DEALLOCATE (qx,qy,auxq)
    !!#ifdef TOR3D
    !!         DEALLOCATE (qt)
    !!#endif
    !!         DEALLOCATE (u)
    !
    !     end subroutine reset_variables








    !***************************************************************
    ! Extract face solution: routine to define a nodal face solution
    ! equal to the elemental solution at face nodes
    !***************************************************************
#ifdef TOR3D
    SUBROUTINE extractFaceSolution
      integer :: neq,Ne,Nf,Nfe,unkF,Np,N2d,Np2d,Nfl,iel,iElem,ifa,iFace,i,itor,ntorloc,nut
      integer :: c,Np1Dpol,Np1Dtor,blk,Nfdir,Fig,Fe,Fi,sh
      integer :: ind_ue(refElTor%Nnodes3D),ind2(refElPol%Nnodes2D)
      integer :: ind3(refElPol%Nnodes2D*refElTor%Nnodes1D),indl(refElPol%Nnodes1D*refElTor%Nnodes1D)
      real*8,allocatable :: u(:,:),u_tilde(:,:)
      integer :: ierr

      sol%u_tilde = 0.d0
      neq = phys%Neq
      N2D = Mesh%Nelems                  ! Number of 2D elements
      Np2D = refElPol%Nnodes2D            ! Number of nodes for each 2D element
#ifdef PARALL
      IF (MPIvar%ntor .gt. 1) THEN
        ntorloc = numer%ntor/MPIvar%ntor + 1
      ELSE
        ntorloc = numer%ntor
      ENDIF
#else
      ntorloc = numer%ntor
#endif
      Ne = N2D*ntorloc                  ! Number of 3D elements
      Nf = Mesh%Nfaces
      Np1Dpol = refElPol%Nnodes1D         ! Number of nodes in the 1D poloidal segments
      Np1Dtor = refElTor%Nnodes1D         ! Number of nodes in the 1D toroidal segments
      Nfl = Np1Dpol*Np1Dtor              ! Number of nodes in the lateral faces
      Nfe = refElPol%Nfaces
      unkF = Mesh%ukf
      Np = Np2D*refElTor%Nnodes1D       ! Number of nodes for each 3D element
#ifdef PARALL
      if (MPIvar%ntor .gt. 1) then
        nut = ntorloc*(Nfl*Nf + Np2d*N2d) + Np2d*N2d ! Size of utilde per equation
      else
        nut = ntorloc*(Nfl*Nf + Np2d*N2d) ! Size of utilde per equation
      endif
#else
      nut = ntorloc*(Nfl*Nf + Np2d*N2d) ! Size of utilde per equation
#endif
      !                                                nut  = ntorloc*(Nfl*Nf + Np2D*N2D)  ! Size of utilde per equation
      Nfdir = Mesh%Ndir

      ! Indices
      indl = (/(i,i=1,Nfl)/)
      ind2 = (/(i,i=1,Np2D)/)
      ind3 = (/(i,i=1,Np)/)

      ALLOCATE (u(Ne*Np,neq))
      ALLOCATE (u_tilde(nut,neq))
      u_tilde = 0.d0
      u = transpose(reshape(sol%u,(/neq,Ne*Np/)))

      ! Loop in elements
      c = 0
      DO itor = 1,ntorloc
        ! Poloidal faces
        DO iel = 1,N2D
          iElem = (itor - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + ind2,:) = u(ind_ue(ind2),:)
          c = c + Np2D
        END DO

        ! Toroidal interior faces
        DO iFace = 1,Mesh%Nintfaces
          iel = Mesh%intfaces(iFace,1)
          ifa = Mesh%intfaces(iFace,2)
          iElem = (itor - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + indl,:) = u(ind_ue(refElTor%faceNodes3(ifa,:)),:)
          c = c + Nfl
        END DO

        ! Toroidal exterior faces
        DO iFace = 1,Mesh%Nextfaces
          iel = Mesh%extfaces(iFace,1)
          ifa = Mesh%extfaces(iFace,2)
          IF (.not. Mesh%Fdir(iel,ifa)) THEN
            iElem = (itor - 1)*N2D+iel
            ind_ue = (iElem - 1)*Np + ind3
            u_tilde(c + indl,:) = u(ind_ue(refElTor%faceNodes3(ifa,:)),:)
            c = c + Nfl
          END IF
        END DO
      END DO

#ifdef PARALL
      ! Add solution on toroidal ghost faces
      IF (MPIvar%ntor .gt. 1) THEN
        sh = (Np1Dtor - 1)*Np2d
        DO iel = 1,N2D
          iElem = (ntorloc - 1)*N2D+iel
          ind_ue = (iElem - 1)*Np + ind3
          u_tilde(c + ind2,:) = u(ind_ue(ind2 + sh),:)
          c = c + Np2D
        END DO
      ENDIF
#endif
      sol%u_tilde = reshape(transpose(u_tilde),(/nut*neq/))

      DEALLOCATE (u,u_tilde)
    END SUBROUTINE extractFaceSolution
#else
    SUBROUTINE extractFaceSolution
      integer :: neq,Ne,Nf,Nfe,unkF,Np,Nfp,iElem,ifa,iFace,i
      integer :: ind_uf(1:Mesh%Nnodesperface),faceNodes(1:Mesh%Nnodesperface)
      integer :: ind_ue(1:Mesh%Nnodesperelem)
      logical :: alreadydone(1:Mesh%ukf)
      real*8,allocatable :: u(:,:),u_tilde(:,:)

      sol%u_tilde = 0.d0
      neq = phys%Neq
      Ne = Mesh%Nelems
      Nf = Mesh%Nfaces
      Nfe = refElPol%Nfaces
      unkF = Mesh%ukf
      Np = Mesh%Nnodesperelem
      Nfp = Mesh%Nnodesperface

      ALLOCATE (u(1:Ne*Np,1:neq))
      ALLOCATE (u_tilde(1:Nf*Nfp,1:neq))
      u_tilde = 0.d0
      u = transpose(reshape(sol%u,(/neq,Ne*Np/)))

      DO iFace = 1,Mesh%Nintfaces
        iElem = Mesh%intfaces(iFace,1)
        ifa = Mesh%intfaces(iFace,2)
        ind_ue = (iElem - 1)*Np + (/(i,i=1,Np)/)
        ind_uf = (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
        faceNodes = refElPol%Face_nodes(ifa,:)
        u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
      END DO

      DO iFace = 1,Mesh%Nextfaces
        iElem = Mesh%extfaces(iFace,1)
        ifa = Mesh%extfaces(iFace,2)
        IF (.not. Mesh%Fdir(iElem,ifa)) THEN
          ind_ue = (iElem - 1)*Np + (/(i,i=1,Np)/)
          ind_uf = Mesh%Nintfaces*Nfp + (iFace - 1)*Nfp + (/(i,i=1,Nfp)/)
          faceNodes = refElPol%Face_nodes(ifa,:)
          u_tilde(ind_uf,:) = u(ind_ue(faceNodes),:)
        END IF
      END DO
      sol%u_tilde = reshape(transpose(u_tilde),(/Nf*Nfp*neq/))

      DEALLOCATE (u,u_tilde)
    END SUBROUTINE extractFaceSolution
#endif



  SUBROUTINE add_perturbation()
    integer             :: Np1d,Np2d,Np
    integer             :: itor,itorg,iel,iel2,i,imod,nmod,ieq,indl,iphi,ntorloc,itheta
    integer,allocatable :: ind(:)
    real*8              :: Xe(Mesh%Nnodesperelem,Mesh%Ndim),pertphi,perttheta
#ifdef TOR3D
    real*8              :: tdiv(numer%ntor + 1)
    real*8              :: htor,tel(refElTor%Nnodes1d)
#endif
    real*8 :: phi,theta,amp
    real*8,allocatable  :: u(:,:)

    allocate(u(size(sol%u)/phys%neq,phys%neq))
    u = transpose(reshape(sol%u,[phys%neq,size(sol%u)/phys%neq]))

    amp = 1e-3
    nmod = 10

#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
      ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
#else
    Np = refElPol%Nnodes2D
#endif

    allocate(ind(Np))

#ifdef TOR3D
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1,numer%ntor
      tdiv(i + 1) = i*htor
    END DO

    DO itor = 1,ntorloc
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
#endif
#endif
      DO iel2 = 1,Mesh%Nelems
        iel = iel2
#ifdef TOR3D
        iel = (itor - 1)*Mesh%Nelems+iel2
#endif
        ind = (iel - 1)*Np + (/(i,i=1,Np)/)
        Xe = Mesh%X(Mesh%T(iel2,:),:)*simpar%refval_length
        perttheta = 1.
        DO imod=1,Nmod
#ifdef TOR3D
          DO itheta =1,refElTor%Nnodes1d
            theta = tel(itheta)
            perttheta = 1+amp*(cos(imod*theta))
#endif
            DO iphi = 1,refElPol%Nnodes2D
              phi = atan2(Xe(iphi,2),Xe(iphi,1)-geom%R0)
              pertphi = (1+amp*(cos(imod*phi)))
              indl = iphi
#ifdef TOR3D
              indl = (itheta-1)*refElPol%Nnodes2D+iphi
#endif
              DO ieq = 1,phys%neq
                u(ind(indl),ieq) = u(ind(indl),ieq)*pertphi*perttheta
              END DO ! ieq
            END DO ! iphi
#ifdef TOR3D
          END DO ! itheta
#endif
        END DO ! modes
      END DO ! elements 2d
#ifdef TOR3D
    END DO ! toroidal loop
#endif
    sol%u = col(transpose(u))
    deallocate(ind,u)

  END SUBROUTINE add_perturbation



  SUBROUTINE add_blob()
    integer             :: itor,itorg,iel,i,imod,nmod
    integer             :: ind(refElPol%Nnodes2D)
    real*8              :: Xe(refElPol%Nnodes2D,Mesh%Ndim)
    real*8              :: xmax,xmin,ymax,ymin,xm,ym,smod,rs,xsource,ysource
    real*8              :: dsource(refElPol%Nnodes2D),aux(refElPol%Nnodes2D)
    real*8,allocatable  :: u(:,:)

#ifdef TOR3D
    write(6,*) "Blob perturbation not implemented yet"
    stop
#endif

    allocate(u(size(sol%u)/phys%neq,phys%neq))
    u = transpose(reshape(sol%u,[phys%neq,size(sol%u)/phys%neq]))
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax+xmin)
    ym = 0.5*(ymax+ymin)

    DO iel = 1,Mesh%Nelems
      ind = (iel - 1)*refElPol%Nnodes2D + (/(i,i=1,refElPol%Nnodes2D)/)
      Xe = Mesh%X(Mesh%T(iel,:),:)
      smod = 0.1
      rs = 0.04/simpar%refval_length
      xsource = xm+0.85*(xmax-xm)
      ysource = ym
      dsource   = sqrt((Xe(:,1)-xsource)**2+(Xe(:,2)-ysource)**2)
      aux = -dsource**2/rs**2
      DO i=1,refElPol%Nnodes2D
        if (aux(i).gt.-30) then
          u(ind(i),1) =  u(ind(i),1)+smod*exp(aux(i))
        endif
      END DO
    END DO ! elements 2d
    sol%u = col(transpose(u))

    deallocate(u)
  END SUBROUTINE add_blob


 SUBROUTINE projectSolutionDifferentMeshes(xs)
 use reference_element, only:orthopoly2d,orthopoly1d,permutation_quads
 use linearAlgebra, only:solve_linear_system_sing,solve_linear_system

 real*8, intent(in)  :: xs(:,:)
 real*8, allocatable:: u_prov(:,:),q_prov(:,:)
 integer :: np,neq,iel,ip,nnod,nnod1
 integer :: ind_u(refElPol%Nnodes2D*phys%Neq),ind_q(refElPol%Nnodes2D*phys%Neq*2)
 integer, allocatable :: correl(:)
 real*8  :: Xe(3,2),d,xieta(2),mat(3,3),bc(3),Xq(4,2)
 real*8  :: V(1:refElPol%Nnodes2D, 1:refElPol%Nnodes2D)
 real*8  :: V1(1:refElPol%Nnodes1D, 1:refElPol%Nnodes1D)
 real*8  :: sf(refElPol%Nnodes2D),sf1(refElPol%Nnodes1D),sf2(refElPol%Nnodes1D)
 real*8  :: b1,b2,b3,dd
 real*8  :: a11,a12,a13,a21,a22,a23,a31,a32,a33
 real*8 :: r, s, p(1:refElPol%Nnodes2D),q1(refElPol%Nnodes1D),q2(refElPol%Nnodes1D)
 integer*4 :: i,j,perm(refElPol%Nnodes2D)
 real*8 :: tol


 IF (MPIvar%glob_id .eq. 0) THEN
    write (6,*) "*** Projecting the solution: find corresponding elements"
 ENDIF

 tol = max(abs(Mesh%xmax),abs(Mesh%xmin),abs(Mesh%ymax),abs(Mesh%ymin))*1e-5

 np = size(xs,1)
 neq = phys%Neq
 nnod = refElPol%Nnodes2D
 nnod1= refElPol%Nnodes1D
 allocate(u_prov(neq,np))
 allocate(q_prov(neq*2,np))
 allocate(correl(np))
 correl = 0

 ! Find corresponding element
 do iel = 1,size(Mesh%T,1)
    Xe = Mesh%X(Mesh%T(iel,1:3),:)
    mat(1,1:3) = Xe(1:3,1)
    mat(2,1:3) = Xe(1:3,2)
    mat(3,1:3) = 1.
				dd = (mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2) - mat(1,2)*mat(2,1)*mat(3,3) + &
							  &mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1))

				a11 =   mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
				a12 = - mat(1,2)*mat(3,3) + mat(1,3)*mat(3,2)
				a13 =   mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
				a21 = - mat(2,1)*mat(3,3) + mat(2,3)*mat(3,1)
				a22 =   mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
				a23 = - mat(1,1)*mat(2,3) + mat(1,3)*mat(2,1)
				a31 =   mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
				a32 = - mat(1,1)*mat(3,2) + mat(1,2)*mat(3,1)
				a33 =   mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
				do i=1,np

							b1 = xs(i,1)
							b2 = xs(i,2)
							b3 = 1.
							bc(1) = (a11*b1+a12*b2+a13*b3)/dd
							bc(2) = (a21*b1+a22*b2+a23*b3)/dd
							bc(3) = (a31*b1+a32*b2+a33*b3)/dd
							if ( bc(1)>=-tol .and. bc(2)>=-tol .and. bc(3)>=-tol .and. bc(1)<=1+tol .and. bc(2)<=1+tol .and. bc(3)<=1+tol) then
			       correl(i) = iel
			    endif
				end do

				if (refElPol%elemType==1) then
						 Xe = Mesh%X(Mesh%T(iel,(/1,3,4/)),:)
						 mat(1,1:3) = Xe(1:3,1)
						 mat(2,1:3) = Xe(1:3,2)
						 mat(3,1:3) = 1.
							dd = (mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2) - mat(1,2)*mat(2,1)*mat(3,3) + &
												&mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1))

							a11 =   mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
							a12 = - mat(1,2)*mat(3,3) + mat(1,3)*mat(3,2)
							a13 =   mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
							a21 = - mat(2,1)*mat(3,3) + mat(2,3)*mat(3,1)
							a22 =   mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
							a23 = - mat(1,1)*mat(2,3) + mat(1,3)*mat(2,1)
							a31 =   mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
							a32 = - mat(1,1)*mat(3,2) + mat(1,2)*mat(3,1)
							a33 =   mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
							do i=1,np

										b1 = xs(i,1)
										b2 = xs(i,2)
										b3 = 1.
										bc(1) = (a11*b1+a12*b2+a13*b3)/dd
										bc(2) = (a21*b1+a22*b2+a23*b3)/dd
										bc(3) = (a31*b1+a32*b2+a33*b3)/dd
										if ( bc(1)>=-tol .and. bc(2)>=-tol .and. bc(3)>=-tol .and. bc(1)<=1+tol .and. bc(2)<=1+tol .and. bc(3)<=1+tol) then
									    correl(i) = iel
									 endif
							end do
				endif
 end do


 IF (MPIvar%glob_id .eq. 0) THEN
    write (6,*) "*** Projecting the solution: interpolation at points"
 ENDIF

! call displayVectorInt(correl(1:1000))
! stop
 do ip = 1,np

    iel = correl(ip)
    if (iel==0) then
     write(6,*) "Element not found for point: ", xs(i,:), " Stopping"
     stop
    endif

    ! Find the corresponding point in the reference element (linear approach so far)
    if (refElPol%elemType==0) then
						 Xe = Mesh%X(Mesh%T(iel,1:3),:)
						 d = 0.5*( (Xe(2,1)-Xe(1,1))*(Xe(3,2)-Xe(1,2))-(Xe(3,1)-Xe(1,1))*(Xe(2,2)-Xe(1,2)))

						 xieta(1)= 1./d*( (Xe(3,2)-Xe(1,2))*(xs(ip,1)-0.5*(Xe(2,1)+Xe(3,1)) ) - (Xe(3,1)-Xe(1,1))*(xs(ip,2)-0.5*(Xe(2,2)+Xe(3,2)) ) )
						 xieta(2)= 1./d*( (Xe(2,1)-Xe(1,1))*(xs(ip,2)-0.5*(Xe(2,2)+Xe(3,2)) ) - (Xe(2,2)-Xe(1,2))*(xs(ip,1)-0.5*(Xe(2,1)+Xe(3,1)) ) )
    else if (refElPol%elemType==1) then
						 Xq = Mesh%X(Mesh%T(iel,1:4),:)
						 d = 0.25*( (Xq(2,1)+Xq(3,1)-Xq(1,1)-Xq(4,1))*(Xq(3,2)+Xq(4,2)-Xq(1,2)-Xq(2,2)) - &
						           &(Xq(2,2)+Xq(3,2)-Xq(1,2)-Xq(4,2))*(Xq(3,1)+Xq(4,1)-Xq(1,1)-Xq(2,1)) )

						 xieta(1)= 1./d*( (Xq(3,2)+Xq(4,2)-Xq(1,2)-Xq(2,2))*(xs(ip,1)-0.25*(Xq(1,1)+Xq(2,1)+Xq(3,1)+Xq(4,1)) ) - &
						                & (Xq(3,1)+Xq(4,1)-Xq(1,1)-Xq(2,1))*(xs(ip,2)-0.25*(Xq(1,2)+Xq(2,2)+Xq(3,2)+Xq(4,2)) ) )
						 xieta(2)= 1./d*( (Xq(2,1)+Xq(3,1)-Xq(1,1)-Xq(4,1))*(xs(ip,2)-0.25*(Xq(1,2)+Xq(2,2)+Xq(3,2)+Xq(4,2)) ) - &
						                 &(Xq(2,2)+Xq(3,2)-Xq(1,2)-Xq(4,2))*(xs(ip,1)-0.25*(Xq(1,1)+Xq(2,1)+Xq(3,1)+Xq(4,1)) ) )

    end if

    if (refElPol%elemType==0) then
    ! Vandermond matrix
						 V = 0.d0
						 DO i = 1, nnod
						   r = refElPol%coord2d(i, 1)
						   s = refElPol%coord2d(i, 2)
						   CALL orthopoly2d(r, s, nnod, refElPol%Ndeg, p)
						   V(i, :) = p(:)
						 END DO
    ! Compute shape functions on triangles
						 CALL orthopoly2d(xieta(1), xieta(2), nnod, refElPol%Ndeg, p)
						 CALL solve_linear_system_sing(transpose(V), p, sf)
    else if (refElPol%elemType==1) then
    ! Vandermond matrix
						 V1 = 0.d0
						 DO i = 1, nnod1
						   r = refElPol%coord1d(i)
						   CALL orthopoly1d(r, refElPol%Ndeg, q1)
						   V1(i, :) = q1(:)
						 END DO
    ! Compute shape functions on quads
       CALL orthopoly1d(xieta(1), refElPol%Ndeg, q1)
       CALL orthopoly1d(xieta(2), refElPol%Ndeg, q2)
       CALL solve_linear_system_sing(transpose(V1), q1, sf1)
       CALL solve_linear_system_sing(transpose(V1), q2, sf2)
       perm = permutation_quads(refElPol)
       sf(perm) = col(tensorProduct(sf1,sf2))
    endif



    ! Interpolation
    ind_u =  (iel-1)*nnod*neq+(/(i,i=1,nnod*neq)/)
    ind_q =  (iel-1)*nnod*neq*2+(/(i,i=1,nnod*neq*2)/)
    u_prov(:,ip) = matmul(sf,transpose(reshape(sol%u(ind_u),[neq,nnod])))
    q_prov(:,ip) = matmul(sf,transpose(reshape(sol%q(ind_q),[neq*2,nnod])))


 end do

 deallocate(sol%u,sol%q,sol%u_tilde)
 allocate(sol%u(np*neq))
 allocate(sol%q(np*neq*2))

! sol%u = col(u_prov)
! sol%q = col(q_prov) ! warning, this may fill the stack

 do j=1,np
   do i=1,neq
   sol%u( i+(j-1)*neq ) = u_prov(i,j)
   end do
 end do

 do j=1,np
   do i=1,neq*2
   sol%q( i+(j-1)*neq*2 ) = q_prov(i,j)
   end do
 end do


 deallocate(u_prov,q_prov,correl)

 END SUBROUTINE projectSolutionDifferentMeshes











END MODULE initialization
