!*****************************************
! project: MHDG
! file: physics.f90
! date: 15/02/2017
! Define the physics of the model
!  ******** N-Gamma system Isothermal ****
!*****************************************
MODULE physics

  USE globals
  USE magnetic_field
  IMPLICIT NONE

CONTAINS

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE initPhys()

    ! number of equation of the problem
    phys%Neq = 2
#ifdef NEUTRAL
    phys%Neq = 3
#endif

    ! number of physical variables
    phys%npv = 3
#ifdef NEUTRAL
    phys%npv = 4
#endif

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho"
    phys%phyVarNam(2) = "u"
    phys%phyVarNam(3) = "Mach"
#ifdef NEUTRAL
    phys%phyVarNam(4)= "rhon"   ! density neutral
#endif

    ! Set the name of the conservative variables
    if (switch%logrho) then
      phys%conVarNam(1) = "log(rho)"
    else
      phys%conVarNam(1) = "rho"
    endif
    phys%conVarNam(2) = "Gamma"
#ifdef NEUTRAL
    phys%conVarNam(3) = "rhon"  ! U3 = rhon
#endif

    simpar%model = 'N-Gamma'
#ifdef NEUTRAL
    simpar%model = 'N-Gamma-Neutral'
#endif
    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq
    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
    simpar%physvar_refval(1) = simpar%refval_density
    simpar%physvar_refval(2) = simpar%refval_speed
    simpar%physvar_refval(3) = 1.
#ifdef NEUTRAL
    simpar%physvar_refval(4) = simpar%refval_neutral
#endif
    if (switch%logrho) then
      simpar%consvar_refval(1) = log(simpar%refval_density)
    else
      simpar%consvar_refval(1) = simpar%refval_density
    endif

    simpar%consvar_refval(2) = simpar%refval_momentum
#ifdef NEUTRAL
    simpar%consvar_refval(3) = simpar%refval_neutral
#endif

  END SUBROUTINE

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up, ua)
    real*8, dimension(:, :), intent(in)  :: up
    real*8, dimension(:, :), intent(out) :: ua

    if (switch%logrho) then
      ua(:, 1) = log(up(:, 1))
    else
      ua(:, 1) = up(:, 1)
    endif
    ua(:, 2) = up(:, 1)*up(:, 2)
#ifdef NEUTRAL
    ua(:,3) = up(:,4)
#endif
  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua, up)
    real*8, dimension(:, :), intent(in)  :: ua
    real*8, dimension(:, :), intent(out) :: up
    real*8 :: dens(size(up,1))
    integer :: i
    
    if (switch%logrho) then
      up(:, 1) = exp(ua(:, 1))                ! density
    else
      up(:, 1) = ua(:, 1)                     ! density
    endif
    dens=up(:, 1)
    if (switch%thresh.ne.0) then
       do i=1,size(up,1)
          dens(i) = max(dens(i),numer%thr)
       end do
    endif
    up(:, 2) = ua(:, 2)/dens              ! u parallel
    up(:, 3) = ua(:, 2)/dens/sqrt(phys%a) ! Mach
#ifdef NEUTRAL
    up(:,4) = abs(ua(:,3))                    ! density neutral
#endif
  END SUBROUTINE cons2phys

   ! ******************************
   ! Split diffusion terms
   ! ******************************
   SUBROUTINE compute_W2(U,W2)
     real*8, intent(IN) :: U(:)
     real*8             :: W2(:)
     W2 = 0.
     W2(1) = (phys%diff_n-phys%diff_u)*U(2)/U(1)
   END SUBROUTINE compute_W2

   SUBROUTINE compute_dW2_dU(U,res)
   real*8, intent(IN) :: U(:)
   real*8             :: res(:,:)
   res = 0.
   res(1,1) = -U(2)*(phys%diff_n-phys%diff_u)/(U(1)**2)
   res(1,2) = 1.*(phys%diff_n-phys%diff_u)/U(1)
   END SUBROUTINE compute_dW2_dU

  SUBROUTINE jacobianMatrices(U, A)
    real*8, intent(in)  :: U(:)
    real*8, intent(out) :: A(:, :)
    real*8 :: dens,gamm

    A = 0.d0
    if (switch%logrho) then
      dens=exp(U(1))
      A(1,1) = -U(2)/dens
      A(1,2) = 1/dens
      A(2,1) = -U(2)**2/dens+phys%a*dens
      A(2,2) = 2*U(2)/dens
    else
      dens=U(1)
!      if (switch%thresh.ne.0) then
!         dens = max(dens,numer%thr)
!      endif
      A(1, 2) = 1.
      A(2, 1) = (-1*U(2)**2/dens**2 + phys%a)
      A(2, 2) = 2*U(2)/dens
    endif

  END SUBROUTINE jacobianMatrices

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U,bn,An)
    real*8,intent(in)  :: U(:),bn
    real*8,intent(out) :: An(:,:)

    An = 0.
    CALL jacobianMatrices(U,An)
    An = bn*An
  END SUBROUTINE jacobianMatricesFace


  SUBROUTINE logrhojacobianVector(U,Up,V)
    real*8, intent(in)  :: U(:),Up(:)
    real*8, intent(out) :: V(:)
    real*8 :: dens

    V = 0.d0
    V(1) = U(2)*U(1)/Up(1)
    V(2) = (U(1)-1)*(U(2)**2/Up(1) - phys%a*Up(1))

  END SUBROUTINE logrhojacobianVector


  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani, Bmod)
    real*8, intent(in)  :: xy(:, :)
    real*8, intent(in)  :: Bmod(:)
    real*8, intent(in)  :: u(:,:) 
    real*8, intent(out) :: d_iso(:, :, :), d_ani(:, :, :)
    real*8              :: iperdiff(size(xy, 1))

    ! d_iso(Neq,Neq,Ngauss),d_ani(Neq,Neq,Ngauss)
    ! first index corresponds to the equation
    ! second index corresponds to the unknown
    ! third index correspond to the gauss point
    ! example d_iso(2,3,ig) is the diffusion in the non-diagonal
    ! diffusion term in the second equation/third variable
    ! d_iso>0,d_ani=0 is an isotropic diffusion
    ! d_iso>0,d_ani=d_iso is a perpendicular diffusion
    ! d_iso=0,d_ani<0 is a (positive) parallel diffusion

    d_iso = 0.
    d_ani = 0.
    !*****************************
    ! Diagonal terms
    !*****************************
    d_iso(1, 1, :) = phys%diff_n
    d_iso(2, 2, :) = phys%diff_u
#ifdef NEUTRAL
    d_iso(3,3,:) = phys%diff_nn
#endif

    d_ani(1, 1, :) = phys%diff_n
    d_ani(2, 2, :) = phys%diff_u
#ifdef NEUTRAL
    d_ani(3,3,:) = 0.
#endif


!write(6,*) "inside d_iso", d_iso(1,1,1),d_iso(2,2,1),d_iso(3,3,1)


    !*****************************
    ! Non diagonal terms
    !*****************************
    ! No non-diagonal terms defined for this model

!    call computeIperDiffusion(xy, u, iperdiff)
!    d_iso(1, 1, :) = d_iso(1, 1, :)*iperdiff
!    d_iso(2, 2, :) = d_iso(2, 2, :)*iperdiff
!#ifdef NEUTRAL    
!    d_iso(3, 3, :) = d_iso(3, 3, :)*iperdiff
!#endif
 
    call computeIperDiffusion(xy, u, iperdiff)
    d_iso(1, 1, :) = d_iso(1, 1, :) + iperdiff
    d_iso(2, 2, :) = d_iso(2, 2, :) + iperdiff
      
#ifdef NEUTRAL    
!    d_iso(3, 3, :) = d_iso(3, 3, :) + iperdiff
#endif

!    if (maxval(iperdiff)>1e-12) then
!    write(6,*) "iperdiff: ", iperdiff
!    endif
!write(6,*) "u: ", u, "iperdiff ", iperdiff

  END SUBROUTINE setLocalDiff

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X, u, ipdiff)
    real*8, intent(IN)  :: X(:, :), u(:,:)
    real*8, intent(OUT) :: ipdiff(:)
    real*8             :: xcorn, ycorn,d,dref,maxamp
    real*8             :: rad(size(X, 1))
    real*8             :: h,rhog
    real*8, parameter   :: tol = 1.e-6
    integer            :: i,g, opt

    ipdiff = 0.
    
    if (switch%difcor .gt. 0) then
						 SELECT CASE (switch%difcor)
						 CASE (1)
						   ! Circular case with infinitely small limiter
						   xcorn = geom%R0
						   ycorn = -0.75
						 CASE (2)
						   ! Circular case with infinitely small limiter
						   xcorn = geom%R0
						   ycorn = -0.287
						 CASE (3)
						   ! West
						   xcorn = 2.7977
						   ycorn = -0.5128
						 CASE DEFAULT
						   WRITE (6, *) "Case not valid"
						   STOP
						 END SELECT

						 !!**********************************************************
						 !! Gaussian around the corner
						 !!**********************************************************
						 h = 10e-3
						 rad = sqrt((X(:, 1)*phys%lscale - xcorn)**2 + (X(:, 2)*phys%lscale - ycorn)**2)
						 ipdiff = 1 + numer%dc_coe*exp(-(2*rad/h)**2)
    end if
    
    maxamp = 4.
    opt = 2
    if (switch%limrho .eq. 2 .and. minval(u(:,1)) .lt. numer%minrho    ) then
       do g = 1,size(u,1)
          rhog = u(g,1)
          if (rhog<0.) rhog=0.
          if (rhog<numer%minrho) then
						       d = numer%minrho-rhog ! 0 < d < minrho
						       if (opt.eq.1) then
						          dref = maxamp * d/numer%minrho  ! 0 < dref < maxamp
!			            ipdiff(g) = exp(dref) ! 1 < ipdiff(g) < exp(maxamp)
						          ipdiff(g) = ipdiff(g) + exp(dref) - 1 ! 0 < ipdiff(g) < exp(maxamp)-1
						       else if (opt==2) then
						          ipdiff(g) = ipdiff(g) +  1./((1.-d/numer%minrho)*2+1./50.) - 1.
						       endif
          endif
       end do
    endif
  END SUBROUTINE computeIperDiffusion

  ! ******************************
  ! Neutral Source terms
  ! ******************************
#ifdef NEUTRAL

  SUBROUTINE compute_niz(U,niz)
    real*8, intent(IN) :: U(:)
    real*8             :: niz,U1,U3
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U3 = U(3)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    niz = U1*U3
  END SUBROUTINE compute_niz


  SUBROUTINE compute_dniz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U3
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U3 = U(3)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    res = 0.
    res(1) = U3
    res(3) = U1
  END SUBROUTINE compute_dniz_dU


  SUBROUTINE compute_nrec(U,nrec)
    real*8, intent(IN) :: U(:)
    real*8             :: nrec,U1
    real,parameter :: tol = 1e-10
    U1 = U(1)
    if (U1<tol) U1=tol
    nrec = U1**2
  END SUBROUTINE compute_nrec


  SUBROUTINE compute_dnrec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1
    real,parameter :: tol = 1e-10
    U1 = U(1)
    if (U1<tol) U1=tol
    res = 0.
    res(1) = 2.*U1
  END SUBROUTINE compute_dnrec_dU


  SUBROUTINE compute_fGammacx(U,fGammacx)
    real*8, intent(IN) :: U(:)
    real*8             :: fGammacx,U2,U3
    real,parameter :: tol = 1e-10
    U2 = U(2)
    U3 = U(3)
    if (U3<tol) U3=tol
    fGammacx = U2*U3
  END SUBROUTINE compute_fGammacx


  SUBROUTINE compute_dfGammacx_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U2,U3
    real,parameter :: tol = 1e-10
    U2 = U(2)
    U3 = U(3)
    if (U3<tol) U3=tol
    res = 0.
    res(2) = U3
    res(3) = U2
  END SUBROUTINE compute_dfGammacx_dU


  SUBROUTINE compute_fGammarec(U,fGammarec)
    real*8, intent(IN) :: U(:)
    real*8             :: fGammarec,U1,U2
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    if (U1<tol) U1=tol
    fGammarec = U1*U2
  END SUBROUTINE compute_fGammarec


  SUBROUTINE compute_dfGammarec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U2
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    if (U1<tol) U1=tol
    res = 0.
    res(1) = U2
    res(2) = U1
  END SUBROUTINE compute_dfGammarec_dU

#endif 
!NEUTRAL




  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up, uc, b, bmod, n, iel, ifa, isext, xy, tau)
    real*8, intent(in)  :: up(:), uc(:), b(:), bmod,n(:), xy(:)
    real, intent(in) :: isext
    integer, intent(in)  :: ifa, iel
    real*8, intent(out) :: tau(:, :)
    integer             :: ndim
#ifdef NEUTRAL
    real*8              :: tau_aux(3),diff_iso(3,3,1),diff_ani(3,3,1)
#else
    real*8              :: tau_aux(2),diff_iso(2,2,1),diff_ani(2,2,1)
#endif
    real*8 :: xc, yc, rad, h, aux, bn, bnorm,xyd(1,size(xy)),uu(1,size(uc)),bmod_d(1)

    ndim = size(n)
    bn = dot_product(b(1:ndim), n)
    bnorm = norm2(b(1:ndim))
    xyd(1,:) = xy(:)
    uu(1,:) = uc(:)
    bmod_d = bmod
    
    call setLocalDiff(xyd, uu, diff_iso, diff_ani, bmod_d)
    
    
!    write(6,*) "diff_iso", diff_iso(1,1,1),diff_iso(2,2,1),diff_iso(3,3,1), " phys%diff_n ", phys%diff_n, " phys%diff_u ",phys%diff_u,  " phys%diff_nn ",phys%diff_nn
    
    if (numer%stab == 2) then
!      tau_aux = abs(up(2)*bn)
!      tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!      tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!#ifdef NEUTRAL
!      tau_aux(3) =              phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!!      tau_aux(3) = tau_aux(3) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!#endif
      tau_aux = abs(up(2)*bn)
      tau_aux(1) = tau_aux(1) + diff_iso(1,1,1)
      tau_aux(2) = tau_aux(2) + diff_iso(2,2,1)
#ifdef NEUTRAL
      tau_aux(3) =              diff_iso(3,3,1)
#endif

    elseif (numer%stab == 3) then
      tau_aux = max(abs((uc(2) + sqrt(phys%a))*bn/up(1)), abs((uc(2) - sqrt(phys%a))*bn/up(1)))
      tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
      tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
      tau_aux(3) = phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif

    elseif (numer%stab == 4) then
      tau_aux = abs((up(2)*bn))
      tau_aux(1) = tau_aux(1) + phys%diff_n
      tau_aux(2) = tau_aux(2) + phys%diff_u
#ifdef NEUTRAL
      tau_aux(3) =  phys%diff_nn
#endif      

    elseif (numer%stab == 5) then
      tau_aux = max(abs((uc(2) + sqrt(phys%a))*bn/up(1)), abs((uc(2) - sqrt(phys%a))*bn/up(1)))
      tau_aux(1) = tau_aux(1) + phys%diff_n
      tau_aux(2) = tau_aux(2) + phys%diff_u
#ifdef NEUTRAL
      tau_aux(3) =  phys%diff_nn
#endif      

    else
      write (6, *) "Wrong stabilization type: ", numer%stab
      stop
    endif
    tau(1, 1) = tau_aux(1)
    tau(2, 2) = tau_aux(2)
#ifdef NEUTRAL
    tau(3,3) = tau_aux(3)
#endif
  END SUBROUTINE computeTauGaussPoints

  SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)
    real*8, intent(in)  :: up(:), uc(:), b(:), n(:), xy(:), isext
    real*8, intent(out) :: tau(:, :)
    integer, intent(in) :: iel
    real*8              :: bn, bnorm
    real*8, parameter :: eps = 1e-12
    real*8 :: U1, U2, U3, U4
    real*8 :: t2, t3, t4, t5, t6, t7, t8, t9
    real*8 :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    real*8 :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    real*8 :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    real*8 :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    real*8 :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
    real*8 :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
    real*8 :: t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80
    real*8 :: x, y, r, h, coef, r0, rc

    x = xy(1)
    y = xy(2)

    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    bn = dot_product(b, n)
    bnorm = norm2(b)
    !************************************
    !
    ! *****     CONVECTIVE PART  ********
    !
    !************************************
    tau(1, 1) = abs((uc(2)*bn)/uc(1))
    tau(2, 2) = abs((uc(2)*bn)/uc(1))
    tau(3, 3) = abs((uc(2)*bn)/uc(1))
    tau(4, 4) = abs((uc(2)*bn)/uc(1))

    !************************************
    !
    ! *****     DIFFUSIVE PART  ********
    !
    !************************************
    tau(1, 1) = tau(1, 1) + phys%diff_n
    tau(2, 2) = tau(2, 2) + phys%diff_u
  END SUBROUTINE computeTauGaussPoints_matrix


END MODULE physics
