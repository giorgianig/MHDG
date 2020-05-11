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
      phys%Neq = 4

      ! number of physical variables
      phys%npv = 4

      ALLOCATE (phys%phyVarNam(phys%npv))
      ALLOCATE (phys%conVarNam(phys%Neq))

      ! Set the name of the physical variables
      phys%phyVarNam(1) = "rho"
      phys%phyVarNam(2) = "Mach"
      phys%phyVarNam(3) = "Vorticity"
      phys%phyVarNam(4) = "Electric potential"

      ! Set the name of the conservative variables
      phys%conVarNam(1) = "rho"
      phys%conVarNam(2) = "Gamma"
      phys%conVarNam(3) = "Vorticity"
      phys%conVarNam(4) = "Electric potential"

      simpar%model = 'N-Gamma-Vorticity'
      simpar%Ndim = 2
#ifdef TOR3D
      simpar%Ndim = 3
#endif
      simpar%Neq = phys%Neq
      ALLOCATE (simpar%physvar_refval(phys%npv))
      ALLOCATE (simpar%consvar_refval(phys%Neq))
      simpar%physvar_refval(1) = simpar%refval_density
      simpar%physvar_refval(2) = 1.
      simpar%physvar_refval(3) = simpar%refval_vorticity
      simpar%physvar_refval(4) = simpar%refval_potential
      simpar%consvar_refval(1) = simpar%refval_density
      simpar%consvar_refval(2) = simpar%refval_momentum
      simpar%consvar_refval(3) = simpar%refval_vorticity
      simpar%consvar_refval(4) = simpar%refval_potential
   END SUBROUTINE

   !*******************************************
   ! Convert physical variable to conservative
   ! variables
   !*******************************************
   SUBROUTINE phys2cons(up,ua)
      real*8,dimension(:,:),intent(in)  :: up
      real*8,dimension(:,:),intent(out) :: ua

      ua(:,1) = up(:,1)
      ua(:,2) = up(:,1)*up(:,2)
      ua(:,3) = up(:,3)
      ua(:,4) = up(:,4)
   END SUBROUTINE phys2cons

   !*******************************************
   ! Convert conservative variable to physical
   ! variables
   !*******************************************
   SUBROUTINE cons2phys(ua,up)
      real*8,dimension(:,:),intent(in)  :: ua
      real*8,dimension(:,:),intent(out) :: up

      up(:,1) = ua(:,1)
      up(:,2) = ua(:,2)/ua(:,1)/sqrt(phys%a)
      up(:,3) = ua(:,3)
      up(:,4) = ua(:,4)
   END SUBROUTINE cons2phys

   SUBROUTINE jacobianMatrices(U,A)
      real*8,intent(in)  :: U(:)
      real*8,intent(out) :: A(:,:)

      A = 0.d0
      A(1,2) = 1.
      A(2,1) = (-1*U(2)**2/U(1)**2 + phys%a)
      A(2,2) = 2*U(2)/U(1)
      IF (switch%convvort) THEN
         A(3,1) = -U(3)*U(2)/U(1)**2
         A(3,2) = U(3)/U(1)
         A(3,3) = U(2)/U(1)
      END IF
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

   !*****************************************
   ! Set the perpendicular diffusion
   !****************************************
   SUBROUTINE setLocalDiff(xy,d_iso,d_ani,Bmod)
      real*8,intent(in)  :: xy(:,:)
      real*8,intent(out) :: d_iso(:,:,:),d_ani(:,:,:)
      real*8              :: Bmod(:)
      real*8              :: iperdiff(size(xy,1))
      real*8              :: coeff(size(Bmod))

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
      coeff = 1./Bmod**2
      !*****************************
      ! Diagonal terms
      !*****************************
      d_iso(1,1,:) = phys%diff_n
      d_iso(2,2,:) = phys%diff_u
      d_iso(3,3,:) = phys%diff_vort
      d_iso(4,4,:) = coeff

      d_ani(1,1,:) = phys%diff_n
      d_ani(2,2,:) = phys%diff_u
      d_ani(3,3,:) = phys%diff_vort
      d_ani(4,4,:) = coeff

      !*****************************
      ! Non diagonal terms
      !*****************************
!   !! term +\Div(1/eta \Grad_par \phi) in the vorticity equation
      d_iso(3,4,:) = 0.
      d_ani(3,4,:) = 1./phys%etapar*phys%c1 ! Negative parallel diffusion
!   !! term -\Div(1/n\Grad_par \n) in the vorticity equation
      d_iso(3,1,:) = 0.
      d_ani(3,1,:) = -1./phys%etapar*phys%c2
!   !! term -\Div(1/B^2\Grad_perp n) in the potential equation
      d_iso(4,1,:) = coeff*phys%Mref
      d_ani(4,1,:) = coeff*phys%Mref

!write(6,*) "diff_pot: ",1./Bmod**2*phys%Mref*0.0232

      if (switch%testcase < 50) then
         d_iso(4,4,:) = phys%diff_pot
         d_ani(4,4,:) = phys%diff_pot
         d_iso(3,4,:) = 0.
         d_ani(3,4,:) = phys%diff_pari ! Negative parallel diffusion
                                                !! term -\Div(1/n\Grad_par \n) in the vorticity equation
         d_iso(3,1,:) = 0.
         d_ani(3,1,:) = -phys%diff_pare
                                                !! term -\Div(1/n 1/B^2\Grad_perp n) in the potential equation
         d_iso(4,1,:) = phys%diff_ee
         d_ani(4,1,:) = phys%diff_ee
      endif

      if (switch%difcor .gt. 0) then
         call computeIperDiffusion(xy,iperdiff)
         d_iso(1,1,:) = d_iso(1,1,:)*iperdiff
         d_iso(2,2,:) = d_iso(2,2,:)*iperdiff
         d_iso(3,3,:) = d_iso(3,3,:)*iperdiff
         d_iso(4,4,:) = d_iso(4,4,:)*iperdiff
      endif

   END SUBROUTINE setLocalDiff

   !*******************************************
   ! Compute local diffusion in points
   !*******************************************
   SUBROUTINE computeIperDiffusion(X,ipdiff)
      real*8,intent(IN)  :: X(:,:)
      real*8,intent(OUT) :: ipdiff(:)
      real*8             :: xcorn,ycorn
      real*8             :: rad(size(X,1))
      real*8             :: h
      real*8,parameter   :: tol = 1.e-6
      integer            :: i

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
         WRITE (6,*) "Case not valid"
         STOP
      END SELECT

                                                                !!**********************************************************
                                                                !! Gaussian around the corner
                                                                !!**********************************************************
      h = 10e-3
      rad = sqrt((X(:,1)*phys%lscale - xcorn)**2 + (X(:,2)*phys%lscale - ycorn)**2)
      ipdiff = 1 + numer%dc_coe*exp(-(2*rad/h)**2)

   END SUBROUTINE computeIperDiffusion

   !***********************************************************************
   !
   !    COMPUTATION OF THE STABILIZATION PARAMETER
   !
   !***********************************************************************

   !*******************************************
   ! Compute the stabilization tensor tau
   !*******************************************
   SUBROUTINE computeTauGaussPoints(up,uc,b,n,iel,ifa,isext,xy,tau)
      real*8,intent(in)  :: up(:),uc(:),b(:),n(:),xy(:)
      real,intent(in) :: isext
      integer,intent(in)  :: ifa,iel
      real*8,intent(out) :: tau(:,:)
      real*8              :: tau_aux(4)
      real*8 :: xc,yc,rad,h,aux,bn,bnorm
      real*8 :: U1,U2,U3,U4
      integer :: ndim
      U1 = uc(1)
      U2 = uc(2)
      U3 = uc(3)
      U4 = uc(4)

      tau = 0.
      ndim = size(n)
      bn = dot_product(b(1:ndim),n)
      bnorm = norm2(b(1:ndim))

      if (numer%stab == 2 .or. numer%stab == 3) then
         tau_aux = abs(uc(2)*bn/uc(1))
#ifdef TOR3D
         if (abs(n(3)) > 0.1) then
            ! Poloidal face
            tau_aux(1) = tau_aux(1) + phys%diff_n*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
            tau_aux(2) = tau_aux(2) + phys%diff_u*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
            tau_aux(3) = tau_aux(3) + 1/phys%etapar*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
            tau_aux(4) = tau_aux(4) + phys%diff_pot*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
            if (abs(isext - 1.) < 1e-8) then
               tau_aux(3) = 1
               tau_aux(4) = phys%etapar
            endif
         else
#endif
            ! Toroidal face
            tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            tau_aux(3) = tau_aux(3) + 1/phys%etapar*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            tau_aux(4) = tau_aux(4) + phys%diff_pot*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            if (abs(isext - 1.) < 1e-8) then
               tau_aux(3) = 1
               tau_aux(4) = phys%etapar
            endif
#ifdef TOR3D
         endif
#endif
      else
         write (6,*) "Wrong stabilization type: ",numer%stab
         stop
      endif
      tau(1,1) = tau_aux(1)
      tau(2,2) = tau_aux(2)
      tau(3,3) = tau_aux(3)
      tau(4,4) = tau_aux(4)
   END SUBROUTINE computeTauGaussPoints

   SUBROUTINE computeTauGaussPoints_matrix(up,uc,b,n,xy,isext,iel,tau)
      real*8,intent(in)  :: up(:),uc(:),b(:),n(:),xy(:),isext
      real*8,intent(out) :: tau(:,:)
      integer,intent(in) :: iel
      real*8              :: bn,bnorm
      real*8,parameter :: eps = 1e-12
      real*8 :: U1,U2,U3,U4
      real*8 :: t2,t3,t4,t5,t6,t7,t8,t9
      real*8 :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
      real*8 :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
      real*8 :: t30,t31,t32,t33,t34,t35,t36,t37,t38,t39
      real*8 :: t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
      real*8 :: t50,t51,t52,t53,t54,t55,t56,t57,t58,t59
      real*8 :: t60,t61,t62,t63,t64,t65,t66,t67,t68,t69
      real*8 :: t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80
      real*8 :: x,y,r,h,coef,r0,rc

      x = xy(1)
      y = xy(2)

      U1 = uc(1)
      U2 = uc(2)
      U3 = uc(3)
      U4 = uc(4)

      bn = dot_product(b,n)
      bnorm = norm2(b)
      !************************************
      !
      ! *****     CONVECTIVE PART  ********
      !
      !************************************
      tau(1,1) = abs((uc(2)*bn)/uc(1))
      tau(2,2) = abs((uc(2)*bn)/uc(1))
      tau(3,3) = abs((uc(2)*bn)/uc(1))
      tau(4,4) = abs((uc(2)*bn)/uc(1))

      !************************************
      !
      ! *****     DIFFUSIVE PART  ********
      !
      !************************************
      tau(1,1) = tau(1,1) + phys%diff_n
      tau(2,2) = tau(2,2) + phys%diff_u
   END SUBROUTINE computeTauGaussPoints_matrix

   !*****************************************
   ! Jacobian matrices
   !****************************************
!                 SUBROUTINE jacobianMatricesBound(u,v,nx,ny,Anp,Anm)
!                 USE LinearAlgebra
!                 real*8,intent(in)  :: u,v,nx,ny
!                 real*8,intent(out) :: Anp(1:3,1:3),Anm(1:3,1:3)
!   real*8              :: An(1:3,1:3),Vm(3,3),Dm(3,3),invVm(3,3)
!
!                        An = 0.d0
!                        An(1,2) = nx
!                        An(1,3) = ny
!                        An(2,1) = phys%a - u**2.d0*nx - u*v*ny
!                        An(2,2) = u*2.d0*nx + v*ny
!                        An(2,3) = u*ny
!                        An(3,1) = (phys%a-v**2)*ny - u*v*nx
!                        An(3,2) = v*nx
!                        An(3,3) = 2*v*ny + u*nx
!
!                        CALL eig(An,Vm,Dm)
!                        CALL invert_matrix(Vm,invVm)
!                        Anp = 0.5*(An+matmul(Vm,matmul(abs(Dm),invVm)) )
!                        Anm = 0.5*(An-matmul(Vm,matmul(abs(Dm),invVm)))

!                 END SUBROUTINE jacobianMatricesBound

   !***********************************************************************
   !
   !                          MAGNETIC FIELD
   !
!                !***********************************************************************
!#ifdef TOR3D
!  SUBROUTINE defineMagneticField(x,y,t,b,divb,drift,Bmod)
!  real*8,intent(in)      :: x(:),y(:),t(:)
!  real*8,intent(out),optional     ::b(:,:),divb(:),drift(:,:),Bmod(:)
!  real*8                  :: xc,yc,R0,q,r
!  real*8                  :: xmax,xmin,ymax,ymin,xm,ym,p,divbp
!  real*8                  :: xx,yy,tt,Bp,Bt,Br,Bz,BB,dmax,B0,xr,yr
!  integer*4               :: i,j,ind,N2D,N1D
!

!  N2d = size(X,1)
!  N1d = size(t,1)
!  xmax = Mesh%xmax
!  xmin = Mesh%xmin
!  ymax = Mesh%ymax
!  ymin = Mesh%ymin
!  xm = 0.5*(xmax+xmin)
!  ym = 0.5*(ymax+ymin)
!
!  xc = -0.5
!  yc = -0.5
!  ! Initialization
!  if (present(b)) then
!     b = 0.
!  endif
!  if (present(divb)) then
!     divb  = 0.
!  endif
!  if (present(divb)) then
!     drift = 0.
!  endif
!  if (present(Bmod)) then
!     Bmod  = 0.
!  endif
!  divbp = 0.
!                        DO i=1,N2d
!                                                 DO j=1,N1d
!                                                     xx = x(i)
!                                                     yy = y(i)
!                                                     tt = t(j)
!                                                     ind = (j-1)*N2d+i
!
!                                                                                        SELECT CASE(switch%testcase)
!                                                                                                        CASE(1)
!                                                                                                        IF (switch%axisym) THEN
!                                                                                                                                WRITE(6,*) "This is NOT an axisymmetric test case!"
!                                                                                                                                stop
!                                                                                                 END IF
!                                                                                                        ! Cartesian case,circular field centered in [xm,ym] in the poloidal plane,Bt = 1
!                           Br = (yy-yc)
!                           Bz = (-xx+xc)
!                           Bt = 1.
!
!                                                                                                CASE(2)
!                                                                                                        IF (.not.switch%axisym) THEN
!                                                                                                                                WRITE(6,*) "This is an axisymmetric test case!"
!                                                                                                                                stop
!                                                                                                 END IF
!                                                                                                        ! Axysimmetric case,circular field centered in [xm,ym] in the poloidal plane,Bt = 1
!                           Br = (yy-ym)/xx
!                           Bz = (-xx+xm)/xx
!                           Bt = 1.
!                           if (present(divb)) then
!                              divbp = -(xm**2*ym+xx**2*ym-xm**2*yy-xx**2*yy+3*ym*yy**2-3*ym**2*yy+ym**3-yy**3-2*xm*xx*ym+2*xm*xx*yy)/(xx**4*((xm-xx)**2/xx**2+(ym-yy)**2/xx**2+1)**(1.5))

!                                                                                                        endif
!                                                                                                        CASE(50:59)
!                                                                                                                  write(6,*) "Error in defineMagneticField: you should not be here!"
!                                                                                                                  STOP
!                                                                                                        CASE(60:69)
!
!                                                                                                                                                ! Circular case with limiter
!                                                                                                                                                R0 = geom%R0
!                                                                                                                                                q  = geom%q
!                                                                                                                                                B0 = 2 ! not influential
!                                                                                                                                                xr = xx*phys%lscale
!                                                                                                                                                yr = yy*phys%lscale
!
!                                                                                                                                                r  = sqrt((xr-R0)**2+yr**2)
!                                                                                                                                                Br = -B0*yr/(xr*q*sqrt(1- (r/R0)**2 ) )
!                                                                                                                                                Bz = B0*(xr-R0)/(xr*q*sqrt(1- (r/R0)**2 ) )
!                                                                                                                                                Bt = B0*R0/xr

!                                                                                                                                                if (present(divb)) then
!                                                                                                                                                   IF (switch%axisym) THEN
!                                                                                                                                                                     divbp = -yy/xx/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
!                                                                                                                                                   ELSE
!                                                                                                                                                                     WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
!                                                                                                                                                                     STOP
!                                                                                                                                                   END IF
!                                                                                                                                                endif
!                                                                                                                                  if (present(divb)) then
!                                                                                                                                                   IF (switch%driftdia) THEN
!                                                                                                                                                                    drift(:,2) =  -1./R0*phys%lscale
!                                                                                                                                                   END IF
!                                                                                                                                                endif
!                                                                                                        CASE DEFAULT
!                                                                                                                                        WRITE(6,*) "Error! Test case not valid"
!                                                                                                                                        STOP
!                                                                                        END SELECT
!
!                                                                                        Bp = sqrt(Br**2+Bz**2)
!                                                     BB = sqrt(Bp**2+Bt**2)
!           if (present(b)) then
!                                                                                                  b(ind,1) = Br/BB
!                                                                                                  b(ind,2) = Bz/BB
!                                                                                                  b(ind,3) = Bt/BB
!           endif
!                                                     if (present(divb)) then
!                                                        divb(ind) = divbp
!                                                     endif
!           if (present(Bmod)) then
!              Bmod(ind) = BB
!           endif
!                                                 END DO
!                        END DO
!  END SUBROUTINE defineMagneticField
!#else
!  SUBROUTINE defineMagneticField(x,y,b,divb,drift,Bmod)
!  real*8,intent(in)      :: x(:),y(:)
!  real*8,intent(out),optional     :: b(:,:),divb(:),drift(:,:),Bmod(:)
!  real*8                  :: xc,yc,R0,q,r(size(x)),B0
!  real*8                  :: Br(size(x)),Bz(size(x)),Bt(size(x)),BB(size(x))
!
!  ! Initialization
!  if (present(b)) then
!     b     = 0.
!  endif
!  if (present(divb)) then
!    divb  = 0.
!  endif
!  if (present(drift)) then
!    drift = 0.
!  endif
!  if (present(Bmod)) then
!     Bmod  = 0.
!  endif
!  BB = 0.
!                                SELECT CASE(switch%testcase)
!                                  CASE(1)
!                                  ! Circular field centered in [xc,yc],n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y)
!                                        xc = 0.
!                                        yc = 0.
!                                 Br = (y-yc)
!                                 Bz = (-x+xc)
!            Bt = 1.
!            BB = sqrt(Br**2+Bz**2+Bt**2)
!            if (present(b)) then
!               b(:,1) = Br/BB
!               b(:,2) = Bz/BB
!            endif
!            if (present(divb)) then
!               divb(:) = 0.
!            endif
!            if (present(Bmod)) then
!               Bmod = BB
!            endif
!!      CASE(2)
!!      ! Axisymmetric case with div(b)~=0
!!                                        xc = 0.
!!                                        yc = 0.
!!            if (present(b)) then
!!                                                         b(:,1) = 1./30.*(x-y**2+2)
!!                                                         b(:,2) = 1./30.*(x*y+y)
!!            endif
!!            if (present(divb)) then
!!               divb(:) = 1./30.+((1./30.)*x-(1./30.)*y**2+1./15.)/x+1./30.*(x+1)
!!            endif
!!                                                                    if (present(Bmod)) then
!!                                                                       Bmod(:) = 1.
!!                                                                    endif
!      CASE(50:59)
!         write(6,*) "Error in defineMagneticField: you should not be here!"
!         STOP
!                                  CASE(60:69)
!
!        ! Circular case with limiter
!        B0 = 1.
!        R0 = geom%R0
!        q  = geom%q
!        r  = phys%lscale*sqrt((x-R0/phys%lscale)**2+y**2)
!        if (present(b)) then
!           b(:,1) = -phys%lscale*y/sqrt(R0**2*q**2+(1-q**2)*r**2)
!           b(:,2) = phys%lscale*(x-R0/phys%lscale)/sqrt(R0**2*q**2+(1-q**2)*r**2)
!        endif
!        if (present(divb)) then
!           IF (switch%axisym) THEN
!               divb(:) = -y/x/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
!           ELSE
!               WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
!               STOP
!           END IF
!        endif
!
!        if (present(drift)) then
!           IF (switch%driftdia) THEN
!              drift(:,2) =  -1./R0*phys%lscale
!           END IF
!        endif
!        if (present(Bmod)) then
!           Bmod(:) = B0/R0
!        endif
!                                  CASE DEFAULT
!                                      WRITE(6,*) "Error! Test case not valid"
!                                      STOP
!                                END SELECT
!
!  END SUBROUTINE defineMagneticField
!#endif

!  SUBROUTINE loadMagneticField()
!                                        USE interpolation
!                                 USE HDF5
!                                 USE HDF5_io_module
!                                        integer        :: i,ierr,ip,jp
!                                        integer(HID_T) :: file_id
!                                        real*8,pointer,dimension(:,:) :: r2D,z2D,flux2D,Br2D,Bz2D,Bphi2D
!                                        real*8,allocatable,dimension(:,:) :: bx,by,bmod,divb,bmodx,bmody,driftx,drifty
!                                        real*8,allocatable,dimension(:)   :: xvec,yvec
!                                        real*8                            :: x,y
!
!                                 WRITE(6,*) "******* Loading magnetic field *******"
!                                        ! Allocate storing space in phys
!                                        ALLOCATE(phys%b(Mesh%Nnodes,Mesh%Ndim))
!                                        ALLOCATE(phys%divb(Mesh%Nnodes))
!                                        ALLOCATE(phys%drift(Mesh%Nnodes,Mesh%Ndim))
!                                        ALLOCATE(phys%flux2d(Mesh%Nnodes))
!
!                                        ! Dimensions of the file storing the magnetic field for West
!                                        ip = 541
!                                        jp = 391
!                                        ALLOCATE(r2D(ip,jp))
!                                        ALLOCATE(z2D(ip,jp))
!                                        ALLOCATE(flux2D(ip,jp))
!                                        ALLOCATE(Br2D(ip,jp))
!                                        ALLOCATE(Bz2D(ip,jp))
!                                        ALLOCATE(Bphi2D(ip,jp))
!                                        ALLOCATE(bx(ip,jp))
!                                        ALLOCATE(by(ip,jp))
!                                        ALLOCATE(bmod(ip,jp))
!                                        ALLOCATE(bmodx(ip,jp))
!                                        ALLOCATE(bmody(ip,jp))
!                                        ALLOCATE(divb(ip,jp))
!                                        ALLOCATE(driftx(ip,jp))
!                                        ALLOCATE(drifty(ip,jp))
!
!                                        ! Read file
!                                        CALL HDF5_open('WEST_far_465.h5',file_id,IERR)
!                                        CALL HDF5_array2D_reading(file_id,r2D,'r2D')
!                                        CALL HDF5_array2D_reading(file_id,z2D,'z2D')
!                                        CALL HDF5_array2D_reading(file_id,flux2D,'flux2D')
!                                        CALL HDF5_array2D_reading(file_id,Br2D,'Br2D')
!                                        CALL HDF5_array2D_reading(file_id,Bz2D,'Bz2D')
!                                        CALL HDF5_array2D_reading(file_id,Bphi2D,'Bphi2D')
!                                        CALL HDF5_close(file_id)
!
!                                 ! Apply length scale
!                                 r2D = r2D/phys%lscale
!                                 z2D = z2D/phys%lscale
!
!                                        ! Compute b
!                                        bmod = sqrt(Br2D**2+Bz2D**2+Bphi2D**2)
!                                        bx = -Br2D/bmod
!                                        by = -Bz2D/bmod
!
!                                        ! Compute divergence of b
!                                        divb = 0.
!                                        IF (switch%axisym) THEN
!                                                ! 1/r*(d r*br/dr)+dbz/dz
!                                                  divb(2:ip-1,2:jp-1) = 1./r2D(2:ip-1,2:jp-1)*(r2D(2:ip-1,3:jp)*bx(2:ip-1,3:jp)- &
!                                                                     r2D(2:ip-1,1:jp-2)*bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-  &
!                                                                     r2D(2:ip-1,1:jp-2))+(by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!
!                                        ELSE
!                                        ! dbr/dr+dbz/dz
!                                                  divb(2:ip-1,2:jp-1) = (bx(2:ip-1,3:jp-1)-bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))+ &
!                                                                       (by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!                                        END IF
!
!                                                ! Compute drift velocity
!                                                driftx = 0.
!                                                drifty = 0.
!                                                IF (switch%driftdia) THEN
!                                                                bmodx = (bmod(2:ip-1,3:jp)-bmod(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))
!                                                                bmody = (bmod(3:ip,2:jp-1)-bmod(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!                                                                driftx(2:ip-1,2:jp-1) =  -Bphi2D(2:ip-1,2:jp-1)*bmody/bmod(2:ip-1,2:jp-1)**3
!                                                                drifty(2:ip-1,2:jp-1) =   Bphi2D(2:ip-1,2:jp-1)*bmodx/bmod(2:ip-1,2:jp-1)**3
!                    END IF
!
!                                                ! Interpolate
!                                                ALLOCATE(xvec(jp))
!                                                ALLOCATE(yvec(ip))
!                                                xvec = r2D(1,:)
!                                                yvec = z2D(:,1)
!                                                DO i = 1,Mesh%Nnodes
!                                                                        x = Mesh%X(i,1)
!                                                                        y = Mesh%X(i,2)
!                                                                        phys%b(i,1) = interpolate( ip,yvec,jp,xvec,bx,y,x,1e-12)
!                                                                        phys%b(i,2) = interpolate( ip,yvec,jp,xvec,by,y,x,1e-12)
!                                                                        phys%divb(i) = interpolate( ip,yvec,jp,xvec,divb,y,x,1e-12)
!                                                                        phys%drift(i,1) = interpolate( ip,yvec,jp,xvec,driftx,y,x,1e-12)
!                                                                        phys%drift(i,2) = interpolate( ip,yvec,jp,xvec,drifty,y,x,1e-12)
!                                                                        phys%flux2D(i) = interpolate( ip,yvec,jp,xvec,flux2D,y,x,1e-12)
!                                                END DO
!
!                                        ! Free memory
!                                        DEALLOCATE(Br2D,Bz2D,Bphi2D,xvec,yvec)
!                                        DEALLOCATE(r2D,z2D,flux2D,bx,by,bmod,bmodx,bmody,divb,driftx,drifty)
!
!  END SUBROUTINE loadMagneticField

!
!
!  SUBROUTINE loadMagneticFieldTemporalEvolution()
!                                 USE HDF5
!                                 USE HDF5_io_module
!                                 USE MPI_OMP
!                                        integer        :: i,ierr,k
!                                 character(LEN=20) :: fname = 'Evolving_equilibrium'
!                                 character(10)  :: npr,nid,nit
!                                 character(len=1000) :: fname_complete
!                                        integer(HID_T) :: file_id
!
!                                 WRITE(6,*) "******* Loading magnetic field *******"
!
!                                        ! Allocate storing space in phys
!                                        ALLOCATE(phys%Br(Mesh%Nnodes))
!                                        ALLOCATE(phys%Bz(Mesh%Nnodes))
!                                        ALLOCATE(phys%Bt(Mesh%Nnodes))
!                                        ALLOCATE(phys%flux2d(Mesh%Nnodes))
!
!     ! File name
!     write(nit,"(i10)") time%it
!     nit = trim(adjustl(nit))
!     k = INDEX(nit," ") -1
!
!                                        IF (MPIvar%glob_size.GT.1) THEN
!                                           write(nid,*) MPIvar%glob_id+1
!                                           write(npr,*) MPIvar%glob_size
!                                           fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0",4 - k)//trim(ADJUSTL(nit))//'.h5'
!                                        ELSE
!                                           fname_complete = trim(adjustl(fname))//'_'//REPEAT("0",4 - k)//trim(ADJUSTL(nit))//'.h5'
!                                        END IF
!
!     write(6,*) 'Magnetic field loaded from file: ',trim(adjustl(fname_complete))
!
!                                        ! Read file
!                                        CALL HDF5_open(fname_complete,file_id,IERR)
!                                        CALL HDF5_array1D_reading(file_id,phys%Br,'Br')
!                                        CALL HDF5_array1D_reading(file_id,phys%Bz,'Bz')
!                                        CALL HDF5_array1D_reading(file_id,phys%Bt,'Bt')
!                                        CALL HDF5_array1D_reading(file_id,phys%flux2D,'flux')
!                                        CALL HDF5_close(file_id)
!
!  END SUBROUTINE loadMagneticFieldTemporalEvolution
!

END MODULE physics
