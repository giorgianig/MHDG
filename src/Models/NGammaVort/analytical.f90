!*****************************************
! project: MHDG
! file: analytical.f90
! date: 15/02/2017
!  ******** N-Gamma system Isothermal ****
! Define the analytical solution and the
! body forces for different cases
!*****************************************
MODULE analytical
  USE prec_const
  USE globals,ONLY: switch
  USE physics
  USE PrintUtils
  IMPLICIT NONE

CONTAINS

#ifdef TOR3D

  !***********************************************************************
  !
  !                            VERSION 3D TOROIDAL
  !
  !***********************************************************************

  !*****************************************
  ! Analytical solution
  !****************************************
  SUBROUTINE analytical_solution(x,y,t,u)
    real*8,dimension(:),intent(IN)        :: x,y,t
    real*8,dimension(:,:),intent(OUT)     :: u
    real*8,dimension(size(u,1),size(u,2))  :: up
    integer:: i,j,k,ind,N2D,N1D
    real*8 :: a,b,r,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
    real*8 :: dsource(refElPol%Nnodes2D),aux(refElPol%Nnodes2D),xsource,ysource,smod

    u = 0.
    up = 0.
    a = 2*pi
    N2D = size(x,1)
    N1D = size(t,1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    DO i = 1,N2d
      DO j = 1,N1d  ! TODO: check if I need to invert i and j
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          IF (switch%axisym) THEN
            WRITE (6,*) "This is NOT an axisymmetric test case!"
            stop
          END IF
          ! Cartesian case,circular field centered in [xm,ym] in the poloidal plane,Bt = 1
          up(ind,1) = 2 + sin(a*xx)*sin(a*yy)
          up(ind,2) = cos(a*xx)*cos(a*yy)
          up(ind,3) = cos(a*xx)*sin(a*yy)
          up(ind,4) = sin(a*xx)*cos(a*yy)

        CASE (50:64)
          up(ind,1) = 1.
          up(ind,2) = 0.
          up(ind,3) = 0.
          up(ind,4) = 0.
        CASE (65)
          up(ind,1) = 1.
          up(ind,2) = 0.
          r = sqrt((xx*phys%lscale - geom%R0)**2 + (yy*phys%lscale - 0.75)**2)
          IF (r .le. 0.05) THEN
            up(ind,2) = 1.
          END IF
        CASE DEFAULT
          WRITE (6,*) "Error! Test case not valid"
          STOP
        END SELECT
      END DO
    END DO
    ! Convert physical variables to conservative variables
    CALL phys2cons(up,u)
  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x,y,t,u,ux,uy,ut)
    real*8,dimension(:),intent(IN)        :: x,y,t
    real*8,dimension(:,:),intent(IN)      :: u
    real*8,dimension(:,:),intent(OUT)     :: ux,uy,ut
    real*8,dimension(size(u,1),size(u,2))  :: upx,upy,upt
    real*8,dimension(size(u,1),phys%npv)   :: up
    integer:: i,j,ind,N1D,N2D
    real*8 :: a,b,r,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
    real*8 :: aux

    N2D = size(x,1)
    N1D = size(t,1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)
    upx = 0.
    upy = 0.
    upt = 0.
    CALL cons2phys(u,up)
    a = 2*pi

    DO i = 1,N2d
      DO j = 1,N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          IF (switch%axisym) THEN
            WRITE (6,*) "This is NOT an axisymmetric test case!"
            stop
          END IF
          ! Circular field centered in [xc,yc],n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y),Ei = 20+cos(wx*x)*sin(wy*y),Ee = 10-sin(wx*x)*cos(wy*y)
          upx(ind,1) = a*cos(a*xx)*sin(a*yy)
          upx(ind,2) = -a*sin(a*xx)*cos(a*yy)
          upx(ind,3) = -a*sin(a*xx)*sin(a*yy)
          upx(ind,4) = a*cos(a*xx)*cos(a*yy)

          upy(ind,1) = a*sin(a*xx)*cos(a*yy)
          upy(ind,2) = -a*cos(a*xx)*sin(a*yy)
          upy(ind,3) = a*cos(a*xx)*cos(a*yy)
          upy(ind,4) = -a*sin(a*xx)*sin(a*yy)

        CASE (2)
          ! Axisimmetric case with div(b)~=0,n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y),Ei = 20+cos(wx*x)*sin(wy*y),Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6,*) "This is an axisymmetric test case!"
            stop
          END IF
          upx(ind,1) = a*cos(a*xx)*sin(a*yy)
          upx(ind,2) = -a*sin(a*xx)*cos(a*yy)
          upx(ind,3) = -a*sin(a*xx)*sin(a*yy)
          upx(ind,4) = -a*cos(a*xx)*cos(a*yy)

          upy(ind,1) = a*sin(a*xx)*cos(a*yy)
          upy(ind,2) = -a*cos(a*xx)*sin(a*yy)
          upy(ind,3) = a*cos(a*xx)*cos(a*yy)
          upy(ind,4) = a*sin(a*xx)*sin(a*yy)
        CASE (3)
          ! Axisimmetric case with div(b)~=0,n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y),Ei = 20+cos(wx*x)*sin(wy*y),Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6,*) "This is an axisymmetric test case!"
            stop
          END IF
          upt(ind,1) = +a*cos(a*tt)
          upt(ind,2) = -a*sin(a*tt)
          upt(ind,3) = -a*sin(a*tt)
          upt(ind,4) = -a*cos(a*tt)
        CASE (5)
          ! Do nothing
        CASE (6)
          ! Do nothing
        CASE (50:)
          ! Do nothing
        CASE DEFAULT
          WRITE (6,*) "Error! Test case not valid"
          STOP
        END SELECT
      END DO
    END DO

    ! Convert physical variables to conservative variables
    ux(:,1) = upx(:,1)
    uy(:,1) = upy(:,1)
    ut(:,1) = upt(:,1)
    ux(:,2) = (upx(:,1)*up(:,2) + up(:,1)*upx(:,2))
    uy(:,2) = (upy(:,1)*up(:,2) + up(:,1)*upy(:,2))
    ut(:,2) = (upt(:,1)*up(:,2) + up(:,1)*upt(:,2))
    ux(:,3) = upx(:,3)
    uy(:,3) = upy(:,3)
    ut(:,3) = upt(:,3)
    ux(:,4) = upx(:,4)
    uy(:,4) = upy(:,4)
    ut(:,4) = upt(:,4)

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x,y,t,f)
    real*8,dimension(:),intent(IN) :: x,y,t
    real*8,dimension(:,:),intent(OUT) :: f
    real*8  :: xc,yc,D,mu,k
    integer:: i,j,N2D,N1D,ind
    real*8 :: a,b,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
    real*8  :: t1,t2,t3,t4,t5,t6,t7,t8,t9,ta,tb,cx2,cy2
    real*8  :: Dw,DparPhi,DPhi,DnDivn,DparN
    real*8 :: r,r2,cx,cy,sx,sy,driftexb

    a = 2*pi
    N2D = size(x,1)
    N1D = size(t,1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    f = 0.
    DO i = 1,N2d
      DO j = 1,N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          ! Circular field centered in [xc,yc],n = 2+sin(2*pi*x)*sin(2*pi*y), u = cos(2*pi*x)*cos(2*pi*y)
          IF (switch%axisym) THEN
            WRITE (6,*) "This is NOT an axisymmetric test case!"
            stop
          END IF
          a = 2*pi
          b = 2*pi
          xc = 0.
          yc = 0.
          D = phys%diff_n
          mu = phys%diff_u
          k = phys%a
          Dw = phys%diff_vort
          DPhi = phys%diff_pot
          DparPhi = phys%diff_pari
          DparN = phys%diff_pare
          DnDivn = phys%diff_ee
          r = ((xx - xc)**2 + (yy - yc)**2 + 1)**0.5
          r2 = ((xx - xc)**2 + (yy - yc)**2 + 1)
          cx = cos(a*xx)
          cy = cos(a*yy)
          sx = sin(a*xx)
          sy = sin(a*yy)

          if (switch%driftexb) then
            driftexb = 1.
          else
            driftexb = 0.
          endif

          f(ind,1) =   driftexb*((a**2*cx*cy**2*sx)/r2 + (a**2*cx*sx*sy**2)/r2 - (a*sx*sy*(2*xx - 2*xc)*(sx*sy + 2))/r2**2 - (a*cx*cy*(2*yy - 2*yc)*(sx*sy + 2))/r2**2) - D*(((yy - yc)*((a*cy*sx)/r + (a**2*sx*sy*(yy - yc))/r + (a**2*cx*cy*(xx - xc))/r - (a*cy*sx*(2*xx - 2*xc)*(xx - xc))/(2*r2**1.5) + (a*cx*sy*(2*xx - 2*xc)*(yy - yc))/(2*r2**1.5)))/r + ((xx - xc)*((a*cx*sy)/r + (a**2*sx*sy*(xx - xc))/r + (a**2*cx*cy*(yy - yc))/r + (a*cy*sx*(2*yy - 2*yc)*(xx - xc))/(2*r2**1.5) - (a*cx*sy*(2*yy - 2*yc)*(yy - yc))/(2*r2**1.5)))/r - 2*a**2*sx*sy - ((2*xx - 2*xc)*((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(yy - yc))/(2*r2**1.5) + ((2*yy - 2*yc)*((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(xx - xc))/(2*r2**1.5)) - ((cx*cy*(2*xx - 2*xc)*(sx*sy + 2)*(yy - yc))/(2*r2**1.5) - (cx*cy*(2*yy - 2*yc)*(sx*sy + 2)*(xx - xc))/(2*r2**1.5) + (a*cx*cy**2*sx*(xx - xc))/r - (a*cx**2*cy*sy*(yy - yc))/r - (a*cx*sy*(sx*sy + 2)*(xx - xc))/r + (a*cy*sx*(sx*sy + 2)*(yy - yc))/r)
          f(ind,2) =   (mu*(2*a**2*xx*yc - 2*a**2*xx*yy + 2*a**2*xc*yy - 2*a**2*xc*yc + 4*a**2*cx*cy + (a*xx*sin(2*a*yy))/2 + (a*yy*sin(2*a*xx))/2 - (a*xc*sin(2*a*yy))/2 - (a*yc*sin(2*a*xx))/2 + 4*a**2*xx*yy*cx**2 - 4*a**2*xx*yc*cx**2 - 4*a**2*xc*yy*cx**2 + 4*a**2*xc*yc*cx**2 + 4*a**2*xx*yy*cy**2 - 4*a**2*xx*yc*cy**2 - 4*a**2*xc*yy*cy**2 + 4*a**2*xc*yc*cy**2 + 2*a**2*xx**2*cx*cy + 2*a**2*xc**2*cx*cy + 2*a**2*yy**2*cx*cy + 2*a**2*yc**2*cx*cy + 2*a*xx*cy*sx - 2*a*xc*cy*sx + 2*a*yy*cx*sy - 2*a*yc*cx*sy - 8*a**2*xx*yy*cx**2*cy**2 + 8*a**2*xx*yc*cx**2*cy**2 + 8*a**2*xc*yy*cx**2*cy**2 - 8*a**2*xc*yc*cx**2*cy**2 - 4*a**2*xx*xc*cx*cy - 4*a**2*yy*yc*cx*cy - 4*a**2*xx*yy*sx*sy + 4*a**2*xx*yc*sx*sy + 4*a**2*xc*yy*sx*sy - 4*a**2*xc*yc*sx*sy - 2*a*xx*cx**2*cy*sy - 2*a*yy*cx*cy**2*sx + 2*a*xc*cx**2*cy*sy + 2*a*yc*cx*cy**2*sx + 8*a**2*cx*cy*sx*sy + 4*a**2*xx**2*cx*cy*sx*sy + 4*a**2*xc**2*cx*cy*sx*sy + 4*a**2*yy**2*cx*cy*sx*sy + 4*a**2*yc**2*cx*cy*sx*sy - 8*a**2*xx*xc*cx*cy*sx*sy - 8*a**2*yy*yc*cx*cy*sx*sy))/(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1) - driftexb*((a**2*cx**2*cy*sy*(sx*sy + 2))/r2 - (a**2*cx**2*cy*sx*sy**2)/r2 - (a**2*cx**2*cy**3*sx)/r2 + (a**2*cy*sx**2*sy*(sx*sy + 2))/r2 + (a*cx**2*cy**2*(2*yy - 2*yc)*(sx*sy + 2))/r2**2 + (a*cx*cy*sx*sy*(2*xx - 2*xc)*(sx*sy + 2))/r2**2) + (a*(2*xx*cx**2*cy*sx - 2*xc*cx**2*cy*sx + 4*xx*cx**2*cy*sy - 4*yy*cx*cy**2*sx - 4*xc*cx**2*cy*sy + 4*yc*cx*cy**2*sx - 2*yy*cx*cy**2*sy + 2*yc*cx*cy**2*sy - 3*xx*cx**2*cy**3*sx + 3*xc*cx**2*cy**3*sx + 3*yy*cx**3*cy**2*sy - 3*yc*cx**3*cy**2*sy - k*xx*cy*sx + k*xc*cy*sx + k*yy*cx*sy - k*yc*cx*sy))/(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1)**0.5
          f(ind,3) =   (Dw*a*(2*a*cx*sy - yy*cx*cy + yc*cx*cy + xx*sx*sy - xc*sx*sy + a*xx**2*cx*sy + a*xc**2*cx*sy + a*yy**2*cx*sy + a*yc**2*cx*sy - 2*a*xx*xc*cx*sy + 2*a*xx*yy*cy*sx - 2*a*xx*yc*cy*sx - 2*a*xc*yy*cy*sx + 2*a*xc*yc*cy*sx - 2*a*yy*yc*cx*sy))/(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1) - (DparPhi*a*(xx*cx*cy - xc*cx*cy - yy*sx*sy + yc*sx*sy + a*xx**2*cy*sx + a*xc**2*cy*sx + a*yy**2*cy*sx + a*yc**2*cy*sx - 2*a*xx*xc*cy*sx - 2*a*xx*yy*cx*sy + 2*a*xx*yc*cx*sy + 2*a*xc*yy*cx*sy - 2*a*xc*yc*cx*sy - 2*a*yy*yc*cy*sx))/(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1) - driftexb*((a**2*sx**2*sy**2)/r2 - (a**2*cx**2*cy**2)/r2 + (a*cx*sx*sy**2*(2*xx - 2*xc))/r2**2 + (a*cx**2*cy*sy*(2*yy - 2*yc))/r2**2) + (a*cx*(xx*cx - xc*cx - 2*xx*cx*cy**2 + 2*xc*cx*cy**2 - 2*yy*cy*sx*sy + 2*yc*cy*sx*sy))/(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1)**0.5 + (DparN*a*((xx*sin(2*a*xx))/2 - (xc*sin(2*a*xx))/2 + (yy*sin(2*a*yy))/2 - (yc*sin(2*a*yy))/2 + a*xx**2 + a*xc**2 + a*yy**2 + a*yc**2 + 2*xx*cx*sy - 2*xc*cx*sy + 2*yy*cy*sx - 2*yc*cy*sx - a*xx**2*cx**2 - a*xc**2*cx**2 - a*yy**2*cy**2 - a*yc**2*cy**2 - 2*a*xx*xc - 2*a*yy*yc - xx*cx*cy**2*sx + xc*cx*cy**2*sx - yy*cx**2*cy*sy + yc*cx**2*cy*sy + 2*a*xx*xc*cx**2 + 2*a*yy*yc*cy**2 + 2*a*xx**2*sx*sy + 2*a*xc**2*sx*sy + 2*a*yy**2*sx*sy + 2*a*yc**2*sx*sy - 4*a*xx*xc*sx*sy - 4*a*yy*yc*sx*sy + 4*a*xx*yy*cx*cy - 4*a*xx*yc*cx*cy - 4*a*xc*yy*cx*cy + 4*a*xc*yc*cx*cy))/((sx*sy + 2)**2*(xx**2 - 2*yy*yc - 2*xx*xc + xc**2 + yy**2 + yc**2 + 1))
          f(ind,4) =   DnDivn*((a**2*sx*sy - ((yy - yc)*((a*cy*sx)/r + (a**2*sx*sy*(yy - yc))/r + (a**2*cx*cy*(xx - xc))/r - (a*cy*sx*(2*xx - 2*xc)*(xx - xc))/(2*r2**1.5) + (a*cx*sy*(2*xx - 2*xc)*(yy - yc))/(2*r2**1.5)))/r + ((2*xx - 2*xc)*((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(yy - yc))/(2*r2**1.5))/(sx*sy + 2) - (((xx - xc)*((a*cx*sy)/r + (a**2*sx*sy*(xx - xc))/r + (a**2*cx*cy*(yy - yc))/r + (a*cy*sx*(2*yy - 2*yc)*(xx - xc))/(2*r2**1.5) - (a*cx*sy*(2*yy - 2*yc)*(yy - yc))/(2*r2**1.5)))/r - a**2*sx*sy + ((2*yy - 2*yc)*((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(xx - xc))/(2*r2**1.5))/(sx*sy + 2) + (a*cy*sx*(a*cy*sx - (((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(xx - xc))/r))/(sx*sy + 2)**2 + (a*cx*sy*(a*cx*sy + (((a*cy*sx*(xx - xc))/r - (a*cx*sy*(yy - yc))/r)*(yy - yc))/r))/(sx*sy + 2)**2) + DPhi*(((xx - xc)*((a**2*cx*sy*(yy - yc))/r - (a**2*cy*sx*(xx - xc))/r - (a*cx*cy)/r + (a*cx*cy*(2*yy - 2*yc)*(yy - yc))/(2*r2**1.5) + (a*sx*sy*(2*yy - 2*yc)*(xx - xc))/(2*r2**1.5)))/r + 2*a**2*cy*sx - ((yy - yc)*((a**2*cy*sx*(yy - yc))/r - (a**2*cx*sy*(xx - xc))/r - (a*sx*sy)/r + (a*cx*cy*(2*xx - 2*xc)*(yy - yc))/(2*r2**1.5) + (a*sx*sy*(2*xx - 2*xc)*(xx - xc))/(2*r2**1.5)))/r - ((2*xx - 2*xc)*((a*cx*cy*(yy - yc))/r + (a*sx*sy*(xx - xc))/r)*(yy - yc))/(2*r2**1.5) + ((2*yy - 2*yc)*((a*cx*cy*(yy - yc))/r + (a*sx*sy*(xx - xc))/r)*(xx - xc))/(2*r2**1.5)) + cx*sy

        CASE (2)
          ! Axisimmetric case with div(b)~=0
          IF (.not. switch%axisym) THEN
            WRITE (6,*) "This is an axisymmetric test case!"
            stop
          END IF
          a = 2*pi
          b = 2*pi
          !                                                                                                                xc = 0.
          !                                                                                                                yc = 0.
          D = phys%diff_n
          mu = phys%diff_u
          k = phys%a

          cx = cos(a*xx)
          cy = cos(a*yy)
          sx = sin(a*xx)
          sy = sin(a*yy)
          cx2 = cx**2
          cy2 = cy**2

          t1 = cos(a*xx)*cos(a*yy)
          t2 = sin(a*xx)*sin(a*yy)
          t3 = cos(a*xx)*sin(a*yy)
          t4 = cos(a*yy)*sin(a*xx)
          t5 = ((xm**2 - 2*ym*yy - 2*xm*xx + 2*xx**2 + ym**2 + yy**2)/xx**2)**(0.5)
          t6 = ((xm**2 - 2*ym*yy - 2*xm*xx + 2*xx**2 + ym**2 + yy**2)/xx**2)**(1.5)

          f(ind,1) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
            &**2)**2*(-D*a*t3*xm**4 - D*a*t3*xx**4*6.0D0 + a*t3*xx**5*(1.0D0/xx**2*&
            &(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - t1*x&
            &x*ym**3*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*&
            &2.0D0 - t1*xx**3*ym*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)*&
            &*2 + 1.0D0)*2.0D0 + t1*xx*yy**3*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**&
            &2*(ym - yy)**2 + 1.0D0)*2.0D0 + t1*xx**3*yy*sqrt(1.0D0/xx**2*(xm - xx)**2 +&
            &1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 + D*a**2*t2*xx**5*6.0D0 + a*cx*t2*&
            &*2*xx**5*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3&
            &.0D0/2.0D0) - D*a*t3*xm**2*xx**2*1.1D+1 - D*a**2*t2*xm*xx**4*1.0D+1 + D*&
            &a**2*t2*xm**4*xx - D*a*t3*xm**2*ym**2 - D*a*t3*xx**2*ym**2*5.0D0 + D*a**&
            &2*t1*xx**4*ym*4.0D0 + D*a**2*t2*xx*ym**4 - D*a*t3*xm**2*yy**2 - D*a*t3*x&
            &x**2*yy**2*5.0D0 - D*a**2*t1*xx**4*yy*4.0D0 + D*a**2*t2*xx*yy**4 - a*sx*&
            &t1**2*xx**5*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*&
            &*(3.0D0/2.0D0) + D*a**2*t2*xm**2*xx**3*9.0D0 - D*a**2*t2*xm**3*xx**2*4&
            &.0D0 + D*a**2*t1*xx**2*ym**3*2.0D0 + D*a**2*t2*xx**3*ym**2*5.0D0 - D*a**&
            &2*t1*xx**2*yy**3*2.0D0 + D*a**2*t2*xx**3*yy**2*5.0D0 + D*a*t3*xm*xx**3&
            &*1.2D+1 + D*a*t3*xm**3*xx*5.0D0 - D*a*t4*xm*ym**3 - D*a*t4*xm**3*ym + D*a*&
            &t4*xx*ym**3*2.0D0 + D*a*t4*xx**3*ym*2.0D0 + D*a*t4*xm*yy**3 + D*a*t4*xm*&
            &*3*yy - D*a*t4*xx*yy**3*2.0D0 - D*a*t4*xx**3*yy*2.0D0 - a*t3*xm*xx**4*(1&
            &.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)&
            &*2.0D0 + a*t4*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**&
            &2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - a*t4*xx**4*yy*(1.0D0/xx**2*(xm - xx)**&
            &2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)*2.0D0 - t1*t2*xx*ym**&
            &3*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t1*t2*&
            &xx**3*ym*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)&
            &+ t1*t2*xx*yy**3*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
            &+ 1.0D0) + t1*t2*xx**3*yy*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym&
            &- yy)**2 + 1.0D0) + t1*xm*xx**2*ym*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx&
            &**2*(ym - yy)**2 + 1.0D0)*4.0D0 - t1*xm**2*xx*ym*sqrt(1.0D0/xx**2*(xm - xx&
            &)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 - t1*xm*xx**2*yy*sqrt(1.0D0&
            &/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*4.0D0 + t1*xm**2*xx*&
            &yy*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0&
            &- t1*xx*ym*yy**2*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
            &+ 1.0D0)*6.0D0 + t1*xx*ym**2*yy*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx*&
            &*2*(ym - yy)**2 + 1.0D0)*6.0D0 + D*a**2*t1*xm**2*xx**2*ym*6.0D0 - D*a**2*t&
            &2*xm*xx**2*ym**2*4.0D0 + D*a**2*t2*xm**2*xx*ym**2*2.0D0 - D*a**2*t1*xm&
            &**2*xx**2*yy*6.0D0 - D*a**2*t2*xm*xx**2*yy**2*4.0D0 + D*a**2*t2*xm**2*&
            &xx*yy**2*2.0D0 + D*a**2*t1*xx**2*ym*yy**2*6.0D0 - D*a**2*t1*xx**2*ym**&
            &2*yy*6.0D0 + D*a**2*t2*xx*ym**2*yy**2*6.0D0 + D*a*t3*xm*xx*ym**2*3.0D0&
            &- D*a*t4*xm*xx**2*ym*4.0D0 + D*a*t4*xm**2*xx*ym*4.0D0 + D*a*t3*xm*xx*yy&
            &**2*3.0D0 + D*a*t4*xm*xx**2*yy*4.0D0 - D*a*t4*xm**2*xx*yy*4.0D0 + D*a*t3&
            &*xm**2*ym*yy*2.0D0 - D*a*t4*xm*ym*yy**2*3.0D0 + D*a*t4*xm*ym**2*yy*3.0&
            &D0 + D*a*t3*xx**2*ym*yy*1.0D+1 + D*a*t4*xx*ym*yy**2*6.0D0 - D*a*t4*xx*ym&
            &**2*yy*6.0D0 + t1*t2*xm*xx**2*ym*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/x&
            &x**2*(ym - yy)**2 + 1.0D0)*2.0D0 - t1*t2*xm**2*xx*ym*sqrt(1.0D0/xx**2*(x&
            &m - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t1*t2*xm*xx**2*yy*sqrt(1.0D&
            &0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)*2.0D0 + t1*t2*xm**2&
            &*xx*yy*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0) - t&
            &1*t2*xx*ym*yy**2*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**&
            &2 + 1.0D0)*3.0D0 + t1*t2*xx*ym**2*yy*sqrt(1.0D0/xx**2*(xm - xx)**2 + 1.0D0&
            &/xx**2*(ym - yy)**2 + 1.0D0)*3.0D0 - a*cx*t2**2*xm*xx**4*(1.0D0/xx**2*(x&
            &m - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) - D*a**2*t1*xm&
            &*xx*ym**3*2.0D0 - D*a**2*t1*xm*xx**3*ym*8.0D0 - D*a**2*t1*xm**3*xx*ym*&
            &2.0D0 + D*a**2*t1*xm*xx*yy**3*2.0D0 + D*a**2*t1*xm*xx**3*yy*8.0D0 + D*a*&
            &*2*t1*xm**3*xx*yy*2.0D0 - D*a**2*t2*xx*ym*yy**3*4.0D0 - D*a**2*t2*xx*y&
            &m**3*yy*4.0D0 - D*a**2*t2*xx**3*ym*yy*1.0D+1 + a*sx*t1**2*xm*xx**4*(1.&
            &0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) +&
            &a*sy*t4**2*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2&
            &+ 1.0D0)**(3.0D0/2.0D0) - a*sy*t4**2*xx**4*yy*(1.0D0/xx**2*(xm - xx)**2&
            &+ 1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0) - D*a**2*t1*xm*xx*ym*y&
            &y**2*6.0D0 + D*a**2*t1*xm*xx*ym**2*yy*6.0D0 + D*a**2*t2*xm*xx**2*ym*yy&
            &*8.0D0 - D*a**2*t2*xm**2*xx*ym*yy*4.0D0 - D*a*t3*xm*xx*ym*yy*6.0D0 - a*c&
            &x2*cy*sy*xx**4*ym*(1.0D0/xx**2*(xm - xx)**2 + 1.0D0/xx**2*(ym - yy)**2 + 1&
            &.0D0)**(3.0D0/2.0D0) + a*cx2*cy*sy*xx**4*yy*(1.0D0/xx**2*(xm - xx)**2 +&
            &1.0D0/xx**2*(ym - yy)**2 + 1.0D0)**(3.0D0/2.0D0)))/xx

          f(ind,2) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
            &**2)**2*(a**2*mu*xx**2*ym**3*2.0D0 - a**2*mu*xx**2*yy**3*2.0D0 + a*mu*&
            &t4*xm**4*2.0D0 + a*mu*t4*xx**4*1.2D+1 - a*t4*t6*xx**5*1.0D+1 + (a*mu*xm*&
            &*4*sin(a*yy*2.0D0))/2.0D0 + a*mu*xx**4*sin(a*yy*2.0D0)*3.0D0 + a**2*mu&
            &*t1*xx**5*1.2D+1 + a**2*mu*xx**4*ym*4.0D0 - a**2*mu*xx**4*yy*4.0D0 + (a*&
            &mu*xm*ym**3*sin(a*xx*2.0D0))/2.0D0 + (a*mu*xm**3*ym*sin(a*xx*2.0D0))&
            &/2.0D0 - a*mu*xx*ym**3*sin(a*xx*2.0D0) - a*mu*xx**3*ym*sin(a*xx*2.0D0)&
            &- a*mu*xm*xx**3*sin(a*yy*2.0D0)*6.0D0 - (a*mu*xm*yy**3*sin(a*xx*2.0D0&
            &))/2.0D0 - a*mu*xm**3*xx*sin(a*yy*2.0D0)*(5.0D0/2.0D0) - (a*mu*xm**3*y&
            &y*sin(a*xx*2.0D0))/2.0D0 + a*mu*xx*yy**3*sin(a*xx*2.0D0) + a*mu*xx**3*&
            &yy*sin(a*xx*2.0D0) + a**2*mu*xx**2*yy**3*cos(a*xx)**2*4.0D0 - a**2*cx2&
            &*mu*xx**4*ym*8.0D0 - a**2*cy2*mu*xx**4*ym*8.0D0 + a**2*cx2*mu*xx**4*yy&
            &*8.0D0 + a**2*cy2*mu*xx**4*yy*8.0D0 + a**2*mu*t1*t2*xx**5*2.4D+1 - a**2*&
            &mu*t1*xm*xx**4*2.0D+1 + a**2*mu*t1*xm**4*xx*2.0D0 + a*mu*t4*xm**2*xx**&
            &2*2.2D+1 + a*mu*t4*xm**2*ym**2*2.0D0 + a**2*mu*t1*xx*ym**4*2.0D0 + a*mu*&
            &t4*xx**2*ym**2*1.0D+1 + a**2*mu*t2*xx**4*ym*8.0D0 + a*mu*t4*xm**2*yy**&
            &2*2.0D0 + a**2*mu*t1*xx*yy**4*2.0D0 + a*mu*t4*xx**2*yy**2*1.0D+1 - a**2*&
            &mu*t2*xx**4*yy*8.0D0 - a**2*mu*xm*xx*ym**3*2.0D0 - a**2*mu*xm*xx**3*ym&
            &*8.0D0 - a**2*mu*xm**3*xx*ym*2.0D0 + a**2*mu*xm*xx*yy**3*2.0D0 + a**2*mu&
            &*xm*xx**3*yy*8.0D0 + a**2*mu*xm**3*xx*yy*2.0D0 + cx2*t5*xx*yy**3*cos(a&
            &*yy)**2*2.0D0 + a*mu*xm**2*xx**2*sin(a*yy*2.0D0)*(1.1D+1/2.0D0) + (a*m&
            &u*xm**2*ym**2*sin(a*yy*2.0D0))/2.0D0 + a*mu*xx**2*ym**2*sin(a*yy*2.0&
            &D0)*(5.0D0/2.0D0) + (a*mu*xm**2*yy**2*sin(a*yy*2.0D0))/2.0D0 + a*mu*xx&
            &**2*yy**2*sin(a*yy*2.0D0)*(5.0D0/2.0D0) - a**2*cx2*mu*xx**2*ym**3*4.&
            &0D0 - a**2*cy2*mu*xx**2*ym**3*4.0D0 + a**2*cy2*mu*xx**2*yy**3*4.0D0 + a*&
            &*2*mu*t1*xm**2*xx**3*1.8D+1 - a**2*mu*t1*xm**3*xx**2*8.0D0 + a**2*mu*t&
            &1*xx**3*ym**2*1.0D+1 + a**2*mu*t2*xx**2*ym**3*4.0D0 + a**2*mu*t1*xx**3&
            &*yy**2*1.0D+1 - a**2*mu*t2*xx**2*yy**3*4.0D0 + a**2*mu*xm**2*xx**2*ym*&
            &6.0D0 - a**2*mu*xm**2*xx**2*yy*6.0D0 + a**2*mu*xx**2*ym*yy**2*6.0D0 - a*&
            &*2*mu*xx**2*ym**2*yy*6.0D0 + a*cx2*t4*xx**5*(1.0D0/xx**2*(xm*xx*(-2.&
            &0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*2.&
            &0D0 - cx2*cy2*xx**3*ym*sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 +&
            &xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*2.0D0 + cx2*cy2*xx**3*yy*sqrt(1.0D0/&
            &xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*&
            &2.0D0 - a*t3*xx**4*ym*(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2&
            &+ xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*1.0D+1 + a*t3*xx**4*yy*(1.&
            &0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**&
            &2))**(3.0D0/2.0D0)*1.0D+1 - cx2*cy2*t5*xx*ym**3*2.0D0 - a*mu*t4*xm*xx*&
            &*3*2.4D+1 - a*mu*t4*xm**3*xx*1.0D+1 + a*mu*t3*xm*ym**3*2.0D0 + a*mu*t3*x&
            &m**3*ym*2.0D0 - a*mu*t3*xx*ym**3*4.0D0 - a*mu*t3*xx**3*ym*4.0D0 - a*mu*t&
            &3*xm*yy**3*2.0D0 - a*mu*t3*xm**3*yy*2.0D0 + a*mu*t3*xx*yy**3*4.0D0 + a*m&
            &u*t3*xx**3*yy*4.0D0 + a*t4*t6*xm*xx**4*1.0D+1 + a**2*cx**2*cy2*mu*xx**&
            &4*ym*1.6D+1 + a**2*cx2*cy2*mu*xx**2*ym**3*8.0D0 - a**2*cx2*cy2*mu*xx**&
            &2*yy**3*8.0D0 - a*t1**2*t6*xx**4*yy*sin(a*xx)*4.0D0 - a**2*cx2*mu*xm**&
            &2*xx**2*ym*1.2D+1 - a**2*cy2*mu*xm**2*xx**2*ym*1.2D+1 + a**2*cx2*mu*xm&
            &**2*xx**2*yy*1.2D+1 + a**2*cy2*mu*xm**2*xx**2*yy*1.2D+1 - a**2*cx2*mu*&
            &xx**2*ym*yy**2*1.2D+1 + a**2*cx2*mu*xx**2*ym**2*yy*1.2D+1 - a**2*cy2*m&
            &u*xx**2*ym*yy**2*1.2D+1 + a**2*cy2*mu*xx**2*ym**2*yy*1.2D+1 + a**2*mu*&
            &t1*t2*xm**2*xx**3*3.6D+1 - a**2*mu*t1*t2*xm**3*xx**2*1.6D+1 + a**2*mu*&
            &t1*t2*xx**3*ym**2*2.0D+1 + a**2*mu*t1*t2*xx**3*yy**2*2.0D+1 - a**2*mu*&
            &t1*xm*xx**2*ym**2*8.0D0 + a**2*mu*t1*xm**2*xx*ym**2*4.0D0 + a**2*mu*t2&
            &*xm**2*xx**2*ym*1.2D+1 - a**2*mu*t1*xm*xx**2*yy**2*8.0D0 + a**2*mu*t1*&
            &xm**2*xx*yy**2*4.0D0 - a**2*mu*t2*xm**2*xx**2*yy*1.2D+1 + a**2*mu*t1*x&
            &x*ym**2*yy**2*1.2D+1 + a**2*mu*t2*xx**2*ym*yy**2*1.2D+1 - a**2*mu*t2*x&
            &x**2*ym**2*yy*1.2D+1 + a*cx2*cy*sy*xx**5*(1.0D0/xx**2*(xm*xx*(-2.0D0&
            &) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*4.0D0&
            &- a*cx2*cy*mu*sy*xx**4*1.2D+1 - cx2*cy2*t2*xx**3*ym*sqrt(1.0D0/xx**2*&
            &(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + cx2*cy&
            &2*t2*xx*yy**3*sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + x&
            &x**2*2.0D0 + ym**2 + yy**2)) + a**2*mu*xm*xx*ym**3*cos(a*yy)**2*4.0D0 + a*&
            &*2*mu*xm**3*xx*ym*cos(a*yy)**2*4.0D0 - a*cx2*t4*t6*xm*xx**4*2.0D0 - cx&
            &2*cy2*t5*xm**2*xx*ym*2.0D0 - cx2*cy2*t5*xm*xx**2*yy*4.0D0 + cx2*cy2*t5&
            &*xm**2*xx*yy*2.0D0 - cx2*cy2*t5*xx*ym*yy**2*6.0D0 + a*mu*t3*xm*xx**2*y&
            &m*8.0D0 - a*mu*t3*xm**2*xx*ym*8.0D0 - a*mu*t4*xm*xx*ym**2*6.0D0 - a*mu*t&
            &3*xm*xx**2*yy*8.0D0 + a*mu*t3*xm**2*xx*yy*8.0D0 - a*mu*t4*xm*xx*yy**2*&
            &6.0D0 + a*mu*t3*xm*ym*yy**2*6.0D0 - a*mu*t3*xm*ym**2*yy*6.0D0 - a*mu*t4*&
            &xm**2*ym*yy*4.0D0 - a*mu*t3*xx*ym*yy**2*1.2D+1 + a*mu*t3*xx*ym**2*yy*1&
            &.2D+1 - a*mu*t4*xx**2*ym*yy*2.0D+1 + a*mu*xm*xx**2*ym*sin(a*xx*2.0D0)*&
            &2.0D0 - a*mu*xm**2*xx*ym*sin(a*xx*2.0D0)*2.0D0 - a*mu*xm*xx**2*yy*sin(&
            &a*xx*2.0D0)*2.0D0 + a*mu*xm**2*xx*yy*sin(a*xx*2.0D0)*2.0D0 - a*mu*xm*x&
            &x*ym**2*sin(a*yy*2.0D0)*(3.0D0/2.0D0) + a*mu*xm*ym*yy**2*sin(a*xx*2.&
            &0D0)*(3.0D0/2.0D0) - a*mu*xm*ym**2*yy*sin(a*xx*2.0D0)*(3.0D0/2.0D0) -&
            &a*mu*xx*ym*yy**2*sin(a*xx*2.0D0)*3.0D0 + a*mu*xx*ym**2*yy*sin(a*xx*2&
            &.0D0)*3.0D0 - a*mu*xm*xx*yy**2*sin(a*yy*2.0D0)*(3.0D0/2.0D0) - a**2*cx&
            &2*cy2*mu*xx**4*yy*1.6D+1 - a*mu*xm**2*ym*yy*sin(a*yy*2.0D0) - a*mu*xx*&
            &*2*ym*yy*sin(a*yy*2.0D0)*5.0D0 - a*cx2*cy**3*sx*t6*xx**5*3.0D0 + a**2*&
            &cx2*mu*xm*xx*ym**3*4.0D0 + a**2*cx2*mu*xm*xx**3*ym*1.6D+1 + a**2*cx2*m&
            &u*xm**3*xx*ym*4.0D0 + a**2*cy2*mu*xm*xx**3*ym*1.6D+1 - a**2*cx2*mu*xm*&
            &xx*yy**3*4.0D0 - a**2*cx2*mu*xm*xx**3*yy*1.6D+1 - a**2*cx2*mu*xm**3*xx&
            &*yy*4.0D0 - a**2*cy2*mu*xm*xx*yy**3*4.0D0 - a**2*cy2*mu*xm*xx**3*yy*1.&
            &6D+1 - a**2*cy2*mu*xm**3*xx*yy*4.0D0 - a**2*mu*t1*t2*xm*xx**4*4.0D+1 + a&
            &**2*mu*t1*t2*xm**4*xx*4.0D0 - a*mu*sx*t1**2*xm*ym**3*2.0D0 - a*mu*sx*t&
            &1**2*xm**3*ym*2.0D0 + a*mu*sx*t1**2*xx**3*ym*4.0D0 + a*mu*sx*t1**2*xm*&
            &yy**3*2.0D0 + a*mu*sx*t1**2*xm**3*yy*2.0D0 - a*mu*sx*t1**2*xx**3*yy*4.&
            &0D0 + a**2*mu*t1*t2*xx*ym**4*4.0D0 + a**2*mu*t1*t2*xx*yy**4*4.0D0 - a**2&
            &*mu*t2*xm*xx*ym**3*4.0D0 - a**2*mu*t2*xm*xx**3*ym*1.6D+1 - a**2*mu*t2*&
            &xm**3*xx*ym*4.0D0 + a**2*mu*t2*xm*xx*yy**3*4.0D0 + a**2*mu*t2*xm*xx**3&
            &*yy*1.6D+1 + a**2*mu*t2*xm**3*xx*yy*4.0D0 - a**2*mu*t1*xx*ym*yy**3*8.0&
            &D0 - a**2*mu*t1*xx*ym**3*yy*8.0D0 - a**2*mu*t1*xx**3*ym*yy*2.0D+1 + a*sx&
            &*t1**2*t6*xx**4*ym*4.0D0 + a*sy*t1**2*t6*xx**4*ym*2.0D0 - a*sy*t1**2*t&
            &6*xx**4*yy*2.0D0 - a**2*mu*xm*xx*ym*yy**2*6.0D0 + a**2*mu*xm*xx*ym**2*&
            &yy*6.0D0 - a*cy*mu*sy*xm**4*cos(a*xx)**2*2.0D0 - cy2*t2*xx*ym**3*cos(a&
            &*xx)**2*sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2&
            &.0D0 + ym**2 + yy**2)) + cy2*t2*xx**3*yy*cos(a*xx)**2*sqrt(1.0D0/xx**2*(&
            &xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + cy2*t5*&
            &xm*xx**2*ym*cos(a*xx)**2*4.0D0 + cx2*t5*xx*ym**2*yy*cos(a*yy)**2*6.0&
            &D0 + a*mu*t1**2*xx*ym**3*sin(a*xx)*4.0D0 - a*mu*t1**2*xx*yy**3*sin(a*x&
            &x)*4.0D0 + a**2*cy2*mu*xx**2*ym*yy**2*cos(a*xx)**2*2.4D+1 - a*cx2*cy*m&
            &u*sy*xm**2*xx**2*2.2D+1 - a*cx2*cy*mu*sy*xm**2*ym**2*2.0D0 - a*cx2*cy*&
            &mu*sy*xx**2*ym**2*1.0D+1 - a*cx2*cy*mu*sy*xm**2*yy**2*2.0D0 - a**2*cx2&
            &*cy2*mu*xm*xx*ym**3*8.0D0 - a**2*cx2*cy2*mu*xm**3*xx*ym*8.0D0 + a**2*c&
            &x2*cy2*mu*xm*xx*yy**3*8.0D0 + a**2*cx2*cy2*mu*xm*xx**3*yy*3.2D+1 + a**&
            &2*cx2*cy2*mu*xm**3*xx*yy*8.0D0 + a*cx**3*cy2*sy*t6*xx**4*yy*3.0D0 - a*&
            &*2*cx2*mu*xm*xx*ym**2*yy*1.2D+1 - a**2*cy2*mu*xm*xx*ym**2*yy*1.2D+1 -&
            &a*cy2*sy*xx**4*ym*cos(a*xx)**3*(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*&
            &2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*3.0D0 - a*mu*sx&
            &*t1**2*xm*xx**2*ym*8.0D0 + a*mu*sx*t1**2*xm*xx**2*yy*8.0D0 - a*mu*sx*t&
            &1**2*xm*ym*yy**2*6.0D0 + a*mu*sx*t1**2*xm*ym**2*yy*6.0D0 + a*mu*sx*t1*&
            &*2*xx*ym*yy**2*1.2D+1 - a*mu*sx*t1**2*xx*ym**2*yy*1.2D+1 - a**2*mu*t1*&
            &t2*xx*ym*yy**3*1.6D+1 - a**2*mu*t1*t2*xx*ym**3*yy*1.6D+1 - a**2*mu*t1*&
            &t2*xx**3*ym*yy*4.0D+1 + a**2*mu*t1*xm*xx**2*ym*yy*1.6D+1 - a**2*mu*t1*&
            &xm**2*xx*ym*yy*8.0D0 - a**2*mu*t2*xm*xx*ym*yy**2*1.2D+1 + a**2*mu*t2*x&
            &m*xx*ym**2*yy*1.2D+1 - cx2*t2*xm**2*xx*ym*cos(a*yy)**2*sqrt(1.0D0/xx&
            &**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) - cy&
            &2*t2*xx*ym*yy**2*cos(a*xx)**2*sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*&
            &yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*3.0D0 + a*mu*t1**2*xm**2*xx&
            &*ym*sin(a*xx)*8.0D0 - a*mu*t1**2*xm**2*xx*yy*sin(a*xx)*8.0D0 + a**2*cx&
            &2*cy2*mu*xm**2*xx**2*ym*2.4D+1 - a**2*cx2*cy2*mu*xm**2*xx**2*yy*2.4D&
            &+1 - a**2*cx2*cy2*mu*xx**2*ym**2*yy*2.4D+1 + a*mu*t4*xm*xx*ym*yy*1.2D+&
            &1 + a**2*cx**2*mu*xm*xx*ym*yy**2*1.2D+1 + a**2*cy**2*mu*xm*xx*ym*yy**2&
            &*1.2D+1 - a**2*mu*t1*t2*xm*xx**2*ym**2*1.6D+1 + a**2*mu*t1*t2*xm**2*xx&
            &*ym**2*8.0D0 - a**2*mu*t1*t2*xm*xx**2*yy**2*1.6D+1 + a**2*mu*t1*t2*xm*&
            &*2*xx*yy**2*8.0D0 + a**2*mu*t1*t2*xx*ym**2*yy**2*2.4D+1 - a*cy*mu*sy*x&
            &x**2*yy**2*cos(a*xx)**2*1.0D+1 - a**2*cy2*mu*xm*xx**3*ym*cos(a*xx)**&
            &2*3.2D+1 + a*cy**3*sx*t6*xm*xx**4*cos(a*xx)**2*3.0D0 + a*cx2*cy*mu*sy*&
            &xm*xx**3*2.4D+1 + cx2*cy2*t2*xm*xx**2*ym*sqrt(1.0D0/xx**2*(xm*xx*(-2&
            &.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*2.0D0 - cx2*cy2*t2&
            &*xm*xx**2*yy*sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx&
            &**2*2.0D0 + ym**2 + yy**2))*2.0D0 + cx2*cy2*t2*xm**2*xx*yy*sqrt(1.0D0/xx&
            &**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + a*&
            &mu*xm*xx*ym*yy*sin(a*yy*2.0D0)*3.0D0 + cx2*cy2*t2*xx*ym**2*yy*sqrt(1&
            &.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy*&
            &*2))*3.0D0 - a*cx2*cy*sy*t6*xm*xx**4*4.0D0 + a*cx2*mu*sy*xm**3*xx*cos(&
            &a*yy)*1.0D+1 - a**2*cx2*mu*xm*xx*ym*yy**2*cos(a*yy)**2*2.4D+1 + a*cx2*&
            &cy*mu*sy*xm*xx*ym**2*6.0D0 + a*cx2*cy*mu*sy*xx**2*ym*yy*2.0D+1 + a*cx2&
            &*cy*mu*xm**2*ym*yy*sin(a*yy)*4.0D0 + a**2*cx2*cy2*mu*xm*xx*ym**2*yy*&
            &2.4D+1 + a**2*mu*t1*t2*xm*xx**2*ym*yy*3.2D+1 - a**2*mu*t1*t2*xm**2*xx*&
            &ym*yy*1.6D+1 + a*cy*mu*sy*xm*xx*yy**2*cos(a*xx)**2*6.0D0 - a*cy*mu*sy*&
            &xm*xx*ym*yy*cos(a*xx)**2*1.2D+1))/xx

        CASE (50:)
          !Do nothing
        CASE DEFAULT
          WRITE (6,*) "Error! Test case not valid"
          STOP

        END SELECT
      END DO
    END DO
  END SUBROUTINE body_force
#else
  !***********************************************************************
  !
  !                            VERSION 2D
  !
  !***********************************************************************
  !*****************************************
  ! Analytical solution
  !****************************************
  SUBROUTINE analytical_solution(x,y,u)
    use mpi_omp
    real*8,dimension(:),intent(IN)        :: x,y
    real*8,dimension(:,:),intent(OUT)     :: u
    real*8,dimension(size(u,1),size(u,2))  :: up
    real*8 :: a,b,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
    real*8 :: dsource(size(x)),aux(size(x)),xsource,ysource,smod
    integer:: i,np
    real*8 :: r(size(x)),rs,umin

    np = size(x)
    u = 0.
    up = 0.
    a = 2*pi
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    SELECT CASE (switch%testcase)
    CASE (1)
      IF (switch%axisym) THEN
        WRITE (6,*) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Circular field centered in [xc,yc],n = 2+sin(a*x)*sin(a*y), u = cos(a*x)*cos(a*y)
      ! Case 9 of the Matlab version: for convergence purpose
      up(:,1) = 2 + sin(a*x)*sin(a*y)
      up(:,2) = cos(a*x)*cos(a*y)
      up(:,3) = cos(a*x)*sin(a*y)
      up(:,4) = sin(a*x)*cos(a*y)

    CASE (5)
      IF (switch%axisym) THEN
        WRITE (6,*) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Cartesian case,circular field centered in [xm,ym] in the poloidal plane,Bt = 1
      smod = 1.
      rs = 0.04/simpar%refval_length
      xsource = xm-0.5*(xmax-xm)
      ysource = ym
      dsource   = sqrt((x-xsource)**2+(y-ysource)**2)
      aux = -dsource**2/rs**2
      up(:,1) = 1e-3
      DO i=1,np
        if (aux(i).gt.-30) then
          up(i,1) =  up(i,1)+smod*exp(aux(i))
        endif
      END DO
    CASE (6)
      IF (.not.switch%axisym) THEN
        WRITE (6,*) "This is an axisymmetric test case!"
        stop
      END IF
      ! Cartesian case,circular field centered in [xm,ym] in the poloidal plane,Bt = 1
      smod = 1.
      rs = 0.04/simpar%refval_length
      xsource = xm-0.5*(xmax-xm)
      ysource = ym
      dsource   = sqrt((x-xsource)**2+(y-ysource)**2)
      aux = -dsource**2/rs**2
      up(:,1) = 1e-3
      DO i=1,np
        if (aux(i).gt.-30) then
          up(i,1) =  up(i,1)+smod*exp(aux(i))
        endif
      END DO
    CASE (7)
      IF (.not.switch%axisym) THEN
        WRITE (6,*) "This is an axisymmetric test case!"
        stop
      END IF
      ! Turbulence interchange case
      umin = 1e-2
      xsource = xm-0.5*(xmax-xm)
      rs = 0.01/simpar%refval_length
      dsource = -(x-xsource)/rs
      aux = 0.5*(1+tanh(dsource))

      up(:,1) = aux+umin*(1-aux)
    CASE (50:64)
      up(:,1) = 1.
      up(:,2) = 0.
    CASE (65)
      up(:,1) = 1.
      up(:,2) = 0.
      r = sqrt((x*phys%lscale - geom%R0)**2 + (y*phys%lscale - 0.75)**2)
      DO i = 1,size(x)
        IF (r(i) .le. 0.05) THEN
          up(i,2) = 1.
        END IF
      END DO
    CASE DEFAULT
      WRITE (6,*) "Error! Test case not valid"
      STOP
    END SELECT
    ! Convert physical variables to conservative variables
    CALL phys2cons(up,u)
  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x,y,u,ux,uy)
    real*8,dimension(:),intent(IN)        :: x,y
    real*8,dimension(:,:),intent(IN)      :: u
    real*8,dimension(:,:),intent(OUT)     :: ux,uy
    real*8,dimension(size(u,1),size(u,2))  :: upx,upy
    real*8,dimension(size(u,1),phys%npv)  :: up
    real*8 :: a

    upx = 0.
    upy = 0.
    CALL cons2phys(u,up)
    a = 2*pi
    SELECT CASE (switch%testcase)
    CASE (1)
      IF (switch%axisym) THEN
        WRITE (6,*) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Circular field centered in [xc,yc],n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y),Ei = 20+cos(wx*x)*sin(wy*y),Ee = 10-sin(wx*x)*cos(wy*y)
      upx(:,1) = a*cos(a*x)*sin(a*y)
      upx(:,2) = -a*sin(a*x)*cos(a*y)
      upx(:,3) = -a*sin(a*x)*sin(a*y)
      upx(:,4) = a*cos(a*x)*cos(a*y)

      upy(:,1) = a*sin(a*x)*cos(a*y)
      upy(:,2) = -a*cos(a*x)*sin(a*y)
      upy(:,3) = a*cos(a*x)*cos(a*y)
      upy(:,4) = -a*sin(a*x)*sin(a*y)

    CASE (2)
      ! Axisimmetric case with div(b)~=0,n = 2+sin(wx*x )*sin(wy*y), u = cos(wx*x)*cos(wy*y),Ei = 20+cos(wx*x)*sin(wy*y),Ee = 10-sin(wx*x)*cos(wy*y)
      IF (.not. switch%axisym) THEN
        WRITE (6,*) "This is an axisymmetric test case!"
        stop
      END IF
      upx(:,1) = a*cos(a*x)*sin(a*y)
      upx(:,2) = -a*sin(a*x)*cos(a*y)
      upx(:,3) = -a*sin(a*x)*sin(a*y)
      upx(:,4) = -a*cos(a*x)*cos(a*y)

      upy(:,1) = a*sin(a*x)*cos(a*y)
      upy(:,2) = -a*cos(a*x)*sin(a*y)
      upy(:,3) = a*cos(a*x)*cos(a*y)
      upy(:,4) = a*sin(a*x)*sin(a*y)
    CASE (5:7)
      ! Do nothing
    CASE (50:)
      ! Do nothing
    CASE DEFAULT
      WRITE (6,*) "Error! Test case not valid"
      STOP
    END SELECT
    ! Convert physical variables to conservative variables
    ux(:,1) = upx(:,1)
    uy(:,1) = upy(:,1)
    ux(:,2) = (upx(:,1)*up(:,2) + up(:,1)*upx(:,2))
    uy(:,2) = (upy(:,1)*up(:,2) + up(:,1)*upy(:,2))
    ux(:,3) = upx(:,3)
    uy(:,3) = upy(:,3)
    ux(:,4) = upx(:,4)
    uy(:,4) = upy(:,4)

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x,y,f)
    real*8,dimension(:),intent(IN) :: x,y
    real*8,dimension(:,:),intent(OUT) :: f
    integer                            :: i,n
    real*8  :: a,b,xc,yc,D,mu,k,source
    real*8  :: Dw,DparPhi,DPhi,DnDivn,DparN
    real*8 :: r(size(x)),r2(size(x)),cx(size(x)),cy(size(x)),sx(size(x)),sy(size(x)),driftexb

    real*8 :: xmax,xmin,ymax,ymin,xm,ym

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    Dw = phys%diff_vort
    DPhi = phys%diff_pot
    DparPhi = phys%diff_pari
    DparN = phys%diff_pare
    DnDivn = phys%diff_ee

    n = size(x)

    if (switch%driftexb) then
      driftexb = 1.
    else
      driftexb = 0.
    endif

    f = 0.
    SELECT CASE (switch%testcase)
    CASE (1)
      ! Circular field centered in [xc,yc],n = 2+sin(2*pi*x)*sin(2*pi*y), u = cos(2*pi*x)*cos(2*pi*y)
      ! Case 9 of the Matlab version: for convergence purpose
      a = 2*pi
      b = 2*pi
      xc = 0.
      yc = 0.
      D = phys%diff_n
      mu = phys%diff_u
      k = phys%a

      r = ((x - xc)**2 + (y - yc)**2 + 1)**0.5
      r2 = ((x - xc)**2 + (y - yc)**2 + 1)
      cx = cos(a*x)
      cy = cos(a*y)
      sx = sin(a*x)
      sy = sin(a*y)

      f(:,1) =   driftexb*((a**2*cx*cy**2*sx)/r2 + (a**2*cx*sx*sy**2)/r2 - (a*sx*sy*(2*x - 2*xc)*(sx*sy + 2))/r2**2 - (a*cx*cy*(2*y - 2*yc)*(sx*sy + 2))/r2**2) - D*(((y - yc)*((a*cy*sx)/r + (a**2*sx*sy*(y - yc))/r + (a**2*cx*cy*(x - xc))/r - (a*cy*sx*(2*x - 2*xc)*(x - xc))/(2*r2**1.5) + (a*cx*sy*(2*x - 2*xc)*(y - yc))/(2*r2**1.5)))/r + ((x - xc)*((a*cx*sy)/r + (a**2*sx*sy*(x - xc))/r + (a**2*cx*cy*(y - yc))/r + (a*cy*sx*(2*y - 2*yc)*(x - xc))/(2*r2**1.5) - (a*cx*sy*(2*y - 2*yc)*(y - yc))/(2*r2**1.5)))/r - 2*a**2*sx*sy - ((2*x - 2*xc)*((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(y - yc))/(2*r2**1.5) + ((2*y - 2*yc)*((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(x - xc))/(2*r2**1.5)) - ((cx*cy*(2*x - 2*xc)*(sx*sy + 2)*(y - yc))/(2*r2**1.5) - (cx*cy*(2*y - 2*yc)*(sx*sy + 2)*(x - xc))/(2*r2**1.5) + (a*cx*cy**2*sx*(x - xc))/r - (a*cx**2*cy*sy*(y - yc))/r - (a*cx*sy*(sx*sy + 2)*(x - xc))/r + (a*cy*sx*(sx*sy + 2)*(y - yc))/r)
      f(:,2) =   (mu*(2*a**2*x*yc - 2*a**2*x*y + 2*a**2*xc*y - 2*a**2*xc*yc + 4*a**2*cx*cy + (a*x*sin(2*a*y))/2 + (a*y*sin(2*a*x))/2 - (a*xc*sin(2*a*y))/2 - (a*yc*sin(2*a*x))/2 + 4*a**2*x*y*cx**2 - 4*a**2*x*yc*cx**2 - 4*a**2*xc*y*cx**2 + 4*a**2*xc*yc*cx**2 + 4*a**2*x*y*cy**2 - 4*a**2*x*yc*cy**2 - 4*a**2*xc*y*cy**2 + 4*a**2*xc*yc*cy**2 + 2*a**2*x**2*cx*cy + 2*a**2*xc**2*cx*cy + 2*a**2*y**2*cx*cy + 2*a**2*yc**2*cx*cy + 2*a*x*cy*sx - 2*a*xc*cy*sx + 2*a*y*cx*sy - 2*a*yc*cx*sy - 8*a**2*x*y*cx**2*cy**2 + 8*a**2*x*yc*cx**2*cy**2 + 8*a**2*xc*y*cx**2*cy**2 - 8*a**2*xc*yc*cx**2*cy**2 - 4*a**2*x*xc*cx*cy - 4*a**2*y*yc*cx*cy - 4*a**2*x*y*sx*sy + 4*a**2*x*yc*sx*sy + 4*a**2*xc*y*sx*sy - 4*a**2*xc*yc*sx*sy - 2*a*x*cx**2*cy*sy - 2*a*y*cx*cy**2*sx + 2*a*xc*cx**2*cy*sy + 2*a*yc*cx*cy**2*sx + 8*a**2*cx*cy*sx*sy + 4*a**2*x**2*cx*cy*sx*sy + 4*a**2*xc**2*cx*cy*sx*sy + 4*a**2*y**2*cx*cy*sx*sy + 4*a**2*yc**2*cx*cy*sx*sy - 8*a**2*x*xc*cx*cy*sx*sy - 8*a**2*y*yc*cx*cy*sx*sy))/(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1) - driftexb*((a**2*cx**2*cy*sy*(sx*sy + 2))/r2 - (a**2*cx**2*cy*sx*sy**2)/r2 - (a**2*cx**2*cy**3*sx)/r2 + (a**2*cy*sx**2*sy*(sx*sy + 2))/r2 + (a*cx**2*cy**2*(2*y - 2*yc)*(sx*sy + 2))/r2**2 + (a*cx*cy*sx*sy*(2*x - 2*xc)*(sx*sy + 2))/r2**2) + (a*(2*x*cx**2*cy*sx - 2*xc*cx**2*cy*sx + 4*x*cx**2*cy*sy - 4*y*cx*cy**2*sx - 4*xc*cx**2*cy*sy + 4*yc*cx*cy**2*sx - 2*y*cx*cy**2*sy + 2*yc*cx*cy**2*sy - 3*x*cx**2*cy**3*sx + 3*xc*cx**2*cy**3*sx + 3*y*cx**3*cy**2*sy - 3*yc*cx**3*cy**2*sy - k*x*cy*sx + k*xc*cy*sx + k*y*cx*sy - k*yc*cx*sy))/(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1)**0.5
      f(:,3) =   (Dw*a*(2*a*cx*sy - y*cx*cy + yc*cx*cy + x*sx*sy - xc*sx*sy + a*x**2*cx*sy + a*xc**2*cx*sy + a*y**2*cx*sy + a*yc**2*cx*sy - 2*a*x*xc*cx*sy + 2*a*x*y*cy*sx - 2*a*x*yc*cy*sx - 2*a*xc*y*cy*sx + 2*a*xc*yc*cy*sx - 2*a*y*yc*cx*sy))/(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1) - (DparPhi*a*(x*cx*cy - xc*cx*cy - y*sx*sy + yc*sx*sy + a*x**2*cy*sx + a*xc**2*cy*sx + a*y**2*cy*sx + a*yc**2*cy*sx - 2*a*x*xc*cy*sx - 2*a*x*y*cx*sy + 2*a*x*yc*cx*sy + 2*a*xc*y*cx*sy - 2*a*xc*yc*cx*sy - 2*a*y*yc*cy*sx))/(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1) - driftexb*((a**2*sx**2*sy**2)/r2 - (a**2*cx**2*cy**2)/r2 + (a*cx*sx*sy**2*(2*x - 2*xc))/r2**2 + (a*cx**2*cy*sy*(2*y - 2*yc))/r2**2) + (a*cx*(x*cx - xc*cx - 2*x*cx*cy**2 + 2*xc*cx*cy**2 - 2*y*cy*sx*sy + 2*yc*cy*sx*sy))/(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1)**0.5 + (DparN*a*((x*sin(2*a*x))/2 - (xc*sin(2*a*x))/2 + (y*sin(2*a*y))/2 - (yc*sin(2*a*y))/2 + a*x**2 + a*xc**2 + a*y**2 + a*yc**2 + 2*x*cx*sy - 2*xc*cx*sy + 2*y*cy*sx - 2*yc*cy*sx - a*x**2*cx**2 - a*xc**2*cx**2 - a*y**2*cy**2 - a*yc**2*cy**2 - 2*a*x*xc - 2*a*y*yc - x*cx*cy**2*sx + xc*cx*cy**2*sx - y*cx**2*cy*sy + yc*cx**2*cy*sy + 2*a*x*xc*cx**2 + 2*a*y*yc*cy**2 + 2*a*x**2*sx*sy + 2*a*xc**2*sx*sy + 2*a*y**2*sx*sy + 2*a*yc**2*sx*sy - 4*a*x*xc*sx*sy - 4*a*y*yc*sx*sy + 4*a*x*y*cx*cy - 4*a*x*yc*cx*cy - 4*a*xc*y*cx*cy + 4*a*xc*yc*cx*cy))/((sx*sy + 2)**2*(x**2 - 2*y*yc - 2*x*xc + xc**2 + y**2 + yc**2 + 1))
      f(:,4) =   DnDivn*((a**2*sx*sy - ((y - yc)*((a*cy*sx)/r + (a**2*sx*sy*(y - yc))/r + (a**2*cx*cy*(x - xc))/r - (a*cy*sx*(2*x - 2*xc)*(x - xc))/(2*r2**1.5) + (a*cx*sy*(2*x - 2*xc)*(y - yc))/(2*r2**1.5)))/r + ((2*x - 2*xc)*((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(y - yc))/(2*r2**1.5))/(sx*sy + 2) - (((x - xc)*((a*cx*sy)/r + (a**2*sx*sy*(x - xc))/r + (a**2*cx*cy*(y - yc))/r + (a*cy*sx*(2*y - 2*yc)*(x - xc))/(2*r2**1.5) - (a*cx*sy*(2*y - 2*yc)*(y - yc))/(2*r2**1.5)))/r - a**2*sx*sy + ((2*y - 2*yc)*((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(x - xc))/(2*r2**1.5))/(sx*sy + 2) + (a*cy*sx*(a*cy*sx - (((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(x - xc))/r))/(sx*sy + 2)**2 + (a*cx*sy*(a*cx*sy + (((a*cy*sx*(x - xc))/r - (a*cx*sy*(y - yc))/r)*(y - yc))/r))/(sx*sy + 2)**2) + DPhi*(((x - xc)*((a**2*cx*sy*(y - yc))/r - (a**2*cy*sx*(x - xc))/r - (a*cx*cy)/r + (a*cx*cy*(2*y - 2*yc)*(y - yc))/(2*r2**1.5) + (a*sx*sy*(2*y - 2*yc)*(x - xc))/(2*r2**1.5)))/r + 2*a**2*cy*sx - ((y - yc)*((a**2*cy*sx*(y - yc))/r - (a**2*cx*sy*(x - xc))/r - (a*sx*sy)/r + (a*cx*cy*(2*x - 2*xc)*(y - yc))/(2*r2**1.5) + (a*sx*sy*(2*x - 2*xc)*(x - xc))/(2*r2**1.5)))/r - ((2*x - 2*xc)*((a*cx*cy*(y - yc))/r + (a*sx*sy*(x - xc))/r)*(y - yc))/(2*r2**1.5) + ((2*y - 2*yc)*((a*cx*cy*(y - yc))/r + (a*sx*sy*(x - xc))/r)*(x - xc))/(2*r2**1.5)) + cx*sy


    CASE (2)
      ! Axisimmetric case with div(b)~=0
      IF (.not. switch%axisym) THEN
        WRITE (6,*) "This is an axisymmetric test case!"
        stop
      END IF
      a = 2*pi
      b = 2*pi
      xc = 0.
      yc = 0.
      D = phys%diff_n
      mu = phys%diff_u
      k = phys%a
      f(:,1) = 0.
      f(:,2) = 0.
    CASE (5:6)
      ! Do nothing
    CASE (7)
      !         source = 1e-5
      !         do i=1,n
      !            if (x(i).gt.xm-0.5*(xmax-xm)) then
      !               f(i,1) = -source
      !               f(i,3) = -source
      !            endif
      !         end do
    CASE (50:)
      !Do nothing
    CASE DEFAULT
      WRITE (6,*) "Error! Test case not valid"
      STOP

    END SELECT
  END SUBROUTINE body_force
#endif

END MODULE analytical
