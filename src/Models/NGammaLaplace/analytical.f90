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
  USE globals, ONLY: switch
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
  SUBROUTINE analytical_solution(x, y, t, u)
    real*8, dimension(:), intent(IN)        :: x, y, t
    real*8, dimension(:, :), intent(OUT)     :: u
    real*8, dimension(size(u, 1), phys%npv)  :: up
    integer:: i, j, ind, N2D, N1D
    real*8 :: a, b, r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    real*8 :: aux

    u = 0.
    a = 2*pi
    N2D = size(x, 1)
    N1D = size(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    DO i = 1, N2d
      DO j = 1, N1d  ! TODO: check if I need to invert i and j
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          IF (switch%axisym) THEN
            WRITE (6, *) "This is NOT an axisymmetric test case!"
            stop
          END IF
          ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
          up(ind, 1) = 2 + sin(a*xx)*sin(a*yy)
          up(ind, 2) = cos(a*xx)*cos(a*yy)
          up(ind, 3) = up(ind, 1)
        CASE (2)
          ! Axisimmetric case with div(b)~=0
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          up(ind, 1) = 2 + sin(a*xx)*sin(a*yy)
          up(ind, 2) = cos(a*xx)*cos(a*yy)
          up(ind, 3) = up(ind, 1)
        CASE (50:64)
          up(ind, 1) = 1.
          up(ind, 2) = 0.
          up(ind, 3) = up(ind, 1)
        CASE (65)
          up(ind, 1) = 1.
          up(ind, 2) = 0.
          up(ind, 3) = up(ind, 1)
          r = sqrt((xx*phys%lscale - geom%R0)**2 + (yy*phys%lscale - 0.75)**2)
          IF (r .le. 0.05) THEN
            up(ind, 2) = 1.
          END IF
        CASE DEFAULT
          WRITE (6, *) "Error! Test case not valid"
          STOP
        END SELECT
      END DO
    END DO
    ! Convert physical variables to conservative variables
    CALL phys2cons(up, u)
  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x, y, t, u, ux, uy, ut)
    real*8, dimension(:), intent(IN)        :: x, y, t
    real*8, dimension(:, :), intent(IN)      :: u
    real*8, dimension(:, :), intent(OUT)     :: ux, uy, ut
    real*8, dimension(size(u, 1), phys%npv)  :: upx, upy, upt
    real*8, dimension(size(u, 1), phys%npv)   :: up
    integer:: i, j, ind, N1D, N2D
    real*8 :: a, b, r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    real*8 :: aux

    N2D = size(x, 1)
    N1D = size(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)
    upx = 0.
    upy = 0.
    upt = 0.
    CALL cons2phys(u, up)
    a = 2*pi

    DO i = 1, N2d
      DO j = 1, N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          IF (switch%axisym) THEN
            WRITE (6, *) "This is NOT an axisymmetric test case!"
            stop
          END IF
          ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          upx(ind, 1) = a*cos(a*xx)*sin(a*yy)
          upx(ind, 2) = -a*sin(a*xx)*cos(a*yy)

          upy(ind, 1) = a*sin(a*xx)*cos(a*yy)
          upy(ind, 2) = -a*cos(a*xx)*sin(a*yy)

          upx(ind, 3) = upx(ind, 1)
          upy(ind, 3) = upy(ind, 1)
        CASE (2)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          upx(ind, 1) = a*cos(a*xx)*sin(a*yy)
          upx(ind, 2) = -a*sin(a*xx)*cos(a*yy)

          upy(ind, 1) = a*sin(a*xx)*cos(a*yy)
          upy(ind, 2) = -a*cos(a*xx)*sin(a*yy)

          upx(ind, 3) = upx(ind, 1)
          upy(ind, 3) = upy(ind, 1)
        CASE (3)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          upt(ind, 1) = +a*cos(a*tt)
          upt(ind, 2) = -a*sin(a*tt)
          upt(ind, 3) = upt(ind, 1)
        CASE (5)
          ! Do nothing
        CASE (6)
          ! Do nothing
        CASE (50:)
          ! Do nothing
        CASE DEFAULT
          WRITE (6, *) "Error! Test case not valid"
          STOP
        END SELECT
      END DO
    END DO

    ! Convert physical variables to conservative variables
    ux(:, 1) = upx(:, 1)
    uy(:, 1) = upy(:, 1)
    ut(:, 1) = upt(:, 1)
    ux(:, 2) = (upx(:, 1)*up(:, 2) + up(:, 1)*upx(:, 2))
    uy(:, 2) = (upy(:, 1)*up(:, 2) + up(:, 1)*upy(:, 2))
    ut(:, 2) = (upt(:, 1)*up(:, 2) + up(:, 1)*upt(:, 2))
    ux(:, 3) = upx(:, 3)
    uy(:, 3) = upy(:, 3)
    ut(:, 3) = upt(:, 3)
  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, t, f)
    real*8, dimension(:), intent(IN) :: x, y, t
    real*8, dimension(:, :), intent(OUT) :: f
    real*8  :: xc, yc, D, mu, k
    integer:: i, j, N2D, N1D, ind
    real*8 :: a, b, r, r2, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    real*8  :: t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb, cx, cy, sx, sy, cx2, cy2
    f = 0.
    a = 2*pi
    N2D = size(x, 1)
    N1D = size(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    f = 0.
    DO i = 1, N2d
      DO j = 1, N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        r2 = ((xm - xx)**2 + (ym - yy)**2)**(1.5)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)
        CASE (1)
          ! Circular field centered in [xc, yc], n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y)
          IF (switch%axisym) THEN
            WRITE (6, *) "This is NOT an axisymmetric test case!"
            stop
          END IF
          a = 2*pi
          b = 2*pi
          !                                                                                                                xc = 0.
          !                                                                                                                yc = 0.
          !                                                                                                                xm = -0.5
          !                                                                                                                ym = -0.5
          D = phys%diff_n
          mu = phys%diff_u
          k = phys%a

          t1 = cos(a*xx)*cos(a*yy)
          t2 = sin(a*xx)*sin(a*yy)
          t3 = cos(a*xx)*sin(a*yy)
          t4 = cos(a*yy)*sin(a*xx)
          t5 = cos(a*yy)*sin(a*yy)
          t8 = cos(a*xx)**2*cos(a*yy)**2
          ta = (xm**2 - 2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1)

          f(ind,1) = (2*D*a**2*sin(a*xx)*sin(a*yy)-2*a*xm*cos(a*xx)*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*xx*cos(a*xx)*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*ym*cos(a*yy)*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*a*yy*cos(a*yy)*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+D*a**2*xm**2*sin(a*xx)*sin(a*yy)+D*a**2*xx**2*sin(a*xx)*sin(a*yy)+D*a**2*ym**2*sin(a*xx)*sin(a*yy)+D*a**2*yy**2*sin(a*xx)*sin(a*yy)+D*a*xm*cos(a*xx)*sin(a*yy)-D*a*xx*cos(a*xx)*sin(a*yy)+D*a*ym*cos(a*yy)*sin(a*xx)-D*a*yy*cos(a*yy)*sin(a*xx)+a*xm*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xx*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xm*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*xx*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*D*a**2*xm*ym*cos(a*xx)*cos(a*yy)+2*D*a**2*xx*ym*cos(a*xx)*cos(a*yy)+2*D*a**2*xm*yy*cos(a*xx)*cos(a*yy)-2*D*a**2*xx*yy*cos(a*xx)*cos(a*yy)-2*D*a**2*xm*xx*sin(a*xx)*sin(a*yy)-2*D*a**2*ym*yy*sin(a*xx)*sin(a*yy))/((xm-xx)**2+(ym-yy)**2+1);
          f(ind,2) = ((a*mu*xx*sin(2*a*yy))/2-(a*mu*xm*sin(2*a*yy))/2-(a*mu*ym*sin(2*a*xx))/2+(a*mu*yy*sin(2*a*xx))/2-2*a**2*mu*xm*ym+2*a**2*mu*xx*ym+2*a**2*mu*xm*yy-2*a**2*mu*xx*yy+4*a**2*mu*cos(a*xx)*cos(a*yy)+4*a**2*mu*xm*ym*cos(a*xx)**2-4*a**2*mu*xx*ym*cos(a*xx)**2-4*a**2*mu*xm*yy*cos(a*xx)**2+4*a**2*mu*xx*yy*cos(a*xx)**2+4*a**2*mu*xm*ym*cos(a*yy)**2-4*a**2*mu*xx*ym*cos(a*yy)**2-4*a**2*mu*xm*yy*cos(a*yy)**2+4*a**2*mu*xx*yy*cos(a*yy)**2+2*a**2*mu*xm**2*cos(a*xx)*cos(a*yy)+2*a**2*mu*xx**2*cos(a*xx)*cos(a*yy)+2*a**2*mu*ym**2*cos(a*xx)*cos(a*yy)+2*a**2*mu*yy**2*cos(a*xx)*cos(a*yy)-2*a*mu*xm*cos(a*yy)*sin(a*xx)+2*a*mu*xx*cos(a*yy)*sin(a*xx)-2*a*mu*ym*cos(a*xx)*sin(a*yy)+2*a*mu*yy*cos(a*xx)*sin(a*yy)+8*a**2*mu*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)-2*a*xm*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+2*a*xx*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+4*a*ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a*xm*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+4*a*xx*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+2*a*ym*cos(a*xx)*cos(a*yy)**2*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-2*a*yy*cos(a*xx)*cos(a*yy)**2*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-8*a**2*mu*xm*ym*cos(a*xx)**2*cos(a*yy)**2+8*a**2*mu*xx*ym*cos(a*xx)**2*cos(a*yy)**2+8*a**2*mu*xm*yy*cos(a*xx)**2*cos(a*yy)**2-8*a**2*mu*xx*yy*cos(a*xx)**2*cos(a*yy)**2+3*a*xm*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-3*a*xx*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-3*a*ym*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+3*a*yy*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+a*k*xm*cos(a*yy)*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-a*k*xx*cos(a*yy)*sin(a*xx)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-a*k*ym*cos(a*xx)*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)+a*k*yy*cos(a*xx)*sin(a*yy)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)**(0.5)-4*a**2*mu*xm*xx*cos(a*xx)*cos(a*yy)-4*a**2*mu*ym*yy*cos(a*xx)*cos(a*yy)-4*a**2*mu*xm*ym*sin(a*xx)*sin(a*yy)+4*a**2*mu*xx*ym*sin(a*xx)*sin(a*yy)+4*a**2*mu*xm*yy*sin(a*xx)*sin(a*yy)-4*a**2*mu*xx*yy*sin(a*xx)*sin(a*yy)+2*a*mu*ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+2*a*mu*xm*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-2*a*mu*xx*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-2*a*mu*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+4*a**2*mu*xm**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)+4*a**2*mu*xx**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)+4*a**2*mu*ym**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)+4*a**2*mu*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)-8*a**2*mu*xm*xx*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)-8*a**2*mu*ym*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy))/(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1);
          f(ind, 3) = f(ind, 1)
          !                                                                                                                f(ind,1) = (2*D*a**2*t2-2*a*xm*t3*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*xx*t3*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+2*a*ym*t4*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*a*yy*t4*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+D*a**2*xm**2*t2+D*a**2*xx**2*t2+D*a**2*ym**2*t2+D*a**2*yy**2*t2+D*a*xm*t3-D*a*xx*t3+D*a*ym*t4-D*a*yy*t4+a*xm*t1**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xx*t1**2*sin(a*xx)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*ym*cos(a*xx)**2*t5*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*yy*cos(a*xx)**2*t5*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*xm*cos(a*xx)*t2**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*xx*cos(a*xx)*t2**2*((xm-xx)**2+(ym-yy)**2+1)**(0.5)+a*ym*t4**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-a*yy*t4**2*sin(a*yy)*((xm-xx)**2+(ym-yy)**2+1)**(0.5)-2*D*a**2*xm*ym*t1+2*D*a**2*xx*ym*t1+2*D*a**2*xm*yy*t1-2*D*a**2*xx*yy*t1-2*D*a**2*xm*xx*t2-2*D*a**2*ym*yy*t2)/((xm-xx)**2+(ym-yy)**2+1)
          !                                                                                                                f(ind,2) = ((a*mu*xx*sin(2*a*yy))/2-(a*mu*xm*sin(2*a*yy))/2-(a*mu*ym*sin(2*a*xx))/2+(a*mu*yy*sin(2*a*xx))/2-2*a**2*mu*xm*ym+2*a**2*mu*xx*ym+2*a**2*mu*xm*yy-2*a**2*mu*xx*yy+4*a**2*mu*t1+4*a**2*mu*xm*ym*cos(a*xx)**2-4*a**2*mu*xx*ym*cos(a*xx)**2-4*a**2*mu*xm*yy*cos(a*xx)**2+4*a**2*mu*xx*yy*cos(a*xx)**2+4*a**2*mu*xm*ym*cos(a*yy)**2-4*a**2*mu*xx*ym*cos(a*yy)**2-4*a**2*mu*xm*yy*cos(a*yy)**2+4*a**2*mu*xx*yy*cos(a*yy)**2+2*a**2*mu*xm**2*t1+2*a**2*mu*xx**2*t1+2*a**2*mu*ym**2*t1+2*a**2*mu*yy**2*t1-2*a*mu*xm*t4+2*a*mu*xx*t4-2*a*mu*ym*t3+2*a*mu*yy*t3+8*a**2*mu*t1*t2-2*a*xm*cos(a*xx)**2*t4*ta**(0.5)+2*a*xx*cos(a*xx)**2*t4*ta**(0.5)+4*a*ym*t1**2*sin(a*xx)*ta**(0.5)-4*a*xm*cos(a*xx)**2*t5*ta**(0.5)+4*a*xx*cos(a*xx)**2*t5*ta**(0.5)-4*a*yy*t1**2*sin(a*xx)*ta**(0.5)+2*a*ym*t1**2*sin(a*yy)*ta**(0.5)-2*a*yy*t1**2*sin(a*yy)*ta**(0.5)-8*a**2*mu*xm*ym*t8+8*a**2*mu*xx*ym*t8+8*a**2*mu*xm*yy*t8-8*a**2*mu*xx*yy*t8+3*a*xm*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*ta**(0.5)-3*a*xx*cos(a*xx)**2*cos(a*yy)**3*sin(a*xx)*ta**(0.5)-3*a*ym*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*ta**(0.5)+3*a*yy*cos(a*xx)**3*cos(a*yy)**2*sin(a*yy)*ta**(0.5)+a*k*xm*t4*ta**(0.5)-a*k*xx*t4*ta**(0.5)-a*k*ym*t3*ta**(0.5)+a*k*yy*t3*ta**(0.5)-4*a**2*mu*xm*xx*t1-4*a**2*mu*ym*yy*t1-4*a**2*mu*xm*ym*t2+4*a**2*mu*xx*ym*t2+4*a**2*mu*xm*yy*t2-4*a**2*mu*xx*yy*t2+2*a*mu*ym*t1**2*sin(a*xx)+2*a*mu*xm*cos(a*xx)**2*t5-2*a*mu*xx*cos(a*xx)**2*t5-2*a*mu*yy*t1**2*sin(a*xx)+4*a**2*mu*xm**2*t1*t2+4*a**2*mu*xx**2*t1*t2+4*a**2*mu*ym**2*t1*t2+4*a**2*mu*yy**2*t1*t2-8*a**2*mu*xm*xx*t1*t2-8*a**2*mu*ym*yy*t1*t2)/ta
        CASE (2)
          ! Axisimmetric case with div(b)~=0
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
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

          f(ind, 1) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
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

          f(ind, 2) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
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

          f(ind, 3) = f(ind, 1)

        CASE (50:)
          !Do nothing
        CASE DEFAULT
          WRITE (6, *) "Error! Test case not valid"
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
  SUBROUTINE analytical_solution(x, y, u)
    real*8, dimension(:), intent(IN)        :: x, y
    real*8, dimension(:, :), intent(OUT)     :: u
    real*8, dimension(size(u, 1), phys%npv)  :: up
    integer:: i
    real*8 :: a, r(size(x))

    up = 0.
    a = 2*pi
    SELECT CASE (switch%testcase)
    CASE (1)
      IF (switch%axisym) THEN
        WRITE (6, *) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Circular field centered in [xc, yc], n = 2+sin(a*x)*sin(a*y),  u = cos(a*x)*cos(a*y)
      ! Case 9 of the Matlab version: for convergence purpose
      up(:, 1) = 2 + sin(a*x)*sin(a*y)
      up(:, 2) = cos(a*x)*cos(a*y)
      up(:, 3) = up(:, 1)
    CASE (2)
      ! Axisimmetric case with div(b)~=0
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      up(:, 1) = 2 + sin(a*x)*sin(a*y)
      up(:, 2) = cos(a*x)*cos(a*y)
      up(:, 3) = up(:, 1)
    CASE (50:64)
      up(:, 1) = 1.
      up(:, 2) = 0.
      up(:, 3) = up(:, 1)
    CASE (65)
      up(:, 1) = 1.
      up(:, 2) = 0.
      up(:, 3) = up(:, 1)
      r = sqrt((x*phys%lscale - geom%R0)**2 + (y*phys%lscale - 0.75)**2)
      DO i = 1, size(x)
        IF (r(i) .le. 0.05) THEN
          up(i, 2) = 1.
        END IF
      END DO
    CASE DEFAULT
      WRITE (6, *) "Error! Test case not valid"
      STOP
    END SELECT
    ! Convert physical variables to conservative variables
    CALL phys2cons(up, u)
  END SUBROUTINE analytical_solution

  !*****************************************
  ! Analytical gradient
  !****************************************
  SUBROUTINE analytical_gradient(x, y, u, ux, uy)
    real*8, dimension(:), intent(IN)        :: x, y
    real*8, dimension(:, :), intent(IN)      :: u
    real*8, dimension(:, :), intent(OUT)     :: ux, uy
    real*8, dimension(size(u, 1), phys%npv)  :: upx, upy
    real*8, dimension(size(u, 1), phys%npv)  :: up
    real*8 :: a

    upx = 0.
    upy = 0.
    CALL cons2phys(u, up)
    a = 2*pi
    SELECT CASE (switch%testcase)
    CASE (1)
      IF (switch%axisym) THEN
        WRITE (6, *) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      upx(:, 1) = a*cos(a*x)*sin(a*y)
      upx(:, 2) = -a*sin(a*x)*cos(a*y)

      upy(:, 1) = a*sin(a*x)*cos(a*y)
      upy(:, 2) = -a*cos(a*x)*sin(a*y)

      upx(:, 3) = upx(:, 1)
      upy(:, 3) = upx(:, 1)
    CASE (2)
      ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      upx(:, 1) = a*cos(a*x)*sin(a*y)
      upx(:, 2) = -a*sin(a*x)*cos(a*y)

      upy(:, 1) = a*sin(a*x)*cos(a*y)
      upy(:, 2) = -a*cos(a*x)*sin(a*y)

      upx(:, 3) = upx(:, 1)
      upy(:, 3) = upx(:, 1)
    CASE (5)
      ! Do nothing
    CASE (6)
      ! Do nothing
    CASE (50:)
      ! Do nothing
    CASE DEFAULT
      WRITE (6, *) "Error! Test case not valid"
      STOP
    END SELECT
    ! Convert physical variables to conservative variables
    ux(:, 1) = upx(:, 1)
    uy(:, 1) = upy(:, 1)
    ux(:, 2) = (upx(:, 1)*up(:, 2) + up(:, 1)*upx(:, 2))
    uy(:, 2) = (upy(:, 1)*up(:, 2) + up(:, 1)*upy(:, 2))
    ux(:, 3) = upx(:, 3)
    uy(:, 3) = upy(:, 3)
  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, f)
    real*8, dimension(:), intent(IN) :: x, y
    real*8, dimension(:, :), intent(OUT) :: f
    integer                            :: i, n
    real*8  :: a, b, xc, yc, D, mu, k

    n = size(x)

    f = 0.
    SELECT CASE (switch%testcase)
    CASE (1)
      ! Circular field centered in [xc, yc], n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y)
      ! Case 9 of the Matlab version: for convergence purpose
      a = 2*pi
      b = 2*pi
      xc = 0.
      yc = 0.
      D = phys%diff_n
      mu = phys%diff_u
      k = phys%a
      f(:,1) = cos(a * x)**2 * a  * sin(b * y) * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+sin(a * x) * cos(b * y)**2 * b * cos(a * x) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-D * (-sin(a * x) * a**2 * sin(b * y)+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * a**2 * sin(b * y) / 0.10e2-(y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * sin(a * x) * cos(b * y) * b * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * cos(b * y) * b / 0.10e2) / 0.10e2-sin(a * x) * sin(b * y) * b**2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * ((y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * cos(b * y) * b / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * sin(a * x) * cos(b * y) * b * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * sin(b * y) * b**2 / 0.10e2) / 0.10e2)
      f(:,2) = cos(a * x)**3 * a  * sin(b * y) * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * a  / 0.5e1-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y)**2 * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+sin(a * x) * cos(b * y)**3 * b * cos(a * x)**2 * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) / 0.10e2-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y) * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(b * y) * b / 0.5e1-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x)**2 * cos(b * y)**2 * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+k * (y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * cos(a * x) * a  * sin(b * y) / 0.10e2+k * (-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * sin(a * x) * cos(b * y) * b / 0.10e2-mu * (-0.3e1 * cos(a * x) * a**2 * sin(b * y) * cos(b * y) * sin(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * a**2 * cos(b * y)+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2) * (0.2e1 * x-0.2e1 * xc) / 0.20e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-0.3e1 * cos(a * x) * a**2 * sin(b * y) * cos(b * y) * sin(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * a**2 * cos(b * y)) / 0.10e2-(y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) * (0.2e1 * x-0.2e1 * xc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * cos(b * y)**2 * b-sin(a * x)**2 * cos(b * y)**2 * b * a -cos(a * x)**2 * a  * sin(b * y)**2 * b+(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * sin(b * y) * b) / 0.10e2) / 0.10e2-0.3e1 * sin(a * x) * cos(b * y) * b**2 * cos(a * x) * sin(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b**2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * ((y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.10e2) * (0.2e1 * y-0.2e1 * yc) / 0.20e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * ((y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) / 0.10e2-(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (cos(a * x)**2 * a  * sin(b * y) * cos(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * cos(b * y)) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(y-yc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (cos(a * x)**2 * a  * cos(b * y)**2 * b-sin(a * x)**2 * cos(b * y)**2 * b * a -cos(a * x)**2 * a  * sin(b * y)**2 * b+(0.2e1+sin(a * x) * sin(b * y)) * sin(a * x) * a  * sin(b * y) * b) / 0.10e2-(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.3e1 / 0.2e1) * (sin(a * x) * cos(b * y)**2 * b * cos(a * x)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) * (0.2e1 * y-0.2e1 * yc) / 0.20e2+(-x+xc) * (y**2-0.2e1 * y * yc+yc**2+x**2-0.2e1 * x * xc+xc**2)**(-0.1e1 / 0.2e1) * (-0.3e1 * sin(a * x) * cos(b * y) * b**2 * cos(a * x) * sin(b * y)-(0.2e1+sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b**2) / 0.10e2) / 0.10e2)
      f(:, 3) = f(:, 1)
    CASE (2)
      ! Axisimmetric case with div(b)~=0
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      a = 2*pi
      b = 2*pi
      xc = 0.
      yc = 0.
      D = phys%diff_n
      mu = phys%diff_u
      k = phys%a
      f(:, 1) = 0.
      f(:, 2) = 0.
      f(:, 3) = f(:, 1)
    CASE (50:)
      !Do nothing
    CASE DEFAULT
      WRITE (6, *) "Error! Test case not valid"
      STOP

    END SELECT
  END SUBROUTINE body_force
#endif

END MODULE analytical
