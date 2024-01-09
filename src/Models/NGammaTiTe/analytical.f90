!*****************************************
! project: MHDG
! file: analytical.f90
! date: 20/09/2017
!  ******** N-Gamma-Energy system     ****
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
    integer:: i, j, ind, N1D, N2D
    real*8 :: a, b, r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    real*8 :: aux

    u = 0.
    up = 0.
    a = 2*pi
    b = 2*pi
    N2D = size(x, 1)
    N1D = size(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

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
          ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
          up(ind, 1) = 2 + sin(a*xx)*sin(a*yy)
          up(ind, 2) = cos(a*xx)*cos(a*yy)
          up(ind, 3) = 20 + cos(a*xx)*sin(a*yy)
          up(ind, 4) = 10 - sin(a*xx)*cos(a*yy)

        CASE (2)
          ! Axisimmetric case with div(b)~=0
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          up(ind, 1) = 2 + sin(a*xx)*sin(a*yy)
          up(ind, 2) = cos(a*xx)*cos(a*yy)
          up(ind, 3) = 20 + cos(a*xx)*sin(a*yy)
          up(ind, 4) = 10 - sin(a*xx)*cos(a*yy)

        CASE (3)
          ! Axisimmetric case with div(b)~=0
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          up(ind, 1) = 2 + sin(a*tt)
          up(ind, 2) = cos(a*tt)
          up(ind, 3) = 20 + cos(a*tt)
          up(ind, 4) = 10 - sin(a*tt)

        CASE (50:63)
          up(ind, 1) = 1.
          up(ind, 3) = 18.
          up(ind, 4) = 18.
#ifdef NEUTRAL
          up(ind,11) = 0.
#endif
        CASE (64)
          up(ind, 1) = 1.
          up(ind, 3) = 18.
          up(ind, 4) = 18.
          r = sqrt((xx - geom%R0)**2 + yy**2)
          !               if (abs(r-1).lt.1e-6) then
          !               ! outer boundary
          !               endif
        CASE (65)
          up(ind, 1) = 1.
          up(ind, 2) = 0.
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
          upx(ind, 3) = -a*sin(a*xx)*sin(a*yy)
          upx(ind, 4) = -a*cos(a*xx)*cos(a*yy)

          upy(ind, 1) = a*sin(a*xx)*cos(a*yy)
          upy(ind, 2) = -a*cos(a*xx)*sin(a*yy)
          upy(ind, 3) = a*cos(a*xx)*cos(a*yy)
          upy(ind, 4) = a*sin(a*xx)*sin(a*yy)

        CASE (2)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          upx(ind, 1) = a*cos(a*xx)*sin(a*yy)
          upx(ind, 2) = -a*sin(a*xx)*cos(a*yy)
          upx(ind, 3) = -a*sin(a*xx)*sin(a*yy)
          upx(ind, 4) = -a*cos(a*xx)*cos(a*yy)

          upy(ind, 1) = a*sin(a*xx)*cos(a*yy)
          upy(ind, 2) = -a*cos(a*xx)*sin(a*yy)
          upy(ind, 3) = a*cos(a*xx)*cos(a*yy)
          upy(ind, 4) = a*sin(a*xx)*sin(a*yy)
        CASE (3)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          upt(ind, 1) = +a*cos(a*tt)
          upt(ind, 2) = -a*sin(a*tt)
          upt(ind, 3) = -a*sin(a*tt)
          upt(ind, 4) = -a*cos(a*tt)
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
    ux(:, 3) = (upx(:, 1)*up(:, 3) + up(:, 1)*upx(:, 3))
    uy(:, 3) = (upy(:, 1)*up(:, 3) + up(:, 1)*upy(:, 3))
    ut(:, 3) = (upt(:, 1)*up(:, 3) + up(:, 1)*upt(:, 3))
    ux(:, 4) = (upx(:, 1)*up(:, 4) + up(:, 1)*upx(:, 4))
    uy(:, 4) = (upy(:, 1)*up(:, 4) + up(:, 1)*upy(:, 4))
    ut(:, 4) = (upt(:, 1)*up(:, 4) + up(:, 1)*upt(:, 4))

  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, t, f)
    real*8, dimension(:), intent(IN) :: x, y, t
    real*8, dimension(:, :), intent(OUT) :: f
    integer                            ::  n
    real*8  :: a, b, D, mu, csii, csie, kpari, kpare, Mref, epn, tie, pcourr
    integer:: i, j, ind, N1D, N2D
    real*8 :: r, xx, yy, tt, xmax, xmin, ymax, ymin, xm, ym
    real*8 :: aux
    real*8 :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
    real*8 :: t11, t12, t13, t14, t15, t16, t17, t18, t19, t20
    real*8 :: t21, t22, t23, t24, t25, t26, t27, t28, t29, t30
    real*8 :: t31, t32, t33, t34, t35, t36, t37, t38, t39, t40
    real*8 :: t41, t42, t43, t44, t45, t46, t47, t48, t49, t50
    real*8 :: t51, t52, t53, t54, t55, t56, t57, t58, t59, t60
    real*8 :: t61, t62, t63, t64, t65, t66, t67, t68, t69, t70
    real*8 :: cx, cy, sx, sy, cx2, sx2, cy2, sy2, r2, cc, ss, ct, st

    N2D = size(x, 1)
    N1D = size(t, 1)

    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    !  xm = -0.5
    !  ym = -0.5

    a = 2*pi
    f = 0.
    DO i = 1, N2d
      DO j = 1, N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        r = sqrt((xx - xm)**2 + (yy - ym)**2)
        ind = (j - 1)*N2d+i
        SELECT CASE (switch%testcase)

        CASE (1)
          ! Cartesian case, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          a = 2*pi
          b = 2*pi
          D = phys%diff_n
          mu = phys%diff_u
          csii = phys%diff_e
          csie = phys%diff_ee
          kpari = phys%diff_pari
          kpare = phys%diff_pare
          Mref = phys%Mref
          epn = phys%epn
          tie = phys%tie
          pcourr = 1.

          cx = cos(a*xx)
          cy = cos(a*yy)
          sx = sin(a*xx)
          sy = sin(a*yy)
          cx2 = cx**2
          cy2 = cy**2
          sx2 = sx**2
          sy2 = sy**2
          r = ((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)
          r2 = ((xm - xx)**2 + (ym - yy)**2 + 1)**(1.5)

          f(ind, 1) = (D*a**2*sx*sy*2.0D0 + D*a**2*sx*sy*xm**2 + D*a**2*sx*sy*xx**2 + D*a&
            &**2*sx*sy*ym**2 + D*a**2*sx*sy*yy**2 + D*a*cx*sy*xm - D*a*cx*sy*xx + D*a*c&
            &y*sx*ym - D*a*cy*sx*yy - a*cx*r*sy*xm*2.0D0 + a*cx*r*sy*xx*2.0D0 + a*cy*r*&
            &sx*ym*2.0D0 - a*cy*r*sx*yy*2.0D0 - D*a**2*cx*cy*xm*ym*2.0D0 + D*a**2*cx*&
            &cy*xx*ym*2.0D0 + D*a**2*cx*cy*xm*yy*2.0D0 - D*a**2*cx*cy*xx*yy*2.0D0 - D&
            &*a**2*sx*sy*xm*xx*2.0D0 - D*a**2*sx*sy*ym*yy*2.0D0 + a*cx*cy2*r*sx*xm -&
            &a*cx*cy2*r*sx*xx - a*cx2*cy*r*sy*ym + a*cx2*cy*r*sy*yy - a*cx*r*sx*sy2*x&
            &m + a*cx*r*sx*sy2*xx + a*cy*r*sx2*sy*ym - a*cy*r*sx2*sy*yy)/((xm - xx)**2 +&
            &(ym - yy)**2 + 1.0D0)

          f(ind, 2) = (a*(cx*cy*r*xm*4.0D0 - cx*cy*r*xx*4.0D0 + cx*cy*r*ym*4.0D0 - cx*cy*&
            &r*yy*4.0D0 - cy*mu*sx*xm*6.0D0 + cy*mu*sx*xx*6.0D0 - cx*mu*sy*ym*6.0D0 + c&
            &x*mu*sy*yy*6.0D0 + cy*r*sx*xm*6.0D+1 - cy2*r*sx2*xm*2.0D0 - cy*r*sx*xx*6&
            &.0D+1 + cy2*r*sx2*xx*2.0D0 - cx*r*sy*ym*6.0D+1 - cx2*r*sy2*ym*2.0D0 + cx*r&
            &*sy*yy*6.0D+1 + cx2*r*sy2*yy*2.0D0 + r*sx*sy*xm*4.0D0 + r*sx2*sy2*xm*2.0&
            &D0 - r*sx*sy*xx*4.0D0 - r*sx2*sy2*xx*2.0D0 + r*sx*sy*ym*4.0D0 + r*sx2*sy2*&
            &ym*2.0D0 - r*sx*sy*yy*4.0D0 - r*sx2*sy2*yy*2.0D0 + a*cx*cy*mu*1.2D+1 + cx2&
            &*cy*mu*sy*xm*3.0D0 - cx2*cy*mu*sy*xx*3.0D0 + cx*cy2*mu*sx*ym*3.0D0 - cx*&
            &cy2*mu*sx*yy*3.0D0 - cx2*cy*r*sy*xm*8.0D0 + cx2*cy*r*sy*xx*8.0D0 + cx*cy&
            &2*r*sx*ym*8.0D0 - cx*cy2*r*sx*yy*8.0D0 - cy*mu*sx2*sy*xm*3.0D0 + cy*mu*s&
            &x2*sy*xx*3.0D0 - cx*mu*sx*sy2*ym*3.0D0 + cx*mu*sx*sy2*yy*3.0D0 + a*cx*cy&
            &*mu*xm**2*6.0D0 + a*cx*cy*mu*xx**2*6.0D0 + a*cx*cy*mu*ym**2*6.0D0 + a*cx&
            &*cy*mu*yy**2*6.0D0 + cx2*cy**3*r*sx*xm*2.0D0 - cx2*cy**3*r*sx*xx*2.0D0&
            &- cx**3*cy2*r*sy*ym*2.0D0 + cx**3*cy2*r*sy*yy*2.0D0 - a*mu*sx*sy*xm*ym*&
            &1.2D+1 - a*mu*sx2*sy2*xm*ym*6.0D0 + a*mu*sx*sy*xx*ym*1.2D+1 + a*mu*sx2*s&
            &y2*xx*ym*6.0D0 + a*mu*sx*sy*xm*yy*1.2D+1 + a*mu*sx2*sy2*xm*yy*6.0D0 - a*&
            &mu*sx*sy*xx*yy*1.2D+1 - a*mu*sx2*sy2*xx*yy*6.0D0 + a*cx*cy*mu*sx*sy*2.&
            &4D+1 - a*cx*cy*mu*xm*xx*1.2D+1 - a*cx2*cy2*mu*xm*ym*6.0D0 + a*cx2*cy2*mu&
            &*xx*ym*6.0D0 + a*cx2*cy2*mu*xm*yy*6.0D0 - a*cx2*cy2*mu*xx*yy*6.0D0 - a*c&
            &x*cy*mu*ym*yy*1.2D+1 + a*cx2*mu*sy2*xm*ym*6.0D0 + a*cy2*mu*sx2*xm*ym*6&
            &.0D0 - a*cx2*mu*sy2*xx*ym*6.0D0 - a*cy2*mu*sx2*xx*ym*6.0D0 - a*cx2*mu*sy&
            &2*xm*yy*6.0D0 - a*cy2*mu*sx2*xm*yy*6.0D0 + a*cx2*mu*sy2*xx*yy*6.0D0 + a*&
            &cy2*mu*sx2*xx*yy*6.0D0 + cx*cy*r*sx*sy*xm*4.0D0 - cx2*cy*r*sx*sy2*xm*4&
            &.0D0 - cx*cy*r*sx*sy*xx*4.0D0 + cx2*cy*r*sx*sy2*xx*4.0D0 + cx*cy*r*sx*sy&
            &*ym*4.0D0 + cx*cy2*r*sx2*sy*ym*4.0D0 - cx*cy*r*sx*sy*yy*4.0D0 - cx*cy2*r&
            &*sx2*sy*yy*4.0D0 + a*cx*cy*mu*sx*sy*xm**2*1.2D+1 + a*cx*cy*mu*sx*sy*xx&
            &**2*1.2D+1 + a*cx*cy*mu*sx*sy*ym**2*1.2D+1 + a*cx*cy*mu*sx*sy*yy**2*1.&
            &2D+1 - a*cx*cy*mu*sx*sy*xm*xx*2.4D+1 - a*cx*cy*mu*sx*sy*ym*yy*2.4D+1))&
            &/(xm*xx*(-6.0D0) - ym*yy*6.0D0 + xm**2*3.0D0 + xx**2*3.0D0 + ym**2*3.0D0 + y&
            &y**2*3.0D0 + 3.0D0)

          f(ind, 3) = (csii*(-a*xm + a*xx + a**2*sin(a*xx*2.0D0)*3.0D0 + a**2*sx*sy*4.0D+&
            &1 + a*cx2*xm*2.0D0 + a*cy2*xm - a*cx2*xx*2.0D0 - a*cy2*xx + a**2*xm**2*sin(a&
            &*xx*2.0D0)*2.0D0 + a**2*xx**2*sin(a*xx*2.0D0)*2.0D0 + a**2*ym**2*sin(a&
            &*xx*2.0D0) + a**2*yy**2*sin(a*xx*2.0D0) + a**2*cx*sy*4.0D0 + a*cx*cy*ym*&
            &2.0D0 - a*cx*cy*yy*2.0D0 + a*cx*sy*xm*2.0D+1 - a*cx*sy*xx*2.0D+1 + a*cy*sx&
            &*ym*2.0D+1 - a*cy*sx*yy*2.0D+1 - a*sx*sy*xm*2.0D0 + a*sx*sy*xx*2.0D0 - a**&
            &2*cx*cy2*sx*8.0D0 - a**2*xm*xx*sin(a*xx*2.0D0)*4.0D0 + a**2*xm*ym*sin(&
            &a*yy*2.0D0)*2.0D0 - a**2*xx*ym*sin(a*yy*2.0D0)*2.0D0 - a**2*ym*yy*sin(&
            &a*xx*2.0D0)*2.0D0 - a**2*xm*yy*sin(a*yy*2.0D0)*2.0D0 + a**2*xx*yy*sin(&
            &a*yy*2.0D0)*2.0D0 + a**2*cx*sy*xm**2*2.0D0 + a**2*cx*sy*xx**2*2.0D0 + a*&
            &*2*cx*sy*ym**2*2.0D0 + a**2*cx*sy*yy**2*2.0D0 + a**2*sx*sy*xm**2*2.0D+&
            &1 + a**2*sx*sy*xx**2*2.0D+1 + a**2*sx*sy*ym**2*2.0D+1 + a**2*sx*sy*yy**2&
            &*2.0D+1 - a*cx2*cy2*xm*2.0D0 + a*cx2*cy2*xx*2.0D0 - a**2*cx*cy2*sx*xm**2&
            &*4.0D0 - a**2*cx*cy2*sx*xx**2*4.0D0 - a**2*cx*cy2*sx*ym**2*4.0D0 - a**2*&
            &cx*cy2*sx*yy**2*4.0D0 - a**2*cx*cy*xm*ym*4.0D+1 + a**2*cx*cy*xx*ym*4.0&
            &D+1 + a**2*cx*cy*xm*yy*4.0D+1 - a**2*cx*cy*xx*yy*4.0D+1 - a**2*cx*sy*xm*&
            &xx*4.0D0 + a**2*cy*sx*xm*ym*4.0D0 - a**2*cy*sx*xx*ym*4.0D0 - a**2*cy*sx*&
            &xm*yy*4.0D0 + a**2*cy*sx*xx*yy*4.0D0 - a**2*cx*sy*ym*yy*4.0D0 - a**2*sx*&
            &sy*xm*xx*4.0D+1 - a**2*sx*sy*ym*yy*4.0D+1 + a**2*cx*cy2*sx*xm*xx*8.0D0&
            &- a**2*cx2*cy*sy*xm*ym*8.0D0 + a**2*cx2*cy*sy*xx*ym*8.0D0 + a**2*cx2*cy&
            &*sy*xm*yy*8.0D0 - a**2*cx2*cy*sy*xx*yy*8.0D0 + a**2*cx*cy2*sx*ym*yy*8.&
            &0D0 + a*cx*cy*sx*sy*ym*2.0D0 - a*cx*cy*sx*sy*yy*2.0D0))/(xm*xx*(-2.0D0&
            &) - ym*yy*2.0D0 + xm**2 + xx**2 + ym**2 + yy**2 + 1.0D0) - (1.0D0/Mref**2*(3.0D0&
            &**(-epn + 1.0D0)*a**2*cx2*cy2*epn*kpari*xm**2*((-cx2*cy2 + cx*sy*2.0D0&
            &+ 4.0D+1)/Mref)**(epn - 1.0D0)*4.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx2*cy2&
            &*epn*kpari*xx**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)&
            &*4.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*epn*kpari*sx2*sy2*ym**2*((-cx2*cy2&
            &+ cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*4.0D0 + 3.0D0**(-epn + 1.0D0)*&
            &a**2*epn*kpari*sx2*sy2*yy**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)*&
            &*(epn - 1.0D0)*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy2*kpari*xm*&
            &*2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.&
            &0D0)*Mref*a**2*cx2*cy2*kpari*xx**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/&
            &Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy2*kpari*ym**2&
            &*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.0D&
            &0)*Mref*a**2*cx2*cy2*kpari*yy**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mr&
            &ef)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpari*sy*xm**2*((-&
            &cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*M&
            &ref*a**2*cx*kpari*sy*xx**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**e&
            &pn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpari*sy*ym**2*((-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*&
            &*2*cx*kpari*sy*yy**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0&
            &D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cx*cy*kpari*ym*((-cx2*cy2 + cx*sy*2.0D&
            &0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cx*cy*kpari*&
            &yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.&
            &0D0)*Mref*a**2*cx2*kpari*xm**2*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 +&
            &4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*kpari*x&
            &x**2*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3&
            &.0D0**(-epn + 1.0D0)*Mref*a**2*cy2*kpari*ym**2*(cx2 - 1.0D0)*((-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*&
            &*2*cy2*kpari*yy**2*(cx2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref&
            &)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*kpari*sx*sy*xm*((-cx2*cy2 +&
            &cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*kp&
            &ari*sx*sy*xx*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0&
            &**(-epn + 1.0D0)*Mref*a**2*cx2*kpari*xm*xx*(cy2 - 1.0D0)*((-cx2*cy2 + cx&
            &*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*c&
            &y2*kpari*ym*yy*(cx2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**e&
            &pn*4.0D0 - 3.0D0**(-epn + 1.0D0)*a**2*cx2*cy2*epn*kpari*xm*xx*((-cx2*c&
            &y2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 - 3.0D0**(-epn + 1.0D0&
            &)*a**2*epn*kpari*sx2*sy2*ym*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref&
            &)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cx*cy2*kpari*sx*xm&
            &*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0**(-epn + 1.0D&
            &0)*Mref*a*cx*cy2*kpari*sx*xx*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)*&
            &*epn*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cx2*cy*kpari*sy*ym*((-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*&
            &cx2*cy*kpari*sy*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*2.0D0&
            &- 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy2*kpari*xm*xx*((-cx2*cy2 + cx*s&
            &y*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2&
            &*cy2*kpari*ym*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 + 3&
            &.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpari*sy*xm*xx*((-cx2*cy2 + cx*sy*2.&
            &0D0 + 4.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpar&
            &i*sx*xm*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 - 3.0D0**&
            &(-epn + 1.0D0)*Mref*a**2*cy*kpari*sx*xx*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.&
            &0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpari*sx*x&
            &m*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn +&
            &1.0D0)*Mref*a**2*cy*kpari*sx*xx*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/&
            &Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpari*sy*ym*yy*(&
            &(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)&
            &*a**2*cx**3*cy2*epn*kpari*sy*xm**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/&
            &Mref)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx**3*cy2*epn*kp&
            &ari*sy*xx**2*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0&
            &D0 - 3.0D0**(-epn + 1.0D0)*a**2*cx**3*cy2*epn*kpari*sy*xm*xx*((-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*1.6D+1 - 3.0D0**(-epn + 1.0D0&
            &)*a**2*cx2*cy**3*epn*kpari*sx*xm*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)&
            &/Mref)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx2*cy**3*epn*k&
            &pari*sx*xx*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.&
            &0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx2*cy**3*epn*kpari*sx*xm*yy*((-cx2*c&
            &y2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 - 3.0D0**(-epn + 1.0D0&
            &)*a**2*cx2*cy**3*epn*kpari*sx*xx*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)&
            &/Mref)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx*cy*epn*kpari&
            &*sx*sy*xm*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0&
            &D0 - 3.0D0**(-epn + 1.0D0)*a**2*cx*cy*epn*kpari*sx*sy*xx*ym*((-cx2*cy2&
            &+ cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 - 3.0D0**(-epn + 1.0D0)*&
            &a**2*cx*cy*epn*kpari*sx*sy*xm*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mr&
            &ef)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*a**2*cx*cy*epn*kpari*sx&
            &*sy*xx*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 -&
            &(3.0D0**(-epn + 1.0D0)*Mref*a**2*cx**4*cy2*epn*kpari*xm**2*(cy2 - 1.0D&
            &0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0)/(-cx2*cy2 + cx*s&
            &y*2.0D0 + 4.0D+1) - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx**4*cy2*epn*kpari&
            &*xx**2*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0&
            &)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2&
            &*cy**4*epn*kpari*ym**2*(cx2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/&
            &Mref)**epn*4.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) - (3.0D0**(-epn + 1.0D&
            &0)*Mref*a**2*cx2*cy**4*epn*kpari*yy**2*(cx2 - 1.0D0)*((-cx2*cy2 + cx*s&
            &y*2.0D0 + 4.0D+1)/Mref)**epn*4.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) + 3.&
            &0D0**(-epn + 1.0D0)*Mref*a**2*cx*cy*kpari*sx*sy*xm*ym*((-cx2*cy2 + cx*&
            &sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx&
            &*cy*kpari*sx*sy*xx*ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.&
            &0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*cy*kpari*sx*sy*xm*yy*((-cx2*c&
            &y2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a&
            &**2*cx*cy*kpari*sx*sy*xx*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**&
            &epn*8.0D0 - 3.0D0**(-epn + 1.0D0)*a**2*cx**3*cy**3*epn*kpari*sx*sy*xm*&
            &ym*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 + 3.0D0**&
            &(-epn + 1.0D0)*a**2*cx**3*cy**3*epn*kpari*sx*sy*xx*ym*((-cx2*cy2 + cx*&
            &sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*8.0D0 + 3.0D0**(-epn + 1.0D0)*a**2&
            &*cx**3*cy**3*epn*kpari*sx*sy*xm*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/&
            &Mref)**(epn - 1.0D0)*8.0D0 - 3.0D0**(-epn + 1.0D0)*a**2*cx**3*cy**3*epn*&
            &kpari*sx*sy*xx*yy*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0&
            &)*8.0D0 + (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*cy2*epn*kpari*sy*ym**2*(&
            &cx2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0)/(-cx2*&
            &cy2 + cx*sy*2.0D0 + 4.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*cy2*epn*&
            &kpari*sy*yy**2*(cx2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**e&
            &pn*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*&
            &a**2*cx**4*cy2*epn*kpari*xm*xx*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 +&
            &4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) + (3.0D0**(-&
            &epn + 1.0D0)*Mref*a**2*cx2*cy**4*epn*kpari*ym*yy*(cx2 - 1.0D0)*((-cx2*&
            &cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.&
            &0D+1) - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy*epn*kpari*sx*xm*ym*(cy&
            &2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy*epn*kp&
            &ari*sx*xx*ym*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn&
            &*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*a*&
            &*2*cx2*cy*epn*kpari*sx*xm*yy*(cy2 - 1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.&
            &0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1) - (3.0D0**(-ep&
            &n + 1.0D0)*Mref*a**2*cx2*cy*epn*kpari*sx*xx*yy*(cy2 - 1.0D0)*((-cx2*cy&
            &2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2 + cx*sy*2.0D0 + 4.0D&
            &+1) - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*cy2*epn*kpari*sy*ym*yy*(cx2 -&
            &1.0D0)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*1.6D+1)/(-cx2*cy2&
            &+ cx*sy*2.0D0 + 4.0D+1)))/((xm - xx)**2*9.0D0 + (ym - yy)**2*9.0D0 + 9.0D0) + M&
            &ref*cx*cy*pcourr*((a*sx*(xm - xx)*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*sx*2.0D&
            &0)*(2.0D0/3.0D0))/(Mref*r) + (a*cx*(ym - yy)*(cy - sy*5.0D0 + cy*sx*sy)*(4&
            &.0D0/3.0D0))/(Mref*r)) - (a*cx*cy*(ym - yy)*(cx*sy*1.0D+2 + cx2*sy2*5.0D&
            &0 - sx*sy*1.0D+1 - sx2*sy2*5.0D0 + cx*cy2*sx*4.0D0 - cx**3*cy2*sy + cx*cy2*s&
            &x2*sy*2.0D0))/(r*3.0D0) + (sqrt(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0 - 2.0D&
            &+1)/Mref)*(sx*sy + 2.0D0)**2*(-cx2*cy2 + cx*sy*2.0D0 + cy*sx*2.0D0 + 2.0D+&
            &1))/(tie*(cy*sx - 1.0D+1)*2.0D0) + (a*cx*cy2*(xm - xx)*(cx*1.0D+1 + sx*1.0&
            &D+2 + cx2*sx*2.0D0 + cx2*sy*4.0D0 - cx2*cy2*sx*3.0D0 + cx*sx*sy*1.0D+1))/(&
            &r*3.0D0) - (cx*cy*(xm*2.0D0 - xx*2.0D0)*(sx*sy + 2.0D0)*(ym - yy)*(-cx2*cy&
            &2 + cx*sy*5.0D0 + 1.0D+2))/(r2*6.0D0) + (cx*cy*(ym*2.0D0 - yy*2.0D0)*(sx*s&
            &y + 2.0D0)*(xm - xx)*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/(r2*6.0D0) - (a*cx*s&
            &y*(sx*sy + 2.0D0)*(xm - xx)*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/(r*3.0D0) + (&
            &a*cy*sx*(sx*sy + 2.0D0)*(ym - yy)*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/(r*3.&
            &0D0)

          f(ind, 4) = (csie*(a*ym - a*yy - a**2*sin(a*yy*2.0D0)*3.0D0 + a**2*sx*sy*2.0D+1&
            &- a*cx2*ym - a*cy2*ym*2.0D0 + a*cx2*yy + a*cy2*yy*2.0D0 - a**2*xm**2*sin(a*&
            &yy*2.0D0) - a**2*xx**2*sin(a*yy*2.0D0) - a**2*ym**2*sin(a*yy*2.0D0)*2.&
            &0D0 - a**2*yy**2*sin(a*yy*2.0D0)*2.0D0 - a**2*cy*sx*4.0D0 + a*cx2*cy2*ym&
            &*2.0D0 - a*cx2*cy2*yy*2.0D0 + a*cx*sy*xm*1.0D+1 - a*cx*sy*xx*1.0D+1 + a*cy&
            &*sx*ym*1.0D+1 - a*cy*sx*yy*1.0D+1 + a*sx*sy*ym*2.0D0 - a*sx*sy*yy*2.0D0 +&
            &a**2*cx2*cy*sy*8.0D0 - a**2*xm*ym*sin(a*xx*2.0D0)*2.0D0 + a**2*xx*ym*s&
            &in(a*xx*2.0D0)*2.0D0 + a**2*xm*xx*sin(a*yy*2.0D0)*2.0D0 + a**2*xm*yy*s&
            &in(a*xx*2.0D0)*2.0D0 - a**2*xx*yy*sin(a*xx*2.0D0)*2.0D0 + a**2*ym*yy*s&
            &in(a*yy*2.0D0)*4.0D0 - a**2*cy*sx*xm**2*2.0D0 - a**2*cy*sx*xx**2*2.0D0&
            &- a**2*cy*sx*ym**2*2.0D0 - a**2*cy*sx*yy**2*2.0D0 + a**2*sx*sy*xm**2*1.&
            &0D+1 + a**2*sx*sy*xx**2*1.0D+1 + a**2*sx*sy*ym**2*1.0D+1 + a**2*sx*sy*yy&
            &**2*1.0D+1 - a*cx*cy*xm*2.0D0 + a*cx*cy*xx*2.0D0 + a**2*cx2*cy*sy*xm**2*&
            &4.0D0 + a**2*cx2*cy*sy*xx**2*4.0D0 + a**2*cx2*cy*sy*ym**2*4.0D0 + a**2*c&
            &x2*cy*sy*yy**2*4.0D0 - a**2*cx*cy*xm*ym*2.0D+1 + a**2*cx*cy*xx*ym*2.0D&
            &+1 + a**2*cx*cy*xm*yy*2.0D+1 - a**2*cx*cy*xx*yy*2.0D+1 + a**2*cy*sx*xm*x&
            &x*4.0D0 - a**2*cx*sy*xm*ym*4.0D0 + a**2*cx*sy*xx*ym*4.0D0 + a**2*cx*sy*x&
            &m*yy*4.0D0 - a**2*cx*sy*xx*yy*4.0D0 + a**2*cy*sx*ym*yy*4.0D0 - a**2*sx*s&
            &y*xm*xx*2.0D+1 - a**2*sx*sy*ym*yy*2.0D+1 - a**2*cx2*cy*sy*xm*xx*8.0D0 +&
            &a**2*cx*cy2*sx*xm*ym*8.0D0 - a**2*cx*cy2*sx*xx*ym*8.0D0 - a**2*cx*cy2*&
            &sx*xm*yy*8.0D0 + a**2*cx*cy2*sx*xx*yy*8.0D0 - a**2*cx2*cy*sy*ym*yy*8.0&
            &D0 - a*cx*cy*sx*sy*xm*2.0D0 + a*cx*cy*sx*sy*xx*2.0D0))/(xm*xx*(-2.0D0)&
            &- ym*yy*2.0D0 + xm**2 + xx**2 + ym**2 + yy**2 + 1.0D0) - Mref*cx*cy*pcourr*((a*&
            &sx*(xm - xx)*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*sx*2.0D0)*(2.0D0/3.0D0))/(Mr&
            &ef*r) + (a*cx*(ym - yy)*(cy - sy*5.0D0 + cy*sx*sy)*(4.0D0/3.0D0))/(Mref*r)&
            &) + (a*cx2*cy*(ym - yy)*(cy - sy*5.0D0 + cy*sx*sy)*(1.0D+1/3.0D0))/r - (sqrt&
            &(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0 - 2.0D+1)/Mref)*(sx*sy + 2.0D0)**2*(-&
            &cx2*cy2 + cx*sy*2.0D0 + cy*sx*2.0D0 + 2.0D+1))/(tie*(cy*sx - 1.0D+1)*2.0D0&
            &) + (a*cx*cy*sx*(xm - xx)*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*sx*2.0D0)*(5.0D0/&
            &3.0D0))/r + (a*cx*sy*(cy*sx - 1.0D+1)*(sx*sy + 2.0D0)*(xm - xx)*(5.0D0/3.0&
            &D0))/r - (a*cy*sx*(cy*sx - 1.0D+1)*(sx*sy + 2.0D0)*(ym - yy)*(5.0D0/3.0D0)&
            &)/r + (cx*cy*(xm*2.0D0 - xx*2.0D0)*(cy*sx - 1.0D+1)*(sx*sy + 2.0D0)*(ym - yy&
            &)*(5.0D0/6.0D0))/r2 - (cx*cy*(ym*2.0D0 - yy*2.0D0)*(cy*sx - 1.0D+1)*(sx*&
            &sy + 2.0D0)*(xm - xx)*(5.0D0/6.0D0))/r2 + (3.0D0**(-epn - 1.0D0)*a*kpare*(&
            &-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(cx*cy*xm*(-1.0D+1) + cx*cy*xx*1.0D&
            &+1 + sx*sy*ym*1.0D+1 - sx*sy*yy*1.0D+1 + cx*cy2*sx*xm - cx*cy2*sx*xx - cy*sx&
            &2*sy*ym + cy*sx2*sy*yy + a*cy*sx*xm**2*1.0D+1 - a*cy2*sx2*xm**2 + a*cy*sx*&
            &xx**2*1.0D+1 - a*cy2*sx2*xx**2 + a*cy*sx*ym**2*1.0D+1 - a*cy2*sx2*ym**2 +&
            &a*cy*sx*yy**2*1.0D+1 - a*cy2*sx2*yy**2 - a*cy*sx*xm*xx*2.0D+1 + a*cy2*sx&
            &2*xm*xx*2.0D0 - a*cx*sy*xm*ym*2.0D+1 + a*cx*sy*xx*ym*2.0D+1 + a*cx*sy*xm&
            &*yy*2.0D+1 - a*cx*sy*xx*yy*2.0D+1 - a*cy*sx*ym*yy*2.0D+1 + a*cy2*sx2*ym*&
            &yy*2.0D0 + a*cx2*cy2*epn*ym**2 + a*cx2*cy2*epn*yy**2 + a*epn*sx2*sy2*xm*&
            &*2 + a*epn*sx2*sy2*xx**2 - a*cx2*cy2*epn*ym*yy*2.0D0 - a*epn*sx2*sy2*xm*&
            &xx*2.0D0 + a*cx*cy*sx*sy*xm*ym*2.0D0 - a*cx*cy*sx*sy*xx*ym*2.0D0 - a*cx*&
            &cy*sx*sy*xm*yy*2.0D0 + a*cx*cy*sx*sy*xx*yy*2.0D0 + a*cx*cy*epn*sx*sy*x&
            &m*ym*2.0D0 - a*cx*cy*epn*sx*sy*xx*ym*2.0D0 - a*cx*cy*epn*sx*sy*xm*yy*2&
            &.0D0 + a*cx*cy*epn*sx*sy*xx*yy*2.0D0)*2.0D0)/(Mref*(cy*sx - 1.0D+1)*(x&
            &m*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2 + ym**2 + yy**2 + 1.0D0))

        CASE (2)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          a = 2*pi
          b = 2*pi

          D = phys%diff_n
          mu = phys%diff_u
          csii = phys%diff_e
          csie = phys%diff_ee
          kpari = phys%diff_pari
          kpare = phys%diff_pare
          Mref = phys%Mref
          epn = phys%epn
          tie = phys%tie
          pcourr = 1.

          cx = cos(a*xx)
          cy = cos(a*yy)
          sx = sin(a*xx)
          sy = sin(a*yy)
          cx2 = cx**2
          cy2 = cy**2
          sx2 = sx**2
          sy2 = sy**2
          r = ((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)
          r2 = ((xm - xx)**2 + (ym - yy)**2 + 1)**(1.5)
          cc = cx*cy
          ss = sx*sy

          f(ind, 1) = (1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy&
            &**2)**2*(cc*r*xx*ym**3*(-2.0D0) - cc*r*xx**3*ym*2.0D0 + cc*r*xx*yy**3*&
            &2.0D0 + D*a**2*ss*xx**5*6.0D0 + D*a**2*cc*xx**4*ym*4.0D0 - D*a**2*cc*xx*&
            &*4*yy*4.0D0 - D*a**2*ss*xm*xx**4*1.0D+1 + D*a**2*ss*xm**4*xx + D*a**2*ss&
            &*xx*ym**4 + D*a**2*ss*xx*yy**4 + D*a**2*cc*xx**2*ym**3*2.0D0 - D*a**2*cc&
            &*xx**2*yy**3*2.0D0 + D*a**2*ss*xm**2*xx**3*9.0D0 - D*a**2*ss*xm**3*xx*&
            &*2*4.0D0 + D*a**2*ss*xx**3*ym**2*5.0D0 + D*a**2*ss*xx**3*yy**2*5.0D0 - D&
            &*a*cx*sy*xm**4 - D*a*cx*sy*xx**4*6.0D0 + a*cx*r2*sy*xx**5*2.0D0 + cx*cy*&
            &r*xx**3*yy*2.0D0 - cc*r*ss*xx*ym**3 - cc*r*ss*xx**3*ym + cc*r*ss*xx*yy**&
            &3 + cc*r*ss*xx**3*yy + cc*r*xm*xx**2*ym*4.0D0 - cc*r*xm**2*xx*ym*2.0D0 - c&
            &c*r*xm*xx**2*yy*4.0D0 + cc*r*xm**2*xx*yy*2.0D0 - cc*r*xx*ym*yy**2*6.0D&
            &0 + cc*r*xx*ym**2*yy*6.0D0 - D*a**2*ss*xm*xx**2*ym**2*4.0D0 + D*a**2*ss*&
            &xm**2*xx*ym**2*2.0D0 - D*a**2*ss*xm*xx**2*yy**2*4.0D0 + D*a**2*ss*xm**&
            &2*xx*yy**2*2.0D0 + D*a**2*ss*xx*ym**2*yy**2*6.0D0 + D*a*cx*sy*xm*xx**3&
            &*1.2D+1 + D*a*cx*sy*xm**3*xx*5.0D0 - D*a*cy*sx*xm*ym**3 - D*a*cy*sx*xm**&
            &3*ym + D*a*cy*sx*xx*ym**3*2.0D0 + D*a*cy*sx*xx**3*ym*2.0D0 + D*a*cy*sx*x&
            &m*yy**3 + D*a*cy*sx*xm**3*yy - D*a*cy*sx*xx*yy**3*2.0D0 - D*a*cy*sx*xx**&
            &3*yy*2.0D0 - a*cx*cy2*r2*sx*xx**5 + a*cx*r2*sx*sy2*xx**5 - a*cx*r2*sy*xm&
            &*xx**4*2.0D0 + a*cy*r2*sx*xx**4*ym*2.0D0 - a*cy*r2*sx*xx**4*yy*2.0D0 + c&
            &c*r*ss*xm*xx**2*ym*2.0D0 - cc*r*ss*xm**2*xx*ym - cc*r*ss*xm*xx**2*yy*2&
            &.0D0 + cc*r*ss*xm**2*xx*yy - cc*r*ss*xx*ym*yy**2*3.0D0 + cc*r*ss*xx*ym**&
            &2*yy*3.0D0 - D*a*cx*sy*xm**2*xx**2*1.1D+1 - D*a*cx*sy*xm**2*ym**2 - D*a*&
            &cx*sy*xx**2*ym**2*5.0D0 - D*a*cx*sy*xm**2*yy**2 - D*a*cx*sy*xx**2*yy**&
            &2*5.0D0 - D*a**2*cc*xm*xx*ym**3*2.0D0 - D*a**2*cc*xm*xx**3*ym*8.0D0 - D*&
            &a**2*cc*xm**3*xx*ym*2.0D0 + D*a**2*cc*xm*xx*yy**3*2.0D0 + D*a**2*cc*xm&
            &*xx**3*yy*8.0D0 - D*a**2*ss*xx*ym*yy**3*4.0D0 - D*a**2*ss*xx*ym**3*yy*&
            &4.0D0 - D*a**2*ss*xx**3*ym*yy*1.0D+1 + D*a**2*cc*xm*xx*ym**2*yy*6.0D0 +&
            &D*a**2*ss*xm*xx**2*ym*yy*8.0D0 - D*a**2*ss*xm**2*xx*ym*yy*4.0D0 + D*a*&
            &*2*cx*cy*xm**2*xx**2*ym*6.0D0 - D*a**2*cx*cy*xm**2*xx**2*yy*6.0D0 + D*&
            &a**2*cx*cy*xx**2*ym*yy**2*6.0D0 - D*a**2*cx*cy*xx**2*ym**2*yy*6.0D0 +&
            &D*a*cx*sy*xm*xx*ym**2*3.0D0 - D*a*cy*sx*xm*xx**2*ym*4.0D0 + D*a*cy*sx*&
            &xm**2*xx*ym*4.0D0 + D*a*cx*sy*xm*xx*yy**2*3.0D0 + D*a*cy*sx*xm*xx**2*y&
            &y*4.0D0 - D*a*cy*sx*xm**2*xx*yy*4.0D0 + D*a*cx*sy*xm**2*ym*yy*2.0D0 - D*&
            &a*cy*sx*xm*ym*yy**2*3.0D0 + D*a*cy*sx*xm*ym**2*yy*3.0D0 + D*a*cx*sy*xx&
            &**2*ym*yy*1.0D+1 + D*a*cy*sx*xx*ym*yy**2*6.0D0 - D*a*cy*sx*xx*ym**2*yy&
            &*6.0D0 + a*cx*cy2*r2*sx*xm*xx**4 - a*cx2*cy*r2*sy*xx**4*ym + a*cx2*cy*r2&
            &*sy*xx**4*yy - a*cx*r2*sx*sy2*xm*xx**4 + a*cy*r2*sx2*sy*xx**4*ym - a*cy*&
            &r2*sx2*sy*xx**4*yy + D*a**2*cx*cy*xm**3*xx*yy*2.0D0 - D*a**2*cx*cy*xm*&
            &xx*ym*yy**2*6.0D0 - D*a*cx*sy*xm*xx*ym*yy*6.0D0))/xx

          f(ind, 2) = mu*((xx*(a**2*cc*(ss + 2.0D0) + a**2*cc*ss*3.0D0 - ((ym - yy)*((a*cx*&
            &sy*(ss + 2.0D0) - a*cx*cy2*sx)/(r*xx) + ((xm - xx)*(a**2*ss*(ss + 2.0D0) + a**&
            &2*cx2*cy2 - a**2*cx2*sy2 - a**2*cy2*sx2))/(r*xx) + ((a**2*cc*(ss + 2.0D0) +&
            &a**2*cc*ss*3.0D0)*(ym - yy))/(r*xx) + (1.0D0/xx**2*(xm - xx)*(a*cx*sy*(s&
            &s + 2.0D0) - a*cx*cy2*sx))/r - (1.0D0/xx**2*(ym - yy)*(a*cy*sx*(ss + 2.0D0) -&
            &a*cx2*cy*sy))/r - ((xm - xx)*(a*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx)*(1.0D0/x&
            &x**2*(xm*2.0D0 - xx*2.0D0) + 1.0D0/xx**3*(xm - xx)**2*2.0D0 + 1.0D0/xx**3*&
            &(ym - yy)**2*2.0D0))/(r2*xx*2.0D0) + ((ym - yy)*(a*cy*sx*(ss + 2.0D0) - a*cx&
            &2*cy*sy)*(1.0D0/xx**2*(xm*2.0D0 - xx*2.0D0) + 1.0D0/xx**3*(xm - xx)**2*2&
            &.0D0 + 1.0D0/xx**3*(ym - yy)**2*2.0D0))/(r2*xx*2.0D0)))/(r*xx) - (1.0D0/&
            &xx**2*(((xm - xx)*(a*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx))/(r*xx) - ((ym - yy)*&
            &(a*cy*sx*(ss + 2.0D0) - a*cx2*cy*sy))/(r*xx))*(ym - yy))/r + ((((xm - xx)*(a&
            &*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx))/(r*xx) - ((ym - yy)*(a*cy*sx*(ss + 2.0D0&
            &) - a*cx2*cy*sy))/(r*xx))*(ym - yy)*(1.0D0/xx**2*(xm*2.0D0 - xx*2.0D0) + 1&
            &.0D0/xx**3*(xm - xx)**2*2.0D0 + 1.0D0/xx**3*(ym - yy)**2*2.0D0))/(r2*xx*&
            &2.0D0)) + ((((xm - xx)*(a*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx))/(r*xx) - ((ym - y&
            &y)*(a*cy*sx*(ss + 2.0D0) - a*cx2*cy*sy))/(r*xx))*(ym - yy))/(r*xx) + a*cy*&
            &sx*(ss + 2.0D0) - a*cx2*cy*sy)/xx + a**2*cc*(ss + 2.0D0) + a**2*cc*ss*3.0D0 -&
            &((xm - xx)*((a*cy*sx*(ss + 2.0D0) - a*cx2*cy*sy)/(r*xx) + ((ym - yy)*(a**2*s&
            &s*(ss + 2.0D0) + a**2*cx2*cy2 - a**2*cx2*sy2 - a**2*cy2*sx2))/(r*xx) + ((a**&
            &2*cc*(ss + 2.0D0) + a**2*cc*ss*3.0D0)*(xm - xx))/(r*xx) + (1.0D0/xx**3*(ym&
            &*2.0D0 - yy*2.0D0)*(xm - xx)*(a*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx))/(r2*2.0&
            &D0) - (1.0D0/xx**3*(ym*2.0D0 - yy*2.0D0)*(ym - yy)*(a*cy*sx*(ss + 2.0D0) - a&
            &*cx2*cy*sy))/(r2*2.0D0)))/(r*xx) - (1.0D0/xx**3*(ym*2.0D0 - yy*2.0D0)*&
            &(((xm - xx)*(a*cx*sy*(ss + 2.0D0) - a*cx*cy2*sx))/(r*xx) - ((ym - yy)*(a*cy*&
            &sx*(ss + 2.0D0) - a*cx2*cy*sy))/(r*xx))*(xm - xx))/(r2*2.0D0)) - ((cx2*cy2&
            &*(ym - yy)*(ss + 2.0D0)*(1.0D0/xx**2*(xm*2.0D0 - xx*2.0D0) + 1.0D0/xx**3*(&
            &xm - xx)**2*2.0D0 + 1.0D0/xx**3*(ym - yy)**2*2.0D0))/(r2*2.0D0) + (a*cx**3&
            &*cy2*sy*(ym - yy))/r - (a*cx*cy2*sx*(ym - yy)*(ss + 2.0D0)*2.0D0)/r)/xx + Mr&
            &ef*(((xm - xx)*(((ss*2.0D0 + 4.0D0)*(a*cx*cy + a*cx2*cy*sy))/(Mref*3.0D0&
            &) + (a*ss*(ss*2.0D0 + 4.0D0))/(Mref*3.0D0) + (a*cy*sx*(cx2*cy2*(-1.0D0/2&
            &.0D0) + cx*sy + 2.0D+1)*(2.0D0/3.0D0))/Mref - (a*cy*sx*(cy*sx - 1.0D+1)*(2&
            &.0D0/3.0D0))/Mref))/(r*xx) + ((ym - yy)*(((ss*2.0D0 + 4.0D0)*(a*ss - a*cx*&
            &cy2*sx))/(Mref*3.0D0) + (a*cc*(ss*2.0D0 + 4.0D0))/(Mref*3.0D0) - (a*cx*s&
            &y*(cx2*cy2*(-1.0D0/2.0D0) + cx*sy + 2.0D+1)*(2.0D0/3.0D0))/Mref + (a*cx*&
            &sy*(cy*sx - 1.0D+1)*(2.0D0/3.0D0))/Mref))/(r*xx)) + (cx2*cy2*1.0D0/xx*&
            &*3*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*(ss + 2.0D0))/(r2*2.0D0) + (a*cx2*cy**3&
            &*sx*(xm - xx))/(r*xx) - (a*cx2*cy*sy*(xm - xx)*(ss + 2.0D0)*2.0D0)/(r*xx)

          f(ind, 3) = -kpari*(((3.0D0**(-epn - 1.0D0)*xx*(ym - yy)*((-cx2*cy2 + cx*sy*2.0&
            &D0 + 4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**&
            &2*2.0D0 + ym**2 + yy**2)**2*(a**2*cy2*ym**3*(-2.0D0) + a**2*cy2*yy**3*2.&
            &0D0 - a*ss*xm*ym*2.0D0 + a*ss*xx*ym*4.0D0 + a*ss*xm*yy*2.0D0 - a*ss*xx*yy*&
            &4.0D0 + a*cx*cy*ym**2*2.0D0 + a*cx*cy*yy**2*2.0D0 + a**2*cx2*cy2*ym**3*4&
            &.0D0 - a**2*cx2*cy2*yy**3*4.0D0 + a**2*cc*ss*xm**3*4.0D0 - a**2*cc*ss*xx&
            &**3*8.0D0 + a**2*cy*sx*xm**3*2.0D0 - a**2*cy*sx*xx**3*4.0D0 - a**2*cx*sy&
            &*ym**3*2.0D0 + a**2*cx*sy*yy**3*2.0D0 - a**2*cy2*xm**2*ym*2.0D0 - a**2*c&
            &y2*xx**2*ym*4.0D0 + a**2*cy2*xm**2*yy*2.0D0 + a**2*cy2*xx**2*yy*4.0D0 -&
            &a**2*cy2*ym*yy**2*6.0D0 + a**2*cy2*ym**2*yy*6.0D0 + a**2*cx2*cy2*xm**2&
            &*ym*4.0D0 + a**2*cx2*cy2*xx**2*ym*8.0D0 - a**2*cx2*cy2*xm**2*yy*4.0D0 -&
            &a**2*cx2*cy2*xx**2*yy*8.0D0 + a**2*cx2*cy2*ym*yy**2*1.2D+1 - a**2*cx2*&
            &cy2*ym**2*yy*1.2D+1 + a**2*cc*ss*xm*xx**2*1.6D+1 - a**2*cc*ss*xm**2*xx&
            &*1.2D+1 + a**2*cy*sx*xm*xx**2*8.0D0 - a**2*cy*sx*xm**2*xx*6.0D0 + a**2*c&
            &c*ss*xm*ym**2*4.0D0 - a**2*cc*ss*xx*ym**2*4.0D0 + a**2*cc*ss*xm*yy**2*&
            &4.0D0 - a**2*cc*ss*xx*yy**2*4.0D0 - a**2*cx*sy*xm**2*ym*2.0D0 + a**2*cy*&
            &sx*xm*ym**2*2.0D0 - a**2*cx*sy*xx**2*ym*4.0D0 - a**2*cy*sx*xx*ym**2*2.&
            &0D0 + a**2*cx*sy*xm**2*yy*2.0D0 + a**2*cy*sx*xm*yy**2*2.0D0 + a**2*cx*sy&
            &*xx**2*yy*4.0D0 - a**2*cy*sx*xx*yy**2*2.0D0 - a**2*cx*sy*ym*yy**2*6.0D&
            &0 + a**2*cx*sy*ym**2*yy*6.0D0 + a*cx*cy*xm*xx*2.0D0 - a*cx*cy*ym*yy*4.0D&
            &0 + a*cx2*cy*sy*ym**2*2.0D0 + a*cx2*cy*sy*yy**2*2.0D0 + a**2*cy2*xm*xx*y&
            &m*4.0D0 - a**2*cy2*xm*xx*yy*4.0D0 - a**2*cx2*cy2*xm*xx*ym*8.0D0 + a**2*c&
            &x2*cy2*xm*xx*yy*8.0D0 + a**2*cx*sy*xm*xx*ym*4.0D0 - a**2*cx*sy*xm*xx*y&
            &y*4.0D0 - a**2*cc*ss*xm*ym*yy*8.0D0 + a**2*cc*ss*xx*ym*yy*8.0D0 - a**2*c&
            &y*sx*xm*ym*yy*4.0D0 + a**2*cy*sx*xx*ym*yy*4.0D0 + a*cx2*cy*sy*xm*xx*2.&
            &0D0 + a*cx*cy2*sx*xm*ym*2.0D0 - a*cx*cy2*sx*xx*ym*4.0D0 - a*cx*cy2*sx*xm&
            &*yy*2.0D0 + a*cx*cy2*sx*xx*yy*4.0D0 - a*cx2*cy*sy*ym*yy*4.0D0))/Mref - (&
            &3.0D0**(-epn - 1.0D0)*a*(ym - yy)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)&
            &**epn*(-xm*xx - ym*yy*2.0D0 + xm**2 + ym**2 + yy**2)*1.0D0/(xm*xx*(-2.0D0)&
            &- ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*(ss*ym - ss*yy + cx*cy*&
            &xm - cx*cy*xx + cx2*cy*sy*xm - cx2*cy*sy*xx - cx*cy2*sx*ym + cx*cy2*sx*yy)*2&
            &.0D0)/Mref + (3.0D0**(-epn - 1.0D0)*1.0D0/Mref**2*a**2*epn*sx*xx*(ym - y&
            &y)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*(sy - cx*cy2)*(&
            &ss*ym - ss*yy + cx*cy*xm - cx*cy*xx + cx2*cy*sy*xm - cx2*cy*sy*xx - cx*cy2*sx*&
            &ym + cx*cy2*sx*yy)*4.0D0)/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.&
            &0D0 + ym**2 + yy**2))/xx - (3.0D0**(-epn - 1.0D0)*(xm - xx)*((-cx2*cy2 + cx*sy&
            &*2.0D0 + 4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 +&
            &xx**2*2.0D0 + ym**2 + yy**2)**2*(a*ss*xm**2*2.0D0 + a*ss*xx**2*4.0D0 + a**&
            &2*cx2*xm**3*2.0D0 - a**2*cx2*xx**3*4.0D0 - a*ss*xm*xx*4.0D0 - a**2*cx2*c&
            &y2*xm**3*4.0D0 + a**2*cx2*cy2*xx**3*8.0D0 + a**2*cx*sy*xm**3*2.0D0 - a**&
            &2*cx*sy*xx**3*4.0D0 - a**2*cc*ss*ym**3*4.0D0 + a**2*cc*ss*yy**3*4.0D0 -&
            &a**2*cy*sx*ym**3*2.0D0 + a**2*cy*sx*yy**3*2.0D0 + a**2*cx2*xm*xx**2*8.&
            &0D0 - a**2*cx2*xm**2*xx*6.0D0 + a**2*cx2*xm*ym**2*2.0D0 - a**2*cx2*xx*ym&
            &**2*2.0D0 + a**2*cx2*xm*yy**2*2.0D0 - a**2*cx2*xx*yy**2*2.0D0 - a**2*cx2&
            &*cy2*xm*xx**2*1.6D+1 + a**2*cx2*cy2*xm**2*xx*1.2D+1 - a**2*cx2*cy2*xm*&
            &ym**2*4.0D0 + a**2*cx2*cy2*xx*ym**2*4.0D0 - a**2*cx2*cy2*xm*yy**2*4.0D&
            &0 + a**2*cx2*cy2*xx*yy**2*4.0D0 + a**2*cx*sy*xm*xx**2*8.0D0 - a**2*cx*sy&
            &*xm**2*xx*6.0D0 - a**2*cc*ss*xm**2*ym*4.0D0 - a**2*cc*ss*xx**2*ym*8.0D&
            &0 + a**2*cc*ss*xm**2*yy*4.0D0 + a**2*cc*ss*xx**2*yy*8.0D0 + a**2*cx*sy*x&
            &m*ym**2*2.0D0 - a**2*cy*sx*xm**2*ym*2.0D0 - a**2*cx*sy*xx*ym**2*2.0D0 -&
            &a**2*cy*sx*xx**2*ym*4.0D0 + a**2*cx*sy*xm*yy**2*2.0D0 + a**2*cy*sx*xm*&
            &*2*yy*2.0D0 - a**2*cx*sy*xx*yy**2*2.0D0 + a**2*cy*sx*xx**2*yy*4.0D0 - a*&
            &*2*cc*ss*ym*yy**2*1.2D+1 + a**2*cc*ss*ym**2*yy*1.2D+1 - a**2*cy*sx*ym*&
            &yy**2*6.0D0 + a**2*cy*sx*ym**2*yy*6.0D0 - a*cx*cy*xm*ym*2.0D0 + a*cx*cy*&
            &xx*ym*2.0D0 + a*cx*cy*xm*yy*2.0D0 - a*cx*cy*xx*yy*2.0D0 - a*cx*cy2*sx*xm&
            &**2*2.0D0 - a*cx*cy2*sx*xx**2*4.0D0 - a**2*cx2*xm*ym*yy*4.0D0 + a**2*cx2&
            &*xx*ym*yy*4.0D0 + a**2*cx2*cy2*xm*ym*yy*8.0D0 - a**2*cx2*cy2*xx*ym*yy*&
            &8.0D0 + a**2*cc*ss*xm*xx*ym*8.0D0 - a**2*cc*ss*xm*xx*yy*8.0D0 + a**2*cy*&
            &sx*xm*xx*ym*4.0D0 - a**2*cy*sx*xm*xx*yy*4.0D0 - a**2*cx*sy*xm*ym*yy*4.&
            &0D0 + a**2*cx*sy*xx*ym*yy*4.0D0 + a*cx*cy2*sx*xm*xx*4.0D0 - a*cx2*cy*sy*&
            &xm*ym*2.0D0 + a*cx2*cy*sy*xx*ym*2.0D0 + a*cx2*cy*sy*xm*yy*2.0D0 - a*cx2*&
            &cy*sy*xx*yy*2.0D0))/Mref + (3.0D0**(-epn - 1.0D0)*a*(ym*2.0D0 - yy*2.0D0&
            &)*(xm - xx)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(&
            &-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*(ss*ym - ss*yy&
            &+ cx*cy*xm - cx*cy*xx + cx2*cy*sy*xm - cx2*cy*sy*xx - cx*cy2*sx*ym + cx*cy2*s&
            &x*yy))/Mref + (3.0D0**(-epn - 1.0D0)*1.0D0/Mref**2*a**2*cc*epn*(cx*sy +&
            &1.0D0)*(xm - xx)*((-cx2*cy2 + cx*sy*2.0D0 + 4.0D+1)/Mref)**(epn - 1.0D0)*(&
            &ss*ym - ss*yy + cx*cy*xm - cx*cy*xx + cx2*cy*sy*xm - cx2*cy*sy*xx - cx*cy2*sx*&
            &ym + cx*cy2*sx*yy)*4.0D0)/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.&
            &0D0 + ym**2 + yy**2)) - ((a*cc*(ym - yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2&
            &.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(ss*(-1.0D+1) + cx&
            &*sy*1.0D+2 + cx2*sy2*5.0D0 - sx2*sy2*5.0D0 + cx*cy2*sx*4.0D0 - cx**3*cy2*s&
            &y + cx*cy2*sx2*sy*2.0D0))/3.0D0 - (a*cy*sx*(ym - yy)*(ss + 2.0D0)*1.0D0/sq&
            &rt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2&
            &+ yy**2))*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/3.0D0 + (cc*1.0D0/xx**3*(ym -&
            &yy)*(ss + 2.0D0)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**&
            &2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(-cx2*cy2 + cx*sy*5.0D0 + 1&
            &.0D+2)*(-xm*xx - ym*yy*2.0D0 + xm**2 + ym**2 + yy**2))/3.0D0)/xx + Mref*cc*p&
            &courr*((a*sx*(xm - xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*&
            &2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*s&
            &x*2.0D0)*(2.0D0/3.0D0))/(Mref*xx) + (a*cx*(ym - yy)*1.0D0/sqrt(1.0D0/x&
            &x**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(&
            &cy - sy*5.0D0 + cy*ss)*(4.0D0/3.0D0))/(Mref*xx)) + (a*csii*1.0D0/(xm*xx*&
            &(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*(xm**2*xx**&
            &2*1.1D+1 + xm**2*ym**2 + xx**2*ym**2*5.0D0 + xm**2*yy**2 + xx**2*yy**2*5.0&
            &D0 - cx2*xm**4*2.0D0 - cy2*xm**4 - cx2*xx**4*1.2D+1 - cy2*xx**4*6.0D0 + ss*x&
            &m**4*2.0D0 + ss*xx**4*1.2D+1 - xm*xx**3*1.2D+1 - xm**3*xx*5.0D0 + xm**4 + xx&
            &**4*6.0D0 + a*ss*xx**5*1.2D+2 - cx*sy*xm**4*2.0D+1 - cx*sy*xx**4*1.2D+2 +&
            &cx2*xm*xx**3*2.4D+1 + cx2*xm**3*xx*1.0D+1 + cy2*xm*xx**3*1.2D+1 + cy2*xm&
            &**3*xx*5.0D0 - ss*xm*xx**3*2.4D+1 - ss*xm**3*xx*1.0D+1 - xm*xx*ym**2*3.0&
            &D0 - xm*xx*yy**2*3.0D0 - xm**2*ym*yy*2.0D0 - xx**2*ym*yy*1.0D+1 + a*xx**5*&
            &sin(a*xx*2.0D0)*1.0D+1 - cx2*xm**2*xx**2*2.2D+1 - cy2*xm**2*xx**2*1.1D&
            &+1 - cx2*xm**2*ym**2*2.0D0 - cy2*xm**2*ym**2 - cx2*xx**2*ym**2*1.0D+1 - cx&
            &2*xm**2*yy**2*2.0D0 - cy2*xx**2*ym**2*5.0D0 - cy2*xm**2*yy**2 - cx2*xx**&
            &2*yy**2*1.0D+1 - cy2*xx**2*yy**2*5.0D0 + ss*xm**2*xx**2*2.2D+1 + ss*xm**&
            &2*ym**2*2.0D0 + ss*xx**2*ym**2*1.0D+1 + ss*xm**2*yy**2*2.0D0 + ss*xx**2*&
            &yy**2*1.0D+1 + cx2*cy2*xm**4*2.0D0 + cx2*cy2*xx**4*1.2D+1 + a*xm**2*xx**&
            &3*sin(a*xx*2.0D0)*1.7D+1 - a*xm**3*xx**2*sin(a*xx*2.0D0)*8.0D0 + a*xx*&
            &*3*ym**2*sin(a*xx*2.0D0)*7.0D0 + a*xx**3*yy**2*sin(a*xx*2.0D0)*7.0D0&
            &- a*xx**2*ym**3*sin(a*yy*2.0D0)*2.0D0 + a*xx**2*yy**3*sin(a*yy*2.0D0)&
            &*2.0D0 + xm*xx*ym*yy*6.0D0 + a*cx*sy*xx**5*1.2D+1 - cx2*cy2*xm*xx**3*2.4&
            &D+1 - cx2*cy2*xm**3*xx*1.0D+1 - cx*cy*xm*ym**3*2.0D0 - cx*cy*xm**3*ym*2.&
            &0D0 + cx*cy*xx*ym**3*4.0D0 + cx*cy*xx**3*ym*4.0D0 + cx*cy*xm*yy**3*2.0D0&
            &+ cx*cy*xm**3*yy*2.0D0 - cx*cy*xx*yy**3*4.0D0 - cx*cy*xx**3*yy*4.0D0 - a*&
            &ss*xm*xx**4*2.0D+2 + a*ss*xm**4*xx*2.0D+1 + a*ss*xx*ym**4*2.0D+1 + a*ss*&
            &xx*yy**4*2.0D+1 + cx*sy*xm*xx**3*2.4D+2 + cx*sy*xm**3*xx*1.0D+2 - cc*ss*&
            &xm*ym**3*2.0D0 - cc*ss*xm**3*ym*2.0D0 + cc*ss*xx*ym**3*4.0D0 + cc*ss*xx*&
            &*3*ym*4.0D0 + cc*ss*xm*yy**3*2.0D0 + cc*ss*xm**3*yy*2.0D0 - cc*ss*xx*yy*&
            &*3*4.0D0 - cc*ss*xx**3*yy*4.0D0 - cy*sx*xm*ym**3*2.0D+1 - cy*sx*xm**3*ym&
            &*2.0D+1 + cy*sx*xx*ym**3*4.0D+1 + cy*sx*xx**3*ym*4.0D+1 + cy*sx*xm*yy**3&
            &*2.0D+1 + cy*sx*xm**3*yy*2.0D+1 - cy*sx*xx*yy**3*4.0D+1 - cy*sx*xx**3*yy&
            &*4.0D+1 + cx2*xm*xx*ym**2*6.0D0 + cy2*xm*xx*ym**2*3.0D0 + cx2*xm*xx*yy**&
            &2*6.0D0 + cy2*xm*xx*yy**2*3.0D0 + cx2*xm**2*ym*yy*4.0D0 + cy2*xm**2*ym*y&
            &y*2.0D0 + cx2*xx**2*ym*yy*2.0D+1 + cy2*xx**2*ym*yy*1.0D+1 - ss*xm*xx*ym*&
            &*2*6.0D0 - ss*xm*xx*yy**2*6.0D0 - ss*xm**2*ym*yy*4.0D0 - ss*xx**2*ym*yy*&
            &2.0D+1 - a*xm*xx**4*sin(a*xx*2.0D0)*1.8D+1 + a*xm**4*xx*sin(a*xx*2.0D0&
            &)*2.0D0 + a*xx*ym**4*sin(a*xx*2.0D0) + a*xx*yy**4*sin(a*xx*2.0D0) - a*xx&
            &**4*ym*sin(a*yy*2.0D0)*4.0D0 + a*xx**4*yy*sin(a*yy*2.0D0)*4.0D0 + cx2*&
            &cy2*xm**2*xx**2*2.2D+1 + cx2*cy2*xm**2*ym**2*2.0D0 + cx2*cy2*xx**2*ym*&
            &*2*1.0D+1 + cx2*cy2*xm**2*yy**2*2.0D0 + cx2*cy2*xx**2*yy**2*1.0D+1 + a*s&
            &s*xm**2*xx**3*1.8D+2 - a*ss*xm**3*xx**2*8.0D+1 + a*ss*xx**3*ym**2*1.0D&
            &+2 + a*ss*xx**3*yy**2*1.0D+2 - cx*sy*xm**2*xx**2*2.2D+2 - cx*sy*xm**2*ym&
            &**2*2.0D+1 - cx*sy*xx**2*ym**2*1.0D+2 - cx*sy*xm**2*yy**2*2.0D+1 - cx*sy&
            &*xx**2*yy**2*1.0D+2 - cc*ss*xm*xx**2*ym*8.0D0 + cc*ss*xm**2*xx*ym*8.0D&
            &0 + cc*ss*xm*xx**2*yy*8.0D0 - cc*ss*xm**2*xx*yy*8.0D0 + cx*sy*xm*xx*ym**&
            &2*6.0D+1 - cy*sx*xm*xx**2*ym*8.0D+1 + cy*sx*xm**2*xx*ym*8.0D+1 + cx*sy*x&
            &m*xx*yy**2*6.0D+1 + cy*sx*xm*xx**2*yy*8.0D+1 - cy*sx*xm**2*xx*yy*8.0D+&
            &1 - cc*ss*xm*ym*yy**2*6.0D0 + cc*ss*xm*ym**2*yy*6.0D0 + cc*ss*xx*ym*yy**&
            &2*1.2D+1 - cc*ss*xx*ym**2*yy*1.2D+1 + cx*sy*xm**2*ym*yy*4.0D+1 - cy*sx*x&
            &m*ym*yy**2*6.0D+1 + cy*sx*xm*ym**2*yy*6.0D+1 + cx*sy*xx**2*ym*yy*2.0D+&
            &2 + cy*sx*xx*ym*yy**2*1.2D+2 - cy*sx*xx*ym**2*yy*1.2D+2 + a*xm*xx*ym**3*&
            &sin(a*yy*2.0D0)*2.0D0 + a*xm*xx**3*ym*sin(a*yy*2.0D0)*8.0D0 + a*xm**3*&
            &xx*ym*sin(a*yy*2.0D0)*2.0D0 - a*xx*ym*yy**3*sin(a*xx*2.0D0)*4.0D0 - a*&
            &xx*ym**3*yy*sin(a*xx*2.0D0)*4.0D0 - a*xx**3*ym*yy*sin(a*xx*2.0D0)*1.&
            &4D+1 - a*xm*xx*yy**3*sin(a*yy*2.0D0)*2.0D0 - a*xm*xx**3*yy*sin(a*yy*2.&
            &0D0)*8.0D0 - a*xm**3*xx*yy*sin(a*yy*2.0D0)*2.0D0 + a*cx*cy*xx**2*ym**3&
            &*4.0D+1 - a*cx*cy*xx**2*yy**3*4.0D+1 + a*cx*sy*xm**2*xx**3*1.8D+1 - a*cx&
            &*sy*xm**3*xx**2*8.0D0 + a*cx*sy*xx**3*ym**2*1.0D+1 - a*cy*sx*xx**2*ym*&
            &*3*4.0D0 + a*cx*sy*xx**3*yy**2*1.0D+1 + a*cy*sx*xx**2*yy**3*4.0D0 - a*ss&
            &*xm*xx**2*ym**2*8.0D+1 + a*ss*xm**2*xx*ym**2*4.0D+1 - a*ss*xm*xx**2*yy&
            &**2*8.0D+1 + a*ss*xm**2*xx*yy**2*4.0D+1 + a*ss*xx*ym**2*yy**2*1.2D+2 - a&
            &*xm*xx**2*ym**2*sin(a*xx*2.0D0)*6.0D0 + a*xm**2*xx*ym**2*sin(a*xx*2.&
            &0D0)*3.0D0 - a*xm*xx**2*yy**2*sin(a*xx*2.0D0)*6.0D0 + a*xm**2*xx*yy**2&
            &*sin(a*xx*2.0D0)*3.0D0 - a*xm**2*xx**2*ym*sin(a*yy*2.0D0)*6.0D0 + a*xx&
            &*ym**2*yy**2*sin(a*xx*2.0D0)*6.0D0 + a*xm**2*xx**2*yy*sin(a*yy*2.0D0&
            &)*6.0D0 - a*xx**2*ym*yy**2*sin(a*yy*2.0D0)*6.0D0 + a*xx**2*ym**2*yy*si&
            &n(a*yy*2.0D0)*6.0D0 - cx2*xm*xx*ym*yy*1.2D+1 - cy2*xm*xx*ym*yy*6.0D0 + s&
            &s*xm*xx*ym*yy*1.2D+1 - a*cx*cy2*sx*xx**5*2.4D+1 + a*cx*cy*xx**4*ym*8.0&
            &D+1 - a*cx*cy*xx**4*yy*8.0D+1 - a*cx*sy*xm*xx**4*2.0D+1 + a*cx*sy*xm**4*&
            &xx*2.0D0 + a*cx*sy*xx*ym**4*2.0D0 - a*cy*sx*xx**4*ym*8.0D0 + a*cx*sy*xx*&
            &yy**4*2.0D0 + a*cy*sx*xx**4*yy*8.0D0 - cx*cy*xm*xx**2*ym*8.0D0 + cx*cy*x&
            &m**2*xx*ym*8.0D0 - cx2*cy2*xm*xx*ym**2*6.0D0 + cx*cy*xm*xx**2*yy*8.0D0&
            &- cx*cy*xm**2*xx*yy*8.0D0 - cx2*cy2*xm*xx*yy**2*6.0D0 - cx*cy*xm*ym*yy*&
            &*2*6.0D0 + cx*cy*xm*ym**2*yy*6.0D0 - cx2*cy2*xm**2*ym*yy*4.0D0 + cx*cy*x&
            &x*ym*yy**2*1.2D+1 - cx*cy*xx*ym**2*yy*1.2D+1 - cx2*cy2*xx**2*ym*yy*2.0&
            &D+1 - a*ss*xx*ym*yy**3*8.0D+1 - a*ss*xx*ym**3*yy*8.0D+1 - a*ss*xx**3*ym*&
            &yy*2.0D+2 + cx2*cy2*xm*xx*ym*yy*1.2D+1 - cx*sy*xm*xx*ym*yy*1.2D+2 + a*cx&
            &*cy2*sx*xm*xx**4*4.0D+1 - a*cx*cy2*sx*xm**4*xx*4.0D0 - a*cx*cy2*sx*xx*&
            &ym**4*4.0D0 + a*cx2*cy*sy*xx**4*ym*1.6D+1 - a*cx*cy2*sx*xx*yy**4*4.0D0&
            &- a*cx2*cy*sy*xx**4*yy*1.6D+1 - a*cx*cy*xm*xx*ym**3*4.0D+1 - a*cx*cy*xm&
            &*xx**3*ym*1.6D+2 - a*cx*cy*xm**3*xx*ym*4.0D+1 + a*cx*cy*xm*xx*yy**3*4.&
            &0D+1 + a*cx*cy*xm*xx**3*yy*1.6D+2 + a*cx*cy*xm**3*xx*yy*4.0D+1 + a*cy*sx&
            &*xm*xx*ym**3*4.0D0 + a*cy*sx*xm*xx**3*ym*1.6D+1 + a*cy*sx*xm**3*xx*ym*&
            &4.0D0 - a*cy*sx*xm*xx*yy**3*4.0D0 - a*cy*sx*xm*xx**3*yy*1.6D+1 - a*cy*sx&
            &*xm**3*xx*yy*4.0D0 - a*cx*sy*xx*ym*yy**3*8.0D0 - a*cx*sy*xx*ym**3*yy*8&
            &.0D0 - a*cx*sy*xx**3*ym*yy*2.0D+1 + a*ss*xm*xx**2*ym*yy*1.6D+2 - a*ss*xm&
            &**2*xx*ym*yy*8.0D+1 - a*cx*cy2*sx*xm**2*xx**3*3.6D+1 + a*cx*cy2*sx*xm*&
            &*3*xx**2*1.6D+1 - a*cx*cy2*sx*xx**3*ym**2*2.0D+1 + a*cx2*cy*sy*xx**2*y&
            &m**3*8.0D0 - a*cx*cy2*sx*xx**3*yy**2*2.0D+1 - a*cx2*cy*sy*xx**2*yy**3*&
            &8.0D0 + a*xm*xx**2*ym*yy*sin(a*xx*2.0D0)*1.2D+1 - a*xm**2*xx*ym*yy*sin&
            &(a*xx*2.0D0)*6.0D0 + a*cx*cy*xm**2*xx**2*ym*1.2D+2 - a*cx*cy*xm**2*xx*&
            &*2*yy*1.2D+2 + a*xm*xx*ym*yy**2*sin(a*yy*2.0D0)*6.0D0 - a*xm*xx*ym**2*&
            &yy*sin(a*yy*2.0D0)*6.0D0 + a*cx*cy*xx**2*ym*yy**2*1.2D+2 - a*cx*cy*xx*&
            &*2*ym**2*yy*1.2D+2 - a*cx*sy*xm*xx**2*ym**2*8.0D0 + a*cx*sy*xm**2*xx*y&
            &m**2*4.0D0 - a*cy*sx*xm**2*xx**2*ym*1.2D+1 - a*cx*sy*xm*xx**2*yy**2*8.&
            &0D0 + a*cx*sy*xm**2*xx*yy**2*4.0D0 + a*cy*sx*xm**2*xx**2*yy*1.2D+1 + a*c&
            &x*sy*xx*ym**2*yy**2*1.2D+1 - a*cy*sx*xx**2*ym*yy**2*1.2D+1 + a*cy*sx*x&
            &x**2*ym**2*yy*1.2D+1 + a*cx*cy2*sx*xm*xx**2*ym**2*1.6D+1 - a*cx*cy2*sx&
            &*xm**2*xx*ym**2*8.0D0 + a*cx2*cy*sy*xm**2*xx**2*ym*2.4D+1 + a*cx*cy2*s&
            &x*xm*xx**2*yy**2*1.6D+1 - a*cx*cy2*sx*xm**2*xx*yy**2*8.0D0 - a*cx2*cy*&
            &sy*xm**2*xx**2*yy*2.4D+1 - a*cx*cy2*sx*xx*ym**2*yy**2*2.4D+1 + a*cx2*c&
            &y*sy*xx**2*ym*yy**2*2.4D+1 - a*cx2*cy*sy*xx**2*ym**2*yy*2.4D+1 - a*cx2&
            &*cy*sy*xm*xx*ym**3*8.0D0 - a*cx2*cy*sy*xm*xx**3*ym*3.2D+1 - a*cx2*cy*s&
            &y*xm**3*xx*ym*8.0D0 + a*cx2*cy*sy*xm*xx*yy**3*8.0D0 + a*cx2*cy*sy*xm*x&
            &x**3*yy*3.2D+1 + a*cx2*cy*sy*xm**3*xx*yy*8.0D0 + a*cx*cy2*sx*xx*ym*yy*&
            &*3*1.6D+1 + a*cx*cy2*sx*xx*ym**3*yy*1.6D+1 + a*cx*cy2*sx*xx**3*ym*yy*4&
            &.0D+1 - a*cx*cy*xm*xx*ym*yy**2*1.2D+2 + a*cx*cy*xm*xx*ym**2*yy*1.2D+2 +&
            &a*cx*sy*xm*xx**2*ym*yy*1.6D+1 - a*cx*sy*xm**2*xx*ym*yy*8.0D0 + a*cy*sx&
            &*xm*xx*ym*yy**2*1.2D+1 - a*cy*sx*xm*xx*ym**2*yy*1.2D+1 - a*cx*cy2*sx*x&
            &m*xx**2*ym*yy*3.2D+1 + a*cx*cy2*sx*xm**2*xx*ym*yy*1.6D+1 - a*cx2*cy*sy&
            &*xm*xx*ym*yy**2*2.4D+1 + a*cx2*cy*sy*xm*xx*ym**2*yy*2.4D+1))/xx + (sqr&
            &t(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0 - 2.0D+1)/Mref)*(ss + 2.0D0)**2*(-cx&
            &2*cy2 + cx*sy*2.0D0 + cy*sx*2.0D0 + 2.0D+1))/(tie*(cy*sx - 1.0D+1)*2.0D0) +&
            &(cc*1.0D0/xx**3*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*(ss + 2.0D0)*1.0D0/(1.0D&
            &0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &)**(3.0D0/2.0D0)*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/6.0D0 + (a*cx*cy2*(x&
            &m - xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx*&
            &*2*2.0D0 + ym**2 + yy**2))*(cx*1.0D+1 + sx*1.0D+2 + cx*ss*1.0D+1 + cx2*sx*2.&
            &0D0 + cx2*sy*4.0D0 - cx2*cy2*sx*3.0D0))/(xx*3.0D0) - (a*cx*sy*(xm - xx)*(s&
            &s + 2.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 +&
            &xx**2*2.0D0 + ym**2 + yy**2))*(-cx2*cy2 + cx*sy*5.0D0 + 1.0D+2))/(xx*3.0D0&
            &)

          f(ind, 4) = (a*cx2*cy*(ym - yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*y&
            &y*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(cy - sy*5.0D0 + cy*ss)*(1.0D+&
            &1/3.0D0) - a*cy*sx*(cy*sx - 1.0D+1)*(ym - yy)*(ss + 2.0D0)*1.0D0/sqrt(1.0D&
            &0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &)*(5.0D0/3.0D0) + cc*1.0D0/xx**3*(cy*sx - 1.0D+1)*(ym - yy)*(ss + 2.0D0)*1&
            &.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym&
            &**2 + yy**2))**(3.0D0/2.0D0)*(-xm*xx - ym*yy*2.0D0 + xm**2 + ym**2 + yy**2)*&
            &(5.0D0/3.0D0))/xx + (1.0D0/Mref**2*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0&
            &+ xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*(3.0D0**(-epn + 1.0D0)*Mref*a*cc*&
            &kpare*ym**4*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1&
            &.0D0)*Mref*a*cc*kpare*yy**4*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D&
            &0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*ym*yy**3*(-(cy*sx*2.0D0 - 2.0D&
            &+1)/Mref)**epn*1.6D+1 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*ym**3*yy&
            &*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*1.6D+1 + 3.0D0**(-epn + 1.0D0)*Mref&
            &*a*kpare*ss*xm*ym**3*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0&
            &**(-epn + 1.0D0)*Mref*a*kpare*ss*xm**3*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mre&
            &f)**epn*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*kpare*ss*xx*ym**3*(-(cy*s&
            &x*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*kpare*&
            &ss*xx**3*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 1&
            &.0D0)*Mref*a*kpare*ss*xm*yy**3*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4&
            &.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*kpare*ss*xm**3*yy*(-(cy*sx*2.0D0 - 2&
            &.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*kpare*ss*xx*yy*&
            &*3*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mre&
            &f*a*kpare*ss*xx**3*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D&
            &0**(-epn + 1.0D0)*Mref*a*cc*kpare*xx**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)&
            &**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*2&
            &.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*ym**2*(-(cy*sx*2.0D0 - 2.0D&
            &+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2&
            &+ yy**2)*2.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*yy**2*(-(cy*sx*2&
            &.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.&
            &0D0 + ym**2 + yy**2)*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*xm**2*y&
            &m**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*M&
            &ref*a*cc*kpare*xx**2*ym**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0&
            &+ 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*xm**2*yy**2*(-(cy*sx*2.0D0 - 2.&
            &0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*xx**2*y&
            &y**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 2.0D0)*M&
            &ref*a*cc*kpare*ym**2*yy**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*8.0D0&
            &- 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpare*sx*xx**3*(-(cy*sx*2.0D0 - 2.&
            &0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym*&
            &*2 + yy**2)*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*xm*xx*(-(cy*sx&
            &*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*&
            &2.0D0 + ym**2 + yy**2)*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*ym*yy&
            &*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm*&
            &*2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*kpare&
            &*ss*xm*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*&
            &2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mre&
            &f*a*kpare*ss*xm*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D&
            &0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*2.0D0 - 3.0D0**(-epn + 1&
            &.0D0)*Mref*a*cc*kpare*xm*xx*ym**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**ep&
            &n*8.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*xm*xx*yy**2*(-(cy*sx*2&
            &.0D0 - 2.0D+1)/Mref)**epn*8.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a*cc*kpare*&
            &xm**2*ym*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*8.0D0 - 3.0D0**(-epn + 1&
            &.0D0)*Mref*a*cc*kpare*xx**2*ym*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**ep&
            &n*8.0D0 + 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*xm*xx**2*ym*(-(cy*sx*2&
            &.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*&
            &xm**2*xx*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 2&
            &.0D0)*Mref*a*kpare*ss*xm*xx**2*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**ep&
            &n*4.0D0 + 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*xm**2*xx*yy*(-(cy*sx*2&
            &.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*&
            &xm*ym*yy**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 - 3.0D0**(-epn + 2&
            &.0D0)*Mref*a*kpare*ss*xm*ym**2*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**ep&
            &n*4.0D0 - 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*xx*ym*yy**2*(-(cy*sx*2&
            &.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 2.0D0)*Mref*a*kpare*ss*&
            &xx*ym**2*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*4.0D0 + 3.0D0**(-epn + 1&
            &.0D0)*Mref*a*cc*kpare*xm*xx*ym*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**ep&
            &n*1.6D+1 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpare*sx*xm*xx**2*(-(cy*&
            &sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**&
            &2*2.0D0 + ym**2 + yy**2)*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpare*&
            &sx*xm**2*xx*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*y&
            &y*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*2.0D0 - 3.0D0**(-epn + 1.0D0)*M&
            &ref*a**2*cx*kpare*sy*xx**2*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(x&
            &m*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0 - 3.0&
            &D0**(-epn + 1.0D0)*Mref*a**2*cy*kpare*sx*xx*ym**2*(-(cy*sx*2.0D0 - 2.0&
            &D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**&
            &2 + yy**2)*2.0D0 + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpare*sy*xx**2*yy*&
            &(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**&
            &2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0 - 3.0D0**(-epn + 1.0D0)*Mref*a**2*cy*&
            &kpare*sx*xx*yy**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0&
            &) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*2.0D0 + (3.0D0**(-epn + 1&
            &.0D0)*Mref*a**2*epn*kpare*sx2*sy2*xx**3*(-(cy*sx*2.0D0 - 2.0D+1)/Mre&
            &f)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &*4.0D0)/(cy*sx*2.0D0 - 2.0D+1) + 3.0D0**(-epn + 1.0D0)*Mref*a**2*cx*kpar&
            &e*sy*xm*xx*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym&
            &*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0 - 3.0D0**(-epn + 1.0D0)&
            &*Mref*a**2*cx*kpare*sy*xm*xx*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*&
            &(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0 + 3&
            &.0D0**(-epn + 1.0D0)*Mref*a**2*cy*kpare*sx*xx*ym*yy*(-(cy*sx*2.0D0 - 2&
            &.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym&
            &**2 + yy**2)*4.0D0 - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cc*epn*kpare*ss*xx&
            &**2*ym*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0&
            &D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*8.0D0)/(cy*sx*2.0D0 - 2.0D+1) + (3.0&
            &D0**(-epn + 1.0D0)*Mref*a**2*cc*epn*kpare*ss*xx**2*yy*(-(cy*sx*2.0D0&
            &- 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 +&
            &ym**2 + yy**2)*8.0D0)/(cy*sx*2.0D0 - 2.0D+1) - (3.0D0**(-epn + 1.0D0)*Mref&
            &*a**2*epn*kpare*sx2*sy2*xm*xx**2*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn&
            &*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*8.0D0)&
            &/(cy*sx*2.0D0 - 2.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*a**2*epn*kpare*sx2&
            &*sy2*xm**2*xx*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym&
            &*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*4.0D0)/(cy*sx*2.0D0 - 2.0D+&
            &1) + (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*cy2*epn*kpare*xx*ym**2*(-(cy&
            &*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx*&
            &*2*2.0D0 + ym**2 + yy**2)*4.0D0)/(cy*sx*2.0D0 - 2.0D+1) + (3.0D0**(-epn + 1.&
            &0D0)*Mref*a**2*cx2*cy2*epn*kpare*xx*yy**2*(-(cy*sx*2.0D0 - 2.0D+1)/M&
            &ref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**&
            &2)*4.0D0)/(cy*sx*2.0D0 - 2.0D+1) - (3.0D0**(-epn + 1.0D0)*Mref*a**2*cx2*&
            &cy2*epn*kpare*xx*ym*yy*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-&
            &2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)*8.0D0)/(cy*sx*2.&
            &0D0 - 2.0D+1) + (3.0D0**(-epn + 1.0D0)*Mref*a**2*cc*epn*kpare*ss*xm*xx*y&
            &m*(-(cy*sx*2.0D0 - 2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm&
            &**2 + xx**2*2.0D0 + ym**2 + yy**2)*8.0D0)/(cy*sx*2.0D0 - 2.0D+1) - (3.0D0**(&
            &-epn + 1.0D0)*Mref*a**2*cc*epn*kpare*ss*xm*xx*yy*(-(cy*sx*2.0D0 - 2.0D&
            &+1)/Mref)**epn*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2&
            &+ yy**2)*8.0D0)/(cy*sx*2.0D0 - 2.0D+1)))/(xx*9.0D0) - Mref*cc*pcourr*((&
            &a*sx*(xm - xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm&
            &**2 + xx**2*2.0D0 + ym**2 + yy**2))*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*sx*2.0D0)&
            &*(2.0D0/3.0D0))/(Mref*xx) + (a*cx*(ym - yy)*1.0D0/sqrt(1.0D0/xx**2*(xm&
            &*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(cy - sy*5.&
            &0D0 + cy*ss)*(4.0D0/3.0D0))/(Mref*xx)) + (a*csie*1.0D0/(xm*xx*(-2.0D0)&
            &- ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*(a*ss*xx**5*6.0D+1 +&
            &cc*ss*xm**4*2.0D0 + cc*ss*xx**4*1.2D+1 - cx*sy*xm**4*1.0D+1 - cx*sy*xx**&
            &4*6.0D+1 - ss*xm*ym**3*2.0D0 - ss*xm**3*ym*2.0D0 + ss*xx*ym**3*4.0D0 + ss*&
            &xx**3*ym*4.0D0 + ss*xm*yy**3*2.0D0 + ss*xm**3*yy*2.0D0 - ss*xx*yy**3*4.0&
            &D0 - ss*xx**3*yy*4.0D0 + cx*cy*xm**4*2.0D0 + cx*cy*xx**4*1.2D+1 - a*cy*sx*&
            &xx**5*1.2D+1 - cx*cy*xm*xx**3*2.4D+1 - cx*cy*xm**3*xx*1.0D+1 - a*ss*xm*x&
            &x**4*1.0D+2 + a*ss*xm**4*xx*1.0D+1 + a*ss*xx*ym**4*1.0D+1 + a*ss*xx*yy**&
            &4*1.0D+1 - cc*ss*xm*xx**3*2.4D+1 - cc*ss*xm**3*xx*1.0D+1 + cx*sy*xm*xx**&
            &3*1.2D+2 + cx*sy*xm**3*xx*5.0D+1 - cy*sx*xm*ym**3*1.0D+1 - cy*sx*xm**3*y&
            &m*1.0D+1 + cy2*sx2*xm*ym**3 + cy2*sx2*xm**3*ym + cy*sx*xx*ym**3*2.0D+1 + c&
            &y*sx*xx**3*ym*2.0D+1 - cy2*sx2*xx*ym**3*2.0D0 - cy2*sx2*xx**3*ym*2.0D0&
            &+ cy*sx*xm*yy**3*1.0D+1 + cy*sx*xm**3*yy*1.0D+1 - cy2*sx2*xm*yy**3 - cy2*&
            &sx2*xm**3*yy - cy*sx*xx*yy**3*2.0D+1 - cy*sx*xx**3*yy*2.0D+1 + cy2*sx2*x&
            &x*yy**3*2.0D0 + cy2*sx2*xx**3*yy*2.0D0 - sx2*sy2*xm*ym**3 - sx2*sy2*xm**&
            &3*ym + sx2*sy2*xx*ym**3*2.0D0 + sx2*sy2*xx**3*ym*2.0D0 + sx2*sy2*xm*yy**&
            &3 + sx2*sy2*xm**3*yy - sx2*sy2*xx*yy**3*2.0D0 - sx2*sy2*xx**3*yy*2.0D0 - s&
            &s*xm*xx**2*ym*8.0D0 + ss*xm**2*xx*ym*8.0D0 + ss*xm*xx**2*yy*8.0D0 - ss*x&
            &m**2*xx*yy*8.0D0 - ss*xm*ym*yy**2*6.0D0 + ss*xm*ym**2*yy*6.0D0 + ss*xx*y&
            &m*yy**2*1.2D+1 - ss*xx*ym**2*yy*1.2D+1 + cx*cy*xm**2*xx**2*2.2D+1 + cx*c&
            &y*xm**2*ym**2*2.0D0 + cx*cy*xx**2*ym**2*1.0D+1 + cx*cy*xm**2*yy**2*2.0&
            &D0 + cx*cy*xx**2*yy**2*1.0D+1 + a*ss*xm**2*xx**3*9.0D+1 - a*ss*xm**3*xx*&
            &*2*4.0D+1 + a*ss*xx**3*ym**2*5.0D+1 + a*ss*xx**3*yy**2*5.0D+1 + cc*ss*xm&
            &**2*xx**2*2.2D+1 - cx*sy*xm**2*xx**2*1.1D+2 + cc*ss*xm**2*ym**2*2.0D0 +&
            &cc*ss*xx**2*ym**2*1.0D+1 + cc*ss*xm**2*yy**2*2.0D0 + cc*ss*xx**2*yy**2&
            &*1.0D+1 - cx*sy*xm**2*ym**2*1.0D+1 - cx*sy*xx**2*ym**2*5.0D+1 - cx*sy*xm&
            &**2*yy**2*1.0D+1 - cx*sy*xx**2*yy**2*5.0D+1 - cc*ss*xm*xx*ym**2*6.0D0 -&
            &cc*ss*xm*xx*yy**2*6.0D0 + cx*sy*xm*xx*ym**2*3.0D+1 - cy*sx*xm*xx**2*ym&
            &*4.0D+1 + cy*sx*xm**2*xx*ym*4.0D+1 + cy2*sx2*xm*xx**2*ym*4.0D0 - cy2*sx2&
            &*xm**2*xx*ym*4.0D0 + cx*sy*xm*xx*yy**2*3.0D+1 + cy*sx*xm*xx**2*yy*4.0D&
            &+1 - cy*sx*xm**2*xx*yy*4.0D+1 - cy2*sx2*xm*xx**2*yy*4.0D0 + cy2*sx2*xm**&
            &2*xx*yy*4.0D0 - cc*ss*xm**2*ym*yy*4.0D0 - cc*ss*xx**2*ym*yy*2.0D+1 + cx*&
            &sy*xm**2*ym*yy*2.0D+1 - cy*sx*xm*ym*yy**2*3.0D+1 + cy*sx*xm*ym**2*yy*3&
            &.0D+1 + cy2*sx2*xm*ym*yy**2*3.0D0 - cy2*sx2*xm*ym**2*yy*3.0D0 + cx*sy*xx&
            &**2*ym*yy*1.0D+2 + cy*sx*xx*ym*yy**2*6.0D+1 - cy*sx*xx*ym**2*yy*6.0D+1&
            &- cy2*sx2*xx*ym*yy**2*6.0D0 + cy2*sx2*xx*ym**2*yy*6.0D0 - sx2*sy2*xm*xx&
            &**2*ym*4.0D0 + sx2*sy2*xm**2*xx*ym*4.0D0 + sx2*sy2*xm*xx**2*yy*4.0D0 - s&
            &x2*sy2*xm**2*xx*yy*4.0D0 - sx2*sy2*xm*ym*yy**2*3.0D0 + sx2*sy2*xm*ym**&
            &2*yy*3.0D0 + sx2*sy2*xx*ym*yy**2*6.0D0 - sx2*sy2*xx*ym**2*yy*6.0D0 + a*c&
            &x*cy*xx**2*ym**3*2.0D+1 - a*cx*cy*xx**2*yy**3*2.0D+1 - a*cy*sx*xm**2*x&
            &x**3*1.8D+1 + a*cy*sx*xm**3*xx**2*8.0D0 + a*cx*sy*xx**2*ym**3*4.0D0 - a*&
            &cy*sx*xx**3*ym**2*1.0D+1 - a*cx*sy*xx**2*yy**3*4.0D0 - a*cy*sx*xx**3*y&
            &y**2*1.0D+1 - a*ss*xm*xx**2*ym**2*4.0D+1 + a*ss*xm**2*xx*ym**2*2.0D+1 -&
            &a*ss*xm*xx**2*yy**2*4.0D+1 + a*ss*xm**2*xx*yy**2*2.0D+1 + a*ss*xx*ym**&
            &2*yy**2*6.0D+1 + a*cx2*cy*sy*xx**5*8.0D0 + a*cx*cy*xx**4*ym*4.0D+1 - a*c&
            &x*cy*xx**4*yy*4.0D+1 - a*cy*sx2*sy*xx**5*1.6D+1 + a*cy*sx*xm*xx**4*2.0&
            &D+1 - a*cy*sx*xm**4*xx*2.0D0 + a*cx*sy*xx**4*ym*8.0D0 - a*cy*sx*xx*ym**4&
            &*2.0D0 - a*cx*sy*xx**4*yy*8.0D0 - a*cy*sx*xx*yy**4*2.0D0 - cx*cy*xm*xx*y&
            &m**2*6.0D0 - cx*cy*xm*xx*yy**2*6.0D0 - cx*cy*xm**2*ym*yy*4.0D0 - cx*cy*x&
            &x**2*ym*yy*2.0D+1 - a*ss*xx*ym*yy**3*4.0D+1 - a*ss*xx*ym**3*yy*4.0D+1 -&
            &a*ss*xx**3*ym*yy*1.0D+2 + cx*cy*xm*xx*ym*yy*1.2D+1 + cc*ss*xm*xx*ym*yy&
            &*1.2D+1 - cx*sy*xm*xx*ym*yy*6.0D+1 - a*cx2*cy*sy*xm*xx**4*1.6D+1 + a*cx2&
            &*cy*sy*xm**4*xx*2.0D0 - a*cx*cy2*sx*xx**4*ym*8.0D0 + a*cx*cy2*sx*xx**4&
            &*yy*8.0D0 - a*cx*cy*xm*xx*ym**3*2.0D+1 - a*cx*cy*xm*xx**3*ym*8.0D+1 - a*&
            &cx*cy*xm**3*xx*ym*2.0D+1 + a*cx*cy*xm*xx*yy**3*2.0D+1 + a*cx*cy*xm*xx*&
            &*3*yy*8.0D+1 + a*cx*cy*xm**3*xx*yy*2.0D+1 + a*cy*sx2*sy*xm*xx**4*2.4D+&
            &1 - a*cy*sx2*sy*xm**4*xx*2.0D0 + a*cx*sx*sy2*xx**4*ym*8.0D0 - a*cy*sx2*s&
            &y*xx*ym**4*4.0D0 - a*cx*sx*sy2*xx**4*yy*8.0D0 - a*cy*sx2*sy*xx*yy**4*4&
            &.0D0 - a*cx*sy*xm*xx*ym**3*4.0D0 - a*cx*sy*xm*xx**3*ym*1.6D+1 - a*cx*sy*&
            &xm**3*xx*ym*4.0D0 + a*cx*sy*xm*xx*yy**3*4.0D0 + a*cx*sy*xm*xx**3*yy*1.&
            &6D+1 + a*cx*sy*xm**3*xx*yy*4.0D0 + a*cy*sx*xx*ym*yy**3*8.0D0 + a*cy*sx*x&
            &x*ym**3*yy*8.0D0 + a*cy*sx*xx**3*ym*yy*2.0D+1 + a*ss*xm*xx**2*ym*yy*8.&
            &0D+1 - a*ss*xm**2*xx*ym*yy*4.0D+1 + a*cx2*cy*sy*xm**2*xx**3*1.6D+1 - a*c&
            &x2*cy*sy*xm**3*xx**2*8.0D0 - a*cx*cy2*sx*xx**2*ym**3*4.0D0 + a*cx2*cy*&
            &sy*xx**3*ym**2*4.0D0 + a*cx*cy2*sx*xx**2*yy**3*4.0D0 + a*cx2*cy*sy*xx*&
            &*3*yy**2*4.0D0 + a*cx*cy*xm**2*xx**2*ym*6.0D+1 - a*cx*cy*xm**2*xx**2*y&
            &y*6.0D+1 + a*cx*cy*xx**2*ym*yy**2*6.0D+1 - a*cx*cy*xx**2*ym**2*yy*6.0D&
            &+1 - a*cy*sx2*sy*xm**2*xx**3*2.0D+1 + a*cy*sx2*sy*xm**3*xx**2*8.0D0 + a*&
            &cx*sx*sy2*xx**2*ym**3*4.0D0 - a*cy*sx2*sy*xx**3*ym**2*1.6D+1 - a*cx*sx&
            &*sy2*xx**2*yy**3*4.0D0 - a*cy*sx2*sy*xx**3*yy**2*1.6D+1 + a*cx*sy*xm**&
            &2*xx**2*ym*1.2D+1 + a*cy*sx*xm*xx**2*ym**2*8.0D0 - a*cy*sx*xm**2*xx*ym&
            &**2*4.0D0 - a*cx*sy*xm**2*xx**2*yy*1.2D+1 + a*cy*sx*xm*xx**2*yy**2*8.0&
            &D0 - a*cy*sx*xm**2*xx*yy**2*4.0D0 + a*cx*sy*xx**2*ym*yy**2*1.2D+1 - a*cx&
            &*sy*xx**2*ym**2*yy*1.2D+1 - a*cy*sx*xx*ym**2*yy**2*1.2D+1 - a*cx*cy2*s&
            &x*xm**2*xx**2*ym*1.2D+1 - a*cx2*cy*sy*xm*xx**2*ym**2*4.0D0 + a*cx2*cy*&
            &sy*xm**2*xx*ym**2*2.0D0 + a*cx*cy2*sx*xm**2*xx**2*yy*1.2D+1 - a*cx2*cy&
            &*sy*xm*xx**2*yy**2*4.0D0 + a*cx2*cy*sy*xm**2*xx*yy**2*2.0D0 - a*cx*cy2&
            &*sx*xx**2*ym*yy**2*1.2D+1 + a*cx*cy2*sx*xx**2*ym**2*yy*1.2D+1 + a*cx*s&
            &x*sy2*xm**2*xx**2*ym*1.2D+1 + a*cy*sx2*sy*xm*xx**2*ym**2*1.2D+1 - a*cy&
            &*sx2*sy*xm**2*xx*ym**2*6.0D0 - a*cx*sx*sy2*xm**2*xx**2*yy*1.2D+1 + a*c&
            &y*sx2*sy*xm*xx**2*yy**2*1.2D+1 - a*cy*sx2*sy*xm**2*xx*yy**2*6.0D0 + a*&
            &cx*sx*sy2*xx**2*ym*yy**2*1.2D+1 - a*cx*sx*sy2*xx**2*ym**2*yy*1.2D+1 -&
            &a*cy*sx2*sy*xx*ym**2*yy**2*2.4D+1 + a*cx*cy2*sx*xm*xx*ym**3*4.0D0 + a*&
            &cx*cy2*sx*xm*xx**3*ym*1.6D+1 + a*cx*cy2*sx*xm**3*xx*ym*4.0D0 - a*cx*cy&
            &2*sx*xm*xx*yy**3*4.0D0 - a*cx*cy2*sx*xm*xx**3*yy*1.6D+1 - a*cx*cy2*sx*&
            &xm**3*xx*yy*4.0D0 - a*cx2*cy*sy*xx**3*ym*yy*8.0D0 - a*cx*cy*xm*xx*ym*y&
            &y**2*6.0D+1 + a*cx*cy*xm*xx*ym**2*yy*6.0D+1 - a*cx*sx*sy2*xm*xx*ym**3*&
            &4.0D0 - a*cx*sx*sy2*xm*xx**3*ym*1.6D+1 - a*cx*sx*sy2*xm**3*xx*ym*4.0D0&
            &+ a*cx*sx*sy2*xm*xx*yy**3*4.0D0 + a*cx*sx*sy2*xm*xx**3*yy*1.6D+1 + a*cx&
            &*sx*sy2*xm**3*xx*yy*4.0D0 + a*cy*sx2*sy*xx*ym*yy**3*1.6D+1 + a*cy*sx2*&
            &sy*xx*ym**3*yy*1.6D+1 + a*cy*sx2*sy*xx**3*ym*yy*3.2D+1 - a*cx*sy*xm*xx&
            &*ym*yy**2*1.2D+1 + a*cx*sy*xm*xx*ym**2*yy*1.2D+1 - a*cy*sx*xm*xx**2*ym&
            &*yy*1.6D+1 + a*cy*sx*xm**2*xx*ym*yy*8.0D0 + a*cx*cy2*sx*xm*xx*ym*yy**2&
            &*1.2D+1 - a*cx*cy2*sx*xm*xx*ym**2*yy*1.2D+1 + a*cx2*cy*sy*xm*xx**2*ym*&
            &yy*8.0D0 - a*cx2*cy*sy*xm**2*xx*ym*yy*4.0D0 - a*cx*sx*sy2*xm*xx*ym*yy*&
            &*2*1.2D+1 + a*cx*sx*sy2*xm*xx*ym**2*yy*1.2D+1 - a*cy*sx2*sy*xm*xx**2*y&
            &m*yy*2.4D+1 + a*cy*sx2*sy*xm**2*xx*ym*yy*1.2D+1))/xx - (sqrt(3.0D0)*1.&
            &0D0/sqrt(-(cy*sx*2.0D0 - 2.0D+1)/Mref)*(ss + 2.0D0)**2*(-cx2*cy2 + cx*sy&
            &*2.0D0 + cy*sx*2.0D0 + 2.0D+1))/(tie*(cy*sx - 1.0D+1)*2.0D0) + (a*cc*sx*(x&
            &m - xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx*&
            &*2*2.0D0 + ym**2 + yy**2))*(cy*1.0D+1 + sx + sy*2.0D0 - cy2*sx*2.0D0)*(5.0D0&
            &/3.0D0))/xx - cc*1.0D0/xx**3*(ym*2.0D0 - yy*2.0D0)*(cy*sx - 1.0D+1)*(xm -&
            &xx)*(ss + 2.0D0)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**&
            &2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(5.0D0/6.0D0) + (a*cx*sy*&
            &(cy*sx - 1.0D+1)*(xm - xx)*(ss + 2.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-&
            &2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(5.0D0/3.0D0))/&
            &xx

        CASE (3)
          ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
          IF (.not. switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          a = 2*pi
          b = 2*pi

          D = phys%diff_n
          mu = phys%diff_u
          csii = phys%diff_e
          csie = phys%diff_ee
          kpari = phys%diff_pari
          kpare = phys%diff_pare
          Mref = phys%Mref
          epn = phys%epn
          tie = phys%tie
          pcourr = 1.

          st = sin(a*tt)
          ct = cos(a*tt)
          r = (xm**2 - 2*ym*yy - 2*xm*xx + 2*xx**2 + ym**2 + yy**2)

          f(ind, 1) = D*((a**2*1.0D0/xx**2*sin(a*tt)*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm&
            &**2 + xx**2 + ym**2 + yy**2))/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.&
            &0D0 + ym**2 + yy**2) + a*cos(a*tt)*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(xm&
            &*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2 - (a*cos(&
            &a*tt)*(ym - yy)*(ym*yy*(-2.0D0) + xm**2 - xx**2*2.0D0 + ym**2 + yy**2)*1.0D0&
            &/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2)/xx&
            &) - (a*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**&
            &2*2.0D0 + ym**2 + yy**2))*(sin(a*tt)*2.0D0 + sin(a*tt)**2*2.0D0 - 1.0D0))/&
            &xx - 1.0D0/xx**4*cos(a*tt)*(sin(a*tt) + 2.0D0)*(ym - yy)*1.0D0/(1.0D0/xx&
            &**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(&
            &3.0D0/2.0D0)*(-xm*xx - ym*yy*2.0D0 + xm**2 + ym**2 + yy**2) + (1.0D0/xx**3*c&
            &os(a*tt)*(sin(a*tt) + 2.0D0)*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(1.0D&
            &0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &)**(3.0D0/2.0D0))/2.0D0

          f(ind, 2) = mu*(-a*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*(sin(a*tt)*2.0D0 + sin(a*tt)&
            &**2*2.0D0 - 1.0D0)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0&
            &D0 + ym**2 + yy**2)**2 + (a**2*1.0D0/xx**2*(cos(a*tt)*4.0D0 + cos(a*tt)*si&
            &n(a*tt)*8.0D0)*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2 + ym**2 + yy**2&
            &))/(xm*xx*(-4.0D0) - ym*yy*4.0D0 + xm**2*2.0D0 + xx**2*4.0D0 + ym**2*2.0D0&
            &+ yy**2*2.0D0) + (a*(ym - yy)*(sin(a*tt)*2.0D0 + sin(a*tt)**2*2.0D0 - 1.0D0&
            &)*(ym*yy*(-2.0D0) + xm**2 - xx**2*2.0D0 + ym**2 + yy**2)*1.0D0/(xm*xx*(-2.&
            &0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2)/xx) - (a*1.0D0/s&
            &qrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**&
            &2 + yy**2))*(cos(a*tt)*(-5.8D+1) + sin(a*tt)*4.0D0 - cos(a*tt)**2*4.0D0 +&
            &cos(a*tt)**3*3.0D0 + 2.0D0))/(xx*3.0D0) - (a*cos(a*tt)*1.0D0/sqrt(1.0D&
            &0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &)*(sin(a*tt)*4.0D0 + sin(a*tt)**2*3.0D0 - 1.0D0))/xx + (1.0D0/xx**3*cos(&
            &a*tt)**2*(sin(a*tt) + 2.0D0)*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(1.0D&
            &0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)&
            &)**(3.0D0/2.0D0))/2.0D0 - 1.0D0/xx**4*cos(a*tt)**2*(sin(a*tt) + 2.0D0)&
            &*(ym - yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**&
            &2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(-xm*xx - ym*yy*2.0D0 + xm**2 + ym*&
            &*2 + yy**2)

          f(ind, 3) = csii*(a*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*(cos(a*tt)*2.0D+1 + cos(a*t&
            &t*2.0D0) - sin(a*tt)*2.0D0)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 +&
            &xx**2*2.0D0 + ym**2 + yy**2)**2 + (a**2*1.0D0/xx**2*(cos(a*tt)*2.0D0 + sin&
            &(a*tt)*2.0D+1 + sin(a*tt*2.0D0)*2.0D0)*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + x&
            &m**2 + xx**2 + ym**2 + yy**2))/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2&
            &.0D0 + ym**2 + yy**2) - (a*(ym - yy)*(cos(a*tt)*2.0D+1 + cos(a*tt*2.0D0) - sin&
            &(a*tt)*2.0D0)*(ym*yy*(-2.0D0) + xm**2 - xx**2*2.0D0 + ym**2 + yy**2)*1.0D0&
            &/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2)/xx&
            &) - kpari*((3.0D0**(-epn - 1.0D0)*a**2*((cos(a*tt)*2.0D0 - cos(a*tt)**2 +&
            &4.0D+1)/Mref)**epn*(cos(a*tt) - 1.0D0)*(epn*(-2.0D0) + cos(a*tt)*8.2D+&
            &1 + epn*cos(a*tt)*2.0D0 + epn*cos(a*tt)**2*2.0D0 - epn*cos(a*tt)**3*2.0D&
            &0 + cos(a*tt)**2*3.0D0 - cos(a*tt)**3*2.0D0 + 4.0D+1)*2.0D0)/(Mref*(cos(&
            &a*tt)*2.0D0 - cos(a*tt)**2 + 4.0D+1)*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2&
            &+ xx**2*2.0D0 + ym**2 + yy**2)) + (3.0D0**(-epn - 1.0D0)*a*sin(a*tt)*(ym*2.&
            &0D0 - yy*2.0D0)*(xm - xx)*((cos(a*tt)*2.0D0 - cos(a*tt)**2 + 4.0D+1)/Mref)&
            &**epn*(cos(a*tt) - 1.0D0)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx&
            &**2*2.0D0 + ym**2 + yy**2)**2*2.0D0)/Mref - (3.0D0**(-epn - 1.0D0)*a*sin(a&
            &*tt)*(ym - yy)*((cos(a*tt)*2.0D0 - cos(a*tt)**2 + 4.0D+1)/Mref)**epn*(co&
            &s(a*tt) - 1.0D0)*(ym*yy*(-2.0D0) + xm**2 - xx**2*2.0D0 + ym**2 + yy**2)*1.0D&
            &0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*2.&
            &0D0)/(Mref*xx)) - (a*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.&
            &0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(cos(a*tt)*1.0D+1 + sin(a*tt)*1.&
            &94D+2 + sin(a*tt*2.0D0)*1.0D+1 + sin(a*tt)**3*6.0D0 - cos(a*tt)**2*2.03D&
            &+2 - cos(a*tt)**3*1.5D+1 + cos(a*tt)**4*4.0D0 + 1.0D+2))/(xx*3.0D0) + (sqr&
            &t(3.0D0)*(sin(a*tt) + 2.0D0)**2*1.0D0/sqrt(-(sin(a*tt)*2.0D0 - 2.0D+1)&
            &/Mref)*(-cos(a*tt*2.0D0) + sqrt(2.0D0)*sin(3.141592653589793D0/4.0D0&
            &+ a*tt)*4.0D0 + 3.9D+1))/(tie*(sin(a*tt) - 1.0D+1)*4.0D0) - (a*pcourr*cos&
            &(a*tt)**2*(sin(a*tt) - 4.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)&
            &- ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(4.0D0/3.0D0))/xx - (1.&
            &0D0/xx**4*cos(a*tt)*(sin(a*tt) + 2.0D0)*(ym - yy)*1.0D0/(1.0D0/xx**2*(&
            &xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0&
            &/2.0D0)*(cos(a*tt)*5.0D0 - cos(a*tt)**2 + 1.0D+2)*(-xm*xx - ym*yy*2.0D0 +&
            &xm**2 + ym**2 + yy**2))/3.0D0 + (1.0D0/xx**3*cos(a*tt)*(sin(a*tt) + 2.0D0)&
            &*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym&
            &*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(cos(a*tt&
            &)*5.0D0 - cos(a*tt)**2 + 1.0D+2))/6.0D0

          f(ind, 4) = csie*((a**2*1.0D0/xx**2*(sin(a*tt)*4.0D0 - sin(a*tt)**2*2.0D0 + 1&
            &.0D0)*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2 + ym**2 + yy**2)*2.0D0)/&
            &(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2) - a*cos(a&
            &*tt)*(sin(a*tt) - 4.0D0)*(ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(xm*xx*(-&
            &2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)**2*2.0D0 + (a*cos(&
            &a*tt)*(sin(a*tt) - 4.0D0)*(ym - yy)*(ym*yy*(-2.0D0) + xm**2 - xx**2*2.0D0 +&
            &ym**2 + yy**2)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + y&
            &m**2 + yy**2)**2*2.0D0)/xx) + kpare*((3.0D0**(-epn - 1.0D0)*a**2*(-(sin(&
            &a*tt)*2.0D0 - 2.0D+1)/Mref)**epn*(epn + sin(a*tt)*1.0D+1 - sin(a*tt)**2 -&
            &epn*sin(a*tt)**2)*2.0D0)/(Mref*(sin(a*tt) - 1.0D+1)*(xm*xx*(-2.0D0) -&
            &ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2)) + (3.0D0**(-epn - 1.0D0)*a&
            &*cos(a*tt)*(ym*2.0D0 - yy*2.0D0)*(-(sin(a*tt)*2.0D0 - 2.0D+1)/Mref)**e&
            &pn*(xm - xx)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym*&
            &*2 + yy**2)**2*2.0D0)/Mref - (3.0D0**(-epn - 1.0D0)*a*cos(a*tt)*(-(sin(a&
            &*tt)*2.0D0 - 2.0D+1)/Mref)**epn*(ym - yy)*(ym*yy*(-2.0D0) + xm**2 - xx**2*&
            &2.0D0 + ym**2 + yy**2)*1.0D0/(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2&
            &.0D0 + ym**2 + yy**2)**2*2.0D0)/(Mref*xx)) - (a*1.0D0/sqrt(1.0D0/xx**2*(&
            &xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))*(sin(a*&
            &tt)*2.2D+1 + sin(a*tt)**2*1.6D+1 - sin(a*tt)**3*3.0D0 - 8.0D0)*(5.0D0/3.&
            &0D0))/xx - (sqrt(3.0D0)*(sin(a*tt) + 2.0D0)**2*1.0D0/sqrt(-(sin(a*tt)*&
            &2.0D0 - 2.0D+1)/Mref)*(-cos(a*tt*2.0D0) + sqrt(2.0D0)*sin(3.1415926535&
            &89793D0/4.0D0 + a*tt)*4.0D0 + 3.9D+1))/(tie*(sin(a*tt) - 1.0D+1)*4.0D0) -&
            &1.0D0/xx**4*cos(a*tt)*(ym - yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - y&
            &m*yy*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(sin(a*t&
            &t)*(4.0D+1/3.0D0) - sin(a*tt)**2*(5.0D0/3.0D0) + 1.0D+2/3.0D0)*(-xm*xx&
            &- ym*yy*2.0D0 + xm**2 + ym**2 + yy**2) + (a*pcourr*cos(a*tt)**2*(sin(a*tt) -&
            &4.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*yy*2.0D0 + xm**2 + xx&
            &**2*2.0D0 + ym**2 + yy**2))*(4.0D0/3.0D0))/xx + (1.0D0/xx**3*cos(a*tt)*(&
            &ym*2.0D0 - yy*2.0D0)*(xm - xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0) - ym*y&
            &y*2.0D0 + xm**2 + xx**2*2.0D0 + ym**2 + yy**2))**(3.0D0/2.0D0)*(sin(a*tt)*&
            &(4.0D+1/3.0D0) - sin(a*tt)**2*(5.0D0/3.0D0) + 1.0D+2/3.0D0))/2.0D0

        CASE (5:6)
          ! Do nothing
        CASE (50:)
          ! Do nothing

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
  SUBROUTINE analytical_solution(iel,x, y, u)
    integer, intent(IN)                      :: iel
    real*8, dimension(:), intent(IN)         :: x, y
    real*8, dimension(:, :), intent(OUT)     :: u
    real*8, dimension(size(u, 1), phys%npv)  :: up
    integer                                  :: i
    real*8                                   :: a, r(size(x)), th(size(x))
    real*8                                   :: sigma,fluxel(refElPol%Nnodes2d)
    integer*4                                :: ierr
    real*8                                   :: minFlux2D, maxFlux2D
    real*8                                   :: xmax, xmin, ymax, ymin, xm, ym
    
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
        WRITE (6, *) "This is NOT an axisymmetric test case!"
        stop
      END IF
      ! Cartesian case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      up(:, 1) = 2 + sin(a*x)*sin(a*y)
      up(:, 2) = cos(a*x)*cos(a*y)
      up(:, 3) = 20 + cos(a*x)*sin(a*y)
      up(:, 4) = 10 - sin(a*x)*cos(a*y)
    CASE (2)
      ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      up(:, 1) = 2 + sin(a*x)*sin(a*y)
      up(:, 2) = cos(a*x)*cos(a*y)
      up(:, 3) = 20 + cos(a*x)*sin(a*y)
      up(:, 4) = 10 - sin(a*x)*cos(a*y)
    CASE (5:6)
      ! Cartesian case, square mesh, horizontal field
      up(:, 1) = 1.
      up(:, 3) = 18.
      up(:, 4) = 18.
    CASE (50:64)
      fluxel = phys%magnetic_flux(Mesh%T(iel,:))  
      fluxel = (fluxel - phys%Flux2Dmin)/(phys%Flux2Dmax - phys%Flux2Dmin)                                            
      sigma = 0.4
      up(:, 1) = 1.*exp(-fluxel**2/(2*sigma**2))
      up(:, 3) = 18.*exp(-fluxel**2/(2*sigma**2))
      up(:, 4) = 18.*exp(-fluxel**2/(2*sigma**2))
#ifdef NEUTRAL
      up(:,11) = 1.e-8
#endif
#ifdef NEUTRALGAMMA
      up(:,12) = 0.
#endif

      !                                      r =         sqrt ( (x*phys%lscale-geom%R0)**2 + (y*phys%lscale)**2 )
      !                                      th = atan(x-geom%R0,y)
      !                                      up(:,1) = 2.+ cos(r*2*pi)*sin(th*4*pi)
      !                                      up(:,2) = 10. + cos(r*2*pi)*sin(th*4*pi)
      !                                      up(:,3) = 36. + cos(r*2*pi)*sin(th*4*pi)
      !                                      up(:,4) = 36. + cos(r*2*pi)*sin(th*4*pi)

    CASE (65)
      up(:, 1) = 1.
      r = sqrt((x*phys%lscale - geom%R0)**2 + (y*phys%lscale - 0.75)**2)
      DO i = 1, size(x)
        IF (r(i) .le. 0.05) THEN
          up(i, 2) = 1.
        END IF
      END DO
    CASE (80:89)
    !Define an anylitical solution with a Gaussian shape respect to the normalized flux surface
      fluxel = phys%magnetic_flux(Mesh%T(iel,:))  
      fluxel = (fluxel - phys%Flux2Dmin)/(phys%Flux2Dmax - phys%Flux2Dmin)                                            
      sigma = 0.3
      up(:, 1) = 1.*exp(-fluxel**2/(2*sigma**2))
      up(:, 3) = 18.*exp(-fluxel**2/(2*sigma**2))
      up(:, 4) = 18.*exp(-fluxel**2/(2*sigma**2))
    !Define an anylitical solution with a Gaussian shape respect to geometrical center of symmetry of the domain 
    !  sigma = 1.5
    !  up(:, 1) = 1.*exp(-((x*phys%lscale - xm*phys%lscale)**2 + (y*phys%lscale - ym*phys%lscale)**2)/(2*sigma**2))
    !  up(:, 3) = 18.*exp(-((x*phys%lscale - xm*phys%lscale)**2 + (y*phys%lscale - ym*phys%lscale)**2)/(2*sigma**2))
    !  up(:, 4) = 18.*exp(-((x*phys%lscale - xm*phys%lscale)**2 + (y*phys%lscale - ym*phys%lscale)**2)/(2*sigma**2))
#ifdef NEUTRAL
      up(:,11) = 1.e-8
#endif
#ifdef NEUTRALGAMMA
      up(:,12) = 0.
#endif
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
    ux = 0.
    uy = 0.
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
      upx(:, 3) = -a*sin(a*x)*sin(a*y)
      upx(:, 4) = -a*cos(a*x)*cos(a*y)

      upy(:, 1) = a*sin(a*x)*cos(a*y)
      upy(:, 2) = -a*cos(a*x)*sin(a*y)
      upy(:, 3) = a*cos(a*x)*cos(a*y)
      upy(:, 4) = a*sin(a*x)*sin(a*y)

    CASE (2)
      ! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      upx(:, 1) = a*cos(a*x)*sin(a*y)
      upx(:, 2) = -a*sin(a*x)*cos(a*y)
      upx(:, 3) = -a*sin(a*x)*sin(a*y)
      upx(:, 4) = -a*cos(a*x)*cos(a*y)

      upy(:, 1) = a*sin(a*x)*cos(a*y)
      upy(:, 2) = -a*cos(a*x)*sin(a*y)
      upy(:, 3) = a*cos(a*x)*cos(a*y)
      upy(:, 4) = a*sin(a*x)*sin(a*y)
    CASE (5)
      ! Do nothing
    CASE (6)
      ! Do nothing
    CASE (50:)
#ifdef NEUTRAL
      upx = 1.e-12
      upy = 1.e-12
#endif
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
    ux(:, 3) = (upx(:, 1)*up(:, 3) + up(:, 1)*upx(:, 3))
    uy(:, 3) = (upy(:, 1)*up(:, 3) + up(:, 1)*upy(:, 3))
    ux(:, 4) = (upx(:, 1)*up(:, 4) + up(:, 1)*upx(:, 4))
    uy(:, 4) = (upy(:, 1)*up(:, 4) + up(:, 1)*upy(:, 4))
#ifdef NEUTRAL
    ux(:,5) = upx(:,5)
    uy(:,5) = upy(:,5)
#endif
#ifdef NEUTRALGAMMA
    ux(:, 6) = (upx(:, 11)*up(:, 12) + up(:, 11)*upx(:, 12))
    uy(:, 6) = (upy(:, 11)*up(:, 12) + up(:, 11)*upy(:, 12))
#endif
  END SUBROUTINE analytical_gradient

  !*****************************************
  ! Body forces
  !****************************************
  SUBROUTINE body_force(x, y, f)
    real*8, dimension(:), intent(IN) :: x, y
    real*8, dimension(:, :), intent(OUT) :: f
    integer                            ::  n, i
    real*8  :: a, b, xc, yc, D, mu, csii, csie, kpari, kpare, Mref, epn, tie, pcourr
    real*8  :: xmax, xmin, ymax, ymin, xm, ym, xx, yy
    real*8 :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
    real*8 :: t11, t12, t13, t14, t15, t16, t17, t18, t19, t20
    real*8 :: t21, t22, t23, t24, t25, t26, t27, t28, t29, t30
    real*8 :: t31, t32, t33, t34, t35, t36, t37, t38, t39, t40
    real*8 :: t41, t42, t43, t44, t45, t46, t47, t48, t49, t50
    real*8 :: t51, t52, t53, t54, t55, t56, t57, t58, t59, t60
    real*8 :: t61, t62, t63, t64, t65, t66, t67, t68, t69, t70
    real*8 :: cx2, sx2, cy2, sy2

    n = size(x)
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)
    f = 0.
    SELECT CASE (switch%testcase)
    CASE (1)
      ! Cartesian case, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      a = 2*pi
      b = 2*pi
      xc = 0.
      yc = 0.
      D = phys%diff_n
      mu = phys%diff_u
      csii = phys%diff_e
      csie = phys%diff_ee
      kpari = phys%diff_pari
      kpare = phys%diff_pare
      Mref = phys%Mref
      epn = phys%epn
      tie = phys%tie
      pcourr = 1.

      !  xm = -0.5
      !  ym = -0.5

      DO i = 1, size(x)
        xx = x(i)
        yy = y(i)

        cx2 = cos(a*xx)**2
        sx2 = sin(a*xx)**2
        cy2 = cos(a*yy)**2
        sy2 = sin(a*yy)**2
        t1 = cos(a*xx)*sin(a*yy)
        t2 = cos(a*yy)*sin(a*xx)
        t3 = ((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)
        t4 = ((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)
        t5 = cx2*cy2
        t6 = cos(a*xx)*cos(a*yy)
        t7 = sin(a*xx)*sin(a*yy)
        f(i,1) = (2*D*a**2*t7-2*a*xm*t1*t3+2*a*xx*t1*t3+2*a*ym*t2*t3-2*a*yy*t2*t3+D*a**2*xm**2*t7+D*a**2*xx**2*t7+D*a**2*ym**2*t7+D*a**2*yy**2*t7+D*a*xm*t1-D*a*xx*t1+&
          &D*a*ym*t2-D*a*yy*t2+a*xm*cos(a*xx)*cy2*sin(a*xx)*t3-a*xx*cos(a*xx)*cy2*sin(a*xx)*t3-a*ym*cx2*cos(a*yy)*sin(a*yy)*t3+a*yy*cx2*cos(a*yy)*sin(a*yy)*t3-&
          &a*xm*cos(a*xx)*sin(a*xx)*sy2*t3+a*xx*cos(a*xx)*sin(a*xx)*sy2*t3+a*ym*t2**2*sin(a*yy)*t3-a*yy*t2**2*sin(a*yy)*t3-2*D*a**2*xm*ym*t6+2*D*a**2*xx*ym*t6+&
          &2*D*a**2*xm*yy*t6 - 2*D*a**2*xx*yy*t6 - 2*D*a**2*xm*xx*t7 - 2*D*a**2*ym*yy*t7)/((xm - xx)**2 + (ym - yy)**2 + 1)

        f(i,2) = (a*(60*xm*t2*t3-60*xx*t2*t3-60*ym*t1*t3+60*yy*t1*t3+4*xm*t7*t3-4*xx*t7*t3+4*ym*t7*t3-4*yy*t7*t3+12*a*mu*t6-6*mu*xm*t2+&
          &6*mu*xx*t2-6*mu*ym*t1+6*mu*yy*t1-2*xm*cy2*sx2*t3+2*xx*cy2*sx2*t3-2*ym*cx2*sy2*t3+2*yy*cx2*sy2*t3+2*xm*sx2*sy2*t3-2*xx*sx2*sy2*t3+&
          &2*ym*sx2*sy2*t3-2*yy*sx2*sy2*t3+4*xm*t6*t3-4*xx*t6*t3+4*ym*t6*t3-4*yy*t6*t3+6*a*mu*xm**2*t6+6*a*mu*xx**2*t6+6*a*mu*ym**2*t6+&
          &6*a*mu*yy**2*t6 + 2*xm*cx2*cos(a*yy)**3*sin(a*xx)*t3 - 2*xx*cx2*cos(a*yy)**3*sin(a*xx)*t3 - 2*ym*cos(a*xx)**3*cy2*sin(a*yy)*t3 +&
          &2*yy*cos(a*xx)**3*cy2*sin(a*yy)*t3+3*mu*ym*cos(a*xx)*cy2*sin(a*xx)+3*mu*xm*cx2*cos(a*yy)*sin(a*yy)-3*mu*xx*cx2*cos(a*yy)*sin(a*yy)-&
          &3*mu*yy*cos(a*xx)*cy2*sin(a*xx)-3*mu*ym*cos(a*xx)*sin(a*xx)*sy2-3*mu*xm*t2**2*sin(a*yy)+3*mu*xx*t2**2*sin(a*yy)+3*mu*yy*cos(a*xx)*sin(a*xx)*sy2+&
          &8*ym*cos(a*xx)*cy2*sin(a*xx)*t3-8*xm*cx2*cos(a*yy)*sin(a*yy)*t3+8*xx*cx2*cos(a*yy)*sin(a*yy)*t3-8*yy*cos(a*xx)*cy2*sin(a*xx)*t3-&
          &4*xm*cx2*t2*sy2*t3 + 4*xx*cx2*t2*sy2*t3 + 4*ym*cos(a*xx)*cy2*sx2*sin(a*yy)*t3 - 4*yy*cos(a*xx)*cy2*sx2*sin(a*yy)*t3 -&
          &6*a*mu*xm*ym*t5 + 6*a*mu*xx*ym*t5 + 6*a*mu*xm*yy*t5 - 6*a*mu*xx*yy*t5 + 6*a*mu*xm*ym*cx2*sy2 + 6*a*mu*xm*ym*cy2*sx2 -&
          &6*a*mu*xx*ym*cx2*sy2 - 6*a*mu*xx*ym*cy2*sx2 - 6*a*mu*xm*yy*cx2*sy2 - 6*a*mu*xm*yy*cy2*sx2 + 6*a*mu*xx*yy*cx2*sy2 +&
          &6*a*mu*xx*yy*cy2*sx2 - 6*a*mu*xm*ym*sx2*sy2 + 6*a*mu*xx*ym*sx2*sy2 + 6*a*mu*xm*yy*sx2*sy2 - 6*a*mu*xx*yy*sx2*sy2 +&
          &4*xm*cos(a*xx)*t2*sin(a*yy)*t3-4*xx*cos(a*xx)*t2*sin(a*yy)*t3+4*ym*cos(a*xx)*t2*sin(a*yy)*t3-4*yy*cos(a*xx)*t2*sin(a*yy)*t3-&
          &12*a*mu*xm*xx*t6 - 12*a*mu*ym*yy*t6 - 12*a*mu*xm*ym*t7 + 12*a*mu*xx*ym*t7 + 12*a*mu*xm*yy*t7 - 12*a*mu*xx*yy*t7 +&
          &24*a*mu*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*xm**2*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*xx**2*cos(a*xx)*t2*sin(a*yy) +&
          &12*a*mu*ym**2*cos(a*xx)*t2*sin(a*yy) + 12*a*mu*yy**2*cos(a*xx)*t2*sin(a*yy) - 24*a*mu*xm*xx*cos(a*xx)*t2*sin(a*yy) -&
          &24*a*mu*ym*yy*cos(a*xx)*t2*sin(a*yy)))/(3*(xm**2 - 2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1))

        f(i, 3) = (csii*(a*xx - a*xm + 3*a**2*sin(2*a*xx) + 2*a*xm*cx2 - 2*a*xx*cx2 + a*xm*cy2 - a*xx*cy2 + 4*a**2*t1 +&
          &40*a**2*t7+2*a**2*xm**2*sin(2*a*xx)+2*a**2*xx**2*sin(2*a*xx)+a**2*ym**2*sin(2*a*xx)+a**2*yy**2*sin(2*a*xx)-2*a*xm*t5+2*a*xx*t5+&
          &2*a**2*xm**2*t1 + 2*a**2*xx**2*t1 + 2*a**2*ym**2*t1 + 2*a**2*yy**2*t1 + 20*a**2*xm**2*t7 +&
          &20*a**2*xx**2*t7 + 20*a**2*ym**2*t7 + 20*a**2*yy**2*t7 - 8*a**2*cos(a*xx)*cy2*sin(a*xx) +&
          &2*a*ym*t6 - 2*a*yy*t6 + 20*a*xm*t1 - 20*a*xx*t1 + 20*a*ym*t2 - 20*a*yy*t2 - 2*a*xm*t7 + 2*a*xx*t7 -&
          &4*a**2*xm*xx*sin(2*a*xx) + 2*a**2*xm*ym*sin(2*a*yy) - 2*a**2*xx*ym*sin(2*a*yy) -&
          &2*a**2*ym*yy*sin(2*a*xx) - 2*a**2*xm*yy*sin(2*a*yy) + 2*a**2*xx*yy*sin(2*a*yy) - 4*a**2*xm**2*cos(a*xx)*cy2*sin(a*xx) -&
          &4*a**2*xx**2*cos(a*xx)*cy2*sin(a*xx) - 4*a**2*ym**2*cos(a*xx)*cy2*sin(a*xx) -&
          &4*a**2*yy**2*cos(a*xx)*cy2*sin(a*xx) - 40*a**2*xm*ym*t6 + 40*a**2*xx*ym*t6 +&
          &40*a**2*xm*yy*t6 - 40*a**2*xx*yy*t6 - 4*a**2*xm*xx*t1 + 4*a**2*xm*ym*t2 -&
          &4*a**2*xx*ym*t2 - 4*a**2*xm*yy*t2 + 4*a**2*xx*yy*t2 - 4*a**2*ym*yy*t1 -&
          &40*a**2*xm*xx*t7 - 40*a**2*ym*yy*t7 + 8*a**2*xm*xx*cos(a*xx)*cy2*sin(a*xx) -&
          &8*a**2*xm*ym*cx2*cos(a*yy)*sin(a*yy) + 8*a**2*xx*ym*cx2*cos(a*yy)*sin(a*yy) +&
          &8*a**2*ym*yy*cos(a*xx)*cy2*sin(a*xx) + 8*a**2*xm*yy*cx2*cos(a*yy)*sin(a*yy) -&
          &8*a**2*xx*yy*cx2*cos(a*yy)*sin(a*yy) + 2*a*ym*cos(a*xx)*t2*sin(a*yy) -&
          &2*a*yy*cos(a*xx)*t2*sin(a*yy)))/(xm**2 - 2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1) -&
          &(kpari*(2*3**(1 - epn)*Mref*a**2*xx**3*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a**2*xx**3*t1*((2*t1 - t5 + 40)/Mref)**epn + 4*3**(1 - epn)*&
          &a**2*epn*xx**3*t5*((2*t1 - t5 + 40)/Mref)**(epn - 1) + 2*3**(1 - epn)*Mref*a**2*&
          &xx**3*t5*((2*t1 - t5 + 40)/Mref)**epn + 2*3**(1 - epn)*Mref*a*xx**2*t7*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a*ym**2*t7*((2*t1 - t5 + 40)/Mref)**epn - 2*3**(1 - epn)*Mref*a*yy**2*t7*&
          &((2*t1 - t5 + 40)/Mref)**epn + 4*3**(1 - epn)*a**2*epn*xx*ym**2*sx2*sy2*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 4*3**(1 - epn)*a**2*epn*xx*yy**2*sx2*sy2*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) - 2*3**(1 - epn)*Mref*a*xm*xx*t7*&
          &((2*t1 - t5 + 40)/Mref)**epn + 4*3**(1 - epn)*Mref*a*ym*yy*t7*((2*t1 - t5 + 40)/Mref)**epn +&
          &8*3**(1 - epn)*a**2*epn*xx**3*cos(a*xx)**3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**(epn - 1) -&
          &4*3**(1 - epn)*Mref*a**2*xm*xx**2*t5*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xm**2*xx*t5*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xx*ym**2*t5*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xx*yy**2*t5*((2*t1 - t5 + 40)/Mref)**epn +&
          &4*3**(1 - epn)*Mref*a**2*xm*xx**2*t1*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a**2*xm**2*xx*t1*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a**2*xx*ym**2*t1*((2*t1 - t5 + 40)/Mref)**epn -&
          &4*3**(1 - epn)*Mref*a**2*xx**2*ym*t2*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a**2*xx*yy**2*t1*((2*t1 - t5 + 40)/Mref)**epn +&
          &4*3**(1 - epn)*Mref*a**2*xx**2*yy*t2*((2*t1 - t5 + 40)/Mref)**epn -&
          &4*3**(1 - epn)*Mref*a**2*xm*xx**2*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xm**2*xx*cx2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xx*ym**2*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a**2*xx*yy**2*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn -&
          &2*3**(1 - epn)*Mref*a*xx**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a*ym**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a*yy**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**epn -&
          &8*3**(1 - epn)*a**2*epn*xm*xx**2*t5*((2*t1 - t5 + 40)/Mref)**(epn - 1) + 4*3**(1 - epn)*&
          &a**2*epn*xm**2*xx*t5*((2*t1 - t5 + 40)/Mref)**(epn - 1) - 2*3**(1 - epn)*Mref*a*xm*ym*&
          &t6*((2*t1 - t5 + 40)/Mref)**epn + 4*3**(1 - epn)*Mref*a*xx*ym*t6*((2*t1 - t5 + 40)/Mref)**epn +&
          &2*3**(1 - epn)*Mref*a*xm*yy*t6*((2*t1 - t5 + 40)/Mref)**epn - 4*3**(1 - epn)*Mref*a*xx*yy*t6*&
          &((2*t1 - t5 + 40)/Mref)**epn - 8*3**(1 - epn)*a**2*epn*xx*ym*yy*sx2*sy2*((2*t1 - t5 + 40)/Mref)**(epn - 1) -&
          &4*3**(1 - epn)*Mref*a**2*xx*ym*yy*t5*((2*t1 - t5 + 40)/Mref)**epn + 8*3**(1 - epn)*a**2*epn*xx**2*ym*&
          &cx2*cos(a*yy)**3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**(epn - 1) - 16*3**(1 - epn)*a**2*epn*xm*xx**2*&
          &cos(a*xx)**3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*a**2*epn*xm**2*xx*&
          &cos(a*xx)**3*cy2*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**(epn - 1) - 8*3**(1 - epn)*a**2*epn*xx**2*yy*&
          &cx2*cos(a*yy)**3*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**(epn - 1) + 4*3**(1 - epn)*Mref*a**2*xm*xx*ym*&
          &t2*((2*t1 - t5 + 40)/Mref)**epn - 4*3**(1 - epn)*Mref*a**2*xm*xx*yy*t2*((2*t1 - t5 + 40)/Mref)**epn +&
          &4*3**(1 - epn)*Mref*a**2*xx*ym*yy*t1*((2*t1 - t5 + 40)/Mref)**epn - 4*3**(1 - epn)*Mref*a**2*xx*ym*&
          &yy*cy2*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn + 2*3**(1 - epn)*Mref*a*xm*xx*cos(a*xx)*cy2*sin(a*xx)*&
          &((2*t1 - t5 + 40)/Mref)**epn - 2*3**(1 - epn)*Mref*a*xm*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**epn +&
          &4*3**(1 - epn)*Mref*a*xx*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**epn - 4*3**(1 - epn)*Mref*a*&
          &ym*yy*cos(a*xx)*cy2*sin(a*xx)*((2*t1 - t5 + 40)/Mref)**epn + 2*3**(1 - epn)*Mref*a*xm*yy*cx2*&
          &cos(a*yy)*sin(a*yy)*((2*t1 - t5 + 40)/Mref)**epn - 4*3**(1 - epn)*Mref*a*xx*yy*cx2*cos(a*yy)*&
          &sin(a*yy)*((2*t1 - t5 + 40)/Mref)**epn - 8*3**(1 - epn)*Mref*a**2*xx**2*ym*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**epn + 8*3**(1 - epn)*Mref*a**2*xx**2*yy*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**epn + 8*3**(1 - epn)*a**2*epn*xx**2*ym*cos(a*xx)**3*cos(a*yy)**3*t7*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) - 8*3**(1 - epn)*a**2*epn*xx**2*yy*cos(a*xx)**3*cos(a*yy)**3*&
          &t7*((2*t1 - t5 + 40)/Mref)**(epn - 1) - (4*3**(1 - epn)*Mref*a**2*epn*xx**3*cos(a*xx)**4*cy2*(cy2 - 1)*&
          &((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) - 8*3**(1 - epn)*a**2*epn*xx**2*ym*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*a**2*epn*xx**2*yy*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) - 8*3**(1 - epn)*a**2*epn*xm*xx*ym*cx2*cos(a*yy)**3*sin(a*xx)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*a**2*epn*xm*xx*yy*cx2*cos(a*yy)**3*sin(a*xx)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*Mref*a**2*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**epn - 8*3**(1 - epn)*Mref*a**2*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**epn - 8*3**(1 - epn)*a**2*epn*xm*xx*ym*cos(a*xx)**3*cos(a*yy)**3*t7*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*a**2*epn*xm*xx*yy*cos(a*xx)**3*cos(a*yy)**3*t7*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) + 8*3**(1 - epn)*a**2*epn*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) - 8*3**(1 - epn)*a**2*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*&
          &((2*t1 - t5 + 40)/Mref)**(epn - 1) - (4*3**(1 - epn)*Mref*a**2*epn*xx*ym**2*cx2*cos(a*yy)**4*&
          &(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) + (8*3**(1 - epn)*Mref*a**2*epn*xm*xx**2*&
          &cos(a*xx)**4*cy2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) - (4*3**(1 - epn)*Mref*&
          &a**2*epn*xm**2*xx*cos(a*xx)**4*cy2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) -&
          &(4*3**(1 - epn)*Mref*a**2*epn*xx*yy**2*cx2*cos(a*yy)**4*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) +&
          &(8*3**(1 - epn)*Mref*a**2*epn*xx*ym*yy*cx2*cos(a*yy)**4*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) +&
          &(8*3**(1 - epn)*Mref*a**2*epn*xx*ym**2*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) +&
          &(8*3**(1 - epn)*Mref*a**2*epn*xx**2*ym*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) + (8*3**(1 - epn)*Mref*&
          &a**2*epn*xx*yy**2*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) - (8*3**(1 - epn)*Mref*&
          &a**2*epn*xx**2*yy*cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) - (8*3**(1 - epn)*Mref*a**2*epn*xm*xx*ym*&
          &cx2*t2*(cy2 - 1)*((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) + (8*3**(1 - epn)*Mref*a**2*epn*xm*xx*yy*cx2*t2*(cy2 - 1)*&
          &((2*t1 - t5 + 40)/Mref)**epn)/(2*t1 - t5 + 40) - (16*3**(1 - epn)*Mref*a**2*epn*xx*ym*yy*cos(a*xx)*cy2*sin(a*yy)*(cx2 - 1)*&
          &((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)))/(9*Mref**2*xx*((xm-xx)**2+(ym-yy)**2+1))+Mref*pcourr*t6*((2*a*sin(a*xx)*(xm-xx)*&
          &(10*cos(a*yy) + sin(a*xx) + 2*sin(a*yy) - 2*cy2*sin(a*xx)))/(3*Mref*t3) + (4*a*cos(a*xx)*(ym - yy)*(cos(a*yy) - 5*sin(a*yy) +&
          &t2*sin(a*yy)))/(3*Mref*t3)) - (a*t6*(ym - yy)*(5*cx2*sy2 - 5*sx2*sy2 + 100*t1 - 10*t7 - cos(a*xx)**3*cy2*sin(a*yy) + 4*cos(a*xx)*&
          &cy2*sin(a*xx) + 2*cos(a*xx)*cy2*sx2*sin(a*yy)))/(3*((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)) + (a*cos(a*xx)*cy2*(xm - xx)*&
          &(10*cos(a*xx) + 100*sin(a*xx) + 2*cx2*sin(a*xx) + 4*cx2*sin(a*yy) - 3*t5*sin(a*xx) + 10*cos(a*xx)*t7))/(3*t3) +&
          &(3**(0.5)*(t7 + 2)**2*(2*t1 - t5 + 2*t2 + 20))/(2*tie*(t2 - 10)*(-(2*t2 - 20)/Mref)**(0.5)) -&
          &(a*t1*(t7 + 2)*(xm - xx)*(5*t1 - t5 + 100))/(3*((xm - xx)**2 + (ym - yy)**2 + 1)**(0.5)) +&
          &(a*t2*(t7+2)*(ym-yy)*(5*t1-t5+100))/(3*t3)-(t6*(2*xm-2*xx)*(t7+2)*(ym-yy)*(5*t1-t5+100))/(6*t4)+(t6*(2*ym-2*yy)*(t7+2)*(xm-xx)*(5*t1-t5+100))/(6*t4)

        f(i, 4) = (csie*(a*ym - a*yy - 3*a**2*sin(2*a*yy) - a*ym*cx2 + a*yy*cx2 - 2*a*ym*cy2 + 2*a*yy*cy2 -&
          &4*a**2*t2 + 20*a**2*t7 - a**2*xm**2*sin(2*a*yy) - a**2*xx**2*sin(2*a*yy) - 2*a**2*ym**2*&
          &sin(2*a*yy) - 2*a**2*yy**2*sin(2*a*yy) + 2*a*ym*t5 - 2*a*yy*t5 - 2*a**2*xm**2*t2 -&
          &2*a**2*xx**2*t2 - 2*a**2*ym**2*t2 - 2*a**2*yy**2*t2 + 10*a**2*xm**2*t7 + 10*a**2*xx**2*t7 +&
          &10*a**2*ym**2*t7 + 10*a**2*yy**2*t7 + 8*a**2*cx2*cos(a*yy)*sin(a*yy) - 2*a*xm*t6 + 2*a*&
          &xx*t6 + 10*a*xm*t1 - 10*a*xx*t1 + 10*a*ym*t2 - 10*a*yy*t2 + 2*a*ym*t7 - 2*a*yy*t7 - 2*a**2*&
          &xm*ym*sin(2*a*xx) + 2*a**2*xx*ym*sin(2*a*xx) + 2*a**2*xm*xx*sin(2*a*yy) + 2*a**2*xm*&
          &yy*sin(2*a*xx) - 2*a**2*xx*yy*sin(2*a*xx) + 4*a**2*ym*yy*sin(2*a*yy) + 4*a**2*xm**2*&
          &cx2*cos(a*yy)*sin(a*yy) + 4*a**2*xx**2*cx2*cos(a*yy)*sin(a*yy) + 4*a**2*ym**2*cx2*cos(a*yy)*&
          &sin(a*yy) + 4*a**2*yy**2*cx2*cos(a*yy)*sin(a*yy) - 20*a**2*xm*ym*t6 + 20*a**2*xx*ym*t6 + 20*a**2*xm*yy*t6 -&
          &20*a**2*xx*yy*t6 + 4*a**2*xm*xx*t2 - 4*a**2*xm*ym*t1 + 4*a**2*xx*ym*t1 + 4*a**2*xm*yy*t1 - 4*&
          &a**2*xx*yy*t1 + 4*a**2*ym*yy*t2 - 20*a**2*xm*xx*t7 - 20*a**2*ym*yy*t7 + 8*a**2*xm*ym*cos(a*xx)*&
          &cy2*sin(a*xx) - 8*a**2*xx*ym*cos(a*xx)*cy2*sin(a*xx) - 8*a**2*xm*xx*cx2*cos(a*yy)*sin(a*yy) -&
          &8*a**2*xm*yy*cos(a*xx)*cy2*sin(a*xx) + 8*a**2*xx*yy*cos(a*xx)*cy2*sin(a*xx) - 8*a**2*ym*yy*&
          &cx2*cos(a*yy)*sin(a*yy) - 2*a*xm*cos(a*xx)*t2*sin(a*yy) + 2*a*xx*cos(a*xx)*t2*sin(a*yy)))/(xm**2 -&
          &2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1) - Mref*pcourr*t6*((2*a*sin(a*xx)*(xm - xx)*(10*cos(a*yy) +&
          &sin(a*xx) + 2*sin(a*yy) - 2*cy2*sin(a*xx)))/(3*Mref*t3) + (4*a*cos(a*xx)*(ym - yy)*(cos(a*yy) -&
          &5*sin(a*yy) + t2*sin(a*yy)))/(3*Mref*t3)) + (10*a*cx2*cos(a*yy)*(ym - yy)*(cos(a*yy) - 5*sin(a*yy) +&
          &t2*sin(a*yy)))/(3*t3) - (3**(0.5)*(t7 + 2)**2*(2*t1 - t5 + 2*t2 + 20))/(2*tie*(t2 - 10)*(-(2*t2 -&
          &20)/Mref)**(0.5)) + (5*a*cos(a*xx)*t2*(xm - xx)*(10*cos(a*yy) + sin(a*xx) + 2*sin(a*yy) - 2*cy2*&
          &sin(a*xx)))/(3*t3) + (5*t6*(t2 - 10)*(2*xm - 2*xx)*(t7 + 2)*(ym - yy))/(6*t4) - (5*t6*(t2 - 10)*&
          &(2*ym - 2*yy)*(t7 + 2)*(xm - xx))/(6*t4) + (5*a*t1*(t2 - 10)*(t7 + 2)*(xm - xx))/(3*t3) - (5*a*t2*(t2 - 10)*&
          &(t7 + 2)*(ym - yy))/(3*t3) + (2/3**(epn + 1)*a*kpare*(-(2*t2 - 20)/Mref)**epn*(10*xx**2*t6 - 10*ym**2*&
          &t6 - 10*yy**2*t6 - xx**2*cos(a*xx)*cy2*sin(a*xx) + ym**2*cos(a*xx)*cy2*sin(a*xx) + yy**2*cos(a*xx)*&
          &cy2*sin(a*xx) - 10*xm*xx*t6 + 20*ym*yy*t6 - a*xx**3*cy2*sx2 - 10*xm*ym*t7 + 20*xx*ym*t7 + 10*xm*yy*t7 -&
          &20*xx*yy*t7 + 10*a*xx**3*t2 + a*epn*xx**3*sx2*sy2 + 2*a*xm*xx**2*cy2*sx2 - a*xm**2*xx*cy2*sx2 - a*xx*&
          &ym**2*cy2*sx2 - a*xx*yy**2*cy2*sx2 - 20*a*xm*xx**2*t2 + 10*a*xm**2*xx*t2 + 10*a*xx*ym**2*t2 + 20*a*&
          &xx**2*ym*t1 + 10*a*xx*yy**2*t2 - 20*a*xx**2*yy*t1 + xm*xx*cos(a*xx)*cy2*sin(a*xx) - 2*ym*yy*&
          &cos(a*xx)*cy2*sin(a*xx) + xm*ym*t2**2*sin(a*yy) - 2*xx*ym*t2**2*sin(a*yy) - xm*yy*t2**2*&
          &sin(a*yy) + 2*xx*yy*t2**2*sin(a*yy) + 2*a*xx*ym*yy*cy2*sx2 + a*epn*xx*ym**2*t5 + a*epn*xx*yy**2*t5 - 20*a*xm*xx*ym*t1 + 20*a*xm*&
          &xx*yy*t1 - 20*a*xx*ym*yy*t2 - 2*a*epn*xm*xx**2*sx2*sy2 + a*epn*xm**2*xx*sx2*sy2 - 2*a*epn*xx*&
          &ym*yy*t5 - 2*a*xx**2*ym*cos(a*xx)*t2*sin(a*yy) + 2*a*xx**2*yy*cos(a*xx)*t2*sin(a*yy) -&
          &2*a*epn*xx**2*ym*cos(a*xx)*t2*sin(a*yy) + 2*a*epn*xx**2*yy*cos(a*xx)*t2*sin(a*yy) +&
          &2*a*xm*xx*ym*cos(a*xx)*t2*sin(a*yy) - 2*a*xm*xx*yy*cos(a*xx)*t2*sin(a*yy) + 2*a*epn*&
          &xm*xx*ym*cos(a*xx)*t2*sin(a*yy) - 2*a*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)))/(Mref*&
          &xx*(t2 - 10)*(xm**2 - 2*ym*yy - 2*xm*xx + xx**2 + ym**2 + yy**2 + 1))

      END DO

      !      f(:,1) =  cos(a * x) ** 2 * a * sin(b * y) * cos(b * y) * (x / 0.30D2 &
      !     &- y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - (0.2D1 + sin(a * x) * sin(b &
      !     &* y)) * sin(a * x) * a * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D&
      !     &2 + 0.1D1 / 0.15D2) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * &
      !     &x) * cos(b * y) / 0.30D2 + sin(a * x) * cos(b * y) ** 2 * b * cos(&
      !     &a * x) * (x * y / 0.30D2 + y / 0.30D2) - (0.2D1 + sin(a * x) * sin&
      !     &(b * y)) * cos(a * x) * sin(b * y) * b * (x * y / 0.30D2 + y / 0.3&
      !     &0D2) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y)&
      !     & * (x / 0.30D2 + 0.1D1 / 0.30D2) - D * (-sin(a * x) * a ** 2 * sin&
      !     &(b * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a &
      !     &* x) * a * sin(b * y) / 0.30D2 - (x * y / 0.30D2 + y / 0.30D2) * s&
      !     &in(a * x) * cos(b * y) * b / 0.30D2 - (x / 0.30D2 - y ** 2 / 0.30D&
      !     &2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) / 0.30D2 - (x /&
      !     & 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a ** 2 &
      !     &* sin(b * y) + y * sin(a * x) * cos(b * y) * b / 0.30D2 + (x * y /&
      !     & 0.30D2 + y / 0.30D2) * cos(a * x) * a * cos(b * y) * b) - sin(a *&
      !     & x) * sin(b * y) * b ** 2 - (x / 0.30D2 + 0.1D1 / 0.30D2) * ((x / &
      !     &0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * sin(&
      !     &b * y) + (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * cos(b * y) *&
      !     & b) - (x * y / 0.30D2 + y / 0.30D2) * (-y * cos(a * x) * a * sin(b&
      !     & * y) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) *&
      !     & cos(a * x) * a * cos(b * y) * b + (x / 0.30D2 + 0.1D1 / 0.30D2) *&
      !     & sin(a * x) * cos(b * y) * b - (x * y / 0.30D2 + y / 0.30D2) * sin&
      !     &(a * x) * sin(b * y) * b ** 2))
      !
      !      f(:,2) = cos(a * x) ** 3 * a * sin(b * y) * cos(b * y) ** 2 * (x / 0.&
      !     &30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - 0.2D1 * (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * cos(a * x) * cos(b * y) ** 2 * (x / 0.30D2 - &
      !     &y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a + (0.2D1 + sin(&
      !     &a * x) * sin(b * y)) * cos(a * x) ** 2 * cos(b * y) ** 2 / 0.30D2 &
      !     &+ sin(a * x) * cos(b * y) ** 3 * b * cos(a * x) ** 2 * (x * y / 0.&
      !     &30D2 + y / 0.30D2) - 0.2D1 * (0.2D1 + sin(a * x) * sin(b * y)) * c&
      !     &os(a * x) ** 2 * cos(b * y) * (x * y / 0.30D2 + y / 0.30D2) * sin(&
      !     &b * y) * b + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) ** 2 *&
      !     & cos(b * y) ** 2 * (x / 0.30D2 + 0.1D1 / 0.30D2) + Mref * ((x / 0.&
      !     &30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (0.2D1 / 0.3D1 * cos(a &
      !     &* x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a &
      !     &* x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref + 0.2D1 / 0.3D1 * (0.2&
      !     &D1 + sin(a * x) * sin(b * y)) * (-sin(a * x) * a * sin(b * y) + co&
      !     &s(a * x) * cos(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D&
      !     &1 * cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y&
      !     &)) / Mref - 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * co&
      !     &s(a * x) * a * cos(b * y) / Mref) + (x * y / 0.30D2 + y / 0.30D2) &
      !     &* (0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a *&
      !     & x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mr&
      !     &ef + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (cos(a * &
      !     &x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b * y) * &
      !     &b) / Mref + 0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.10D2 &
      !     &- sin(a * x) * cos(b * y)) / Mref + 0.2D1 / 0.3D1 * (0.2D1 + sin(a&
      !     & * x) * sin(b * y)) * sin(a * x) * sin(b * y) * b / Mref)) - mu * &
      !     &(-0.3D1 * cos(a * x) * a ** 2 * sin(b * y) * cos(b * y) * sin(a * &
      !     &x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a ** 2 * cos&
      !     &(b * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a&
      !     & * x) ** 2 * a * sin(b * y) * cos(b * y) - (0.2D1 + sin(a * x) * s&
      !     &in(b * y)) * sin(a * x) * a * cos(b * y)) / 0.30D2 - (x * y / 0.30&
      !     &D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) &
      !     &- (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b)&
      !     & / 0.30D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos&
      !     &(a * x) ** 2 * a * sin(b * y) * cos(b * y) / 0.30D2 - (0.2D1 + sin&
      !     &(a * x) * sin(b * y)) * sin(a * x) * a * cos(b * y) / 0.30D2 + (x &
      !     &/ 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-0.3D1 * cos(a * x&
      !     &) * a ** 2 * sin(b * y) * cos(b * y) * sin(a * x) - (0.2D1 + sin(a&
      !     & * x) * sin(b * y)) * cos(a * x) * a ** 2 * cos(b * y)) + y * (sin&
      !     &(a * x) * cos(b * y) ** 2 * b * cos(a * x) - (0.2D1 + sin(a * x) *&
      !     & sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.30D2 + (x * y / 0.&
      !     &30D2 + y / 0.30D2) * (cos(a * x) ** 2 * a * cos(b * y) ** 2 * b - &
      !     &sin(a * x) ** 2 * cos(b * y) ** 2 * b * a - cos(a * x) ** 2 * a * &
      !     &sin(b * y) ** 2 * b + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * &
      !     &x) * a * sin(b * y) * b)) - 0.3D1 * sin(a * x) * cos(b * y) * b **&
      !     & 2 * cos(a * x) * sin(b * y) - (0.2D1 + sin(a * x) * sin(b * y)) *&
      !     & cos(a * x) * cos(b * y) * b ** 2 - (x / 0.30D2 + 0.1D1 / 0.30D2) &
      !     &* ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) *&
      !     &* 2 * a * sin(b * y) * cos(b * y) - (0.2D1 + sin(a * x) * sin(b * &
      !     &y)) * sin(a * x) * a * cos(b * y)) + (x * y / 0.30D2 + y / 0.30D2)&
      !     & * (sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) - (0.2D1 + sin(a&
      !     & * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b)) - (x * y / 0.3&
      !     &0D2 + y / 0.30D2) * (-y * (cos(a * x) ** 2 * a * sin(b * y) * cos(&
      !     &b * y) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * cos(&
      !     &b * y)) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2)&
      !     & * (cos(a * x) ** 2 * a * cos(b * y) ** 2 * b - sin(a * x) ** 2 * &
      !     &cos(b * y) ** 2 * b * a - cos(a * x) ** 2 * a * sin(b * y) ** 2 * &
      !     &b + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b * y&
      !     &) * b) + (x / 0.30D2 + 0.1D1 / 0.30D2) * (sin(a * x) * cos(b * y) &
      !     &** 2 * b * cos(a * x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a &
      !     &* x) * sin(b * y) * b) + (x * y / 0.30D2 + y / 0.30D2) * (-0.3D1 *&
      !     & sin(a * x) * cos(b * y) * b ** 2 * cos(a * x) * sin(b * y) - (0.2&
      !     &D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b ** 2))&
      !     &)

      !      f(:,3) = (cos(a * x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b &
      !     &* y)) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b&
      !     & * y) + 0.2D1 / 0.3D1 * cos(a * x) * a * sin(b * y) * (0.20D2 + co&
      !     &s(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1)&
      !     & + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (-sin(a * x&
      !     &) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) * a&
      !     &)) * cos(a * x) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1&
      !     &D1 / 0.15D2) - ((0.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(&
      !     &a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b&
      !     & * y)) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos&
      !     &(b * y) ** 2 / 0.2D1)) * sin(a * x) * a * cos(b * y) * (x / 0.30D2&
      !     & - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) + ((0.2D1 + sin(a * x) * sin(&
      !     &b * y)) * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.&
      !     &2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(a * x) * sin(b * y)&
      !     & - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1)) * cos(a * x) * cos(&
      !     &b * y) / 0.30D2 + (sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a *&
      !     & x) * sin(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x)&
      !     & * cos(b * y) * b + 0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * &
      !     &(0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) *&
      !     &* 2 / 0.2D1) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) *&
      !     & (cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin&
      !     &(b * y) * b)) * cos(a * x) * cos(b * y) * (x * y / 0.30D2 + y / 0.&
      !     &30D2) - ((0.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(a * x) &
      !     &* sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) &
      !     &* (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y)&
      !     & ** 2 / 0.2D1)) * cos(a * x) * sin(b * y) * b * (x * y / 0.30D2 + &
      !     &y / 0.30D2) + ((0.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(a&
      !     & * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b &
      !     &* y)) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(&
      !     &b * y) ** 2 / 0.2D1)) * cos(a * x) * cos(b * y) * (x / 0.30D2 + 0.&
      !     &1D1 / 0.30D2) - csii * (-sin(a * x) * a ** 2 * sin(b * y) * (0.20D&
      !     &2 + cos(a * x) * sin(b * y)) - 0.2D1 * cos(a * x) * a ** 2 * sin(b&
      !     & * y) ** 2 * sin(a * x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(&
      !     &a * x) * a ** 2 * sin(b * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1&
      !     &D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) * (0.20D2 + cos(a * x)&
      !     & * sin(b * y)) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * &
      !     &a * sin(b * y)) / 0.30D2 - (x * y / 0.30D2 + y / 0.30D2) * (sin(a &
      !     &* x) * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)) + (0.2D&
      !     &1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b) / 0.30&
      !     &D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x)&
      !     & * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)) / 0.30D2 - &
      !     &(0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b * y) / &
      !     &0.30D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a&
      !     & * x) * a ** 2 * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)) -&
      !     & 0.2D1 * cos(a * x) * a ** 2 * sin(b * y) ** 2 * sin(a * x) - (0.2&
      !     &D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a ** 2 * sin(b * y)) &
      !     &+ y * (sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b &
      !     &* y)) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y&
      !     &) * b) / 0.30D2 + (x * y / 0.30D2 + y / 0.30D2) * (cos(a * x) * a &
      !     &* cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)) - sin(a * x)&
      !     & ** 2 * cos(b * y) * b * a * sin(b * y) + cos(a * x) ** 2 * a * si&
      !     &n(b * y) * cos(b * y) * b - (0.2D1 + sin(a * x) * sin(b * y)) * si&
      !     &n(a * x) * a * cos(b * y) * b)) - sin(a * x) * sin(b * y) * b ** 2&
      !     & * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 * sin(a * x) * cos(b&
      !     & * y) ** 2 * b ** 2 * cos(a * x) - (0.2D1 + sin(a * x) * sin(b * y&
      !     &)) * cos(a * x) * sin(b * y) * b ** 2 - (x / 0.30D2 + 0.1D1 / 0.30&
      !     &D2) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * &
      !     &x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)) - (0.2D1 &
      !     &+ sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b * y)) + (x * y&
      !     & / 0.30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.20D2 +&
      !     & cos(a * x) * sin(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * co&
      !     &s(a * x) * cos(b * y) * b)) - (x * y / 0.30D2 + y / 0.30D2) * (-y &
      !     &* (cos(a * x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)&
      !     &) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b * y&
      !     &)) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (c&
      !     &os(a * x) * a * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)&
      !     &) - sin(a * x) ** 2 * cos(b * y) * b * a * sin(b * y) + cos(a * x)&
      !     & ** 2 * a * sin(b * y) * cos(b * y) * b - (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * sin(a * x) * a * cos(b * y) * b) + (x / 0.30D2 + 0.1D1&
      !     & / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a * x) *&
      !     & sin(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * co&
      !     &s(b * y) * b) + (x * y / 0.30D2 + y / 0.30D2) * (-sin(a * x) * sin&
      !     &(b * y) * b ** 2 * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 * si&
      !     &n(a * x) * cos(b * y) ** 2 * b ** 2 * cos(a * x) - (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * cos(a * x) * sin(b * y) * b ** 2))) - kpari *&
      !     & ((0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) &
      !     &** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * epn * (-sin(a * x&
      !     &) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) * a&
      !     &) / (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * &
      !     &y) ** 2 / 0.2D1) * (0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 &
      !     &+ 0.1D1 / 0.15D2) * (-sin(a * x) * a * sin(b * y) + cos(a * x) * c&
      !     &os(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1 * (x * y /&
      !     & 0.30D2 + y / 0.30D2) * (cos(a * x) * cos(b * y) * b + cos(a * x) &
      !     &** 2 * cos(b * y) * sin(b * y) * b) / Mref) * (x / 0.30D2 - y ** 2&
      !     & / 0.30D2 + 0.1D1 / 0.15D2) + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x&
      !     &) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref&
      !     &) ** epn * ((-sin(a * x) * a * sin(b * y) + cos(a * x) * cos(b * y&
      !     &) ** 2 * sin(a * x) * a) / Mref / 0.45D2 + 0.2D1 / 0.3D1 * (x / 0.&
      !     &30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-cos(a * x) * a ** 2 *&
      !     & sin(b * y) - sin(a * x) ** 2 * a ** 2 * cos(b * y) ** 2 + cos(a *&
      !     & x) ** 2 * cos(b * y) ** 2 * a ** 2) / Mref + y * (cos(a * x) * co&
      !     &s(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) / Mr&
      !     &ef / 0.45D2 + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * (-si&
      !     &n(a * x) * a * cos(b * y) * b - 0.2D1 * cos(a * x) * cos(b * y) * &
      !     &sin(b * y) * b * sin(a * x) * a) / Mref) * (x / 0.30D2 - y ** 2 / &
      !     &0.30D2 + 0.1D1 / 0.15D2) + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) *&
      !     & sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) *&
      !     &* epn * (0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0&
      !     &.15D2) * (-sin(a * x) * a * sin(b * y) + cos(a * x) * cos(b * y) *&
      !     &* 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y&
      !     & / 0.30D2) * (cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(&
      !     &b * y) * sin(b * y) * b) / Mref) / 0.30D2 + (0.2D1 / 0.3D1 * (0.20&
      !     &D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 /&
      !     & 0.2D1) / Mref) ** epn * epn * (cos(a * x) * cos(b * y) * b + cos(&
      !     &a * x) ** 2 * cos(b * y) * sin(b * y) * b) / (0.20D2 + cos(a * x) &
      !     &* sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) * (0.2D1&
      !     & / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin&
      !     &(a * x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * &
      !     &x) * a) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * (&
      !     &cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b&
      !     & * y) * b) / Mref) * (x * y / 0.30D2 + y / 0.30D2) + (0.2D1 / 0.3D&
      !     &1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * &
      !     &y) ** 2 / 0.2D1) / Mref) ** epn * (-0.2D1 / 0.45D2 * y * (-sin(a *&
      !     & x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) *&
      !     & a) / Mref + 0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1&
      !     & / 0.15D2) * (-sin(a * x) * a * cos(b * y) * b - 0.2D1 * cos(a * x&
      !     &) * cos(b * y) * sin(b * y) * b * sin(a * x) * a) / Mref + 0.2D1 /&
      !     & 0.3D1 * (x / 0.30D2 + 0.1D1 / 0.30D2) * (cos(a * x) * cos(b * y) &
      !     &* b + cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) / Mref + 0.2D&
      !     &1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * (-cos(a * x) * sin(b *&
      !     & y) * b ** 2 - cos(a * x) ** 2 * sin(b * y) ** 2 * b ** 2 + cos(a &
      !     &* x) ** 2 * cos(b * y) ** 2 * b ** 2) / Mref) * (x * y / 0.30D2 + &
      !     &y / 0.30D2) + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) -&
      !     & cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * (0.2D&
      !     &1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-si&
      !     &n(a * x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a *&
      !     & x) * a) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * &
      !     &(cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(&
      !     &b * y) * b) / Mref) * (x / 0.30D2 + 0.1D1 / 0.30D2)) + pcourr * Mr&
      !     &ef * cos(a * x) * cos(b * y) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.&
      !     &1D1 / 0.15D2) * (0.2D1 / 0.3D1 * cos(a * x) * a * sin(b * y) * (0.&
      !     &10D2 - sin(a * x) * cos(b * y)) / Mref - 0.2D1 / 0.3D1 * (0.2D1 + &
      !     &sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y) / Mref) + (&
      !     &x * y / 0.30D2 + y / 0.30D2) * (0.2D1 / 0.3D1 * sin(a * x) * cos(b&
      !     & * y) * b * (0.10D2 - sin(a * x) * cos(b * y)) / Mref + 0.2D1 / 0.&
      !     &3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * sin(b * y) &
      !     &* b / Mref)) + 0.3D1 / 0.4D1 * (0.2D1 + sin(a * x) * sin(b * y)) *&
      !     &* 2 / tie * (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / &
      !     &Mref - 0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a *&
      !     & x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) * sqrt(0.2D1) * sqrt(0&
      !     &.3D1) * ((0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** (-0.3D1 / 0&
      !     &.2D1)

      !      f(:,4) = 0.5D1 / 0.3D1 * cos(a * x) ** 2 * a * sin(b * y) * (0.10D2 -&
      !     & sin(a * x) * cos(b * y)) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.&
      !     &30D2 + 0.1D1 / 0.15D2) - 0.5D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin&
      !     &(b * y)) * cos(a * x) ** 2 * a * cos(b * y) ** 2 * (x / 0.30D2 - y&
      !     & ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - 0.5D1 / 0.3D1 * (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * (0.10D2 - sin(a * x) * cos(b * y)) * sin(a * &
      !     &x) * a * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15&
      !     &D2) + (0.2D1 + sin(a * x) * sin(b * y)) * (0.10D2 - sin(a * x) * c&
      !     &os(b * y)) * cos(a * x) * cos(b * y) / 0.18D2 + 0.5D1 / 0.3D1 * si&
      !     &n(a * x) * cos(b * y) ** 2 * b * (0.10D2 - sin(a * x) * cos(b * y)&
      !     &) * cos(a * x) * (x * y / 0.30D2 + y / 0.30D2) + 0.5D1 / 0.3D1 * (&
      !     &0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * sin(b * y) * b * c&
      !     &os(a * x) * cos(b * y) * (x * y / 0.30D2 + y / 0.30D2) - 0.5D1 / 0&
      !     &.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (0.10D2 - sin(a * x) * &
      !     &cos(b * y)) * cos(a * x) * sin(b * y) * b * (x * y / 0.30D2 + y / &
      !     &0.30D2) + 0.5D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (0.1&
      !     &0D2 - sin(a * x) * cos(b * y)) * cos(a * x) * cos(b * y) * (x / 0.&
      !     &30D2 + 0.1D1 / 0.30D2) - csie * (-sin(a * x) * a ** 2 * sin(b * y)&
      !     & * (0.10D2 - sin(a * x) * cos(b * y)) - 0.2D1 * cos(a * x) ** 2 * &
      !     &a ** 2 * sin(b * y) * cos(b * y) + (0.2D1 + sin(a * x) * sin(b * y&
      !     &)) * sin(a * x) * a ** 2 * cos(b * y) - (x / 0.30D2 - y ** 2 / 0.3&
      !     &0D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) * (0.10D2 - s&
      !     &in(a * x) * cos(b * y)) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(&
      !     &a * x) * a * cos(b * y)) / 0.30D2 - (x * y / 0.30D2 + y / 0.30D2) &
      !     &* (sin(a * x) * cos(b * y) * b * (0.10D2 - sin(a * x) * cos(b * y)&
      !     &) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * sin(b * y) * &
      !     &b) / 0.30D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (c&
      !     &os(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) / &
      !     &0.30D2 - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(&
      !     &b * y) / 0.30D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) &
      !     &* (-sin(a * x) * a ** 2 * sin(b * y) * (0.10D2 - sin(a * x) * cos(&
      !     &b * y)) - 0.2D1 * cos(a * x) ** 2 * a ** 2 * sin(b * y) * cos(b * &
      !     &y) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a ** 2 * cos&
      !     &(b * y)) + y * (sin(a * x) * cos(b * y) * b * (0.10D2 - sin(a * x)&
      !     & * cos(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * &
      !     &sin(b * y) * b) / 0.30D2 + (x * y / 0.30D2 + y / 0.30D2) * (cos(a &
      !     &* x) * a * cos(b * y) * b * (0.10D2 - sin(a * x) * cos(b * y)) - s&
      !     &in(a * x) * cos(b * y) ** 2 * b * cos(a * x) * a + cos(a * x) * a &
      !     &* sin(b * y) ** 2 * sin(a * x) * b + (0.2D1 + sin(a * x) * sin(b *&
      !     & y)) * cos(a * x) * a * sin(b * y) * b)) - sin(a * x) * sin(b * y)&
      !     & * b ** 2 * (0.10D2 - sin(a * x) * cos(b * y)) + 0.2D1 * sin(a * x&
      !     &) ** 2 * cos(b * y) * b ** 2 * sin(b * y) + (0.2D1 + sin(a * x) * &
      !     &sin(b * y)) * sin(a * x) * cos(b * y) * b ** 2 - (x / 0.30D2 + 0.1&
      !     &D1 / 0.30D2) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * &
      !     &(cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) &
      !     &- (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y))&
      !     & + (x * y / 0.30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) * b * &
      !     &(0.10D2 - sin(a * x) * cos(b * y)) + (0.2D1 + sin(a * x) * sin(b *&
      !     & y)) * sin(a * x) * sin(b * y) * b)) - (x * y / 0.30D2 + y / 0.30D&
      !     &2) * (-y * (cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * c&
      !     &os(b * y)) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * &
      !     &cos(b * y)) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.1&
      !     &5D2) * (cos(a * x) * a * cos(b * y) * b * (0.10D2 - sin(a * x) * c&
      !     &os(b * y)) - sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) * a + c&
      !     &os(a * x) * a * sin(b * y) ** 2 * sin(a * x) * b + (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * cos(a * x) * a * sin(b * y) * b) + (x / 0.30D&
      !     &2 + 0.1D1 / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.10D2 - sin&
      !     &(a * x) * cos(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a &
      !     &* x) * sin(b * y) * b) + (x * y / 0.30D2 + y / 0.30D2) * (-sin(a *&
      !     & x) * sin(b * y) * b ** 2 * (0.10D2 - sin(a * x) * cos(b * y)) + 0&
      !     &.2D1 * sin(a * x) ** 2 * cos(b * y) * b ** 2 * sin(b * y) + (0.2D1&
      !     & + sin(a * x) * sin(b * y)) * sin(a * x) * cos(b * y) * b ** 2))) &
      !     &- kpare * (-(0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / &
      !     &Mref) ** epn * epn * cos(a * x) * a * cos(b * y) / (0.10D2 - sin(a&
      !     & * x) * cos(b * y)) * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.3&
      !     &0D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1&
      !     & / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y)&
      !     & * b / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) + (&
      !     &0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** epn &
      !     &* (-cos(a * x) * a * cos(b * y) / Mref / 0.45D2 + 0.2D1 / 0.3D1 * &
      !     &(x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a *&
      !     &* 2 * cos(b * y) / Mref + y * sin(a * x) * sin(b * y) * b / Mref /&
      !     & 0.45D2 + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * cos(a * &
      !     &x) * a * sin(b * y) * b / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + &
      !     &0.1D1 / 0.15D2) + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * &
      !     &y)) / Mref) ** epn * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30&
      !     &D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 &
      !     &/ 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) &
      !     &* b / Mref) / 0.30D2 + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos&
      !     &(b * y)) / Mref) ** epn * epn * sin(a * x) * sin(b * y) * b / (0.1&
      !     &0D2 - sin(a * x) * cos(b * y)) * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y&
      !     & ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * cos(b * y) / M&
      !     &ref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) *&
      !     & sin(b * y) * b / Mref) * (x * y / 0.30D2 + y / 0.30D2) + (0.2D1 /&
      !     & 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** epn * (0.2D&
      !     &1 / 0.45D2 * y * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D&
      !     &1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) *&
      !     & a * sin(b * y) * b / Mref + 0.2D1 / 0.3D1 * (x / 0.30D2 + 0.1D1 /&
      !     & 0.30D2) * sin(a * x) * sin(b * y) * b / Mref + 0.2D1 / 0.3D1 * (x&
      !     & * y / 0.30D2 + y / 0.30D2) * sin(a * x) * cos(b * y) * b ** 2 / M&
      !     &ref) * (x * y / 0.30D2 + y / 0.30D2) + (0.2D1 / 0.3D1 * (0.10D2 - &
      !     &sin(a * x) * cos(b * y)) / Mref) ** epn * (-0.2D1 / 0.3D1 * (x / 0&
      !     &.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * cos(b&
      !     & * y) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * sin&
      !     &(a * x) * sin(b * y) * b / Mref) * (x / 0.30D2 + 0.1D1 / 0.30D2)) &
      !     &- pcourr * Mref * cos(a * x) * cos(b * y) * ((x / 0.30D2 - y ** 2 &
      !     &/ 0.30D2 + 0.1D1 / 0.15D2) * (0.2D1 / 0.3D1 * cos(a * x) * a * sin&
      !     &(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) / Mref - 0.2D1 / 0.3D&
      !     &1 * (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y&
      !     &) / Mref) + (x * y / 0.30D2 + y / 0.30D2) * (0.2D1 / 0.3D1 * sin(a&
      !     & * x) * cos(b * y) * b * (0.10D2 - sin(a * x) * cos(b * y)) / Mref&
      !     & + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) &
      !     &* sin(b * y) * b / Mref)) - 0.3D1 / 0.4D1 * (0.2D1 + sin(a * x) * &
      !     &sin(b * y)) ** 2 / tie * (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * c&
      !     &os(b * y)) / Mref - 0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b *&
      !     & y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) * sqrt(0.&
      !     &2D1) * sqrt(0.3D1) * ((0.10D2 - sin(a * x) * cos(b * y)) / Mref) *&
      !     &* (-0.3D1 / 0.2D1)

    CASE (2)

      ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      IF (.not. switch%axisym) THEN
        WRITE (6, *) "This is an axisymmetric test case!"
        stop
      END IF
      a = 2*pi
      b = 2*pi

      D = phys%diff_n
      mu = phys%diff_u
      csii = phys%diff_e
      csie = phys%diff_ee
      kpari = phys%diff_pari
      kpare = phys%diff_pare
      Mref = phys%Mref
      epn = phys%epn
      tie = phys%tie
      pcourr = 1.

      DO i = 1, size(x)
        xx = x(i)
        yy = y(i)

            f(i, 1) = (1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy&
     &**2)**2*(-D*a*xm**4*cos(a*xx)*sin(a*yy)-D*a*xx**4*cos(a*xx)*sin(a*&
     &yy)*6.0D0+a*xx**5*cos(a*xx)*sin(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D&
     &0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0-xx*ym**3*cos(a*xx)*&
     &cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
     &)*2.0D0-xx**3*ym*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1&
     &.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xx*yy**3*cos(a*xx)*cos(a*yy)*sq&
     &rt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xx**&
     &3*yy*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(&
     &ym-yy)**2+1.0D0)*2.0D0+D*a**2*xx**5*sin(a*xx)*sin(a*yy)*6.0D0-D*a*&
     &*2*xm*xx**4*sin(a*xx)*sin(a*yy)*1.0D+1+D*a**2*xm**4*xx*sin(a*xx)*s&
     &in(a*yy)+D*a**2*xx*ym**4*sin(a*xx)*sin(a*yy)+D*a**2*xx*yy**4*sin(a&
     &*xx)*sin(a*yy)-a*xx**5*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*(1.0D0/xx*&
     &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a*xx**5&
     &*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx&
     &**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+D*a**2*xx**2*ym**3*cos(a*xx)*&
     &cos(a*yy)*2.0D0-D*a**2*xx**2*yy**3*cos(a*xx)*cos(a*yy)*2.0D0+D*a**&
     &2*xm**2*xx**3*sin(a*xx)*sin(a*yy)*9.0D0-D*a**2*xm**3*xx**2*sin(a*x&
     &x)*sin(a*yy)*4.0D0+D*a**2*xx**3*ym**2*sin(a*xx)*sin(a*yy)*5.0D0+D*&
     &a**2*xx**3*yy**2*sin(a*xx)*sin(a*yy)*5.0D0+D*a*xm*xx**3*cos(a*xx)*&
     &sin(a*yy)*1.2D+1+D*a*xm**3*xx*cos(a*xx)*sin(a*yy)*5.0D0-D*a*xm*ym*&
     &*3*cos(a*yy)*sin(a*xx)-D*a*xm**3*ym*cos(a*yy)*sin(a*xx)+D*a*xx*ym*&
     &*3*cos(a*yy)*sin(a*xx)*2.0D0+D*a*xx**3*ym*cos(a*yy)*sin(a*xx)*2.0D&
     &0+D*a*xm*yy**3*cos(a*yy)*sin(a*xx)+D*a*xm**3*yy*cos(a*yy)*sin(a*xx&
     &)-D*a*xx*yy**3*cos(a*yy)*sin(a*xx)*2.0D0-D*a*xx**3*yy*cos(a*yy)*si&
     &n(a*xx)*2.0D0-a*xm*xx**4*cos(a*xx)*sin(a*yy)*(1.0D0/xx**2*(xm-xx)*&
     &*2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D0+a*xx**4*ym*c&
     &os(a*yy)*sin(a*xx)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+&
     &1.0D0)**(3.0D0/2.0D0)*2.0D0-a*xx**4*yy*cos(a*yy)*sin(a*xx)*(1.0D0/&
     &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)*2.0D&
     &0+xm*xx**2*ym*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D&
     &0/xx**2*(ym-yy)**2+1.0D0)*4.0D0-xm**2*xx*ym*cos(a*xx)*cos(a*yy)*sq&
     &rt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0-xm*x&
     &x**2*yy*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**&
     &2*(ym-yy)**2+1.0D0)*4.0D0+xm**2*xx*yy*cos(a*xx)*cos(a*yy)*sqrt(1.0&
     &D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0-xx*ym*yy**&
     &2*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
     &yy)**2+1.0D0)*6.0D0+xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sqrt(1.0D0/xx*&
     &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*6.0D0+D*a**2*xx**4*ym*&
     &cos(a*xx)*cos(a*yy)*4.0D0-D*a**2*xx**4*yy*cos(a*xx)*cos(a*yy)*4.0D&
     &0-D*a*xm**2*xx**2*cos(a*xx)*sin(a*yy)*1.1D+1-D*a*xm**2*ym**2*cos(a&
     &*xx)*sin(a*yy)-D*a*xx**2*ym**2*cos(a*xx)*sin(a*yy)*5.0D0-D*a*xm**2&
     &*yy**2*cos(a*xx)*sin(a*yy)-D*a*xx**2*yy**2*cos(a*xx)*sin(a*yy)*5.0&
     &D0+D*a*xm*xx*ym**2*cos(a*xx)*sin(a*yy)*3.0D0-D*a*xm*xx**2*ym*cos(a&
     &*yy)*sin(a*xx)*4.0D0+D*a*xm**2*xx*ym*cos(a*yy)*sin(a*xx)*4.0D0+D*a&
     &*xm*xx*yy**2*cos(a*xx)*sin(a*yy)*3.0D0+D*a*xm*xx**2*yy*cos(a*yy)*s&
     &in(a*xx)*4.0D0-D*a*xm**2*xx*yy*cos(a*yy)*sin(a*xx)*4.0D0-D*a*xm*ym&
     &*yy**2*cos(a*yy)*sin(a*xx)*3.0D0+D*a*xm*ym**2*yy*cos(a*yy)*sin(a*x&
     &x)*3.0D0+D*a*xm**2*ym*yy*cos(a*xx)*sin(a*yy)*2.0D0+D*a*xx*ym*yy**2&
     &*cos(a*yy)*sin(a*xx)*6.0D0-D*a*xx*ym**2*yy*cos(a*yy)*sin(a*xx)*6.0&
     &D0+D*a*xx**2*ym*yy*cos(a*xx)*sin(a*yy)*1.0D+1-D*a**2*xm*xx*ym**3*c&
     &os(a*xx)*cos(a*yy)*2.0D0-D*a**2*xm*xx**3*ym*cos(a*xx)*cos(a*yy)*8.&
     &0D0-D*a**2*xm**3*xx*ym*cos(a*xx)*cos(a*yy)*2.0D0+D*a**2*xm*xx*yy**&
     &3*cos(a*xx)*cos(a*yy)*2.0D0+D*a**2*xm*xx**3*yy*cos(a*xx)*cos(a*yy)&
     &*8.0D0+D*a**2*xm**3*xx*yy*cos(a*xx)*cos(a*yy)*2.0D0-D*a**2*xx*ym*y&
     &y**3*sin(a*xx)*sin(a*yy)*4.0D0-D*a**2*xx*ym**3*yy*sin(a*xx)*sin(a*&
     &yy)*4.0D0-D*a**2*xx**3*ym*yy*sin(a*xx)*sin(a*yy)*1.0D+1-xx*ym**3*c&
     &os(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2&
     &+1.0D0/xx**2*(ym-yy)**2+1.0D0)-xx**3*ym*cos(a*xx)*cos(a*yy)*sin(a*&
     &xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1&
     &.0D0)+xx*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/&
     &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+xx**3*yy*cos(a*xx)*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx&
     &**2*(ym-yy)**2+1.0D0)+a*xm*xx**4*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*&
     &(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D&
     &0)-a*xx**4*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*(1.0D0/xx**2*(xm-xx&
     &)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+a*xx**4*yy*cos(a&
     &*xx)**2*cos(a*yy)*sin(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(y&
     &m-yy)**2+1.0D0)**(3.0D0/2.0D0)-a*xm*xx**4*cos(a*xx)*sin(a*xx)*sin(&
     &a*yy)**2*(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3&
     &.0D0/2.0D0)+a*xx**4*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*(1.0D0/xx*&
     &*2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)-a*xx**4&
     &*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*(1.0D0/xx**2*(xm-xx)**2+1.0D0&
     &/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0)+D*a**2*xm**2*xx**2*ym*cos(&
     &a*xx)*cos(a*yy)*6.0D0-D*a**2*xm**2*xx**2*yy*cos(a*xx)*cos(a*yy)*6.&
     &0D0+D*a**2*xx**2*ym*yy**2*cos(a*xx)*cos(a*yy)*6.0D0-D*a**2*xx**2*y&
     &m**2*yy*cos(a*xx)*cos(a*yy)*6.0D0-D*a**2*xm*xx**2*ym**2*sin(a*xx)*&
     &sin(a*yy)*4.0D0+D*a**2*xm**2*xx*ym**2*sin(a*xx)*sin(a*yy)*2.0D0-D*&
     &a**2*xm*xx**2*yy**2*sin(a*xx)*sin(a*yy)*4.0D0+D*a**2*xm**2*xx*yy**&
     &2*sin(a*xx)*sin(a*yy)*2.0D0+D*a**2*xx*ym**2*yy**2*sin(a*xx)*sin(a*&
     &yy)*6.0D0+D*a**2*xm*xx**2*ym*yy*sin(a*xx)*sin(a*yy)*8.0D0-D*a**2*x&
     &m**2*xx*ym*yy*sin(a*xx)*sin(a*yy)*4.0D0+xm*xx**2*ym*cos(a*xx)*cos(&
     &a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*&
     &(ym-yy)**2+1.0D0)*2.0D0-xm**2*xx*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)*&
     &sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
     &)-xm*xx**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/x&
     &x**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0+xm**2*xx*yy*co&
     &s(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+&
     &1.0D0/xx**2*(ym-yy)**2+1.0D0)-xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(&
     &a*xx)*sin(a*yy)*sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2&
     &+1.0D0)*3.0D0+xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*&
     &sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*3.0D0-D*&
     &a*xm*xx*ym*yy*cos(a*xx)*sin(a*yy)*6.0D0-D*a**2*xm*xx*ym*yy**2*cos(&
     &a*xx)*cos(a*yy)*6.0D0+D*a**2*xm*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*6.&
     &0D0))/xx

    f(i, 2) = Mref*(((xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*&
     &(ym-yy)**2+1.0D0)*(((sin(a*xx)*sin(a*yy)*2.0D0+4.0D0)*(a*cos(a*xx)&
     &*cos(a*yy)+a*cos(a*xx)**2*cos(a*yy)*sin(a*yy)))/(Mref*3.0D0)+(a*co&
     &s(a*yy)*sin(a*xx)*(cos(a*xx)**2*cos(a*yy)**2*(-1.0D0/2.0D0)+cos(a*&
     &xx)*sin(a*yy)+2.0D+1)*(2.0D0/3.0D0))/Mref-(a*cos(a*yy)*sin(a*xx)*(&
     &cos(a*yy)*sin(a*xx)-1.0D+1)*(2.0D0/3.0D0))/Mref+(a*sin(a*xx)*sin(a&
     &*yy)*(sin(a*xx)*sin(a*yy)*2.0D0+4.0D0))/(Mref*3.0D0)))/xx+((ym-yy)&
     &*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*(&
     &((sin(a*xx)*sin(a*yy)*2.0D0+4.0D0)*(a*sin(a*xx)*sin(a*yy)-a*cos(a*&
     &xx)*cos(a*yy)**2*sin(a*xx)))/(Mref*3.0D0)-(a*cos(a*xx)*sin(a*yy)*(&
     &cos(a*xx)**2*cos(a*yy)**2*(-1.0D0/2.0D0)+cos(a*xx)*sin(a*yy)+2.0D+&
     &1)*(2.0D0/3.0D0))/Mref+(a*cos(a*xx)*cos(a*yy)*(sin(a*xx)*sin(a*yy)&
     &*2.0D0+4.0D0))/(Mref*3.0D0)+(a*cos(a*xx)*sin(a*yy)*(cos(a*yy)*sin(&
     &a*xx)-1.0D+1)*(2.0D0/3.0D0))/Mref))/xx)-(a*cos(a*xx)**3*cos(a*yy)*&
     &*2*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2&
     &*(ym-yy)**2+1.0D0)+(cos(a*xx)**2*cos(a*yy)**2*(sin(a*xx)*sin(a*yy)&
     &+2.0D0)*(ym-yy)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(xm-x&
     &x)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0)*1.0D0/(1.0D0/xx**2*(xm-x&
     &x)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))/2.0D0-a*cos(a*&
     &xx)*cos(a*yy)**2*sin(a*xx)*(sin(a*xx)*sin(a*yy)+2.0D0)*(ym-yy)*1.0&
     &D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0&
     &)/xx+(xx*(mu*(sin(a*xx)*sin(a*yy)+2.0D0)*(a**2*cos(a*xx)*cos(a*yy)&
     &-((ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2&
     &+1.0D0)*((a*cos(a*xx)*sin(a*yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+&
     &1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(a**2*cos(a*xx)*cos(a*yy)*(ym-yy&
     &)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))&
     &/xx+(a**2*sin(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-x&
     &x)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+a*1.0D0/xx**2*cos(a*xx)*si&
     &n(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
     &yy)**2+1.0D0)-a*1.0D0/xx**2*cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt&
     &(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)-(a*cos(a*xx)&
     &*sin(a*yy)*(xm-xx)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(x&
     &m-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0)*1.0D0/(1.0D0/xx**2*(x&
     &m-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))/(xx*2.0D0)+&
     &(a*cos(a*yy)*sin(a*xx)*(ym-yy)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.&
     &0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0)*1.0D0/(1.&
     &0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))&
     &/(xx*2.0D0)))/xx-1.0D0/xx**2*((a*cos(a*xx)*sin(a*yy)*(xm-xx)*1.0D0&
     &/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*&
     &cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D&
     &0/xx**2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-&
     &xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+(((a*cos(a*xx)*sin(a*yy)*(xm-&
     &xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
     &))/xx-(a*cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx&
     &)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*(1.0D0/xx**2*(xm*2&
     &.0D0-xx*2.0D0)+1.0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2&
     &*2.0D0)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
     &)**(3.0D0/2.0D0))/(xx*2.0D0))+D*cos(a*xx)*cos(a*yy)*(a**2*sin(a*xx&
     &)*sin(a*yy)-((ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2&
     &*(ym-yy)**2+1.0D0)*(-(a*cos(a*yy)*sin(a*xx)*1.0D0/sqrt(1.0D0/xx**2&
     &*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(a**2*cos(a*xx)*cos(&
     &a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy&
     &)**2+1.0D0))/xx+(a**2*sin(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0&
     &/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-a*1.0D0/xx**2*&
     &cos(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D&
     &0/xx**2*(ym-yy)**2+1.0D0)+a*1.0D0/xx**2*cos(a*xx)*sin(a*yy)*(ym-yy&
     &)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)+&
     &(a*cos(a*yy)*sin(a*xx)*(xm-xx)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.&
     &0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0)*1.0D0/(1.&
     &0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))&
     &/(xx*2.0D0)-(a*cos(a*xx)*sin(a*yy)*(ym-yy)*(1.0D0/xx**2*(xm*2.0D0-&
     &xx*2.0D0)+1.0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D&
     &0)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3&
     &.0D0/2.0D0))/(xx*2.0D0)))/xx+1.0D0/xx**2*((a*cos(a*yy)*sin(a*xx)*(&
     &xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.&
     &0D0))/xx-(a*cos(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm&
     &-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D&
     &0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)-(((a*cos(a*yy)*si&
     &n(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
     &yy)**2+1.0D0))/xx-(a*cos(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/&
     &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*(1.0D0&
     &/xx**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**&
     &3*(ym-yy)**2*2.0D0)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-&
     &yy)**2+1.0D0)**(3.0D0/2.0D0))/(xx*2.0D0))+D*a*cos(a*yy)*sin(a*xx)*&
     &(a*cos(a*xx)*sin(a*yy)+(((a*cos(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/sqrt&
     &(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*cos(a&
     &*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx*&
     &*2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**&
     &2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)+a*mu*cos(a*xx)*sin(a*yy)*(a*c&
     &os(a*yy)*sin(a*xx)+(((a*cos(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0&
     &D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*cos(a*yy)&
     &*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(&
     &ym-yy)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.&
     &0D0/xx**2*(ym-yy)**2+1.0D0))/xx))+mu*(sin(a*xx)*sin(a*yy)+2.0D0)*(&
     &a*cos(a*yy)*sin(a*xx)+(((a*cos(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(&
     &1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*cos(a*&
     &yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**&
     &2*(ym-yy)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2&
     &+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)-D*cos(a*xx)*cos(a*yy)*(a*cos(a&
     &*xx)*sin(a*yy)+(((a*cos(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/x&
     &x**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*cos(a*xx)*sin&
     &(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-y&
     &y)**2+1.0D0))/xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/&
     &xx**2*(ym-yy)**2+1.0D0))/xx))/xx-mu*(sin(a*xx)*sin(a*yy)+2.0D0)*(-&
     &a**2*cos(a*xx)*cos(a*yy)+((xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**&
     &2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*((a*cos(a*yy)*sin(a*xx)*1.0D0/sqrt&
     &(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(a**2*co&
     &s(a*xx)*cos(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/&
     &xx**2*(ym-yy)**2+1.0D0))/xx+(a**2*sin(a*xx)*sin(a*yy)*(ym-yy)*1.0D&
     &0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(a&
     &*1.0D0/xx**3*cos(a*xx)*sin(a*yy)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0&
     &/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0&
     &D0))/2.0D0-(a*1.0D0/xx**3*cos(a*yy)*sin(a*xx)*(ym*2.0D0-yy*2.0D0)*&
     &(ym-yy)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0&
     &)**(3.0D0/2.0D0))/2.0D0))/xx+(1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*((a*&
     &cos(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D&
     &0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0&
     &/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(x&
     &m-xx)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)*&
     &*(3.0D0/2.0D0))/2.0D0)+D*cos(a*xx)*cos(a*yy)*(a**2*sin(a*xx)*sin(a&
     &*yy)-((xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy&
     &)**2+1.0D0)*(-(a*cos(a*xx)*sin(a*yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx&
     &)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(a**2*cos(a*xx)*cos(a*yy)*(&
     &ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.&
     &0D0))/xx+(a**2*sin(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*&
     &(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx-(a*1.0D0/xx**3*cos(a*&
     &yy)*sin(a*xx)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm-x&
     &x)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))/2.0D0+(a*1.0D0&
     &/xx**3*cos(a*xx)*sin(a*yy)*(ym*2.0D0-yy*2.0D0)*(ym-yy)*1.0D0/(1.0D&
     &0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))/2&
     &.0D0))/xx+(1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*((a*cos(a*yy)*sin(a*xx)&
     &*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+&
     &1.0D0))/xx-(a*cos(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(&
     &xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(xm-xx)*1.0D0/(1.0D0/&
     &xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)**(3.0D0/2.0D0))/2.0&
     &D0)+D*a*cos(a*xx)*sin(a*yy)*(a*cos(a*yy)*sin(a*xx)-(((a*cos(a*yy)*&
     &sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(y&
     &m-yy)**2+1.0D0))/xx-(a*cos(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D&
     &0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(xm-xx)*1.0D&
     &0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)+a&
     &*mu*cos(a*yy)*sin(a*xx)*(a*cos(a*xx)*sin(a*yy)-(((a*cos(a*xx)*sin(&
     &a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy&
     &)**2+1.0D0))/xx-(a*cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx&
     &**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)*(xm-xx)*1.0D0/sq&
     &rt(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx)+(a*co&
     &s(a*xx)**2*cos(a*yy)**3*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(&
     &xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0))/xx+(1.0D0/xx**3*cos(a*xx)&
     &**2*cos(a*yy)**2*(ym*2.0D0-yy*2.0D0)*(sin(a*xx)*sin(a*yy)+2.0D0)*(&
     &xm-xx)*1.0D0/(1.0D0/xx**2*(xm-xx)**2+1.0D0/xx**2*(ym-yy)**2+1.0D0)&
     &**(3.0D0/2.0D0))/2.0D0-(a*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*(sin(a*&
     &xx)*sin(a*yy)+2.0D0)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm-xx)**2+1.0&
     &D0/xx**2*(ym-yy)**2+1.0D0)*2.0D0)/xx

    f(i, 3) = -kpari*(((3.0D0**(-epn-1.0D0)*xx*((-cos(a*xx)**2*cos(a*yy)**2&
     &+cos(a*xx)*sin(a*yy)*2.0D0+4.0D+1)/Mref)**epn*(ym-yy)*1.0D0/(xm*xx&
     &*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(a**2*ym**&
     &3*cos(a*yy)**2*(-2.0D0)+a**2*yy**3*cos(a*yy)**2*2.0D0+a**2*xm**3*c&
     &os(a*yy)*sin(a*xx)*2.0D0-a**2*xx**3*cos(a*yy)*sin(a*xx)*4.0D0-a**2&
     &*ym**3*cos(a*xx)*sin(a*yy)*2.0D0+a**2*yy**3*cos(a*xx)*sin(a*yy)*2.&
     &0D0-a**2*xm**2*ym*cos(a*yy)**2*2.0D0-a**2*xx**2*ym*cos(a*yy)**2*4.&
     &0D0+a**2*xm**2*yy*cos(a*yy)**2*2.0D0+a**2*xx**2*yy*cos(a*yy)**2*4.&
     &0D0-a**2*ym*yy**2*cos(a*yy)**2*6.0D0+a**2*ym**2*yy*cos(a*yy)**2*6.&
     &0D0+a*ym**2*cos(a*xx)*cos(a*yy)*2.0D0+a*yy**2*cos(a*xx)*cos(a*yy)*&
     &2.0D0+a**2*ym**3*cos(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*yy**3*cos(a*&
     &xx)**2*cos(a*yy)**2*4.0D0-a*xm*ym*sin(a*xx)*sin(a*yy)*2.0D0+a*xx*y&
     &m*sin(a*xx)*sin(a*yy)*4.0D0+a*xm*yy*sin(a*xx)*sin(a*yy)*2.0D0-a*xx&
     &*yy*sin(a*xx)*sin(a*yy)*4.0D0+a**2*xm**2*ym*cos(a*xx)**2*cos(a*yy)&
     &**2*4.0D0+a**2*xx**2*ym*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*xm**2&
     &*yy*cos(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*xx**2*yy*cos(a*xx)**2*cos&
     &(a*yy)**2*8.0D0+a**2*ym*yy**2*cos(a*xx)**2*cos(a*yy)**2*1.2D+1-a**&
     &2*ym**2*yy*cos(a*xx)**2*cos(a*yy)**2*1.2D+1+a**2*xm*xx*ym*cos(a*yy&
     &)**2*4.0D0-a**2*xm*xx*yy*cos(a*yy)**2*4.0D0+a**2*xm*xx**2*cos(a*yy&
     &)*sin(a*xx)*8.0D0-a**2*xm**2*xx*cos(a*yy)*sin(a*xx)*6.0D0+a**2*xm*&
     &ym**2*cos(a*yy)*sin(a*xx)*2.0D0-a**2*xm**2*ym*cos(a*xx)*sin(a*yy)*&
     &2.0D0-a**2*xx*ym**2*cos(a*yy)*sin(a*xx)*2.0D0-a**2*xx**2*ym*cos(a*&
     &xx)*sin(a*yy)*4.0D0+a**2*xm*yy**2*cos(a*yy)*sin(a*xx)*2.0D0+a**2*x&
     &m**2*yy*cos(a*xx)*sin(a*yy)*2.0D0-a**2*xx*yy**2*cos(a*yy)*sin(a*xx&
     &)*2.0D0+a**2*xx**2*yy*cos(a*xx)*sin(a*yy)*4.0D0-a**2*ym*yy**2*cos(&
     &a*xx)*sin(a*yy)*6.0D0+a**2*ym**2*yy*cos(a*xx)*sin(a*yy)*6.0D0+a*ym&
     &**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0+a*yy**2*cos(a*xx)**2*co&
     &s(a*yy)*sin(a*yy)*2.0D0+a*xm*xx*cos(a*xx)*cos(a*yy)*2.0D0-a*ym*yy*&
     &cos(a*xx)*cos(a*yy)*4.0D0-a**2*xm*xx*ym*cos(a*xx)**2*cos(a*yy)**2*&
     &8.0D0+a**2*xm*xx*yy*cos(a*xx)**2*cos(a*yy)**2*8.0D0+a**2*xm**3*cos&
     &(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*xx**3*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*xm*xx*ym*cos(a*xx)*sin(a*yy&
     &)*4.0D0-a**2*xm*xx*yy*cos(a*xx)*sin(a*yy)*4.0D0-a**2*xm*ym*yy*cos(&
     &a*yy)*sin(a*xx)*4.0D0+a**2*xx*ym*yy*cos(a*yy)*sin(a*xx)*4.0D0+a*xm&
     &*ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*xx*ym*cos(a*xx)*cos(a&
     &*yy)**2*sin(a*xx)*4.0D0+a*xm*xx*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2&
     &.0D0-a*xm*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0+a*xx*yy*cos(a*&
     &xx)*cos(a*yy)**2*sin(a*xx)*4.0D0-a*ym*yy*cos(a*xx)**2*cos(a*yy)*si&
     &n(a*yy)*4.0D0+a**2*xm*xx**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy&
     &)*1.6D+1-a**2*xm**2*xx*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.2&
     &D+1+a**2*xm*ym**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a*&
     &*2*xx*ym**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*xm*&
     &yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*xx*yy**2*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*xm*ym*yy*cos(a*&
     &xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*xx*ym*yy*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0))/Mref-(3.0D0**(-epn-1.0D0)*a*((&
     &-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2.0D0+4.0D+1)/Mref)&
     &**epn*(ym-yy)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)*1.0D0/(xm*xx*&
     &(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(xm*cos(a*x&
     &x)*cos(a*yy)-xx*cos(a*xx)*cos(a*yy)+ym*sin(a*xx)*sin(a*yy)-yy*sin(&
     &a*xx)*sin(a*yy)-ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+xm*cos(a*xx)**&
     &2*cos(a*yy)*sin(a*yy)-xx*cos(a*xx)**2*cos(a*yy)*sin(a*yy)+yy*cos(a&
     &*xx)*cos(a*yy)**2*sin(a*xx))*2.0D0)/Mref+(3.0D0**(-epn-1.0D0)*1.0D&
     &0/Mref**2*a**2*epn*xx*sin(a*xx)*(sin(a*yy)-cos(a*xx)*cos(a*yy)**2)&
     &*((-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2.0D0+4.0D+1)/Mr&
     &ef)**(epn-1.0D0)*(ym-yy)*(xm*cos(a*xx)*cos(a*yy)-xx*cos(a*xx)*cos(&
     &a*yy)+ym*sin(a*xx)*sin(a*yy)-yy*sin(a*xx)*sin(a*yy)-ym*cos(a*xx)*c&
     &os(a*yy)**2*sin(a*xx)+xm*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-xx*cos(a&
     &*xx)**2*cos(a*yy)*sin(a*yy)+yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx))*4&
     &.0D0)/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))/&
     &xx-(3.0D0**(-epn-1.0D0)*((-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin&
     &(a*yy)*2.0D0+4.0D+1)/Mref)**epn*(xm-xx)*1.0D0/(xm*xx*(-2.0D0)-ym*y&
     &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(a**2*xm**3*cos(a*xx)**2&
     &*2.0D0-a**2*xx**3*cos(a*xx)**2*4.0D0+a**2*xm**3*cos(a*xx)*sin(a*yy&
     &)*2.0D0-a**2*xx**3*cos(a*xx)*sin(a*yy)*4.0D0-a**2*ym**3*cos(a*yy)*&
     &sin(a*xx)*2.0D0+a**2*yy**3*cos(a*yy)*sin(a*xx)*2.0D0+a**2*xm*xx**2&
     &*cos(a*xx)**2*8.0D0-a**2*xm**2*xx*cos(a*xx)**2*6.0D0+a**2*xm*ym**2&
     &*cos(a*xx)**2*2.0D0-a**2*xx*ym**2*cos(a*xx)**2*2.0D0+a**2*xm*yy**2&
     &*cos(a*xx)**2*2.0D0-a**2*xx*yy**2*cos(a*xx)**2*2.0D0-a**2*xm**3*co&
     &s(a*xx)**2*cos(a*yy)**2*4.0D0+a**2*xx**3*cos(a*xx)**2*cos(a*yy)**2&
     &*8.0D0+a*xm**2*sin(a*xx)*sin(a*yy)*2.0D0+a*xx**2*sin(a*xx)*sin(a*y&
     &y)*4.0D0-a*xm*xx*sin(a*xx)*sin(a*yy)*4.0D0-a**2*xm*xx**2*cos(a*xx)&
     &**2*cos(a*yy)**2*1.6D+1+a**2*xm**2*xx*cos(a*xx)**2*cos(a*yy)**2*1.&
     &2D+1-a**2*xm*ym**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a**2*xx*ym**2*c&
     &os(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*xm*yy**2*cos(a*xx)**2*cos(a*yy&
     &)**2*4.0D0+a**2*xx*yy**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*xm*y&
     &m*yy*cos(a*xx)**2*4.0D0+a**2*xx*ym*yy*cos(a*xx)**2*4.0D0+a**2*xm*x&
     &x**2*cos(a*xx)*sin(a*yy)*8.0D0-a**2*xm**2*xx*cos(a*xx)*sin(a*yy)*6&
     &.0D0+a**2*xm*ym**2*cos(a*xx)*sin(a*yy)*2.0D0-a**2*xm**2*ym*cos(a*y&
     &y)*sin(a*xx)*2.0D0-a**2*xx*ym**2*cos(a*xx)*sin(a*yy)*2.0D0-a**2*xx&
     &**2*ym*cos(a*yy)*sin(a*xx)*4.0D0+a**2*xm*yy**2*cos(a*xx)*sin(a*yy)&
     &*2.0D0+a**2*xm**2*yy*cos(a*yy)*sin(a*xx)*2.0D0-a**2*xx*yy**2*cos(a&
     &*xx)*sin(a*yy)*2.0D0+a**2*xx**2*yy*cos(a*yy)*sin(a*xx)*4.0D0-a**2*&
     &ym*yy**2*cos(a*yy)*sin(a*xx)*6.0D0+a**2*ym**2*yy*cos(a*yy)*sin(a*x&
     &x)*6.0D0-a*xm**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*xx**2*co&
     &s(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0-a*xm*ym*cos(a*xx)*cos(a*yy)*2&
     &.0D0+a*xx*ym*cos(a*xx)*cos(a*yy)*2.0D0+a*xm*yy*cos(a*xx)*cos(a*yy)&
     &*2.0D0-a*xx*yy*cos(a*xx)*cos(a*yy)*2.0D0+a**2*xm*ym*yy*cos(a*xx)**&
     &2*cos(a*yy)**2*8.0D0-a**2*xx*ym*yy*cos(a*xx)**2*cos(a*yy)**2*8.0D0&
     &-a**2*ym**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*yy*&
     &*3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*xm*xx*ym*cos&
     &(a*yy)*sin(a*xx)*4.0D0-a**2*xm*xx*yy*cos(a*yy)*sin(a*xx)*4.0D0-a**&
     &2*xm*ym*yy*cos(a*xx)*sin(a*yy)*4.0D0+a**2*xx*ym*yy*cos(a*xx)*sin(a&
     &*yy)*4.0D0+a*xm*xx*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0-a*xm*ym*&
     &cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0+a*xx*ym*cos(a*xx)**2*cos(a*&
     &yy)*sin(a*yy)*2.0D0+a*xm*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0&
     &-a*xx*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a**2*xm**2*ym*cos(&
     &a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*xx**2*ym*cos(a*xx)*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*xm**2*yy*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*xx**2*yy*cos(a*xx)*cos(a*yy)*si&
     &n(a*xx)*sin(a*yy)*8.0D0-a**2*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx&
     &)*sin(a*yy)*1.2D+1+a**2*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin&
     &(a*yy)*1.2D+1+a**2*xm*xx*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy&
     &)*8.0D0-a**2*xm*xx*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D&
     &0))/Mref+(3.0D0**(-epn-1.0D0)*a*(ym*2.0D0-yy*2.0D0)*((-cos(a*xx)**&
     &2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2.0D0+4.0D+1)/Mref)**epn*(xm-xx&
     &)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
     &**2*(xm*cos(a*xx)*cos(a*yy)-xx*cos(a*xx)*cos(a*yy)+ym*sin(a*xx)*si&
     &n(a*yy)-yy*sin(a*xx)*sin(a*yy)-ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)&
     &+xm*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-xx*cos(a*xx)**2*cos(a*yy)*sin&
     &(a*yy)+yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)))/Mref+(3.0D0**(-epn-1.&
     &0D0)*1.0D0/Mref**2*a**2*epn*cos(a*xx)*cos(a*yy)*(cos(a*xx)*sin(a*y&
     &y)+1.0D0)*((-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2.0D0+4&
     &.0D+1)/Mref)**(epn-1.0D0)*(xm-xx)*(xm*cos(a*xx)*cos(a*yy)-xx*cos(a&
     &*xx)*cos(a*yy)+ym*sin(a*xx)*sin(a*yy)-yy*sin(a*xx)*sin(a*yy)-ym*co&
     &s(a*xx)*cos(a*yy)**2*sin(a*xx)+xm*cos(a*xx)**2*cos(a*yy)*sin(a*yy)&
     &-xx*cos(a*xx)**2*cos(a*yy)*sin(a*yy)+yy*cos(a*xx)*cos(a*yy)**2*sin&
     &(a*xx))*4.0D0)/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2&
     &+yy**2))-((a*cos(a*xx)*cos(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(x&
     &m*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*x&
     &x)**2*sin(a*yy)**2*5.0D0-sin(a*xx)**2*sin(a*yy)**2*5.0D0+cos(a*xx)&
     &*sin(a*yy)*1.0D+2-sin(a*xx)*sin(a*yy)*1.0D+1-cos(a*xx)**3*cos(a*yy&
     &)**2*sin(a*yy)+cos(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0+cos(a*xx)*co&
     &s(a*yy)**2*sin(a*xx)**2*sin(a*yy)*2.0D0))/3.0D0-(a*cos(a*yy)*sin(a&
     &*xx)*(sin(a*xx)*sin(a*yy)+2.0D0)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(x&
     &m*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(-cos(a*&
     &xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*5.0D0+1.0D+2))/3.0D0+(1.0D&
     &0/xx**3*cos(a*xx)*cos(a*yy)*(sin(a*xx)*sin(a*yy)+2.0D0)*(ym-yy)*1.&
     &0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym*&
     &*2+yy**2))**(3.0D0/2.0D0)*(-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*si&
     &n(a*yy)*5.0D0+1.0D+2)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2))/3.0D&
     &0)/xx+(sin(a*xx)*sin(a*yy)+2.0D0)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D&
     &0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(a**2*csii*xx**4*cos(a*xx)**2*&
     &2.0D0+a**2*csii*ym**4*cos(a*xx)**2+a**2*csii*yy**4*cos(a*xx)**2-a*&
     &*2*mu*xx**4*cos(a*xx)**2*2.0D0-a**2*mu*ym**4*cos(a*xx)**2-a**2*mu*&
     &yy**4*cos(a*xx)**2+a**2*csii*xm**2*xx**2*cos(a*xx)**2+a**2*csii*xm&
     &**2*ym**2*cos(a*xx)**2+a**2*csii*xx**2*ym**2*cos(a*xx)**2*3.0D0+a*&
     &*2*csii*xm**2*yy**2*cos(a*xx)**2+a**2*csii*xx**2*yy**2*cos(a*xx)**&
     &2*3.0D0+a**2*csii*ym**2*yy**2*cos(a*xx)**2*6.0D0-a**2*mu*xm**2*xx*&
     &*2*cos(a*xx)**2-a**2*mu*xm**2*ym**2*cos(a*xx)**2-a**2*mu*xx**2*ym*&
     &*2*cos(a*xx)**2*3.0D0-a**2*mu*xm**2*yy**2*cos(a*xx)**2-a**2*mu*xx*&
     &*2*yy**2*cos(a*xx)**2*3.0D0-a**2*mu*ym**2*yy**2*cos(a*xx)**2*6.0D0&
     &-a**2*csii*xx**4*cos(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*csii*ym**4*c&
     &os(a*xx)**2*cos(a*yy)**2*2.0D0-a**2*csii*yy**4*cos(a*xx)**2*cos(a*&
     &yy)**2*2.0D0-a*csii*xm**3*sin(a*xx)*sin(a*yy)+a*csii*xx**3*sin(a*x&
     &x)*sin(a*yy)*2.0D0+a**2*mu*xx**4*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a&
     &**2*mu*ym**4*cos(a*xx)**2*cos(a*yy)**2*2.0D0+a**2*mu*yy**4*cos(a*x&
     &x)**2*cos(a*yy)**2*2.0D0+a**2*csii*xx**4*cos(a*xx)*sin(a*yy)*2.0D0&
     &+a**2*csii*ym**4*cos(a*xx)*sin(a*yy)+a**2*csii*yy**4*cos(a*xx)*sin&
     &(a*yy)-a**2*csii*xm*xx**3*cos(a*xx)**2*2.0D0-a**2*csii*ym*yy**3*co&
     &s(a*xx)**2*4.0D0-a**2*csii*ym**3*yy*cos(a*xx)**2*4.0D0+a**2*mu*xm*&
     &xx**3*cos(a*xx)**2*2.0D0+a**2*mu*ym*yy**3*cos(a*xx)**2*4.0D0+a**2*&
     &mu*ym**3*yy*cos(a*xx)**2*4.0D0-a**2*csii*xm*xx**3*cos(a*xx)*sin(a*&
     &yy)*2.0D0+a**2*csii*xm*ym**3*cos(a*yy)*sin(a*xx)+a**2*csii*xm**3*y&
     &m*cos(a*yy)*sin(a*xx)-a**2*csii*xx*ym**3*cos(a*yy)*sin(a*xx)-a**2*&
     &csii*xx**3*ym*cos(a*yy)*sin(a*xx)*2.0D0-a**2*csii*xm*yy**3*cos(a*y&
     &y)*sin(a*xx)-a**2*csii*xm**3*yy*cos(a*yy)*sin(a*xx)+a**2*csii*xx*y&
     &y**3*cos(a*yy)*sin(a*xx)+a**2*csii*xx**3*yy*cos(a*yy)*sin(a*xx)*2.&
     &0D0-a**2*csii*ym*yy**3*cos(a*xx)*sin(a*yy)*4.0D0-a**2*csii*ym**3*y&
     &y*cos(a*xx)*sin(a*yy)*4.0D0-a**2*csii*xm**2*xx**2*cos(a*xx)**2*cos&
     &(a*yy)**2*2.0D0-a**2*csii*xm**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*2.&
     &0D0-a**2*csii*xx**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*6.0D0-a**2*csi&
     &i*xm**2*yy**2*cos(a*xx)**2*cos(a*yy)**2*2.0D0-a**2*csii*xx**2*yy**&
     &2*cos(a*xx)**2*cos(a*yy)**2*6.0D0-a**2*csii*ym**2*yy**2*cos(a*xx)*&
     &*2*cos(a*yy)**2*1.2D+1+a**2*mu*xm**2*xx**2*cos(a*xx)**2*cos(a*yy)*&
     &*2*2.0D0+a**2*mu*xm**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*2.0D0+a**2*&
     &mu*xx**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*6.0D0+a**2*mu*xm**2*yy**2&
     &*cos(a*xx)**2*cos(a*yy)**2*2.0D0+a**2*mu*xx**2*yy**2*cos(a*xx)**2*&
     &cos(a*yy)**2*6.0D0+a**2*mu*ym**2*yy**2*cos(a*xx)**2*cos(a*yy)**2*1&
     &.2D+1+a*csii*xm**3*cos(a*xx)*cos(a*yy)**2*sin(a*xx)-a*csii*xx**3*c&
     &os(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a**2*csii*xm*xx*ym**2*cos(a*&
     &xx)**2*2.0D0-a**2*csii*xm*xx*yy**2*cos(a*xx)**2*2.0D0-a**2*csii*xm&
     &**2*ym*yy*cos(a*xx)**2*2.0D0-a**2*csii*xx**2*ym*yy*cos(a*xx)**2*6.&
     &0D0-a*mu*xm**3*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+a*mu*xx**3*cos(a*x&
     &x)*cos(a*yy)**2*sin(a*xx)*2.0D0+a**2*mu*xm*xx*ym**2*cos(a*xx)**2*2&
     &.0D0+a**2*mu*xm*xx*yy**2*cos(a*xx)**2*2.0D0+a**2*mu*xm**2*ym*yy*co&
     &s(a*xx)**2*2.0D0+a**2*mu*xx**2*ym*yy*cos(a*xx)**2*6.0D0+a**2*csii*&
     &xm**2*xx**2*cos(a*xx)*sin(a*yy)+a**2*csii*xm**2*ym**2*cos(a*xx)*si&
     &n(a*yy)+a**2*csii*xx**2*ym**2*cos(a*xx)*sin(a*yy)*3.0D0+a**2*csii*&
     &xm**2*yy**2*cos(a*xx)*sin(a*yy)+a**2*csii*xx**2*yy**2*cos(a*xx)*si&
     &n(a*yy)*3.0D0+a**2*csii*ym**2*yy**2*cos(a*xx)*sin(a*yy)*6.0D0+a*cs&
     &ii*xm**2*ym*cos(a*xx)*cos(a*yy)*2.0D0+a*csii*xx**2*ym*cos(a*xx)*co&
     &s(a*yy)*2.0D0-a*csii*xm**2*yy*cos(a*xx)*cos(a*yy)*2.0D0-a*csii*xx*&
     &*2*yy*cos(a*xx)*cos(a*yy)*2.0D0+a**2*csii*xm*xx**3*cos(a*xx)**2*co&
     &s(a*yy)**2*4.0D0-a*csii*xm*xx**2*sin(a*xx)*sin(a*yy)*4.0D0+a*csii*&
     &xm**2*xx*sin(a*xx)*sin(a*yy)*3.0D0+a**2*csii*ym*yy**3*cos(a*xx)**2&
     &*cos(a*yy)**2*8.0D0+a**2*csii*ym**3*yy*cos(a*xx)**2*cos(a*yy)**2*8&
     &.0D0+a*csii*xm*ym**2*sin(a*xx)*sin(a*yy)-a*csii*xx*ym**2*sin(a*xx)&
     &*sin(a*yy)+a*csii*xm*yy**2*sin(a*xx)*sin(a*yy)-a*csii*xx*yy**2*sin&
     &(a*xx)*sin(a*yy)-a**2*mu*xm*xx**3*cos(a*xx)**2*cos(a*yy)**2*4.0D0-&
     &a**2*mu*ym*yy**3*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*mu*ym**3*yy*&
     &cos(a*xx)**2*cos(a*yy)**2*8.0D0+a**2*csii*xm*xx*ym**2*cos(a*xx)**2&
     &*cos(a*yy)**2*4.0D0+a**2*csii*xm*xx*yy**2*cos(a*xx)**2*cos(a*yy)**&
     &2*4.0D0+a**2*csii*xm**2*ym*yy*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a**2&
     &*csii*xx**2*ym*yy*cos(a*xx)**2*cos(a*yy)**2*1.2D+1-a**2*mu*xm*xx*y&
     &m**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0-a**2*mu*xm*xx*yy**2*cos(a*xx)&
     &**2*cos(a*yy)**2*4.0D0-a**2*mu*xm**2*ym*yy*cos(a*xx)**2*cos(a*yy)*&
     &*2*4.0D0-a**2*mu*xx**2*ym*yy*cos(a*xx)**2*cos(a*yy)**2*1.2D+1+a**2&
     &*csii*xm*xx*ym*yy*cos(a*xx)**2*4.0D0-a**2*mu*xm*xx*ym*yy*cos(a*xx)&
     &**2*4.0D0-a**2*csii*xm*xx*ym**2*cos(a*xx)*sin(a*yy)*2.0D0+a**2*csi&
     &i*xm*xx**2*ym*cos(a*yy)*sin(a*xx)*4.0D0-a**2*csii*xm**2*xx*ym*cos(&
     &a*yy)*sin(a*xx)*3.0D0-a**2*csii*xm*xx*yy**2*cos(a*xx)*sin(a*yy)*2.&
     &0D0-a**2*csii*xm*xx**2*yy*cos(a*yy)*sin(a*xx)*4.0D0+a**2*csii*xm**&
     &2*xx*yy*cos(a*yy)*sin(a*xx)*3.0D0+a**2*csii*xm*ym*yy**2*cos(a*yy)*&
     &sin(a*xx)*3.0D0-a**2*csii*xm*ym**2*yy*cos(a*yy)*sin(a*xx)*3.0D0-a*&
     &*2*csii*xm**2*ym*yy*cos(a*xx)*sin(a*yy)*2.0D0-a**2*csii*xx*ym*yy**&
     &2*cos(a*yy)*sin(a*xx)*3.0D0+a**2*csii*xx*ym**2*yy*cos(a*yy)*sin(a*&
     &xx)*3.0D0-a**2*csii*xx**2*ym*yy*cos(a*xx)*sin(a*yy)*6.0D0+a*csii*x&
     &m*xx**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0-a*csii*xm**2*xx*cos&
     &(a*xx)*cos(a*yy)**2*sin(a*xx)*3.0D0-a*csii*xm*ym**2*cos(a*xx)*cos(&
     &a*yy)**2*sin(a*xx)+a*csii*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx&
     &)-a*csii*xm*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+a*csii*xx*yy**2&
     &*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+a*csii*xm**2*ym*cos(a*xx)**2*cos&
     &(a*yy)*sin(a*yy)*2.0D0+a*csii*xx**2*ym*cos(a*xx)**2*cos(a*yy)*sin(&
     &a*yy)*2.0D0-a*csii*xm**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0&
     &-a*csii*xx**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a*mu*xm*xx&
     &**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0+a*mu*xm**2*xx*cos(a*xx)&
     &*cos(a*yy)**2*sin(a*xx)*3.0D0+a*mu*xm*ym**2*cos(a*xx)*cos(a*yy)**2&
     &*sin(a*xx)-a*mu*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+a*mu*xm*&
     &yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)-a*mu*xx*yy**2*cos(a*xx)*cos&
     &(a*yy)**2*sin(a*xx)-a*mu*xm**2*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)&
     &*2.0D0-a*mu*xx**2*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0+a*mu*x&
     &m**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0+a*mu*xx**2*yy*cos(a&
     &*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a*csii*xm*xx*ym*cos(a*xx)*cos(a*&
     &yy)*4.0D0+a*csii*xm*xx*yy*cos(a*xx)*cos(a*yy)*4.0D0-a*csii*xm*ym*y&
     &y*sin(a*xx)*sin(a*yy)*2.0D0+a*csii*xx*ym*yy*sin(a*xx)*sin(a*yy)*2.&
     &0D0-a**2*csii*xm*xx*ym*yy*cos(a*xx)**2*cos(a*yy)**2*8.0D0+a**2*csi&
     &i*xm*ym**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0+a**2*csii&
     &*xm**3*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-a**2*csii*&
     &xx*ym**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-a**2*csii*x&
     &x**3*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*csii*xm&
     &*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-a**2*csii*xm*&
     &*3*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0+a**2*csii*xx*y&
     &y**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0+a**2*csii*xx**3&
     &*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*mu*xm*xx*ym&
     &*yy*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*mu*xm*ym**3*cos(a*xx)*cos&
     &(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-a**2*mu*xm**3*ym*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*2.0D0+a**2*mu*xx*ym**3*cos(a*xx)*cos(a*yy)&
     &*sin(a*xx)*sin(a*yy)*2.0D0+a**2*mu*xx**3*ym*cos(a*xx)*cos(a*yy)*si&
     &n(a*xx)*sin(a*yy)*4.0D0+a**2*mu*xm*yy**3*cos(a*xx)*cos(a*yy)*sin(a&
     &*xx)*sin(a*yy)*2.0D0+a**2*mu*xm**3*yy*cos(a*xx)*cos(a*yy)*sin(a*xx&
     &)*sin(a*yy)*2.0D0-a**2*mu*xx*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*s&
     &in(a*yy)*2.0D0-a**2*mu*xx**3*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(&
     &a*yy)*4.0D0+a**2*csii*xm*xx*ym*yy*cos(a*xx)*sin(a*yy)*4.0D0-a*csii&
     &*xm*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0+a*csii*xm*ym*yy*c&
     &os(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*csii*xx*ym*yy*cos(a*xx)*co&
     &s(a*yy)**2*sin(a*xx)*2.0D0+a*csii*xm*xx*yy*cos(a*xx)**2*cos(a*yy)*&
     &sin(a*yy)*4.0D0+a*mu*xm*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0&
     &D0-a*mu*xm*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0+a*mu*xx*ym&
     &*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*mu*xm*xx*yy*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*yy)*4.0D0+a**2*csii*xm*xx**2*ym*cos(a*xx)*cos(&
     &a*yy)*sin(a*xx)*sin(a*yy)*8.0D0-a**2*csii*xm**2*xx*ym*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*csii*xm*xx**2*yy*cos(a*xx)*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*csii*xm**2*xx*yy*cos(a*xx&
     &)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0+a**2*csii*xm*ym*yy**2*cos(a*&
     &xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*csii*xm*ym**2*yy*cos(&
     &a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*csii*xx*ym*yy**2*co&
     &s(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0+a**2*csii*xx*ym**2*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*mu*xm*xx**2*ym*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*mu*xm**2*xx*ym*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0+a**2*mu*xm*xx**2*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0-a**2*mu*xm**2*xx*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*mu*xm*ym*yy**2*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0+a**2*mu*xm*ym**2*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0+a**2*mu*xx*ym*yy**2*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0-a**2*mu*xx*ym**2*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*6.0D0)+D*(cos(a*xx)*sin(a*&
     &yy)+2.0D+1)*(a**2*sin(a*xx)*sin(a*yy)-((xm-xx)*1.0D0/sqrt(1.0D0/xx&
     &**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(-&
     &(a*cos(a*xx)*sin(a*yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*y&
     &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a**2*cos(a*xx)*cos(a*&
     &yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm*&
     &*2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a**2*sin(a*xx)*sin(a*yy)*(xm-xx)&
     &*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
     &0D0+ym**2+yy**2)))/xx-(a*1.0D0/xx**3*cos(a*yy)*sin(a*xx)*(ym*2.0D0&
     &-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+&
     &xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0))/2.0D0+(a*1.0D0/xx*&
     &*3*cos(a*xx)*sin(a*yy)*(ym*2.0D0-yy*2.0D0)*(ym-yy)*1.0D0/(1.0D0/xx&
     &**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(&
     &3.0D0/2.0D0))/2.0D0))/xx+(a*(ym*2.0D0-yy*2.0D0)*(xm-xx)*(xm*cos(a*&
     &yy)*sin(a*xx)-xx*cos(a*yy)*sin(a*xx)-ym*cos(a*xx)*sin(a*yy)+yy*cos&
     &(a*xx)*sin(a*yy))*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
     &0D0+ym**2+yy**2)**2)/2.0D0)+(1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm*&
     &*2+xx**2*2.0D0+ym**2+yy**2)**2*(a**2*csii*xx**5*cos(a*yy)**2*8.0D0&
     &-a**2*mu*xx**5*cos(a*yy)**2*8.0D0+a*csii*xm**4*sin(a*xx)**2*sin(a*&
     &yy)**2+a*csii*xx**4*sin(a*xx)**2*sin(a*yy)**2*4.0D0-D*a*xm**4*cos(&
     &a*xx)*sin(a*yy)*2.0D+1-D*a*xx**4*cos(a*xx)*sin(a*yy)*8.0D+1+a**2*c&
     &sii*xm**2*xx**3*cos(a*yy)**2*1.6D+1-a**2*csii*xm**3*xx**2*cos(a*yy&
     &)**2*8.0D0+a**2*csii*xx**3*ym**2*cos(a*yy)**2*4.0D0+a**2*csii*xx**&
     &3*yy**2*cos(a*yy)**2*4.0D0-a**2*mu*xm**2*xx**3*cos(a*yy)**2*1.6D+1&
     &+a**2*mu*xm**3*xx**2*cos(a*yy)**2*8.0D0-a**2*mu*xx**3*ym**2*cos(a*&
     &yy)**2*4.0D0-a**2*mu*xx**3*yy**2*cos(a*yy)**2*4.0D0-a**2*csii*xx**&
     &5*cos(a*xx)**2*cos(a*yy)**2*1.6D+1+a*csii*xm**4*sin(a*xx)*sin(a*yy&
     &)*2.0D0+a*csii*xx**4*sin(a*xx)*sin(a*yy)*8.0D0+a**2*mu*xx**5*cos(a&
     &*xx)**2*cos(a*yy)**2*1.6D+1+D*a**2*xx**5*sin(a*xx)*sin(a*yy)*8.0D+&
     &1+a**2*csii*xx**5*cos(a*xx)*sin(a*yy)*8.0D0-a**2*csii*xm*xx**4*cos&
     &(a*yy)**2*1.6D+1+a**2*csii*xm**4*xx*cos(a*yy)**2*2.0D0+a**2*mu*xm*&
     &xx**4*cos(a*yy)**2*1.6D+1-a**2*mu*xm**4*xx*cos(a*yy)**2*2.0D0-D*a*&
     &xm**4*cos(a*xx)**2*sin(a*yy)**2-D*a*xx**4*cos(a*xx)**2*sin(a*yy)**&
     &2*4.0D0-D*a**2*xm*xx**4*sin(a*xx)*sin(a*yy)*1.6D+2+D*a**2*xm**4*xx&
     &*sin(a*xx)*sin(a*yy)*2.0D+1-a**2*csii*xm*xx**4*cos(a*xx)*sin(a*yy)&
     &*1.6D+1+a**2*csii*xm**4*xx*cos(a*xx)*sin(a*yy)*2.0D0-a**2*csii*xx*&
     &*4*ym*cos(a*yy)*sin(a*xx)*4.0D0+a**2*csii*xx**4*yy*cos(a*yy)*sin(a&
     &*xx)*4.0D0-a**2*csii*xx**5*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0&
     &)*4.0D0+a**2*mu*xx**5*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*4.0&
     &D0-a**2*csii*xm**2*xx**3*cos(a*xx)**2*cos(a*yy)**2*3.2D+1+a**2*csi&
     &i*xm**3*xx**2*cos(a*xx)**2*cos(a*yy)**2*1.6D+1-a**2*csii*xx**3*ym*&
     &*2*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*csii*xx**3*yy**2*cos(a*xx)&
     &**2*cos(a*yy)**2*8.0D0+a*csii*xm**2*xx**2*sin(a*xx)*sin(a*yy)*1.6D&
     &+1+a*csii*xm**2*ym**2*sin(a*xx)*sin(a*yy)*2.0D0+a*csii*xx**2*ym**2&
     &*sin(a*xx)*sin(a*yy)*1.2D+1+a*csii*xm**2*yy**2*sin(a*xx)*sin(a*yy)&
     &*2.0D0+a*csii*xx**2*yy**2*sin(a*xx)*sin(a*yy)*1.2D+1+a**2*mu*xm**2&
     &*xx**3*cos(a*xx)**2*cos(a*yy)**2*3.2D+1-a**2*mu*xm**3*xx**2*cos(a*&
     &xx)**2*cos(a*yy)**2*1.6D+1+a**2*mu*xx**3*ym**2*cos(a*xx)**2*cos(a*&
     &yy)**2*8.0D0+a**2*mu*xx**3*yy**2*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a&
     &*csii*xm**4*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*csii*xx**4*co&
     &s(a*xx)*cos(a*yy)**2*sin(a*xx)*8.0D0-a**2*csii*xx**3*ym*yy*cos(a*y&
     &y)**2*8.0D0+a*mu*xm**4*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0+a*mu&
     &*xx**4*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*8.0D0+D*a**2*xx**2*ym**3*c&
     &os(a*xx)*cos(a*yy)*2.0D+1-D*a**2*xx**2*yy**3*cos(a*xx)*cos(a*yy)*2&
     &.0D+1+a**2*mu*xx**3*ym*yy*cos(a*yy)**2*8.0D0+D*a*xm*xx**3*cos(a*xx&
     &)**2*sin(a*yy)**2*8.0D0+D*a*xm**3*xx*cos(a*xx)**2*sin(a*yy)**2*4.0&
     &D0+D*a**2*xm**2*xx**3*sin(a*xx)*sin(a*yy)*1.6D+2-D*a**2*xm**3*xx**&
     &2*sin(a*xx)*sin(a*yy)*8.0D+1+D*a**2*xx**3*ym**2*sin(a*xx)*sin(a*yy&
     &)*4.0D+1+D*a**2*xx**3*yy**2*sin(a*xx)*sin(a*yy)*4.0D+1+a**2*csii*x&
     &m**2*xx**3*cos(a*xx)*sin(a*yy)*1.6D+1-a**2*csii*xm**3*xx**2*cos(a*&
     &xx)*sin(a*yy)*8.0D0-a**2*csii*xx**2*ym**3*cos(a*yy)*sin(a*xx)*2.0D&
     &0+a**2*csii*xx**3*ym**2*cos(a*xx)*sin(a*yy)*4.0D0+a**2*csii*xx**2*&
     &yy**3*cos(a*yy)*sin(a*xx)*2.0D0+a**2*csii*xx**3*yy**2*cos(a*xx)*si&
     &n(a*yy)*4.0D0-a*csii*xm*xx**3*sin(a*xx)**2*sin(a*yy)**2*8.0D0-a*cs&
     &ii*xm**3*xx*sin(a*xx)**2*sin(a*yy)**2*4.0D0+D*a**2*xx**5*cos(a*xx)&
     &*sin(a*xx)*sin(a*yy)**2*8.0D0+D*a*xm*xx**3*cos(a*xx)*sin(a*yy)*1.6&
     &D+2+D*a*xm**3*xx*cos(a*xx)*sin(a*yy)*8.0D+1-D*a*xm*ym**3*cos(a*yy)&
     &*sin(a*xx)*2.0D+1-D*a*xm**3*ym*cos(a*yy)*sin(a*xx)*2.0D+1+D*a*xx*y&
     &m**3*cos(a*yy)*sin(a*xx)*4.0D+1+D*a*xm*yy**3*cos(a*yy)*sin(a*xx)*2&
     &.0D+1+D*a*xm**3*yy*cos(a*yy)*sin(a*xx)*2.0D+1-D*a*xx*yy**3*cos(a*y&
     &y)*sin(a*xx)*4.0D+1-a**2*csii*xm*xx**2*ym**2*cos(a*yy)**2*4.0D0+a*&
     &*2*csii*xm**2*xx*ym**2*cos(a*yy)**2*2.0D0-a**2*csii*xm*xx**2*yy**2&
     &*cos(a*yy)**2*4.0D0+a**2*csii*xm**2*xx*yy**2*cos(a*yy)**2*2.0D0-a*&
     &csii*xm*ym**3*cos(a*xx)*cos(a*yy)*2.0D0-a*csii*xm**3*ym*cos(a*xx)*&
     &cos(a*yy)*2.0D0+a*csii*xx*ym**3*cos(a*xx)*cos(a*yy)*4.0D0+a*csii*x&
     &m*yy**3*cos(a*xx)*cos(a*yy)*2.0D0+a*csii*xm**3*yy*cos(a*xx)*cos(a*&
     &yy)*2.0D0-a*csii*xx*yy**3*cos(a*xx)*cos(a*yy)*4.0D0+a**2*csii*xx**&
     &5*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*8.0D0+a**2*mu*xm*xx**2*ym**2*co&
     &s(a*yy)**2*4.0D0-a**2*mu*xm**2*xx*ym**2*cos(a*yy)**2*2.0D0+a**2*mu&
     &*xm*xx**2*yy**2*cos(a*yy)**2*4.0D0-a**2*mu*xm**2*xx*yy**2*cos(a*yy&
     &)**2*2.0D0-D*a*xm**2*xx**2*cos(a*xx)**2*sin(a*yy)**2*8.0D0-D*a*xm*&
     &*2*ym**2*cos(a*xx)**2*sin(a*yy)**2-D*a*xx**2*ym**2*cos(a*xx)**2*si&
     &n(a*yy)**2*6.0D0-D*a*xm**2*yy**2*cos(a*xx)**2*sin(a*yy)**2-D*a*xx*&
     &*2*yy**2*cos(a*xx)**2*sin(a*yy)**2*6.0D0+a**2*csii*xm*xx**4*cos(a*&
     &xx)**2*cos(a*yy)**2*3.2D+1-a**2*csii*xm**4*xx*cos(a*xx)**2*cos(a*y&
     &y)**2*4.0D0-a*csii*xm*xx**3*sin(a*xx)*sin(a*yy)*1.6D+1-a*csii*xm**&
     &3*xx*sin(a*xx)*sin(a*yy)*8.0D0-a**2*mu*xm*xx**4*cos(a*xx)**2*cos(a&
     &*yy)**2*3.2D+1+a**2*mu*xm**4*xx*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a*&
     &csii*xm**2*xx**2*sin(a*xx)**2*sin(a*yy)**2*8.0D0+a*csii*xm**2*ym**&
     &2*sin(a*xx)**2*sin(a*yy)**2+a*csii*xx**2*ym**2*sin(a*xx)**2*sin(a*&
     &yy)**2*6.0D0+a*csii*xm**2*yy**2*sin(a*xx)**2*sin(a*yy)**2+a*csii*x&
     &x**2*yy**2*sin(a*xx)**2*sin(a*yy)**2*6.0D0+D*a**2*xx**4*ym*cos(a*x&
     &x)*cos(a*yy)*4.0D+1-D*a**2*xx**4*yy*cos(a*xx)*cos(a*yy)*4.0D+1-D*a&
     &*xm**2*xx**2*cos(a*xx)*sin(a*yy)*1.6D+2-D*a*xm**2*ym**2*cos(a*xx)*&
     &sin(a*yy)*2.0D+1-D*a*xx**2*ym**2*cos(a*xx)*sin(a*yy)*1.2D+2-D*a*xm&
     &**2*yy**2*cos(a*xx)*sin(a*yy)*2.0D+1-D*a*xx**2*yy**2*cos(a*xx)*sin&
     &(a*yy)*1.2D+2+a**2*mu*xm**2*xx**3*sin(a*xx)*sin(a*yy)*(sin(a*yy)**&
     &2-1.0D0)*8.0D0-a**2*mu*xm**3*xx**2*sin(a*xx)*sin(a*yy)*(sin(a*yy)*&
     &*2-1.0D0)*4.0D0+a**2*mu*xx**3*ym**2*sin(a*xx)*sin(a*yy)*(sin(a*yy)&
     &**2-1.0D0)*2.0D0+a**2*mu*xx**3*yy**2*sin(a*xx)*sin(a*yy)*(sin(a*yy&
     &)**2-1.0D0)*2.0D0-a*csii*xm*xx*ym**2*sin(a*xx)**2*sin(a*yy)**2*4.0&
     &D0-a*csii*xm*xx*yy**2*sin(a*xx)**2*sin(a*yy)**2*4.0D0-a*csii*xm**2&
     &*ym*yy*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a*csii*xx**2*ym*yy*sin(a*xx&
     &)**2*sin(a*yy)**2*1.2D+1-D*a**2*xm*xx**4*cos(a*xx)*sin(a*xx)*sin(a&
     &*yy)**2*1.6D+1+D*a**2*xm**4*xx*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.&
     &0D0-D*a**2*xx**4*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*2.0D0+D*a**2*&
     &xx**4*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*2.0D0+D*a*xm*xx*ym**2*co&
     &s(a*xx)*sin(a*yy)*8.0D+1+D*a*xm**2*xx*ym*cos(a*yy)*sin(a*xx)*4.0D+&
     &1+D*a*xm*xx*yy**2*cos(a*xx)*sin(a*yy)*8.0D+1-D*a*xm**2*xx*yy*cos(a&
     &*yy)*sin(a*xx)*4.0D+1-a*csii*xm**2*xx**2*cos(a*xx)*cos(a*yy)**2*si&
     &n(a*xx)*1.6D+1-D*a*xm*ym*yy**2*cos(a*yy)*sin(a*xx)*6.0D+1+D*a*xm*y&
     &m**2*yy*cos(a*yy)*sin(a*xx)*6.0D+1+D*a*xm**2*ym*yy*cos(a*xx)*sin(a&
     &*yy)*4.0D+1+D*a*xx*ym*yy**2*cos(a*yy)*sin(a*xx)*1.2D+2-D*a*xx*ym**&
     &2*yy*cos(a*yy)*sin(a*xx)*1.2D+2+D*a*xx**2*ym*yy*cos(a*xx)*sin(a*yy&
     &)*2.4D+2-a*csii*xm**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0&
     &-a*csii*xx**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.2D+1-a*csii&
     &*xm**2*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0-a*csii*xx**2*y&
     &y**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.2D+1+a**2*csii*xx**4*ym*co&
     &s(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a**2*csii*xx**4*yy*cos(a*xx)*&
     &*2*cos(a*yy)*sin(a*yy)*2.0D0+a*mu*xm**2*xx**2*cos(a*xx)*cos(a*yy)*&
     &*2*sin(a*xx)*1.6D+1+a*mu*xm**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*&
     &xx)*2.0D0+a*mu*xx**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.2D+1&
     &+a*mu*xm**2*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.0D0+a*mu*xx**&
     &2*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.2D+1+a*csii*xm**2*xx*ym&
     &*cos(a*xx)*cos(a*yy)*4.0D0-a*csii*xm**2*xx*yy*cos(a*xx)*cos(a*yy)*&
     &4.0D0-a*csii*xm*ym*yy**2*cos(a*xx)*cos(a*yy)*6.0D0+a*csii*xm*ym**2&
     &*yy*cos(a*xx)*cos(a*yy)*6.0D0+a*csii*xx*ym*yy**2*cos(a*xx)*cos(a*y&
     &y)*1.2D+1-a*csii*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*1.2D+1-a**2*csii*&
     &xm*xx**4*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*1.6D+1+a**2*csii*xm**4*x&
     &x*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.0D0-a**2*csii*xx**2*ym**3*cos&
     &(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)+a**2*csii*xx**2*yy**3*cos&
     &(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)-a**2*csii*xx**4*ym*cos(a*&
     &yy)*sin(a*xx)**2*sin(a*yy)*2.0D0+a**2*csii*xx**4*yy*cos(a*yy)*sin(&
     &a*xx)**2*sin(a*yy)*2.0D0+a**2*mu*xx**2*ym**3*cos(a*xx)**3*cos(a*yy&
     &)*(cos(a*yy)**2-1.0D0)-a**2*mu*xx**2*yy**3*cos(a*xx)**3*cos(a*yy)*&
     &(cos(a*yy)**2-1.0D0)+a**2*csii*xx**3*ym*yy*cos(a*xx)**2*cos(a*yy)*&
     &*2*1.6D+1-a*csii*xm*xx*ym**2*sin(a*xx)*sin(a*yy)*8.0D0-a*csii*xm*x&
     &x*yy**2*sin(a*xx)*sin(a*yy)*8.0D0-a*csii*xm**2*ym*yy*sin(a*xx)*sin&
     &(a*yy)*4.0D0-a*csii*xx**2*ym*yy*sin(a*xx)*sin(a*yy)*2.4D+1-a*csii*&
     &xm**4*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)-a*csii*xx**4*c&
     &os(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*4.0D0-a**2*mu*xx**3*y&
     &m*yy*cos(a*xx)**2*cos(a*yy)**2*1.6D+1+a*mu*xm**4*cos(a*xx)*cos(a*y&
     &y)**2*sin(a*xx)**2*sin(a*yy)+a*mu*xx**4*cos(a*xx)*cos(a*yy)**2*sin&
     &(a*xx)**2*sin(a*yy)*4.0D0+D*a**2*xx**2*ym**3*cos(a*xx)**2*cos(a*yy&
     &)*sin(a*yy)-D*a**2*xx**2*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-D*&
     &a**2*xm*xx*ym**3*cos(a*xx)*cos(a*yy)*2.0D+1-D*a**2*xm*xx**3*ym*cos&
     &(a*xx)*cos(a*yy)*8.0D+1-D*a**2*xm**3*xx*ym*cos(a*xx)*cos(a*yy)*2.0&
     &D+1+D*a**2*xm*xx*yy**3*cos(a*xx)*cos(a*yy)*2.0D+1+D*a**2*xm*xx**3*&
     &yy*cos(a*xx)*cos(a*yy)*8.0D+1+D*a**2*xm**3*xx*yy*cos(a*xx)*cos(a*y&
     &y)*2.0D+1+D*a**2*xm**2*xx**3*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*1.6D&
     &+1-D*a**2*xm**3*xx**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*8.0D0+D*a**&
     &2*xx**3*ym**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0D0+D*a**2*xx**3*&
     &yy**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0D0-D*a**2*xx**2*ym**3*co&
     &s(a*yy)*sin(a*xx)**2*sin(a*yy)+D*a**2*xx**2*yy**3*cos(a*yy)*sin(a*&
     &xx)**2*sin(a*yy)+a**2*csii*xx**2*ym**3*cos(a*xx)**2*cos(a*yy)*sin(&
     &a*yy)-a**2*csii*xx**2*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)+a**2*&
     &csii*xm**2*xx**3*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*1.6D+1-a**2*csii&
     &*xm**3*xx**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*8.0D0-D*a**2*xx**3*y&
     &m*yy*sin(a*xx)*sin(a*yy)*8.0D+1+a**2*csii*xx**3*ym**2*cos(a*xx)*si&
     &n(a*xx)*sin(a*yy)**2*4.0D0+a**2*csii*xx**3*yy**2*cos(a*xx)*sin(a*x&
     &x)*sin(a*yy)**2*4.0D0-a**2*csii*xx**2*ym**3*cos(a*yy)*sin(a*xx)**2&
     &*sin(a*yy)+a**2*csii*xx**2*yy**3*cos(a*yy)*sin(a*xx)**2*sin(a*yy)+&
     &a**2*csii*xm*xx*ym**3*cos(a*yy)*sin(a*xx)*2.0D0+a**2*csii*xm*xx**3&
     &*ym*cos(a*yy)*sin(a*xx)*8.0D0+a**2*csii*xm**3*xx*ym*cos(a*yy)*sin(&
     &a*xx)*2.0D0-a**2*csii*xm*xx*yy**3*cos(a*yy)*sin(a*xx)*2.0D0-a**2*c&
     &sii*xm*xx**3*yy*cos(a*yy)*sin(a*xx)*8.0D0-a**2*csii*xm**3*xx*yy*co&
     &s(a*yy)*sin(a*xx)*2.0D0-a**2*csii*xx**3*ym*yy*cos(a*xx)*sin(a*yy)*&
     &8.0D0+a**2*csii*xm*xx**4*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*&
     &8.0D0-a**2*csii*xm**4*xx*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)-&
     &a**2*mu*xm*xx**4*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*8.0D0+a*&
     &*2*mu*xm**4*xx*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)+a**2*csii*&
     &xm*xx**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*csii*xm**2*xx*&
     &ym**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a**2*csii*xm*xx**2*yy**2*cos&
     &(a*xx)**2*cos(a*yy)**2*8.0D0-a**2*csii*xm**2*xx*yy**2*cos(a*xx)**2&
     &*cos(a*yy)**2*4.0D0-a**2*csii*xx**5*cos(a*xx)**2*cos(a*yy)**2*sin(&
     &a*xx)*sin(a*yy)*1.2D+1-a**2*mu*xm*xx**2*ym**2*cos(a*xx)**2*cos(a*y&
     &y)**2*8.0D0+a**2*mu*xm**2*xx*ym**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0&
     &-a**2*mu*xm*xx**2*yy**2*cos(a*xx)**2*cos(a*yy)**2*8.0D0+a**2*mu*xm&
     &**2*xx*yy**2*cos(a*xx)**2*cos(a*yy)**2*4.0D0+a**2*mu*xx**5*cos(a*x&
     &x)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*1.2D+1+a*csii*xm*xx**3*cos(&
     &a*xx)*cos(a*yy)**2*sin(a*xx)*1.6D+1+a*csii*xm**3*xx*cos(a*xx)*cos(&
     &a*yy)**2*sin(a*xx)*8.0D0-a*csii*xm*ym**3*cos(a*xx)**2*cos(a*yy)*si&
     &n(a*yy)*2.0D0-a*csii*xm**3*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0&
     &D0+a*csii*xx*ym**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0+a*csii*x&
     &m*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0+a*csii*xm**3*yy*cos&
     &(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a*csii*xx*yy**3*cos(a*xx)**2*c&
     &os(a*yy)*sin(a*yy)*4.0D0+a**2*csii*xm*xx**2*ym*yy*cos(a*yy)**2*8.0&
     &D0-a**2*csii*xm**2*xx*ym*yy*cos(a*yy)**2*4.0D0-a*mu*xm*xx**3*cos(a&
     &*xx)*cos(a*yy)**2*sin(a*xx)*1.6D+1-a*mu*xm**3*xx*cos(a*xx)*cos(a*y&
     &y)**2*sin(a*xx)*8.0D0+a*mu*xm*ym**3*cos(a*xx)**2*cos(a*yy)*sin(a*y&
     &y)*2.0D0+a*mu*xm**3*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a*mu&
     &*xx*ym**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0-a*mu*xm*yy**3*cos&
     &(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-a*mu*xm**3*yy*cos(a*xx)**2*cos&
     &(a*yy)*sin(a*yy)*2.0D0+a*mu*xx*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*&
     &yy)*4.0D0+D*a**2*xm**2*xx**2*ym*cos(a*xx)*cos(a*yy)*6.0D+1-D*a**2*&
     &xm**2*xx**2*yy*cos(a*xx)*cos(a*yy)*6.0D+1+D*a**2*xx**2*ym*yy**2*co&
     &s(a*xx)*cos(a*yy)*6.0D+1-D*a**2*xx**2*ym**2*yy*cos(a*xx)*cos(a*yy)&
     &*6.0D+1-a**2*csii*xx**4*ym*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.&
     &0D0)*2.0D0+a**2*csii*xx**4*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2&
     &-1.0D0)*2.0D0-a**2*mu*xm*xx**2*ym*yy*cos(a*yy)**2*8.0D0+a**2*mu*xm&
     &**2*xx*ym*yy*cos(a*yy)**2*4.0D0+a**2*mu*xx**4*ym*cos(a*xx)**3*cos(&
     &a*yy)*(cos(a*yy)**2-1.0D0)*2.0D0-a**2*mu*xx**4*yy*cos(a*xx)**3*cos&
     &(a*yy)*(cos(a*yy)**2-1.0D0)*2.0D0+D*a*xm*xx*ym**2*cos(a*xx)**2*sin&
     &(a*yy)**2*4.0D0+D*a*xm*xx*yy**2*cos(a*xx)**2*sin(a*yy)**2*4.0D0+D*&
     &a*xm**2*ym*yy*cos(a*xx)**2*sin(a*yy)**2*2.0D0+D*a*xx**2*ym*yy*cos(&
     &a*xx)**2*sin(a*yy)**2*1.2D+1-D*a**2*xm*xx**2*ym**2*sin(a*xx)*sin(a&
     &*yy)*4.0D+1+D*a**2*xm**2*xx*ym**2*sin(a*xx)*sin(a*yy)*2.0D+1-D*a**&
     &2*xm*xx**2*yy**2*sin(a*xx)*sin(a*yy)*4.0D+1+D*a**2*xm**2*xx*yy**2*&
     &sin(a*xx)*sin(a*yy)*2.0D+1-a**2*csii*xm*xx**2*ym**2*cos(a*xx)*sin(&
     &a*yy)*4.0D0+a**2*csii*xm**2*xx*ym**2*cos(a*xx)*sin(a*yy)*2.0D0-a**&
     &2*csii*xm**2*xx**2*ym*cos(a*yy)*sin(a*xx)*6.0D0-a**2*csii*xm*xx**2&
     &*yy**2*cos(a*xx)*sin(a*yy)*4.0D0+a**2*csii*xm**2*xx*yy**2*cos(a*xx&
     &)*sin(a*yy)*2.0D0+a**2*csii*xm**2*xx**2*yy*cos(a*yy)*sin(a*xx)*6.0&
     &D0-a**2*csii*xx**2*ym*yy**2*cos(a*yy)*sin(a*xx)*6.0D0+a**2*csii*xx&
     &**2*ym**2*yy*cos(a*yy)*sin(a*xx)*6.0D0-a**2*csii*xm**2*xx**3*sin(a&
     &*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*8.0D0+a**2*csii*xm**3*xx**2*si&
     &n(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*4.0D0-a**2*csii*xx**3*ym**2&
     &*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0-a**2*csii*xx**3*yy&
     &**2*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0+D*a**2*xx**4*ym&
     &*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0-D*a**2*xx**4*yy*cos(a*xx)*&
     &*2*cos(a*yy)*sin(a*yy)*2.0D0+D*a**2*xm*xx**2*ym*yy*sin(a*xx)*sin(a&
     &*yy)*8.0D+1-D*a**2*xm**2*xx*ym*yy*sin(a*xx)*sin(a*yy)*4.0D+1-a**2*&
     &csii*xm*xx**2*ym**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0D0+a**2*cs&
     &ii*xm**2*xx*ym**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.0D0-a**2*csii&
     &*xm*xx**2*yy**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0D0+a**2*csii*x&
     &m**2*xx*yy**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.0D0-a**2*csii*xm*&
     &*2*xx**2*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0+a**2*csii*xm**2&
     &*xx**2*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0-a**2*csii*xx**2*y&
     &m*yy**2*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0+a**2*csii*xx**2*ym*&
     &*2*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0-a*csii*xm*ym**3*cos(a&
     &*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)-a*csii*xm**3*ym*cos(a*xx)*cos(a&
     &*yy)*sin(a*xx)*sin(a*yy)+a*csii*xx*ym**3*cos(a*xx)*cos(a*yy)*sin(a&
     &*xx)*sin(a*yy)*2.0D0+a*csii*xm*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)&
     &*sin(a*yy)+a*csii*xm**3*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)&
     &-a*csii*xx*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0+a**&
     &2*csii*xm*xx*ym*yy**2*cos(a*yy)*sin(a*xx)*6.0D0-a**2*csii*xm*xx*ym&
     &**2*yy*cos(a*yy)*sin(a*xx)*6.0D0+a**2*csii*xm*xx**2*ym*yy*cos(a*xx&
     &)*sin(a*yy)*8.0D0-a**2*csii*xm**2*xx*ym*yy*cos(a*xx)*sin(a*yy)*4.0&
     &D0+a**2*csii*xx**3*ym*yy*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*&
     &4.0D0-a**2*mu*xx**3*ym*yy*sin(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)&
     &*4.0D0+a*csii*xm*xx*ym*yy*sin(a*xx)**2*sin(a*yy)**2*8.0D0-a*csii*x&
     &m**2*xx**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*8.0D0+a**&
     &2*csii*xm*xx**4*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*2.4D&
     &+1-a**2*csii*xm**4*xx*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy&
     &)*3.0D0-a*csii*xm**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin&
     &(a*yy)-a*csii*xx**2*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(&
     &a*yy)*6.0D0-a**2*csii*xx**4*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*si&
     &n(a*yy)**2*4.0D0-a*csii*xm**2*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*x&
     &x)**2*sin(a*yy)-a*csii*xx**2*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx&
     &)**2*sin(a*yy)*6.0D0+a**2*csii*xx**4*yy*cos(a*xx)*cos(a*yy)*sin(a*&
     &xx)**2*sin(a*yy)**2*4.0D0+a*mu*xm**2*xx**2*cos(a*xx)*cos(a*yy)**2*&
     &sin(a*xx)**2*sin(a*yy)*8.0D0-a**2*mu*xm*xx**4*cos(a*xx)**2*cos(a*y&
     &y)**2*sin(a*xx)*sin(a*yy)*2.4D+1+a**2*mu*xm**4*xx*cos(a*xx)**2*cos&
     &(a*yy)**2*sin(a*xx)*sin(a*yy)*3.0D0+a*mu*xm**2*ym**2*cos(a*xx)*cos&
     &(a*yy)**2*sin(a*xx)**2*sin(a*yy)+a*mu*xx**2*ym**2*cos(a*xx)*cos(a*&
     &yy)**2*sin(a*xx)**2*sin(a*yy)*6.0D0+a**2*mu*xx**4*ym*cos(a*xx)*cos&
     &(a*yy)*sin(a*xx)**2*sin(a*yy)**2*4.0D0+a*mu*xm**2*yy**2*cos(a*xx)*&
     &cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)+a*mu*xx**2*yy**2*cos(a*xx)*cos&
     &(a*yy)**2*sin(a*xx)**2*sin(a*yy)*6.0D0-a**2*mu*xx**4*yy*cos(a*xx)*&
     &cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*4.0D0-D*a*xm*xx*ym*yy*cos(a*xx&
     &)*sin(a*yy)*1.6D+2+a*csii*xm*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a&
     &*xx)*8.0D0+a*csii*xm*xx*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*8.0&
     &D0+a*csii*xm**2*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0+a*csi&
     &i*xm**2*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0+a*csii*xx**2*&
     &ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*2.4D+1-a*csii*xm**2*xx*yy*c&
     &os(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0-a*csii*xm*ym*yy**2*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*yy)*6.0D0+a*csii*xm*ym**2*yy*cos(a*xx)**2*cos(&
     &a*yy)*sin(a*yy)*6.0D0+a*csii*xx*ym*yy**2*cos(a*xx)**2*cos(a*yy)*si&
     &n(a*yy)*1.2D+1-a*csii*xx*ym**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)&
     &*1.2D+1-a*mu*xm*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*8.0D0-a*&
     &mu*xm*xx*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*8.0D0-a*mu*xm**2*x&
     &x*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0-a*mu*xm**2*ym*yy*cos(a&
     &*xx)*cos(a*yy)**2*sin(a*xx)*4.0D0-a*mu*xx**2*ym*yy*cos(a*xx)*cos(a&
     &*yy)**2*sin(a*xx)*2.4D+1+a*mu*xm**2*xx*yy*cos(a*xx)**2*cos(a*yy)*s&
     &in(a*yy)*4.0D0+a*mu*xm*ym*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*6&
     &.0D0-a*mu*xm*ym**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*6.0D0-a*mu*&
     &xx*ym*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*1.2D+1+a*mu*xx*ym**2*&
     &yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*1.2D+1+a**2*csii*xm*xx*ym**3*c&
     &os(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)+a**2*csii*xm*xx**3*ym*c&
     &os(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*4.0D0+a**2*csii*xm**3*x&
     &x*ym*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)-a**2*csii*xm*xx*y&
     &y**3*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)-a**2*csii*xm*xx**&
     &3*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*4.0D0-a**2*csii*x&
     &m**3*xx*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)-a**2*mu*xm*&
     &xx*ym**3*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)-a**2*mu*xm*xx&
     &**3*ym*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*4.0D0-a**2*mu*x&
     &m**3*xx*ym*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)+a**2*mu*xm*&
     &xx*yy**3*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)+a**2*mu*xm*xx&
     &**3*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*4.0D0+a**2*mu*x&
     &m**3*xx*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)+a*csii*xm*x&
     &x*ym*yy*sin(a*xx)*sin(a*yy)*1.6D+1-a**2*csii*xx**4*ym*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*csii*xx**4*yy*cos(a*xx)*cos&
     &(a*yy)*sin(a*xx)*sin(a*yy)*8.0D0+a**2*mu*xx**4*ym*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*8.0D0-a**2*mu*xx**4*yy*cos(a*xx)*cos(a*yy)&
     &*sin(a*xx)*sin(a*yy)*8.0D0+a**2*csii*xm*xx**2*ym**2*sin(a*xx)*sin(&
     &a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0-a**2*csii*xm**2*xx*ym**2*sin(a*xx&
     &)*sin(a*yy)*(sin(a*yy)**2-1.0D0)+a**2*csii*xm*xx**2*yy**2*sin(a*xx&
     &)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0-a**2*csii*xm**2*xx*yy**2*si&
     &n(a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)-D*a**2*xm*xx*ym**3*cos(a*xx&
     &)**2*cos(a*yy)*sin(a*yy)-D*a**2*xm*xx**3*ym*cos(a*xx)**2*cos(a*yy)&
     &*sin(a*yy)*4.0D0-D*a**2*xm**3*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*y&
     &y)+D*a**2*xm*xx*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*yy)+D*a**2*xm*x&
     &x**3*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D0+D*a**2*xm**3*xx*yy*&
     &cos(a*xx)**2*cos(a*yy)*sin(a*yy)-a**2*mu*xm*xx**2*ym**2*sin(a*xx)*&
     &sin(a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0+a**2*mu*xm**2*xx*ym**2*sin(a*&
     &xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)-a**2*mu*xm*xx**2*yy**2*sin(a*xx&
     &)*sin(a*yy)*(sin(a*yy)**2-1.0D0)*2.0D0+a**2*mu*xm**2*xx*yy**2*sin(&
     &a*xx)*sin(a*yy)*(sin(a*yy)**2-1.0D0)+D*a**2*xm*xx*ym**3*cos(a*yy)*&
     &sin(a*xx)**2*sin(a*yy)+D*a**2*xm*xx**3*ym*cos(a*yy)*sin(a*xx)**2*s&
     &in(a*yy)*4.0D0+D*a**2*xm**3*xx*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)&
     &-D*a**2*xx**3*ym*yy*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*8.0D0-D*a**2*&
     &xm*xx*yy**3*cos(a*yy)*sin(a*xx)**2*sin(a*yy)-D*a**2*xm*xx**3*yy*co&
     &s(a*yy)*sin(a*xx)**2*sin(a*yy)*4.0D0-D*a**2*xm**3*xx*yy*cos(a*yy)*&
     &sin(a*xx)**2*sin(a*yy)-a**2*csii*xm**2*xx**3*cos(a*xx)**2*cos(a*yy&
     &)**2*sin(a*xx)*sin(a*yy)*2.4D+1+a**2*csii*xm**3*xx**2*cos(a*xx)**2&
     &*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*1.2D+1-a**2*csii*xx**2*ym**3*cos&
     &(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a**2*csii*xx**3*y&
     &m**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*6.0D0+a**2*csii&
     &*xx**2*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a&
     &**2*csii*xx**3*yy**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)&
     &*6.0D0+a**2*mu*xm**2*xx**3*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin&
     &(a*yy)*2.4D+1-a**2*mu*xm**3*xx**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*&
     &xx)*sin(a*yy)*1.2D+1+a**2*mu*xx**2*ym**3*cos(a*xx)*cos(a*yy)*sin(a&
     &*xx)**2*sin(a*yy)**2*2.0D0+a**2*mu*xx**3*ym**2*cos(a*xx)**2*cos(a*&
     &yy)**2*sin(a*xx)*sin(a*yy)*6.0D0-a**2*mu*xx**2*yy**3*cos(a*xx)*cos&
     &(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0+a**2*mu*xx**3*yy**2*cos(a*x&
     &x)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*6.0D0-a**2*csii*xm*xx*ym**3&
     &*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-a**2*csii*xm*xx**3*ym*cos(a*xx)*&
     &*2*cos(a*yy)*sin(a*yy)*4.0D0-a**2*csii*xm**3*xx*ym*cos(a*xx)**2*co&
     &s(a*yy)*sin(a*yy)+a**2*csii*xm*xx*yy**3*cos(a*xx)**2*cos(a*yy)*sin&
     &(a*yy)+a**2*csii*xm*xx**3*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*4.0D&
     &0+a**2*csii*xm**3*xx*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-a**2*csii&
     &*xm**2*xx**2*ym*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*3.0D0+&
     &a**2*csii*xm**2*xx**2*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D&
     &0)*3.0D0-a**2*csii*xx**2*ym*yy**2*cos(a*xx)**3*cos(a*yy)*(cos(a*yy&
     &)**2-1.0D0)*3.0D0+a**2*csii*xx**2*ym**2*yy*cos(a*xx)**3*cos(a*yy)*&
     &(cos(a*yy)**2-1.0D0)*3.0D0+a**2*csii*xm*xx*ym**3*cos(a*yy)*sin(a*x&
     &x)**2*sin(a*yy)+a**2*csii*xm*xx**3*ym*cos(a*yy)*sin(a*xx)**2*sin(a&
     &*yy)*4.0D0+a**2*csii*xm**3*xx*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)-&
     &a**2*csii*xx**3*ym*yy*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*8.0D0-a**2*&
     &csii*xm*xx*yy**3*cos(a*yy)*sin(a*xx)**2*sin(a*yy)-a**2*csii*xm*xx*&
     &*3*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*4.0D0-a**2*csii*xm**3*xx*yy&
     &*cos(a*yy)*sin(a*xx)**2*sin(a*yy)+a**2*mu*xm**2*xx**2*ym*cos(a*xx)&
     &**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*3.0D0-a**2*mu*xm**2*xx**2*yy*co&
     &s(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*3.0D0+a**2*mu*xx**2*ym*y&
     &y**2*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*3.0D0-a**2*mu*xx*&
     &*2*ym**2*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D0)*3.0D0-a**2&
     &*csii*xm*xx**2*ym*yy*cos(a*xx)**2*cos(a*yy)**2*1.6D+1+a**2*csii*xm&
     &**2*xx*ym*yy*cos(a*xx)**2*cos(a*yy)**2*8.0D0+a*csii*xm*xx**3*cos(a&
     &*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*8.0D0+a*csii*xm**3*xx*cos&
     &(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*4.0D0-a*csii*xm*ym**3*c&
     &os(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2-a*csii*xm**3*ym*cos(a&
     &*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2+a*csii*xx*ym**3*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*2.0D0-a**2*csii*xx**2*ym**3*c&
     &os(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a*csii*xm*yy**3*cos(a&
     &*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2+a*csii*xm**3*yy*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2-a*csii*xx*yy**3*cos(a*xx)**2*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)**2*2.0D0+a**2*csii*xx**2*yy**3*cos(a&
     &*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*mu*xm*xx**2*ym*yy*co&
     &s(a*xx)**2*cos(a*yy)**2*1.6D+1-a**2*mu*xm**2*xx*ym*yy*cos(a*xx)**2&
     &*cos(a*yy)**2*8.0D0-a*mu*xm*xx**3*cos(a*xx)*cos(a*yy)**2*sin(a*xx)&
     &**2*sin(a*yy)*8.0D0-a*mu*xm**3*xx*cos(a*xx)*cos(a*yy)**2*sin(a*xx)&
     &**2*sin(a*yy)*4.0D0+a*mu*xm*ym**3*cos(a*xx)**2*cos(a*yy)*sin(a*xx)&
     &*sin(a*yy)**2+a*mu*xm**3*ym*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a&
     &*yy)**2-a*mu*xx*ym**3*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**&
     &2*2.0D0+a**2*mu*xx**2*ym**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy&
     &)*4.0D0-a*mu*xm*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**&
     &2-a*mu*xm**3*yy*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2+a*mu&
     &*xx*yy**3*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*2.0D0-a**2&
     &*mu*xx**2*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+D*a*&
     &*2*xm**2*xx**2*ym*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0-D*a**2*xm&
     &**2*xx**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0+D*a**2*xx**2*y&
     &m*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0-D*a**2*xx**2*ym**2*&
     &yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0-D*a**2*xm*xx*ym*yy**2*co&
     &s(a*xx)*cos(a*yy)*6.0D+1+D*a**2*xm*xx*ym**2*yy*cos(a*xx)*cos(a*yy)&
     &*6.0D+1-D*a**2*xm*xx**2*ym**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0&
     &D0+D*a**2*xm**2*xx*ym**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.0D0-D*&
     &a**2*xm*xx**2*yy**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*4.0D0+D*a**2*&
     &xm**2*xx*yy**2*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*2.0D0-D*a**2*xm**2&
     &*xx**2*ym*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0+D*a**2*xm**2*xx**&
     &2*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0-D*a**2*xx**2*ym*yy**2*&
     &cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0+D*a**2*xx**2*ym**2*yy*cos(a&
     &*yy)*sin(a*xx)**2*sin(a*yy)*3.0D0-D*a*xm*ym**3*cos(a*xx)*cos(a*yy)&
     &*sin(a*xx)*sin(a*yy)-D*a*xm**3*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)*si&
     &n(a*yy)+D*a*xx*ym**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0&
     &+D*a*xm*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)+D*a*xm**3*yy&
     &*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)-D*a*xx*yy**3*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-D*a*xm*xx*ym*yy*cos(a*xx)**2*sin&
     &(a*yy)**2*8.0D0+a**2*csii*xm**2*xx**2*ym*cos(a*xx)**2*cos(a*yy)*si&
     &n(a*yy)*3.0D0-a**2*csii*xm**2*xx**2*yy*cos(a*xx)**2*cos(a*yy)*sin(&
     &a*yy)*3.0D0+a**2*csii*xx**2*ym*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*&
     &yy)*3.0D0-a**2*csii*xx**2*ym**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy&
     &)*3.0D0-a**2*mu*xm*xx**2*ym**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)&
     &*sin(a*yy)*6.0D0+a**2*mu*xm**2*xx*ym**2*cos(a*xx)**2*cos(a*yy)**2*&
     &sin(a*xx)*sin(a*yy)*3.0D0+a**2*mu*xm**2*xx**2*ym*cos(a*xx)*cos(a*y&
     &y)*sin(a*xx)**2*sin(a*yy)**2*6.0D0-a**2*mu*xm*xx**2*yy**2*cos(a*xx&
     &)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*6.0D0+a**2*mu*xm**2*xx*yy**2&
     &*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*3.0D0-a**2*mu*xm**2&
     &*xx**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0+a**2&
     &*mu*xx**2*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6&
     &.0D0-a**2*mu*xx**2*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a&
     &*yy)**2*6.0D0-a**2*csii*xm*xx*ym*yy**2*cos(a*xx)**2*cos(a*yy)*sin(&
     &a*yy)*3.0D0+a**2*csii*xm*xx*ym**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*&
     &yy)*3.0D0+a**2*csii*xm*xx**2*ym*yy*cos(a*xx)*sin(a*xx)*sin(a*yy)**&
     &2*8.0D0-a**2*csii*xm**2*xx*ym*yy*cos(a*xx)*sin(a*xx)*sin(a*yy)**2*&
     &4.0D0+a**2*csii*xm*xx*ym*yy**2*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.&
     &0D0-a**2*csii*xm*xx*ym**2*yy*cos(a*yy)*sin(a*xx)**2*sin(a*yy)*3.0D&
     &0+a*csii*xm*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)&
     &*4.0D0+a*csii*xm**2*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*y&
     &y)**2*2.0D0-a**2*csii*xm**2*xx**2*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)&
     &*sin(a*yy)*1.2D+1+a*csii*xm*xx*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a*&
     &xx)**2*sin(a*yy)*4.0D0-a*csii*xm**2*xx*yy*cos(a*xx)**2*cos(a*yy)*s&
     &in(a*xx)*sin(a*yy)**2*2.0D0+a**2*csii*xm**2*xx**2*yy*cos(a*xx)*cos&
     &(a*yy)*sin(a*xx)*sin(a*yy)*1.2D+1-a*csii*xm*ym*yy**2*cos(a*xx)**2*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)**2*3.0D0+a*csii*xm*ym**2*yy*cos(a*xx&
     &)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*3.0D0+a*csii*xm**2*ym*yy*cos&
     &(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*2.0D0+a*csii*xx*ym*yy**&
     &2*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*6.0D0-a*csii*xx*ym&
     &**2*yy*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*6.0D0+a*csii*&
     &xx**2*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*1.2D+1-a&
     &**2*csii*xx**2*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.&
     &2D+1+a**2*csii*xx**2*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*&
     &yy)*1.2D+1-a*mu*xm*xx*ym**2*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*si&
     &n(a*yy)*4.0D0-a*mu*xm**2*xx*ym*cos(a*xx)**2*cos(a*yy)*sin(a*xx)*si&
     &n(a*yy)**2*2.0D0+a**2*mu*xm**2*xx**2*ym*cos(a*xx)*cos(a*yy)*sin(a*&
     &xx)*sin(a*yy)*1.2D+1-a*mu*xm*xx*yy**2*cos(a*xx)*cos(a*yy)**2*sin(a&
     &*xx)**2*sin(a*yy)*4.0D0+a*mu*xm**2*xx*yy*cos(a*xx)**2*cos(a*yy)*si&
     &n(a*xx)*sin(a*yy)**2*2.0D0-a**2*mu*xm**2*xx**2*yy*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*1.2D+1+a*mu*xm*ym*yy**2*cos(a*xx)**2*cos(a&
     &*yy)*sin(a*xx)*sin(a*yy)**2*3.0D0-a*mu*xm*ym**2*yy*cos(a*xx)**2*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)**2*3.0D0-a*mu*xm**2*ym*yy*cos(a*xx)*co&
     &s(a*yy)**2*sin(a*xx)**2*sin(a*yy)*2.0D0-a*mu*xx*ym*yy**2*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*6.0D0+a*mu*xx*ym**2*yy*cos(a*&
     &xx)**2*cos(a*yy)*sin(a*xx)*sin(a*yy)**2*6.0D0-a*mu*xx**2*ym*yy*cos&
     &(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*1.2D+1+a**2*mu*xx**2*ym&
     &*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.2D+1-a**2*mu*xx**&
     &2*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.2D+1-a*csii*x&
     &m*xx*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.6D+1+a*mu*xm*xx*ym*y&
     &y*cos(a*xx)*cos(a*yy)**2*sin(a*xx)*1.6D+1+D*a*xm**2*xx*ym*cos(a*xx&
     &)*cos(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-D*a*xm**2*xx*yy*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*2.0D0-D*a*xm*ym*yy**2*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*3.0D0+D*a*xm*ym**2*yy*cos(a*xx)*cos(a*yy)*&
     &sin(a*xx)*sin(a*yy)*3.0D0+D*a*xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(&
     &a*xx)*sin(a*yy)*6.0D0-D*a*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx&
     &)*sin(a*yy)*6.0D0+a*csii*xm**2*xx*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)&
     &*sin(a*yy)*2.0D0-a*csii*xm**2*xx*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*&
     &sin(a*yy)*2.0D0-a*csii*xm*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*s&
     &in(a*yy)*3.0D0+a*csii*xm*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*si&
     &n(a*yy)*3.0D0+a*csii*xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin&
     &(a*yy)*6.0D0-a*csii*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(&
     &a*yy)*6.0D0-a**2*csii*xm*xx**2*ym*yy*sin(a*xx)*sin(a*yy)*(sin(a*yy&
     &)**2-1.0D0)*4.0D0+a**2*csii*xm**2*xx*ym*yy*sin(a*xx)*sin(a*yy)*(si&
     &n(a*yy)**2-1.0D0)*2.0D0+a**2*mu*xm*xx**2*ym*yy*sin(a*xx)*sin(a*yy)&
     &*(sin(a*yy)**2-1.0D0)*4.0D0-a**2*mu*xm**2*xx*ym*yy*sin(a*xx)*sin(a&
     &*yy)*(sin(a*yy)**2-1.0D0)*2.0D0+a**2*csii*xm*xx*ym**3*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0+a**2*csii*xm*xx**3*ym*cos(&
     &a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*8.0D0+a**2*csii*xm**3*xx&
     &*ym*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a**2*csii*&
     &xm*xx*yy**3*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a*&
     &*2*csii*xm*xx**3*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*&
     &8.0D0-a**2*csii*xm**3*xx*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a&
     &*yy)**2*2.0D0+a**2*csii*xx**3*ym*yy*cos(a*xx)**2*cos(a*yy)**2*sin(&
     &a*xx)*sin(a*yy)*1.2D+1-a**2*mu*xm*xx*ym**3*cos(a*xx)*cos(a*yy)*sin&
     &(a*xx)**2*sin(a*yy)**2*2.0D0-a**2*mu*xm*xx**3*ym*cos(a*xx)*cos(a*y&
     &y)*sin(a*xx)**2*sin(a*yy)**2*8.0D0-a**2*mu*xm**3*xx*ym*cos(a*xx)*c&
     &os(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0+a**2*mu*xm*xx*yy**3*cos(a&
     &*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0+a**2*mu*xm*xx**3*yy&
     &*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*8.0D0+a**2*mu*xm**3&
     &*xx*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*2.0D0-a**2*mu&
     &*xx**3*ym*yy*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*1.2D+1+&
     &a**2*csii*xm*xx*ym*yy**2*cos(a*xx)**3*cos(a*yy)*(cos(a*yy)**2-1.0D&
     &0)*3.0D0-a**2*csii*xm*xx*ym**2*yy*cos(a*xx)**3*cos(a*yy)*(cos(a*yy&
     &)**2-1.0D0)*3.0D0-a**2*mu*xm*xx*ym*yy**2*cos(a*xx)**3*cos(a*yy)*(c&
     &os(a*yy)**2-1.0D0)*3.0D0+a**2*mu*xm*xx*ym**2*yy*cos(a*xx)**3*cos(a&
     &*yy)*(cos(a*yy)**2-1.0D0)*3.0D0+a**2*csii*xm*xx*ym**3*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*csii*xm*xx**3*ym*cos(a*xx)*&
     &cos(a*yy)*sin(a*xx)*sin(a*yy)*1.6D+1+a**2*csii*xm**3*xx*ym*cos(a*x&
     &x)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*csii*xm*xx*yy**3*cos(a&
     &*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*csii*xm*xx**3*yy*cos&
     &(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.6D+1-a**2*csii*xm**3*xx*yy*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*mu*xm*xx*ym**3*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-a**2*mu*xm*xx**3*ym*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.6D+1-a**2*mu*xm**3*xx*ym&
     &*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*mu*xm*xx*yy**3&
     &*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0+a**2*mu*xm*xx**3*yy&
     &*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.6D+1+a**2*mu*xm**3*xx*y&
     &y*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*4.0D0-D*a**2*xm*xx*ym*yy&
     &**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0+D*a**2*xm*xx*ym**2*yy*c&
     &os(a*xx)**2*cos(a*yy)*sin(a*yy)*3.0D0+D*a**2*xm*xx**2*ym*yy*cos(a*&
     &xx)*sin(a*xx)*sin(a*yy)**2*8.0D0-D*a**2*xm**2*xx*ym*yy*cos(a*xx)*s&
     &in(a*xx)*sin(a*yy)**2*4.0D0+D*a**2*xm*xx*ym*yy**2*cos(a*yy)*sin(a*&
     &xx)**2*sin(a*yy)*3.0D0-D*a**2*xm*xx*ym**2*yy*cos(a*yy)*sin(a*xx)**&
     &2*sin(a*yy)*3.0D0+a**2*csii*xm*xx**2*ym**2*cos(a*xx)**2*cos(a*yy)*&
     &*2*sin(a*xx)*sin(a*yy)*6.0D0-a**2*csii*xm**2*xx*ym**2*cos(a*xx)**2&
     &*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*3.0D0-a**2*csii*xm**2*xx**2*ym*c&
     &os(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0+a**2*csii*xm*xx&
     &**2*yy**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*6.0D0-a**2&
     &*csii*xm**2*xx*yy**2*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)&
     &*3.0D0+a**2*csii*xm**2*xx**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*s&
     &in(a*yy)**2*6.0D0-a**2*csii*xx**2*ym*yy**2*cos(a*xx)*cos(a*yy)*sin&
     &(a*xx)**2*sin(a*yy)**2*6.0D0+a**2*csii*xx**2*ym**2*yy*cos(a*xx)*co&
     &s(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0+a**2*csii*xm*xx*ym*yy**2*c&
     &os(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0-a**2*csii*xm*xx&
     &*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0-a**2&
     &*csii*xm*xx**2*ym*yy*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)&
     &*1.2D+1+a**2*csii*xm**2*xx*ym*yy*cos(a*xx)**2*cos(a*yy)**2*sin(a*x&
     &x)*sin(a*yy)*6.0D0-a**2*mu*xm*xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(&
     &a*xx)**2*sin(a*yy)**2*6.0D0+a**2*mu*xm*xx*ym**2*yy*cos(a*xx)*cos(a&
     &*yy)*sin(a*xx)**2*sin(a*yy)**2*6.0D0+a**2*mu*xm*xx**2*ym*yy*cos(a*&
     &xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*1.2D+1-a**2*mu*xm**2*xx*ym&
     &*yy*cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*sin(a*yy)*6.0D0-a*csii*xm*&
     &xx*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*yy)*8.0D0+a**2*&
     &csii*xm*xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*1.2D+1&
     &-a**2*csii*xm*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*&
     &1.2D+1+a*mu*xm*xx*ym*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)**2*sin(a*&
     &yy)*8.0D0-a**2*mu*xm*xx*ym*yy**2*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin&
     &(a*yy)*1.2D+1+a**2*mu*xm*xx*ym**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)&
     &*sin(a*yy)*1.2D+1))/xx+Mref*pcourr*cos(a*xx)*cos(a*yy)*((a*sin(a*x&
     &x)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**&
     &2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)*1.0D+1+sin(a*xx)+sin(a*yy)*&
     &2.0D0-cos(a*yy)**2*sin(a*xx)*2.0D0)*(2.0D0/3.0D0))/(Mref*xx)+(a*co&
     &s(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0&
     &+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)-sin(a*yy)*5.0D0+cos(a*&
     &yy)*sin(a*xx)*sin(a*yy))*(4.0D0/3.0D0))/(Mref*xx))-D*a*cos(a*xx)*c&
     &os(a*yy)*(a*cos(a*yy)*sin(a*xx)-(a*(xm-xx)*(xm*cos(a*yy)*sin(a*xx)&
     &-xx*cos(a*yy)*sin(a*xx)-ym*cos(a*xx)*sin(a*yy)+yy*cos(a*xx)*sin(a*&
     &yy)))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))-&
     &(a**2*cos(a*yy)*sin(a*xx)*(csii*xx**2*cos(a*xx)*cos(a*yy)+csii*ym*&
     &*2*cos(a*xx)*cos(a*yy)+csii*yy**2*cos(a*xx)*cos(a*yy)-csii*xm*ym*s&
     &in(a*xx)*sin(a*yy)+csii*xx*ym*sin(a*xx)*sin(a*yy)+csii*xm*yy*sin(a&
     &*xx)*sin(a*yy)-csii*xx*yy*sin(a*xx)*sin(a*yy)+csii*xx**2*cos(a*xx)&
     &**2*cos(a*yy)*sin(a*yy)+csii*ym**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy&
     &)+csii*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-mu*xx**2*cos(a*xx)**&
     &2*cos(a*yy)*sin(a*yy)-mu*ym**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-mu&
     &*yy**2*cos(a*xx)**2*cos(a*yy)*sin(a*yy)-csii*ym*yy*cos(a*xx)*cos(a&
     &*yy)*2.0D0-mu*xm*ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+mu*xx*ym*cos(&
     &a*xx)*cos(a*yy)**2*sin(a*xx)+mu*xm*yy*cos(a*xx)*cos(a*yy)**2*sin(a&
     &*xx)-mu*xx*yy*cos(a*xx)*cos(a*yy)**2*sin(a*xx)+mu*ym*yy*cos(a*xx)*&
     &*2*cos(a*yy)*sin(a*yy)*2.0D0+csii*xm*ym*cos(a*xx)*cos(a*yy)**2*sin&
     &(a*xx)-csii*xx*ym*cos(a*xx)*cos(a*yy)**2*sin(a*xx)-csii*xm*yy*cos(&
     &a*xx)*cos(a*yy)**2*sin(a*xx)+csii*xx*yy*cos(a*xx)*cos(a*yy)**2*sin&
     &(a*xx)-csii*ym*yy*cos(a*xx)**2*cos(a*yy)*sin(a*yy)*2.0D0))/(xm*xx*&
     &(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)+(sqrt(3.0D0)*1&
     &.0D0/sqrt(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)*(sin(a*xx)*sin&
     &(a*yy)+2.0D0)**2*(-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2&
     &.0D0+cos(a*yy)*sin(a*xx)*2.0D0+2.0D+1))/(tie*(cos(a*yy)*sin(a*xx)-&
     &1.0D+1)*2.0D0)+(a*cos(a*xx)*cos(a*yy)**2*(xm-xx)*1.0D0/sqrt(1.0D0/&
     &xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*&
     &(cos(a*xx)*1.0D+1+sin(a*xx)*1.0D+2+cos(a*xx)**2*sin(a*xx)*2.0D0+co&
     &s(a*xx)**2*sin(a*yy)*4.0D0-cos(a*xx)**2*cos(a*yy)**2*sin(a*xx)*3.0&
     &D0+cos(a*xx)*sin(a*xx)*sin(a*yy)*1.0D+1))/(xx*3.0D0)+(1.0D0/xx**3*&
     &cos(a*xx)*cos(a*yy)*(ym*2.0D0-yy*2.0D0)*(sin(a*xx)*sin(a*yy)+2.0D0&
     &)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
     &*2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(-cos(a*xx)**2*cos(a*yy)**2+&
     &cos(a*xx)*sin(a*yy)*5.0D0+1.0D+2))/6.0D0-(a*cos(a*xx)*sin(a*yy)*(s&
     &in(a*xx)*sin(a*yy)+2.0D0)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-&
     &2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(-cos(a*xx)**2*&
     &cos(a*yy)**2+cos(a*xx)*sin(a*yy)*5.0D0+1.0D+2))/(xx*3.0D0)


            f(i, 4) = (-xx*(csie*(sin(a*xx)*sin(a*yy)+2.0D0)*(a**2*cos(a*yy)*sin(a*&
     &xx)-((ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm&
     &**2+xx**2*2.0D0+ym**2+yy**2))*((a*sin(a*xx)*sin(a*yy)*1.0D0/sqrt(1&
     &.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy*&
     &*2)))/xx-(a**2*cos(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*&
     &(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a&
     &**2*cos(a*yy)*sin(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.&
     &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+a*1.0D0/xx**2&
     &*cos(a*xx)*cos(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0&
     &)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))+a*1.0D0/xx**2*sin(a*&
     &xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy&
     &*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))-a*1.0D0/xx**4*cos(a*xx)*cos&
     &(a*yy)*(ym-yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**&
     &2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm*&
     &*2+ym**2+yy**2)-a*1.0D0/xx**4*sin(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/(1&
     &.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy*&
     &*2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)))/xx-(a&
     &*(ym-yy)*(ym*cos(a*xx)*cos(a*yy)-yy*cos(a*xx)*cos(a*yy)+xm*sin(a*x&
     &x)*sin(a*yy)-xx*sin(a*xx)*sin(a*yy)))/(xx*(xm*xx*(-2.0D0)-ym*yy*2.&
     &0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))+(a*(ym-yy)*(ym*cos(a*xx)*cos(a&
     &*yy)-yy*cos(a*xx)*cos(a*yy)+xm*sin(a*xx)*sin(a*yy)-xx*sin(a*xx)*si&
     &n(a*yy))*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)*1.0D0/(xm*xx*(-2.0&
     &D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/xx)+D*(cos(a*yy&
     &)*sin(a*xx)-1.0D+1)*(a**2*sin(a*xx)*sin(a*yy)-((ym-yy)*1.0D0/sqrt(&
     &1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy&
     &**2))*(-(a*cos(a*yy)*sin(a*xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0&
     &D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a**2*cos(a*xx&
     &)*cos(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2&
     &.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a**2*sin(a*xx)*sin(a*yy)&
     &*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+&
     &xx**2*2.0D0+ym**2+yy**2)))/xx-a*1.0D0/xx**2*cos(a*yy)*sin(a*xx)*(x&
     &m-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
     &*2*2.0D0+ym**2+yy**2))+a*1.0D0/xx**2*cos(a*xx)*sin(a*yy)*(ym-yy)*1&
     &.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D&
     &0+ym**2+yy**2))+a*1.0D0/xx**4*cos(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/(1&
     &.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy*&
     &*2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)-a*1.0D0&
     &/xx**4*cos(a*xx)*sin(a*yy)*(ym-yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0&
     &D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(-x&
     &m*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)))/xx+(a*(ym-yy)*(xm*cos(a*yy)*&
     &sin(a*xx)-xx*cos(a*yy)*sin(a*xx)-ym*cos(a*xx)*sin(a*yy)+yy*cos(a*x&
     &x)*sin(a*yy)))/(xx*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+y&
     &m**2+yy**2))-(a*(ym-yy)*(xm*cos(a*yy)*sin(a*xx)-xx*cos(a*yy)*sin(a&
     &*xx)-ym*cos(a*xx)*sin(a*yy)+yy*cos(a*xx)*sin(a*yy))*(-xm*xx-ym*yy*&
     &2.0D0+xm**2+ym**2+yy**2)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+x&
     &x**2*2.0D0+ym**2+yy**2)**2)/xx)-D*a*cos(a*xx)*cos(a*yy)*(a*cos(a*x&
     &x)*sin(a*yy)+(a*(ym-yy)*(xm*cos(a*yy)*sin(a*xx)-xx*cos(a*yy)*sin(a&
     &*xx)-ym*cos(a*xx)*sin(a*yy)+yy*cos(a*xx)*sin(a*yy)))/(xm*xx*(-2.0D&
     &0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))-a*csie*cos(a*xx)*si&
     &n(a*yy)*(a*cos(a*xx)*cos(a*yy)-(a*(ym-yy)*(ym*cos(a*xx)*cos(a*yy)-&
     &yy*cos(a*xx)*cos(a*yy)+xm*sin(a*xx)*sin(a*yy)-xx*sin(a*xx)*sin(a*y&
     &y)))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))+&
     &D*(cos(a*yy)*sin(a*xx)-1.0D+1)*(a*cos(a*xx)*sin(a*yy)+(a*(ym-yy)*(&
     &xm*cos(a*yy)*sin(a*xx)-xx*cos(a*yy)*sin(a*xx)-ym*cos(a*xx)*sin(a*y&
     &y)+yy*cos(a*xx)*sin(a*yy)))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
     &2*2.0D0+ym**2+yy**2))+csie*(sin(a*xx)*sin(a*yy)+2.0D0)*(a*cos(a*xx&
     &)*cos(a*yy)-(a*(ym-yy)*(ym*cos(a*xx)*cos(a*yy)-yy*cos(a*xx)*cos(a*&
     &yy)+xm*sin(a*xx)*sin(a*yy)-xx*sin(a*xx)*sin(a*yy)))/(xm*xx*(-2.0D0&
     &)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a*cos(a*xx)**2*&
     &cos(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0&
     &D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)-sin(a*yy)*5.0D0+cos(&
     &a*yy)*sin(a*xx)*sin(a*yy))*(1.0D+1/3.0D0)-a*cos(a*yy)*sin(a*xx)*(c&
     &os(a*yy)*sin(a*xx)-1.0D+1)*(sin(a*xx)*sin(a*yy)+2.0D0)*(ym-yy)*1.0&
     &D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+&
     &ym**2+yy**2))*(5.0D0/3.0D0)+1.0D0/xx**3*cos(a*xx)*cos(a*yy)*(cos(a&
     &*yy)*sin(a*xx)-1.0D+1)*(sin(a*xx)*sin(a*yy)+2.0D0)*(ym-yy)*1.0D0/(&
     &1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy&
     &**2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)*(5.0D0&
     &/3.0D0))/xx-D*(cos(a*yy)*sin(a*xx)-1.0D+1)*(a**2*sin(a*xx)*sin(a*y&
     &y)-((xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm*&
     &*2+xx**2*2.0D0+ym**2+yy**2))*(-(a*cos(a*xx)*sin(a*yy)*1.0D0/sqrt(1&
     &.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy*&
     &*2)))/xx+(a**2*cos(a*xx)*cos(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*&
     &(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a&
     &**2*sin(a*xx)*sin(a*yy)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.&
     &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx-(a*1.0D0/xx**&
     &3*cos(a*yy)*sin(a*xx)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx*&
     &*2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3&
     &.0D0/2.0D0))/2.0D0+(a*1.0D0/xx**3*cos(a*xx)*sin(a*yy)*(ym*2.0D0-yy&
     &*2.0D0)*(ym-yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm*&
     &*2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0))/2.0D0))/xx+(a*(ym*2.0&
     &D0-yy*2.0D0)*(xm-xx)*(xm*cos(a*yy)*sin(a*xx)-xx*cos(a*yy)*sin(a*xx&
     &)-ym*cos(a*xx)*sin(a*yy)+yy*cos(a*xx)*sin(a*yy))*1.0D0/(xm*xx*(-2.&
     &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/2.0D0)+csie*(s&
     &in(a*xx)*sin(a*yy)+2.0D0)*(-a**2*cos(a*yy)*sin(a*xx)+((xm-xx)*1.0D&
     &0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+y&
     &m**2+yy**2))*(-(a*cos(a*xx)*cos(a*yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*x&
     &x*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a**2*c&
     &os(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-&
     &ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx-(a**2*cos(a*xx)*si&
     &n(a*yy)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0&
     &+xm**2+xx**2*2.0D0+ym**2+yy**2)))/xx+(a*1.0D0/xx**3*sin(a*xx)*sin(&
     &a*yy)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D&
     &0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0))/2.0&
     &D0+(a*1.0D0/xx**3*cos(a*xx)*cos(a*yy)*(ym*2.0D0-yy*2.0D0)*(ym-yy)*&
     &1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+y&
     &m**2+yy**2))**(3.0D0/2.0D0))/2.0D0))/xx+(a*(ym*2.0D0-yy*2.0D0)*(xm&
     &-xx)*(ym*cos(a*xx)*cos(a*yy)-yy*cos(a*xx)*cos(a*yy)+xm*sin(a*xx)*s&
     &in(a*yy)-xx*sin(a*xx)*sin(a*yy))*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0&
     &+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/2.0D0)+(1.0D0/Mref**2*1.0D0/(x&
     &m*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(3.0D0&
     &**(-epn+1.0D0)*Mref*a*kpare*ym**4*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)&
     &*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref&
     &*a*kpare*yy**4*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.&
     &0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*xx**3*c&
     &os(a*yy)*sin(a*xx)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn&
     &*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-&
     &3.0D0**(-epn+1.0D0)*Mref*a*kpare*ym*yy**3*cos(a*xx)*cos(a*yy)*(-(c&
     &os(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*1.6D+1-3.0D0**(-epn+1.&
     &0D0)*Mref*a*kpare*ym**3*yy*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*&
     &xx)*2.0D0-2.0D+1)/Mref)**epn*1.6D+1+3.0D0**(-epn+1.0D0)*Mref*a*kpa&
     &re*xm*ym**3*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+&
     &1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm**3*ym*sin(&
     &a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.&
     &0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*xx*ym**3*sin(a*xx)*sin(a*yy)*&
     &(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn&
     &+1.0D0)*Mref*a*kpare*xx**3*ym*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin&
     &(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a*k&
     &pare*xm*yy**3*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0&
     &D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm**3*yy*si&
     &n(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*&
     &4.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xx*yy**3*sin(a*xx)*sin(a*yy&
     &)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-e&
     &pn+1.0D0)*Mref*a*kpare*xx**3*yy*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*s&
     &in(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a&
     &*kpare*xx**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D&
     &+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2&
     &+yy**2)*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*ym**2*cos(a*xx)*cos&
     &(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.&
     &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-3.0D0**(-epn&
     &+1.0D0)*Mref*a*kpare*yy**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*&
     &xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
     &*2*2.0D0+ym**2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm**2&
     &*ym**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mr&
     &ef)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xx**2*ym**2*cos(a*&
     &xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D&
     &0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm**2*yy**2*cos(a*xx)*cos(a*yy)&
     &*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-ep&
     &n+1.0D0)*Mref*a*kpare*xx**2*yy**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)&
     &*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*Mref&
     &*a*kpare*ym**2*yy**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.&
     &0D0-2.0D+1)/Mref)**epn*8.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*x&
     &m*xx**2*cos(a*yy)*sin(a*xx)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/M&
     &ref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**&
     &2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*xm**2*xx*cos(a*yy)*si&
     &n(a*xx)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2&
     &.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-3.0D0**(-ep&
     &n+1.0D0)*Mref*a**2*kpare*xx*ym**2*cos(a*yy)*sin(a*xx)*(-(cos(a*yy)&
     &*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm&
     &**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*k&
     &pare*xx**2*ym*cos(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0&
     &D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
     &2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*xx*yy**2*cos(a*&
     &yy)*sin(a*xx)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*&
     &xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0+3.0D0&
     &**(-epn+1.0D0)*Mref*a**2*kpare*xx**2*yy*cos(a*xx)*sin(a*yy)*(-(cos&
     &(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.&
     &0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0+3.0D0**(-epn+1.0D0)*Mref*&
     &a*kpare*xm*xx*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0&
     &D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
     &2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*ym*yy*cos(a*xx)*co&
     &s(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2&
     &.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-ep&
     &n+1.0D0)*Mref*a*kpare*xm*xx*ym**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)&
     &*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*8.0D0-3.0D0**(-epn+1.0D0)*Mref&
     &*a*kpare*xm*xx*yy**2*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.&
     &0D0-2.0D+1)/Mref)**epn*8.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm**&
     &2*ym*yy*cos(a*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/M&
     &ref)**epn*8.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*xx**2*ym*yy*cos(a&
     &*xx)*cos(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*8.0&
     &D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm*ym*sin(a*xx)*sin(a*yy)*(-(c&
     &os(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*&
     &2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mre&
     &f*a*kpare*xm*yy*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2&
     &.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym&
     &**2+yy**2)*2.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*xm*xx**2*ym*sin(&
     &a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.&
     &0D0-3.0D0**(-epn+2.0D0)*Mref*a*kpare*xm**2*xx*ym*sin(a*xx)*sin(a*y&
     &y)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-&
     &epn+2.0D0)*Mref*a*kpare*xm*xx**2*yy*sin(a*xx)*sin(a*yy)*(-(cos(a*y&
     &y)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*Mr&
     &ef*a*kpare*xm**2*xx*yy*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*&
     &2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*xm&
     &*ym*yy**2*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)&
     &/Mref)**epn*4.0D0-3.0D0**(-epn+2.0D0)*Mref*a*kpare*xm*ym**2*yy*sin&
     &(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4&
     &.0D0-3.0D0**(-epn+2.0D0)*Mref*a*kpare*xx*ym*yy**2*sin(a*xx)*sin(a*&
     &yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(&
     &-epn+2.0D0)*Mref*a*kpare*xx*ym**2*yy*sin(a*xx)*sin(a*yy)*(-(cos(a*&
     &yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*M&
     &ref*a**2*kpare*xm*xx*ym*cos(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)&
     &*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*&
     &2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*xm*xx&
     &*yy*cos(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)&
     &**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4&
     &.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*kpare*xx*ym*yy*cos(a*yy)*sin(a*&
     &xx)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0&
     &)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0+(3.0D0**(-epn+1&
     &.0D0)*Mref*a**2*epn*kpare*xx**3*sin(a*xx)**2*sin(a*yy)**2*(-(cos(a&
     &*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D&
     &0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0)/(cos(a*yy)*sin(a*xx)*2.0D0&
     &-2.0D+1)+3.0D0**(-epn+1.0D0)*Mref*a*kpare*xm*xx*ym*yy*cos(a*xx)*co&
     &s(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*1.6D+1+(3.&
     &0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*xx*ym**2*cos(a*xx)**2*cos(a*&
     &yy)**2*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.&
     &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0)/(cos(a*yy)*&
     &sin(a*xx)*2.0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*x&
     &x*yy**2*cos(a*xx)**2*cos(a*yy)**2*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0&
     &D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
     &2+yy**2)*4.0D0)/(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)-(3.0D0**(-epn+1&
     &.0D0)*Mref*a**2*epn*kpare*xm*xx**2*sin(a*xx)**2*sin(a*yy)**2*(-(co&
     &s(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2&
     &.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(cos(a*yy)*sin(a*xx)*2.&
     &0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*xm**2*xx*sin(&
     &a*xx)**2*sin(a*yy)**2*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**&
     &epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0&
     &D0)/(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)-(3.0D0**(-epn+1.0D0)*Mref*a&
     &**2*epn*kpare*xx*ym*yy*cos(a*xx)**2*cos(a*yy)**2*(-(cos(a*yy)*sin(&
     &a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+x&
     &x**2*2.0D0+ym**2+yy**2)*8.0D0)/(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)-&
     &(3.0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*xx**2*ym*cos(a*xx)*cos(a*&
     &yy)*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)&
     &**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8&
     &.0D0)/(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref&
     &*a**2*epn*kpare*xx**2*yy*cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*(&
     &-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*&
     &yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(cos(a*yy)*sin(a*xx&
     &)*2.0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*xm*xx*ym*&
     &cos(a*xx)*cos(a*yy)*sin(a*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0&
     &D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D&
     &0+ym**2+yy**2)*8.0D0)/(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)-(3.0D0**(&
     &-epn+1.0D0)*Mref*a**2*epn*kpare*xm*xx*yy*cos(a*xx)*cos(a*yy)*sin(a&
     &*xx)*sin(a*yy)*(-(cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)**epn*(xm&
     &*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(co&
     &s(a*yy)*sin(a*xx)*2.0D0-2.0D+1)))/(xx*9.0D0)-Mref*pcourr*cos(a*xx)&
     &*cos(a*yy)*((a*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2&
     &.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)*1.0D+&
     &1+sin(a*xx)+sin(a*yy)*2.0D0-cos(a*yy)**2*sin(a*xx)*2.0D0)*(2.0D0/3&
     &.0D0))/(Mref*xx)+(a*cos(a*xx)*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*x&
     &x*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)-&
     &sin(a*yy)*5.0D0+cos(a*yy)*sin(a*xx)*sin(a*yy))*(4.0D0/3.0D0))/(Mre&
     &f*xx))-D*a*sin(a*xx)*sin(a*yy)*(a*cos(a*yy)*sin(a*xx)-(a*(xm-xx)*(&
     &xm*cos(a*yy)*sin(a*xx)-xx*cos(a*yy)*sin(a*xx)-ym*cos(a*xx)*sin(a*y&
     &y)+yy*cos(a*xx)*sin(a*yy)))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
     &2*2.0D0+ym**2+yy**2))-a*csie*cos(a*yy)*sin(a*xx)*(a*sin(a*xx)*sin(&
     &a*yy)-(a*(xm-xx)*(ym*cos(a*xx)*cos(a*yy)-yy*cos(a*xx)*cos(a*yy)+xm&
     &*sin(a*xx)*sin(a*yy)-xx*sin(a*xx)*sin(a*yy)))/(xm*xx*(-2.0D0)-ym*y&
     &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))-(sqrt(3.0D0)*1.0D0/sqrt(-(&
     &cos(a*yy)*sin(a*xx)*2.0D0-2.0D+1)/Mref)*(sin(a*xx)*sin(a*yy)+2.0D0&
     &)**2*(-cos(a*xx)**2*cos(a*yy)**2+cos(a*xx)*sin(a*yy)*2.0D0+cos(a*y&
     &y)*sin(a*xx)*2.0D0+2.0D+1))/(tie*(cos(a*yy)*sin(a*xx)-1.0D+1)*2.0D&
     &0)+(a*cos(a*xx)*sin(a*yy)*(cos(a*yy)*sin(a*xx)-1.0D+1)*(sin(a*xx)*&
     &sin(a*yy)+2.0D0)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym&
     &*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(5.0D0/3.0D0))/xx+(a*cos&
     &(a*xx)*cos(a*yy)*sin(a*xx)*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(&
     &-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*yy)*1.0&
     &D+1+sin(a*xx)+sin(a*yy)*2.0D0-cos(a*yy)**2*sin(a*xx)*2.0D0)*(5.0D0&
     &/3.0D0))/xx-1.0D0/xx**3*cos(a*xx)*cos(a*yy)*(cos(a*yy)*sin(a*xx)-1&
     &.0D+1)*(ym*2.0D0-yy*2.0D0)*(sin(a*xx)*sin(a*yy)+2.0D0)*(xm-xx)*1.0&
     &D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
     &2+yy**2))**(3.0D0/2.0D0)*(5.0D0/6.0D0)
      END DO

      !                                                ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
      !                                                   IF (.not.switch%axisym) THEN
      !                                                      WRITE(6,*) "This is an axisymmetric test case!"
      !                                                      stop
      !                                                   END IF
      !                                                                        a = 2*pi
      !                                                                        b = 2*pi
      !                                                                        xc = 0.
      !                                                                        yc = 0.
      !                                                                        D     = phys%diff_n
      !                                                                        mu    = phys%diff_u
      !                                                                        csii  = phys%diff_e
      !                                                                        csie  = phys%diff_ee
      !                                                                        kpari = phys%diff_pari
      !                                                                        kpare = phys%diff_pare
      !                                                                        Mref  = phys%Mref
      !                                                                        epn   = phys%epn
      !                                                                        tie   = phys%tie
      !                                                                        pcourr = 1.

      !                                                                        f(:,1) = cos(a * x) ** 2 * a * sin(b * y) * cos(b * y) * (x / 0.30D2 &
      !     &- y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - (0.2D1 + sin(a * x) * sin(b &
      !     &* y)) * sin(a * x) * a * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D&
      !     &2 + 0.1D1 / 0.15D2) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * &
      !     &x) * cos(b * y) / 0.30D2 + (0.2D1 + sin(a * x) * sin(b * y)) * cos&
      !     &(a * x) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.1&
      !     &5D2) / x + sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) * (x * y &
      !     &/ 0.30D2 + y / 0.30D2) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a&
      !     & * x) * sin(b * y) * b * (x * y / 0.30D2 + y / 0.30D2) + (0.2D1 + &
      !     &sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * (x / 0.30D2 +&
      !     & 0.1D1 / 0.30D2) - D * (-sin(a * x) * a ** 2 * sin(b * y) - (x / 0&
      !     &.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * sin(b&
      !     & * y) / 0.30D2 - (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * cos(&
      !     &b * y) * b / 0.30D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15&
      !     &D2) * (cos(a * x) * a * sin(b * y) / 0.30D2 - (x / 0.30D2 - y ** 2&
      !     & / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a ** 2 * sin(b * y) + y&
      !     & * sin(a * x) * cos(b * y) * b / 0.30D2 + (x * y / 0.30D2 + y / 0.&
      !     &30D2) * cos(a * x) * a * cos(b * y) * b) + (cos(a * x) * a * sin(b&
      !     & * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * ((x / 0.&
      !     &30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * sin(b &
      !     &* y) + (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * cos(b * y) * b&
      !     &)) / x - sin(a * x) * sin(b * y) * b ** 2 - (x / 0.30D2 + 0.1D1 / &
      !     &0.30D2) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * cos(a&
      !     & * x) * a * sin(b * y) + (x * y / 0.30D2 + y / 0.30D2) * sin(a * x&
      !     &) * cos(b * y) * b) - (x * y / 0.30D2 + y / 0.30D2) * (-y * cos(a &
      !     &* x) * a * sin(b * y) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0&
      !     &.1D1 / 0.15D2) * cos(a * x) * a * cos(b * y) * b + (x / 0.30D2 + 0&
      !     &.1D1 / 0.30D2) * sin(a * x) * cos(b * y) * b - (x * y / 0.30D2 + y&
      !     & / 0.30D2) * sin(a * x) * sin(b * y) * b ** 2))

      !                                                                        f(:,2) =  cos(a * x) ** 3 * a * sin(b * y) * cos(b * y) ** 2 * (x / 0.&
      !     &30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - 0.2D1 * (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * cos(a * x) * cos(b * y) ** 2 * (x / 0.30D2 - &
      !     &y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a + (0.2D1 + sin(&
      !     &a * x) * sin(b * y)) * cos(a * x) ** 2 * cos(b * y) ** 2 / 0.30D2 &
      !     &+ (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) ** 2 * cos(b * y)&
      !     & ** 2 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) / x + sin(&
      !     &a * x) * cos(b * y) ** 3 * b * cos(a * x) ** 2 * (x * y / 0.30D2 +&
      !     & y / 0.30D2) - 0.2D1 * (0.2D1 + sin(a * x) * sin(b * y)) * cos(a *&
      !     & x) ** 2 * cos(b * y) * (x * y / 0.30D2 + y / 0.30D2) * sin(b * y)&
      !     & * b + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) ** 2 * cos(b&
      !     & * y) ** 2 * (x / 0.30D2 + 0.1D1 / 0.30D2) + Mref * ((x / 0.30D2 -&
      !     & y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (0.2D1 / 0.3D1 * cos(a * x) *&
      !     & a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) *&
      !     &* 2 * cos(b * y) ** 2 / 0.2D1) / Mref + 0.2D1 / 0.3D1 * (0.2D1 + s&
      !     &in(a * x) * sin(b * y)) * (-sin(a * x) * a * sin(b * y) + cos(a * &
      !     &x) * cos(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1 * co&
      !     &s(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) / M&
      !     &ref - 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * &
      !     &x) * a * cos(b * y) / Mref) + (x * y / 0.30D2 + y / 0.30D2) * (0.2&
      !     &D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a * x) * &
      !     &sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref + 0&
      !     &.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (cos(a * x) * c&
      !     &os(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) / M&
      !     &ref + 0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.10D2 - sin(&
      !     &a * x) * cos(b * y)) / Mref + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) &
      !     &* sin(b * y)) * sin(a * x) * sin(b * y) * b / Mref)) - mu * (-0.3D&
      !     &1 * cos(a * x) * a ** 2 * sin(b * y) * cos(b * y) * sin(a * x) - (&
      !     &0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a ** 2 * cos(b * y&
      !     &) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) &
      !     &** 2 * a * sin(b * y) * cos(b * y) - (0.2D1 + sin(a * x) * sin(b *&
      !     & y)) * sin(a * x) * a * cos(b * y)) / 0.30D2 - (x * y / 0.30D2 + y&
      !     & / 0.30D2) * (sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) - (0.2&
      !     &D1 + sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * b) / 0.3&
      !     &0D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x&
      !     &) ** 2 * a * sin(b * y) * cos(b * y) / 0.30D2 - (0.2D1 + sin(a * x&
      !     &) * sin(b * y)) * sin(a * x) * a * cos(b * y) / 0.30D2 + (x / 0.30&
      !     &D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-0.3D1 * cos(a * x) * a &
      !     &** 2 * sin(b * y) * cos(b * y) * sin(a * x) - (0.2D1 + sin(a * x) &
      !     &* sin(b * y)) * cos(a * x) * a ** 2 * cos(b * y)) + y * (sin(a * x&
      !     &) * cos(b * y) ** 2 * b * cos(a * x) - (0.2D1 + sin(a * x) * sin(b&
      !     & * y)) * cos(a * x) * sin(b * y) * b) / 0.30D2 + (x * y / 0.30D2 +&
      !     & y / 0.30D2) * (cos(a * x) ** 2 * a * cos(b * y) ** 2 * b - sin(a &
      !     &* x) ** 2 * cos(b * y) ** 2 * b * a - cos(a * x) ** 2 * a * sin(b &
      !     &* y) ** 2 * b + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a&
      !     & * sin(b * y) * b)) + (cos(a * x) ** 2 * a * sin(b * y) * cos(b * &
      !     &y) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * cos(b * &
      !     &y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * ((x / 0.30D&
      !     &2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) ** 2 * a * sin&
      !     &(b * y) * cos(b * y) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a *&
      !     & x) * a * cos(b * y)) + (x * y / 0.30D2 + y / 0.30D2) * (sin(a * x&
      !     &) * cos(b * y) ** 2 * b * cos(a * x) - (0.2D1 + sin(a * x) * sin(b&
      !     & * y)) * cos(a * x) * sin(b * y) * b))) / x - 0.3D1 * sin(a * x) *&
      !     & cos(b * y) * b ** 2 * cos(a * x) * sin(b * y) - (0.2D1 + sin(a * &
      !     &x) * sin(b * y)) * cos(a * x) * cos(b * y) * b ** 2 - (x / 0.30D2 &
      !     &+ 0.1D1 / 0.30D2) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D&
      !     &2) * (cos(a * x) ** 2 * a * sin(b * y) * cos(b * y) - (0.2D1 + sin&
      !     &(a * x) * sin(b * y)) * sin(a * x) * a * cos(b * y)) + (x * y / 0.&
      !     &30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) ** 2 * b * cos(a * x&
      !     &) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * y) * &
      !     &b)) - (x * y / 0.30D2 + y / 0.30D2) * (-y * (cos(a * x) ** 2 * a *&
      !     & sin(b * y) * cos(b * y) - (0.2D1 + sin(a * x) * sin(b * y)) * sin&
      !     &(a * x) * a * cos(b * y)) / 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2&
      !     & + 0.1D1 / 0.15D2) * (cos(a * x) ** 2 * a * cos(b * y) ** 2 * b - &
      !     &sin(a * x) ** 2 * cos(b * y) ** 2 * b * a - cos(a * x) ** 2 * a * &
      !     &sin(b * y) ** 2 * b + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * &
      !     &x) * a * sin(b * y) * b) + (x / 0.30D2 + 0.1D1 / 0.30D2) * (sin(a &
      !     &* x) * cos(b * y) ** 2 * b * cos(a * x) - (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * cos(a * x) * sin(b * y) * b) + (x * y / 0.30D2 + y / 0&
      !     &.30D2) * (-0.3D1 * sin(a * x) * cos(b * y) * b ** 2 * cos(a * x) *&
      !     & sin(b * y) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos&
      !     &(b * y) * b ** 2)))

      !       f(:,3) =  (cos(a * x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b&
      !     & * y)) - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(&
      !     &b * y) + 0.2D1 / 0.3D1 * cos(a * x) * a * sin(b * y) * (0.20D2 + c&
      !     &os(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1&
      !     &) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * (-sin(a * &
      !     &x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) * &
      !     &a)) * cos(a * x) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.&
      !     &1D1 / 0.15D2) - ((0.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos&
      !     &(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(&
      !     &b * y)) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * co&
      !     &s(b * y) ** 2 / 0.2D1)) * sin(a * x) * a * cos(b * y) * (x / 0.30D&
      !     &2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) + ((0.2D1 + sin(a * x) * sin&
      !     &(b * y)) * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0&
      !     &.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(a * x) * sin(b * y&
      !     &) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1)) * cos(a * x) * cos&
      !     &(b * y) / 0.30D2 + ((0.2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + &
      !     &cos(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * s&
      !     &in(b * y)) * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 *&
      !     & cos(b * y) ** 2 / 0.2D1)) * cos(a * x) * cos(b * y) * (x / 0.30D2&
      !     & - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) / x + (sin(a * x) * cos(b * y&
      !     &) * b * (0.20D2 + cos(a * x) * sin(b * y)) + (0.2D1 + sin(a * x) *&
      !     & sin(b * y)) * cos(a * x) * cos(b * y) * b + 0.2D1 / 0.3D1 * sin(a&
      !     & * x) * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y) - cos(a&
      !     & * x) ** 2 * cos(b * y) ** 2 / 0.2D1) + 0.2D1 / 0.3D1 * (0.2D1 + s&
      !     &in(a * x) * sin(b * y)) * (cos(a * x) * cos(b * y) * b + cos(a * x&
      !     &) ** 2 * cos(b * y) * sin(b * y) * b)) * cos(a * x) * cos(b * y) *&
      !     & (x * y / 0.30D2 + y / 0.30D2) - ((0.2D1 + sin(a * x) * sin(b * y)&
      !     &) * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.2D1 + &
      !     &sin(a * x) * sin(b * y)) * (0.20D2 + cos(a * x) * sin(b * y) - cos&
      !     &(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1)) * cos(a * x) * sin(b * y)&
      !     & * b * (x * y / 0.30D2 + y / 0.30D2) + ((0.2D1 + sin(a * x) * sin(&
      !     &b * y)) * (0.20D2 + cos(a * x) * sin(b * y)) + 0.2D1 / 0.3D1 * (0.&
      !     &2D1 + sin(a * x) * sin(b * y)) * (0.20D2 + cos(a * x) * sin(b * y)&
      !     & - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1)) * cos(a * x) * cos(&
      !     &b * y) * (x / 0.30D2 + 0.1D1 / 0.30D2) - csii * (-sin(a * x) * a *&
      !     &* 2 * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)) - 0.2D1 * co&
      !     &s(a * x) * a ** 2 * sin(b * y) ** 2 * sin(a * x) - (0.2D1 + sin(a &
      !     &* x) * sin(b * y)) * cos(a * x) * a ** 2 * sin(b * y) - (x / 0.30D&
      !     &2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * &
      !     &y) * (0.20D2 + cos(a * x) * sin(b * y)) - (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * sin(a * x) * a * sin(b * y)) / 0.30D2 - (x * y / 0.30D&
      !     &2 + y / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a *&
      !     & x) * sin(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x)&
      !     & * cos(b * y) * b) / 0.30D2 - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D&
      !     &1 / 0.15D2) * (cos(a * x) * a * sin(b * y) * (0.20D2 + cos(a * x) &
      !     &* sin(b * y)) / 0.30D2 - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a&
      !     & * x) * a * sin(b * y) / 0.30D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + &
      !     &0.1D1 / 0.15D2) * (-sin(a * x) * a ** 2 * sin(b * y) * (0.20D2 + c&
      !     &os(a * x) * sin(b * y)) - 0.2D1 * cos(a * x) * a ** 2 * sin(b * y)&
      !     & ** 2 * sin(a * x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x&
      !     &) * a ** 2 * sin(b * y)) + y * (sin(a * x) * cos(b * y) * b * (0.2&
      !     &0D2 + cos(a * x) * sin(b * y)) + (0.2D1 + sin(a * x) * sin(b * y))&
      !     & * cos(a * x) * cos(b * y) * b) / 0.30D2 + (x * y / 0.30D2 + y / 0&
      !     &.30D2) * (cos(a * x) * a * cos(b * y) * b * (0.20D2 + cos(a * x) *&
      !     & sin(b * y)) - sin(a * x) ** 2 * cos(b * y) * b * a * sin(b * y) +&
      !     & cos(a * x) ** 2 * a * sin(b * y) * cos(b * y) * b - (0.2D1 + sin(&
      !     &a * x) * sin(b * y)) * sin(a * x) * a * cos(b * y) * b)) + (cos(a &
      !     &* x) * a * sin(b * y) * (0.20D2 + cos(a * x) * sin(b * y)) - (0.2D&
      !     &1 + sin(a * x) * sin(b * y)) * sin(a * x) * a * sin(b * y) - (x / &
      !     &0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * ((x / 0.30D2 - y ** 2&
      !     & / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) * (0.20&
      !     &D2 + cos(a * x) * sin(b * y)) - (0.2D1 + sin(a * x) * sin(b * y)) &
      !     &* sin(a * x) * a * sin(b * y)) + (x * y / 0.30D2 + y / 0.30D2) * (&
      !     &sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)) +&
      !     & (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b))&
      !     &) / x - sin(a * x) * sin(b * y) * b ** 2 * (0.20D2 + cos(a * x) * &
      !     &sin(b * y)) + 0.2D1 * sin(a * x) * cos(b * y) ** 2 * b ** 2 * cos(&
      !     &a * x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * sin(b * &
      !     &y) * b ** 2 - (x / 0.30D2 + 0.1D1 / 0.30D2) * ((x / 0.30D2 - y ** &
      !     &2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) * (0.2&
      !     &0D2 + cos(a * x) * sin(b * y)) - (0.2D1 + sin(a * x) * sin(b * y))&
      !     & * sin(a * x) * a * sin(b * y)) + (x * y / 0.30D2 + y / 0.30D2) * &
      !     &(sin(a * x) * cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)) &
      !     &+ (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b)&
      !     &) - (x * y / 0.30D2 + y / 0.30D2) * (-y * (cos(a * x) * a * sin(b &
      !     &* y) * (0.20D2 + cos(a * x) * sin(b * y)) - (0.2D1 + sin(a * x) * &
      !     &sin(b * y)) * sin(a * x) * a * sin(b * y)) / 0.15D2 + (x / 0.30D2 &
      !     &- y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * cos(b * y)&
      !     & * b * (0.20D2 + cos(a * x) * sin(b * y)) - sin(a * x) ** 2 * cos(&
      !     &b * y) * b * a * sin(b * y) + cos(a * x) ** 2 * a * sin(b * y) * c&
      !     &os(b * y) * b - (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a&
      !     & * cos(b * y) * b) + (x / 0.30D2 + 0.1D1 / 0.30D2) * (sin(a * x) *&
      !     & cos(b * y) * b * (0.20D2 + cos(a * x) * sin(b * y)) + (0.2D1 + si&
      !     &n(a * x) * sin(b * y)) * cos(a * x) * cos(b * y) * b) + (x * y / 0&
      !     &.30D2 + y / 0.30D2) * (-sin(a * x) * sin(b * y) * b ** 2 * (0.20D2&
      !     & + cos(a * x) * sin(b * y)) + 0.2D1 * sin(a * x) * cos(b * y) ** 2&
      !     & * b ** 2 * cos(a * x) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a&
      !     & * x) * sin(b * y) * b ** 2))) - kpari * ((0.2D1 / 0.3D1 * (0.20D2&
      !     & + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0&
      !     &.2D1) / Mref) ** epn * epn * (-sin(a * x) * a * sin(b * y) + cos(a&
      !     & * x) * cos(b * y) ** 2 * sin(a * x) * a) / (0.20D2 + cos(a * x) *&
      !     & sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) * (0.2D1 &
      !     &/ 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(&
      !     &a * x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x&
      !     &) * a) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * (c&
      !     &os(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b &
      !     &* y) * b) / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2&
      !     &) + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x&
      !     &) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * ((-sin(a * x) *&
      !     & a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) * a) /&
      !     & Mref / 0.45D2 + 0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0&
      !     &.1D1 / 0.15D2) * (-cos(a * x) * a ** 2 * sin(b * y) - sin(a * x) *&
      !     &* 2 * a ** 2 * cos(b * y) ** 2 + cos(a * x) ** 2 * cos(b * y) ** 2&
      !     & * a ** 2) / Mref + y * (cos(a * x) * cos(b * y) * b + cos(a * x) &
      !     &** 2 * cos(b * y) * sin(b * y) * b) / Mref / 0.45D2 + 0.2D1 / 0.3D&
      !     &1 * (x * y / 0.30D2 + y / 0.30D2) * (-sin(a * x) * a * cos(b * y) &
      !     &* b - 0.2D1 * cos(a * x) * cos(b * y) * sin(b * y) * b * sin(a * x&
      !     &) * a) / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) +&
      !     & (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) *&
      !     &* 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * (0.2D1 / 0.3D1 * (&
      !     &x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a * x) * a &
      !     &* sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin(a * x) * a) / Mr&
      !     &ef + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * (cos(a * x) *&
      !     & cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) /&
      !     & Mref) / 0.30D2 + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * &
      !     &y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * (&
      !     &0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * &
      !     &(-sin(a * x) * a * sin(b * y) + cos(a * x) * cos(b * y) ** 2 * sin&
      !     &(a * x) * a) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2&
      !     &) * (cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * &
      !     &sin(b * y) * b) / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / &
      !     &0.15D2) / x + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin(b * y) -&
      !     & cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** epn * epn *&
      !     & (cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b * y) * sin&
      !     &(b * y) * b) / (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2&
      !     & * cos(b * y) ** 2 / 0.2D1) * (0.2D1 / 0.3D1 * (x / 0.30D2 - y ** &
      !     &2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a * x) * a * sin(b * y) + cos&
      !     &(a * x) * cos(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1&
      !     & * (x * y / 0.30D2 + y / 0.30D2) * (cos(a * x) * cos(b * y) * b + &
      !     &cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) / Mref) * (x * y / &
      !     &0.30D2 + y / 0.30D2) + (0.2D1 / 0.3D1 * (0.20D2 + cos(a * x) * sin&
      !     &(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0.2D1) / Mref) ** ep&
      !     &n * (-0.2D1 / 0.45D2 * y * (-sin(a * x) * a * sin(b * y) + cos(a *&
      !     & x) * cos(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1 * (&
      !     &x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a * x) * a &
      !     &* cos(b * y) * b - 0.2D1 * cos(a * x) * cos(b * y) * sin(b * y) * &
      !     &b * sin(a * x) * a) / Mref + 0.2D1 / 0.3D1 * (x / 0.30D2 + 0.1D1 /&
      !     & 0.30D2) * (cos(a * x) * cos(b * y) * b + cos(a * x) ** 2 * cos(b &
      !     &* y) * sin(b * y) * b) / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + &
      !     &y / 0.30D2) * (-cos(a * x) * sin(b * y) * b ** 2 - cos(a * x) ** 2&
      !     & * sin(b * y) ** 2 * b ** 2 + cos(a * x) ** 2 * cos(b * y) ** 2 * &
      !     &b ** 2) / Mref) * (x * y / 0.30D2 + y / 0.30D2) + (0.2D1 / 0.3D1 *&
      !     & (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) &
      !     &** 2 / 0.2D1) / Mref) ** epn * (0.2D1 / 0.3D1 * (x / 0.30D2 - y **&
      !     & 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a * x) * a * sin(b * y) + co&
      !     &s(a * x) * cos(b * y) ** 2 * sin(a * x) * a) / Mref + 0.2D1 / 0.3D&
      !     &1 * (x * y / 0.30D2 + y / 0.30D2) * (cos(a * x) * cos(b * y) * b +&
      !     & cos(a * x) ** 2 * cos(b * y) * sin(b * y) * b) / Mref) * (x / 0.3&
      !     &0D2 + 0.1D1 / 0.30D2)) + pcourr * Mref * cos(a * x) * cos(b * y) *&
      !     & ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (0.2D1 / 0.3D1&
      !     & * cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)&
      !     &) / Mref - 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b * y)) * cos&
      !     &(a * x) * a * cos(b * y) / Mref) + (x * y / 0.30D2 + y / 0.30D2) *&
      !     & (0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.10D2 - sin(a * &
      !     &x) * cos(b * y)) / Mref + 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * sin(a * x) * sin(b * y) * b / Mref)) + 0.3D1 / 0.4D1 *&
      !     & (0.2D1 + sin(a * x) * sin(b * y)) ** 2 / tie * (0.2D1 / 0.3D1 * (&
      !     &0.10D2 - sin(a * x) * cos(b * y)) / Mref - 0.2D1 / 0.3D1 * (0.20D2&
      !     & + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b * y) ** 2 / 0&
      !     &.2D1) / Mref) * sqrt(0.2D1) * sqrt(0.3D1) * ((0.10D2 - sin(a * x) &
      !     &* cos(b * y)) / Mref) ** (-0.3D1 / 0.2D1)
      !
      !       f(:,4) = 0.5D1 / 0.3D1 * cos(a * x) ** 2 * a * sin(b * y) * (0.10D2 &
      !     &- sin(a * x) * cos(b * y)) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0&
      !     &.30D2 + 0.1D1 / 0.15D2) - 0.5D1 / 0.3D1 * (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * cos(a * x) ** 2 * a * cos(b * y) ** 2 * (x / 0.30D2 - &
      !     &y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) - 0.5D1 / 0.3D1 * (0.2D1 + sin(a&
      !     & * x) * sin(b * y)) * (0.10D2 - sin(a * x) * cos(b * y)) * sin(a *&
      !     & x) * a * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.1&
      !     &5D2) + (0.2D1 + sin(a * x) * sin(b * y)) * (0.10D2 - sin(a * x) * &
      !     &cos(b * y)) * cos(a * x) * cos(b * y) / 0.18D2 + 0.5D1 / 0.3D1 * (&
      !     &0.2D1 + sin(a * x) * sin(b * y)) * (0.10D2 - sin(a * x) * cos(b * &
      !     &y)) * cos(a * x) * cos(b * y) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.&
      !     &1D1 / 0.15D2) / x + 0.5D1 / 0.3D1 * sin(a * x) * cos(b * y) ** 2 *&
      !     & b * (0.10D2 - sin(a * x) * cos(b * y)) * cos(a * x) * (x * y / 0.&
      !     &30D2 + y / 0.30D2) + 0.5D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin(b *&
      !     & y)) * sin(a * x) * sin(b * y) * b * cos(a * x) * cos(b * y) * (x &
      !     &* y / 0.30D2 + y / 0.30D2) - 0.5D1 / 0.3D1 * (0.2D1 + sin(a * x) *&
      !     & sin(b * y)) * (0.10D2 - sin(a * x) * cos(b * y)) * cos(a * x) * s&
      !     &in(b * y) * b * (x * y / 0.30D2 + y / 0.30D2) + 0.5D1 / 0.3D1 * (0&
      !     &.2D1 + sin(a * x) * sin(b * y)) * (0.10D2 - sin(a * x) * cos(b * y&
      !     &)) * cos(a * x) * cos(b * y) * (x / 0.30D2 + 0.1D1 / 0.30D2) - csi&
      !     &e * (-sin(a * x) * a ** 2 * sin(b * y) * (0.10D2 - sin(a * x) * co&
      !     &s(b * y)) - 0.2D1 * cos(a * x) ** 2 * a ** 2 * sin(b * y) * cos(b &
      !     &* y) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * a ** 2 * c&
      !     &os(b * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos&
      !     &(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) - (0&
      !     &.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y)) / 0&
      !     &.30D2 - (x * y / 0.30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) *&
      !     & b * (0.10D2 - sin(a * x) * cos(b * y)) + (0.2D1 + sin(a * x) * si&
      !     &n(b * y)) * sin(a * x) * sin(b * y) * b) / 0.30D2 - (x / 0.30D2 - &
      !     &y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * a * sin(b * y) *&
      !     & (0.10D2 - sin(a * x) * cos(b * y)) / 0.30D2 - (0.2D1 + sin(a * x)&
      !     & * sin(b * y)) * cos(a * x) * a * cos(b * y) / 0.30D2 + (x / 0.30D&
      !     &2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (-sin(a * x) * a ** 2 * si&
      !     &n(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) - 0.2D1 * cos(a * x)&
      !     & ** 2 * a ** 2 * sin(b * y) * cos(b * y) + (0.2D1 + sin(a * x) * s&
      !     &in(b * y)) * sin(a * x) * a ** 2 * cos(b * y)) + y * (sin(a * x) *&
      !     & cos(b * y) * b * (0.10D2 - sin(a * x) * cos(b * y)) + (0.2D1 + si&
      !     &n(a * x) * sin(b * y)) * sin(a * x) * sin(b * y) * b) / 0.30D2 + (&
      !     &x * y / 0.30D2 + y / 0.30D2) * (cos(a * x) * a * cos(b * y) * b * &
      !     &(0.10D2 - sin(a * x) * cos(b * y)) - sin(a * x) * cos(b * y) ** 2 &
      !     &* b * cos(a * x) * a + cos(a * x) * a * sin(b * y) ** 2 * sin(a * &
      !     &x) * b + (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * sin(&
      !     &b * y) * b)) + (cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x)&
      !     & * cos(b * y)) - (0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * &
      !     &a * cos(b * y) - (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) *&
      !     & ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) * &
      !     &a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) - (0.2D1 + sin&
      !     &(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y)) + (x * y / 0.&
      !     &30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.10D2 - sin(&
      !     &a * x) * cos(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a *&
      !     & x) * sin(b * y) * b))) / x - sin(a * x) * sin(b * y) * b ** 2 * (&
      !     &0.10D2 - sin(a * x) * cos(b * y)) + 0.2D1 * sin(a * x) ** 2 * cos(&
      !     &b * y) * b ** 2 * sin(b * y) + (0.2D1 + sin(a * x) * sin(b * y)) *&
      !     & sin(a * x) * cos(b * y) * b ** 2 - (x / 0.30D2 + 0.1D1 / 0.30D2) &
      !     &* ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a * x) *&
      !     & a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) - (0.2D1 + si&
      !     &n(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y)) + (x * y / 0&
      !     &.30D2 + y / 0.30D2) * (sin(a * x) * cos(b * y) * b * (0.10D2 - sin&
      !     &(a * x) * cos(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a &
      !     &* x) * sin(b * y) * b)) - (x * y / 0.30D2 + y / 0.30D2) * (-y * (c&
      !     &os(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x) * cos(b * y)) - &
      !     &(0.2D1 + sin(a * x) * sin(b * y)) * cos(a * x) * a * cos(b * y)) /&
      !     & 0.15D2 + (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * (cos(a&
      !     & * x) * a * cos(b * y) * b * (0.10D2 - sin(a * x) * cos(b * y)) - &
      !     &sin(a * x) * cos(b * y) ** 2 * b * cos(a * x) * a + cos(a * x) * a&
      !     & * sin(b * y) ** 2 * sin(a * x) * b + (0.2D1 + sin(a * x) * sin(b &
      !     &* y)) * cos(a * x) * a * sin(b * y) * b) + (x / 0.30D2 + 0.1D1 / 0&
      !     &.30D2) * (sin(a * x) * cos(b * y) * b * (0.10D2 - sin(a * x) * cos&
      !     &(b * y)) + (0.2D1 + sin(a * x) * sin(b * y)) * sin(a * x) * sin(b &
      !     &* y) * b) + (x * y / 0.30D2 + y / 0.30D2) * (-sin(a * x) * sin(b *&
      !     & y) * b ** 2 * (0.10D2 - sin(a * x) * cos(b * y)) + 0.2D1 * sin(a &
      !     &* x) ** 2 * cos(b * y) * b ** 2 * sin(b * y) + (0.2D1 + sin(a * x)&
      !     & * sin(b * y)) * sin(a * x) * cos(b * y) * b ** 2))) - kpare * (-(&
      !     &0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** epn &
      !     &* epn * cos(a * x) * a * cos(b * y) / (0.10D2 - sin(a * x) * cos(b&
      !     & * y)) * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 /&
      !     & 0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D1 * (x&
      !     & * y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) * b / Mref) &
      !     &* (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) + (0.2D1 / 0.3D1&
      !     & * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** epn * (-cos(a * x&
      !     &) * a * cos(b * y) / Mref / 0.45D2 + 0.2D1 / 0.3D1 * (x / 0.30D2 -&
      !     & y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * sin(a * x) * a ** 2 * cos(b *&
      !     & y) / Mref + y * sin(a * x) * sin(b * y) * b / Mref / 0.45D2 + 0.2&
      !     &D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) * cos(a * x) * a * sin(&
      !     &b * y) * b / Mref) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D&
      !     &2) + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) *&
      !     &* epn * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / &
      !     &0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D1 * (x &
      !     &* y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) * b / Mref) /&
      !     & 0.30D2 + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mr&
      !     &ef) ** epn * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1&
      !     &D1 / 0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D1 &
      !     &* (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) * b / Mr&
      !     &ef) * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) / x + (0.2D1&
      !     & / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref) ** epn * epn&
      !     & * sin(a * x) * sin(b * y) * b / (0.10D2 - sin(a * x) * cos(b * y)&
      !     &) * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15&
      !     &D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D1 * (x * y &
      !     &/ 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) * b / Mref) * (x &
      !     &* y / 0.30D2 + y / 0.30D2) + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x)&
      !     & * cos(b * y)) / Mref) ** epn * (0.2D1 / 0.45D2 * y * cos(a * x) *&
      !     & a * cos(b * y) / Mref + 0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.&
      !     &30D2 + 0.1D1 / 0.15D2) * cos(a * x) * a * sin(b * y) * b / Mref + &
      !     &0.2D1 / 0.3D1 * (x / 0.30D2 + 0.1D1 / 0.30D2) * sin(a * x) * sin(b&
      !     & * y) * b / Mref + 0.2D1 / 0.3D1 * (x * y / 0.30D2 + y / 0.30D2) *&
      !     & sin(a * x) * cos(b * y) * b ** 2 / Mref) * (x * y / 0.30D2 + y / &
      !     &0.30D2) + (0.2D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mr&
      !     &ef) ** epn * (-0.2D1 / 0.3D1 * (x / 0.30D2 - y ** 2 / 0.30D2 + 0.1&
      !     &D1 / 0.15D2) * cos(a * x) * a * cos(b * y) / Mref + 0.2D1 / 0.3D1 &
      !     &* (x * y / 0.30D2 + y / 0.30D2) * sin(a * x) * sin(b * y) * b / Mr&
      !     &ef) * (x / 0.30D2 + 0.1D1 / 0.30D2)) - pcourr * Mref * cos(a * x) &
      !     &* cos(b * y) * ((x / 0.30D2 - y ** 2 / 0.30D2 + 0.1D1 / 0.15D2) * &
      !     &(0.2D1 / 0.3D1 * cos(a * x) * a * sin(b * y) * (0.10D2 - sin(a * x&
      !     &) * cos(b * y)) / Mref - 0.2D1 / 0.3D1 * (0.2D1 + sin(a * x) * sin&
      !     &(b * y)) * cos(a * x) * a * cos(b * y) / Mref) + (x * y / 0.30D2 +&
      !     & y / 0.30D2) * (0.2D1 / 0.3D1 * sin(a * x) * cos(b * y) * b * (0.1&
      !     &0D2 - sin(a * x) * cos(b * y)) / Mref + 0.2D1 / 0.3D1 * (0.2D1 + s&
      !     &in(a * x) * sin(b * y)) * sin(a * x) * sin(b * y) * b / Mref)) - 0&
      !     &.3D1 / 0.4D1 * (0.2D1 + sin(a * x) * sin(b * y)) ** 2 / tie * (0.2&
      !     &D1 / 0.3D1 * (0.10D2 - sin(a * x) * cos(b * y)) / Mref - 0.2D1 / 0&
      !     &.3D1 * (0.20D2 + cos(a * x) * sin(b * y) - cos(a * x) ** 2 * cos(b&
      !     & * y) ** 2 / 0.2D1) / Mref) * sqrt(0.2D1) * sqrt(0.3D1) * ((0.10D2&
      !     & - sin(a * x) * cos(b * y)) / Mref) ** (-0.3D1 / 0.2D1)

    CASE (5:6)
      ! Do nothing
    CASE (50:)
      ! Do nothing

    CASE DEFAULT
      WRITE (6, *) "Error! Test case not valid"
      STOP

    END SELECT
  END SUBROUTINE body_force
#endif
END MODULE analytical
