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
  SUBROUTINE analytical_solution(x,y,t,u)
  real*8, dimension(:), intent(IN)        :: x,y,t
  real*8, dimension(:,:), intent(OUT)     :: u
  real*8, dimension(size(u,1),phys%npv)  :: up
  integer:: i,j,ind,N1D,N2D
  real*8 :: a,b,r,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
  real*8 :: aux
    
  u = 0.
  up = 0.
  a = 2*pi
  b = 2*pi
  N2D = size(x,1)
  N1D = size(t,1)
  
  xmax = Mesh%xmax
  xmin = Mesh%xmin
  ymax = Mesh%ymax
  ymin = Mesh%ymin
  xm = 0.5*(xmax+xmin)
  ym = 0.5*(ymax+ymin)  
  
			DO i=1,N2d
						 DO j=1,N1d
						     xx = x(i)
						     yy = y(i)
						     tt = t(j)
						     r = sqrt((xx-xm)**2+(yy-ym)**2)
						     ind = (j-1)*N2d + i
						     SELECT CASE(switch%testcase)
											CASE(1)
															IF (switch%axisym) THEN
																 WRITE(6,*) "This is NOT an axisymmetric test case!"
																 stop
														END IF
               ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
															up(ind,1) = 2+sin(a* xx)*sin(a* yy)
															up(ind,2) = cos(a* xx)*cos(a* yy)
															up(ind,3) = 20+cos(a*xx)*sin(a*yy)
															up(ind,4) = 10-sin(a*xx)*cos(a*yy)												
			
											CASE(2)
											! Axisimmetric case with div(b)~=0
														IF (.not.switch%axisym) THEN
																 WRITE(6,*) "This is an axisymmetric test case!"
																 stop
														END IF 
															up(ind,1) = 2+sin(a* xx)*sin(a* yy)
															up(ind,2) = cos(a* xx)*cos(a* yy)
															up(ind,3) = 20+cos(a*xx)*sin(a*yy)
															up(ind,4) = 10-sin(a*xx)*cos(a*yy)			
															
											CASE(3)
											! Axisimmetric case with div(b)~=0
														IF (.not.switch%axisym) THEN
																 WRITE(6,*) "This is an axisymmetric test case!"
																 stop
														END IF 
															up(ind,1) = 2+sin(a* tt)
															up(ind,2) = cos(a* tt)
															up(ind,3) = 20+cos(a*tt)
															up(ind,4) = 10-sin(a*tt)
																																										
											CASE(50:63)
														 up(ind,1) = 1.
														 up(ind,3) = 18.
														 up(ind,4) = 18.
           CASE(64)
														 up(ind,1) = 1.
														 up(ind,3) = 18.
														 up(ind,4) = 18.
               r = sqrt((xx-geom%R0)**2+yy**2)
!               if (abs(r-1).lt.1e-6) then
!               ! outer boundary
!               endif 
											CASE(65)
														 up(ind,1) = 1.
														 up(ind,2) = 0.
														 r = 	sqrt ( (xx*phys%lscale-geom%R0)**2 + (yy*phys%lscale-0.75)**2 )
											    IF (r.le. 0.05) THEN
											       up(ind,2) = 1.
											    END IF     
											CASE DEFAULT
														 WRITE(6,*) "Error! Test case not valid"
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
  real*8, dimension(:), intent(IN)        :: x,y,t
  real*8, dimension(:,:), intent(IN)      :: u
  real*8, dimension(:,:), intent(OUT)     :: ux,uy,ut
  real*8, dimension(size(u,1),size(u,2))  :: upx,upy,upt
  real*8, dimension(size(u,1),phys%npv)   :: up
  integer:: i,j,ind,N1D,N2D
  real*8 :: a,b,r,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
  real*8 :: aux
    
  N2D = size(x,1)
  N1D = size(t,1)
  
  xmax = Mesh%xmax
  xmin = Mesh%xmin
  ymax = Mesh%ymax
  ymin = Mesh%ymin
  xm = 0.5*(xmax+xmin)
  ym = 0.5*(ymax+ymin)  
  upx = 0.
  upy = 0.
  upt = 0.
  CALL cons2phys(u,up)
  a = 2*pi
    
			DO i=1,N2d
						 DO j=1,N1d
						     xx = x(i)
						     yy = y(i)
						     tt = t(j)
						     r = sqrt((xx-xm)**2+(yy-ym)**2)
						     ind = (j-1)*N2d + i
											SELECT CASE(switch%testcase)
													CASE(1)
																	IF (switch%axisym) THEN
																			WRITE(6,*) "This is NOT an axisymmetric test case!"
																			stop
																END IF
													! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
														 upx(ind,1) =  a*cos(a*xx)*sin(a*yy)
														 upx(ind,2) = -a*sin(a*xx)*cos(a*yy)
														 upx(ind,3) = -a*sin(a*xx)*sin(a*yy)
														 upx(ind,4) = -a*cos(a*xx)*cos(a*yy)

														 upy(ind,1) =  a*sin(a*xx)*cos(a*yy)
														 upy(ind,2) = -a*cos(a*xx)*sin(a*yy)
														 upy(ind,3) =  a*cos(a*xx)*cos(a*yy)
														 upy(ind,4) =  a*sin(a*xx)*sin(a*yy)
												
													CASE(2)
													! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
																IF (.not.switch%axisym) THEN
																			WRITE(6,*) "This is an axisymmetric test case!"
																			stop
																END IF 
														 upx(ind,1) =  a*cos(a*xx)*sin(a*yy)
														 upx(ind,2) = -a*sin(a*xx)*cos(a*yy)
														 upx(ind,3) = -a*sin(a*xx)*sin(a*yy)
														 upx(ind,4) = -a*cos(a*xx)*cos(a*yy)

														 upy(ind,1) =  a*sin(a*xx)*cos(a*yy)
														 upy(ind,2) = -a*cos(a*xx)*sin(a*yy)
														 upy(ind,3) =  a*cos(a*xx)*cos(a*yy)
														 upy(ind,4) =  a*sin(a*xx)*sin(a*yy)
													CASE(3)
													! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
																IF (.not.switch%axisym) THEN
																			WRITE(6,*) "This is an axisymmetric test case!"
																			stop
																END IF 														 
															upt(ind,1) = +a*cos(a* tt)
															upt(ind,2) = -a*sin(a* tt)
															upt(ind,3) = -a*sin(a*tt)
															upt(ind,4) = -a*cos(a*tt)														 
													CASE(5)
														   ! Do nothing
													CASE(6)
														   ! Do nothing
													CASE(50:)
																	! Do nothing   
													CASE DEFAULT
																	WRITE(6,*) "Error! Test case not valid"
																	STOP
						     END SELECT
						 END DO
			END DO
			
			
				! Convert physical variables to conservative variables
				ux(:,1) = upx(:,1) 
				uy(:,1) = upy(:,1)
				ut(:,1) = upt(:,1)
				ux(:,2) = ( upx(:,1)*up(:,2)+up(:,1)*upx(:,2) )
				uy(:,2) = ( upy(:,1)*up(:,2)+up(:,1)*upy(:,2) )
				ut(:,2) = ( upt(:,1)*up(:,2)+up(:,1)*upt(:,2) )				
				ux(:,3) = ( upx(:,1)*up(:,3)+up(:,1)*upx(:,3) )
				uy(:,3) = ( upy(:,1)*up(:,3)+up(:,1)*upy(:,3) )
				ut(:,3) = ( upt(:,1)*up(:,3)+up(:,1)*upt(:,3) )				
				ux(:,4) = ( upx(:,1)*up(:,4)+up(:,1)*upx(:,4) )
				uy(:,4) = ( upy(:,1)*up(:,4)+up(:,1)*upy(:,4) )
				ut(:,4) = ( upt(:,1)*up(:,4)+up(:,1)*upt(:,4) )				

  END SUBROUTINE analytical_gradient
  
  
		!*****************************************
		! Body forces
		!**************************************** 
  SUBROUTINE body_force(x,y,t,f)
  real*8, dimension(:), intent(IN) :: x,y,t
  real*8, dimension(:,:), intent(OUT) :: f
  integer                            ::  n
  real*8  :: a,b,D,mu,csii,csie,kpari,kpare,Mref,epn,tie,pcourr
  integer:: i,j,ind,N1D,N2D
  real*8 :: r,xx,yy,tt,xmax,xmin,ymax,ymin,xm,ym
  real*8 :: aux
  real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
  real*8 :: t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
  real*8 :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t30
  real*8 :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t40
  real*8 :: t41,t42,t43,t44,t45,t46,t47,t48,t49,t50
  real*8 :: t51,t52,t53,t54,t55,t56,t57,t58,t59,t60
  real*8 :: t61,t62,t63,t64,t65,t66,t67,t68,t69,t70
  real*8 :: cx,cy,sx,sy,cx2,sx2,cy2,sy2,r2,cc,ss,ct,st
 
  N2D = size(x,1)
  N1D = size(t,1)
  
  xmax = Mesh%xmax
  xmin = Mesh%xmin
  ymax = Mesh%ymax
  ymin = Mesh%ymin
  xm = 0.5*(xmax+xmin)
  ym = 0.5*(ymax+ymin)  
  
!  xm = -0.5
!  ym = -0.5

  a = 2*pi
    f = 0.    
			DO i=1,N2d
						 DO j=1,N1d
						     xx = x(i)
						     yy = y(i)
						     tt = t(j)
						     r = sqrt((xx-xm)**2+(yy-ym)**2)
						     ind = (j-1)*N2d + i
											SELECT CASE(switch%testcase)  
  
				       CASE(1)
				       ! Cartesian case, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
				       a = 2*pi
				       b = 2*pi
				       D     = phys%diff_n
				       mu    = phys%diff_u
				       csii  = phys%diff_e
				       csie  = phys%diff_ee
				       kpari = phys%diff_pari
				       kpare = phys%diff_pare
				       Mref  = phys%Mref
				       epn   = phys%epn
				       tie   = phys%tie
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
	          
				       
           f(ind,1) = (D*a**2*sx*sy*2.0D0+D*a**2*sx*sy*xm**2+D*a**2*sx*sy*xx**2+D*a&
                       &**2*sx*sy*ym**2+D*a**2*sx*sy*yy**2+D*a*cx*sy*xm-D*a*cx*sy*xx+D*a*c&
                       &y*sx*ym-D*a*cy*sx*yy-a*cx*r*sy*xm*2.0D0+a*cx*r*sy*xx*2.0D0+a*cy*r*&
                       &sx*ym*2.0D0-a*cy*r*sx*yy*2.0D0-D*a**2*cx*cy*xm*ym*2.0D0+D*a**2*cx*&
                       &cy*xx*ym*2.0D0+D*a**2*cx*cy*xm*yy*2.0D0-D*a**2*cx*cy*xx*yy*2.0D0-D&
                       &*a**2*sx*sy*xm*xx*2.0D0-D*a**2*sx*sy*ym*yy*2.0D0+a*cx*cy2*r*sx*xm-&
                       &a*cx*cy2*r*sx*xx-a*cx2*cy*r*sy*ym+a*cx2*cy*r*sy*yy-a*cx*r*sx*sy2*x&
                       &m+a*cx*r*sx*sy2*xx+a*cy*r*sx2*sy*ym-a*cy*r*sx2*sy*yy)/((xm-xx)**2+&
                       &(ym-yy)**2+1.0D0)
	          
           f(ind,2) = (a*(cx*cy*r*xm*4.0D0-cx*cy*r*xx*4.0D0+cx*cy*r*ym*4.0D0-cx*cy*&
                       &r*yy*4.0D0-cy*mu*sx*xm*6.0D0+cy*mu*sx*xx*6.0D0-cx*mu*sy*ym*6.0D0+c&
                       &x*mu*sy*yy*6.0D0+cy*r*sx*xm*6.0D+1-cy2*r*sx2*xm*2.0D0-cy*r*sx*xx*6&
                       &.0D+1+cy2*r*sx2*xx*2.0D0-cx*r*sy*ym*6.0D+1-cx2*r*sy2*ym*2.0D0+cx*r&
                       &*sy*yy*6.0D+1+cx2*r*sy2*yy*2.0D0+r*sx*sy*xm*4.0D0+r*sx2*sy2*xm*2.0&
                       &D0-r*sx*sy*xx*4.0D0-r*sx2*sy2*xx*2.0D0+r*sx*sy*ym*4.0D0+r*sx2*sy2*&
                       &ym*2.0D0-r*sx*sy*yy*4.0D0-r*sx2*sy2*yy*2.0D0+a*cx*cy*mu*1.2D+1+cx2&
                       &*cy*mu*sy*xm*3.0D0-cx2*cy*mu*sy*xx*3.0D0+cx*cy2*mu*sx*ym*3.0D0-cx*&
                       &cy2*mu*sx*yy*3.0D0-cx2*cy*r*sy*xm*8.0D0+cx2*cy*r*sy*xx*8.0D0+cx*cy&
                       &2*r*sx*ym*8.0D0-cx*cy2*r*sx*yy*8.0D0-cy*mu*sx2*sy*xm*3.0D0+cy*mu*s&
                       &x2*sy*xx*3.0D0-cx*mu*sx*sy2*ym*3.0D0+cx*mu*sx*sy2*yy*3.0D0+a*cx*cy&
                       &*mu*xm**2*6.0D0+a*cx*cy*mu*xx**2*6.0D0+a*cx*cy*mu*ym**2*6.0D0+a*cx&
                       &*cy*mu*yy**2*6.0D0+cx2*cy**3*r*sx*xm*2.0D0-cx2*cy**3*r*sx*xx*2.0D0&
                       &-cx**3*cy2*r*sy*ym*2.0D0+cx**3*cy2*r*sy*yy*2.0D0-a*mu*sx*sy*xm*ym*&
                       &1.2D+1-a*mu*sx2*sy2*xm*ym*6.0D0+a*mu*sx*sy*xx*ym*1.2D+1+a*mu*sx2*s&
                       &y2*xx*ym*6.0D0+a*mu*sx*sy*xm*yy*1.2D+1+a*mu*sx2*sy2*xm*yy*6.0D0-a*&
                       &mu*sx*sy*xx*yy*1.2D+1-a*mu*sx2*sy2*xx*yy*6.0D0+a*cx*cy*mu*sx*sy*2.&
                       &4D+1-a*cx*cy*mu*xm*xx*1.2D+1-a*cx2*cy2*mu*xm*ym*6.0D0+a*cx2*cy2*mu&
                       &*xx*ym*6.0D0+a*cx2*cy2*mu*xm*yy*6.0D0-a*cx2*cy2*mu*xx*yy*6.0D0-a*c&
                       &x*cy*mu*ym*yy*1.2D+1+a*cx2*mu*sy2*xm*ym*6.0D0+a*cy2*mu*sx2*xm*ym*6&
                       &.0D0-a*cx2*mu*sy2*xx*ym*6.0D0-a*cy2*mu*sx2*xx*ym*6.0D0-a*cx2*mu*sy&
                       &2*xm*yy*6.0D0-a*cy2*mu*sx2*xm*yy*6.0D0+a*cx2*mu*sy2*xx*yy*6.0D0+a*&
                       &cy2*mu*sx2*xx*yy*6.0D0+cx*cy*r*sx*sy*xm*4.0D0-cx2*cy*r*sx*sy2*xm*4&
                       &.0D0-cx*cy*r*sx*sy*xx*4.0D0+cx2*cy*r*sx*sy2*xx*4.0D0+cx*cy*r*sx*sy&
                       &*ym*4.0D0+cx*cy2*r*sx2*sy*ym*4.0D0-cx*cy*r*sx*sy*yy*4.0D0-cx*cy2*r&
                       &*sx2*sy*yy*4.0D0+a*cx*cy*mu*sx*sy*xm**2*1.2D+1+a*cx*cy*mu*sx*sy*xx&
                       &**2*1.2D+1+a*cx*cy*mu*sx*sy*ym**2*1.2D+1+a*cx*cy*mu*sx*sy*yy**2*1.&
                       &2D+1-a*cx*cy*mu*sx*sy*xm*xx*2.4D+1-a*cx*cy*mu*sx*sy*ym*yy*2.4D+1))&
                       &/(xm*xx*(-6.0D0)-ym*yy*6.0D0+xm**2*3.0D0+xx**2*3.0D0+ym**2*3.0D0+y&
                       &y**2*3.0D0+3.0D0)

           f(ind,3) = (csii*(-a*xm+a*xx+a**2*sin(a*xx*2.0D0)*3.0D0+a**2*sx*sy*4.0D+&
                       &1+a*cx2*xm*2.0D0+a*cy2*xm-a*cx2*xx*2.0D0-a*cy2*xx+a**2*xm**2*sin(a&
                       &*xx*2.0D0)*2.0D0+a**2*xx**2*sin(a*xx*2.0D0)*2.0D0+a**2*ym**2*sin(a&
                       &*xx*2.0D0)+a**2*yy**2*sin(a*xx*2.0D0)+a**2*cx*sy*4.0D0+a*cx*cy*ym*&
                       &2.0D0-a*cx*cy*yy*2.0D0+a*cx*sy*xm*2.0D+1-a*cx*sy*xx*2.0D+1+a*cy*sx&
                       &*ym*2.0D+1-a*cy*sx*yy*2.0D+1-a*sx*sy*xm*2.0D0+a*sx*sy*xx*2.0D0-a**&
                       &2*cx*cy2*sx*8.0D0-a**2*xm*xx*sin(a*xx*2.0D0)*4.0D0+a**2*xm*ym*sin(&
                       &a*yy*2.0D0)*2.0D0-a**2*xx*ym*sin(a*yy*2.0D0)*2.0D0-a**2*ym*yy*sin(&
                       &a*xx*2.0D0)*2.0D0-a**2*xm*yy*sin(a*yy*2.0D0)*2.0D0+a**2*xx*yy*sin(&
                       &a*yy*2.0D0)*2.0D0+a**2*cx*sy*xm**2*2.0D0+a**2*cx*sy*xx**2*2.0D0+a*&
                       &*2*cx*sy*ym**2*2.0D0+a**2*cx*sy*yy**2*2.0D0+a**2*sx*sy*xm**2*2.0D+&
                       &1+a**2*sx*sy*xx**2*2.0D+1+a**2*sx*sy*ym**2*2.0D+1+a**2*sx*sy*yy**2&
                       &*2.0D+1-a*cx2*cy2*xm*2.0D0+a*cx2*cy2*xx*2.0D0-a**2*cx*cy2*sx*xm**2&
                       &*4.0D0-a**2*cx*cy2*sx*xx**2*4.0D0-a**2*cx*cy2*sx*ym**2*4.0D0-a**2*&
                       &cx*cy2*sx*yy**2*4.0D0-a**2*cx*cy*xm*ym*4.0D+1+a**2*cx*cy*xx*ym*4.0&
                       &D+1+a**2*cx*cy*xm*yy*4.0D+1-a**2*cx*cy*xx*yy*4.0D+1-a**2*cx*sy*xm*&
                       &xx*4.0D0+a**2*cy*sx*xm*ym*4.0D0-a**2*cy*sx*xx*ym*4.0D0-a**2*cy*sx*&
                       &xm*yy*4.0D0+a**2*cy*sx*xx*yy*4.0D0-a**2*cx*sy*ym*yy*4.0D0-a**2*sx*&
                       &sy*xm*xx*4.0D+1-a**2*sx*sy*ym*yy*4.0D+1+a**2*cx*cy2*sx*xm*xx*8.0D0&
                       &-a**2*cx2*cy*sy*xm*ym*8.0D0+a**2*cx2*cy*sy*xx*ym*8.0D0+a**2*cx2*cy&
                       &*sy*xm*yy*8.0D0-a**2*cx2*cy*sy*xx*yy*8.0D0+a**2*cx*cy2*sx*ym*yy*8.&
                       &0D0+a*cx*cy*sx*sy*ym*2.0D0-a*cx*cy*sx*sy*yy*2.0D0))/(xm*xx*(-2.0D0&
                       &)-ym*yy*2.0D0+xm**2+xx**2+ym**2+yy**2+1.0D0)-(1.0D0/Mref**2*(3.0D0&
                       &**(-epn+1.0D0)*a**2*cx2*cy2*epn*kpari*xm**2*((-cx2*cy2+cx*sy*2.0D0&
                       &+4.0D+1)/Mref)**(epn-1.0D0)*4.0D0+3.0D0**(-epn+1.0D0)*a**2*cx2*cy2&
                       &*epn*kpari*xx**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)&
                       &*4.0D0+3.0D0**(-epn+1.0D0)*a**2*epn*kpari*sx2*sy2*ym**2*((-cx2*cy2&
                       &+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*4.0D0+3.0D0**(-epn+1.0D0)*&
                       &a**2*epn*kpari*sx2*sy2*yy**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)*&
                       &*(epn-1.0D0)*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy2*kpari*xm*&
                       &*2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.&
                       &0D0)*Mref*a**2*cx2*cy2*kpari*xx**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/&
                       &Mref)**epn*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy2*kpari*ym**2&
                       &*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.0D&
                       &0)*Mref*a**2*cx2*cy2*kpari*yy**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mr&
                       &ef)**epn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpari*sy*xm**2*((-&
                       &cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0-3.0D0**(-epn+1.0D0)*M&
                       &ref*a**2*cx*kpari*sy*xx**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**e&
                       &pn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpari*sy*ym**2*((-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*&
                       &*2*cx*kpari*sy*yy**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0&
                       &D0+3.0D0**(-epn+1.0D0)*Mref*a*cx*cy*kpari*ym*((-cx2*cy2+cx*sy*2.0D&
                       &0+4.0D+1)/Mref)**epn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*cx*cy*kpari*&
                       &yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.&
                       &0D0)*Mref*a**2*cx2*kpari*xm**2*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+&
                       &4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*kpari*x&
                       &x**2*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3&
                       &.0D0**(-epn+1.0D0)*Mref*a**2*cy2*kpari*ym**2*(cx2-1.0D0)*((-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*&
                       &*2*cy2*kpari*yy**2*(cx2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref&
                       &)**epn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpari*sx*sy*xm*((-cx2*cy2+&
                       &cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kp&
                       &ari*sx*sy*xx*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0-3.0D0&
                       &**(-epn+1.0D0)*Mref*a**2*cx2*kpari*xm*xx*(cy2-1.0D0)*((-cx2*cy2+cx&
                       &*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*c&
                       &y2*kpari*ym*yy*(cx2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**e&
                       &pn*4.0D0-3.0D0**(-epn+1.0D0)*a**2*cx2*cy2*epn*kpari*xm*xx*((-cx2*c&
                       &y2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0-3.0D0**(-epn+1.0D0&
                       &)*a**2*epn*kpari*sx2*sy2*ym*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref&
                       &)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cx*cy2*kpari*sx*xm&
                       &*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0-3.0D0**(-epn+1.0D&
                       &0)*Mref*a*cx*cy2*kpari*sx*xx*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)*&
                       &*epn*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cx2*cy*kpari*sy*ym*((-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*&
                       &cx2*cy*kpari*sy*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*2.0D0&
                       &-3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy2*kpari*xm*xx*((-cx2*cy2+cx*s&
                       &y*2.0D0+4.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cx2&
                       &*cy2*kpari*ym*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0+3&
                       &.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpari*sy*xm*xx*((-cx2*cy2+cx*sy*2.&
                       &0D0+4.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpar&
                       &i*sx*xm*ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0-3.0D0**&
                       &(-epn+1.0D0)*Mref*a**2*cy*kpari*sx*xx*ym*((-cx2*cy2+cx*sy*2.0D0+4.&
                       &0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpari*sx*x&
                       &m*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+&
                       &1.0D0)*Mref*a**2*cy*kpari*sx*xx*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/&
                       &Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpari*sy*ym*yy*(&
                       &(-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)&
                       &*a**2*cx**3*cy2*epn*kpari*sy*xm**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/&
                       &Mref)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*a**2*cx**3*cy2*epn*kp&
                       &ari*sy*xx**2*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0&
                       &D0-3.0D0**(-epn+1.0D0)*a**2*cx**3*cy2*epn*kpari*sy*xm*xx*((-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*1.6D+1-3.0D0**(-epn+1.0D0&
                       &)*a**2*cx2*cy**3*epn*kpari*sx*xm*ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)&
                       &/Mref)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*a**2*cx2*cy**3*epn*k&
                       &pari*sx*xx*ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.&
                       &0D0+3.0D0**(-epn+1.0D0)*a**2*cx2*cy**3*epn*kpari*sx*xm*yy*((-cx2*c&
                       &y2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0-3.0D0**(-epn+1.0D0&
                       &)*a**2*cx2*cy**3*epn*kpari*sx*xx*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)&
                       &/Mref)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*a**2*cx*cy*epn*kpari&
                       &*sx*sy*xm*ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0&
                       &D0-3.0D0**(-epn+1.0D0)*a**2*cx*cy*epn*kpari*sx*sy*xx*ym*((-cx2*cy2&
                       &+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0-3.0D0**(-epn+1.0D0)*&
                       &a**2*cx*cy*epn*kpari*sx*sy*xm*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mr&
                       &ef)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*a**2*cx*cy*epn*kpari*sx&
                       &*sy*xx*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0-&
                       &(3.0D0**(-epn+1.0D0)*Mref*a**2*cx**4*cy2*epn*kpari*xm**2*(cy2-1.0D&
                       &0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0)/(-cx2*cy2+cx*s&
                       &y*2.0D0+4.0D+1)-(3.0D0**(-epn+1.0D0)*Mref*a**2*cx**4*cy2*epn*kpari&
                       &*xx**2*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*4.0D0&
                       &)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)-(3.0D0**(-epn+1.0D0)*Mref*a**2*cx2&
                       &*cy**4*epn*kpari*ym**2*(cx2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/&
                       &Mref)**epn*4.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)-(3.0D0**(-epn+1.0D&
                       &0)*Mref*a**2*cx2*cy**4*epn*kpari*yy**2*(cx2-1.0D0)*((-cx2*cy2+cx*s&
                       &y*2.0D0+4.0D+1)/Mref)**epn*4.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)+3.&
                       &0D0**(-epn+1.0D0)*Mref*a**2*cx*cy*kpari*sx*sy*xm*ym*((-cx2*cy2+cx*&
                       &sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cx&
                       &*cy*kpari*sx*sy*xx*ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.&
                       &0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cx*cy*kpari*sx*sy*xm*yy*((-cx2*c&
                       &y2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0+3.0D0**(-epn+1.0D0)*Mref*a&
                       &**2*cx*cy*kpari*sx*sy*xx*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**&
                       &epn*8.0D0-3.0D0**(-epn+1.0D0)*a**2*cx**3*cy**3*epn*kpari*sx*sy*xm*&
                       &ym*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0+3.0D0**&
                       &(-epn+1.0D0)*a**2*cx**3*cy**3*epn*kpari*sx*sy*xx*ym*((-cx2*cy2+cx*&
                       &sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*8.0D0+3.0D0**(-epn+1.0D0)*a**2&
                       &*cx**3*cy**3*epn*kpari*sx*sy*xm*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/&
                       &Mref)**(epn-1.0D0)*8.0D0-3.0D0**(-epn+1.0D0)*a**2*cx**3*cy**3*epn*&
                       &kpari*sx*sy*xx*yy*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0&
                       &)*8.0D0+(3.0D0**(-epn+1.0D0)*Mref*a**2*cx*cy2*epn*kpari*sy*ym**2*(&
                       &cx2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0)/(-cx2*&
                       &cy2+cx*sy*2.0D0+4.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*cx*cy2*epn*&
                       &kpari*sy*yy**2*(cx2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**e&
                       &pn*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*&
                       &a**2*cx**4*cy2*epn*kpari*xm*xx*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+&
                       &4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)+(3.0D0**(-&
                       &epn+1.0D0)*Mref*a**2*cx2*cy**4*epn*kpari*ym*yy*(cx2-1.0D0)*((-cx2*&
                       &cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.&
                       &0D+1)-(3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy*epn*kpari*sx*xm*ym*(cy&
                       &2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy*epn*kp&
                       &ari*sx*xx*ym*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn&
                       &*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a*&
                       &*2*cx2*cy*epn*kpari*sx*xm*yy*(cy2-1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.&
                       &0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D+1)-(3.0D0**(-ep&
                       &n+1.0D0)*Mref*a**2*cx2*cy*epn*kpari*sx*xx*yy*(cy2-1.0D0)*((-cx2*cy&
                       &2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*8.0D0)/(-cx2*cy2+cx*sy*2.0D0+4.0D&
                       &+1)-(3.0D0**(-epn+1.0D0)*Mref*a**2*cx*cy2*epn*kpari*sy*ym*yy*(cx2-&
                       &1.0D0)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*1.6D+1)/(-cx2*cy2&
                       &+cx*sy*2.0D0+4.0D+1)))/((xm-xx)**2*9.0D0+(ym-yy)**2*9.0D0+9.0D0)+M&
                       &ref*cx*cy*pcourr*((a*sx*(xm-xx)*(cy*1.0D+1+sx+sy*2.0D0-cy2*sx*2.0D&
                       &0)*(2.0D0/3.0D0))/(Mref*r)+(a*cx*(ym-yy)*(cy-sy*5.0D0+cy*sx*sy)*(4&
                       &.0D0/3.0D0))/(Mref*r))-(a*cx*cy*(ym-yy)*(cx*sy*1.0D+2+cx2*sy2*5.0D&
                       &0-sx*sy*1.0D+1-sx2*sy2*5.0D0+cx*cy2*sx*4.0D0-cx**3*cy2*sy+cx*cy2*s&
                       &x2*sy*2.0D0))/(r*3.0D0)+(sqrt(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0-2.0D&
                       &+1)/Mref)*(sx*sy+2.0D0)**2*(-cx2*cy2+cx*sy*2.0D0+cy*sx*2.0D0+2.0D+&
                       &1))/(tie*(cy*sx-1.0D+1)*2.0D0)+(a*cx*cy2*(xm-xx)*(cx*1.0D+1+sx*1.0&
                       &D+2+cx2*sx*2.0D0+cx2*sy*4.0D0-cx2*cy2*sx*3.0D0+cx*sx*sy*1.0D+1))/(&
                       &r*3.0D0)-(cx*cy*(xm*2.0D0-xx*2.0D0)*(sx*sy+2.0D0)*(ym-yy)*(-cx2*cy&
                       &2+cx*sy*5.0D0+1.0D+2))/(r2*6.0D0)+(cx*cy*(ym*2.0D0-yy*2.0D0)*(sx*s&
                       &y+2.0D0)*(xm-xx)*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/(r2*6.0D0)-(a*cx*s&
                       &y*(sx*sy+2.0D0)*(xm-xx)*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/(r*3.0D0)+(&
                       &a*cy*sx*(sx*sy+2.0D0)*(ym-yy)*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/(r*3.&
                       &0D0)
                      
          f(ind,4) =  (csie*(a*ym-a*yy-a**2*sin(a*yy*2.0D0)*3.0D0+a**2*sx*sy*2.0D+1&
                       &-a*cx2*ym-a*cy2*ym*2.0D0+a*cx2*yy+a*cy2*yy*2.0D0-a**2*xm**2*sin(a*&
                       &yy*2.0D0)-a**2*xx**2*sin(a*yy*2.0D0)-a**2*ym**2*sin(a*yy*2.0D0)*2.&
                       &0D0-a**2*yy**2*sin(a*yy*2.0D0)*2.0D0-a**2*cy*sx*4.0D0+a*cx2*cy2*ym&
                       &*2.0D0-a*cx2*cy2*yy*2.0D0+a*cx*sy*xm*1.0D+1-a*cx*sy*xx*1.0D+1+a*cy&
                       &*sx*ym*1.0D+1-a*cy*sx*yy*1.0D+1+a*sx*sy*ym*2.0D0-a*sx*sy*yy*2.0D0+&
                       &a**2*cx2*cy*sy*8.0D0-a**2*xm*ym*sin(a*xx*2.0D0)*2.0D0+a**2*xx*ym*s&
                       &in(a*xx*2.0D0)*2.0D0+a**2*xm*xx*sin(a*yy*2.0D0)*2.0D0+a**2*xm*yy*s&
                       &in(a*xx*2.0D0)*2.0D0-a**2*xx*yy*sin(a*xx*2.0D0)*2.0D0+a**2*ym*yy*s&
                       &in(a*yy*2.0D0)*4.0D0-a**2*cy*sx*xm**2*2.0D0-a**2*cy*sx*xx**2*2.0D0&
                       &-a**2*cy*sx*ym**2*2.0D0-a**2*cy*sx*yy**2*2.0D0+a**2*sx*sy*xm**2*1.&
                       &0D+1+a**2*sx*sy*xx**2*1.0D+1+a**2*sx*sy*ym**2*1.0D+1+a**2*sx*sy*yy&
                       &**2*1.0D+1-a*cx*cy*xm*2.0D0+a*cx*cy*xx*2.0D0+a**2*cx2*cy*sy*xm**2*&
                       &4.0D0+a**2*cx2*cy*sy*xx**2*4.0D0+a**2*cx2*cy*sy*ym**2*4.0D0+a**2*c&
                       &x2*cy*sy*yy**2*4.0D0-a**2*cx*cy*xm*ym*2.0D+1+a**2*cx*cy*xx*ym*2.0D&
                       &+1+a**2*cx*cy*xm*yy*2.0D+1-a**2*cx*cy*xx*yy*2.0D+1+a**2*cy*sx*xm*x&
                       &x*4.0D0-a**2*cx*sy*xm*ym*4.0D0+a**2*cx*sy*xx*ym*4.0D0+a**2*cx*sy*x&
                       &m*yy*4.0D0-a**2*cx*sy*xx*yy*4.0D0+a**2*cy*sx*ym*yy*4.0D0-a**2*sx*s&
                       &y*xm*xx*2.0D+1-a**2*sx*sy*ym*yy*2.0D+1-a**2*cx2*cy*sy*xm*xx*8.0D0+&
                       &a**2*cx*cy2*sx*xm*ym*8.0D0-a**2*cx*cy2*sx*xx*ym*8.0D0-a**2*cx*cy2*&
                       &sx*xm*yy*8.0D0+a**2*cx*cy2*sx*xx*yy*8.0D0-a**2*cx2*cy*sy*ym*yy*8.0&
                       &D0-a*cx*cy*sx*sy*xm*2.0D0+a*cx*cy*sx*sy*xx*2.0D0))/(xm*xx*(-2.0D0)&
                       &-ym*yy*2.0D0+xm**2+xx**2+ym**2+yy**2+1.0D0)-Mref*cx*cy*pcourr*((a*&
                       &sx*(xm-xx)*(cy*1.0D+1+sx+sy*2.0D0-cy2*sx*2.0D0)*(2.0D0/3.0D0))/(Mr&
                       &ef*r)+(a*cx*(ym-yy)*(cy-sy*5.0D0+cy*sx*sy)*(4.0D0/3.0D0))/(Mref*r)&
                       &)+(a*cx2*cy*(ym-yy)*(cy-sy*5.0D0+cy*sx*sy)*(1.0D+1/3.0D0))/r-(sqrt&
                       &(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0-2.0D+1)/Mref)*(sx*sy+2.0D0)**2*(-&
                       &cx2*cy2+cx*sy*2.0D0+cy*sx*2.0D0+2.0D+1))/(tie*(cy*sx-1.0D+1)*2.0D0&
                       &)+(a*cx*cy*sx*(xm-xx)*(cy*1.0D+1+sx+sy*2.0D0-cy2*sx*2.0D0)*(5.0D0/&
                       &3.0D0))/r+(a*cx*sy*(cy*sx-1.0D+1)*(sx*sy+2.0D0)*(xm-xx)*(5.0D0/3.0&
                       &D0))/r-(a*cy*sx*(cy*sx-1.0D+1)*(sx*sy+2.0D0)*(ym-yy)*(5.0D0/3.0D0)&
                       &)/r+(cx*cy*(xm*2.0D0-xx*2.0D0)*(cy*sx-1.0D+1)*(sx*sy+2.0D0)*(ym-yy&
                       &)*(5.0D0/6.0D0))/r2-(cx*cy*(ym*2.0D0-yy*2.0D0)*(cy*sx-1.0D+1)*(sx*&
                       &sy+2.0D0)*(xm-xx)*(5.0D0/6.0D0))/r2+(3.0D0**(-epn-1.0D0)*a*kpare*(&
                       &-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(cx*cy*xm*(-1.0D+1)+cx*cy*xx*1.0D&
                       &+1+sx*sy*ym*1.0D+1-sx*sy*yy*1.0D+1+cx*cy2*sx*xm-cx*cy2*sx*xx-cy*sx&
                       &2*sy*ym+cy*sx2*sy*yy+a*cy*sx*xm**2*1.0D+1-a*cy2*sx2*xm**2+a*cy*sx*&
                       &xx**2*1.0D+1-a*cy2*sx2*xx**2+a*cy*sx*ym**2*1.0D+1-a*cy2*sx2*ym**2+&
                       &a*cy*sx*yy**2*1.0D+1-a*cy2*sx2*yy**2-a*cy*sx*xm*xx*2.0D+1+a*cy2*sx&
                       &2*xm*xx*2.0D0-a*cx*sy*xm*ym*2.0D+1+a*cx*sy*xx*ym*2.0D+1+a*cx*sy*xm&
                       &*yy*2.0D+1-a*cx*sy*xx*yy*2.0D+1-a*cy*sx*ym*yy*2.0D+1+a*cy2*sx2*ym*&
                       &yy*2.0D0+a*cx2*cy2*epn*ym**2+a*cx2*cy2*epn*yy**2+a*epn*sx2*sy2*xm*&
                       &*2+a*epn*sx2*sy2*xx**2-a*cx2*cy2*epn*ym*yy*2.0D0-a*epn*sx2*sy2*xm*&
                       &xx*2.0D0+a*cx*cy*sx*sy*xm*ym*2.0D0-a*cx*cy*sx*sy*xx*ym*2.0D0-a*cx*&
                       &cy*sx*sy*xm*yy*2.0D0+a*cx*cy*sx*sy*xx*yy*2.0D0+a*cx*cy*epn*sx*sy*x&
                       &m*ym*2.0D0-a*cx*cy*epn*sx*sy*xx*ym*2.0D0-a*cx*cy*epn*sx*sy*xm*yy*2&
                       &.0D0+a*cx*cy*epn*sx*sy*xx*yy*2.0D0)*2.0D0)/(Mref*(cy*sx-1.0D+1)*(x&
                       &m*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2+ym**2+yy**2+1.0D0))
           
        
           
           
           CASE(2)
						     ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
						        IF (.not.switch%axisym) THEN
						           WRITE(6,*) "This is an axisymmetric test case!"
						           stop
						        END IF
									     a = 2*pi
									     b = 2*pi

									     D     = phys%diff_n
									     mu    = phys%diff_u
									     csii  = phys%diff_e
									     csie  = phys%diff_ee
									     kpari = phys%diff_pari
									     kpare = phys%diff_pare
									     Mref  = phys%Mref
									     epn   = phys%epn
									     tie   = phys%tie
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


              f(ind,1) = (1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy&
                         &**2)**2*(cc*r*xx*ym**3*(-2.0D0)-cc*r*xx**3*ym*2.0D0+cc*r*xx*yy**3*&
                         &2.0D0+D*a**2*ss*xx**5*6.0D0+D*a**2*cc*xx**4*ym*4.0D0-D*a**2*cc*xx*&
                         &*4*yy*4.0D0-D*a**2*ss*xm*xx**4*1.0D+1+D*a**2*ss*xm**4*xx+D*a**2*ss&
                         &*xx*ym**4+D*a**2*ss*xx*yy**4+D*a**2*cc*xx**2*ym**3*2.0D0-D*a**2*cc&
                         &*xx**2*yy**3*2.0D0+D*a**2*ss*xm**2*xx**3*9.0D0-D*a**2*ss*xm**3*xx*&
                         &*2*4.0D0+D*a**2*ss*xx**3*ym**2*5.0D0+D*a**2*ss*xx**3*yy**2*5.0D0-D&
                         &*a*cx*sy*xm**4-D*a*cx*sy*xx**4*6.0D0+a*cx*r2*sy*xx**5*2.0D0+cx*cy*&
                         &r*xx**3*yy*2.0D0-cc*r*ss*xx*ym**3-cc*r*ss*xx**3*ym+cc*r*ss*xx*yy**&
                         &3+cc*r*ss*xx**3*yy+cc*r*xm*xx**2*ym*4.0D0-cc*r*xm**2*xx*ym*2.0D0-c&
                         &c*r*xm*xx**2*yy*4.0D0+cc*r*xm**2*xx*yy*2.0D0-cc*r*xx*ym*yy**2*6.0D&
                         &0+cc*r*xx*ym**2*yy*6.0D0-D*a**2*ss*xm*xx**2*ym**2*4.0D0+D*a**2*ss*&
                         &xm**2*xx*ym**2*2.0D0-D*a**2*ss*xm*xx**2*yy**2*4.0D0+D*a**2*ss*xm**&
                         &2*xx*yy**2*2.0D0+D*a**2*ss*xx*ym**2*yy**2*6.0D0+D*a*cx*sy*xm*xx**3&
                         &*1.2D+1+D*a*cx*sy*xm**3*xx*5.0D0-D*a*cy*sx*xm*ym**3-D*a*cy*sx*xm**&
                         &3*ym+D*a*cy*sx*xx*ym**3*2.0D0+D*a*cy*sx*xx**3*ym*2.0D0+D*a*cy*sx*x&
                         &m*yy**3+D*a*cy*sx*xm**3*yy-D*a*cy*sx*xx*yy**3*2.0D0-D*a*cy*sx*xx**&
                         &3*yy*2.0D0-a*cx*cy2*r2*sx*xx**5+a*cx*r2*sx*sy2*xx**5-a*cx*r2*sy*xm&
                         &*xx**4*2.0D0+a*cy*r2*sx*xx**4*ym*2.0D0-a*cy*r2*sx*xx**4*yy*2.0D0+c&
                         &c*r*ss*xm*xx**2*ym*2.0D0-cc*r*ss*xm**2*xx*ym-cc*r*ss*xm*xx**2*yy*2&
                         &.0D0+cc*r*ss*xm**2*xx*yy-cc*r*ss*xx*ym*yy**2*3.0D0+cc*r*ss*xx*ym**&
                         &2*yy*3.0D0-D*a*cx*sy*xm**2*xx**2*1.1D+1-D*a*cx*sy*xm**2*ym**2-D*a*&
                         &cx*sy*xx**2*ym**2*5.0D0-D*a*cx*sy*xm**2*yy**2-D*a*cx*sy*xx**2*yy**&
                         &2*5.0D0-D*a**2*cc*xm*xx*ym**3*2.0D0-D*a**2*cc*xm*xx**3*ym*8.0D0-D*&
                         &a**2*cc*xm**3*xx*ym*2.0D0+D*a**2*cc*xm*xx*yy**3*2.0D0+D*a**2*cc*xm&
                         &*xx**3*yy*8.0D0-D*a**2*ss*xx*ym*yy**3*4.0D0-D*a**2*ss*xx*ym**3*yy*&
                         &4.0D0-D*a**2*ss*xx**3*ym*yy*1.0D+1+D*a**2*cc*xm*xx*ym**2*yy*6.0D0+&
                         &D*a**2*ss*xm*xx**2*ym*yy*8.0D0-D*a**2*ss*xm**2*xx*ym*yy*4.0D0+D*a*&
                         &*2*cx*cy*xm**2*xx**2*ym*6.0D0-D*a**2*cx*cy*xm**2*xx**2*yy*6.0D0+D*&
                         &a**2*cx*cy*xx**2*ym*yy**2*6.0D0-D*a**2*cx*cy*xx**2*ym**2*yy*6.0D0+&
                         &D*a*cx*sy*xm*xx*ym**2*3.0D0-D*a*cy*sx*xm*xx**2*ym*4.0D0+D*a*cy*sx*&
                         &xm**2*xx*ym*4.0D0+D*a*cx*sy*xm*xx*yy**2*3.0D0+D*a*cy*sx*xm*xx**2*y&
                         &y*4.0D0-D*a*cy*sx*xm**2*xx*yy*4.0D0+D*a*cx*sy*xm**2*ym*yy*2.0D0-D*&
                         &a*cy*sx*xm*ym*yy**2*3.0D0+D*a*cy*sx*xm*ym**2*yy*3.0D0+D*a*cx*sy*xx&
                         &**2*ym*yy*1.0D+1+D*a*cy*sx*xx*ym*yy**2*6.0D0-D*a*cy*sx*xx*ym**2*yy&
                         &*6.0D0+a*cx*cy2*r2*sx*xm*xx**4-a*cx2*cy*r2*sy*xx**4*ym+a*cx2*cy*r2&
                         &*sy*xx**4*yy-a*cx*r2*sx*sy2*xm*xx**4+a*cy*r2*sx2*sy*xx**4*ym-a*cy*&
                         &r2*sx2*sy*xx**4*yy+D*a**2*cx*cy*xm**3*xx*yy*2.0D0-D*a**2*cx*cy*xm*&
                         &xx*ym*yy**2*6.0D0-D*a*cx*sy*xm*xx*ym*yy*6.0D0))/xx
             
									             


              f(ind,2) =  mu*((xx*(a**2*cc*(ss+2.0D0)+a**2*cc*ss*3.0D0-((ym-yy)*((a*cx*&
                         &sy*(ss+2.0D0)-a*cx*cy2*sx)/(r*xx)+((xm-xx)*(a**2*ss*(ss+2.0D0)+a**&
                         &2*cx2*cy2-a**2*cx2*sy2-a**2*cy2*sx2))/(r*xx)+((a**2*cc*(ss+2.0D0)+&
                         &a**2*cc*ss*3.0D0)*(ym-yy))/(r*xx)+(1.0D0/xx**2*(xm-xx)*(a*cx*sy*(s&
                         &s+2.0D0)-a*cx*cy2*sx))/r-(1.0D0/xx**2*(ym-yy)*(a*cy*sx*(ss+2.0D0)-&
                         &a*cx2*cy*sy))/r-((xm-xx)*(a*cx*sy*(ss+2.0D0)-a*cx*cy2*sx)*(1.0D0/x&
                         &x**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*&
                         &(ym-yy)**2*2.0D0))/(r2*xx*2.0D0)+((ym-yy)*(a*cy*sx*(ss+2.0D0)-a*cx&
                         &2*cy*sy)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(xm-xx)**2*2&
                         &.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0))/(r2*xx*2.0D0)))/(r*xx)-(1.0D0/&
                         &xx**2*(((xm-xx)*(a*cx*sy*(ss+2.0D0)-a*cx*cy2*sx))/(r*xx)-((ym-yy)*&
                         &(a*cy*sx*(ss+2.0D0)-a*cx2*cy*sy))/(r*xx))*(ym-yy))/r+((((xm-xx)*(a&
                         &*cx*sy*(ss+2.0D0)-a*cx*cy2*sx))/(r*xx)-((ym-yy)*(a*cy*sx*(ss+2.0D0&
                         &)-a*cx2*cy*sy))/(r*xx))*(ym-yy)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1&
                         &.0D0/xx**3*(xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0))/(r2*xx*&
                         &2.0D0))+((((xm-xx)*(a*cx*sy*(ss+2.0D0)-a*cx*cy2*sx))/(r*xx)-((ym-y&
                         &y)*(a*cy*sx*(ss+2.0D0)-a*cx2*cy*sy))/(r*xx))*(ym-yy))/(r*xx)+a*cy*&
                         &sx*(ss+2.0D0)-a*cx2*cy*sy)/xx+a**2*cc*(ss+2.0D0)+a**2*cc*ss*3.0D0-&
                         &((xm-xx)*((a*cy*sx*(ss+2.0D0)-a*cx2*cy*sy)/(r*xx)+((ym-yy)*(a**2*s&
                         &s*(ss+2.0D0)+a**2*cx2*cy2-a**2*cx2*sy2-a**2*cy2*sx2))/(r*xx)+((a**&
                         &2*cc*(ss+2.0D0)+a**2*cc*ss*3.0D0)*(xm-xx))/(r*xx)+(1.0D0/xx**3*(ym&
                         &*2.0D0-yy*2.0D0)*(xm-xx)*(a*cx*sy*(ss+2.0D0)-a*cx*cy2*sx))/(r2*2.0&
                         &D0)-(1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*(ym-yy)*(a*cy*sx*(ss+2.0D0)-a&
                         &*cx2*cy*sy))/(r2*2.0D0)))/(r*xx)-(1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*&
                         &(((xm-xx)*(a*cx*sy*(ss+2.0D0)-a*cx*cy2*sx))/(r*xx)-((ym-yy)*(a*cy*&
                         &sx*(ss+2.0D0)-a*cx2*cy*sy))/(r*xx))*(xm-xx))/(r2*2.0D0))-((cx2*cy2&
                         &*(ym-yy)*(ss+2.0D0)*(1.0D0/xx**2*(xm*2.0D0-xx*2.0D0)+1.0D0/xx**3*(&
                         &xm-xx)**2*2.0D0+1.0D0/xx**3*(ym-yy)**2*2.0D0))/(r2*2.0D0)+(a*cx**3&
                         &*cy2*sy*(ym-yy))/r-(a*cx*cy2*sx*(ym-yy)*(ss+2.0D0)*2.0D0)/r)/xx+Mr&
                         &ef*(((xm-xx)*(((ss*2.0D0+4.0D0)*(a*cx*cy+a*cx2*cy*sy))/(Mref*3.0D0&
                         &)+(a*ss*(ss*2.0D0+4.0D0))/(Mref*3.0D0)+(a*cy*sx*(cx2*cy2*(-1.0D0/2&
                         &.0D0)+cx*sy+2.0D+1)*(2.0D0/3.0D0))/Mref-(a*cy*sx*(cy*sx-1.0D+1)*(2&
                         &.0D0/3.0D0))/Mref))/(r*xx)+((ym-yy)*(((ss*2.0D0+4.0D0)*(a*ss-a*cx*&
                         &cy2*sx))/(Mref*3.0D0)+(a*cc*(ss*2.0D0+4.0D0))/(Mref*3.0D0)-(a*cx*s&
                         &y*(cx2*cy2*(-1.0D0/2.0D0)+cx*sy+2.0D+1)*(2.0D0/3.0D0))/Mref+(a*cx*&
                         &sy*(cy*sx-1.0D+1)*(2.0D0/3.0D0))/Mref))/(r*xx))+(cx2*cy2*1.0D0/xx*&
                         &*3*(ym*2.0D0-yy*2.0D0)*(xm-xx)*(ss+2.0D0))/(r2*2.0D0)+(a*cx2*cy**3&
                         &*sx*(xm-xx))/(r*xx)-(a*cx2*cy*sy*(xm-xx)*(ss+2.0D0)*2.0D0)/(r*xx)
									              
             
              f(ind,3) = -kpari*(((3.0D0**(-epn-1.0D0)*xx*(ym-yy)*((-cx2*cy2+cx*sy*2.0&
                         &D0+4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
                         &2*2.0D0+ym**2+yy**2)**2*(a**2*cy2*ym**3*(-2.0D0)+a**2*cy2*yy**3*2.&
                         &0D0-a*ss*xm*ym*2.0D0+a*ss*xx*ym*4.0D0+a*ss*xm*yy*2.0D0-a*ss*xx*yy*&
                         &4.0D0+a*cx*cy*ym**2*2.0D0+a*cx*cy*yy**2*2.0D0+a**2*cx2*cy2*ym**3*4&
                         &.0D0-a**2*cx2*cy2*yy**3*4.0D0+a**2*cc*ss*xm**3*4.0D0-a**2*cc*ss*xx&
                         &**3*8.0D0+a**2*cy*sx*xm**3*2.0D0-a**2*cy*sx*xx**3*4.0D0-a**2*cx*sy&
                         &*ym**3*2.0D0+a**2*cx*sy*yy**3*2.0D0-a**2*cy2*xm**2*ym*2.0D0-a**2*c&
                         &y2*xx**2*ym*4.0D0+a**2*cy2*xm**2*yy*2.0D0+a**2*cy2*xx**2*yy*4.0D0-&
                         &a**2*cy2*ym*yy**2*6.0D0+a**2*cy2*ym**2*yy*6.0D0+a**2*cx2*cy2*xm**2&
                         &*ym*4.0D0+a**2*cx2*cy2*xx**2*ym*8.0D0-a**2*cx2*cy2*xm**2*yy*4.0D0-&
                         &a**2*cx2*cy2*xx**2*yy*8.0D0+a**2*cx2*cy2*ym*yy**2*1.2D+1-a**2*cx2*&
                         &cy2*ym**2*yy*1.2D+1+a**2*cc*ss*xm*xx**2*1.6D+1-a**2*cc*ss*xm**2*xx&
                         &*1.2D+1+a**2*cy*sx*xm*xx**2*8.0D0-a**2*cy*sx*xm**2*xx*6.0D0+a**2*c&
                         &c*ss*xm*ym**2*4.0D0-a**2*cc*ss*xx*ym**2*4.0D0+a**2*cc*ss*xm*yy**2*&
                         &4.0D0-a**2*cc*ss*xx*yy**2*4.0D0-a**2*cx*sy*xm**2*ym*2.0D0+a**2*cy*&
                         &sx*xm*ym**2*2.0D0-a**2*cx*sy*xx**2*ym*4.0D0-a**2*cy*sx*xx*ym**2*2.&
                         &0D0+a**2*cx*sy*xm**2*yy*2.0D0+a**2*cy*sx*xm*yy**2*2.0D0+a**2*cx*sy&
                         &*xx**2*yy*4.0D0-a**2*cy*sx*xx*yy**2*2.0D0-a**2*cx*sy*ym*yy**2*6.0D&
                         &0+a**2*cx*sy*ym**2*yy*6.0D0+a*cx*cy*xm*xx*2.0D0-a*cx*cy*ym*yy*4.0D&
                         &0+a*cx2*cy*sy*ym**2*2.0D0+a*cx2*cy*sy*yy**2*2.0D0+a**2*cy2*xm*xx*y&
                         &m*4.0D0-a**2*cy2*xm*xx*yy*4.0D0-a**2*cx2*cy2*xm*xx*ym*8.0D0+a**2*c&
                         &x2*cy2*xm*xx*yy*8.0D0+a**2*cx*sy*xm*xx*ym*4.0D0-a**2*cx*sy*xm*xx*y&
                         &y*4.0D0-a**2*cc*ss*xm*ym*yy*8.0D0+a**2*cc*ss*xx*ym*yy*8.0D0-a**2*c&
                         &y*sx*xm*ym*yy*4.0D0+a**2*cy*sx*xx*ym*yy*4.0D0+a*cx2*cy*sy*xm*xx*2.&
                         &0D0+a*cx*cy2*sx*xm*ym*2.0D0-a*cx*cy2*sx*xx*ym*4.0D0-a*cx*cy2*sx*xm&
                         &*yy*2.0D0+a*cx*cy2*sx*xx*yy*4.0D0-a*cx2*cy*sy*ym*yy*4.0D0))/Mref-(&
                         &3.0D0**(-epn-1.0D0)*a*(ym-yy)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)&
                         &**epn*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)*1.0D0/(xm*xx*(-2.0D0)&
                         &-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(ss*ym-ss*yy+cx*cy*&
                         &xm-cx*cy*xx+cx2*cy*sy*xm-cx2*cy*sy*xx-cx*cy2*sx*ym+cx*cy2*sx*yy)*2&
                         &.0D0)/Mref+(3.0D0**(-epn-1.0D0)*1.0D0/Mref**2*a**2*epn*sx*xx*(ym-y&
                         &y)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*(sy-cx*cy2)*(&
                         &ss*ym-ss*yy+cx*cy*xm-cx*cy*xx+cx2*cy*sy*xm-cx2*cy*sy*xx-cx*cy2*sx*&
                         &ym+cx*cy2*sx*yy)*4.0D0)/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
                         &0D0+ym**2+yy**2))/xx-(3.0D0**(-epn-1.0D0)*(xm-xx)*((-cx2*cy2+cx*sy&
                         &*2.0D0+4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+&
                         &xx**2*2.0D0+ym**2+yy**2)**2*(a*ss*xm**2*2.0D0+a*ss*xx**2*4.0D0+a**&
                         &2*cx2*xm**3*2.0D0-a**2*cx2*xx**3*4.0D0-a*ss*xm*xx*4.0D0-a**2*cx2*c&
                         &y2*xm**3*4.0D0+a**2*cx2*cy2*xx**3*8.0D0+a**2*cx*sy*xm**3*2.0D0-a**&
                         &2*cx*sy*xx**3*4.0D0-a**2*cc*ss*ym**3*4.0D0+a**2*cc*ss*yy**3*4.0D0-&
                         &a**2*cy*sx*ym**3*2.0D0+a**2*cy*sx*yy**3*2.0D0+a**2*cx2*xm*xx**2*8.&
                         &0D0-a**2*cx2*xm**2*xx*6.0D0+a**2*cx2*xm*ym**2*2.0D0-a**2*cx2*xx*ym&
                         &**2*2.0D0+a**2*cx2*xm*yy**2*2.0D0-a**2*cx2*xx*yy**2*2.0D0-a**2*cx2&
                         &*cy2*xm*xx**2*1.6D+1+a**2*cx2*cy2*xm**2*xx*1.2D+1-a**2*cx2*cy2*xm*&
                         &ym**2*4.0D0+a**2*cx2*cy2*xx*ym**2*4.0D0-a**2*cx2*cy2*xm*yy**2*4.0D&
                         &0+a**2*cx2*cy2*xx*yy**2*4.0D0+a**2*cx*sy*xm*xx**2*8.0D0-a**2*cx*sy&
                         &*xm**2*xx*6.0D0-a**2*cc*ss*xm**2*ym*4.0D0-a**2*cc*ss*xx**2*ym*8.0D&
                         &0+a**2*cc*ss*xm**2*yy*4.0D0+a**2*cc*ss*xx**2*yy*8.0D0+a**2*cx*sy*x&
                         &m*ym**2*2.0D0-a**2*cy*sx*xm**2*ym*2.0D0-a**2*cx*sy*xx*ym**2*2.0D0-&
                         &a**2*cy*sx*xx**2*ym*4.0D0+a**2*cx*sy*xm*yy**2*2.0D0+a**2*cy*sx*xm*&
                         &*2*yy*2.0D0-a**2*cx*sy*xx*yy**2*2.0D0+a**2*cy*sx*xx**2*yy*4.0D0-a*&
                         &*2*cc*ss*ym*yy**2*1.2D+1+a**2*cc*ss*ym**2*yy*1.2D+1-a**2*cy*sx*ym*&
                         &yy**2*6.0D0+a**2*cy*sx*ym**2*yy*6.0D0-a*cx*cy*xm*ym*2.0D0+a*cx*cy*&
                         &xx*ym*2.0D0+a*cx*cy*xm*yy*2.0D0-a*cx*cy*xx*yy*2.0D0-a*cx*cy2*sx*xm&
                         &**2*2.0D0-a*cx*cy2*sx*xx**2*4.0D0-a**2*cx2*xm*ym*yy*4.0D0+a**2*cx2&
                         &*xx*ym*yy*4.0D0+a**2*cx2*cy2*xm*ym*yy*8.0D0-a**2*cx2*cy2*xx*ym*yy*&
                         &8.0D0+a**2*cc*ss*xm*xx*ym*8.0D0-a**2*cc*ss*xm*xx*yy*8.0D0+a**2*cy*&
                         &sx*xm*xx*ym*4.0D0-a**2*cy*sx*xm*xx*yy*4.0D0-a**2*cx*sy*xm*ym*yy*4.&
                         &0D0+a**2*cx*sy*xx*ym*yy*4.0D0+a*cx*cy2*sx*xm*xx*4.0D0-a*cx2*cy*sy*&
                         &xm*ym*2.0D0+a*cx2*cy*sy*xx*ym*2.0D0+a*cx2*cy*sy*xm*yy*2.0D0-a*cx2*&
                         &cy*sy*xx*yy*2.0D0))/Mref+(3.0D0**(-epn-1.0D0)*a*(ym*2.0D0-yy*2.0D0&
                         &)*(xm-xx)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**epn*1.0D0/(xm*xx*(&
                         &-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(ss*ym-ss*yy&
                         &+cx*cy*xm-cx*cy*xx+cx2*cy*sy*xm-cx2*cy*sy*xx-cx*cy2*sx*ym+cx*cy2*s&
                         &x*yy))/Mref+(3.0D0**(-epn-1.0D0)*1.0D0/Mref**2*a**2*cc*epn*(cx*sy+&
                         &1.0D0)*(xm-xx)*((-cx2*cy2+cx*sy*2.0D0+4.0D+1)/Mref)**(epn-1.0D0)*(&
                         &ss*ym-ss*yy+cx*cy*xm-cx*cy*xx+cx2*cy*sy*xm-cx2*cy*sy*xx-cx*cy2*sx*&
                         &ym+cx*cy2*sx*yy)*4.0D0)/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
                         &0D0+ym**2+yy**2))-((a*cc*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2&
                         &.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(ss*(-1.0D+1)+cx&
                         &*sy*1.0D+2+cx2*sy2*5.0D0-sx2*sy2*5.0D0+cx*cy2*sx*4.0D0-cx**3*cy2*s&
                         &y+cx*cy2*sx2*sy*2.0D0))/3.0D0-(a*cy*sx*(ym-yy)*(ss+2.0D0)*1.0D0/sq&
                         &rt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2&
                         &+yy**2))*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/3.0D0+(cc*1.0D0/xx**3*(ym-&
                         &yy)*(ss+2.0D0)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**&
                         &2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(-cx2*cy2+cx*sy*5.0D0+1&
                         &.0D+2)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2))/3.0D0)/xx+Mref*cc*p&
                         &courr*((a*sx*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*&
                         &2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cy*1.0D+1+sx+sy*2.0D0-cy2*s&
                         &x*2.0D0)*(2.0D0/3.0D0))/(Mref*xx)+(a*cx*(ym-yy)*1.0D0/sqrt(1.0D0/x&
                         &x**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(&
                         &cy-sy*5.0D0+cy*ss)*(4.0D0/3.0D0))/(Mref*xx))+(a*csii*1.0D0/(xm*xx*&
                         &(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(xm**2*xx**&
                         &2*1.1D+1+xm**2*ym**2+xx**2*ym**2*5.0D0+xm**2*yy**2+xx**2*yy**2*5.0&
                         &D0-cx2*xm**4*2.0D0-cy2*xm**4-cx2*xx**4*1.2D+1-cy2*xx**4*6.0D0+ss*x&
                         &m**4*2.0D0+ss*xx**4*1.2D+1-xm*xx**3*1.2D+1-xm**3*xx*5.0D0+xm**4+xx&
                         &**4*6.0D0+a*ss*xx**5*1.2D+2-cx*sy*xm**4*2.0D+1-cx*sy*xx**4*1.2D+2+&
                         &cx2*xm*xx**3*2.4D+1+cx2*xm**3*xx*1.0D+1+cy2*xm*xx**3*1.2D+1+cy2*xm&
                         &**3*xx*5.0D0-ss*xm*xx**3*2.4D+1-ss*xm**3*xx*1.0D+1-xm*xx*ym**2*3.0&
                         &D0-xm*xx*yy**2*3.0D0-xm**2*ym*yy*2.0D0-xx**2*ym*yy*1.0D+1+a*xx**5*&
                         &sin(a*xx*2.0D0)*1.0D+1-cx2*xm**2*xx**2*2.2D+1-cy2*xm**2*xx**2*1.1D&
                         &+1-cx2*xm**2*ym**2*2.0D0-cy2*xm**2*ym**2-cx2*xx**2*ym**2*1.0D+1-cx&
                         &2*xm**2*yy**2*2.0D0-cy2*xx**2*ym**2*5.0D0-cy2*xm**2*yy**2-cx2*xx**&
                         &2*yy**2*1.0D+1-cy2*xx**2*yy**2*5.0D0+ss*xm**2*xx**2*2.2D+1+ss*xm**&
                         &2*ym**2*2.0D0+ss*xx**2*ym**2*1.0D+1+ss*xm**2*yy**2*2.0D0+ss*xx**2*&
                         &yy**2*1.0D+1+cx2*cy2*xm**4*2.0D0+cx2*cy2*xx**4*1.2D+1+a*xm**2*xx**&
                         &3*sin(a*xx*2.0D0)*1.7D+1-a*xm**3*xx**2*sin(a*xx*2.0D0)*8.0D0+a*xx*&
                         &*3*ym**2*sin(a*xx*2.0D0)*7.0D0+a*xx**3*yy**2*sin(a*xx*2.0D0)*7.0D0&
                         &-a*xx**2*ym**3*sin(a*yy*2.0D0)*2.0D0+a*xx**2*yy**3*sin(a*yy*2.0D0)&
                         &*2.0D0+xm*xx*ym*yy*6.0D0+a*cx*sy*xx**5*1.2D+1-cx2*cy2*xm*xx**3*2.4&
                         &D+1-cx2*cy2*xm**3*xx*1.0D+1-cx*cy*xm*ym**3*2.0D0-cx*cy*xm**3*ym*2.&
                         &0D0+cx*cy*xx*ym**3*4.0D0+cx*cy*xx**3*ym*4.0D0+cx*cy*xm*yy**3*2.0D0&
                         &+cx*cy*xm**3*yy*2.0D0-cx*cy*xx*yy**3*4.0D0-cx*cy*xx**3*yy*4.0D0-a*&
                         &ss*xm*xx**4*2.0D+2+a*ss*xm**4*xx*2.0D+1+a*ss*xx*ym**4*2.0D+1+a*ss*&
                         &xx*yy**4*2.0D+1+cx*sy*xm*xx**3*2.4D+2+cx*sy*xm**3*xx*1.0D+2-cc*ss*&
                         &xm*ym**3*2.0D0-cc*ss*xm**3*ym*2.0D0+cc*ss*xx*ym**3*4.0D0+cc*ss*xx*&
                         &*3*ym*4.0D0+cc*ss*xm*yy**3*2.0D0+cc*ss*xm**3*yy*2.0D0-cc*ss*xx*yy*&
                         &*3*4.0D0-cc*ss*xx**3*yy*4.0D0-cy*sx*xm*ym**3*2.0D+1-cy*sx*xm**3*ym&
                         &*2.0D+1+cy*sx*xx*ym**3*4.0D+1+cy*sx*xx**3*ym*4.0D+1+cy*sx*xm*yy**3&
                         &*2.0D+1+cy*sx*xm**3*yy*2.0D+1-cy*sx*xx*yy**3*4.0D+1-cy*sx*xx**3*yy&
                         &*4.0D+1+cx2*xm*xx*ym**2*6.0D0+cy2*xm*xx*ym**2*3.0D0+cx2*xm*xx*yy**&
                         &2*6.0D0+cy2*xm*xx*yy**2*3.0D0+cx2*xm**2*ym*yy*4.0D0+cy2*xm**2*ym*y&
                         &y*2.0D0+cx2*xx**2*ym*yy*2.0D+1+cy2*xx**2*ym*yy*1.0D+1-ss*xm*xx*ym*&
                         &*2*6.0D0-ss*xm*xx*yy**2*6.0D0-ss*xm**2*ym*yy*4.0D0-ss*xx**2*ym*yy*&
                         &2.0D+1-a*xm*xx**4*sin(a*xx*2.0D0)*1.8D+1+a*xm**4*xx*sin(a*xx*2.0D0&
                         &)*2.0D0+a*xx*ym**4*sin(a*xx*2.0D0)+a*xx*yy**4*sin(a*xx*2.0D0)-a*xx&
                         &**4*ym*sin(a*yy*2.0D0)*4.0D0+a*xx**4*yy*sin(a*yy*2.0D0)*4.0D0+cx2*&
                         &cy2*xm**2*xx**2*2.2D+1+cx2*cy2*xm**2*ym**2*2.0D0+cx2*cy2*xx**2*ym*&
                         &*2*1.0D+1+cx2*cy2*xm**2*yy**2*2.0D0+cx2*cy2*xx**2*yy**2*1.0D+1+a*s&
                         &s*xm**2*xx**3*1.8D+2-a*ss*xm**3*xx**2*8.0D+1+a*ss*xx**3*ym**2*1.0D&
                         &+2+a*ss*xx**3*yy**2*1.0D+2-cx*sy*xm**2*xx**2*2.2D+2-cx*sy*xm**2*ym&
                         &**2*2.0D+1-cx*sy*xx**2*ym**2*1.0D+2-cx*sy*xm**2*yy**2*2.0D+1-cx*sy&
                         &*xx**2*yy**2*1.0D+2-cc*ss*xm*xx**2*ym*8.0D0+cc*ss*xm**2*xx*ym*8.0D&
                         &0+cc*ss*xm*xx**2*yy*8.0D0-cc*ss*xm**2*xx*yy*8.0D0+cx*sy*xm*xx*ym**&
                         &2*6.0D+1-cy*sx*xm*xx**2*ym*8.0D+1+cy*sx*xm**2*xx*ym*8.0D+1+cx*sy*x&
                         &m*xx*yy**2*6.0D+1+cy*sx*xm*xx**2*yy*8.0D+1-cy*sx*xm**2*xx*yy*8.0D+&
                         &1-cc*ss*xm*ym*yy**2*6.0D0+cc*ss*xm*ym**2*yy*6.0D0+cc*ss*xx*ym*yy**&
                         &2*1.2D+1-cc*ss*xx*ym**2*yy*1.2D+1+cx*sy*xm**2*ym*yy*4.0D+1-cy*sx*x&
                         &m*ym*yy**2*6.0D+1+cy*sx*xm*ym**2*yy*6.0D+1+cx*sy*xx**2*ym*yy*2.0D+&
                         &2+cy*sx*xx*ym*yy**2*1.2D+2-cy*sx*xx*ym**2*yy*1.2D+2+a*xm*xx*ym**3*&
                         &sin(a*yy*2.0D0)*2.0D0+a*xm*xx**3*ym*sin(a*yy*2.0D0)*8.0D0+a*xm**3*&
                         &xx*ym*sin(a*yy*2.0D0)*2.0D0-a*xx*ym*yy**3*sin(a*xx*2.0D0)*4.0D0-a*&
                         &xx*ym**3*yy*sin(a*xx*2.0D0)*4.0D0-a*xx**3*ym*yy*sin(a*xx*2.0D0)*1.&
                         &4D+1-a*xm*xx*yy**3*sin(a*yy*2.0D0)*2.0D0-a*xm*xx**3*yy*sin(a*yy*2.&
                         &0D0)*8.0D0-a*xm**3*xx*yy*sin(a*yy*2.0D0)*2.0D0+a*cx*cy*xx**2*ym**3&
                         &*4.0D+1-a*cx*cy*xx**2*yy**3*4.0D+1+a*cx*sy*xm**2*xx**3*1.8D+1-a*cx&
                         &*sy*xm**3*xx**2*8.0D0+a*cx*sy*xx**3*ym**2*1.0D+1-a*cy*sx*xx**2*ym*&
                         &*3*4.0D0+a*cx*sy*xx**3*yy**2*1.0D+1+a*cy*sx*xx**2*yy**3*4.0D0-a*ss&
                         &*xm*xx**2*ym**2*8.0D+1+a*ss*xm**2*xx*ym**2*4.0D+1-a*ss*xm*xx**2*yy&
                         &**2*8.0D+1+a*ss*xm**2*xx*yy**2*4.0D+1+a*ss*xx*ym**2*yy**2*1.2D+2-a&
                         &*xm*xx**2*ym**2*sin(a*xx*2.0D0)*6.0D0+a*xm**2*xx*ym**2*sin(a*xx*2.&
                         &0D0)*3.0D0-a*xm*xx**2*yy**2*sin(a*xx*2.0D0)*6.0D0+a*xm**2*xx*yy**2&
                         &*sin(a*xx*2.0D0)*3.0D0-a*xm**2*xx**2*ym*sin(a*yy*2.0D0)*6.0D0+a*xx&
                         &*ym**2*yy**2*sin(a*xx*2.0D0)*6.0D0+a*xm**2*xx**2*yy*sin(a*yy*2.0D0&
                         &)*6.0D0-a*xx**2*ym*yy**2*sin(a*yy*2.0D0)*6.0D0+a*xx**2*ym**2*yy*si&
                         &n(a*yy*2.0D0)*6.0D0-cx2*xm*xx*ym*yy*1.2D+1-cy2*xm*xx*ym*yy*6.0D0+s&
                         &s*xm*xx*ym*yy*1.2D+1-a*cx*cy2*sx*xx**5*2.4D+1+a*cx*cy*xx**4*ym*8.0&
                         &D+1-a*cx*cy*xx**4*yy*8.0D+1-a*cx*sy*xm*xx**4*2.0D+1+a*cx*sy*xm**4*&
                         &xx*2.0D0+a*cx*sy*xx*ym**4*2.0D0-a*cy*sx*xx**4*ym*8.0D0+a*cx*sy*xx*&
                         &yy**4*2.0D0+a*cy*sx*xx**4*yy*8.0D0-cx*cy*xm*xx**2*ym*8.0D0+cx*cy*x&
                         &m**2*xx*ym*8.0D0-cx2*cy2*xm*xx*ym**2*6.0D0+cx*cy*xm*xx**2*yy*8.0D0&
                         &-cx*cy*xm**2*xx*yy*8.0D0-cx2*cy2*xm*xx*yy**2*6.0D0-cx*cy*xm*ym*yy*&
                         &*2*6.0D0+cx*cy*xm*ym**2*yy*6.0D0-cx2*cy2*xm**2*ym*yy*4.0D0+cx*cy*x&
                         &x*ym*yy**2*1.2D+1-cx*cy*xx*ym**2*yy*1.2D+1-cx2*cy2*xx**2*ym*yy*2.0&
                         &D+1-a*ss*xx*ym*yy**3*8.0D+1-a*ss*xx*ym**3*yy*8.0D+1-a*ss*xx**3*ym*&
                         &yy*2.0D+2+cx2*cy2*xm*xx*ym*yy*1.2D+1-cx*sy*xm*xx*ym*yy*1.2D+2+a*cx&
                         &*cy2*sx*xm*xx**4*4.0D+1-a*cx*cy2*sx*xm**4*xx*4.0D0-a*cx*cy2*sx*xx*&
                         &ym**4*4.0D0+a*cx2*cy*sy*xx**4*ym*1.6D+1-a*cx*cy2*sx*xx*yy**4*4.0D0&
                         &-a*cx2*cy*sy*xx**4*yy*1.6D+1-a*cx*cy*xm*xx*ym**3*4.0D+1-a*cx*cy*xm&
                         &*xx**3*ym*1.6D+2-a*cx*cy*xm**3*xx*ym*4.0D+1+a*cx*cy*xm*xx*yy**3*4.&
                         &0D+1+a*cx*cy*xm*xx**3*yy*1.6D+2+a*cx*cy*xm**3*xx*yy*4.0D+1+a*cy*sx&
                         &*xm*xx*ym**3*4.0D0+a*cy*sx*xm*xx**3*ym*1.6D+1+a*cy*sx*xm**3*xx*ym*&
                         &4.0D0-a*cy*sx*xm*xx*yy**3*4.0D0-a*cy*sx*xm*xx**3*yy*1.6D+1-a*cy*sx&
                         &*xm**3*xx*yy*4.0D0-a*cx*sy*xx*ym*yy**3*8.0D0-a*cx*sy*xx*ym**3*yy*8&
                         &.0D0-a*cx*sy*xx**3*ym*yy*2.0D+1+a*ss*xm*xx**2*ym*yy*1.6D+2-a*ss*xm&
                         &**2*xx*ym*yy*8.0D+1-a*cx*cy2*sx*xm**2*xx**3*3.6D+1+a*cx*cy2*sx*xm*&
                         &*3*xx**2*1.6D+1-a*cx*cy2*sx*xx**3*ym**2*2.0D+1+a*cx2*cy*sy*xx**2*y&
                         &m**3*8.0D0-a*cx*cy2*sx*xx**3*yy**2*2.0D+1-a*cx2*cy*sy*xx**2*yy**3*&
                         &8.0D0+a*xm*xx**2*ym*yy*sin(a*xx*2.0D0)*1.2D+1-a*xm**2*xx*ym*yy*sin&
                         &(a*xx*2.0D0)*6.0D0+a*cx*cy*xm**2*xx**2*ym*1.2D+2-a*cx*cy*xm**2*xx*&
                         &*2*yy*1.2D+2+a*xm*xx*ym*yy**2*sin(a*yy*2.0D0)*6.0D0-a*xm*xx*ym**2*&
                         &yy*sin(a*yy*2.0D0)*6.0D0+a*cx*cy*xx**2*ym*yy**2*1.2D+2-a*cx*cy*xx*&
                         &*2*ym**2*yy*1.2D+2-a*cx*sy*xm*xx**2*ym**2*8.0D0+a*cx*sy*xm**2*xx*y&
                         &m**2*4.0D0-a*cy*sx*xm**2*xx**2*ym*1.2D+1-a*cx*sy*xm*xx**2*yy**2*8.&
                         &0D0+a*cx*sy*xm**2*xx*yy**2*4.0D0+a*cy*sx*xm**2*xx**2*yy*1.2D+1+a*c&
                         &x*sy*xx*ym**2*yy**2*1.2D+1-a*cy*sx*xx**2*ym*yy**2*1.2D+1+a*cy*sx*x&
                         &x**2*ym**2*yy*1.2D+1+a*cx*cy2*sx*xm*xx**2*ym**2*1.6D+1-a*cx*cy2*sx&
                         &*xm**2*xx*ym**2*8.0D0+a*cx2*cy*sy*xm**2*xx**2*ym*2.4D+1+a*cx*cy2*s&
                         &x*xm*xx**2*yy**2*1.6D+1-a*cx*cy2*sx*xm**2*xx*yy**2*8.0D0-a*cx2*cy*&
                         &sy*xm**2*xx**2*yy*2.4D+1-a*cx*cy2*sx*xx*ym**2*yy**2*2.4D+1+a*cx2*c&
                         &y*sy*xx**2*ym*yy**2*2.4D+1-a*cx2*cy*sy*xx**2*ym**2*yy*2.4D+1-a*cx2&
                         &*cy*sy*xm*xx*ym**3*8.0D0-a*cx2*cy*sy*xm*xx**3*ym*3.2D+1-a*cx2*cy*s&
                         &y*xm**3*xx*ym*8.0D0+a*cx2*cy*sy*xm*xx*yy**3*8.0D0+a*cx2*cy*sy*xm*x&
                         &x**3*yy*3.2D+1+a*cx2*cy*sy*xm**3*xx*yy*8.0D0+a*cx*cy2*sx*xx*ym*yy*&
                         &*3*1.6D+1+a*cx*cy2*sx*xx*ym**3*yy*1.6D+1+a*cx*cy2*sx*xx**3*ym*yy*4&
                         &.0D+1-a*cx*cy*xm*xx*ym*yy**2*1.2D+2+a*cx*cy*xm*xx*ym**2*yy*1.2D+2+&
                         &a*cx*sy*xm*xx**2*ym*yy*1.6D+1-a*cx*sy*xm**2*xx*ym*yy*8.0D0+a*cy*sx&
                         &*xm*xx*ym*yy**2*1.2D+1-a*cy*sx*xm*xx*ym**2*yy*1.2D+1-a*cx*cy2*sx*x&
                         &m*xx**2*ym*yy*3.2D+1+a*cx*cy2*sx*xm**2*xx*ym*yy*1.6D+1-a*cx2*cy*sy&
                         &*xm*xx*ym*yy**2*2.4D+1+a*cx2*cy*sy*xm*xx*ym**2*yy*2.4D+1))/xx+(sqr&
                         &t(3.0D0)*1.0D0/sqrt(-(cy*sx*2.0D0-2.0D+1)/Mref)*(ss+2.0D0)**2*(-cx&
                         &2*cy2+cx*sy*2.0D0+cy*sx*2.0D0+2.0D+1))/(tie*(cy*sx-1.0D+1)*2.0D0)+&
                         &(cc*1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*(xm-xx)*(ss+2.0D0)*1.0D0/(1.0D&
                         &0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
                         &)**(3.0D0/2.0D0)*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/6.0D0+(a*cx*cy2*(x&
                         &m-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
                         &*2*2.0D0+ym**2+yy**2))*(cx*1.0D+1+sx*1.0D+2+cx*ss*1.0D+1+cx2*sx*2.&
                         &0D0+cx2*sy*4.0D0-cx2*cy2*sx*3.0D0))/(xx*3.0D0)-(a*cx*sy*(xm-xx)*(s&
                         &s+2.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+&
                         &xx**2*2.0D0+ym**2+yy**2))*(-cx2*cy2+cx*sy*5.0D0+1.0D+2))/(xx*3.0D0&
                         &)		      

                     

              f(ind,4) =  (a*cx2*cy*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*y&
                         &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cy-sy*5.0D0+cy*ss)*(1.0D+&
                         &1/3.0D0)-a*cy*sx*(cy*sx-1.0D+1)*(ym-yy)*(ss+2.0D0)*1.0D0/sqrt(1.0D&
                         &0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
                         &)*(5.0D0/3.0D0)+cc*1.0D0/xx**3*(cy*sx-1.0D+1)*(ym-yy)*(ss+2.0D0)*1&
                         &.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym&
                         &**2+yy**2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)*&
                         &(5.0D0/3.0D0))/xx+(1.0D0/Mref**2*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0&
                         &+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(3.0D0**(-epn+1.0D0)*Mref*a*cc*&
                         &kpare*ym**4*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1&
                         &.0D0)*Mref*a*cc*kpare*yy**4*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D&
                         &0-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*ym*yy**3*(-(cy*sx*2.0D0-2.0D&
                         &+1)/Mref)**epn*1.6D+1-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*ym**3*yy&
                         &*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*1.6D+1+3.0D0**(-epn+1.0D0)*Mref&
                         &*a*kpare*ss*xm*ym**3*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0&
                         &**(-epn+1.0D0)*Mref*a*kpare*ss*xm**3*ym*(-(cy*sx*2.0D0-2.0D+1)/Mre&
                         &f)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*ss*xx*ym**3*(-(cy*s&
                         &x*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*&
                         &ss*xx**3*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+1&
                         &.0D0)*Mref*a*kpare*ss*xm*yy**3*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4&
                         &.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare*ss*xm**3*yy*(-(cy*sx*2.0D0-2&
                         &.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a*kpare*ss*xx*yy*&
                         &*3*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mre&
                         &f*a*kpare*ss*xx**3*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D&
                         &0**(-epn+1.0D0)*Mref*a*cc*kpare*xx**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)&
                         &**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2&
                         &.0D0-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*ym**2*(-(cy*sx*2.0D0-2.0D&
                         &+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2&
                         &+yy**2)*2.0D0-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*yy**2*(-(cy*sx*2&
                         &.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
                         &0D0+ym**2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*xm**2*y&
                         &m**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*M&
                         &ref*a*cc*kpare*xx**2*ym**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0&
                         &+3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*xm**2*yy**2*(-(cy*sx*2.0D0-2.&
                         &0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*xx**2*y&
                         &y**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*M&
                         &ref*a*cc*kpare*ym**2*yy**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*8.0D0&
                         &-3.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpare*sx*xx**3*(-(cy*sx*2.0D0-2.&
                         &0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym*&
                         &*2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*xm*xx*(-(cy*sx&
                         &*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*&
                         &2.0D0+ym**2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*ym*yy&
                         &*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm*&
                         &*2+xx**2*2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a*kpare&
                         &*ss*xm*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*&
                         &2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mre&
                         &f*a*kpare*ss*xm*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D&
                         &0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-3.0D0**(-epn+1&
                         &.0D0)*Mref*a*cc*kpare*xm*xx*ym**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**ep&
                         &n*8.0D0-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*xm*xx*yy**2*(-(cy*sx*2&
                         &.0D0-2.0D+1)/Mref)**epn*8.0D0-3.0D0**(-epn+1.0D0)*Mref*a*cc*kpare*&
                         &xm**2*ym*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*8.0D0-3.0D0**(-epn+1&
                         &.0D0)*Mref*a*cc*kpare*xx**2*ym*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**ep&
                         &n*8.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*xm*xx**2*ym*(-(cy*sx*2&
                         &.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*&
                         &xm**2*xx*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+2&
                         &.0D0)*Mref*a*kpare*ss*xm*xx**2*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**ep&
                         &n*4.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*xm**2*xx*yy*(-(cy*sx*2&
                         &.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*&
                         &xm*ym*yy**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0-3.0D0**(-epn+2&
                         &.0D0)*Mref*a*kpare*ss*xm*ym**2*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**ep&
                         &n*4.0D0-3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*xx*ym*yy**2*(-(cy*sx*2&
                         &.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+2.0D0)*Mref*a*kpare*ss*&
                         &xx*ym**2*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*4.0D0+3.0D0**(-epn+1&
                         &.0D0)*Mref*a*cc*kpare*xm*xx*ym*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**ep&
                         &n*1.6D+1+3.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpare*sx*xm*xx**2*(-(cy*&
                         &sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
                         &2*2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpare*&
                         &sx*xm**2*xx*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*y&
                         &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0-3.0D0**(-epn+1.0D0)*M&
                         &ref*a**2*cx*kpare*sy*xx**2*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(x&
                         &m*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0-3.0&
                         &D0**(-epn+1.0D0)*Mref*a**2*cy*kpare*sx*xx*ym**2*(-(cy*sx*2.0D0-2.0&
                         &D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
                         &2+yy**2)*2.0D0+3.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpare*sy*xx**2*yy*&
                         &(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**&
                         &2+xx**2*2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)*Mref*a**2*cy*&
                         &kpare*sx*xx*yy**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0&
                         &)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*2.0D0+(3.0D0**(-epn+1&
                         &.0D0)*Mref*a**2*epn*kpare*sx2*sy2*xx**3*(-(cy*sx*2.0D0-2.0D+1)/Mre&
                         &f)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
                         &*4.0D0)/(cy*sx*2.0D0-2.0D+1)+3.0D0**(-epn+1.0D0)*Mref*a**2*cx*kpar&
                         &e*sy*xm*xx*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym&
                         &*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0-3.0D0**(-epn+1.0D0)&
                         &*Mref*a**2*cx*kpare*sy*xm*xx*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*&
                         &(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0+3&
                         &.0D0**(-epn+1.0D0)*Mref*a**2*cy*kpare*sx*xx*ym*yy*(-(cy*sx*2.0D0-2&
                         &.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym&
                         &**2+yy**2)*4.0D0-(3.0D0**(-epn+1.0D0)*Mref*a**2*cc*epn*kpare*ss*xx&
                         &**2*ym*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0&
                         &D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(cy*sx*2.0D0-2.0D+1)+(3.0&
                         &D0**(-epn+1.0D0)*Mref*a**2*cc*epn*kpare*ss*xx**2*yy*(-(cy*sx*2.0D0&
                         &-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+&
                         &ym**2+yy**2)*8.0D0)/(cy*sx*2.0D0-2.0D+1)-(3.0D0**(-epn+1.0D0)*Mref&
                         &*a**2*epn*kpare*sx2*sy2*xm*xx**2*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn&
                         &*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)&
                         &/(cy*sx*2.0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*epn*kpare*sx2&
                         &*sy2*xm**2*xx*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym&
                         &*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*4.0D0)/(cy*sx*2.0D0-2.0D+&
                         &1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*cy2*epn*kpare*xx*ym**2*(-(cy&
                         &*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
                         &*2*2.0D0+ym**2+yy**2)*4.0D0)/(cy*sx*2.0D0-2.0D+1)+(3.0D0**(-epn+1.&
                         &0D0)*Mref*a**2*cx2*cy2*epn*kpare*xx*yy**2*(-(cy*sx*2.0D0-2.0D+1)/M&
                         &ref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**&
                         &2)*4.0D0)/(cy*sx*2.0D0-2.0D+1)-(3.0D0**(-epn+1.0D0)*Mref*a**2*cx2*&
                         &cy2*epn*kpare*xx*ym*yy*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-&
                         &2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(cy*sx*2.&
                         &0D0-2.0D+1)+(3.0D0**(-epn+1.0D0)*Mref*a**2*cc*epn*kpare*ss*xm*xx*y&
                         &m*(-(cy*sx*2.0D0-2.0D+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm&
                         &**2+xx**2*2.0D0+ym**2+yy**2)*8.0D0)/(cy*sx*2.0D0-2.0D+1)-(3.0D0**(&
                         &-epn+1.0D0)*Mref*a**2*cc*epn*kpare*ss*xm*xx*yy*(-(cy*sx*2.0D0-2.0D&
                         &+1)/Mref)**epn*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2&
                         &+yy**2)*8.0D0)/(cy*sx*2.0D0-2.0D+1)))/(xx*9.0D0)-Mref*cc*pcourr*((&
                         &a*sx*(xm-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm&
                         &**2+xx**2*2.0D0+ym**2+yy**2))*(cy*1.0D+1+sx+sy*2.0D0-cy2*sx*2.0D0)&
                         &*(2.0D0/3.0D0))/(Mref*xx)+(a*cx*(ym-yy)*1.0D0/sqrt(1.0D0/xx**2*(xm&
                         &*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cy-sy*5.&
                         &0D0+cy*ss)*(4.0D0/3.0D0))/(Mref*xx))+(a*csie*1.0D0/(xm*xx*(-2.0D0)&
                         &-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*(a*ss*xx**5*6.0D+1+&
                         &cc*ss*xm**4*2.0D0+cc*ss*xx**4*1.2D+1-cx*sy*xm**4*1.0D+1-cx*sy*xx**&
                         &4*6.0D+1-ss*xm*ym**3*2.0D0-ss*xm**3*ym*2.0D0+ss*xx*ym**3*4.0D0+ss*&
                         &xx**3*ym*4.0D0+ss*xm*yy**3*2.0D0+ss*xm**3*yy*2.0D0-ss*xx*yy**3*4.0&
                         &D0-ss*xx**3*yy*4.0D0+cx*cy*xm**4*2.0D0+cx*cy*xx**4*1.2D+1-a*cy*sx*&
                         &xx**5*1.2D+1-cx*cy*xm*xx**3*2.4D+1-cx*cy*xm**3*xx*1.0D+1-a*ss*xm*x&
                         &x**4*1.0D+2+a*ss*xm**4*xx*1.0D+1+a*ss*xx*ym**4*1.0D+1+a*ss*xx*yy**&
                         &4*1.0D+1-cc*ss*xm*xx**3*2.4D+1-cc*ss*xm**3*xx*1.0D+1+cx*sy*xm*xx**&
                         &3*1.2D+2+cx*sy*xm**3*xx*5.0D+1-cy*sx*xm*ym**3*1.0D+1-cy*sx*xm**3*y&
                         &m*1.0D+1+cy2*sx2*xm*ym**3+cy2*sx2*xm**3*ym+cy*sx*xx*ym**3*2.0D+1+c&
                         &y*sx*xx**3*ym*2.0D+1-cy2*sx2*xx*ym**3*2.0D0-cy2*sx2*xx**3*ym*2.0D0&
                         &+cy*sx*xm*yy**3*1.0D+1+cy*sx*xm**3*yy*1.0D+1-cy2*sx2*xm*yy**3-cy2*&
                         &sx2*xm**3*yy-cy*sx*xx*yy**3*2.0D+1-cy*sx*xx**3*yy*2.0D+1+cy2*sx2*x&
                         &x*yy**3*2.0D0+cy2*sx2*xx**3*yy*2.0D0-sx2*sy2*xm*ym**3-sx2*sy2*xm**&
                         &3*ym+sx2*sy2*xx*ym**3*2.0D0+sx2*sy2*xx**3*ym*2.0D0+sx2*sy2*xm*yy**&
                         &3+sx2*sy2*xm**3*yy-sx2*sy2*xx*yy**3*2.0D0-sx2*sy2*xx**3*yy*2.0D0-s&
                         &s*xm*xx**2*ym*8.0D0+ss*xm**2*xx*ym*8.0D0+ss*xm*xx**2*yy*8.0D0-ss*x&
                         &m**2*xx*yy*8.0D0-ss*xm*ym*yy**2*6.0D0+ss*xm*ym**2*yy*6.0D0+ss*xx*y&
                         &m*yy**2*1.2D+1-ss*xx*ym**2*yy*1.2D+1+cx*cy*xm**2*xx**2*2.2D+1+cx*c&
                         &y*xm**2*ym**2*2.0D0+cx*cy*xx**2*ym**2*1.0D+1+cx*cy*xm**2*yy**2*2.0&
                         &D0+cx*cy*xx**2*yy**2*1.0D+1+a*ss*xm**2*xx**3*9.0D+1-a*ss*xm**3*xx*&
                         &*2*4.0D+1+a*ss*xx**3*ym**2*5.0D+1+a*ss*xx**3*yy**2*5.0D+1+cc*ss*xm&
                         &**2*xx**2*2.2D+1-cx*sy*xm**2*xx**2*1.1D+2+cc*ss*xm**2*ym**2*2.0D0+&
                         &cc*ss*xx**2*ym**2*1.0D+1+cc*ss*xm**2*yy**2*2.0D0+cc*ss*xx**2*yy**2&
                         &*1.0D+1-cx*sy*xm**2*ym**2*1.0D+1-cx*sy*xx**2*ym**2*5.0D+1-cx*sy*xm&
                         &**2*yy**2*1.0D+1-cx*sy*xx**2*yy**2*5.0D+1-cc*ss*xm*xx*ym**2*6.0D0-&
                         &cc*ss*xm*xx*yy**2*6.0D0+cx*sy*xm*xx*ym**2*3.0D+1-cy*sx*xm*xx**2*ym&
                         &*4.0D+1+cy*sx*xm**2*xx*ym*4.0D+1+cy2*sx2*xm*xx**2*ym*4.0D0-cy2*sx2&
                         &*xm**2*xx*ym*4.0D0+cx*sy*xm*xx*yy**2*3.0D+1+cy*sx*xm*xx**2*yy*4.0D&
                         &+1-cy*sx*xm**2*xx*yy*4.0D+1-cy2*sx2*xm*xx**2*yy*4.0D0+cy2*sx2*xm**&
                         &2*xx*yy*4.0D0-cc*ss*xm**2*ym*yy*4.0D0-cc*ss*xx**2*ym*yy*2.0D+1+cx*&
                         &sy*xm**2*ym*yy*2.0D+1-cy*sx*xm*ym*yy**2*3.0D+1+cy*sx*xm*ym**2*yy*3&
                         &.0D+1+cy2*sx2*xm*ym*yy**2*3.0D0-cy2*sx2*xm*ym**2*yy*3.0D0+cx*sy*xx&
                         &**2*ym*yy*1.0D+2+cy*sx*xx*ym*yy**2*6.0D+1-cy*sx*xx*ym**2*yy*6.0D+1&
                         &-cy2*sx2*xx*ym*yy**2*6.0D0+cy2*sx2*xx*ym**2*yy*6.0D0-sx2*sy2*xm*xx&
                         &**2*ym*4.0D0+sx2*sy2*xm**2*xx*ym*4.0D0+sx2*sy2*xm*xx**2*yy*4.0D0-s&
                         &x2*sy2*xm**2*xx*yy*4.0D0-sx2*sy2*xm*ym*yy**2*3.0D0+sx2*sy2*xm*ym**&
                         &2*yy*3.0D0+sx2*sy2*xx*ym*yy**2*6.0D0-sx2*sy2*xx*ym**2*yy*6.0D0+a*c&
                         &x*cy*xx**2*ym**3*2.0D+1-a*cx*cy*xx**2*yy**3*2.0D+1-a*cy*sx*xm**2*x&
                         &x**3*1.8D+1+a*cy*sx*xm**3*xx**2*8.0D0+a*cx*sy*xx**2*ym**3*4.0D0-a*&
                         &cy*sx*xx**3*ym**2*1.0D+1-a*cx*sy*xx**2*yy**3*4.0D0-a*cy*sx*xx**3*y&
                         &y**2*1.0D+1-a*ss*xm*xx**2*ym**2*4.0D+1+a*ss*xm**2*xx*ym**2*2.0D+1-&
                         &a*ss*xm*xx**2*yy**2*4.0D+1+a*ss*xm**2*xx*yy**2*2.0D+1+a*ss*xx*ym**&
                         &2*yy**2*6.0D+1+a*cx2*cy*sy*xx**5*8.0D0+a*cx*cy*xx**4*ym*4.0D+1-a*c&
                         &x*cy*xx**4*yy*4.0D+1-a*cy*sx2*sy*xx**5*1.6D+1+a*cy*sx*xm*xx**4*2.0&
                         &D+1-a*cy*sx*xm**4*xx*2.0D0+a*cx*sy*xx**4*ym*8.0D0-a*cy*sx*xx*ym**4&
                         &*2.0D0-a*cx*sy*xx**4*yy*8.0D0-a*cy*sx*xx*yy**4*2.0D0-cx*cy*xm*xx*y&
                         &m**2*6.0D0-cx*cy*xm*xx*yy**2*6.0D0-cx*cy*xm**2*ym*yy*4.0D0-cx*cy*x&
                         &x**2*ym*yy*2.0D+1-a*ss*xx*ym*yy**3*4.0D+1-a*ss*xx*ym**3*yy*4.0D+1-&
                         &a*ss*xx**3*ym*yy*1.0D+2+cx*cy*xm*xx*ym*yy*1.2D+1+cc*ss*xm*xx*ym*yy&
                         &*1.2D+1-cx*sy*xm*xx*ym*yy*6.0D+1-a*cx2*cy*sy*xm*xx**4*1.6D+1+a*cx2&
                         &*cy*sy*xm**4*xx*2.0D0-a*cx*cy2*sx*xx**4*ym*8.0D0+a*cx*cy2*sx*xx**4&
                         &*yy*8.0D0-a*cx*cy*xm*xx*ym**3*2.0D+1-a*cx*cy*xm*xx**3*ym*8.0D+1-a*&
                         &cx*cy*xm**3*xx*ym*2.0D+1+a*cx*cy*xm*xx*yy**3*2.0D+1+a*cx*cy*xm*xx*&
                         &*3*yy*8.0D+1+a*cx*cy*xm**3*xx*yy*2.0D+1+a*cy*sx2*sy*xm*xx**4*2.4D+&
                         &1-a*cy*sx2*sy*xm**4*xx*2.0D0+a*cx*sx*sy2*xx**4*ym*8.0D0-a*cy*sx2*s&
                         &y*xx*ym**4*4.0D0-a*cx*sx*sy2*xx**4*yy*8.0D0-a*cy*sx2*sy*xx*yy**4*4&
                         &.0D0-a*cx*sy*xm*xx*ym**3*4.0D0-a*cx*sy*xm*xx**3*ym*1.6D+1-a*cx*sy*&
                         &xm**3*xx*ym*4.0D0+a*cx*sy*xm*xx*yy**3*4.0D0+a*cx*sy*xm*xx**3*yy*1.&
                         &6D+1+a*cx*sy*xm**3*xx*yy*4.0D0+a*cy*sx*xx*ym*yy**3*8.0D0+a*cy*sx*x&
                         &x*ym**3*yy*8.0D0+a*cy*sx*xx**3*ym*yy*2.0D+1+a*ss*xm*xx**2*ym*yy*8.&
                         &0D+1-a*ss*xm**2*xx*ym*yy*4.0D+1+a*cx2*cy*sy*xm**2*xx**3*1.6D+1-a*c&
                         &x2*cy*sy*xm**3*xx**2*8.0D0-a*cx*cy2*sx*xx**2*ym**3*4.0D0+a*cx2*cy*&
                         &sy*xx**3*ym**2*4.0D0+a*cx*cy2*sx*xx**2*yy**3*4.0D0+a*cx2*cy*sy*xx*&
                         &*3*yy**2*4.0D0+a*cx*cy*xm**2*xx**2*ym*6.0D+1-a*cx*cy*xm**2*xx**2*y&
                         &y*6.0D+1+a*cx*cy*xx**2*ym*yy**2*6.0D+1-a*cx*cy*xx**2*ym**2*yy*6.0D&
                         &+1-a*cy*sx2*sy*xm**2*xx**3*2.0D+1+a*cy*sx2*sy*xm**3*xx**2*8.0D0+a*&
                         &cx*sx*sy2*xx**2*ym**3*4.0D0-a*cy*sx2*sy*xx**3*ym**2*1.6D+1-a*cx*sx&
                         &*sy2*xx**2*yy**3*4.0D0-a*cy*sx2*sy*xx**3*yy**2*1.6D+1+a*cx*sy*xm**&
                         &2*xx**2*ym*1.2D+1+a*cy*sx*xm*xx**2*ym**2*8.0D0-a*cy*sx*xm**2*xx*ym&
                         &**2*4.0D0-a*cx*sy*xm**2*xx**2*yy*1.2D+1+a*cy*sx*xm*xx**2*yy**2*8.0&
                         &D0-a*cy*sx*xm**2*xx*yy**2*4.0D0+a*cx*sy*xx**2*ym*yy**2*1.2D+1-a*cx&
                         &*sy*xx**2*ym**2*yy*1.2D+1-a*cy*sx*xx*ym**2*yy**2*1.2D+1-a*cx*cy2*s&
                         &x*xm**2*xx**2*ym*1.2D+1-a*cx2*cy*sy*xm*xx**2*ym**2*4.0D0+a*cx2*cy*&
                         &sy*xm**2*xx*ym**2*2.0D0+a*cx*cy2*sx*xm**2*xx**2*yy*1.2D+1-a*cx2*cy&
                         &*sy*xm*xx**2*yy**2*4.0D0+a*cx2*cy*sy*xm**2*xx*yy**2*2.0D0-a*cx*cy2&
                         &*sx*xx**2*ym*yy**2*1.2D+1+a*cx*cy2*sx*xx**2*ym**2*yy*1.2D+1+a*cx*s&
                         &x*sy2*xm**2*xx**2*ym*1.2D+1+a*cy*sx2*sy*xm*xx**2*ym**2*1.2D+1-a*cy&
                         &*sx2*sy*xm**2*xx*ym**2*6.0D0-a*cx*sx*sy2*xm**2*xx**2*yy*1.2D+1+a*c&
                         &y*sx2*sy*xm*xx**2*yy**2*1.2D+1-a*cy*sx2*sy*xm**2*xx*yy**2*6.0D0+a*&
                         &cx*sx*sy2*xx**2*ym*yy**2*1.2D+1-a*cx*sx*sy2*xx**2*ym**2*yy*1.2D+1-&
                         &a*cy*sx2*sy*xx*ym**2*yy**2*2.4D+1+a*cx*cy2*sx*xm*xx*ym**3*4.0D0+a*&
                         &cx*cy2*sx*xm*xx**3*ym*1.6D+1+a*cx*cy2*sx*xm**3*xx*ym*4.0D0-a*cx*cy&
                         &2*sx*xm*xx*yy**3*4.0D0-a*cx*cy2*sx*xm*xx**3*yy*1.6D+1-a*cx*cy2*sx*&
                         &xm**3*xx*yy*4.0D0-a*cx2*cy*sy*xx**3*ym*yy*8.0D0-a*cx*cy*xm*xx*ym*y&
                         &y**2*6.0D+1+a*cx*cy*xm*xx*ym**2*yy*6.0D+1-a*cx*sx*sy2*xm*xx*ym**3*&
                         &4.0D0-a*cx*sx*sy2*xm*xx**3*ym*1.6D+1-a*cx*sx*sy2*xm**3*xx*ym*4.0D0&
                         &+a*cx*sx*sy2*xm*xx*yy**3*4.0D0+a*cx*sx*sy2*xm*xx**3*yy*1.6D+1+a*cx&
                         &*sx*sy2*xm**3*xx*yy*4.0D0+a*cy*sx2*sy*xx*ym*yy**3*1.6D+1+a*cy*sx2*&
                         &sy*xx*ym**3*yy*1.6D+1+a*cy*sx2*sy*xx**3*ym*yy*3.2D+1-a*cx*sy*xm*xx&
                         &*ym*yy**2*1.2D+1+a*cx*sy*xm*xx*ym**2*yy*1.2D+1-a*cy*sx*xm*xx**2*ym&
                         &*yy*1.6D+1+a*cy*sx*xm**2*xx*ym*yy*8.0D0+a*cx*cy2*sx*xm*xx*ym*yy**2&
                         &*1.2D+1-a*cx*cy2*sx*xm*xx*ym**2*yy*1.2D+1+a*cx2*cy*sy*xm*xx**2*ym*&
                         &yy*8.0D0-a*cx2*cy*sy*xm**2*xx*ym*yy*4.0D0-a*cx*sx*sy2*xm*xx*ym*yy*&
                         &*2*1.2D+1+a*cx*sx*sy2*xm*xx*ym**2*yy*1.2D+1-a*cy*sx2*sy*xm*xx**2*y&
                         &m*yy*2.4D+1+a*cy*sx2*sy*xm**2*xx*ym*yy*1.2D+1))/xx-(sqrt(3.0D0)*1.&
                         &0D0/sqrt(-(cy*sx*2.0D0-2.0D+1)/Mref)*(ss+2.0D0)**2*(-cx2*cy2+cx*sy&
                         &*2.0D0+cy*sx*2.0D0+2.0D+1))/(tie*(cy*sx-1.0D+1)*2.0D0)+(a*cc*sx*(x&
                         &m-xx)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx*&
                         &*2*2.0D0+ym**2+yy**2))*(cy*1.0D+1+sx+sy*2.0D0-cy2*sx*2.0D0)*(5.0D0&
                         &/3.0D0))/xx-cc*1.0D0/xx**3*(ym*2.0D0-yy*2.0D0)*(cy*sx-1.0D+1)*(xm-&
                         &xx)*(ss+2.0D0)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**&
                         &2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(5.0D0/6.0D0)+(a*cx*sy*&
                         &(cy*sx-1.0D+1)*(xm-xx)*(ss+2.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-&
                         &2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(5.0D0/3.0D0))/&
                         &xx    


           CASE(3)
						     ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
						        IF (.not.switch%axisym) THEN
						           WRITE(6,*) "This is an axisymmetric test case!"
						           stop
						        END IF
									     a = 2*pi
									     b = 2*pi

									     D     = phys%diff_n
									     mu    = phys%diff_u
									     csii  = phys%diff_e
									     csie  = phys%diff_ee
									     kpari = phys%diff_pari
									     kpare = phys%diff_pare
									     Mref  = phys%Mref
									     epn   = phys%epn
									     tie   = phys%tie
									     pcourr = 1.
									     
              st = sin(a*tt)
              ct = cos(a*tt)
              r = (xm**2-2*ym*yy-2*xm*xx+2*xx**2+ym**2+yy**2)
              
              
              f(ind,1) = D*((a**2*1.0D0/xx**2*sin(a*tt)*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm&
          &**2+xx**2+ym**2+yy**2))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.&
          &0D0+ym**2+yy**2)+a*cos(a*tt)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(xm&
          &*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2-(a*cos(&
          &a*tt)*(ym-yy)*(ym*yy*(-2.0D0)+xm**2-xx**2*2.0D0+ym**2+yy**2)*1.0D0&
          &/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/xx&
          &)-(a*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
          &2*2.0D0+ym**2+yy**2))*(sin(a*tt)*2.0D0+sin(a*tt)**2*2.0D0-1.0D0))/&
          &xx-1.0D0/xx**4*cos(a*tt)*(sin(a*tt)+2.0D0)*(ym-yy)*1.0D0/(1.0D0/xx&
          &**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(&
          &3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym**2+yy**2)+(1.0D0/xx**3*c&
          &os(a*tt)*(sin(a*tt)+2.0D0)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D&
          &0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
          &)**(3.0D0/2.0D0))/2.0D0
          
          f(ind,2) = mu*(-a*(ym*2.0D0-yy*2.0D0)*(xm-xx)*(sin(a*tt)*2.0D0+sin(a*tt)&
          &**2*2.0D0-1.0D0)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0&
          &D0+ym**2+yy**2)**2+(a**2*1.0D0/xx**2*(cos(a*tt)*4.0D0+cos(a*tt)*si&
          &n(a*tt)*8.0D0)*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2+ym**2+yy**2&
          &))/(xm*xx*(-4.0D0)-ym*yy*4.0D0+xm**2*2.0D0+xx**2*4.0D0+ym**2*2.0D0&
          &+yy**2*2.0D0)+(a*(ym-yy)*(sin(a*tt)*2.0D0+sin(a*tt)**2*2.0D0-1.0D0&
          &)*(ym*yy*(-2.0D0)+xm**2-xx**2*2.0D0+ym**2+yy**2)*1.0D0/(xm*xx*(-2.&
          &0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/xx)-(a*1.0D0/s&
          &qrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**&
          &2+yy**2))*(cos(a*tt)*(-5.8D+1)+sin(a*tt)*4.0D0-cos(a*tt)**2*4.0D0+&
          &cos(a*tt)**3*3.0D0+2.0D0))/(xx*3.0D0)-(a*cos(a*tt)*1.0D0/sqrt(1.0D&
          &0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
          &)*(sin(a*tt)*4.0D0+sin(a*tt)**2*3.0D0-1.0D0))/xx+(1.0D0/xx**3*cos(&
          &a*tt)**2*(sin(a*tt)+2.0D0)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D&
          &0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)&
          &)**(3.0D0/2.0D0))/2.0D0-1.0D0/xx**4*cos(a*tt)**2*(sin(a*tt)+2.0D0)&
          &*(ym-yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**&
          &2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(-xm*xx-ym*yy*2.0D0+xm**2+ym*&
          &*2+yy**2)
                     
                     
          f(ind,3) = csii*(a*(ym*2.0D0-yy*2.0D0)*(xm-xx)*(cos(a*tt)*2.0D+1+cos(a*t&
          &t*2.0D0)-sin(a*tt)*2.0D0)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+&
          &xx**2*2.0D0+ym**2+yy**2)**2+(a**2*1.0D0/xx**2*(cos(a*tt)*2.0D0+sin&
          &(a*tt)*2.0D+1+sin(a*tt*2.0D0)*2.0D0)*(xm*xx*(-2.0D0)-ym*yy*2.0D0+x&
          &m**2+xx**2+ym**2+yy**2))/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2&
          &.0D0+ym**2+yy**2)-(a*(ym-yy)*(cos(a*tt)*2.0D+1+cos(a*tt*2.0D0)-sin&
          &(a*tt)*2.0D0)*(ym*yy*(-2.0D0)+xm**2-xx**2*2.0D0+ym**2+yy**2)*1.0D0&
          &/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2)/xx&
          &)-kpari*((3.0D0**(-epn-1.0D0)*a**2*((cos(a*tt)*2.0D0-cos(a*tt)**2+&
          &4.0D+1)/Mref)**epn*(cos(a*tt)-1.0D0)*(epn*(-2.0D0)+cos(a*tt)*8.2D+&
          &1+epn*cos(a*tt)*2.0D0+epn*cos(a*tt)**2*2.0D0-epn*cos(a*tt)**3*2.0D&
          &0+cos(a*tt)**2*3.0D0-cos(a*tt)**3*2.0D0+4.0D+1)*2.0D0)/(Mref*(cos(&
          &a*tt)*2.0D0-cos(a*tt)**2+4.0D+1)*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2&
          &+xx**2*2.0D0+ym**2+yy**2))+(3.0D0**(-epn-1.0D0)*a*sin(a*tt)*(ym*2.&
          &0D0-yy*2.0D0)*(xm-xx)*((cos(a*tt)*2.0D0-cos(a*tt)**2+4.0D+1)/Mref)&
          &**epn*(cos(a*tt)-1.0D0)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx&
          &**2*2.0D0+ym**2+yy**2)**2*2.0D0)/Mref-(3.0D0**(-epn-1.0D0)*a*sin(a&
          &*tt)*(ym-yy)*((cos(a*tt)*2.0D0-cos(a*tt)**2+4.0D+1)/Mref)**epn*(co&
          &s(a*tt)-1.0D0)*(ym*yy*(-2.0D0)+xm**2-xx**2*2.0D0+ym**2+yy**2)*1.0D&
          &0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*2.&
          &0D0)/(Mref*xx))-(a*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.&
          &0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(cos(a*tt)*1.0D+1+sin(a*tt)*1.&
          &94D+2+sin(a*tt*2.0D0)*1.0D+1+sin(a*tt)**3*6.0D0-cos(a*tt)**2*2.03D&
          &+2-cos(a*tt)**3*1.5D+1+cos(a*tt)**4*4.0D0+1.0D+2))/(xx*3.0D0)+(sqr&
          &t(3.0D0)*(sin(a*tt)+2.0D0)**2*1.0D0/sqrt(-(sin(a*tt)*2.0D0-2.0D+1)&
          &/Mref)*(-cos(a*tt*2.0D0)+sqrt(2.0D0)*sin(3.141592653589793D0/4.0D0&
          &+a*tt)*4.0D0+3.9D+1))/(tie*(sin(a*tt)-1.0D+1)*4.0D0)-(a*pcourr*cos&
          &(a*tt)**2*(sin(a*tt)-4.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)&
          &-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(4.0D0/3.0D0))/xx-(1.&
          &0D0/xx**4*cos(a*tt)*(sin(a*tt)+2.0D0)*(ym-yy)*1.0D0/(1.0D0/xx**2*(&
          &xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0&
          &/2.0D0)*(cos(a*tt)*5.0D0-cos(a*tt)**2+1.0D+2)*(-xm*xx-ym*yy*2.0D0+&
          &xm**2+ym**2+yy**2))/3.0D0+(1.0D0/xx**3*cos(a*tt)*(sin(a*tt)+2.0D0)&
          &*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym&
          &*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(cos(a*tt&
          &)*5.0D0-cos(a*tt)**2+1.0D+2))/6.0D0
                     
                     
                     
          f(ind,4) = csie*((a**2*1.0D0/xx**2*(sin(a*tt)*4.0D0-sin(a*tt)**2*2.0D0+1&
          &.0D0)*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2+ym**2+yy**2)*2.0D0)/&
          &(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)-a*cos(a&
          &*tt)*(sin(a*tt)-4.0D0)*(ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(xm*xx*(-&
          &2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2)**2*2.0D0+(a*cos(&
          &a*tt)*(sin(a*tt)-4.0D0)*(ym-yy)*(ym*yy*(-2.0D0)+xm**2-xx**2*2.0D0+&
          &ym**2+yy**2)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+y&
          &m**2+yy**2)**2*2.0D0)/xx)+kpare*((3.0D0**(-epn-1.0D0)*a**2*(-(sin(&
          &a*tt)*2.0D0-2.0D+1)/Mref)**epn*(epn+sin(a*tt)*1.0D+1-sin(a*tt)**2-&
          &epn*sin(a*tt)**2)*2.0D0)/(Mref*(sin(a*tt)-1.0D+1)*(xm*xx*(-2.0D0)-&
          &ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))+(3.0D0**(-epn-1.0D0)*a&
          &*cos(a*tt)*(ym*2.0D0-yy*2.0D0)*(-(sin(a*tt)*2.0D0-2.0D+1)/Mref)**e&
          &pn*(xm-xx)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym*&
          &*2+yy**2)**2*2.0D0)/Mref-(3.0D0**(-epn-1.0D0)*a*cos(a*tt)*(-(sin(a&
          &*tt)*2.0D0-2.0D+1)/Mref)**epn*(ym-yy)*(ym*yy*(-2.0D0)+xm**2-xx**2*&
          &2.0D0+ym**2+yy**2)*1.0D0/(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2&
          &.0D0+ym**2+yy**2)**2*2.0D0)/(Mref*xx))-(a*1.0D0/sqrt(1.0D0/xx**2*(&
          &xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))*(sin(a*&
          &tt)*2.2D+1+sin(a*tt)**2*1.6D+1-sin(a*tt)**3*3.0D0-8.0D0)*(5.0D0/3.&
          &0D0))/xx-(sqrt(3.0D0)*(sin(a*tt)+2.0D0)**2*1.0D0/sqrt(-(sin(a*tt)*&
          &2.0D0-2.0D+1)/Mref)*(-cos(a*tt*2.0D0)+sqrt(2.0D0)*sin(3.1415926535&
          &89793D0/4.0D0+a*tt)*4.0D0+3.9D+1))/(tie*(sin(a*tt)-1.0D+1)*4.0D0)-&
          &1.0D0/xx**4*cos(a*tt)*(ym-yy)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-y&
          &m*yy*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(sin(a*t&
          &t)*(4.0D+1/3.0D0)-sin(a*tt)**2*(5.0D0/3.0D0)+1.0D+2/3.0D0)*(-xm*xx&
          &-ym*yy*2.0D0+xm**2+ym**2+yy**2)+(a*pcourr*cos(a*tt)**2*(sin(a*tt)-&
          &4.0D0)*1.0D0/sqrt(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*yy*2.0D0+xm**2+xx&
          &**2*2.0D0+ym**2+yy**2))*(4.0D0/3.0D0))/xx+(1.0D0/xx**3*cos(a*tt)*(&
          &ym*2.0D0-yy*2.0D0)*(xm-xx)*1.0D0/(1.0D0/xx**2*(xm*xx*(-2.0D0)-ym*y&
          &y*2.0D0+xm**2+xx**2*2.0D0+ym**2+yy**2))**(3.0D0/2.0D0)*(sin(a*tt)*&
          &(4.0D+1/3.0D0)-sin(a*tt)**2*(5.0D0/3.0D0)+1.0D+2/3.0D0))/2.0D0

  
   
   
   
               

           CASE(5:6)
           ! Do nothing
				       CASE(50:)
				       ! Do nothing

				       CASE DEFAULT
				           WRITE(6,*) "Error! Test case not valid"
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
  real*8, dimension(:), intent(IN)        :: x,y
  real*8, dimension(:,:), intent(OUT)     :: u
  real*8, dimension(size(u,1),phys%npv)  :: up
  integer:: i
  real*8 :: a,r(size(x)),th(size(x))
  
  up = 0.
  a = 2*pi
				SELECT CASE(switch%testcase)
				  CASE(1)
				  				IF (switch%axisym) THEN
						      WRITE(6,*) "This is NOT an axisymmetric test case!"
						      stop
						   END IF
				  ! Cartesian case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
								  up(:,1) = 2+sin(a*x)*sin(a*y)
								  up(:,2) = cos(a*x)*cos(a*y)
								  up(:,3) = 20+cos(a*x)*sin(a*y)
								  up(:,4) = 10-sin(a*x)*cos(a*y)
						CASE(2)
						! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
						   IF (.not.switch%axisym) THEN
						      WRITE(6,*) "This is an axisymmetric test case!"
						      stop
						   END IF 
								  up(:,1) = 2+sin(a*x)*sin(a*y)
								  up(:,2) = cos(a*x)*cos(a*y)
								  up(:,3) = 20+cos(a*x)*sin(a*y)
								  up(:,4) = 10-sin(a*x)*cos(a*y)
						CASE(5:6)
						! Cartesian case, square mesh, horizontal field 
				      up(:,1) = 1.
				      up(:,3) = 18.
				      up(:,4) = 18.     
				  CASE(50:64)
				      up(:,1) = 1.
				      up(:,3) = 18.
				      up(:,4) = 18.
				      
				      
!				      r = 	sqrt ( (x*phys%lscale-geom%R0)**2 + (y*phys%lscale)**2 )
!				      th = atan(x-geom%R0,y)
!				      up(:,1) = 2.+ cos(r*2*pi)*sin(th*4*pi)
!				      up(:,2) = 10. + cos(r*2*pi)*sin(th*4*pi)
!				      up(:,3) = 36. + cos(r*2*pi)*sin(th*4*pi)
!				      up(:,4) = 36. + cos(r*2*pi)*sin(th*4*pi)
				      
				      
				      
				  CASE(65)
				      up(:,1) = 1.
				      r = 	sqrt ( (x*phys%lscale-geom%R0)**2 + (y*phys%lscale-0.75)**2 )
				      DO i=1,size(x)
				         IF (r(i).le. 0.05) THEN
				            up(i,2) = 1.
				         END IF
				      END DO			      
				  CASE DEFAULT
				      WRITE(6,*) "Error! Test case not valid"
				      STOP
				END SELECT
				! Convert physical variables to conservative variables
				CALL phys2cons(up,u)
  END SUBROUTINE analytical_solution



		!*****************************************
		! Analytical gradient  
		!**************************************** 
  SUBROUTINE analytical_gradient(x,y,u,ux,uy)
  real*8, dimension(:), intent(IN)        :: x,y
  real*8, dimension(:,:), intent(IN)      :: u
  real*8, dimension(:,:), intent(OUT)     :: ux,uy
  real*8, dimension(size(u,1),size(u,2))  :: upx,upy
  real*8, dimension(size(u,1),phys%npv)  :: up
  real*8 :: a
  
  upx = 0.
  upy = 0.
  CALL cons2phys(u,up)
  a = 2*pi
				SELECT CASE(switch%testcase)
				  CASE(1)
				  				IF (switch%axisym) THEN
						      WRITE(6,*) "This is NOT an axisymmetric test case!"
						      stop
						   END IF
				  ! Circular field centered in [xc, yc], n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
        upx(:,1) =  a*cos(a*x)*sin(a*y)
        upx(:,2) = -a*sin(a*x)*cos(a*y)
        upx(:,3) = -a*sin(a*x)*sin(a*y)
        upx(:,4) = -a*cos(a*x)*cos(a*y)

        upy(:,1) =  a*sin(a*x)*cos(a*y)
        upy(:,2) = -a*cos(a*x)*sin(a*y)
        upy(:,3) =  a*cos(a*x)*cos(a*y)
        upy(:,4) =  a*sin(a*x)*sin(a*y)
								  
						CASE(2)
						! Axisimmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
						   IF (.not.switch%axisym) THEN
						      WRITE(6,*) "This is an axisymmetric test case!"
						      stop
						   END IF 
        upx(:,1) =  a*cos(a*x)*sin(a*y)
        upx(:,2) = -a*sin(a*x)*cos(a*y)
        upx(:,3) = -a*sin(a*x)*sin(a*y)
        upx(:,4) = -a*cos(a*x)*cos(a*y)

        upy(:,1) =  a*sin(a*x)*cos(a*y)
        upy(:,2) = -a*cos(a*x)*sin(a*y)
        upy(:,3) =  a*cos(a*x)*cos(a*y)
        upy(:,4) =  a*sin(a*x)*sin(a*y)
      CASE(5)
          ! Do nothing
      CASE(6)
          ! Do nothing
				  CASE(50:)
				      ! Do nothing   
				  CASE DEFAULT
				      WRITE(6,*) "Error! Test case not valid"
				      STOP
				END SELECT
				! Convert physical variables to conservative variables
				ux(:,1) = upx(:,1) 
				uy(:,1) = upy(:,1)
				ux(:,2) = ( upx(:,1)*up(:,2)+up(:,1)*upx(:,2) )
				uy(:,2) = ( upy(:,1)*up(:,2)+up(:,1)*upy(:,2) )
				ux(:,3) = ( upx(:,1)*up(:,3)+up(:,1)*upx(:,3) )
				uy(:,3) = ( upy(:,1)*up(:,3)+up(:,1)*upy(:,3) )
				ux(:,4) = ( upx(:,1)*up(:,4)+up(:,1)*upx(:,4) )
				uy(:,4) = ( upy(:,1)*up(:,4)+up(:,1)*upy(:,4) )

  END SUBROUTINE analytical_gradient
  
  

		!*****************************************
		! Body forces
		!**************************************** 
  SUBROUTINE body_force(x,y,f)
  real*8, dimension(:), intent(IN) :: x,y
  real*8, dimension(:,:), intent(OUT) :: f
  integer                            ::  n,i
  real*8  :: a,b,xc,yc,D,mu,csii,csie,kpari,kpare,Mref,epn,tie,pcourr
  real*8  :: xmax,xmin,ymax,ymin,xm,ym,xx,yy
  real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
  real*8 :: t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
  real*8 :: t21,t22,t23,t24,t25,t26,t27,t28,t29,t30
  real*8 :: t31,t32,t33,t34,t35,t36,t37,t38,t39,t40
  real*8 :: t41,t42,t43,t44,t45,t46,t47,t48,t49,t50
  real*8 :: t51,t52,t53,t54,t55,t56,t57,t58,t59,t60
  real*8 :: t61,t62,t63,t64,t65,t66,t67,t68,t69,t70
  real*8 :: cx2,sx2,cy2,sy2
  
  n = size(x)
  xmax = Mesh%xmax
  xmin = Mesh%xmin
  ymax = Mesh%ymax
  ymin = Mesh%ymin
  xm = 0.5*(xmax+xmin)
  ym = 0.5*(ymax+ymin)    
    f = 0.
				SELECT CASE(switch%testcase)
				  CASE(1)
				  ! Cartesian case, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
				  a = 2*pi
				  b = 2*pi
				  xc = 0.
				  yc = 0.
				  D     = phys%diff_n
				  mu    = phys%diff_u
				  csii  = phys%diff_e
				  csie  = phys%diff_ee
				  kpari = phys%diff_pari
				  kpare = phys%diff_pare
				  Mref  = phys%Mref
				  epn   = phys%epn
				  tie   = phys%tie
				  pcourr = 1.
				  


  
!  xm = -0.5
!  ym = -0.5

      DO i=1,size(x)
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
	          &2*D*a**2*xm*yy*t6-2*D*a**2*xx*yy*t6-2*D*a**2*xm*xx*t7-2*D*a**2*ym*yy*t7)/((xm-xx)**2+(ym-yy)**2+1)
	          
           f(i,2) = (a*(60*xm*t2*t3-60*xx*t2*t3-60*ym*t1*t3+60*yy*t1*t3+4*xm*t7*t3-4*xx*t7*t3+4*ym*t7*t3-4*yy*t7*t3+12*a*mu*t6-6*mu*xm*t2+&
           &6*mu*xx*t2-6*mu*ym*t1+6*mu*yy*t1-2*xm*cy2*sx2*t3+2*xx*cy2*sx2*t3-2*ym*cx2*sy2*t3+2*yy*cx2*sy2*t3+2*xm*sx2*sy2*t3-2*xx*sx2*sy2*t3+&
           &2*ym*sx2*sy2*t3-2*yy*sx2*sy2*t3+4*xm*t6*t3-4*xx*t6*t3+4*ym*t6*t3-4*yy*t6*t3+6*a*mu*xm**2*t6+6*a*mu*xx**2*t6+6*a*mu*ym**2*t6+&
           &6*a*mu*yy**2*t6+2*xm*cx2*cos(a*yy)**3*sin(a*xx)*t3-2*xx*cx2*cos(a*yy)**3*sin(a*xx)*t3-2*ym*cos(a*xx)**3*cy2*sin(a*yy)*t3+&
           &2*yy*cos(a*xx)**3*cy2*sin(a*yy)*t3+3*mu*ym*cos(a*xx)*cy2*sin(a*xx)+3*mu*xm*cx2*cos(a*yy)*sin(a*yy)-3*mu*xx*cx2*cos(a*yy)*sin(a*yy)-&
           &3*mu*yy*cos(a*xx)*cy2*sin(a*xx)-3*mu*ym*cos(a*xx)*sin(a*xx)*sy2-3*mu*xm*t2**2*sin(a*yy)+3*mu*xx*t2**2*sin(a*yy)+3*mu*yy*cos(a*xx)*sin(a*xx)*sy2+&
           &8*ym*cos(a*xx)*cy2*sin(a*xx)*t3-8*xm*cx2*cos(a*yy)*sin(a*yy)*t3+8*xx*cx2*cos(a*yy)*sin(a*yy)*t3-8*yy*cos(a*xx)*cy2*sin(a*xx)*t3-&
           &4*xm*cx2*t2*sy2*t3+4*xx*cx2*t2*sy2*t3+4*ym*cos(a*xx)*cy2*sx2*sin(a*yy)*t3-4*yy*cos(a*xx)*cy2*sx2*sin(a*yy)*t3-&
           &6*a*mu*xm*ym*t5+6*a*mu*xx*ym*t5+6*a*mu*xm*yy*t5-6*a*mu*xx*yy*t5+6*a*mu*xm*ym*cx2*sy2+6*a*mu*xm*ym*cy2*sx2-&
           &6*a*mu*xx*ym*cx2*sy2-6*a*mu*xx*ym*cy2*sx2-6*a*mu*xm*yy*cx2*sy2-6*a*mu*xm*yy*cy2*sx2+6*a*mu*xx*yy*cx2*sy2+&
           &6*a*mu*xx*yy*cy2*sx2-6*a*mu*xm*ym*sx2*sy2+6*a*mu*xx*ym*sx2*sy2+6*a*mu*xm*yy*sx2*sy2-6*a*mu*xx*yy*sx2*sy2+&
           &4*xm*cos(a*xx)*t2*sin(a*yy)*t3-4*xx*cos(a*xx)*t2*sin(a*yy)*t3+4*ym*cos(a*xx)*t2*sin(a*yy)*t3-4*yy*cos(a*xx)*t2*sin(a*yy)*t3-&
           &12*a*mu*xm*xx*t6-12*a*mu*ym*yy*t6-12*a*mu*xm*ym*t7+12*a*mu*xx*ym*t7+12*a*mu*xm*yy*t7-12*a*mu*xx*yy*t7+&
           &24*a*mu*cos(a*xx)*t2*sin(a*yy)+12*a*mu*xm**2*cos(a*xx)*t2*sin(a*yy)+12*a*mu*xx**2*cos(a*xx)*t2*sin(a*yy)+&
           &12*a*mu*ym**2*cos(a*xx)*t2*sin(a*yy)+12*a*mu*yy**2*cos(a*xx)*t2*sin(a*yy)-24*a*mu*xm*xx*cos(a*xx)*t2*sin(a*yy)-&
           &24*a*mu*ym*yy*cos(a*xx)*t2*sin(a*yy)))/(3*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1))

           f(i,3) = (csii*(a*xx-a*xm+3*a**2*sin(2*a*xx)+2*a*xm*cx2-2*a*xx*cx2+a*xm*cy2-a*xx*cy2+4*a**2*t1+&
           &40*a**2*t7+2*a**2*xm**2*sin(2*a*xx)+2*a**2*xx**2*sin(2*a*xx)+a**2*ym**2*sin(2*a*xx)+a**2*yy**2*sin(2*a*xx)-2*a*xm*t5+2*a*xx*t5+&
           &2*a**2*xm**2*t1+2*a**2*xx**2*t1+2*a**2*ym**2*t1+2*a**2*yy**2*t1+20*a**2*xm**2*t7+&
           &20*a**2*xx**2*t7+20*a**2*ym**2*t7+20*a**2*yy**2*t7-8*a**2*cos(a*xx)*cy2*sin(a*xx)+&
           &2*a*ym*t6-2*a*yy*t6+20*a*xm*t1-20*a*xx*t1+20*a*ym*t2-20*a*yy*t2-2*a*xm*t7+2*a*xx*t7-&
           &4*a**2*xm*xx*sin(2*a*xx)+2*a**2*xm*ym*sin(2*a*yy)-2*a**2*xx*ym*sin(2*a*yy)-&
           &2*a**2*ym*yy*sin(2*a*xx)-2*a**2*xm*yy*sin(2*a*yy)+2*a**2*xx*yy*sin(2*a*yy)-4*a**2*xm**2*cos(a*xx)*cy2*sin(a*xx)-&
           &4*a**2*xx**2*cos(a*xx)*cy2*sin(a*xx)-4*a**2*ym**2*cos(a*xx)*cy2*sin(a*xx)-&
           &4*a**2*yy**2*cos(a*xx)*cy2*sin(a*xx)-40*a**2*xm*ym*t6+40*a**2*xx*ym*t6+&
           &40*a**2*xm*yy*t6-40*a**2*xx*yy*t6-4*a**2*xm*xx*t1+4*a**2*xm*ym*t2-&
           &4*a**2*xx*ym*t2-4*a**2*xm*yy*t2+4*a**2*xx*yy*t2-4*a**2*ym*yy*t1-&
           &40*a**2*xm*xx*t7-40*a**2*ym*yy*t7+8*a**2*xm*xx*cos(a*xx)*cy2*sin(a*xx)-&
           &8*a**2*xm*ym*cx2*cos(a*yy)*sin(a*yy)+8*a**2*xx*ym*cx2*cos(a*yy)*sin(a*yy)+&
           &8*a**2*ym*yy*cos(a*xx)*cy2*sin(a*xx)+8*a**2*xm*yy*cx2*cos(a*yy)*sin(a*yy)-&
           &8*a**2*xx*yy*cx2*cos(a*yy)*sin(a*yy)+2*a*ym*cos(a*xx)*t2*sin(a*yy)-&
           &2*a*yy*cos(a*xx)*t2*sin(a*yy)))/(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)-&
           &(kpari*(2*3**(1-epn)*Mref*a**2*xx**3*cx2*(cy2-1)*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a**2*xx**3*t1*((2*t1-t5+40)/Mref)**epn+4*3**(1-epn)*&
           &a**2*epn*xx**3*t5*((2*t1-t5+40)/Mref)**(epn-1)+2*3**(1-epn)*Mref*a**2*&
           &xx**3*t5*((2*t1-t5+40)/Mref)**epn+2*3**(1-epn)*Mref*a*xx**2*t7*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a*ym**2*t7*((2*t1-t5+40)/Mref)**epn-2*3**(1-epn)*Mref*a*yy**2*t7*&
           &((2*t1-t5+40)/Mref)**epn+4*3**(1-epn)*a**2*epn*xx*ym**2*sx2*sy2*&
           &((2*t1-t5+40)/Mref)**(epn-1)+4*3**(1-epn)*a**2*epn*xx*yy**2*sx2*sy2*&
           &((2*t1-t5+40)/Mref)**(epn-1)-2*3**(1-epn)*Mref*a*xm*xx*t7*&
           &((2*t1-t5+40)/Mref)**epn+4*3**(1-epn)*Mref*a*ym*yy*t7*((2*t1-t5+40)/Mref)**epn+&
           &8*3**(1-epn)*a**2*epn*xx**3*cos(a*xx)**3*cy2*sin(a*yy)*((2*t1-t5+40)/Mref)**(epn-1)-&
           &4*3**(1-epn)*Mref*a**2*xm*xx**2*t5*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xm**2*xx*t5*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xx*ym**2*t5*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xx*yy**2*t5*((2*t1-t5+40)/Mref)**epn+&
           &4*3**(1-epn)*Mref*a**2*xm*xx**2*t1*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a**2*xm**2*xx*t1*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a**2*xx*ym**2*t1*((2*t1-t5+40)/Mref)**epn-&
           &4*3**(1-epn)*Mref*a**2*xx**2*ym*t2*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a**2*xx*yy**2*t1*((2*t1-t5+40)/Mref)**epn+&
           &4*3**(1-epn)*Mref*a**2*xx**2*yy*t2*((2*t1-t5+40)/Mref)**epn-&
           &4*3**(1-epn)*Mref*a**2*xm*xx**2*cx2*(cy2-1)*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xm**2*xx*cx2*(cy2-1)*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xx*ym**2*cy2*(cx2-1)*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a**2*xx*yy**2*cy2*(cx2-1)*((2*t1-t5+40)/Mref)**epn-&
           &2*3**(1-epn)*Mref*a*xx**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a*ym**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a*yy**2*cos(a*xx)*cy2*sin(a*xx)*((2*t1-t5+40)/Mref)**epn-&
           &8*3**(1-epn)*a**2*epn*xm*xx**2*t5*((2*t1-t5+40)/Mref)**(epn-1)+4*3**(1-epn)*&
           &a**2*epn*xm**2*xx*t5*((2*t1-t5+40)/Mref)**(epn-1)-2*3**(1-epn)*Mref*a*xm*ym*&
           &t6*((2*t1-t5+40)/Mref)**epn+4*3**(1-epn)*Mref*a*xx*ym*t6*((2*t1-t5+40)/Mref)**epn+&
           &2*3**(1-epn)*Mref*a*xm*yy*t6*((2*t1-t5+40)/Mref)**epn-4*3**(1-epn)*Mref*a*xx*yy*t6*&
           &((2*t1-t5+40)/Mref)**epn-8*3**(1-epn)*a**2*epn*xx*ym*yy*sx2*sy2*((2*t1-t5+40)/Mref)**(epn-1)-&
           &4*3**(1-epn)*Mref*a**2*xx*ym*yy*t5*((2*t1-t5+40)/Mref)**epn+8*3**(1-epn)*a**2*epn*xx**2*ym*&
           &cx2*cos(a*yy)**3*sin(a*xx)*((2*t1-t5+40)/Mref)**(epn-1)-16*3**(1-epn)*a**2*epn*xm*xx**2*&
           &cos(a*xx)**3*cy2*sin(a*yy)*((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*a**2*epn*xm**2*xx*&
           &cos(a*xx)**3*cy2*sin(a*yy)*((2*t1-t5+40)/Mref)**(epn-1)-8*3**(1-epn)*a**2*epn*xx**2*yy*&
           &cx2*cos(a*yy)**3*sin(a*xx)*((2*t1-t5+40)/Mref)**(epn-1)+4*3**(1-epn)*Mref*a**2*xm*xx*ym*&
           &t2*((2*t1-t5+40)/Mref)**epn-4*3**(1-epn)*Mref*a**2*xm*xx*yy*t2*((2*t1-t5+40)/Mref)**epn+&
           &4*3**(1-epn)*Mref*a**2*xx*ym*yy*t1*((2*t1-t5+40)/Mref)**epn-4*3**(1-epn)*Mref*a**2*xx*ym*&
           &yy*cy2*(cx2-1)*((2*t1-t5+40)/Mref)**epn+2*3**(1-epn)*Mref*a*xm*xx*cos(a*xx)*cy2*sin(a*xx)*&
           &((2*t1-t5+40)/Mref)**epn-2*3**(1-epn)*Mref*a*xm*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1-t5+40)/Mref)**epn+&
           &4*3**(1-epn)*Mref*a*xx*ym*cx2*cos(a*yy)*sin(a*yy)*((2*t1-t5+40)/Mref)**epn-4*3**(1-epn)*Mref*a*&
           &ym*yy*cos(a*xx)*cy2*sin(a*xx)*((2*t1-t5+40)/Mref)**epn+2*3**(1-epn)*Mref*a*xm*yy*cx2*&
           &cos(a*yy)*sin(a*yy)*((2*t1-t5+40)/Mref)**epn-4*3**(1-epn)*Mref*a*xx*yy*cx2*cos(a*yy)*&
           &sin(a*yy)*((2*t1-t5+40)/Mref)**epn-8*3**(1-epn)*Mref*a**2*xx**2*ym*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**epn+8*3**(1-epn)*Mref*a**2*xx**2*yy*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**epn+8*3**(1-epn)*a**2*epn*xx**2*ym*cos(a*xx)**3*cos(a*yy)**3*t7*&
           &((2*t1-t5+40)/Mref)**(epn-1)-8*3**(1-epn)*a**2*epn*xx**2*yy*cos(a*xx)**3*cos(a*yy)**3*&
           &t7*((2*t1-t5+40)/Mref)**(epn-1)-(4*3**(1-epn)*Mref*a**2*epn*xx**3*cos(a*xx)**4*cy2*(cy2-1)*&
           &((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)-8*3**(1-epn)*a**2*epn*xx**2*ym*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*a**2*epn*xx**2*yy*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**(epn-1)-8*3**(1-epn)*a**2*epn*xm*xx*ym*cx2*cos(a*yy)**3*sin(a*xx)*&
           &((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*a**2*epn*xm*xx*yy*cx2*cos(a*yy)**3*sin(a*xx)*&
           &((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*Mref*a**2*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**epn-8*3**(1-epn)*Mref*a**2*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**epn-8*3**(1-epn)*a**2*epn*xm*xx*ym*cos(a*xx)**3*cos(a*yy)**3*t7*&
           &((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*a**2*epn*xm*xx*yy*cos(a*xx)**3*cos(a*yy)**3*t7*&
           &((2*t1-t5+40)/Mref)**(epn-1)+8*3**(1-epn)*a**2*epn*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**(epn-1)-8*3**(1-epn)*a**2*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)*&
           &((2*t1-t5+40)/Mref)**(epn-1)-(4*3**(1-epn)*Mref*a**2*epn*xx*ym**2*cx2*cos(a*yy)**4*&
           &(cx2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+(8*3**(1-epn)*Mref*a**2*epn*xm*xx**2*&
           &cos(a*xx)**4*cy2*(cy2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)-(4*3**(1-epn)*Mref*&
           &a**2*epn*xm**2*xx*cos(a*xx)**4*cy2*(cy2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)-&
           &(4*3**(1-epn)*Mref*a**2*epn*xx*yy**2*cx2*cos(a*yy)**4*(cx2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+&
           &(8*3**(1-epn)*Mref*a**2*epn*xx*ym*yy*cx2*cos(a*yy)**4*(cx2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+&
           &(8*3**(1-epn)*Mref*a**2*epn*xx*ym**2*cos(a*xx)*cy2*sin(a*yy)*(cx2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+&
           &(8*3**(1-epn)*Mref*a**2*epn*xx**2*ym*cx2*t2*(cy2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+(8*3**(1-epn)*Mref*&
           &a**2*epn*xx*yy**2*cos(a*xx)*cy2*sin(a*yy)*(cx2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)-(8*3**(1-epn)*Mref*&
           &a**2*epn*xx**2*yy*cx2*t2*(cy2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)-(8*3**(1-epn)*Mref*a**2*epn*xm*xx*ym*&
           &cx2*t2*(cy2-1)*((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)+(8*3**(1-epn)*Mref*a**2*epn*xm*xx*yy*cx2*t2*(cy2-1)*&
           &((2*t1-t5+40)/Mref)**epn)/(2*t1-  t5+40)-(16*3**(1-epn)*Mref*a**2*epn*xx*ym*yy*cos(a*xx)*cy2*sin(a*yy)*(cx2-1)*&
           &((2*t1-t5+40)/Mref)**epn)/(2*t1-t5+40)))/(9*Mref**2*xx*((xm-xx)**2+(ym-yy)**2+1))+Mref*pcourr*t6*((2*a*sin(a*xx)*(xm-xx)*&
           &(10*cos(a*yy)+sin(a*xx)+2*sin(a*yy)-2*cy2*sin(a*xx)))/(3*Mref*t3)+(4*a*cos(a*xx)*(ym-yy)*(cos(a*yy)-5*sin(a*yy)+&
           &t2*sin(a*yy)))/(3*Mref*t3))-(a*t6*(ym-yy)*(5*cx2*sy2-5*sx2*sy2 +100*t1-10*t7-cos(a*xx)**3*cy2*sin(a*yy)+4*cos(a*xx)*&
           &cy2*sin(a*xx)+2*cos(a*xx)*cy2*sx2*sin(a*yy)))/(3*((xm-xx)**2+(ym-yy)**2+1)**(0.5))+(a*cos(a*xx)*cy2*(xm-xx)*&
           &(10*cos(a*xx)+100*sin(a*xx)+2*cx2*sin(a*xx)+4*cx2*sin(a*yy)-3*t5*sin(a*xx)+10*cos(a*xx)*t7))/(3*t3)+&
           &(3**(0.5)*(t7+2)**2*(2*t1 -t5+2*t2+20))/(2*tie*(t2-10)*(-(2*t2-20)/Mref)**(0.5))-&
           &(a*t1*(t7+2)*(xm-xx)*(5*t1-t5+100))/(3*((xm-xx)**2+(ym-yy)**2+1)**(0.5))+&
           &(a*t2*(t7+2)*(ym-yy)*(5*t1-t5+100))/(3*t3)-(t6*(2*xm-2*xx)*(t7+2)*(ym-yy)*(5*t1-t5+100))/(6*t4)+(t6*(2*ym-2*yy)*(t7+2)*(xm-xx)*(5*t1-t5+100))/(6*t4)
                      
          f(i,4) = (csie*(a*ym-a*yy-3*a**2*sin(2*a*yy)-a*ym*cx2+a*yy*cx2-2*a*ym*cy2+2*a*yy*cy2-&
          &4*a**2*t2+20*a**2*t7-a**2*xm**2*sin(2*a*yy)-a**2*xx**2*sin(2*a*yy)-2*a**2*ym**2*&
          &sin(2*a*yy)-2*a**2*yy**2*sin(2*a*yy)+2*a*ym*t5-2*a*yy*t5-2*a**2*xm**2*t2-&
          &2*a**2*xx**2*t2-2*a**2*ym**2*t2-2*a**2*yy**2*t2+10*a**2*xm**2*t7+10*a**2*xx**2*t7+&
          &10*a**2*ym**2*t7+10*a**2*yy**2*t7+8*a**2*cx2*cos(a*yy)*sin(a*yy)-2*a*xm*t6+2*a*&
          &xx*t6+10*a*xm*t1-10*a*xx*t1+10*a*ym*t2-10*a*yy*t2+2*a*ym*t7-2*a*yy*t7-2*a**2*&
          &xm*ym*sin(2*a*xx)+2*a**2*xx*ym*sin(2*a*xx)+2*a**2*xm*xx*sin(2*a*yy)+2*a**2*xm*&
          &yy*sin(2*a*xx)-2*a**2*xx*yy*sin(2*a*xx)+4*a**2*ym*yy*sin(2*a*yy)+4*a**2*xm**2*&
          &cx2*cos(a*yy)*sin(a*yy)+4*a**2*xx**2*cx2*cos(a*yy)*sin(a*yy)+4*a**2*ym**2*cx2*cos(a*yy)*&
          &sin(a*yy)+4*a**2*yy**2*cx2*cos(a*yy)*sin(a*yy)-20*a**2*xm*ym*t6+20*a**2*xx*ym*t6+20*a**2*xm*yy*t6-&
          &20*a**2*xx*yy*t6+4*a**2*xm*xx*t2-4*a**2*xm*ym*t1+4*a**2*xx*ym*t1+4*a**2*xm*yy*t1-4*&
          &a**2*xx*yy*t1+4*a**2*ym*yy*t2-20*a**2*xm*xx*t7-20*a**2*ym*yy*t7+8*a**2*xm*ym*cos(a*xx)*&
          &cy2*sin(a*xx)-8*a**2*xx*ym*cos(a*xx)*cy2*sin(a*xx)-8*a**2*xm*xx*cx2*cos(a*yy)*sin(a*yy)-&
          &8*a**2*xm*yy*cos(a*xx)*cy2*sin(a*xx)+8*a**2*xx*yy*cos(a*xx)*cy2*sin(a*xx)-8*a**2*ym*yy*&
          &cx2*cos(a*yy)*sin(a*yy)-2*a*xm*cos(a*xx)*t2*sin(a*yy)+2*a*xx*cos(a*xx)*t2*sin(a*yy)))/(xm**2-&
          &2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1)-Mref*pcourr*t6*((2*a*sin(a*xx)*(xm-xx)*(10*cos(a*yy)+&
          &sin(a*xx)+2*sin(a*yy)-2*cy2*sin(a*xx)))/(3*Mref*t3)+(4*a*cos(a*xx)*(ym-yy)*(cos(a*yy)-&
          &5*sin(a*yy)+t2*sin(a*yy)))/(3*Mref*t3))+(10*a*cx2*cos(a*yy)*(ym-yy)*(cos(a*yy)-5*sin(a*yy)+&
          &t2*sin(a*yy)))/(3*t3)-(3**(0.5)*(t7+2)**2*(2*t1-t5+2*t2+20))/(2*tie*(t2-10)*(-(2*t2-&
          &20)/Mref)**(0.5))+(5*a*cos(a*xx)*t2*(xm-xx)*(10*cos(a*yy)+sin(a*xx)+2*sin(a*yy)-2*cy2*&
          &sin(a*xx)))/(3*t3)+(5*t6*(t2-10)*(2*xm-2*xx)*(t7+2)*(ym-yy))/(6*t4)-(5*t6*(t2-10)*&
          &(2*ym-2*yy)*(t7+2)*(xm-xx))/(6*t4)+(5*a*t1*(t2-10)*(t7+2)*(xm-xx))/(3*t3)-(5*a*t2*(t2-10)*&
          &(t7+2)*(ym-yy))/(3*t3)+(2/3**(epn+1)*a*kpare*(-(2*t2-20)/Mref)**epn*(10*xx**2*t6-10*ym**2*&
          &t6-10*yy**2*t6-xx**2*cos(a*xx)*cy2*sin(a*xx)+ym**2*cos(a*xx)*cy2*sin(a*xx)+yy**2*cos(a*xx)*&
          &cy2*sin(a*xx)-10*xm*xx*t6+20*ym*yy*t6-a*xx**3*cy2*sx2-10*xm*ym*t7+20*xx*ym*t7+10*xm*yy*t7-&
          &20*xx*yy*t7+10*a*xx**3*t2+a*epn*xx**3*sx2*sy2+2*a*xm*xx**2*cy2*sx2-a*xm**2*xx*cy2*sx2-a*xx*&
          &ym**2*cy2*sx2-a*xx*yy**2*cy2*sx2-20*a*xm*xx**2*t2+10*a*xm**2*xx*t2+10*a*xx*ym**2*t2+20*a*&
          &xx**2*ym*t1+10*a*xx*yy**2*t2-20*a*xx**2*yy*t1+xm*xx*cos(a*xx)*cy2*sin(a*xx)-2*ym*yy*&
          &cos(a*xx)*cy2*sin(a*xx)+xm*ym*t2**2*sin(a*yy)-2*xx*ym*t2**2*sin(a*yy)-xm*yy*t2**2*&
          &sin(a*yy)+2*xx*yy*t2**2*sin(a*yy)+2*a*xx*ym*yy*cy2*sx2+a*epn*xx*ym**2*t5+a*epn*xx*yy**2*t5-20*a*xm*xx*ym*t1+20*a*xm*&
          &xx*yy*t1-20*a*xx*ym*yy*t2-2*a*epn*xm*xx**2*sx2*sy2+a*epn*xm**2*xx*sx2*sy2-2*a*epn*xx*&
          &ym*yy*t5-2*a*xx**2*ym*cos(a*xx)*t2*sin(a*yy)+2*a*xx**2*yy*cos(a*xx)*t2*sin(a*yy)-&
          &2*a*epn*xx**2*ym*cos(a*xx)*t2*sin(a*yy)+2*a*epn*xx**2*yy*cos(a*xx)*t2*sin(a*yy)+&
          &2*a*xm*xx*ym*cos(a*xx)*t2*sin(a*yy)-2*a*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)+2*a*epn*&
          &xm*xx*ym*cos(a*xx)*t2*sin(a*yy)-2*a*epn*xm*xx*yy*cos(a*xx)*t2*sin(a*yy)))/(Mref*&
          &xx*(t2-10)*(xm**2-2*ym*yy-2*xm*xx+xx**2+ym**2+yy**2+1))
           



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

      
      CASE(2)
      
      
      
      
      
      
      
      
      
      						     ! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
						        IF (.not.switch%axisym) THEN
						           WRITE(6,*) "This is an axisymmetric test case!"
						           stop
						        END IF
									     a = 2*pi
									     b = 2*pi

									     D     = phys%diff_n
									     mu    = phys%diff_u
									     csii  = phys%diff_e
									     csie  = phys%diff_ee
									     kpari = phys%diff_pari
									     kpare = phys%diff_pare
									     Mref  = phys%Mref
									     epn   = phys%epn
									     tie   = phys%tie
									     pcourr = 1.

      DO i=1,size(x)
         xx = x(i)
         yy = y(i)  
         
              t2 = xm**2
              t3 = xx**2
              t4 = a*xx
              t5 = cos(t4)
              t6 = a*yy
              t7 = sin(t6)
              t8 = t3**2
              t9 = xm-xx
              t10 = 1.0D0/xx**2
              t11 = ym-yy
              t12 = ym**2
              t13 = t9**2
              t14 = t10*t13
              t15 = t11**2
              t16 = t10*t15
              t17 = t14+t16+1.0D0
              t18 = cos(t6)
              t19 = sqrt(t17)
              t20 = yy**2
              t21 = a**2
              t22 = sin(t4)
              t23 = t2**2
              t24 = t17**(3.0D0/2.0D0)
              t25 = t18**2
              t26 = t5**2
              t27 = t7**2
              t28 = t22**2
              f(i,1) = (1.0D0/(t2+t3*2.0D0+t12+t20-xm*xx*2.0D0-ym*yy*2.0D0)**2*(D*a*&
              &t5*t7*t8*(-6.0D0)-D*a*t5*t7*t23-a*t5*t7*t8*t24*xm*2.0D0+a*t5*t7*t8&
              &*t24*xx*2.0D0+a*t8*t18*t22*t24*ym*2.0D0-a*t8*t18*t22*t24*yy*2.0D0+&
              &t3*t5*t18*t19*xm*ym*4.0D0-t2*t5*t18*t19*xx*ym*2.0D0-t3*t5*t18*t19*&
              &xx*ym*2.0D0-t5*t12*t18*t19*xx*ym*2.0D0-t5*t18*t19*t20*xx*ym*6.0D0-&
              &t3*t5*t18*t19*xm*yy*4.0D0+t2*t5*t18*t19*xx*yy*2.0D0+t3*t5*t18*t19*&
              &xx*yy*2.0D0+t5*t12*t18*t19*xx*yy*6.0D0+t5*t18*t19*t20*xx*yy*2.0D0+&
              &D*t7*t12**2*t21*t22*xx+D*t7*t20**2*t21*t22*xx-D*a*t2*t3*t5*t7*1.1D&
              &1-D*a*t2*t5*t7*t12-D*a*t3*t5*t7*t12*5.0D0-D*a*t2*t5*t7*t20-D*a*t3*&
              &t5*t7*t20*5.0D0-D*t7*t8*t21*t22*xm*1.0D1+D*t7*t8*t21*t22*xx*6.0D0+&
              &D*t7*t21*t22*t23*xx+D*t5*t8*t18*t21*ym*4.0D0-D*t5*t8*t18*t21*yy*4.&
              &0D0+D*a*t2*t5*t7*xm*xx*5.0D0+D*a*t3*t5*t7*xm*xx*1.2D1+D*a*t5*t7*t1&
              &2*xm*xx*3.0D0+D*a*t5*t7*t20*xm*xx*3.0D0-D*a*t2*t18*t22*xm*ym-D*a*t&
              &3*t18*t22*xm*ym*4.0D0-D*a*t12*t18*t22*xm*ym-D*a*t18*t20*t22*xm*ym*&
              &3.0D0+D*a*t2*t18*t22*xx*ym*4.0D0+D*a*t3*t18*t22*xx*ym*2.0D0+D*a*t1&
              &2*t18*t22*xx*ym*2.0D0+D*a*t18*t20*t22*xx*ym*6.0D0+D*a*t2*t18*t22*x&
              &m*yy+D*a*t3*t18*t22*xm*yy*4.0D0+D*a*t12*t18*t22*xm*yy*3.0D0+D*a*t1&
              &8*t20*t22*xm*yy-D*a*t2*t18*t22*xx*yy*4.0D0-D*a*t3*t18*t22*xx*yy*2.&
              &0D0-D*a*t12*t18*t22*xx*yy*6.0D0-D*a*t18*t20*t22*xx*yy*2.0D0+D*a*t2&
              &*t5*t7*ym*yy*2.0D0+D*a*t3*t5*t7*ym*yy*1.0D1-D*t2*t3*t7*t21*t22*xm*&
              &4.0D0-D*t3*t7*t12*t21*t22*xm*4.0D0-D*t3*t7*t20*t21*t22*xm*4.0D0+D*&
              &t2*t3*t7*t21*t22*xx*9.0D0+D*t2*t7*t12*t21*t22*xx*2.0D0+D*t3*t7*t12&
              &*t21*t22*xx*5.0D0+D*t2*t7*t20*t21*t22*xx*2.0D0+D*t3*t7*t20*t21*t22&
              &*xx*5.0D0+D*t7*t12*t20*t21*t22*xx*6.0D0+D*t2*t3*t5*t18*t21*ym*6.0D&
              &0+D*t3*t5*t12*t18*t21*ym*2.0D0+D*t3*t5*t18*t20*t21*ym*6.0D0-D*t2*t&
              &3*t5*t18*t21*yy*6.0D0-D*t3*t5*t12*t18*t21*yy*6.0D0-D*t3*t5*t18*t20&
              &*t21*yy*2.0D0+a*t5*t8*t22*t24*t25*xm-a*t5*t8*t22*t24*t27*xm-a*t5*t&
              &8*t22*t24*t25*xx+a*t5*t8*t22*t24*t27*xx-a*t7*t8*t18*t24*t26*ym+a*t&
              &7*t8*t18*t24*t28*ym+a*t7*t8*t18*t24*t26*yy-a*t7*t8*t18*t24*t28*yy+&
              &t3*t5*t7*t18*t19*t22*xm*ym*2.0D0-t2*t5*t7*t18*t19*t22*xx*ym-t3*t5*&
              &t7*t18*t19*t22*xx*ym-t5*t7*t12*t18*t19*t22*xx*ym-t5*t7*t18*t19*t20&
              &*t22*xx*ym*3.0D0-t3*t5*t7*t18*t19*t22*xm*yy*2.0D0+t2*t5*t7*t18*t19&
              &*t22*xx*yy+t3*t5*t7*t18*t19*t22*xx*yy+t5*t7*t12*t18*t19*t22*xx*yy*&
              &3.0D0+t5*t7*t18*t19*t20*t22*xx*yy-D*a*t5*t7*xm*xx*ym*yy*6.0D0-D*t2&
              &*t5*t18*t21*xm*xx*ym*2.0D0-D*t3*t5*t18*t21*xm*xx*ym*8.0D0-D*t5*t12&
              &*t18*t21*xm*xx*ym*2.0D0-D*t5*t18*t20*t21*xm*xx*ym*6.0D0+D*t2*t5*t1&
              &8*t21*xm*xx*yy*2.0D0+D*t3*t5*t18*t21*xm*xx*yy*8.0D0+D*t5*t12*t18*t&
              &21*xm*xx*yy*6.0D0+D*t5*t18*t20*t21*xm*xx*yy*2.0D0+D*t3*t7*t21*t22*&
              &xm*ym*yy*8.0D0-D*t2*t7*t21*t22*xx*ym*yy*4.0D0-D*t3*t7*t21*t22*xx*y&
              &m*yy*1.0D1-D*t7*t12*t21*t22*xx*ym*yy*4.0D0-D*t7*t20*t21*t22*xx*ym*&
              &yy*4.0D0))/xx
             
									             

              t2 = xm-xx
              t3 = 1.0D0/xx**2
              t4 = ym-yy
              t5 = a*xx
              t6 = a*yy
              t7 = cos(t5)
              t8 = cos(t6)
              t9 = sin(t6)
              t10 = 1.0D0/Mref
              t11 = sin(t5)
              t12 = t7**2
              t13 = t9*t11*2.0D0
              t14 = t13+4.0D0
              t15 = 1.0D0/xx
              t16 = t2**2
              t17 = t3*t16
              t18 = t4**2
              t19 = t3*t18
              t20 = t17+t19+1.0D0
              t21 = 1.0D0/sqrt(t20)
              t22 = t8**2
              t23 = t7*t9
              t24 = t23-(t12*t22)/2.0D0+2.0D1
              t25 = t8*t11
              t26 = t25-1.0D1
              t27 = 1.0D0/xx**3
              t28 = t9*t11
              t29 = t28+2.0D0
              t30 = a*t8*t9*t12
              t31 = a*t7*t11*t22
              t34 = a*t7*t9*t29
              t32 = t31-t34
              t33 = a**2
              t43 = a*t8*t11*t29
              t35 = t30-t43
              t36 = xm*2.0D0
              t44 = xx*2.0D0
              t37 = t36-t44
              t38 = t3*t37
              t39 = t16*t27*2.0D0
              t40 = t18*t27*2.0D0
              t41 = t38+t39+t40
              t42 = 1.0D0/t20**(3.0D0/2.0D0)
              t45 = t7*t8*t29*t33
              t46 = t2*t15*t21*t32
              t49 = t4*t15*t21*t35
              t47 = t46-t49
              t48 = t7*t8*t9*t11*t33*3.0D0
              t50 = t9**2
              t51 = t12*t33*t50
              t52 = t11**2
              t53 = t22*t33*t52
              t54 = t51+t53-t12*t22*t33-t9*t11*t29*t33
              t55 = t45+t48
              t56 = ym*2.0D0
              t58 = yy*2.0D0
              t57 = t56-t58
              f(i,2) =  -t15*((t4*t12*t22*t29*t41*t42)/2.0D0+a*t4*t7*t9*t12*t21*t22-a&
              &*t4*t7*t11*t21*t22*t29*2.0D0)+mu*(t45+t48-t15*(t30-t43-xx*(t45+t48&
              &+t3*t4*t21*t47+t4*t15*t21*(t15*t21*t32+t2*t3*t21*t32-t3*t4*t21*t35&
              &+t2*t15*t21*t54-t4*t15*t21*t55-(t2*t15*t32*t41*t42)/2.0D0+(t4*t15*&
              &t35*t41*t42)/2.0D0)-(t4*t15*t41*t42*t47)/2.0D0)+t4*t15*t21*t47)+t2&
              &*t15*t21*(t15*t21*t35-t2*t15*t21*t55+t4*t15*t21*t54+(t2*t27*t32*t4&
              &2*t57)/2.0D0-(t4*t27*t35*t42*t57)/2.0D0)+(t2*t27*t42*t47*t57)/2.0D&
              &0)+Mref*(t4*t15*t21*((t10*t14*(a*t9*t11-a*t7*t11*t22))/3.0D0+(a*t7&
              &*t8*t10*t14)/3.0D0-a*t7*t9*t10*t24*(2.0D0/3.0D0)+a*t7*t9*t10*t26*(&
              &2.0D0/3.0D0))+t2*t15*t21*((t10*t14*(t30+a*t7*t8))/3.0D0+(a*t9*t10*&
              &t11*t14)/3.0D0+a*t8*t10*t11*t24*(2.0D0/3.0D0)-a*t8*t10*t11*t26*(2.&
              &0D0/3.0D0)))+(t2*t12*t22*t27*t29*t42*t57)/2.0D0+a*t2*t8*t11*t12*t1&
              &5*t21*t22-a*t2*t8*t9*t12*t15*t21*t29*2.0D0
									              
									              
              t2 = 1.0D0/Mref
              t5 = a*xx
              t3 = cos(t5)
              t6 = a*yy
              t4 = cos(t6)
              t7 = ym**2
              t8 = t4**2
              t9 = a**2
              t10 = yy**2
              t11 = xm**2
              t12 = xx**2
              t13 = sin(t5)
              t14 = sin(t6)
              t15 = t3**2
              t16 = -epn-1.0D0
              t17 = 3.0D0**t16
              t18 = t3*t14*2.0D0
              t27 = t8*t15
              t19 = t18-t27+4.0D1
              t20 = t2*t19
              t21 = t20**epn
              t22 = ym-yy
              t23 = t12*2.0D0
              t26 = ym*yy*2.0D0
              t28 = xm*xx*2.0D0
              t24 = t7+t10+t11+t23-t26-t28
              t25 = 1.0D0/t24**2
              t29 = t3*t4*xm
              t30 = t13*t14*ym
              t31 = t4*t14*t15*xm
              t32 = t3*t8*t13*yy
              t35 = t3*t4*xx
              t36 = t13*t14*yy
              t37 = t3*t8*t13*ym
              t38 = t4*t14*t15*xx
              t33 = t29+t30+t31+t32-t35-t36-t37-t38
              t34 = xm-xx
              t39 = 1.0D0/Mref**2
              t40 = epn-1.0D0
              t41 = t20**t40
              t42 = 1.0D0/t24
              t43 = 1.0D0/xx
              t44 = t14**2
              t45 = t13**2
              t46 = 1.0D0/xx**2
              t47 = t24*t46
              t48 = 1.0D0/sqrt(t47)
              t49 = t13*t14
              t50 = t49+2.0D0
              t51 = t3*t14*5.0D0
              t52 = -t27+t51+1.0D2
              t53 = t7+t10+t11-t26-xm*xx
              t54 = t11**2
              t55 = t12**2
              t56 = a*xx*2.0D0
              t57 = sin(t56)
              t58 = a*yy*2.0D0
              t59 = sin(t58)
              t60 = t7**2
              t61 = t10**2
              t62 = t4*t13*2.0D0
              t63 = 1.0D0/xx**3
              t64 = ym*2.0D0
              t65 = t64-yy*2.0D0
              t66 = 1.0D0/t47**(3.0D0/2.0D0)
              f(i,3) =  -kpari*(t43*(t2*t17*t21*t22*t25*xx*(a*t3*t4*t7*2.0D0+a*t3*t4*&
              &t10*2.0D0-t7*t8*t9*ym*2.0D0-t8*t9*t10*ym*6.0D0-t8*t9*t11*ym*2.0D0-&
              &t8*t9*t12*ym*4.0D0+t7*t8*t9*yy*6.0D0+t8*t9*t10*yy*2.0D0+t8*t9*t11*&
              &yy*2.0D0+t8*t9*t12*yy*4.0D0+a*t4*t7*t14*t15*2.0D0+a*t4*t10*t14*t15&
              &*2.0D0+a*t3*t4*xm*xx*2.0D0-a*t13*t14*xm*ym*2.0D0+a*t13*t14*xx*ym*4&
              &.0D0+a*t13*t14*xm*yy*2.0D0-a*t13*t14*xx*yy*4.0D0-a*t3*t4*ym*yy*4.0&
              &D0+t4*t7*t9*t13*xm*2.0D0+t4*t9*t10*t13*xm*2.0D0+t4*t9*t11*t13*xm*2&
              &.0D0+t4*t9*t12*t13*xm*8.0D0-t4*t7*t9*t13*xx*2.0D0-t4*t9*t10*t13*xx&
              &*2.0D0-t4*t9*t11*t13*xx*6.0D0-t4*t9*t12*t13*xx*4.0D0-t3*t7*t9*t14*&
              &ym*2.0D0-t3*t9*t10*t14*ym*6.0D0-t3*t9*t11*t14*ym*2.0D0-t3*t9*t12*t&
              &14*ym*4.0D0+t7*t8*t9*t15*ym*4.0D0+t8*t9*t10*t15*ym*1.2D1+t8*t9*t11&
              &*t15*ym*4.0D0+t8*t9*t12*t15*ym*8.0D0+t3*t7*t9*t14*yy*6.0D0+t3*t9*t&
              &10*t14*yy*2.0D0+t3*t9*t11*t14*yy*2.0D0+t3*t9*t12*t14*yy*4.0D0-t7*t&
              &8*t9*t15*yy*1.2D1-t8*t9*t10*t15*yy*4.0D0-t8*t9*t11*t15*yy*4.0D0-t8&
              &*t9*t12*t15*yy*8.0D0+t8*t9*xm*xx*ym*4.0D0-t8*t9*xm*xx*yy*4.0D0+a*t&
              &4*t14*t15*xm*xx*2.0D0+a*t3*t8*t13*xm*ym*2.0D0-a*t3*t8*t13*xx*ym*4.&
              &0D0-a*t3*t8*t13*xm*yy*2.0D0+a*t3*t8*t13*xx*yy*4.0D0-a*t4*t14*t15*y&
              &m*yy*4.0D0+t3*t9*t14*xm*xx*ym*4.0D0-t8*t9*t15*xm*xx*ym*8.0D0-t3*t9&
              &*t14*xm*xx*yy*4.0D0+t8*t9*t15*xm*xx*yy*8.0D0-t4*t9*t13*xm*ym*yy*4.&
              &0D0+t4*t9*t13*xx*ym*yy*4.0D0+t3*t4*t7*t9*t13*t14*xm*4.0D0+t3*t4*t9&
              &*t10*t13*t14*xm*4.0D0+t3*t4*t9*t11*t13*t14*xm*4.0D0+t3*t4*t9*t12*t&
              &13*t14*xm*1.6D1-t3*t4*t7*t9*t13*t14*xx*4.0D0-t3*t4*t9*t10*t13*t14*&
              &xx*4.0D0-t3*t4*t9*t11*t13*t14*xx*1.2D1-t3*t4*t9*t12*t13*t14*xx*8.0&
              &D0-t3*t4*t9*t13*t14*xm*ym*yy*8.0D0+t3*t4*t9*t13*t14*xx*ym*yy*8.0D0&
              &)-a*t2*t17*t21*t22*t25*t33*t53*2.0D0+epn*t9*t13*t17*t22*t33*t39*t4&
              &1*t42*xx*(t14-t3*t8)*4.0D0)-t2*t17*t21*t25*t34*(a*t11*t13*t14*2.0D&
              &0+a*t12*t13*t14*4.0D0+t7*t9*t15*xm*2.0D0+t9*t10*t15*xm*2.0D0+t9*t1&
              &1*t15*xm*2.0D0+t9*t12*t15*xm*8.0D0-t7*t9*t15*xx*2.0D0-t9*t10*t15*x&
              &x*2.0D0-t9*t11*t15*xx*6.0D0-t9*t12*t15*xx*4.0D0-a*t3*t8*t11*t13*2.&
              &0D0-a*t3*t8*t12*t13*4.0D0-a*t13*t14*xm*xx*4.0D0-a*t3*t4*xm*ym*2.0D&
              &0+a*t3*t4*xx*ym*2.0D0+a*t3*t4*xm*yy*2.0D0-a*t3*t4*xx*yy*2.0D0+t3*t&
              &7*t9*t14*xm*2.0D0+t3*t9*t10*t14*xm*2.0D0+t3*t9*t11*t14*xm*2.0D0+t3&
              &*t9*t12*t14*xm*8.0D0-t7*t8*t9*t15*xm*4.0D0-t8*t9*t10*t15*xm*4.0D0-&
              &t8*t9*t11*t15*xm*4.0D0-t8*t9*t12*t15*xm*1.6D1-t3*t7*t9*t14*xx*2.0D&
              &0-t3*t9*t10*t14*xx*2.0D0-t3*t9*t11*t14*xx*6.0D0-t3*t9*t12*t14*xx*4&
              &.0D0+t7*t8*t9*t15*xx*4.0D0+t8*t9*t10*t15*xx*4.0D0+t8*t9*t11*t15*xx&
              &*1.2D1+t8*t9*t12*t15*xx*8.0D0-t4*t7*t9*t13*ym*2.0D0-t4*t9*t10*t13*&
              &ym*6.0D0-t4*t9*t11*t13*ym*2.0D0-t4*t9*t12*t13*ym*4.0D0+t4*t7*t9*t1&
              &3*yy*6.0D0+t4*t9*t10*t13*yy*2.0D0+t4*t9*t11*t13*yy*2.0D0+t4*t9*t12&
              &*t13*yy*4.0D0-t9*t15*xm*ym*yy*4.0D0+t9*t15*xx*ym*yy*4.0D0+a*t3*t8*&
              &t13*xm*xx*4.0D0-a*t4*t14*t15*xm*ym*2.0D0+a*t4*t14*t15*xx*ym*2.0D0+&
              &a*t4*t14*t15*xm*yy*2.0D0-a*t4*t14*t15*xx*yy*2.0D0+t4*t9*t13*xm*xx*&
              &ym*4.0D0-t4*t9*t13*xm*xx*yy*4.0D0-t3*t9*t14*xm*ym*yy*4.0D0+t8*t9*t&
              &15*xm*ym*yy*8.0D0+t3*t9*t14*xx*ym*yy*4.0D0-t8*t9*t15*xx*ym*yy*8.0D&
              &0-t3*t4*t7*t9*t13*t14*ym*4.0D0-t3*t4*t9*t10*t13*t14*ym*1.2D1-t3*t4&
              &*t9*t11*t13*t14*ym*4.0D0-t3*t4*t9*t12*t13*t14*ym*8.0D0+t3*t4*t7*t9&
              &*t13*t14*yy*1.2D1+t3*t4*t9*t10*t13*t14*yy*4.0D0+t3*t4*t9*t11*t13*t&
              &14*yy*4.0D0+t3*t4*t9*t12*t13*t14*yy*8.0D0+t3*t4*t9*t13*t14*xm*xx*y&
              &m*8.0D0-t3*t4*t9*t13*t14*xm*xx*yy*8.0D0)+a*t2*t17*t21*t25*t33*t34*&
              &t65+epn*t3*t4*t9*t17*t33*t34*t39*t41*t42*(t3*t14+1.0D0)*4.0D0)-t43&
              &*((a*t3*t4*t22*t48*(t3*t14*1.0D2-t13*t14*1.0D1+t15*t44*5.0D0-t44*t&
              &45*5.0D0+t3*t8*t13*4.0D0-t3*t8*t14*t15+t3*t8*t14*t45*2.0D0))/3.0D0&
              &-(a*t4*t13*t22*t48*t50*t52)/3.0D0+(t3*t4*t22*t50*t52*t53*t63*t66)/&
              &3.0D0)+a*csii*t25*t43*(t54+t55*6.0D0+t7*t11+t7*t12*5.0D0+t10*t11+t&
              &10*t12*5.0D0+t11*t12*1.1D1-t8*t54-t8*t55*6.0D0-t15*t54*2.0D0-t15*t&
              &55*1.2D1-t7*t8*t11-t7*t8*t12*5.0D0-t8*t10*t11-t8*t10*t12*5.0D0-t8*&
              &t11*t12*1.1D1-t7*t11*t15*2.0D0-t7*t12*t15*1.0D1-t10*t11*t15*2.0D0-&
              &t10*t12*t15*1.0D1-t11*t12*t15*2.2D1-t3*t14*t54*2.0D1-t3*t14*t55*1.&
              &2D2+t8*t15*t54*2.0D0+t8*t15*t55*1.2D1+t13*t14*t54*2.0D0+t13*t14*t5&
              &5*1.2D1-t7*xm*xx*3.0D0-t10*xm*xx*3.0D0-t11*xm*xx*5.0D0-t12*xm*xx*1&
              &.2D1-t11*ym*yy*2.0D0-t12*ym*yy*1.0D1-a*t55*t57*xm*1.8D1+a*t54*t57*&
              &xx*2.0D0+a*t55*t57*xx*1.0D1+a*t57*t60*xx+a*t57*t61*xx-a*t55*t59*ym&
              &*4.0D0+a*t55*t59*yy*4.0D0-t3*t7*t11*t14*2.0D1-t3*t7*t12*t14*1.0D2-&
              &t3*t10*t11*t14*2.0D1-t3*t10*t12*t14*1.0D2-t3*t11*t12*t14*2.2D2+t7*&
              &t8*t11*t15*2.0D0+t7*t8*t12*t15*1.0D1+t8*t10*t11*t15*2.0D0+t7*t11*t&
              &13*t14*2.0D0+t8*t10*t12*t15*1.0D1+t7*t12*t13*t14*1.0D1+t8*t11*t12*&
              &t15*2.2D1+t10*t11*t13*t14*2.0D0+t10*t12*t13*t14*1.0D1+t11*t12*t13*&
              &t14*2.2D1+t7*t8*xm*xx*3.0D0+t8*t10*xm*xx*3.0D0+t8*t11*xm*xx*5.0D0+&
              &t8*t12*xm*xx*1.2D1+t7*t15*xm*xx*6.0D0+t10*t15*xm*xx*6.0D0+t11*t15*&
              &xm*xx*1.0D1+t12*t15*xm*xx*2.4D1+t8*t11*ym*yy*2.0D0+t8*t12*ym*yy*1.&
              &0D1+t11*t15*ym*yy*4.0D0+t12*t15*ym*yy*2.0D1+xm*xx*ym*yy*6.0D0-a*t3&
              &*t14*t55*xm*2.0D1-a*t7*t12*t57*xm*6.0D0-a*t10*t12*t57*xm*6.0D0-a*t&
              &11*t12*t57*xm*8.0D0-a*t13*t14*t55*xm*2.0D2+a*t3*t14*t54*xx*2.0D0+a&
              &*t3*t14*t55*xx*1.2D1+a*t7*t10*t57*xx*6.0D0+a*t7*t11*t57*xx*3.0D0+a&
              &*t7*t12*t57*xx*7.0D0+a*t3*t14*t60*xx*2.0D0+a*t3*t14*t61*xx*2.0D0+a&
              &*t10*t11*t57*xx*3.0D0+a*t10*t12*t57*xx*7.0D0+a*t11*t12*t57*xx*1.7D&
              &1+a*t13*t14*t54*xx*2.0D1+a*t13*t14*t55*xx*1.2D2+a*t13*t14*t60*xx*2&
              &.0D1+a*t13*t14*t61*xx*2.0D1+a*t3*t4*t55*ym*8.0D1-a*t4*t13*t55*ym*8&
              &.0D0-a*t7*t12*t59*ym*2.0D0-a*t10*t12*t59*ym*6.0D0-a*t11*t12*t59*ym&
              &*6.0D0-a*t3*t4*t55*yy*8.0D1+a*t4*t13*t55*yy*8.0D0+a*t7*t12*t59*yy*&
              &6.0D0+a*t10*t12*t59*yy*2.0D0+a*t11*t12*t59*yy*6.0D0+t3*t7*t14*xm*x&
              &x*6.0D1+t3*t10*t14*xm*xx*6.0D1+t3*t11*t14*xm*xx*1.0D2+t3*t12*t14*x&
              &m*xx*2.4D2-t7*t8*t15*xm*xx*6.0D0-t8*t10*t15*xm*xx*6.0D0-t7*t13*t14&
              &*xm*xx*6.0D0-t8*t11*t15*xm*xx*1.0D1-t8*t12*t15*xm*xx*2.4D1-t10*t13&
              &*t14*xm*xx*6.0D0-t11*t13*t14*xm*xx*1.0D1-t12*t13*t14*xm*xx*2.4D1-t&
              &3*t4*t7*xm*ym*2.0D0-t3*t4*t10*xm*ym*6.0D0-t3*t4*t11*xm*ym*2.0D0-t3&
              &*t4*t12*xm*ym*8.0D0-t4*t7*t13*xm*ym*2.0D1-t4*t10*t13*xm*ym*6.0D1-t&
              &4*t11*t13*xm*ym*2.0D1-t4*t12*t13*xm*ym*8.0D1+t3*t4*t7*xx*ym*4.0D0+&
              &t3*t4*t10*xx*ym*1.2D1+t3*t4*t11*xx*ym*8.0D0+t3*t4*t12*xx*ym*4.0D0+&
              &t4*t7*t13*xx*ym*4.0D1+t4*t10*t13*xx*ym*1.2D2+t4*t11*t13*xx*ym*8.0D&
              &1+t4*t12*t13*xx*ym*4.0D1+t3*t4*t7*xm*yy*6.0D0+t3*t4*t10*xm*yy*2.0D&
              &0+t3*t4*t11*xm*yy*2.0D0+t3*t4*t12*xm*yy*8.0D0+t4*t7*t13*xm*yy*6.0D&
              &1+t4*t10*t13*xm*yy*2.0D1+t4*t11*t13*xm*yy*2.0D1+t4*t12*t13*xm*yy*8&
              &.0D1-t3*t4*t7*xx*yy*1.2D1-t3*t4*t10*xx*yy*4.0D0-t3*t4*t11*xx*yy*8.&
              &0D0-t3*t4*t12*xx*yy*4.0D0-t4*t7*t13*xx*yy*1.2D2-t4*t10*t13*xx*yy*4&
              &.0D1-t4*t11*t13*xx*yy*8.0D1-t4*t12*t13*xx*yy*4.0D1+t3*t11*t14*ym*y&
              &y*4.0D1+t3*t12*t14*ym*yy*2.0D2-t8*t11*t15*ym*yy*4.0D0-t8*t12*t15*y&
              &m*yy*2.0D1-t11*t13*t14*ym*yy*4.0D0-t12*t13*t14*ym*yy*2.0D1-t8*xm*x&
              &x*ym*yy*6.0D0-t15*xm*xx*ym*yy*1.2D1-a*t3*t7*t12*t14*xm*8.0D0-a*t3*&
              &t10*t12*t14*xm*8.0D0-a*t3*t11*t12*t14*xm*8.0D0-a*t7*t12*t13*t14*xm&
              &*8.0D1-a*t10*t12*t13*t14*xm*8.0D1-a*t11*t12*t13*t14*xm*8.0D1+a*t3*&
              &t8*t13*t55*xm*4.0D1+a*t3*t7*t10*t14*xx*1.2D1+a*t3*t7*t11*t14*xx*4.&
              &0D0+a*t3*t7*t12*t14*xx*1.0D1+a*t3*t10*t11*t14*xx*4.0D0+a*t3*t10*t1&
              &2*t14*xx*1.0D1+a*t3*t11*t12*t14*xx*1.8D1+a*t7*t10*t13*t14*xx*1.2D2&
              &+a*t7*t11*t13*t14*xx*4.0D1+a*t7*t12*t13*t14*xx*1.0D2+a*t10*t11*t13&
              &*t14*xx*4.0D1+a*t10*t12*t13*t14*xx*1.0D2+a*t11*t12*t13*t14*xx*1.8D&
              &2-a*t3*t8*t13*t54*xx*4.0D0-a*t3*t8*t13*t55*xx*2.4D1-a*t3*t8*t13*t6&
              &0*xx*4.0D0-a*t3*t8*t13*t61*xx*4.0D0+a*t3*t4*t7*t12*ym*4.0D1+a*t3*t&
              &4*t10*t12*ym*1.2D2+a*t3*t4*t11*t12*ym*1.2D2-a*t4*t7*t12*t13*ym*4.0&
              &D0-a*t4*t10*t12*t13*ym*1.2D1-a*t4*t11*t12*t13*ym*1.2D1+a*t4*t14*t1&
              &5*t55*ym*1.6D1-a*t3*t4*t7*t12*yy*1.2D2-a*t3*t4*t10*t12*yy*4.0D1-a*&
              &t3*t4*t11*t12*yy*1.2D2+a*t4*t7*t12*t13*yy*1.2D1+a*t4*t10*t12*t13*y&
              &y*4.0D0+a*t4*t11*t12*t13*yy*1.2D1-a*t4*t14*t15*t55*yy*1.6D1+a*t7*t&
              &59*xm*xx*ym*2.0D0+a*t10*t59*xm*xx*ym*6.0D0+a*t11*t59*xm*xx*ym*2.0D&
              &0+a*t12*t59*xm*xx*ym*8.0D0-a*t7*t59*xm*xx*yy*6.0D0-a*t10*t59*xm*xx&
              &*yy*2.0D0-a*t11*t59*xm*xx*yy*2.0D0-a*t12*t59*xm*xx*yy*8.0D0+a*t12*&
              &t57*xm*ym*yy*1.2D1-a*t7*t57*xx*ym*yy*4.0D0-a*t10*t57*xx*ym*yy*4.0D&
              &0-a*t11*t57*xx*ym*yy*6.0D0-a*t12*t57*xx*ym*yy*1.4D1-t3*t14*xm*xx*y&
              &m*yy*1.2D2+t8*t15*xm*xx*ym*yy*1.2D1+t13*t14*xm*xx*ym*yy*1.2D1+a*t3&
              &*t7*t8*t12*t13*xm*1.6D1+a*t3*t8*t10*t12*t13*xm*1.6D1+a*t3*t8*t11*t&
              &12*t13*xm*1.6D1-a*t3*t7*t8*t10*t13*xx*2.4D1-a*t3*t7*t8*t11*t13*xx*&
              &8.0D0-a*t3*t7*t8*t12*t13*xx*2.0D1-a*t3*t8*t10*t11*t13*xx*8.0D0-a*t&
              &3*t8*t10*t12*t13*xx*2.0D1-a*t3*t8*t11*t12*t13*xx*3.6D1+a*t4*t7*t12&
              &*t14*t15*ym*8.0D0+a*t4*t10*t12*t14*t15*ym*2.4D1+a*t4*t11*t12*t14*t&
              &15*ym*2.4D1-a*t4*t7*t12*t14*t15*yy*2.4D1-a*t4*t10*t12*t14*t15*yy*8&
              &.0D0-a*t4*t11*t12*t14*t15*yy*2.4D1-a*t3*t4*t7*xm*xx*ym*4.0D1-a*t3*&
              &t4*t10*xm*xx*ym*1.2D2-a*t3*t4*t11*xm*xx*ym*4.0D1-a*t3*t4*t12*xm*xx&
              &*ym*1.6D2+a*t4*t7*t13*xm*xx*ym*4.0D0+a*t4*t10*t13*xm*xx*ym*1.2D1+a&
              &*t4*t11*t13*xm*xx*ym*4.0D0+a*t4*t12*t13*xm*xx*ym*1.6D1+a*t3*t4*t7*&
              &xm*xx*yy*1.2D2+a*t3*t4*t10*xm*xx*yy*4.0D1+a*t3*t4*t11*xm*xx*yy*4.0&
              &D1+a*t3*t4*t12*xm*xx*yy*1.6D2-a*t4*t7*t13*xm*xx*yy*1.2D1-a*t4*t10*&
              &t13*xm*xx*yy*4.0D0-a*t4*t11*t13*xm*xx*yy*4.0D0-a*t4*t12*t13*xm*xx*&
              &yy*1.6D1+a*t3*t12*t14*xm*ym*yy*1.6D1+a*t12*t13*t14*xm*ym*yy*1.6D2-&
              &a*t3*t7*t14*xx*ym*yy*8.0D0-a*t3*t10*t14*xx*ym*yy*8.0D0-a*t3*t11*t1&
              &4*xx*ym*yy*8.0D0-a*t3*t12*t14*xx*ym*yy*2.0D1-a*t7*t13*t14*xx*ym*yy&
              &*8.0D1-a*t10*t13*t14*xx*ym*yy*8.0D1-a*t11*t13*t14*xx*ym*yy*8.0D1-a&
              &*t12*t13*t14*xx*ym*yy*2.0D2-t3*t4*t7*t13*t14*xm*ym*2.0D0-t3*t4*t10&
              &*t13*t14*xm*ym*6.0D0-t3*t4*t11*t13*t14*xm*ym*2.0D0-t3*t4*t12*t13*t&
              &14*xm*ym*8.0D0+t3*t4*t7*t13*t14*xx*ym*4.0D0+t3*t4*t10*t13*t14*xx*y&
              &m*1.2D1+t3*t4*t11*t13*t14*xx*ym*8.0D0+t3*t4*t12*t13*t14*xx*ym*4.0D&
              &0+t3*t4*t7*t13*t14*xm*yy*6.0D0+t3*t4*t10*t13*t14*xm*yy*2.0D0+t3*t4&
              &*t11*t13*t14*xm*yy*2.0D0+t3*t4*t12*t13*t14*xm*yy*8.0D0-t3*t4*t7*t1&
              &3*t14*xx*yy*1.2D1-t3*t4*t10*t13*t14*xx*yy*4.0D0-t3*t4*t11*t13*t14*&
              &xx*yy*8.0D0-t3*t4*t12*t13*t14*xx*yy*4.0D0-a*t4*t7*t14*t15*xm*xx*ym&
              &*8.0D0-a*t4*t10*t14*t15*xm*xx*ym*2.4D1-a*t4*t11*t14*t15*xm*xx*ym*8&
              &.0D0-a*t4*t12*t14*t15*xm*xx*ym*3.2D1+a*t4*t7*t14*t15*xm*xx*yy*2.4D&
              &1+a*t4*t10*t14*t15*xm*xx*yy*8.0D0+a*t4*t11*t14*t15*xm*xx*yy*8.0D0+&
              &a*t4*t12*t14*t15*xm*xx*yy*3.2D1-a*t3*t8*t12*t13*xm*ym*yy*3.2D1+a*t&
              &3*t7*t8*t13*xx*ym*yy*1.6D1+a*t3*t8*t10*t13*xx*ym*yy*1.6D1+a*t3*t8*&
              &t11*t13*xx*ym*yy*1.6D1+a*t3*t8*t12*t13*xx*ym*yy*4.0D1)+Mref*pcourr&
              &*t3*t4*(a*t2*t3*t22*t43*t48*(t4-t14*5.0D0+t4*t13*t14)*(4.0D0/3.0D0&
              &)+a*t2*t13*t34*t43*t48*(t4*1.0D1+t13+t14*2.0D0-t8*t13*2.0D0)*(2.0D&
              &0/3.0D0))+(sqrt(3.0D0)*t50**2*1.0D0/sqrt(-t2*(t62-2.0D1))*(t18-t27&
              &+t62+2.0D1))/(tie*(t4*t13-1.0D1)*2.0D0)+(a*t3*t8*t34*t43*t48*(t3*1&
              &.0D1+t13*1.0D2+t13*t15*2.0D0+t14*t15*4.0D0+t3*t13*t14*1.0D1-t8*t13&
              &*t15*3.0D0))/3.0D0-(a*t3*t14*t34*t43*t48*t50*t52)/3.0D0+(t3*t4*t34&
              &*t50*t52*t63*t65*t66)/6.0D0									      

                     
               t5 = a*xx
              t2 = cos(t5)
              t3 = a*yy
              t4 = cos(t3)
              t6 = sin(t3)
              t7 = sin(t5)
              t8 = ym-yy
              t9 = 1.0D0/xx**2
              t10 = xm**2
              t11 = xx**2
              t12 = t11*2.0D0
              t13 = ym**2
              t14 = yy**2
              t22 = xm*xx*2.0D0
              t23 = ym*yy*2.0D0
              t15 = t10+t12+t13+t14-t22-t23
              t16 = t9*t15
              t17 = 1.0D0/sqrt(t16)
              t18 = t4*t7
              t19 = t18-1.0D1
              t20 = t6*t7
              t21 = t20+2.0D0
              t24 = 1.0D0/xx
              t25 = -epn+1.0D0
              t26 = 3.0D0**t25
              t27 = 1.0D0/Mref
              t28 = t4*t7*2.0D0
              t29 = t28-2.0D1
              t31 = t27*t29
              t30 = (-t31)**epn
              t32 = a**2
              t33 = -epn+2.0D0
              t34 = 3.0D0**t33
              t35 = t2**2
              t36 = 1.0D0/t29
              t37 = t4**2
              t38 = t7**2
              t39 = t6**2
              t40 = t4*t6*t7
              t41 = t4-t6*5.0D0+t40
              t42 = 1.0D0/t15**2
              t43 = t10**2
              t44 = t11**2
              t45 = t13**2
              t46 = t14**2
              t47 = xm-xx
              t48 = t4*1.0D1
              t49 = t6*2.0D0
              t50 = t7+t48+t49-t7*t37*2.0D0
              t51 = 1.0D0/xx**3
              t52 = 1.0D0/t16**(3.0D0/2.0D0)
              f(i,4) = t24*(a*t4*t8*t17*t35*t41*(1.0D1/3.0D0)-a*t4*t7*t8*t17*t19*t21&
              &*(5.0D0/3.0D0)+t2*t4*t8*t19*t21*t51*t52*(t10+t13+t14-t23-xm*xx)*(5&
              &.0D0/3.0D0))+(1.0D0/Mref**2*t24*t42*(Mref*a*kpare*t2*t4*t26*t30*t4&
              &5*4.0D0+Mref*a*kpare*t2*t4*t26*t30*t46*4.0D0+Mref*a*kpare*t2*t4*t1&
              &0*t13*t26*t30*4.0D0+Mref*a*kpare*t2*t4*t10*t14*t26*t30*4.0D0+Mref*&
              &a*kpare*t2*t4*t11*t13*t26*t30*4.0D0+Mref*a*kpare*t2*t4*t11*t14*t26&
              &*t30*4.0D0-Mref*a*kpare*t2*t4*t11*t15*t26*t30*2.0D0-Mref*a*kpare*t&
              &2*t4*t13*t15*t26*t30*2.0D0-Mref*a*kpare*t2*t4*t14*t15*t26*t30*2.0D&
              &0+Mref*a*kpare*t2*t4*t13*t14*t30*t34*8.0D0-Mref*a*kpare*t2*t4*t13*&
              &t26*t30*xm*xx*8.0D0-Mref*a*kpare*t2*t4*t14*t26*t30*xm*xx*8.0D0+Mre&
              &f*a*kpare*t2*t4*t15*t26*t30*xm*xx*2.0D0+Mref*a*kpare*t6*t7*t10*t26&
              &*t30*xm*ym*4.0D0+Mref*a*kpare*t6*t7*t13*t26*t30*xm*ym*4.0D0-Mref*a&
              &*kpare*t6*t7*t15*t26*t30*xm*ym*2.0D0+Mref*a*kpare*t6*t7*t11*t30*t3&
              &4*xm*ym*4.0D0+Mref*a*kpare*t6*t7*t14*t30*t34*xm*ym*4.0D0-Mref*a*kp&
              &are*t6*t7*t11*t26*t30*xx*ym*4.0D0-Mref*a*kpare*t6*t7*t13*t26*t30*x&
              &x*ym*4.0D0-Mref*a*kpare*t6*t7*t10*t30*t34*xx*ym*4.0D0-Mref*a*kpare&
              &*t6*t7*t14*t30*t34*xx*ym*4.0D0-Mref*a*kpare*t6*t7*t10*t26*t30*xm*y&
              &y*4.0D0-Mref*a*kpare*t6*t7*t14*t26*t30*xm*yy*4.0D0+Mref*a*kpare*t6&
              &*t7*t15*t26*t30*xm*yy*2.0D0-Mref*a*kpare*t6*t7*t11*t30*t34*xm*yy*4&
              &.0D0-Mref*a*kpare*t6*t7*t13*t30*t34*xm*yy*4.0D0+Mref*a*kpare*t6*t7&
              &*t11*t26*t30*xx*yy*4.0D0+Mref*a*kpare*t6*t7*t14*t26*t30*xx*yy*4.0D&
              &0+Mref*a*kpare*t6*t7*t10*t30*t34*xx*yy*4.0D0+Mref*a*kpare*t6*t7*t1&
              &3*t30*t34*xx*yy*4.0D0-Mref*a*kpare*t2*t4*t10*t26*t30*ym*yy*8.0D0-M&
              &ref*a*kpare*t2*t4*t11*t26*t30*ym*yy*8.0D0-Mref*a*kpare*t2*t4*t13*t&
              &26*t30*ym*yy*1.6D1-Mref*a*kpare*t2*t4*t14*t26*t30*ym*yy*1.6D1+Mref&
              &*a*kpare*t2*t4*t15*t26*t30*ym*yy*4.0D0+Mref*kpare*t4*t7*t11*t15*t2&
              &6*t30*t32*xm*4.0D0-Mref*kpare*t4*t7*t10*t15*t26*t30*t32*xx*2.0D0-M&
              &ref*kpare*t4*t7*t11*t15*t26*t30*t32*xx*2.0D0-Mref*kpare*t4*t7*t13*&
              &t15*t26*t30*t32*xx*2.0D0-Mref*kpare*t4*t7*t14*t15*t26*t30*t32*xx*2&
              &.0D0-Mref*kpare*t2*t6*t11*t15*t26*t30*t32*ym*4.0D0+Mref*kpare*t2*t&
              &6*t11*t15*t26*t30*t32*yy*4.0D0+Mref*a*kpare*t2*t4*t26*t30*xm*xx*ym&
              &*yy*1.6D1+Mref*kpare*t2*t6*t15*t26*t30*t32*xm*xx*ym*4.0D0-Mref*kpa&
              &re*t2*t6*t15*t26*t30*t32*xm*xx*yy*4.0D0+Mref*kpare*t4*t7*t15*t26*t&
              &30*t32*xx*ym*yy*4.0D0-Mref*epn*kpare*t11*t15*t26*t30*t32*t36*t38*t&
              &39*xm*8.0D0+Mref*epn*kpare*t13*t15*t26*t30*t32*t35*t36*t37*xx*4.0D&
              &0+Mref*epn*kpare*t14*t15*t26*t30*t32*t35*t36*t37*xx*4.0D0+Mref*epn&
              &*kpare*t10*t15*t26*t30*t32*t36*t38*t39*xx*4.0D0+Mref*epn*kpare*t11&
              &*t15*t26*t30*t32*t36*t38*t39*xx*4.0D0-Mref*epn*kpare*t15*t26*t30*t&
              &32*t35*t36*t37*xx*ym*yy*8.0D0-Mref*epn*kpare*t2*t4*t6*t7*t11*t15*t&
              &26*t30*t32*t36*ym*8.0D0+Mref*epn*kpare*t2*t4*t6*t7*t11*t15*t26*t30&
              &*t32*t36*yy*8.0D0+Mref*epn*kpare*t2*t4*t6*t7*t15*t26*t30*t32*t36*x&
              &m*xx*ym*8.0D0-Mref*epn*kpare*t2*t4*t6*t7*t15*t26*t30*t32*t36*xm*xx&
              &*yy*8.0D0))/9.0D0-Mref*pcourr*t2*t4*(a*t2*t8*t17*t24*t27*t41*(4.0D&
              &0/3.0D0)+a*t7*t17*t24*t27*t47*t50*(2.0D0/3.0D0))+a*csie*t24*t42*(t&
              &2*t4*t43*2.0D0+t2*t4*t44*1.2D1-t2*t6*t43*1.0D1-t2*t6*t44*6.0D1+t2*&
              &t4*t10*t11*2.2D1+t2*t4*t10*t13*2.0D0-t2*t6*t10*t11*1.1D2+t2*t4*t10&
              &*t14*2.0D0+t2*t4*t11*t13*1.0D1+t2*t4*t11*t14*1.0D1-t2*t6*t10*t13*1&
              &.0D1-t2*t6*t10*t14*1.0D1-t2*t6*t11*t13*5.0D1-t2*t6*t11*t14*5.0D1+a&
              &*t4*t7*t44*xm*2.0D1-a*t6*t7*t44*xm*1.0D2-a*t4*t7*t43*xx*2.0D0-a*t4&
              &*t7*t44*xx*1.2D1-a*t4*t7*t45*xx*2.0D0+a*t6*t7*t43*xx*1.0D1-a*t4*t7&
              &*t46*xx*2.0D0+a*t6*t7*t44*xx*6.0D1+a*t6*t7*t45*xx*1.0D1+a*t6*t7*t4&
              &6*xx*1.0D1+a*t2*t4*t44*ym*4.0D1+a*t2*t6*t44*ym*8.0D0-a*t2*t4*t44*y&
              &y*4.0D1-a*t2*t6*t44*yy*8.0D0+t2*t4*t6*t7*t43*2.0D0+t2*t4*t6*t7*t44&
              &*1.2D1-t2*t4*t10*xm*xx*1.0D1-t2*t4*t11*xm*xx*2.4D1+t2*t6*t10*xm*xx&
              &*5.0D1-t2*t4*t13*xm*xx*6.0D0+t2*t6*t11*xm*xx*1.2D2-t2*t4*t14*xm*xx&
              &*6.0D0+t2*t6*t13*xm*xx*3.0D1+t2*t6*t14*xm*xx*3.0D1-t4*t7*t10*xm*ym&
              &*1.0D1-t4*t7*t11*xm*ym*4.0D1-t6*t7*t10*xm*ym*2.0D0-t4*t7*t13*xm*ym&
              &*1.0D1-t6*t7*t11*xm*ym*8.0D0-t4*t7*t14*xm*ym*3.0D1-t6*t7*t13*xm*ym&
              &*2.0D0-t6*t7*t14*xm*ym*6.0D0+t10*t37*t38*xm*ym+t11*t37*t38*xm*ym*4&
              &.0D0-t10*t38*t39*xm*ym-t11*t38*t39*xm*ym*4.0D0+t13*t37*t38*xm*ym+t&
              &14*t37*t38*xm*ym*3.0D0-t13*t38*t39*xm*ym-t14*t38*t39*xm*ym*3.0D0+t&
              &4*t7*t10*xx*ym*4.0D1+t4*t7*t11*xx*ym*2.0D1+t6*t7*t10*xx*ym*8.0D0+t&
              &4*t7*t13*xx*ym*2.0D1+t6*t7*t11*xx*ym*4.0D0+t4*t7*t14*xx*ym*6.0D1+t&
              &6*t7*t13*xx*ym*4.0D0+t6*t7*t14*xx*ym*1.2D1-t10*t37*t38*xx*ym*4.0D0&
              &-t11*t37*t38*xx*ym*2.0D0+t10*t38*t39*xx*ym*4.0D0+t11*t38*t39*xx*ym&
              &*2.0D0-t13*t37*t38*xx*ym*2.0D0-t14*t37*t38*xx*ym*6.0D0+t13*t38*t39&
              &*xx*ym*2.0D0+t14*t38*t39*xx*ym*6.0D0+t4*t7*t10*xm*yy*1.0D1+t4*t7*t&
              &11*xm*yy*4.0D1+t6*t7*t10*xm*yy*2.0D0+t4*t7*t13*xm*yy*3.0D1+t6*t7*t&
              &11*xm*yy*8.0D0+t4*t7*t14*xm*yy*1.0D1+t6*t7*t13*xm*yy*6.0D0+t6*t7*t&
              &14*xm*yy*2.0D0-t10*t37*t38*xm*yy-t11*t37*t38*xm*yy*4.0D0+t10*t38*t&
              &39*xm*yy+t11*t38*t39*xm*yy*4.0D0-t13*t37*t38*xm*yy*3.0D0-t14*t37*t&
              &38*xm*yy+t13*t38*t39*xm*yy*3.0D0+t14*t38*t39*xm*yy-t4*t7*t10*xx*yy&
              &*4.0D1-t4*t7*t11*xx*yy*2.0D1-t6*t7*t10*xx*yy*8.0D0-t4*t7*t13*xx*yy&
              &*6.0D1-t6*t7*t11*xx*yy*4.0D0-t4*t7*t14*xx*yy*2.0D1-t6*t7*t13*xx*yy&
              &*1.2D1-t6*t7*t14*xx*yy*4.0D0+t10*t37*t38*xx*yy*4.0D0+t11*t37*t38*x&
              &x*yy*2.0D0-t10*t38*t39*xx*yy*4.0D0-t11*t38*t39*xx*yy*2.0D0+t13*t37&
              &*t38*xx*yy*6.0D0+t14*t37*t38*xx*yy*2.0D0-t13*t38*t39*xx*yy*6.0D0-t&
              &14*t38*t39*xx*yy*2.0D0-t2*t4*t10*ym*yy*4.0D0-t2*t4*t11*ym*yy*2.0D1&
              &+t2*t6*t10*ym*yy*2.0D1+t2*t6*t11*ym*yy*1.0D2+a*t4*t7*t10*t11*xm*8.&
              &0D0-a*t6*t7*t10*t11*xm*4.0D1+a*t4*t7*t11*t13*xm*8.0D0+a*t4*t7*t11*&
              &t14*xm*8.0D0-a*t6*t7*t11*t13*xm*4.0D1-a*t6*t7*t11*t14*xm*4.0D1-a*t&
              &4*t6*t35*t44*xm*1.6D1+a*t4*t6*t38*t44*xm*2.4D1-a*t4*t7*t10*t11*xx*&
              &1.8D1-a*t4*t7*t10*t13*xx*4.0D0+a*t6*t7*t10*t11*xx*9.0D1-a*t4*t7*t1&
              &0*t14*xx*4.0D0-a*t4*t7*t11*t13*xx*1.0D1-a*t4*t7*t11*t14*xx*1.0D1+a&
              &*t6*t7*t10*t13*xx*2.0D1+a*t6*t7*t10*t14*xx*2.0D1+a*t6*t7*t11*t13*x&
              &x*5.0D1-a*t4*t7*t13*t14*xx*1.2D1+a*t6*t7*t11*t14*xx*5.0D1+a*t6*t7*&
              &t13*t14*xx*6.0D1+a*t4*t6*t35*t43*xx*2.0D0+a*t4*t6*t35*t44*xx*8.0D0&
              &-a*t4*t6*t38*t43*xx*2.0D0-a*t4*t6*t38*t44*xx*1.6D1-a*t4*t6*t38*t45&
              &*xx*4.0D0-a*t4*t6*t38*t46*xx*4.0D0+a*t2*t4*t10*t11*ym*6.0D1+a*t2*t&
              &6*t10*t11*ym*1.2D1+a*t2*t4*t11*t13*ym*2.0D1+a*t2*t4*t11*t14*ym*6.0&
              &D1+a*t2*t6*t11*t13*ym*4.0D0+a*t2*t6*t11*t14*ym*1.2D1-a*t2*t7*t37*t&
              &44*ym*8.0D0+a*t2*t7*t39*t44*ym*8.0D0-a*t2*t4*t10*t11*yy*6.0D1-a*t2&
              &*t6*t10*t11*yy*1.2D1-a*t2*t4*t11*t13*yy*6.0D1-a*t2*t4*t11*t14*yy*2&
              &.0D1-a*t2*t6*t11*t13*yy*1.2D1-a*t2*t6*t11*t14*yy*4.0D0+a*t2*t7*t37&
              &*t44*yy*8.0D0-a*t2*t7*t39*t44*yy*8.0D0+t2*t4*t6*t7*t10*t11*2.2D1+t&
              &2*t4*t6*t7*t10*t13*2.0D0+t2*t4*t6*t7*t10*t14*2.0D0+t2*t4*t6*t7*t11&
              &*t13*1.0D1+t2*t4*t6*t7*t11*t14*1.0D1+t2*t4*xm*xx*ym*yy*1.2D1-t2*t6&
              &*xm*xx*ym*yy*6.0D1-a*t4*t6*t10*t11*t35*xm*8.0D0+a*t4*t6*t10*t11*t3&
              &8*xm*8.0D0-a*t4*t6*t11*t13*t35*xm*4.0D0-a*t4*t6*t11*t14*t35*xm*4.0&
              &D0+a*t4*t6*t11*t13*t38*xm*1.2D1+a*t4*t6*t11*t14*t38*xm*1.2D1+a*t4*&
              &t6*t10*t11*t35*xx*1.6D1+a*t4*t6*t10*t13*t35*xx*2.0D0-a*t4*t6*t10*t&
              &11*t38*xx*2.0D1+a*t4*t6*t10*t14*t35*xx*2.0D0+a*t4*t6*t11*t13*t35*x&
              &x*4.0D0+a*t4*t6*t11*t14*t35*xx*4.0D0-a*t4*t6*t10*t13*t38*xx*6.0D0-&
              &a*t4*t6*t10*t14*t38*xx*6.0D0-a*t4*t6*t11*t13*t38*xx*1.6D1-a*t4*t6*&
              &t11*t14*t38*xx*1.6D1-a*t4*t6*t13*t14*t38*xx*2.4D1-a*t2*t7*t10*t11*&
              &t37*ym*1.2D1+a*t2*t7*t10*t11*t39*ym*1.2D1-a*t2*t7*t11*t13*t37*ym*4&
              &.0D0-a*t2*t7*t11*t14*t37*ym*1.2D1+a*t2*t7*t11*t13*t39*ym*4.0D0+a*t&
              &2*t7*t11*t14*t39*ym*1.2D1+a*t2*t7*t10*t11*t37*yy*1.2D1-a*t2*t7*t10&
              &*t11*t39*yy*1.2D1+a*t2*t7*t11*t13*t37*yy*1.2D1+a*t2*t7*t11*t14*t37&
              &*yy*4.0D0-a*t2*t7*t11*t13*t39*yy*1.2D1-a*t2*t7*t11*t14*t39*yy*4.0D&
              &0-a*t2*t4*t10*xm*xx*ym*2.0D1-a*t2*t4*t11*xm*xx*ym*8.0D1-a*t2*t6*t1&
              &0*xm*xx*ym*4.0D0-a*t2*t4*t13*xm*xx*ym*2.0D1-a*t2*t6*t11*xm*xx*ym*1&
              &.6D1-a*t2*t4*t14*xm*xx*ym*6.0D1-a*t2*t6*t13*xm*xx*ym*4.0D0-a*t2*t6&
              &*t14*xm*xx*ym*1.2D1+a*t2*t4*t10*xm*xx*yy*2.0D1+a*t2*t4*t11*xm*xx*y&
              &y*8.0D1+a*t2*t6*t10*xm*xx*yy*4.0D0+a*t2*t4*t13*xm*xx*yy*6.0D1+a*t2&
              &*t6*t11*xm*xx*yy*1.6D1+a*t2*t4*t14*xm*xx*yy*2.0D1+a*t2*t6*t13*xm*x&
              &x*yy*1.2D1+a*t2*t6*t14*xm*xx*yy*4.0D0-a*t4*t7*t11*xm*ym*yy*1.6D1+a&
              &*t6*t7*t11*xm*ym*yy*8.0D1+a*t4*t7*t10*xx*ym*yy*8.0D0+a*t4*t7*t11*x&
              &x*ym*yy*2.0D1-a*t6*t7*t10*xx*ym*yy*4.0D1+a*t4*t7*t13*xx*ym*yy*8.0D&
              &0-a*t6*t7*t11*xx*ym*yy*1.0D2+a*t4*t7*t14*xx*ym*yy*8.0D0-a*t6*t7*t1&
              &3*xx*ym*yy*4.0D1-a*t6*t7*t14*xx*ym*yy*4.0D1-t2*t4*t6*t7*t10*xm*xx*&
              &1.0D1-t2*t4*t6*t7*t11*xm*xx*2.4D1-t2*t4*t6*t7*t13*xm*xx*6.0D0-t2*t&
              &4*t6*t7*t14*xm*xx*6.0D0-t2*t4*t6*t7*t10*ym*yy*4.0D0-t2*t4*t6*t7*t1&
              &1*ym*yy*2.0D1+t2*t4*t6*t7*xm*xx*ym*yy*1.2D1+a*t2*t7*t10*t37*xm*xx*&
              &ym*4.0D0+a*t2*t7*t11*t37*xm*xx*ym*1.6D1-a*t2*t7*t10*t39*xm*xx*ym*4&
              &.0D0-a*t2*t7*t11*t39*xm*xx*ym*1.6D1+a*t2*t7*t13*t37*xm*xx*ym*4.0D0&
              &+a*t2*t7*t14*t37*xm*xx*ym*1.2D1-a*t2*t7*t13*t39*xm*xx*ym*4.0D0-a*t&
              &2*t7*t14*t39*xm*xx*ym*1.2D1-a*t2*t7*t10*t37*xm*xx*yy*4.0D0-a*t2*t7&
              &*t11*t37*xm*xx*yy*1.6D1+a*t2*t7*t10*t39*xm*xx*yy*4.0D0+a*t2*t7*t11&
              &*t39*xm*xx*yy*1.6D1-a*t2*t7*t13*t37*xm*xx*yy*1.2D1-a*t2*t7*t14*t37&
              &*xm*xx*yy*4.0D0+a*t2*t7*t13*t39*xm*xx*yy*1.2D1+a*t2*t7*t14*t39*xm*&
              &xx*yy*4.0D0+a*t4*t6*t11*t35*xm*ym*yy*8.0D0-a*t4*t6*t11*t38*xm*ym*y&
              &y*2.4D1-a*t4*t6*t10*t35*xx*ym*yy*4.0D0-a*t4*t6*t11*t35*xx*ym*yy*8.&
              &0D0+a*t4*t6*t10*t38*xx*ym*yy*1.2D1+a*t4*t6*t11*t38*xx*ym*yy*3.2D1+&
              &a*t4*t6*t13*t38*xx*ym*yy*1.6D1+a*t4*t6*t14*t38*xx*ym*yy*1.6D1)-(sq&
              &rt(3.0D0)*t21**2*1.0D0/sqrt(-t31)*(t28+t2*t6*2.0D0-t35*t37+2.0D1))&
              &/(t19*tie*2.0D0)-t2*t4*t19*t21*t47*t51*t52*(ym*2.0D0-yy*2.0D0)*(5.&
              &0D0/6.0D0)+a*t2*t6*t17*t19*t21*t24*t47*(5.0D0/3.0D0)+a*t2*t4*t7*t1&
              &7*t24*t47*t50*(5.0D0/3.0D0)           
           END DO
              
              
              
      
      
      
      
!						! Axisimmetric case with div(b)~=0, n = 2+sin(2*pi*x)*sin(2*pi*y),  u = cos(2*pi*x)*cos(2*pi*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 10-sin(wx*x)*cos(wy*y)
!						   IF (.not.switch%axisym) THEN
!						      WRITE(6,*) "This is an axisymmetric test case!"
!						      stop
!						   END IF
!									a = 2*pi
!									b = 2*pi
!									xc = 0.
!									yc = 0.
!									D     = phys%diff_n
!									mu    = phys%diff_u
!									csii  = phys%diff_e
!									csie  = phys%diff_ee
!									kpari = phys%diff_pari
!									kpare = phys%diff_pare
!									Mref  = phys%Mref
!									epn   = phys%epn
!									tie   = phys%tie
!									pcourr = 1.

!									f(:,1) = cos(a * x) ** 2 * a * sin(b * y) * cos(b * y) * (x / 0.30D2 &
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

!									f(:,2) =  cos(a * x) ** 3 * a * sin(b * y) * cos(b * y) ** 2 * (x / 0.&
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

      CASE(5:6)
      ! Do nothing
				  CASE(50:)
				  ! Do nothing

				  CASE DEFAULT
				      WRITE(6,*) "Error! Test case not valid"
				      STOP
			
				END SELECT
  END SUBROUTINE body_force
#endif  
END MODULE analytical
