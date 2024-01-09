!*****************************************
! project: MHDG
! file: physics.f90
! date: 20/09/2017
! Define the physics of the model
!  ******** N-Gamma-Ti-Te system     ****
!*****************************************
MODULE physics

  USE globals
  USE printUtils
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
#ifdef NEUTRAL
    phys%Neq = 5
#endif
#ifdef NEUTRALGAMMA
    phys%Neq = 6
#endif

    ! number of physical variables
    phys%npv = 10
#ifdef NEUTRAL
    phys%npv = 11
#endif
#ifdef NEUTRALGAMMA
    phys%npv = 12
#endif

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho" ! density
    phys%phyVarNam(2) = "u"   ! parallel velocity
    phys%phyVarNam(3) = "Ei"  ! total energy of ions
    phys%phyVarNam(4) = "Ee"  ! total energy of electrons
    phys%phyVarNam(5) = "pi"  ! pressure of ions
    phys%phyVarNam(6) = "pe"  ! pressure of electrons
    phys%phyVarNam(7) = "Ti"  ! temperature of ions
    phys%phyVarNam(8) = "Te"  ! temperature of electrons
    phys%phyVarNam(9) = "Csi" ! sound speed
    phys%phyVarNam(10)= "M"   ! Mach
#ifdef NEUTRAL
    phys%phyVarNam(11)= "rhon"   ! density neutral
#endif
#ifdef NEUTRALGAMMA
    phys%phyVarNam(12)= "un"   ! parallel velocity neutral
#endif

    ! Set the name of the conservative variables
    phys%conVarNam(1) = "rho"   ! U1 = rho
    phys%conVarNam(2) = "Gamma" ! U2 = rho*u
    phys%conVarNam(3) = "nEi"   ! U3 = rho*Ei
    phys%conVarNam(4) = "nEe"   ! U4 = rho*Ee
#ifdef NEUTRAL
    phys%conVarNam(5) = "rhon"  ! U5 = rhon
#endif
#ifdef NEUTRALGAMMA
    phys%conVarNam(6) = "Gamman"! U6 = rhon*un
#endif

    simpar%model = 'N-Gamma-Ti-Te'
#ifdef NEUTRAL
    simpar%model = 'N-Gamma-Ti-Te-Neutral'
#endif
#ifdef NEUTRALGAMMA
    simpar%model = 'N-Gamma-Ti-Te-NeutralGamma'
#endif

    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq

    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
    simpar%physvar_refval(1) = simpar%refval_density
    simpar%physvar_refval(2) = 1.
    simpar%physvar_refval(3) = simpar%refval_specenergy
    simpar%physvar_refval(4) = simpar%refval_specenergy
    simpar%physvar_refval(5) = simpar%refval_specpress
    simpar%physvar_refval(6) = simpar%refval_specpress
    simpar%physvar_refval(7) = simpar%refval_temperature
    simpar%physvar_refval(8) = simpar%refval_temperature
    simpar%physvar_refval(9) = simpar%refval_speed
    simpar%physvar_refval(10) = 1.
#ifdef NEUTRAL
    simpar%physvar_refval(11) = simpar%refval_neutral
#endif
#ifdef NEUTRALGAMMA
    simpar%physvar_refval(12) = simpar%refval_speed
#endif
    simpar%consvar_refval(1) = simpar%refval_density
    simpar%consvar_refval(2) = simpar%refval_momentum
    simpar%consvar_refval(3) = simpar%refval_specenergydens
    simpar%consvar_refval(4) = simpar%refval_specenergydens
#ifdef NEUTRAL
    simpar%consvar_refval(5) = simpar%refval_neutral
#endif
#ifdef NEUTRALGAMMA
    simpar%consvar_refval(6) = simpar%refval_momentum
#endif
  END SUBROUTINE

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up, ua)
    real*8, dimension(:, :), intent(in)  :: up
    real*8, dimension(:, :), intent(out) :: ua

    ua(:, 1) = abs(up(:, 1))
    ua(:, 2) = up(:, 1)*up(:, 2)
    ua(:, 3) = up(:, 1)*up(:, 3)
    ua(:, 4) = up(:, 1)*up(:, 4)
#ifdef NEUTRAL
    ua(:,5) = abs(up(:,11))
#endif
#ifdef NEUTRALGAMMA
    ua(:,6) = up(:,11)*up(:,12)
#endif

  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua, up)
    real*8, dimension(:, :), intent(in)  :: ua
    real*8, dimension(:, :), intent(out) :: up
    real*8,  dimension(size(ua,1))  :: U1,U5
    real, parameter :: tol = 1e-12
    integer :: i

    U1 = abs(ua(:,1))
#ifdef NEUTRAL
    U5 = abs(ua(:,5))
#endif
!    do i = 1, size(ua,1)
!      if (U1(i)<tol) U1(i)=tol
!#ifdef NEUTRAL
!      if (U5(i)<tol) U5(i)=tol
!#endif
!    end do

    up(:, 1) = abs(U1)                                                ! density
    up(:, 2) = ua(:, 2)/U1                                            ! u parallel
    up(:, 3) = ua(:, 3)/U1                                            ! total energy of ions
    up(:, 4) = ua(:, 4)/U1                                            ! total energy of electrons
    up(:, 5) = abs(2./(3.*phys%Mref)*(ua(:, 3) - 0.5*ua(:, 2)**2/U1)) ! pressure of ions
    up(:, 6) = abs(2./(3.*phys%Mref)*ua(:, 4))                        ! pressure of electrons
    up(:, 7) = up(:, 5)/U1                                            ! temperature of ions
    up(:, 8) = up(:, 6)/U1                                            ! temperature of electrons
    up(:, 9) = sqrt((up(:, 7) + up(:, 8))*phys%Mref)                  ! sound speed
    up(:, 10) = up(:, 2)/up(:, 9)                                     ! Mach
#ifdef NEUTRAL
    up(:,11) = abs(U5)                                                ! density neutral
#endif
#ifdef NEUTRALGAMMA
    up(:,12) = ua(:,6)/U5                                             ! u parallel neutral
#endif

    ! Set threshold for low density and temperature
    !DO i = 1, size(ua,1) 
    !   IF (up(i,1) < 1.e-6) up(i,1) = 1.e-6
    !   IF (up(i,7) < 1.e-3) up(i,7) = 1.e-3
    !   IF (up(i,8) < 1.e-3) up(i,8) = 1.e-3
    !END DO

  END SUBROUTINE cons2phys

    ! ******************************
    ! Split diffusion terms
    ! ******************************
    SUBROUTINE computeVu(U, V)
      real*8, intent(IN)  :: U(:)
      real*8, intent(OUT) :: V(:)

      V = 0.d0

      V(1) = -U(2)/U(1)**2
      V(2) = 1./U(1)
    END SUBROUTINE computeVu

    SUBROUTINE compute_dVu_dU(U,dV_dU)
      real*8, intent(IN)  :: U(:)
      real*8, intent(OUT) :: dV_dU(:, :)

      dV_dU = 0.

      dV_dU(1, 1) = 2*U(2)/U(1)**3
      dV_dU(1, 2) = -1/U(1)**2

      dV_dU(2, 1) = -1/U(1)**2
    END SUBROUTINE compute_dVu_dU

    SUBROUTINE compute_W2(U,W2)
      real*8, intent(IN) :: U(:)
      real*8             :: diff_n,diff_u
      real*8             :: D(4),W2(:)

      W2 = 0.
      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_u = D(2)

      W2(1) = (diff_n - diff_u)*U(2)/U(1)
    END SUBROUTINE compute_W2

    SUBROUTINE compute_W3(U,W3)
      real*8, intent(IN) :: U(:)
      real*8             :: D(4),W3(:)
      real*8 		 :: diff_n,diff_u,diff_e
      real*8 		 :: rhovar,sigmavar

      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_u = D(2)
      diff_e = D(3)
      rhovar = (diff_e - diff_u)
      sigmavar = (diff_n - diff_e)

      W3 = 0.
      W3(1) = sigmavar*U(3)/U(1) + rhovar*(U(2)/U(1))**2
      W3(2) = -rhovar*U(2)/U(1)
    END SUBROUTINE compute_W3

    SUBROUTINE compute_W4(U,W4)
      real*8, intent(IN) :: U(:)
      real*8             :: D(4),W4(:)
      real*8             :: diff_n,diff_ee

      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_ee = D(4)

      W4 = 0.
      W4(1) = (diff_n - diff_ee)*U(4)/U(1)
    END SUBROUTINE compute_W4

    SUBROUTINE compute_dW2_dU(U,dW2_dU)
      real*8, intent(IN) :: U(:)
      real*8             :: D(4),dW2_dU(:,:)
      real*8             :: diff_n,diff_u

      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_u = D(2)

      dW2_dU = 0.
      dW2_dU(1,1) = -U(2)/(U(1)**2)
      dW2_dU(1,2) = 1./U(1)

      dW2_dU = (diff_n - diff_u)*dW2_dU
    END SUBROUTINE compute_dW2_dU

    SUBROUTINE compute_dW3_dU(U,res)
      real*8, intent(IN) :: U(:)
      real*8             :: D(4),res(:,:)
      real*8 		 :: diff_n,diff_u,diff_e
      real*8 		 :: rhovar,sigmavar

      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_u = D(2)
      diff_e = D(3)

      rhovar = (diff_e - diff_u)
      sigmavar = (diff_n - diff_e)

      res = 0.
      res(1,1) = -sigmavar*U(3)/(U(1)**2)-2*rhovar*(U(2)**2)/(U(1)**3)
      res(1,2) = 2*rhovar*U(2)/(U(1)**2)
      res(1,3) = sigmavar*1./U(1)

      res(2,1) = rhovar*U(2)/(U(1)**2)
      res(2,2) = -rhovar*1./U(1)
    END SUBROUTINE compute_dW3_dU

    SUBROUTINE compute_dW4_dU(U,res)
      real*8, intent(IN) :: U(:)
      real*8             :: D(4),res(:,:)
      real*8             :: diff_n,diff_ee

      CALL setLocalDiffSplitTerms(U,D)
      diff_n = D(1)
      diff_ee = D(4)

      res = 0.
      res(1,1) = -U(4)*(diff_n - diff_ee)/(U(1)**2)
      res(1,4) = 1.*(diff_n - diff_ee)/U(1)
    END SUBROUTINE compute_dW4_dU


  !*****************************************
  ! Jacobian matrices
  !****************************************
  SUBROUTINE jacobianMatrices(U, A)
    real*8, intent(in)  :: U(:)
    real*8, intent(out) :: A(:, :)
    ![ 0,                                           1,                              0,                0; ...
    ! -2/3*U(2)**2/U(1)**2                          4/3*U(2)/U(1)                   2/3,              2/3; ...
    ! -5/3*U(2)*U(3)/U(1)**2+2/3*U(2)**3/U(1)**3    5/3*U(3)/U(1)-U(2)**2/U(1)**2   5/3*U(2)/U(1),    0 ;   ...
    ! -5/3*U(4)*U(2)/U(1)**2,                       5/3*U(4)/U(1),                  0,                5/3*U(2)/U(1)]

    A = 0.d0
    if (switch%decoup) then
      A(1, 2) = 1.

      A(2, 1) = (-U(2)**2/U(1)**2 + phys%Mref)
      A(2, 2) = 2.*U(2)/U(1)

      A(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      A(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      A(3, 3) = 5./3.*U(2)/U(1)

      A(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      A(4, 2) = 5./3.*U(4)/U(1)
      A(4, 4) = 5./3.*U(2)/U(1)
    else

      A(1, 2) = 1.

      A(2, 1) = -2./3.*U(2)**2/U(1)**2
      A(2, 2) = 4./3.*U(2)/U(1)
      A(2, 3) = 2./3.
      A(2, 4) = 2./3.

      A(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      A(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      A(3, 3) = 5./3.*U(2)/U(1)

      A(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      A(4, 2) = 5./3.*U(4)/U(1)
      A(4, 4) = 5./3.*U(2)/U(1)

#ifdef NEUTRAL
#ifndef NEUTRALGAMMA
      !A(5, 1) = -U(5)*U(2)/U(1)**2
      !A(5, 2) = U(5)/U(1)
      !A(5, 5) = U(2)/U(1)
      !A(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*A(5,:)
#else
      A(5, 6) = 1.
      
      A(6, 1) = 2./3.*U(5)*(- U(3)/U(1)**2 + U(2)**2/U(1)**3)
      A(6, 2) = - 2./3.*U(5)*U(2)/U(1)**2
      A(6, 3) = 2./3.*U(5)/U(1)
      A(6, 5) = - U(6)**2/U(5)**2 + 2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2)
      A(6, 6) = 2.*U(6)/U(5)
#endif
#endif 
    end if
  END SUBROUTINE jacobianMatrices

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesNP(U, Anp)
    real*8, intent(in)  :: U(:)
    real*8, intent(out) :: Anp(:)
    real*8              :: cs_n

    Anp = 0.d0
    cs_n = sqrt(abs(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2)))

    Anp(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anp(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anp(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anp(5) = cs_n
  END SUBROUTINE jacobianMatricesNP
#endif

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U, bn, An)
    real*8, intent(in)  :: U(:), bn
    real*8, intent(out) :: An(:, :)
    An = 0.d0
    if (switch%decoup) then
      An(1, 2) = 1.

      An(2, 1) = (-U(2)**2/U(1)**2 + phys%Mref)
      An(2, 2) = 2.*U(2)/U(1)

      An(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      An(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      An(3, 3) = 5./3.*U(2)/U(1)

      An(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      An(4, 2) = 5./3.*U(4)/U(1)
      An(4, 4) = 5./3.*U(2)/U(1)
    else
      An(1, 2) = 1.

      An(2, 1) = -2./3.*U(2)**2/U(1)**2
      An(2, 2) = 4./3.*U(2)/U(1)
      An(2, 3) = 2./3.
      An(2, 4) = 2./3.

      An(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      An(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      An(3, 3) = 5./3.*U(2)/U(1)

      An(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      An(4, 2) = 5./3.*U(4)/U(1)
      An(4, 4) = 5./3.*U(2)/U(1)
      
#ifdef NEUTRAL
      !An(5, 1) = -U(5)*U(2)/U(1)**2
      !An(5, 2) = U(5)/U(1)
      !An(5, 5) = U(2)/U(1)
      !An(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*An(5,:)
#endif      
    endif
    An = bn*An
  END SUBROUTINE jacobianMatricesFace

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesFaceNP(U, bn, Anpn)
    real*8, intent(in)  :: U(:), bn
    real*8, intent(out) :: Anpn(:)
    real*8              :: cs_n

    Anpn = 0.d0
    cs_n = sqrt(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2))

    Anpn(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anpn(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anpn(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anpn(5) = cs_n

    Anpn = bn*Anpn
  END SUBROUTINE jacobianMatricesFaceNP
#endif

  !*****************************************
  ! Jacobian matrix for the Bohm BC
  !****************************************
  SUBROUTINE jacobianMatricesBohm(U, A)
    real*8, intent(in)  :: U(:)
    real*8, intent(out) :: A(:, :)
    real*8              :: auxi, auxe, RE

    A = 0.
    RE = 0.3
    
    auxi = (5.-2.*phys%Gmbohm)/3.
    auxe = (5.-2.*phys%Gmbohme)/3.
    
    A(1, 1) = -auxi*(U(2)**3/U(1)**3 - U(2)*U(3)/U(1)**2)
    A(1, 2) = auxi*(U(3)/U(1) - 3./2.*U(2)**2/U(1)**2)
    A(1, 3) = auxi*U(2)/U(1)

    A(2, 1) = -auxe*U(2)*U(4)/U(1)**2
    A(2, 2) = auxe*U(4)/U(1)
    A(2, 4) = auxe*U(2)/U(1)
    
#ifdef NEUTRAL
    ! Energy recycled by neutrals
    !A(1, 1) = A(1,1) - RE*(U(3) + U(4))*U(2)/U(1)**2 !+ RE*(U(2)/U(1))**3
    !A(1, 2) = A(1,2) + RE*(U(3) + U(4))/U(1) !- RE*3./2.*(U(2)/U(1))**2
    !A(1, 3) = A(1,3) + RE*U(2)/U(1)
    !A(1, 4) = A(1,4) + RE*U(2)/U(1)
    ! Energy recycled by neutrals: convective neutral flux u_n = u, Tn = Ti
    !A(1, 1) = A(1,1) - RE*U(5)*U(3)*U(2)*2./U(1)**3 + RE*3./2.*U(5)*U(2)**3/U(1)**4
    !A(1, 2) = A(1,2) + RE*U(5)*U(3)/U(1)**2 - 3./2.*U(2)**2/U(1)**3
    !A(1, 3) = A(1,3) + RE*U(5)*U(2)/U(1)**2
    !A(1, 5) = A(1,5) + RE*U(3)*U(2)/U(1)**2 - 1./2.*(U(2)/U(1))**3

    ! Pinch flux
    !A(5, 1) = -U(5)*U(2)/U(1)**2
    !A(5, 2) = U(5)/U(1)
    !A(5, 5) = U(2)/U(1) 
    !A(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*A(5,:)
#endif
  END SUBROUTINE jacobianMatricesBohm

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesBohmNP(U, Anp)
    real*8, intent(in)  :: U(:)
    real*8, intent(out) :: Anp(:)
    real*8              :: cs_n

    Anp = 0.d0
    cs_n = sqrt(abs(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2)))
    
    Anp(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anp(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anp(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anp(5) = cs_n
  END SUBROUTINE jacobianMatricesBohmNP
#endif
  
#ifdef NEUTRAL
  !*****************************************
  ! Jacobian matrices for Neutrals
  !****************************************
  SUBROUTINE jacobianMatricesN(U, Up, Q, b, sigmaviz, sigmavcx, Ax, Ay)
  real*8, intent(in)  :: U(:), Up(:), Q(:,:), b(:), sigmaviz, sigmavcx
  real*8, intent(out) :: Ax(:,:),Ay(:,:)
  real*8              :: GradTi(simpar%Ndim), GradTimod, Vnn, Csnn
  
  GradTi(1) = 2./(3*phys%Mref)*( (U(2)**2/U(1)**3 - U(3)/U(1)**2)*Q(1,1) - (U(2)/U(1)**2)*Q(1,2) + 1./U(1)*Q(1,3) )
  GradTi(2) = 2./(3*phys%Mref)*( (U(2)**2/U(1)**3 - U(3)/U(1)**2)*Q(2,1) - (U(2)/U(1)**2)*Q(2,2) + 1./U(1)*Q(2,3) )
  GradTi = simpar%refval_temperature/simpar%refval_length*GradTi
  
  GradTimod = sqrt(GradTi(1)**2 + GradTi(2)**2)
  
  Vnn = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTimod/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
  Csnn = sqrt(simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass*Up(7))
   
  Ax = 0.d0
  !Neutral convective velocity
!  if (Csnn .ge. Vnn) then
!     Ax(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(1)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
!  else
!     Ax(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(1)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))*Csnn/Vnn
!  endif
!  Ax = Ax/simpar%refval_speed
  !Neutral velocity parallel to magnetic field
  Ax(5,5) = Ax(5,5) - Up(2)*b(1)
  
  Ay = 0.d0
!  if (Csnn .ge. Vnn) then
!     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
!  else
!     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))*Csnn/Vnn
!  endif
!  Ay = Ay/simpar%refval_speed
  !Neutral velocity parallel to magnetic field
  Ay(5,5) = Ay(5,5) - Up(2)*b(2)
  
  Ax = 0.
  Ay = 0.
  
  END SUBROUTINE jacobianMatricesN
#endif

  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy, u, q, d_iso, d_ani)
    real*8, intent(in)  		:: xy(:, :)
    real*8, intent(in)  		:: u(:,:), q(:,:)
    real*8, intent(out)		 :: d_iso(:, :, :), d_ani(:, :, :)
    real*8		              :: iperdiff(size(xy, 1))
#ifdef NEUTRAL
    integer             		:: i
    real*8				            :: Ery = 13.6, cs_n, DnnTh
    real*8, dimension(size(u,1))	:: U1, U2, U3, U4, U5, E0iz, E0cx, sigmaviz, sigmavcx, Dnn
#endif

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
    d_iso(3, 3, :) = phys%diff_e
    d_iso(4, 4, :) = phys%diff_ee
    if (switch%ME .eq. .TRUE.  .AND. switch%testcase .gt. 84) then
       if (switch%testcase .eq. 85) then !Iter core-edge with evolving equilibria plus diffusion decrease
          d_iso(1, 1, :) = phys%diff_n - (phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(2, 2, :) = phys%diff_u - (phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(3, 3, :) = phys%diff_e - (phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(4, 4, :) = phys%diff_ee - (phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
      else if (switch%testcase .eq. 86) then
          d_iso(1,1,:) = max(switch%diffmin,phys%diff_n - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(2,2,:) = max(switch%diffmin,phys%diff_u - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(3,3,:) = max(switch%diffmin,phys%diff_e - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(4,4,:) = max(switch%diffmin,phys%diff_ee - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
      else if (switch%testcase .eq. 87) then
          d_iso(1,1,:) = max(switch%diffmin,phys%diff_n - 2.38*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(2,2,:) = max(switch%diffmin,phys%diff_u - 2.38*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(3,3,:) = max(switch%diffmin,phys%diff_e - 2.0134*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_iso(4,4,:) = max(switch%diffmin,phys%diff_ee - 2.0134*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
      endif
      phys%ME_diff_n = d_iso(1,1,1)
      phys%ME_diff_u = d_iso(2,2,1)
      phys%ME_diff_e = d_iso(3,3,1)
      phys%ME_diff_ee = d_iso(4,4,1)
    endif  
#ifndef NEUTRALP
#ifdef NEUTRAL
    !d_iso(5, 5, :) = phys%diff_nn
    U1 = u(:,1)
    U2 = u(:,2)
    U3 = u(:,3)
    U4 = u(:,4)
    U5 = u(:,5)
    DO i=1,size(u,1)
       CALL compute_sigmaviz(u(i,:),sigmaviz(i))
       CALL compute_sigmavcx(u(i,:),sigmavcx(i))
    END DO
    Dnn = simpar%refval_charge*max(simpar%refval_temperature*2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2),0.1)/(simpar%refval_mass*simpar%refval_density*U1*(sigmaviz + sigmavcx))   
    Dnn = Dnn*simpar%refval_time/simpar%refval_length**2
    Dnn(:) = Dnn(:)*(1 + (Dnn(:)/phys%diff_nn)**20)**(-1./20)
    d_iso(5,5,:) = Dnn
    !Set a threshold on Dnn
    !DO i=1,size(Dnn,1)
       !if (Dnn(i) .gt.  phys%diff_nn) d_iso(5,5,i) = phys%diff_nn
       !d_iso(5,5,i)= d_iso(5,5,i)*(1 + (d_iso(5,5,i)/phys%diff_nn)**1)**(-1./1)
       !if (Dnn(i) .lt.  19) d_iso(5,5,i) = 19
       !if (Dnn(i) .lt. 50*simpar%refval_time/simpar%refval_length**2) d_iso(5,5,i) = 50*simpar%refval_time/simpar%refval_length**2
       !cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3(i)/U1(i) - 1./2.*(U2(i)/U1(i))**2)))
       !DnnTh = abs(cs_n*U5(i)/(sqrt(q(i,9)**2 + q(i,10)**2)))
       !!Max Value
       !if (Dnn(i) .gt.  DnnTh) then 
	      !   if (DnnTh .gt. 20.*phys%diff_n) then 
	      !      d_iso(5,5,i) = DnnTh
       !      if (DnnTh .gt. phys%diff_nn) d_iso(5,5,i) = phys%diff_nn
	      !   end if 
	      !!Min Value
	      !!if (Dnn(i) .lt.  phys%diff_n) d_iso(5,5,i) = phys%diff_n
       !end if
    !END DO
    !Dnn = sum(Dnn)/size(u,1)
#endif
#else
    d_iso(5,5,:) = 0.
#endif
    d_ani(1, 1, :) = phys%diff_n
    d_ani(2, 2, :) = phys%diff_u
    d_ani(3, 3, :) = phys%diff_e
    d_ani(4, 4, :) = phys%diff_ee
    if (switch%ME .eq. .TRUE. .AND. switch%testcase .gt. 84) then !Iter core-edge with evolving equilibria plus diffusion decrease
       if (switch%testcase .eq. 85) then
         d_ani(1, 1, :) = phys%diff_n - (phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(2, 2, :) = phys%diff_u - (phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(3, 3, :) = phys%diff_e - (phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(4, 4, :) = phys%diff_ee - (phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
       else if (switch%testcase .eq. 86) then
          d_ani(1,1,:) = max(switch%diffmin,phys%diff_n - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(2,2,:) = max(switch%diffmin,phys%diff_u - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(3,3,:) = max(switch%diffmin,phys%diff_e - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(4,4,:) = max(switch%diffmin,phys%diff_ee - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
       else if (switch%testcase .eq. 87) then
          d_ani(1,1,:) = max(switch%diffmin,phys%diff_n - 2.38*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(2,2,:) = max(switch%diffmin,phys%diff_u - 2.38*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(3,3,:) = max(switch%diffmin,phys%diff_e - 2.0134*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          d_ani(4,4,:) = max(switch%diffmin,phys%diff_ee - 2.0134*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
       endif
    endif  
#ifdef NEUTRAL
    d_ani(5,5,:) = 0.
#endif

    !*****************************
    ! Non diagonal terms
    !*****************************
    ! No non-diagonal terms defined for this model
    call computeIperDiffusion(xy,u, iperdiff)
    d_iso(1, 1, :) = d_iso(1, 1, :)+iperdiff
    d_iso(2, 2, :) = d_iso(2, 2, :)+iperdiff
    d_iso(3, 3, :) = d_iso(3, 3, :)+iperdiff
    d_iso(4, 4, :) = d_iso(4, 4, :)+iperdiff
    d_ani(1, 1, :) = d_ani(1, 1, :) + iperdiff
    d_ani(2, 2, :) = d_ani(2, 2, :) + iperdiff
    d_ani(3, 3, :) = d_ani(3, 3, :) + iperdiff
    d_ani(4, 4, :) = d_ani(4, 4, :) + iperdiff
  END SUBROUTINE setLocalDiff

  !*****************************************
  ! Compute local perpendicular diffusion 
  ! for split diffusion terms for plasma
  !*****************************************
  SUBROUTINE setLocalDiffSplitTerms(u, D)
    real*8, intent(IN)    :: u(:)
    real*8, intent(OUT)   :: D(:)

    ! Reference case: diffusion from param file
    D(1) = phys%diff_n
    D(2) = phys%diff_u
    D(3) = phys%diff_e
    D(4) = phys%diff_ee

    ! Diffusion evlution for ME cases
    ! Iter core-edge with evolving equilibria plus diffusion decrease
    if (switch%ME .eq. .TRUE.  .AND. switch%testcase .gt. 84) then
       if (switch%testcase .eq. 85) then 
          D(1) = phys%diff_n - (phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          D(2) = D(1)
          D(3) = D(1)
          D(4) = D(1)
      else if (switch%testcase .eq. 86) then
          D(1) = max(switch%diffmin,phys%diff_n - 2.18*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          D(2) = D(1)
          D(3) = D(1)
          D(4) = D(1)
      else if (switch%testcase .eq. 87) then
          D(1) = max(switch%diffmin,phys%diff_n - 2.38*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          D(2) = D(1)
          D(3) = max(switch%diffmin,phys%diff_e - 2.0134*simpar%refval_time/simpar%refval_length**2*(tanh((phys%I_p - 0.35)/5.)))
          D(4) = D(3)
      endif
    endif

  END SUBROUTINE

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X, u, ipdiff)
    real*8, intent(IN)  :: X(:, :), u(:,:)
    real*8, intent(OUT) :: ipdiff(:)
    real*8             :: d,dref,maxamp
    real*8,allocatable :: xcorn(:), ycorn(:)
    real*8             :: rad(size(X, 1))
    real*8             :: h,rhog,up(size(u,1),phys%npv),Ti(size(u,1)),Te(size(u,1))
    real*8, parameter   :: tol = 1.e-6
    integer            :: i,g, opt

    ipdiff = 0.

    if (switch%difcor .gt. 0) then

						 SELECT CASE (switch%difcor)
										CASE (1)
												! Circular case with infinitely small limiter
												allocate(xcorn(1),ycorn(1))
												xcorn = geom%R0
												ycorn = -0.75
										CASE (2)
												! Circular case with infinitely small limiter
												allocate(xcorn(1),ycorn(1))
												xcorn = geom%R0
												ycorn = -0.287
										CASE (3)
												! West
												allocate(xcorn(1),ycorn(1))
												xcorn = 2.7977
												ycorn = -0.5128
										case(4)
										 ! ITER
										 allocate(xcorn(4),ycorn(4))
												xcorn(1) =  4.023150000000000
												ycorn(1) = -2.544840000000000
												xcorn(2) =  6.23648
												ycorn(2) = -3.23689
												xcorn(3) =  4.02837
												ycorn(3) =  3.588
												xcorn(4) =  5.74627
												ycorn(4) =  4.51401
										CASE DEFAULT
												WRITE (6, *) "Case not valid"
												STOP
						 END SELECT

						 !!**********************************************************
						 !! Gaussian around the corner
						 !!**********************************************************
						 h = 1e-1
						 do i=1,size(xcorn)
										rad = sqrt((X(:, 1)*phys%lscale - xcorn(i))**2 + (X(:, 2)*phys%lscale - ycorn(i))**2)
										ipdiff = ipdiff + numer%dc_coe*exp(-(2*rad/h)**2)

       end do
       do i=1,size(ipdiff)
          if (ipdiff(i)<1e-5) ipdiff(i)=0.
       end do
       deallocate(xcorn,ycorn)
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

    !Temperature IperDiff
    !CALL cons2phys(u,up)
    !Ti = up(:,7)
    !Te = up(:,8)
                                                                                                                   
    !IF (minval(Te) .le. 6.e-4) THEN                                                                                                                      !    !do g = 1, size(Ei,1)
    !   ipdiff = 2*tanh(6.e-4/minval(abs(Te)) - 1.)*simpar%refval_time/simpar%refval_length**2
    !   !WRITE(6,*) 'size ipdiff = ', size(ipdiff,1)                                                                                                    
       !WRITE(6,*) 'ipdiff = ', ipdiff                                                                                                                 
    !END IF

  END SUBROUTINE computeIperDiffusion
   
  !*****************************************
  ! Pinch term 
  !***************************************
!  SUBROUTINE setPinch(u)
!    real*8,intent(IN)                    :: u(:,:)
!    real*8                               :: nu,nu0 = 0.
!    real*8,dimension(refElPol%Nnodes2D)  :: psiel,Yel,nel,Tel
!    integer*4                            :: iel,N2D,i,inde(refElPol%Nnodes2D),Npel,inode,ierr

!    !u = transpose(reshape(sol%u,[phys%Neq,size(sol%u)/phys%Neq]))
!    N2D = Mesh%Nelems
!    Npel = refElPol%Nnodes2D

    

!    ! Loop on the number of elements
!    DO iel = 1,N2D
!       ! Z coordinates of the nodes of the element
!       Yel = Mesh%X(Mesh%T(iel,:),2)

!       ! Normalazied poloidal flux of the nodes of the element
!       psiel = phys%magnetic_psi(Mesh%T(iel,:))
       
!       ! Indices and elemental solution
!       inde = (iel - 1)*Npel + (/(i,i=1,Npel)/)
!       nel = u(inde,1)
!       Tel = (2./(3.*phys%Mref))*(u(inde,4)/u(inde,1))

!       ! Loop on the number of nodes
!#ifdef PARALL
!       IF (Mesh%ghostElems(iel) .eq. 0) THEN
!#endif
!          DO inode = 1,Npel
!             IF (switch%testcase .ge. 80 .and. switch%testcase .lt. 90) THEN !Machine dependent: ITER testcase
!                IF (psiel(inode) .le. 0.95 .and. Yel(inode) .gt. -3.2/simpar%refval_length) THEN
!                   nu = nel(inode)/(Tel(inode)**(3./2.))
!                   IF (nu .gt. nu0) nu0 = nu
!                END IF
!             END IF
!          END DO
!#ifdef PARALL
!       END IF
!#endif
!    END DO

!    phys%nu0 = nu0

!#ifdef PARALL
!    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%nu0, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
!#endif

!  END SUBROUTINE
  
  SUBROUTINE computeNuth(nu_th)
    real*8, intent(OUT)   :: nu_th
    
    nu_th = 0.

    ! Decrease collisionality threshold for pinch term 
    ! ITER: nu_th = 0.04 at flat-top phase as E. Militello-Asp NF 2022
    if (switch%testcase .gt. 80) then
       nu_th = max(0.04,1. - 0.96*tanh((phys%I_p - 0.35)/5.))
    end if

  END SUBROUTINE

  SUBROUTINE computeCollisionality(U,q,nu)
    ! Compute normalized electron collisionality
    ! Reference O.Sauter PoP 1999
    real*8, intent(IN)   :: U(:),q
    real*8, intent(OUT)  :: nu
    real*8               :: R0,a,epsilon,n_e,T_e,lambda_e,Zeff
    
    nu = 0.

    ! ITER testcase
    IF (switch%testcase .ge. 80 .and. switch%testcase .lt. 90) THEN
       R0 = 6.2
       a = 2.
       Zeff = 1.6
       epsilon = a/R0
       n_e = U(1)*simpar%refval_density
       T_e = 2./(3.*phys%Mref)*U(4)/U(1)*simpar%refval_temperature
       lambda_e = 31.3 - log(sqrt(abs(n_e))/T_e)
       nu = 6.921e-18*(R0*q*n_e*Zeff*lambda_e)/((epsilon**(1.5))*(T_e**(2.)))
    END IF

  END SUBROUTINE

  SUBROUTINE computePinch(U,b,psi,q,APinch)
    real*8, intent(IN)     :: U(:),b(:),psi,q
    real*8, intent(OUT)    :: APinch(:,:)
    real*8                 :: v_p,bnorm(2)
    real*8                 :: nu,nu_star,f_nu,nu_th,D(4)
  
    ! Inintialize matrix and parameters
    APinch = 0.
    v_p = 0.
  
    ! Normalized poloidal magnetic field vector
    bnorm = b(:)/norm2(b)
  
    ! IMPLEMENTATION: v_p = v_p(psi) ad hoc profile  
    IF (psi .le. 0.95) THEN
       v_p = phys%v_p*(psi**2 + psi**2*tanh((0.95 - psi)/0.02))
       !IF (v_p .lt. 1.e-4/simpar%refval_speed) v_p = 0.   
    END IF

    ! IMPLEMENTATION: v_p = v_p(nu) self consistent profile as a function of 
    ! normalized electorn collisionality 
    ! Tuned as in JINTRAC (see E. Militello-Asp NF 2022)
    !f_nu = 0.
    !IF (psi .le. 0.95) THEN ! Apply just in the closed field line region (MUST balance BC for PFR )
    !   CALL computeNuth(nu_th)
    !   CALL computeCollisionality(U,q,nu) 
    !   CALL setLocalDiffSplitTerms(U, D)
    !   D = D*simpar%refval_length**2/simpar%refval_time
    !   nu_Star = nu/nu_th
    !   !nu_star = min(nu/nu_th,10.)
    !   f_nu = exp(1 - nu_star)
    !   v_p = f_nu*0.5*D(1)*sqrt(abs(psi))
    !   v_p = phys%v_p*v_p
    !END IF
    
    phys%v_pmax = min(phys%v_pmax,v_p)

    ! Assembly Pinch Matrix
    APinch(1,1) = v_p*bnorm(2) 
    APinch(1,2) = v_p*(-bnorm(1))

  END SUBROUTINE
    
  !*****************************************
  ! Curvature term matrix
  !****************************************
  SUBROUTINE GimpMatrix(U, divb, G)
    real*8, intent(in)  :: U(:), divb
    real*8, intent(out) :: G(:, :)

    !G = divb*[              0,                            0,                            0                              0; ...
    !                  1/3*U(2)**2/U(1)**2            -2/3*U(2)/U(1)                      2/3,                           2/3;...
    !                        0                             0                             0                              0; ...
    !                        0                             0                             0                              0];

    G = 0.d0
    if (switch%decoup) then
      G(2, 1) = phys%Mref
    else
      G(2, 1) = 1./3.*U(2)**2/U(1)**2
      G(2, 2) = -2./3.*U(2)/U(1)
      G(2, 3) = 2./3.
      G(2, 4) = 2./3.
    end if
    G = divb*G
  END SUBROUTINE GimpMatrix

  !*****************************************
  ! Ion pressure gradient
  !****************************************
  SUBROUTINE compute_dPi_dU(U,dPi_dU)
    real*8, intent(IN)   :: U(:)
    real*8, intent(OUT)  :: dPi_dU(:)

    dPi_dU = 0.
    
    dPi_dU(1) = 1./2.*(U(2)/U(1))**2
    dPi_dU(2) = -U(2)/U(1)
    dPi_dU(3) = 1.

  END SUBROUTINE
  
  !*****************************************
  ! Parallel diffusion terms
  !****************************************
  SUBROUTINE computeVi(U, V)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: V(:)
    V = 0.d0
    V(1) = U(2)**2/U(1)**3 - U(3)/U(1)**2
    V(2) = -U(2)/U(1)**2
    V(3) = 1./U(1)
  END SUBROUTINE computeVi

  SUBROUTINE computeVe(U, V)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: V(:)
    V = 0.d0
    V(1) = -U(4)/U(1)**2
    V(4) = 1./U(1)
  END SUBROUTINE computeVe

  !                                SUBROUTINE computeVe(U,V)
  !                                real*8, intent(IN)  :: U(:)
  !                                real*8, intent(OUT) :: V(:)
  !                                V = 0.d0
  !                                V(1) = U(2)**2/U(1)**3 - U(4)/U(1)**2
  !                                V(2) = -U(2)/U(1)**2
  !                                V(4) = 1./U(1)
  !                                END SUBROUTINE computeVe

  SUBROUTINE compute_dV_dUi(U, dV_dU)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: dV_dU(:, :)
    dV_dU = 0.
    dV_dU(1, 1) = 2*U(3)/U(1)**3 - 3*U(2)**2/U(1)**4
    dV_dU(1, 2) = 2*U(2)/U(1)**3
    dV_dU(1, 3) = -1/U(1)**2
    dV_dU(2, 1) = 2*U(2)/U(1)**3
    dV_dU(2, 2) = -1/U(1)**2
    dV_dU(3, 1) = -1/U(1)**2
  END SUBROUTINE compute_dV_dUi

  SUBROUTINE compute_dV_dUe(U, dV_dU)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: dV_dU(:, :)
    dV_dU = 0.
    dV_dU(1, 1) = 2*U(4)/U(1)**3
    dV_dU(1, 4) = -1/U(1)**2
    dV_dU(4, 1) = -1/U(1)**2
  END SUBROUTINE compute_dV_dUe

  FUNCTION computeAlphai(U) RESULT(res)
    real*8 :: U(:)
    real*8 :: res, aux
    real, parameter :: tol = 1e-5
    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
    if ((2./(3.*phys%Mref)*aux > 1.) .and. (switch%testcase .ne. 2)) then
      res = (3.*phys%Mref/2)**(phys%epn)
    else
      if (aux<tol) aux = tol
      res = aux**phys%epn
    endif
  END FUNCTION computeAlphai

  FUNCTION computeAlphae(U) RESULT(res)
    real*8 :: U(:)
    real*8 :: res, aux
    real, parameter :: tol = 1e-5
    aux = U(4)/U(1)
    if ((2./(3.*phys%Mref)*aux > 1.) .and. (switch%testcase .ne. 2)) then
      res = (3.*phys%Mref/2)**(phys%epn)
    else
      if (aux<tol) aux = tol
      res = aux**phys%epn
    endif
  END FUNCTION computeAlphae

  SUBROUTINE compute_dAlpha_dUi(U, res)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: res(:)
    real*8             :: aux
    real, parameter :: tol = 1e-5
    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
    if ((2./(3.*phys%Mref)*aux > 1.) .and. (switch%testcase .ne. 2)) then !! don't apply flux limiter if it is a convergence test
      res = 0.
    else
      if (aux<0) aux = tol
      res = 0.d0
      res(1) = -U(3)/U(1)**2+U(2)**2/U(1)**3
      res(2) = -U(2)/U(1)**2
      res(3) = 1./U(1)
      res=phys%epn*aux**(phys%epn-1)*res
    endif
  END SUBROUTINE compute_dAlpha_dUi

  SUBROUTINE compute_dAlpha_dUe(U, res)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: res(:)
    real*8             :: aux
    real, parameter :: tol = 1e-5
    aux = U(4)/U(1)
    if ((2./(3.*phys%Mref)*aux > 1.) .and. (switch%testcase .ne. 2)) then !! don't apply flux limiter if it is a convergence test
      res = 0.
    else
      if (aux<0) aux = tol
      res = 0.d0
      res(1) = -U(4)/U(1)**2
      res(4) = 1./U(1)
      res=phys%epn*aux**(phys%epn-1)*res
    endif
  END SUBROUTINE compute_dAlpha_dUe

  ! ******************************
  ! Parallel electric field terms
  ! ******************************
  SUBROUTINE compute_W(U, W)
    real*8, intent(IN) :: U(:)
    real*8             :: W

    W = 2./3.*U(2)/U(1)
  END SUBROUTINE compute_W

  SUBROUTINE compute_dW_dU(U, res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:, :)
    res = 0.
    res(4, 1) = -2./3.*U(2)/U(1)**2
    res(4, 2) = 2./3./U(1)
  END SUBROUTINE compute_dW_dU

  ! ******************************
  ! Temperature exchange terms
  ! ******************************
  SUBROUTINE compute_s(U, s)
    real*8, intent(IN) :: U(:)
    real*8             :: s, U1, U4, U3
    real, parameter :: tol = 1e-5
    U1 = U(1)
    U4 = U(4)
    U3 = U(3)
    if (U4 < tol) U4 = tol
    if (U1 < tol) U1 = tol
    if (U3 < tol) U3 = tol
    if (switch%ME .eq. .FALSE.) then
       if ((phys%diff_ee .gt. 0.0380) .and. (switch%testcase .ne. 2)) then
          s = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
       else
          s = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
       endif
    else
       s = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4 - U3 + 0.5*(U(2)**2/U1))
    endif
  END SUBROUTINE compute_s

  SUBROUTINE compute_ds_dU(U, res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:), U1, U4, U3
    real, parameter :: tol = 1e-5
    U1 = U(1)
    U4 = U(4)
    U3 = U(3)
    if (U4 < tol) U4 = tol
    if (U1 < tol) U1 = tol
    if (U3 < tol) U3 = tol
    res = 0.
    res(1) = 2.5*(U1/U4)**1.5*(U4 - U3 + 0.5*(U(2)**2/U1)) - 0.5*U1**0.5*U(2)**2/U4**1.5
    res(2) = U(2)*(U1/U4)**1.5
    res(3) = -U1**2.5/U4**1.5
    res(4) = -1.5*(U1/U4)**2.5*(U4 - U3 + 0.5*U(2)**2/U1) + U1**2.5/U4**1.5
    if (switch%ME .eq. .FALSE.) then
       if ((phys%diff_ee .gt. 0.0380) .and. (switch%testcase .ne. 2)) then
          res = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*res
       else
          res = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*res
       endif
    else
       res = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*res
    endif
  END SUBROUTINE compute_ds_dU


  !*****************************
  !Ohmic Heating Source
  !*****************************
  SUBROUTINE compute_Sohmic(U,Sohmic)
    real*8, intent(IN) :: U(:)
    real*8             :: Sohmic,U1,U4
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    if (U4<tol) U4=tol
    if (U1<tol) U1=tol
    Sohmic = phys%Pohmic*(((3*phys%Mref)/2)**1.5)*((U1/U4)**1.5)
  END SUBROUTINE compute_Sohmic


  SUBROUTINE compute_dSohmic_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    if (U4<tol) U4=tol
    if (U1<tol) U1=tol
    res = 0.
    res(1) = (1.5*U1**0.5)/(U4**1.5)
    res(4) = -(1.5*U1**1.5)/(U4**2.5)
    res = phys%Pohmic*(((3*phys%Mref)/2)**1.5)*res
  END SUBROUTINE compute_dSohmic_dU

  ! ******************************
  ! Neutral Source terms
  ! ******************************
#ifdef NEUTRAL

  SUBROUTINE compute_niz(U,niz)
    real*8, intent(IN) :: U(:)
    real*8             :: niz,U1,U5
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    niz = U1*U5
  END SUBROUTINE compute_niz


  SUBROUTINE compute_dniz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U5
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 0.
    res(1) = U5
    res(5) = U1
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
    real*8             :: fGammacx,U2,U5
    real,parameter :: tol = 1e-10
    U2 = U(2)
    U5 = U(5)
    if (U5<tol) U5=tol
    fGammacx = U2*U5
  END SUBROUTINE compute_fGammacx


  SUBROUTINE compute_dfGammacx_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U2,U5
    real,parameter :: tol = 1e-10
    U2 = U(2)
    U5 = U(5)
    if (U5<tol) U5=tol
    res = 0.
    res(2) = U5
    res(5) = U2
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

#ifdef NEUTRALGAMMA
    SUBROUTINE compute_fGammaN(U,fGammaN)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: fGammaN
    real*8             :: U1,U6

    fGammaN = 0.

    U1 = U(1)
    U6 = U(6)

    fGammaN = U1*U6
  END SUBROUTINE compute_fGammaN


  SUBROUTINE compute_dfGammaN_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: res(:)
    real*8             :: U1,U6

    U1 = U(1)
    U6 = U(6)

    res = 0.

    res(1) = U6
    res(6) = U1
  END SUBROUTINE compute_dfGammaN_dU
#endif

#ifdef TEMPERATURE
  !SUBROUTINE compute_sigmaviz(U,sigmaviz)
  ! real*8, intent(IN) :: U(:)
  ! real*8             :: sigmaviz,U1,U4,T0,Ery,E0
  ! real, parameter    :: tol = 1e-10
  ! U1 = U(1)
  ! U4 = U(4)
  ! T0 = 50.
  ! Ery = 13.6
  ! if (U1<tol) U1=tol
  ! if (U4<tol) U4=tol
  ! E0 = (2.*T0*U4)/(3.*phys%Mref*Ery*U1)
  ! !Threshold on Te >= 0.2 eV
  ! if (E0*Ery .le. 0.05) E0 = 0.05/Ery
  ! sigmaviz = 1.e-11*sqrt(E0)*(1./((Ery**1.5)*(6. + E0)))*exp(-1./E0)
  !END SUBROUTINE compute_sigmaviz


  !SUBROUTINE compute_dsigmaviz_dU(U,res)
  ! real*8, intent(IN) :: U(:)
  ! real*8             :: res(:),U1,U4,T0,Ery,E0
  ! real, parameter    :: tol = 1e-10
  ! U1 = U(1)
  ! U4 = U(4)
  ! T0 = 50.
  ! Ery = 13.6
  ! if (U1<tol) U1=tol
  ! if (U4<tol) U4=tol
  ! res = 0.
  ! E0 = (2.*T0*U4)/(3.*phys%Mref*Ery*U1)
  ! !Threshold on Te >= 0.2 eV
  ! !if (E0*Ery .le. 0.05)  E0 = 0.05/Ery
  ! if (E0*Ery .ge. 0.05) then
  !    res(1) = (1./U1)*( - 1./2 + 1./((6./E0) + 1)   - 1./E0 )
  !    res(4) = 1./(2.*U4) - 1./(6.*(U4/E0) + U4) + 1./(E0*U4)
  !    res = 1.e-11*sqrt(E0)*(1./((Ery**1.5)*(6. + E0)))*exp(-1./E0)*res
  ! end if
  !END SUBROUTINE compute_dsigmaviz_dU


  !SUBROUTINE compute_sigmavrec(U,sigmavrec)
  ! real*8, intent(IN) :: U(:)
  ! real*8             :: sigmavrec,U1,U4,T0,Ery,E0
  ! real, parameter    :: tol = 1e-10
  ! U1 = U(1)
  ! U4 = U(4)
  ! T0 = 50.
  ! Ery = 13.6
  ! if (U1<tol) U1=tol
  ! if (U4<tol) U4=tol
  ! E0 = (Ery*3.*phys%Mref*U1)/(T0*2.*U4)
  ! !Threshold on Te >= 0.1 eV
  ! if (E0/Ery .ge. 10.) E0 = 10.*Ery
  ! sigmavrec = 5.2e-20*sqrt(E0)*(0.43 + 0.5*log(E0) + 0.469*(E0**(-1./3)))
  !END SUBROUTINE compute_sigmavrec


  !SUBROUTINE compute_dsigmavrec_dU(U,res)
  ! real*8, intent(IN) :: U(:)
  ! real*8             :: res(:),U1,U4,T0,Ery,E0
  ! real, parameter    :: tol = 1e-10
  ! U1 = U(1)
  ! U4 = U(4)
  ! T0 = 50.
  ! Ery = 13.6
  ! if (U1<tol) U1=tol
  ! if (U4<tol) U4=tol
  ! E0 = (Ery*3.*phys%Mref*U1)/(T0*2.*U4)
  ! res = 0.
  ! !Threshold on Te >= 0.1 eV
  ! if (E0/Ery .le. 10.) then
  !    res(1) = 5.2e-20*sqrt(E0)*(0.43 + 0.5*log(E0) + 0.469*(E0**(-1./3)))/(2.*U1) + 5.2e-20*sqrt(E0)*(0.5/U1 +&
  !    &0.469*(E0**(-1./3))*(-1./(3.*U1)))
  !    res(4) = - 5.2e-20*sqrt(E0)*(0.43 + 0.5*log(E0) + 0.469*(E0**(-1./3)))/(2.*U4) + 5.2e-20*sqrt(E0)*( - 0.5/U4 +&
  !    &0.469*(E0**(-1./3))*(1./(3.*U4)))
  ! endif 
  !END SUBROUTINE compute_dsigmavrec_dU

  SUBROUTINE compute_sigmaviz(U,sigmaviz)
   real*8, intent(IN) :: U(:)
   real*8 :: sigmaviz,U1,U4,T0,Ery,E0,te,ne,n0
   real*8, dimension(9,9) :: alpha
   real, parameter :: tol = 1e-10
   integer :: i, j
   U1 = U(1)
   U4 = U(4)
   T0 = 50.
   n0 = 1.e19
   ! In EIRENE the density is scaled to 1.e14
   ne = n0*U1/1.e14
   te = T0*2/3/phys%Mref*U4/U1
   !if (U1<tol) U1=tol
   !if (U4<tol) U4=tol
   if (ne<1.1e0) ne=1.1e0
   if (ne>2.e7) ne = 2.e7
   if (te<0.1) te=0.1
   if (te>1.e4) te = 1.e4
  
   sigmaviz = 0.

   ! EIRENE documentation
   !Reaction 2.1.5JH AMJUEL
   alpha(1,:) = (/-3.29264710e+01,1.29348138e-02,5.51756251e-03,&
    -7.85381632e-04,1.43612850e-04,-3.88375028e-07,&
    -1.48977436e-06,1.41636143e-07,-3.89093208e-09/)
   alpha(2,:) = (/1.42397767e+01, -1.17314396e-02, 1.06344011e-03,&
    -1.60005353e-03, 1.13655464e-05, 5.17766228e-05,&
    -7.94799990e-06, 4.50850568e-07, -8.95261409e-09/)
   alpha(3,:) = (/-6.51943873e+00, -7.18982575e-03, 9.24737741e-04,&
    2.03702675e-03, -3.66871720e-04, 5.36863032e-06,&
    3.71395891e-06, -3.12576437e-07, 7.45121322e-09/)
   alpha(4,:) = (/2.00999615e+00, 1.27597974e-02, -4.69347962e-03,&
    -2.38922414e-05, 1.35806992e-04, -1.45489756e-05,&
    4.21203150e-08, 5.50604467e-08, -1.85267764e-09/)
   alpha(5,:) = (/-4.28959442e-01, -5.34086632e-03, 2.32458236e-03,&
    -3.21722808e-04, 6.66058141e-06, 2.39653187e-06,&
    -1.78520832e-07, 6.09564957e-10, 1.47020423e-10/)
   alpha(6,:) = (/6.04783461e-02, 9.62490059e-04, -4.18298118e-04,&
    7.95723018e-05, -7.44704256e-06, 1.84915526e-07,&
    1.61823364e-08, -8.18292830e-10, 4.83578962e-12/)
   alpha(7,:) = (/-5.30473797e-03, -7.85487245e-05, 2.73582380e-05,&
    -5.91534856e-06, 8.66630287e-07, -6.11551482e-08,&
    1.07547317e-09, 4.11800067e-11, -1.08932309e-12/)
   alpha(8,:) = (/2.60694695e-04, 2.31744225e-06, 5.14889078e-08,&
    -7.14418252e-09, -2.54019475e-08, 4.09785784e-09,&
    -2.04865734e-10, 3.02791637e-12, 1.15585402e-14/)
   alpha(9,:) = (/-5.46790307e-06, 6.07738004e-09, -4.71289307e-08,&
    1.08685876e-08, -3.44841725e-10, -8.71418322e-11,&
    8.02366070e-12, -2.39651850e-13, 2.17364528e-15/) 

   !Reaction 2.1.5 AMJUEL                                                                                            
   !alpha(1,:) = (/-3.248025330340e+01, -5.440669186583e-02, 9.048888225109e-02,&
   !      -4.054078993576e-02, 8.976513750477e-03, -1.060334011186e-03,&
   !       6.846238436472e-05, -2.242955329604e-06, 2.890437688072e-08/)
   ! alpha(2,:) = (/ 1.425332391510e+01, -3.594347160760e-02, -2.014729121556e-02,&
   !       1.039773615730e-02, -1.771792153042e-03, 1.237467264294e-04,&
   !       -3.130184159149e-06, -3.051994601527e-08, 1.888148175469e-09/)
   ! alpha(3,:) = (/-6.632235026785e+00, 9.255558353174e-02, -5.580210154625e-03,&
   !        -5.902218748238e-03, 1.295609806553e-03, -1.056721622588e-04,&
   !        4.646310029498e-06, -1.479612391848e-07, 2.852251258320e-09/)
   ! alpha(4,:) = (/ 2.059544135448e+00, -7.562462086943e-02, 1.519595967433e-02,&
   !        5.803498098354e-04, -3.527285012725e-04, 3.201533740322e-05,&
   !        -1.835196889733e-06, 9.474014343303e-08, -2.342505583774e-09/)
   ! alpha(5,:) = (/-4.425370331410e-01, 2.882634019199e-02, -7.285771485050e-03,&
   !        4.643389885987e-04, 1.145700685235e-06, 8.493662724988e-07,&
   !        -1.001032516512e-08, -1.476839184318e-08, 6.047700368169e-10/)
   ! alpha(6,:) = (/ 6.309381861496e-02, -5.788686535780e-03, 1.507382955250e-03,&
   !       -1.201550548662e-04, 6.574487543511e-06, -9.678782818849e-07,&
   !         5.176265845225e-08, 1.291551676860e-09, -9.685157340473e-11/)
   ! alpha(7,:) = (/-5.620091829261e-03, 6.329105568040e-04, -1.527777697951e-04,&
   !        8.270124691336e-06, 3.224101773605e-08, 4.377402649057e-08,&
   !        -2.622921686955e-09, -2.259663431436e-10, 1.161438990709e-11/)
   ! alpha(8,:) = (/ 2.812016578355e-04, -3.564132950345e-05, 7.222726811078e-06,&
   !        1.433018694347e-07, -1.097431215601e-07, 7.789031791949e-09,&
   !       -4.197728680251e-10, 3.032260338723e-11, -8.911076930014e-13/)
   ! alpha(9,:) = (/-6.011143453374e-06, 8.089651265488e-07, -1.186212683668e-07,&
   !        -2.381080756307e-08, 6.271173694534e-09, -5.483010244930e-10,&
   !        3.064611702159e-11, -1.355903284487e-12, 2.935080031599e-14/)

   call compute_eirene_rate(te,ne,alpha,sigmaviz)
  
  END SUBROUTINE compute_sigmaviz 
  
  SUBROUTINE compute_dsigmaviz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8 :: res(:),sigmaviz,U1,U4,T0,Ery,E0,te,ne,n0
    real*8, dimension(9,9) :: alpha
    real, parameter :: tol = 1e-10
    integer :: i, j
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19
    ne = n0*U1/1.e14
    te = T0*2/3/phys%Mref*U4/U1
    !if (U1<tol) U1=tol
    !if (U4<tol) U4=tol
    if (ne<1.1e0) ne=1.1e0
    if (ne>2.e7) ne = 2.e7
    if (te<0.1) te=0.1
    if (te>1.e4) te = 1.e4
  
    !sigmaviz = 0.
  
    !call compute_sigmaviz(U,sigmaviz)

    ! EIRENE documentation
    !Reaction 2.1.5JH AMJUEL
    alpha(1,:) = (/-3.29264710e+01,1.29348138e-02,5.51756251e-03,&
     -7.85381632e-04,1.43612850e-04,-3.88375028e-07,&
     -1.48977436e-06,1.41636143e-07,-3.89093208e-09/)
    alpha(2,:) = (/1.42397767e+01, -1.17314396e-02, 1.06344011e-03,&
     -1.60005353e-03, 1.13655464e-05, 5.17766228e-05,&
     -7.94799990e-06, 4.50850568e-07, -8.95261409e-09/)
    alpha(3,:) = (/-6.51943873e+00, -7.18982575e-03, 9.24737741e-04,&
     2.03702675e-03, -3.66871720e-04, 5.36863032e-06,&
     3.71395891e-06, -3.12576437e-07, 7.45121322e-09/)
    alpha(4,:) = (/2.00999615e+00, 1.27597974e-02, -4.69347962e-03,&
     -2.38922414e-05, 1.35806992e-04, -1.45489756e-05,&
     4.21203150e-08, 5.50604467e-08, -1.85267764e-09/)
    alpha(5,:) = (/-4.28959442e-01, -5.34086632e-03, 2.32458236e-03,&
     -3.21722808e-04, 6.66058141e-06, 2.39653187e-06,&
     -1.78520832e-07, 6.09564957e-10, 1.47020423e-10/)
    alpha(6,:) = (/6.04783461e-02, 9.62490059e-04, -4.18298118e-04,&
     7.95723018e-05, -7.44704256e-06, 1.84915526e-07,&
     1.61823364e-08, -8.18292830e-10, 4.83578962e-12/)
    alpha(7,:) = (/-5.30473797e-03, -7.85487245e-05, 2.73582380e-05,&
     -5.91534856e-06, 8.66630287e-07, -6.11551482e-08,&
     1.07547317e-09, 4.11800067e-11, -1.08932309e-12/)
    alpha(8,:) = (/2.60694695e-04, 2.31744225e-06, 5.14889078e-08,&
     -7.14418252e-09, -2.54019475e-08, 4.09785784e-09,&
     -2.04865734e-10, 3.02791637e-12, 1.15585402e-14/)
    alpha(9,:) = (/-5.46790307e-06, 6.07738004e-09, -4.71289307e-08,&
     1.08685876e-08, -3.44841725e-10, -8.71418322e-11,&
     8.02366070e-12, -2.39651850e-13, 2.17364528e-15/)
  
    !Reaction 2.1.5
    !alpha(1,:) = (/-3.248025330340e+01, -5.440669186583e-02, 9.048888225109e-02,&
    !     -4.054078993576e-02, 8.976513750477e-03, -1.060334011186e-03,&
    !      6.846238436472e-05, -2.242955329604e-06, 2.890437688072e-08/)
    !alpha(2,:) = (/ 1.425332391510e+01, -3.594347160760e-02, -2.014729121556e-02,&
    !      1.039773615730e-02, -1.771792153042e-03, 1.237467264294e-04,&
    !      -3.130184159149e-06, -3.051994601527e-08, 1.888148175469e-09/)
    !alpha(3,:) = (/-6.632235026785e+00, 9.255558353174e-02, -5.580210154625e-03,&
    !       -5.902218748238e-03, 1.295609806553e-03, -1.056721622588e-04,&
    !       4.646310029498e-06, -1.479612391848e-07, 2.852251258320e-09/)
    !alpha(4,:) = (/ 2.059544135448e+00, -7.562462086943e-02, 1.519595967433e-02,&
    !       5.803498098354e-04, -3.527285012725e-04, 3.201533740322e-05,&
    !       -1.835196889733e-06, 9.474014343303e-08, -2.342505583774e-09/)
    !alpha(5,:) = (/-4.425370331410e-01, 2.882634019199e-02, -7.285771485050e-03,&
    !       4.643389885987e-04, 1.145700685235e-06, 8.493662724988e-07,&
    !       -1.001032516512e-08, -1.476839184318e-08, 6.047700368169e-10/)
    !alpha(6,:) = (/ 6.309381861496e-02, -5.788686535780e-03, 1.507382955250e-03,&
    !      -1.201550548662e-04, 6.574487543511e-06, -9.678782818849e-07,&
    !        5.176265845225e-08, 1.291551676860e-09, -9.685157340473e-11/)
    !alpha(7,:) = (/-5.620091829261e-03, 6.329105568040e-04, -1.527777697951e-04,&
    !       8.270124691336e-06, 3.224101773605e-08, 4.377402649057e-08,&
    !       -2.622921686955e-09, -2.259663431436e-10, 1.161438990709e-11/)
    !alpha(8,:) = (/ 2.812016578355e-04, -3.564132950345e-05, 7.222726811078e-06,&
    !       1.433018694347e-07, -1.097431215601e-07, 7.789031791949e-09,&
    !      -4.197728680251e-10, 3.032260338723e-11, -8.911076930014e-13/)
    !alpha(9,:) = (/-6.011143453374e-06, 8.089651265488e-07, -1.186212683668e-07,&
    !       -2.381080756307e-08, 6.271173694534e-09, -5.483010244930e-10,&
    !       3.064611702159e-11, -1.355903284487e-12, 2.935080031599e-14/)

    res = 0.
    
    call compute_eirene_rate_du(U1,U4,te,ne,alpha,res)
  
  END SUBROUTINE compute_dsigmaviz_dU
  
    !SUBROUTINE compute_dsigmaviz_dne(U,res)
    !real*8, intent(IN) :: U(:)
    !real*8 :: res,sigmaviz,U1,U4,T0,Ery,E0,te,ne,n0
    !real*8, dimension(9,9) :: alpha
    !real, parameter :: tol = 1e-10
    !integer :: i, j
    !U1 = U(1)
    !U4 = U(4)
    !T0 = 50.
    !n0 = 1.e19
    !ne = n0*U1/1.e14
    !te = T0*2/3/phys%Mref*U4/U1
    !if (U1<tol) U1=tol
    !if (U4<tol) U4=tol
    !if (ne<1.1e0) ne=1.1e0
    !if (ne>2.e7) ne = 2.e7
    !if (te<0.1) te=0.1
    !if (te>1.e4) te = 1.e4

    !sigmaviz = 0.                                                                                                        

    !call compute_sigmaviz(U,sigmaviz)                                                                                    

    ! EIRENE documentation                                                                                                
    !alpha(1,:) = (/-3.29264710e+01,1.29348138e-02,5.51756251e-03,&
    ! -7.85381632e-04,1.43612850e-04,-3.88375028e-07,&
    ! -1.48977436e-06,1.41636143e-07,-3.89093208e-09/)
    !alpha(2,:) = (/1.42397767e+01, -1.17314396e-02, 1.06344011e-03,&
    ! -1.60005353e-03, 1.13655464e-05, 5.17766228e-05,&
    ! -7.94799990e-06, 4.50850568e-07, -8.95261409e-09/)
    !alpha(3,:) = (/-6.51943873e+00, -7.18982575e-03, 9.24737741e-04,&
    ! 2.03702675e-03, -3.66871720e-04, 5.36863032e-06,&
    ! 3.71395891e-06, -3.12576437e-07, 7.45121322e-09/)
    !alpha(4,:) = (/2.00999615e+00, 1.27597974e-02, -4.69347962e-03,&
    ! -2.38922414e-05, 1.35806992e-04, -1.45489756e-05,&
    ! 4.21203150e-08, 5.50604467e-08, -1.85267764e-09/)
    !alpha(5,:) = (/-4.28959442e-01, -5.34086632e-03, 2.32458236e-03,&
    ! -3.21722808e-04, 6.66058141e-06, 2.39653187e-06,&
    ! -1.78520832e-07, 6.09564957e-10, 1.47020423e-10/)
    !alpha(6,:) = (/6.04783461e-02, 9.62490059e-04, -4.18298118e-04,&
    ! 7.95723018e-05, -7.44704256e-06, 1.84915526e-07,&
    ! 1.61823364e-08, -8.18292830e-10, 4.83578962e-12/)
    !alpha(7,:) = (/-5.30473797e-03, -7.85487245e-05, 2.73582380e-05,&
    ! -5.91534856e-06, 8.66630287e-07, -6.11551482e-08,&
    ! 1.07547317e-09, 4.11800067e-11, -1.08932309e-12/)
    !alpha(8,:) = (/2.60694695e-04, 2.31744225e-06, 5.14889078e-08,&
    ! -7.14418252e-09, -2.54019475e-08, 4.09785784e-09,&
    ! -2.04865734e-10, 3.02791637e-12, 1.15585402e-14/)
    !alpha(9,:) = (/-5.46790307e-06, 6.07738004e-09, -4.71289307e-08,&
    ! 1.08685876e-08, -3.44841725e-10, -8.71418322e-11,&
    ! 8.02366070e-12, -2.39651850e-13, 2.17364528e-15/)

    !res = 0.

    !call compute_eirene_rate_dne(U1,U4,te,ne,alpha,res)

  !END SUBROUTINE compute_dsigmaviz_dne  

  SUBROUTINE compute_sigmavrec(U,sigmavrec)
    real*8, intent(IN) :: U(:)
    real*8 :: sigmavrec,U1,U4,T0,Ery,E0,te,ne,n0
    real*8, dimension(9,9) :: alpha
    real, parameter :: tol = 1e-10
    integer :: i, j
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19
    ne = n0*U1/1.e14
    te = T0*2/3/phys%Mref*U4/U1
    !if (U1<tol) U1=tol
    !if (U4<tol) U4=tol
    if (ne<1.1e0) ne=1.1e0
    if (ne>2.e7) ne = 2.e7
    if (te<0.1) te=0.1
    if (te>1.e4) te = 1.e4
  
    sigmavrec = 0.
  
    ! EIRENE documentation
    !Reaction 2.1.8JH
    alpha(1,:) = (/-2.85572848e+01, 3.48856323e-02, -2.79964439e-02,&
     1.20954532e-02, -2.43663080e-03, 2.83789372e-04,&
     -1.88651117e-05, 6.75215560e-07, -1.00589386e-08/)
    alpha(2,:) = (/-7.66404261e-01, -3.58323337e-03, -7.45251429e-03,&
     2.70929976e-03, -7.74512977e-04, 1.14244470e-04,&
     -9.38278352e-06, 3.90280010e-07, -6.38741159e-09/)
    alpha(3,:) = (/-4.93042400e-03, -3.62024535e-03, 6.95871196e-03,&
     -2.13925730e-03, 4.60388371e-04, -5.99163684e-05,&
     4.72926255e-06, -1.99348540e-07, 3.35258987e-09/)
    alpha(4,:) = (/-5.38683098e-03, -9.53284048e-04, 4.63175381e-04,&
     -5.37117970e-04, 1.54335050e-04, -2.25756584e-05,&
     1.73078295e-06, -6.61824078e-08, 1.01336428e-09/)
    alpha(5,:) = (/-1.62603924e-04, 1.88804863e-04, 1.28857769e-04,&
     -1.63458052e-05, -9.60103695e-06, 3.42526239e-06,&
     -4.07701994e-07, 2.04204110e-08, -3.70797772e-10/)
    alpha(6,:) = (/6.08090765e-06, -1.01489068e-05, -1.14502889e-04,&
     5.94219398e-05, -1.21185172e-05, 1.11896550e-06,&
     -4.27532157e-08, 3.70861611e-10, 7.06845011e-12/)
    alpha(7,:) = (/2.10110205e-05, 2.24567656e-05, -2.24562427e-06,&
     -2.94487376e-06, 1.00210510e-06, -1.29132080e-07,&
     7.78615546e-09, -2.44112778e-10, 3.77320848e-12/)
    alpha(8,:) = (/-2.77071760e-06, -4.69598237e-06, 3.25087887e-06,&
     -9.38729079e-07, 1.39239163e-07, -1.13909329e-08,&
     5.17850560e-10, -9.45240216e-12, -4.67272402e-14/)
    alpha(9,:) = (/1.03823594e-07, 2.52316661e-07, -2.14539040e-07,&
     7.38143524e-08, -1.29971368e-08, 1.26518958e-09,&
     -6.85420397e-11, 1.83661503e-12, -1.64049236e-14/)

    !Reaction 2.1.8a AMJUEL
    !alpha(1,:) = (/-2.861779556590e+01, -1.786166918005e-02, 6.391553337864e-04,&     
    ! -4.509415260040e-04, 7.095459017274e-05, -5.660309928918e-06,&
    ! 1.160186631232e-07, 7.564986067995e-09, -2.969815025786e-10/)
    !alpha(2,:) = (/-7.251997071478e-01, 3.210966054964e-03, 4.550251497787e-03,&
    !  -1.882306456891e-03, 3.983133042462e-04, -4.851835293564e-05,&
    !  3.404834497087e-06, -1.280839994482e-07, 1.982839967575e-09/)
    !alpha(3,:) = (/ -1.735023322687e-02, -3.112517426840e-03, 1.077863345492e-03,&
    !  -2.616958968739e-04, 5.459332810644e-05, -8.635308675130e-06,&
    !  8.383106368091e-07, -4.133352004945e-08, 7.872491728981e-10/)
    !alpha(4,:) = (/-3.557752804131e-03, 1.558966107388e-03, -1.037331531958e-03,&
    !  2.817237174744e-04, -4.407815167942e-05, 4.646017350681e-06,&
    !  -3.365654551356e-07, 1.428350791171e-08, -2.522153346435e-10/)
    !alpha(5,:) = (/-2.777882255016e-04, -9.329932857673e-05, 1.096331766957e-04,&
    !  -4.567488387292e-05, 8.495787235165e-06, -7.261076273040e-07,&
    !   2.326992940046e-08, 2.208089550616e-10, -1.989979386039e-11/)
    !alpha(6,:) = (/ 2.060295404466e-05, -1.283711654633e-04, 7.312311894769e-05,&
    !   -1.064805149480e-05, -1.498776433806e-07, 1.199087596048e-07,&
    !   -5.668079133507e-09, -1.018554043516e-10, 7.766578964142e-12/)
    !alpha(7,:) = (/ 1.593238392469e-05, 3.705503401064e-05, -2.407235857913e-05,&
    !   4.915213917257e-06, -3.346609397503e-07, -4.912753691671e-09,&
    !   1.302393677822e-09, -3.169013613822e-11, -1.783762758524e-13/)
    !alpha(8,:) = (/-2.116580756634e-06, -3.854172456142e-06, 2.662392026941e-06,&
    !   -6.120846201882e-07, 5.663728215333e-08, -1.474221162308e-09,&
    !   -7.373095178045e-11, 4.314457229158e-12, -4.791677504810e-14/)
    !alpha(9,:) = (/7.665990100168e-08, 1.400789118322e-07, -1.008951470934e-07,&
    !    2.495214914834e-08, -2.678484130657e-09, 1.170138331019e-10,&
    !    -1.588254701759e-13, -1.226345218681e-13, 2.329402447113e-15/)
  
     call compute_eirene_rate(te,ne,alpha,sigmavrec)
  
  END SUBROUTINE compute_sigmavrec

  SUBROUTINE compute_dsigmavrec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8 :: res(:),sigmavrec,U1,U4,T0,Ery,E0,te,ne,n0
    real*8, dimension(9,9) :: alpha
    real, parameter :: tol = 1e-10
    integer :: i, j
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19
    ne = n0*U1/1.e14
    te = T0*2/3/phys%Mref*U4/U1
    !if (U1<tol) U1=tol
    !if (U4<tol) U4=tol
    if (ne<1.1e0) ne=1.1e0
    if (ne>2.e7) ne = 2.e7
    if (te<0.1) te=0.1
    if (te>1.e4) te = 1.e4
  
    !sigmavrec = 0.
  
    !call compute_sigmavrec(U,sigmavrec)
  
    ! EIRENE documentation
    !Reaction 2.1.8JH AMJUEL
    alpha(1,:) = (/-2.85572848e+01, 3.48856323e-02, -2.79964439e-02,&
     1.20954532e-02, -2.43663080e-03, 2.83789372e-04,&
     -1.88651117e-05, 6.75215560e-07, -1.00589386e-08/)
    alpha(2,:) = (/-7.66404261e-01, -3.58323337e-03, -7.45251429e-03,&
     2.70929976e-03, -7.74512977e-04, 1.14244470e-04,&
     -9.38278352e-06, 3.90280010e-07, -6.38741159e-09/)
    alpha(3,:) = (/-4.93042400e-03, -3.62024535e-03, 6.95871196e-03,&
     -2.13925730e-03, 4.60388371e-04, -5.99163684e-05,&
     4.72926255e-06, -1.99348540e-07, 3.35258987e-09/)
    alpha(4,:) = (/-5.38683098e-03, -9.53284048e-04, 4.63175381e-04,&
     -5.37117970e-04, 1.54335050e-04, -2.25756584e-05,&
     1.73078295e-06, -6.61824078e-08, 1.01336428e-09/)
    alpha(5,:) = (/-1.62603924e-04, 1.88804863e-04, 1.28857769e-04,&
     -1.63458052e-05, -9.60103695e-06, 3.42526239e-06,&
     -4.07701994e-07, 2.04204110e-08, -3.70797772e-10/)
    alpha(6,:) = (/6.08090765e-06, -1.01489068e-05, -1.14502889e-04,&
     5.94219398e-05, -1.21185172e-05, 1.11896550e-06,&
     -4.27532157e-08, 3.70861611e-10, 7.06845011e-12/)
    alpha(7,:) = (/2.10110205e-05, 2.24567656e-05, -2.24562427e-06,&
     -2.94487376e-06, 1.00210510e-06, -1.29132080e-07,&
     7.78615546e-09, -2.44112778e-10, 3.77320848e-12/)
    alpha(8,:) = (/-2.77071760e-06, -4.69598237e-06, 3.25087887e-06,&
     -9.38729079e-07, 1.39239163e-07, -1.13909329e-08,&
     5.17850560e-10, -9.45240216e-12, -4.67272402e-14/)
    alpha(9,:) = (/1.03823594e-07, 2.52316661e-07, -2.14539040e-07,&
     7.38143524e-08, -1.29971368e-08, 1.26518958e-09,&
     -6.85420397e-11, 1.83661503e-12, -1.64049236e-14/)

    !Reaction 2.1.8a AMJUEL 
    !alpha(1,:) = (/-2.861779556590e+01, -1.786166918005e-02, 6.391553337864e-04,&
    ! -4.509415260040e-04, 7.095459017274e-05, -5.660309928918e-06,&
    ! 1.160186631232e-07, 7.564986067995e-09, -2.969815025786e-10/)
    !alpha(2,:) = (/-7.251997071478e-01, 3.210966054964e-03, 4.550251497787e-03,&
    !  -1.882306456891e-03, 3.983133042462e-04, -4.851835293564e-05,&
    !  3.404834497087e-06, -1.280839994482e-07, 1.982839967575e-09/)
    !alpha(3,:) = (/ -1.735023322687e-02, -3.112517426840e-03, 1.077863345492e-03,&
    !  -2.616958968739e-04, 5.459332810644e-05, -8.635308675130e-06,&
    !  8.383106368091e-07, -4.133352004945e-08, 7.872491728981e-10/)
    !alpha(4,:) = (/-3.557752804131e-03, 1.558966107388e-03, -1.037331531958e-03,&
    !  2.817237174744e-04, -4.407815167942e-05, 4.646017350681e-06,&
    !  -3.365654551356e-07, 1.428350791171e-08, -2.522153346435e-10/)
    !alpha(5,:) = (/-2.777882255016e-04, -9.329932857673e-05, 1.096331766957e-04,&
    !  -4.567488387292e-05, 8.495787235165e-06, -7.261076273040e-07,&
    !   2.326992940046e-08, 2.208089550616e-10, -1.989979386039e-11/)
    !alpha(6,:) = (/ 2.060295404466e-05, -1.283711654633e-04, 7.312311894769e-05,&
    !   -1.064805149480e-05, -1.498776433806e-07, 1.199087596048e-07,&
    !   -5.668079133507e-09, -1.018554043516e-10, 7.766578964142e-12/)
    !alpha(7,:) = (/ 1.593238392469e-05, 3.705503401064e-05, -2.407235857913e-05,&
    !   4.915213917257e-06, -3.346609397503e-07, -4.912753691671e-09,&
    !   1.302393677822e-09, -3.169013613822e-11, -1.783762758524e-13/)
    !alpha(8,:) = (/-2.116580756634e-06, -3.854172456142e-06, 2.662392026941e-06,&
    !   -6.120846201882e-07, 5.663728215333e-08, -1.474221162308e-09,&
    !   -7.373095178045e-11, 4.314457229158e-12, -4.791677504810e-14/)
    !alpha(9,:) = (/7.665990100168e-08, 1.400789118322e-07, -1.008951470934e-07,&
    !    2.495214914834e-08, -2.678484130657e-09, 1.170138331019e-10,&
    !    -1.588254701759e-13, -1.226345218681e-13, 2.329402447113e-15/)
  
     res = 0.
  
    call compute_eirene_rate_du(U1,U4,te,ne,alpha,res)
  
  END SUBROUTINE compute_dsigmavrec_dU

  !SUBROUTINE compute_dsigmavrec_dne(U,res)
  !  real*8, intent(IN) :: U(:)
  !  real*8 :: res,sigmavrec,U1,U4,T0,Ery,E0,te,ne,n0
  !  real*8, dimension(9,9) :: alpha
  !  real, parameter :: tol = 1e-10
  !  integer :: i, j
  !  U1 = U(1)
  !  U4 = U(4)
  !  T0 = 50.
  !  n0 = 1.e19
  !  ne = n0*U1/1.e14
  !  te = T0*2/3/phys%Mref*U4/U1
  !  if (U1<tol) U1=tol
  !  if (U4<tol) U4=tol
  !  if (ne<1.1e0) ne=1.1e0
  !  if (ne>2.e7) ne = 2.e7
  !  if (te<0.1) te=0.1
  !  if (te>1.e4) te = 1.e4

    !sigmavrec = 0.                                                                                                       

    !call compute_sigmavrec(U,sigmavrec)                                                                                  

    ! EIRENE documentation                                                                                                
  !  alpha(1,:) = (/-2.85572848e+01, 3.48856323e-02, -2.79964439e-02,&
  !   1.20954532e-02, -2.43663080e-03, 2.83789372e-04,&
  !   -1.88651117e-05, 6.75215560e-07, -1.00589386e-08/)
  !  alpha(2,:) = (/-7.66404261e-01, -3.58323337e-03, -7.45251429e-03,&
  !   2.70929976e-03, -7.74512977e-04, 1.14244470e-04,&
  !   -9.38278352e-06, 3.90280010e-07, -6.38741159e-09/)
  !  alpha(3,:) = (/-4.93042400e-03, -3.62024535e-03, 6.95871196e-03,&
  !   -2.13925730e-03, 4.60388371e-04, -5.99163684e-05,&
  !   4.72926255e-06, -1.99348540e-07, 3.35258987e-09/)
  !  alpha(4,:) = (/-5.38683098e-03, -9.53284048e-04, 4.63175381e-04,&
  !   -5.37117970e-04, 1.54335050e-04, -2.25756584e-05,&
  !   1.73078295e-06, -6.61824078e-08, 1.01336428e-09/)
  !  alpha(5,:) = (/-1.62603924e-04, 1.88804863e-04, 1.28857769e-04,&
  !   -1.63458052e-05, -9.60103695e-06, 3.42526239e-06,&
  !   -4.07701994e-07, 2.04204110e-08, -3.70797772e-10/)
  !  alpha(6,:) = (/6.08090765e-06, -1.01489068e-05, -1.14502889e-04,&
  !   5.94219398e-05, -1.21185172e-05, 1.11896550e-06,&
  !   -4.27532157e-08, 3.70861611e-10, 7.06845011e-12/)
  !  alpha(7,:) = (/2.10110205e-05, 2.24567656e-05, -2.24562427e-06,&
  !   -2.94487376e-06, 1.00210510e-06, -1.29132080e-07,&
  !   7.78615546e-09, -2.44112778e-10, 3.77320848e-12/)
  !  alpha(8,:) = (/-2.77071760e-06, -4.69598237e-06, 3.25087887e-06,&
  !   -9.38729079e-07, 1.39239163e-07, -1.13909329e-08,&
  !   5.17850560e-10, -9.45240216e-12, -4.67272402e-14/)
  !  alpha(9,:) = (/1.03823594e-07, 2.52316661e-07, -2.14539040e-07,&
  !   7.38143524e-08, -1.29971368e-08, 1.26518958e-09,&
  !   -6.85420397e-11, 1.83661503e-12, -1.64049236e-14/)

  !   res = 0.

  !  call compute_eirene_rate_dne(U1,U4,te,ne,alpha,res)

  !END SUBROUTINE compute_dsigmavrec_dne


  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavcx,U1,U4,T0,E0
    real*8              :: p1,p2,p3,p4,p5
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    !if (U1<tol) U1=tol
    !if (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand
    !E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    !sigmavcx = (2.5e-15/exp(-0.5))*exp(-E0)
    !Fit log log from ADAS database
    p1 = -0.0006837
    p2 = 0.008004
    p3 = -0.03581
    p4 = 0.4518
    p5 = -32.59
    E0 = T0*2./(3.*phys%Mref)*U4/U1

    !Threshold on Te >= 0.2 eV
    if (E0 .le. 0.05) E0 = 0.05
    sigmavcx = exp(p1*log(E0)**4 + p2*log(E0)**3 + p3*log(E0)**2 + p4*log(E0) + p5)   
  END SUBROUTINE compute_sigmavcx


  SUBROUTINE compute_dsigmavcx_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,T0,E0
    real*8             :: p1,p2,p3,p4,p5
    real, parameter    :: tol = 1e-10
    T0 = 50.
    U1 = U(1)
    U4 = U(4)
    if (U1<tol) U1=tol
    if (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand
    !E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    !res = 0.
    !res(1) = -1./U4
    !res(4) = U1/(U4**2)
    !res = (1.5*phys%Mref/(2.*T0))*(2.5e-15/exp(-0.5))*exp(-E0)*res
    !Fit log log from ADAS database
    p1 = -0.0006837
    p2 = 0.008004
    p3 = -0.03581
    p4 = 0.4518
    p5 = -32.59
    E0 = T0*2./(3.*phys%Mref)*U4/U1
    res = 0.
    !Threshold on Te >= 0.2 eV
    !if (E0 .le. 0.05) E0 = 0.05
    if (E0 .ge. 0.05) then
       res(1) = (4.*p1*log(E0)**3 + 3.*p2*log(E0)**2 + 2.*p3*log(E0) + p4)*(-U4/U1**2) 
       res(4) = (4.*p1*log(E0)**3 + 3.*p2*log(E0)**2 + 2.*p3*log(E0) + p4)*1./U1
       res = exp(p1*log(E0)**4 + p2*log(E0)**3 + p3*log(E0)**2 + p4*log(E0) + p5)*U1/U4*res
    end if
  END SUBROUTINE compute_dsigmavcx_dU

  SUBROUTINE compute_eirene_rate(te,ne,alpha,rate)
    real*8, intent(IN) :: te,ne,alpha(:,:) 
    real*8, intent(OUT):: rate
    integer :: i,j 
 
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    do i=1,size(alpha,1)
       do j = 1,size(alpha,2) 
          rate = rate + alpha(i,j)*log(ne)**(j-1)*log(te)**(i-1)
      end do
    end do
    ! rate is in cm^3/s in EIRENE
    rate = exp(rate)/1.e6
  
  END SUBROUTINE compute_eirene_rate
  
  SUBROUTINE compute_eirene_rate_du(U1,U4,te,ne,alpha,rate_du)
    real*8, intent(IN) :: U1,U4,te,ne,alpha(:,:)
    real*8 :: rate
    real*8, intent(OUT):: rate_du(:)
    integer :: i,j
  
    call compute_eirene_rate(te,ne,alpha,rate)
  
    ! In EIRENE the density is scaled to 1.e14
    rate_du = 0.
    
    do i=1,size(alpha,1)
       do j = 1,size(alpha,2)
          if (ne > 1.1e0 .and.  ne < 2.e7) then ! Limit on ne
             rate_du(1) = rate_du(1)+alpha(i,j)*log(ne)**(j-2)*(j-1)*log(te)**(i-1)*(1./U1)
          end if
          if (te > 0.1 .and. te < 1.e4) then ! Limit on Te
             rate_du(1) = rate_du(1)+alpha(i,j)*log(ne)**(j-1)*(i-1)*log(te)**(i-2)*(-1./U1)
             rate_du(4) = rate_du(4)+alpha(i,j)*log(ne)**(j-1)*(i-1)*log(te)**(i-2)*(1./U4)
          end if      
       end do
    end do 
    rate_du = rate_du*rate
  
  END SUBROUTINE compute_eirene_rate_du  

  !SUBROUTINE compute_eirene_rate_dne(U1,U4,te,ne,alpha,rate_du)
  !  real*8, intent(IN) :: U1,U4,te,ne,alpha(:,:)
  !  real*8 :: rate
  !  real*8, intent(OUT):: rate_du
  !  integer :: i,j

  !  call compute_eirene_rate(te,ne,alpha,rate)

    ! In EIRENE the density is scaled to 1.e14                                                                            
  !  rate_du = 0.

  !  do i=1,size(alpha,1)
  !     do j = 1,size(alpha,2)
  !        if (ne > 1.1e0 .and.  ne < 2.e7) then ! Limit on ne                                        !                   
  !           rate_du = rate_du + alpha(i,j)*log(ne)**(j-2)*(j-1)*log(te)**(i-1)*(1./U1)
  !        end if
  !     end do
  !  end do
  !  rate_du = rate_du*rate*U1

  !END SUBROUTINE compute_eirene_rate_dne
  

  SUBROUTINE compute_Tloss(U,Tloss)
    real*8, intent(IN) :: U(:)
    real*8             :: Tloss,U1,U4,T0
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    if (U1<tol) U1=tol
    if (U4<tol) U4=tol
    Tloss = 25. + 170.*exp(-(T0*U4)/(3.*phys%Mref*U1))
  END SUBROUTINE compute_Tloss


  SUBROUTINE compute_dTloss_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,T0
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    if (U1<tol) U1=tol
    if (U4<tol) U4=tol
    res = 0.
    res(1) = U4/(U1**2)
    res(4) = -1./U1
    res = 170.*exp(-(T0*U4)/(3.*phys%Mref*U1))*(T0/(3.*phys%Mref))*res
  END SUBROUTINE compute_dTloss_dU


  SUBROUTINE compute_Tlossrec(U,Tlossrec)
    real*8, intent(IN) :: U(:)
    real*8             :: Tlossrec,U1,U4,T0
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    if (U1<tol) U1=tol
    if (U4<tol) U4=tol
    !Tlossrec = min(250.,8.*exp((2.*U4*T0)/(3.*phys%Mref*U1*9.)))
    Tlossrec = (2.*U4*T0)/(3.*phys%Mref*U1*9.)
    if (Tlossrec .le. log(250./8.)) then
       Tlossrec = 8.*exp(Tlossrec)
    else
       Tlossrec = 250.
    endif
  END SUBROUTINE compute_Tlossrec


  SUBROUTINE compute_dTlossrec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,T0,Tlossrec
    real, parameter    :: tol = 1e-10
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    if (U1<tol) U1=tol
    if (U4<tol) U4=tol
    res = 0.
    !Tlossrec = 8.*exp((2.*U4*T0)/(3.*phys%Mref*U1*9.))
    !if (Tlossrec .le. 250.) then
    !  res(1) = -U4/(U1**2)
    !  res(4) = 1./U1
    !  res = Tlossrec*((2.*T0)/(27.*phys%Mref))*res
    !endif
    Tlossrec = (2.*U4*T0)/(3.*phys%Mref*U1*9.)
    if (Tlossrec .le. log(250./8.)) then
       res(1) = -U4/(U1**2)
       res(4) = 1./U1
       res = 8.*exp(Tlossrec)*((2.*T0)/(27.*phys%Mref))*res
    endif 
  END SUBROUTINE compute_dTlossrec_dU


  SUBROUTINE compute_fEiiz(U,fEiiz)
    real*8, intent(IN) :: U(:)
    real*8             :: fEiiz,U3,U5
    real,parameter :: tol = 1e-10
    U3 = U(3)
    U5 = U(5)
    if (U3<tol) U3=tol
    if (U5<tol) U5=tol
    fEiiz = U3*U5
  END SUBROUTINE compute_fEiiz


  SUBROUTINE compute_dfEiiz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U3,U5
    real,parameter :: tol = 1e-10
    U3 = U(3)
    U5 = U(5)
    if (U3<tol) U3=tol
    if (U5<tol) U5=tol
    res = 0.
    res(3) = U5
    res(5) = U3
  END SUBROUTINE compute_dfEiiz_dU


  SUBROUTINE compute_fEirec(U,fEirec)
    real*8, intent(IN) :: U(:)
    real*8             :: fEirec,U1,U3
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U3 = U(3)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    fEirec = U1*U3
  END SUBROUTINE compute_fEirec


  SUBROUTINE compute_dfEirec_dU(U,res)
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
  END SUBROUTINE compute_dfEirec_dU


  SUBROUTINE compute_fEicx(U,fEicx)
    real*8, intent(IN) :: U(:)
    real*8             :: fEicx,U1,U2,U5
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    fEicx = (U5*U2**2)/U1
  END SUBROUTINE compute_fEicx


  SUBROUTINE compute_dfEicx_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U2,U5
    real,parameter :: tol = 1e-10
    U1 = U(1)
    U2 = U(2)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 0.
    res(1) = -U5*(U2/U1)**2
    res(2) = 2.*U5*U2/U1
    res(5) = (U2**2)/U1
  END SUBROUTINE compute_dfEicx_dU
  
#ifdef NEUTRALP
  SUBROUTINE computeDpn(U,Q,Vpn,Dpn)
    real*8, intent(IN) :: U(:),Q(:,:),Vpn(:)
	   real*8, intent(OUT):: Dpn
    real*8             :: U1,U2,U3,U4,U5
    real*8             ::sigmaviz,sigmavcx,cs_n,Dpn_th
    real*8             :: Grad_Pn(simpar%Ndim)
    real,parameter :: tol = 1e-10
	   Dpn = 0.
	   U1 = U(1)
	   U2 = U(2)
	   U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U2<tol) U2=tol
    if (U3<tol) U3=tol
    if (U4<tol) U4=tol
    if (U5<tol) U5=tol
    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavcx(U,sigmavcx)
    Dpn = 1./(simpar%refval_time*simpar%refval_density*U1*(sigmaviz + sigmavcx))
    Dpn = 2./3.*Dpn   
    !Set a threshold Dpn*|grad(Pn)| <= cs_n*n_n
    !cs_n = 0.
    !Grad_Pn = 0. 
    !cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    !Grad_Pn = 2./(3.*phys%Mref)*matmul(Q,Vpn)
    !Dpn_th = abs(cs_n*U5/norm2(Grad_Pn))
    !if (Dpn .gt. Dpn_th) then 
    if (Dpn .gt. phys%diff_nn) then
	Dpn = phys%diff_nn 
	!fth(3) = 5.e-7*abs(U5)!*abs(U3/U1)
	!fth(4) = 5.e-7*abs(U5)!*abs(U4/U1)
    end if
    !Dpn*exp(abs(Dpn - Dpn_th)/Dpn_th) + Dpn_th 
    END SUBROUTINE computeDpn
	 
    SUBROUTINE compute_dDpn_dU(U,Q,Vpn,res)
    real*8, intent(IN) :: U(:),Q(:,:),Vpn(:)
    real*8, intent(OUT):: res(:)
    real*8             :: sigmaviz,sigmavcx,Dpn
    real*8             :: dsigmaviz_dU(phys%Neq),dsigmavcx_dU(phys%Neq)
    real*8             :: U1,U2,U3,U4,U5
    real,parameter     :: tol = 1e-12
	   U1 = U(1)
    U2 = U(2)
	   U3 = U(3)
	   U4 = U(4)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 0.
    
    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavcx(U,sigmavcx)
    call compute_dsigmaviz_dU(U,dsigmaviz_dU)
    call compute_dsigmavcx_dU(U,dsigmavcx_dU)
    call computeDpn(U,Q,Vpn,Dpn)
    
    IF (Dpn .lt. phys%diff_nn) THEN
       res(1) = sigmaviz + sigmavcx + U1*(dsigmaviz_dU(1) + dsigmavcx_dU(1))
       res(4) = U1*(dsigmaviz_dU(4) + dsigmavcx_dU(4))
       res = -3./2*Dpn**2*res
    END IF
  END SUBROUTINE
  
  SUBROUTINE computeVpn(U,Vpn)
    real*8, intent(IN) :: U(:)
	   real*8, intent(OUT):: Vpn(:)
    real*8             :: U1,U2,U3,U5
    real,parameter :: tol = 1e-12
	   Vpn = 0.
    U1 = U(1)
    U2 = U(2)
	   U3 = U(3)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    Vpn(1) = U5*(U2**2/U1**3 - U3/U1**2)
	   Vpn(2) = - U5*U2/U1**2
	   Vpn(3) = U5/U1
	   Vpn(5) = U3/U1 - 1./2.*U2**2/U1**2
  END SUBROUTINE computeVpn
  
  SUBROUTINE compute_dVpn_dU(U,dVpn_dU)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: dVpn_dU(:, :)
	   real*8              :: U1,U2,U3,U5
    real,parameter :: tol = 1e-12
	   U1 = U(1)
    U2 = U(2)
	   U3 = U(3)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    dVpn_dU = 0.
	
    dVpn_dU(1, 1) = (2.*U3/U1**3 - 3.*U2**2/U1**4)*U5
    dVpn_dU(1, 2) = (2.*U2/U1**3)*U5
	   dVpn_dU(1, 3) = -U5/U1**2 
   	dVpn_dU(1, 5) = U2**2/U1**3 - U3/U1**2
	
	   dVpn_dU(2, 1) = (2.*U2/U1**3)*U5
	   dVpn_dU(2, 2) = - U5/U1**2
	   dVpn_dU(2, 5) = - U2/U1**2
	
	   dVpn_dU(3, 1) = -U5/U1**2 
	   dVpn_dU(3, 5) = 1./U1
	
	   dVpn_dU(5, 1) = U2**2/U1**3 - U3/U1**2
    dVpn_dU(5, 2) = -U2/U1**2
    dVpn_dU(5, 3) = 1./U1 	
  END SUBROUTINE compute_dVpn_dU
  
  SUBROUTINE computeGammared(U,res)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: res
	   real*8              :: U1,U2,U3,U5,Tmin,T
    real,parameter      :: tol = 1e-12
	   U1 = U(1)
    U2 = U(2)
	   U3 = U(3)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    
    res = 0.
    Tmin = 0.2/simpar%refval_temperature  ! Threshold at 0.2 eV
    T = 2./(3.*phys%Mref)*(U3/U1 - 1./2.*U2**2/U1**2)
    
    !WRITE(6,*) 'Tmin/T = ', Tmin/T
    
    IF (Tmin/T .le. 1) THEN
       res = 0.*3./2.*(phys%Mref**2)*(Tmin/T)*U5
    ELSE
       res = 0.*3./2.*(phys%Mref**2)*U5
    END IF
    
  END SUBROUTINE
  
SUBROUTINE computeAlphaCoeff(U,Q,Vpn,res)
    real*8, intent(IN)    :: U(:),Q(:,:),Vpn(:)
    real*8, intent(OUT)   :: res
    real*8                :: U1,U2,U3,U4,U5
    real*8                :: cs_n,Dpn,t
    real*8                :: Grad_Pn(simpar%Ndim)
    real,parameter        :: gamma = 4.,tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 1.

    cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*matmul(Q,Vpn)
    call computeDpn(U,Q,Vpn,Dpn)
    !t = abs(Dpn*norm2(Grad_Pn))/abs(cs_n*U5)
    t = abs(Dpn)/phys%diff_nn
    
    IF (t .ge. 1) THEN
       res = exp((1-t)/gamma)
    END IF
  END SUBROUTINE
  
  SUBROUTINE computeBetaCoeff(U,Q,Vpn,res)
    real*8, intent(IN)    :: U(:),Q(:,:),Vpn(:)
    real*8, intent(OUT)   :: res
    real*8                :: U1,U2,U3,U4,U5
    real*8                :: cs_n,Dpn,t
    real*8                :: Grad_Pn(simpar%Ndim)
    real,parameter        :: gamma = 4.,tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 0.

    cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*matmul(Q,Vpn)
    call computeDpn(U,Q,Vpn,Dpn)
    !t = abs(Dpn*norm2(Grad_Pn))/abs(cs_n*U5)
    t = abs(Dpn)/phys%diff_nn
    
    IF (t .ge. 1) THEN
       res = 1 - exp((1-t)/gamma)
    END IF
  END SUBROUTINE

  SUBROUTINE computeGammaLim(U,Q,Vpn,res)
    real*8, intent(IN)    :: U(:),Q(:,:),Vpn(:)
    real*8, intent(OUT)   :: res
    real*8                :: U1,U2,U3,U4,U5
    real*8                :: cs_n,Dpn,GammaDpn,GammaLim
    real*8                :: Grad_Pn(simpar%Ndim)
    real,parameter        :: tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    if (U1<tol) U1=tol
    if (U5<tol) U5=tol
    res = 0.

    cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*matmul(Q,Vpn)
    call computeDpn(U,Q,Vpn,Dpn)
    GammaDpn = abs(Dpn*norm2(Grad_Pn))
    GammaLim = abs(cs_n*U5)

    res = 1./(1. + GammaDpn/GammaLim)
  END SUBROUTINE
#endif 
!NEUTRAL PRESSURE

! NEUTRAL GAMMA Section
#ifdef NEUTRALGAMMA
  !**************************************** 
  ! Curvature term matrix
  !****************************************
  SUBROUTINE GimpMatrixN(U, divb, Gn)
    real*8, intent(in)  :: U(:), divb
    real*8, intent(out) :: Gn(:, :)
    real*8              :: U1,U2,U3,U5
    
    Gn = 0.d0
   
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U5 = U(5)

    Gn(6, 1) = 2./3.*U5*(- U3/U1**2 + U2**2/U1**3)
    Gn(6, 2) = - 2./3.*U5*U2/U1**2
    Gn(6, 3) = 2./3.*U5/U1
    Gn(6, 5) = 2./3.*(U3/U1 - 1./2.*U2**2/U1**2) 

    Gn = divb*Gn

  END SUBROUTINE GimpMatrixN
  
  SUBROUTINE computeEtan(U,Etan)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: Etan
    real*8             :: U1,U2,U3,U4,U5,U6
    real*8             :: Tn,sigmavcx
    
    Etan = 0.

    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    U6 = U(6)

    CALL compute_sigmavcx(U,sigmavcx)
    
    Tn = simpar%refval_temperature*2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2) 
    Etan = simpar%refval_density*U5*simpar%refval_charge*max(Tn,0.1)/(simpar%refval_mass*simpar%refval_density*U1*sigmavcx)
    Etan = Etan*simpar%refval_time/simpar%refval_length**2/simpar%refval_density

    if (Etan/U5 .gt. phys%diff_nn) then
        Etan = U5*phys%diff_nn
    end if
  
  END SUBROUTINE computeEtan

  SUBROUTINE computeVun(U,Vun)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: Vun(:)
    real*8             :: U5,U6
 
    Vun = 0.
    
    U5 = U(5)
    U6 = U(6)

    Vun(5) = - U6/U5**2
    Vun(6) = 1./U5

  END SUBROUTINE computeVun

  SUBROUTINE compute_dVun_dU(U,dVun_dU)
    real*8, intent(IN)  :: U(:)
    real*8, intent(OUT) :: dVun_dU(:, :)
    real*8              :: U5,U6

    dVun_dU = 0.

    U5 = U(5)
    U6 = U(6)

    dVun_dU(5, 5) = 2.*U6/U5**3
    dVun_dU(5, 6) = - 1./U5**2
    
    dVun_dU(6, 5) = - 1./U5**2
  END SUBROUTINE compute_dVun_dU

#endif
! NEUTRAL GAMMA

#endif
! TEMPERATURE

#endif
! NEUTRAL

  !***********************************************************************
  !
  !    COMPUTATION OF THE STABILIZATION PARAMETER
  !
  !***********************************************************************

  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up, uc, q, b, n, iel, ifa, isext, xy, tau)
    real*8, intent(in)  :: up(:), uc(:), q(:), b(:), n(:), xy(:)
    real, intent(in)    :: isext
    integer, intent(in) :: ifa, iel
    real*8, intent(out) :: tau(:, :)
    real*8              :: tau_aux(phys%Neq),diff_iso(phys%Neq,phys%Neq,1),diff_ani(phys%Neq,phys%Neq,1)
#ifdef NEUTRALP
    real*8              :: Dpn
    real*8              :: Vpn(simpar%Neq),Qpr(simpar%Ndim,simpar%Neq)
#endif
#ifdef NEUTRALGAMMA
    real*8              :: Etan
#endif
    integer             :: ndim
    real*8              :: xc, yc, rad, h, aux, bn, bnorm,xyd(1,size(xy)),uu(1,size(uc)),qq(1,size(q))
    real*8              :: U1, U2, U3, U4
    
    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    tau = 0.
    ndim = size(n)
    bn = dot_product(b(1:ndim), n)
    bnorm = norm2(b(1:ndim))
    xyd(1,:) = xy(:)
    uu(1,:) = uc(:)
    qq(1,:) = q(:)

    call setLocalDiff(xyd, uu, qq, diff_iso, diff_ani)
    
#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(uc,Vpn)
    ! Compute Dpn(U^(k-1))
    Qpr = reshape(q,(/simpar%Ndim,simpar%Neq/))
    !CALL computeDpn(uc,Qpr,Vpn,Dpn)
#endif

#ifdef NEUTRALGAMMA
    CALL computeEtan(uc,Etan)
#endif

    if (numer%stab == 2) then
      if (abs(isext - 1.) .lt. 1e-12) then
        ! exterior faces
        tau_aux = abs((4*uc(2)*bn)/uc(1))
      else
        tau_aux = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(1) + sqrt(abs(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)))
      endif
#ifdef TOR3D
      if (abs(n(3)) > 0.1) then
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifdef NEUTRAL
        tau_aux(5) = tau_aux(5) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
      else
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
        tau_aux(5) = tau_aux(5) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif
#ifdef TOR3D
      endif
#endif

    elseif (numer%stab == 3) then
      if (abs(isext - 1.) .lt. 1e-12) then
        ! exterior faces
        tau_aux = max(abs((5*U2 - 2*U2*phys%Gmbohme)/(3*U1)), abs((5*U2 - 2*U2*phys%Gmbohm)/(3*U1)))
      else
        tau_aux = max(abs((3*U2*bn + 5**(0.5)*bn*(-U2**2 + 2*U1*U3 + 2*U1*U4)**(0.5))/(3*U1)),&
          &abs((3*U2*bn - 5**(0.5)*bn*(-U2**2 + 2*U1*U3 + 2*U1*U4)**(0.5))/(3*U1)),&
          &abs((U2*bn)/U1),&
          &abs((5*U2*bn)/(3*U1)))
      endif
#ifdef TOR3D
      if (abs(n(3)) > 0.1) then
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifdef NEUTRAL
        tau_aux(5) = tau_aux(5) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
      else
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
        tau_aux(5) = tau_aux(5) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif
#ifdef TOR3D
      endif
#endif

    elseif (numer%stab == 4) then
      tau_aux = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(1) + sqrt(abs(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)), &
        phys%lscale/geom%R0*abs(bn)*phys%diff_pari*up(7)**2.5, phys%lscale/geom%R0*abs(bn)*phys%diff_pare*up(8)**2.5)

    elseif (numer%stab == 5) then
      if (abs(isext - 1.) .lt. 1e-12) then
        ! exterior faces
        tau_aux = abs((4*uc(2)*bn)/uc(1))
      else
        tau_aux = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(1) + sqrt(abs(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)))
!#ifdef NEUTRALP 
!        tau_aux(5) = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(5) + sqrt(abs(10*uc(3)/uc(1)*uc(5)**2 - 5*(uc(5)*uc(2)/uc(1))**2)))/uc(5)))
!#endif      
      endif
#ifdef TOR3D
      if (abs(n(3)) > 0.1) then
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n
        tau_aux(2) = tau_aux(2) + phys%diff_u
        tau_aux(3) = tau_aux(3) + phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(4) = tau_aux(4) + phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifndef NEUTRALP
#ifdef NEUTRAL
        tau_aux(5) = numer%tau(5) !tau_aux(5) + diff_iso(5,5,1)
#endif
#else
        tau_aux(5) = tau_aux(5) + phys%diff_nn 
#endif
      else
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + 6.*diff_iso(1,1,1)
        tau_aux(2) = tau_aux(2) + 6.*diff_iso(2,2,1)
        tau_aux(3) = tau_aux(3) + 6.*diff_iso(3,3,1) + abs(bn)*phys%diff_pari*(min(1.,up(7)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(4) = tau_aux(4) + 6.*diff_iso(4,4,1) + abs(bn)*phys%diff_pare*(min(1.,up(8)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifndef NEUTRALP
#ifdef NEUTRAL
        tau_aux(5) = tau_aux(5) + diff_iso(5,5,1) !phys%diff_nn !numer%tau(5) 
#endif
#else
        tau_aux(5) = tau_aux(5) + numer%tau(5) !Dpn
#endif
#ifdef NEUTRALGAMMA
        tau_aux(6) = tau_aux(6) + Etan/uc(5)         ! diff_iso(5,5,1)
#endif
!        ! Toroidal face
!        tau_aux(1) = tau_aux(1) +  diff_iso(1,1,1)
!        tau_aux(2) = tau_aux(2) +  diff_iso(2,2,1)
!        tau_aux(3) = tau_aux(3) +  diff_iso(3,3,1) + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!        tau_aux(4) = tau_aux(4) +  diff_iso(4,4,1) + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!#ifdef NEUTRAL
!        tau_aux(5) = tau_aux(5) +  diff_iso(5,5,1)
!#endif
#ifdef TOR3D
      endif
#endif

    else
      write (6, *) "Wrong stabilization type: ", numer%stab
      stop
    endif
    tau(1, 1) = tau_aux(1)
    tau(2, 2) = tau_aux(2)
    tau(3, 3) = tau_aux(3)
    tau(4, 4) = tau_aux(4)
#ifdef NEUTRAL
    tau(5, 5) = tau_aux(5)
#endif
#ifdef NEUTRALGAMMA
    tau(6,6) = tau_aux(6)
#endif
  END SUBROUTINE computeTauGaussPoints

  !!

  SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)

    real*8, intent(in)  :: up(:), uc(:), b(:), n(:), xy(:), isext
    real*8, intent(out) :: tau(:, :)
    integer, intent(in) :: iel
    real*8, parameter :: eps = 1e-12
    real*8              :: bn, bnorm

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

    tau = 0.
    bn = dot_product(b, n)
    bnorm = norm2(b)

    x = xy(1)
    y = xy(2)

    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    !************************************
    !
    ! *****     CONVECTIVE PART  ********
    !
    !************************************
    IF (abs(isext - 1.) .lt. 1e-12) THEN

      !************************************
      !   EXTERIOR FACES
      !************************************
      tau(3, 1) = (1.0D0/U1**2*abs(U2*bn)*(U1*U3 - U2**2)*(-4.0D0))/abs(U1)
      tau(3, 2) = (abs(U2*bn)*(U1*U3*2.0D0 - U2**2*3.0D0)*2.0D0)/(U1*(U2 + eps)*abs(U1))
      tau(3, 3) = (abs(U2*bn)*4.0D0)/abs(U1)

    ELSE
      !************************************
      !   INTERIOR FACES
      !************************************
      t2 = abs(U1)
      t3 = U1*U3*2.0D0
      t4 = U1*U4*2.0D0
      t5 = U2**2
      !          t6 = t3+t4-t5
      t6 = abs(t3 + t4 - t5)
      t7 = U2*bn
      t8 = abs(t7)
      t9 = t5**2
      t10 = U2*bn*3.0D0
      t11 = sqrt(5.0D0)
      t12 = sqrt(t6)
      t13 = bn*t11*t12
      t14 = t10 + t13
      t15 = abs(t14)
      t16 = t10 - t13
      t17 = abs(t16)
      t18 = t6**(3.0D0/2.0D0)
      t19 = 1.0D0/t6
      t20 = 1.0D0/t2
      t21 = U2*3.0D0
      t22 = t11*t12
      t23 = U1**2
      t24 = t21 + t22
      t25 = bn*t24
      t26 = abs(t25)
      t27 = t21 - t22
      t28 = bn*t27
      t29 = abs(t28)
      t30 = t8*(-6.0D0) + t26 + t29
      t31 = t19*t20*t23*t30*(1.0D0/5.0D0)
      t32 = 1.0D0/U1
      t33 = U1*U3*1.0D1
      t34 = U1*U4*1.0D1
      t51 = t5*9.0D0
      t35 = t33 + t34 - t51
      t36 = 1.0D0/t35
      t37 = U2*t12*5.0D0
      t38 = U1*U3*t11*1.0D1
      t39 = U1*U4*t11*1.0D1
      t40 = U2*t11*t12
      t41 = t3 + t4
      t42 = t2*t6*1.35D2
      t43 = t42 - t2*t41*6.0D1
      t44 = 1.0D0/t43
      t45 = U1*U2*t5*t8*7.2D1
      t46 = U1*U2*t6*t15*1.5D1
      t47 = U1*U2*t6*t17*1.5D1
      t48 = U1*t11*t15*t18*5.0D0
      t49 = U1*t5*t11*t12*t17*4.0D0
      t50 = t19*t44*(t45+t46+t47+t48+t49-U1*t11*t17*t18*5.0D0-U1*U2*t6*t8*9.0D1-U1*U2*t5*t15*1.2D1-U1*U2*t5*t17*1.2D1-U1*t5*t11*t12*t15*4.0D0)
      t52 = U1*U3*5.0D0
      t53 = U1*U4*5.0D0
      t54 = t5*(-4.0D0) + t52 + t53
      t55 = 1.0D0/U1**2
      t56 = t5*7.0D0
      t57 = 1.0D0/t6**(3.0D0/2.0D0)
      t58 = t33 + t34 - t40 - t56
      t59 = t9*1.5D1
      t60 = U3**2
      t61 = t23*t60*5.0D1
      t62 = U3*U4*t23*5.0D1
      t63 = U2*t5*t11*t12*3.0D0
      t65 = U1*U3*t5*5.5D1
      t66 = U1*U4*t5*3.0D1
      t64 = t59 + t61 + t62 + t63 - t65 - t66
      t67 = U2*2.0D0
      t68 = U1*U3*5.0D1
      t69 = U1*U4*5.0D1
      t73 = t5*4.5D1
      t70 = t68 + t69 - t73
      t71 = 1.0D0/t70
      t72 = t22 + t67
      t74 = t59 + t61 + t62 - t63 - t65 - t66
      t75 = t11*t20*t26*t57*t71*t72*t74*(1.0D0/1.5D1)
      t76 = t11*t20*t29*t57*t64*t71*(t22 - t67)*(1.0D0/1.5D1)
      t77 = 1.0D0/sqrt(t6)
      t78 = U1*U4*t12*t26*5.0D0
      t79 = U1*U4*t12*t29*5.0D0
      t80 = U1*U2*U4*t11*t26*2.0D0
      tau(1,1) = (t19*(t8*t9*2.4D1-t9*t15*4.0D0-t9*t17*4.0D0+t6**2*t8*5.0D1-t5*t6*t8*7.0D1+t5*t6*t15*5.0D0+t5*t6*t17*5.0D0-U2*t11*t15*t18*5.0D0+U2*t11*t17*t18*5.0D0+U2*t5*t11*t12*t15*4.0D0-U2*t5*t11*t12*t17*4.0D0))/(t2*t6*9.0D1-t2*t41*4.0D1)
      tau(1, 2) = t19*t20*(U1*U2*t8*(-1.2D1) + U1*U2*t15*2.0D0 + U1*U2*t17*2.0D0 - U1*t11*t12*t15 + U1*t11*t12*t17)*(-1.0D0/1.0D1)
      tau(1, 3) = t31
      tau(1, 4) = t31
      tau(2,1) = U2*t8*t19*t20*t32*t54*(2.0D0/5.0D0)+U2*t11*t19*t20*t29*t32*t36*t58*(t37-t38-t39+t5*t11*1.1D1)*(1.0D0/1.5D2)-U2*t11*t19*t20*t26*t32*t36*(t37+t38+t39-t5*t11*1.1D1)*(t5*(-7.0D0)+t33+t34+t40)*(1.0D0/1.5D2)
      tau(2,2) = t19*t20*(t5*t8*3.6D1-t5*t15*6.0D0+t6*t15*5.0D0-t5*t17*6.0D0+t6*t17*5.0D0+U2*t11*t12*t15-U2*t11*t12*t17)*(1.0D0/3.0D1)
      tau(2, 3) = t50
      tau(2, 4) = t50
      tau(3,1) = t5*t8*t19*t20*t54*t55*(1.0D0/5.0D0)-U4*t5*t8*t20*t32*t36*(2.5D1/3.0D0)+U2*t11*t20*t29*t36*t55*t57*t58*t64*(1.0D0/1.5D2)-U2*t11*t20*t26*t36*t55*t57*(t33+t34+t40-t56)*(t59+t61+t62-U1*U3*t5*5.5D1-U1*U4*t5*3.0D1-U2*t5*t11*t12*3.0D0)*(1.0D0/1.5D2)
      tau(3,2) = t11*t20*t29*t32*t57*t64*(-1.0D0/1.5D2)+t11*t20*t26*t32*t57*(t59+t61+t62-t63-U1*U3*t5*5.5D1-U1*U4*t5*3.0D1)*(1.0D0/1.5D2)+U2*t5*t8*t19*t20*t32*(3.0D0/5.0D0)
      tau(3, 3) = t75 + t76 - t5*t8*t19*t20*(3.0D0/5.0D0) + U1*U4*t8*t20*t36*(5.0D1/3.0D0)
      tau(3, 4) = t75 + t76 - t5*t8*t19*t20*(3.0D0/5.0D0) - t8*t20*t36*(t33 - t51)*(5.0D0/3.0D0)
      tau(4,1) = (t11*t77*(U2*U4*t5*t15*(-1.0D1)+U2*U4*t6*t15*2.5D1+U2*U4*t5*t17*1.0D1-U2*U4*t6*t17*2.5D1-U4*t5*t8*t11*t12*5.0D1+U4*t5*t11*t12*t15*5.0D0+U4*t5*t11*t12*t17*5.0D0)*(-1.0D0/5.0D0))/(U1*t2*t6*5.4D1-U1*t2*t41*2.4D1)
      tau(4, 2) = t20*t77*(U4*t11*t15 - U4*t11*t17)*(1.0D0/6.0D0)
      tau(4, 3) = t20*t36*t77*(t78 + t79 + t80 - U1*U4*t8*t12*5.0D1 - U1*U2*U4*t11*t29*2.0D0)*(1.0D0/3.0D0)
      tau(4, 4) = t20*t36*t77*(t78 + t79 + t80 - t5*t8*t12*4.5D1 + U1*U3*t8*t12*5.0D1 - U1*U2*U4*t11*t29*2.0D0)*(1.0D0/3.0D0)
    END IF

    !************************************
    !
    ! *****     DIFFUSIVE PART  ********
    !
    !************************************
    tau(1, 1) = tau(1, 1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(2, 2) = tau(2, 2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(3, 3) = tau(3, 3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(4, 4) = tau(4, 4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale

  END SUBROUTINE computeTauGaussPoints_matrix

!  !***********************************************************************
!  !
!  !                           MAGNETIC FIELD
!  !
!  !***********************************************************************
!  #ifdef TOR3D
!  SUBROUTINE defineMagneticField(x,y,t,b,divb,drift)
!    real*8, intent(in)      :: x(:),y(:),t(:)
!    real*8, intent(out),optional     :::: b(:,:),divb(:),drift(:,:)
!    real*8                  :: R0,q,r
!    real*8                  :: xmax,xmin,ymax,ymin,xm,ym,p,divbp
!    real*8                  :: xx,yy,tt,Bp,Bt,Br,Bz,BB,dmax,B0,xr,yr
!    integer*4               :: i,j,ind,N2D,N1D
!
!
!    N2d = size(X,1)
!    N1d = size(t,1)
!    xmax = Mesh%xmax
!    xmin = Mesh%xmin
!    ymax = Mesh%ymax
!    ymin = Mesh%ymin
!    xm = 0.5*(xmax+xmin)
!    ym = 0.5*(ymax+ymin)
!    !  xm = -0.5
!    !  ym = -0.5
!    ! Initialization
!    ! Initialization
!    if (present(b)) then
!      b     = 0.
!    endif
!    if (present(divb)) then
!      divb  = 0.
!    endif
!    if (present(drift)) then
!      drift = 0.
!    endif
!    if (present(Bmod)) then
!      Bmod  = 0.
!    endif
!
!    DO i=1,N2d
!      DO j=1,N1d
!        xx = x(i)
!        yy = y(i)
!        tt = t(j)
!        ind = (j-1)*N2d+i
!
!        SELECT CASE(switch%testcase)
!        CASE(1)
!          IF (switch%axisym) THEN
!            WRITE(6,*) "This is NOT an axisymmetric test case!"
!            stop
!          END IF
!          ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!          Br = (yy-ym)
!          Bz = (-xx+xm)
!          Bt = 1.
!          divbp = 0.
!        CASE(2)
!          IF (.not.switch%axisym) THEN
!            WRITE(6,*) "This is an axisymmetric test case!"
!            stop
!          END IF
!          ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!          Br = (yy-ym)/xx
!          Bz = (-xx+xm)/xx
!          Bt = 1.
!          divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!        CASE(3)
!          IF (.not.switch%axisym) THEN
!            WRITE(6,*) "This is an axisymmetric test case!"
!            stop
!          END IF
!          ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!          Br = (yy-ym)/xx
!          Bz = (-xx+xm)/xx
!          Bt = 1.
!          divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!
!        CASE(50:59)
!          write(6,*) "Error in defineMagneticField: you should not be here!"
!          STOP
!        CASE(60:69)
!
!          ! Circular case with limiter
!          R0 = geom%R0
!          q  = geom%q
!          B0 = 2 ! not influential
!          xr = xx*phys%lscale
!          yr = yy*phys%lscale
!
!          r  = sqrt((xr-R0)**2+yr**2)
!          Br = -B0*yr/(xr*q*sqrt(1- (r/R0)**2 ) )
!          Bz = B0*(xr-R0)/(xr*q*sqrt(1- (r/R0)**2 ) )
!          Bt = B0*R0/xr
!
!          if (present(divb)) then
!            IF (switch%axisym) THEN
!              divbp = -yy/xx/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
!            ELSE
!              WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
!              STOP
!            END IF
!          endif
!          if (present(divb)) then
!            IF (switch%driftdia) THEN
!              drift(:,2) =  -1./R0*phys%lscale
!            END IF
!          endif
!
!
!
!
!
!
!        CASE DEFAULT
!          WRITE(6,*) "Error! Test case not valid"
!          STOP
!        END SELECT
!
!        Bp = sqrt(Br**2+Bz**2)
!        BB = sqrt(Bp**2+Bt**2)
!        b(ind,1) = Br/BB
!        b(ind,2) = Bz/BB
!        b(ind,3) = Bt/BB
!        if (present(divb)) then
!          divb(ind) = divbp
!        endif
!      END DO
!    END DO
!  END SUBROUTINE defineMagneticField
!  #else
!  SUBROUTINE defineMagneticField(x,y,b,divb,drift,Bmod)
!    real*8, intent(in)      :: x(:),y(:)
!    real*8, intent(out)     :: b(:,:)
!    real*8, intent(out),optional::divb(:),drift(:,:),Bmod(:)
!    real*8                  :: R0,q,r(size(x)),auxdivb(size(x)),auxdrift(size(b,1),size(b,2))
!    real*8                  :: xmax,xmin,ymax,ymin,xm,ym,xx,yy
!    real*8                  :: Br,Bz,Bt,BB,Bp
!    integer :: i
!    ! Initialization
!    b     = 0.
!    auxdivb(:) = 0.
!    xmax = Mesh%xmax
!    xmin = Mesh%xmin
!    ymax = Mesh%ymax
!    ymin = Mesh%ymin
!    xm = 0.5*(xmax+xmin)
!    ym = 0.5*(ymax+ymin)
!    if (present(Bmod)) then
!      Bmod  = 0.
!    endif
!    !  xm = -0.5
!    !  ym = -0.5
!    SELECT CASE(switch%testcase)
!    CASE(1)
!      ! Cartesian case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
!      !            b(:,1) = 1./30.*(x-y**2+2)
!      !            b(:,2) = 1./30.*(x*y+y)
!      !            auxdivb(:) = 1./15.+(1./30.)*x
!      DO i=1,size(x)
!        xx = x(i)
!        yy = y(i)
!        Br = (yy-ym)
!        Bz = (-xx+xm)
!        Bt = 1.
!        Bp = sqrt(Br**2+Bz**2)
!        BB = sqrt(Bp**2+Bt**2)
!        b(i,1) = Br/BB
!        b(i,2) = Bz/BB
!        auxdivb(i) = 0.
!      END DO
!
!    CASE(2)
!      ! Axisymmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
!      !            b(:,1) = 1./30.*(x-y**2+2)
!      !            b(:,2) = 1./30.*(x*y+y)
!      !            auxdivb(:) = 1./30.+((1./30.)*x-(1./30.)*y**2+1./15.)/x+1./30.*(x+1)
!      DO i=1,size(x)
!        xx = x(i)
!        yy = y(i)
!        Br = (yy-ym)/xx
!        Bz = (-xx+xm)/xx
!        Bt = 1.
!        Bp = sqrt(Br**2+Bz**2)
!        BB = sqrt(Bp**2+Bt**2)
!        b(i,1) = Br/BB
!        b(i,2) = Bz/BB
!        auxdivb(i) = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!      END DO
!    CASE(5)
!      ! Cartesian case, square mesh, horizontal field
!      b(:,1) = 0.1
!      b(:,2) = 0.
!    CASE(6)
!      ! Cartesian case, square mesh, horizontal field
!      DO i=1,size(x)
!        IF (y(i).ge.0.5) THEN
!          b(i,1) = 0.1
!        ELSE
!          b(i,1) = -0.1
!        END IF
!      END DO
!    CASE(50:59)
!      write(6,*) "Error in defineMagneticField: you should not be here!"
!      STOP
!    CASE(60:69)
!
!      ! Circular case with limiter
!      R0 = geom%R0
!      q  = geom%q
!      r  = phys%lscale*sqrt((x-R0/phys%lscale)**2+y**2)
!      b(:,1) = -phys%lscale*y/sqrt(R0**2*q**2+(1-q**2)*r**2)
!      b(:,2) = phys%lscale*(x-R0/phys%lscale)/sqrt(R0**2*q**2+(1-q**2)*r**2)
!
!      IF (switch%axisym) THEN
!        auxdivb(:) = -y/x/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
!      ELSE
!        WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
!        STOP
!      END IF
!
!      IF (switch%driftdia) THEN
!        auxdrift(:,2) =  -1./R0*phys%lscale
!      END IF
!    CASE DEFAULT
!      WRITE(6,*) "Error! Test case not valid"
!      STOP
!    END SELECT
!
!    IF (present(divb)) THEN
!      divb = auxdivb
!    ENDIF
!    IF (present(drift)) THEN
!      drift = auxdrift
!    ENDIF
!  END SUBROUTINE defineMagneticField
!  #endif
!
!  SUBROUTINE loadMagneticField()
!    USE interpolation
!    USE HDF5
!    USE HDF5_io_module
!    integer        :: i,ierr,ip,jp
!    integer(HID_T) :: file_id
!    real*8,pointer,dimension(:,:) :: r2D,z2D,flux2D,Br2D,Bz2D,Bphi2D
!    real*8,allocatable,dimension(:,:) :: bx,by,bmod,divb,bmodx,bmody,driftx,drifty
!    real*8,allocatable,dimension(:)   :: xvec,yvec
!    real*8                            :: x,y
!
!    WRITE(6,*) "******* Loading magnetic field *******"
!    ! Allocate storing space in phys
!    ALLOCATE(phys%b(Mesh%Nnodes,Mesh%Ndim))
!    ALLOCATE(phys%divb(Mesh%Nnodes))
!    ALLOCATE(phys%drift(Mesh%Nnodes,Mesh%Ndim))
!    ALLOCATE(phys%flux2d(Mesh%Nnodes))
!
!    ! Dimensions of the file storing the magnetic field for West
!    ip = 541
!    jp = 391
!    ALLOCATE(r2D(ip,jp))
!    ALLOCATE(z2D(ip,jp))
!    ALLOCATE(flux2D(ip,jp))
!    ALLOCATE(Br2D(ip,jp))
!    ALLOCATE(Bz2D(ip,jp))
!    ALLOCATE(Bphi2D(ip,jp))
!    ALLOCATE(bx(ip,jp))
!    ALLOCATE(by(ip,jp))
!    ALLOCATE(bmod(ip,jp))
!    ALLOCATE(bmodx(ip,jp))
!    ALLOCATE(bmody(ip,jp))
!    ALLOCATE(divb(ip,jp))
!    ALLOCATE(driftx(ip,jp))
!    ALLOCATE(drifty(ip,jp))
!
!    ! Read file
!    CALL HDF5_open('WEST_far_465.h5',file_id,IERR)
!    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
!    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
!    CALL HDF5_array2D_reading(file_id,flux2D,'flux2D')
!    CALL HDF5_array2D_reading(file_id,Br2D,'Br2D')
!    CALL HDF5_array2D_reading(file_id,Bz2D,'Bz2D')
!    CALL HDF5_array2D_reading(file_id,Bphi2D,'Bphi2D')
!    CALL HDF5_close(file_id)
!
!    ! Apply length scale
!    r2D = r2D/phys%lscale
!    z2D = z2D/phys%lscale
!
!    ! Compute b
!    bmod = sqrt(Br2D**2+Bz2D**2+Bphi2D**2)
!    bx = -Br2D/bmod
!    by = -Bz2D/bmod
!
!    ! Compute divergence of b
!    divb = 0.
!    IF (switch%axisym) THEN
!      ! 1/r*(d r*br/dr)+dbz/dz
!      divb(2:ip-1,2:jp-1) = 1./r2D(2:ip-1,2:jp-1)*(r2D(2:ip-1,3:jp)*bx(2:ip-1,3:jp)- &
!        r2D(2:ip-1,1:jp-2)*bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-  &
!        r2D(2:ip-1,1:jp-2))+(by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!
!    ELSE
!      ! dbr/dr+dbz/dz
!      divb(2:ip-1,2:jp-1) = (bx(2:ip-1,3:jp-1)-bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))+ &
!        (by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!    END IF
!
!    ! Compute drift velocity
!    driftx = 0.
!    drifty = 0.
!    IF (switch%driftdia) THEN
!      bmodx = (bmod(2:ip-1,3:jp)-bmod(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))
!      bmody = (bmod(3:ip,2:jp-1)-bmod(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
!      driftx(2:ip-1,2:jp-1) =  -Bphi2D(2:ip-1,2:jp-1)*bmody/bmod(2:ip-1,2:jp-1)**3
!      drifty(2:ip-1,2:jp-1) =   Bphi2D(2:ip-1,2:jp-1)*bmodx/bmod(2:ip-1,2:jp-1)**3
!    END IF
!
!    ! Interpolate
!    ALLOCATE(xvec(jp))
!    ALLOCATE(yvec(ip))
!    xvec = r2D(1,:)
!    yvec = z2D(:,1)
!    DO i = 1,Mesh%Nnodes
!      x = Mesh%X(i,1)
!      y = Mesh%X(i,2)
!      phys%b(i,1) = interpolate( ip, yvec,jp, xvec, bx, y,x, 1e-12)
!      phys%b(i,2) = interpolate( ip, yvec,jp, xvec, by, y,x, 1e-12)
!      phys%divb(i) = interpolate( ip, yvec,jp, xvec, divb, y,x, 1e-12)
!      phys%drift(i,1) = interpolate( ip, yvec,jp, xvec,driftx, y,x, 1e-12)
!      phys%drift(i,2) = interpolate( ip, yvec,jp, xvec,drifty, y,x, 1e-12)
!      phys%flux2D(i) = interpolate( ip, yvec,jp, xvec,flux2D, y,x, 1e-12)
!    END DO
!
!    ! Free memory
!    DEALLOCATE(Br2D,Bz2D,Bphi2D,xvec,yvec)
!    DEALLOCATE(r2D,z2D,flux2D,bx,by,bmod,bmodx,bmody,divb,driftx,drifty)
!
!  END SUBROUTINE loadMagneticField
!
!
!
!  SUBROUTINE loadMagneticFieldTemporalEvolution()
!    USE HDF5
!    USE HDF5_io_module
!    USE MPI_OMP
!    integer        :: ierr,k
!    character(LEN=20) :: fname = 'Evolving_equilibrium'
!    character(10)  :: npr,nid,nit
!    character(len=1000) :: fname_complete
!    integer(HID_T) :: file_id
!
!    WRITE(6,*) "******* Loading magnetic field *******"
!
!    ! Allocate storing space in phys
!    ALLOCATE(phys%Br(Mesh%Nnodes))
!    ALLOCATE(phys%Bz(Mesh%Nnodes))
!    ALLOCATE(phys%Bt(Mesh%Nnodes))
!    ALLOCATE(phys%flux2d(Mesh%Nnodes))
!
!    ! File name
!    write(nit, "(i10)") time%it
!    nit = trim(adjustl(nit))
!    k = INDEX(nit, " ") -1
!
!    IF (MPIvar%glob_size.GT.1) THEN
!      write(nid,*) MPIvar%glob_id+1
!      write(npr,*) MPIvar%glob_size
!      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
!    ELSE
!      fname_complete = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
!    END IF
!
!    write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))
!
!    ! Read file
!    CALL HDF5_open(fname_complete,file_id,IERR)
!    CALL HDF5_array1D_reading(file_id,phys%Br,'Br')
!    CALL HDF5_array1D_reading(file_id,phys%Bz,'Bz')
!    CALL HDF5_array1D_reading(file_id,phys%Bt,'Bt')
!    CALL HDF5_array1D_reading(file_id,phys%flux2D,'flux')
!    CALL HDF5_close(file_id)
!
!  END SUBROUTINE loadMagneticFieldTemporalEvolution


END MODULE physics
