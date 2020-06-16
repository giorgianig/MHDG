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

      ! number of physical variables
      phys%npv = 10

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
      phys%phyVarNam(10) = "M"   ! Mach

      ! Set the name of the conservative variables
      phys%conVarNam(1) = "rho"   ! U1 = rho
      phys%conVarNam(2) = "Gamma" ! U2 = rho*u
      phys%conVarNam(3) = "nEi"   ! U3 = rho*Ei
      phys%conVarNam(4) = "nEe"   ! U4 = rho*Ee

      simpar%model = 'N-Gamma-Ti-Te'
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
      simpar%consvar_refval(1) = simpar%refval_density
      simpar%consvar_refval(2) = simpar%refval_momentum
      simpar%consvar_refval(3) = simpar%refval_specenergydens
      simpar%consvar_refval(4) = simpar%refval_specenergydens
   END SUBROUTINE

   !*******************************************
   ! Convert physical variable to conservative
   ! variables
   !*******************************************
   SUBROUTINE phys2cons(up, ua)
      real*8, dimension(:, :), intent(in)  :: up
      real*8, dimension(:, :), intent(out) :: ua

      ua(:, 1) = up(:, 1)
      ua(:, 2) = up(:, 1)*up(:, 2)
      ua(:, 3) = up(:, 1)*up(:, 3)
      ua(:, 4) = up(:, 1)*up(:, 4)

   END SUBROUTINE phys2cons

   !*******************************************
   ! Convert conservative variable to physical
   ! variables
   !*******************************************
   SUBROUTINE cons2phys(ua, up)
      real*8, dimension(:, :), intent(in)  :: ua
      real*8, dimension(:, :), intent(out) :: up
      integer :: i

      up(:, 1) = abs(ua(:, 1))                                                 ! density
      up(:, 2) = ua(:, 2)/ua(:, 1)                                         ! u parallel
      up(:, 3) = ua(:, 3)/ua(:, 1)                                         ! total energy of ions
      up(:, 4) = ua(:, 4)/ua(:, 1)                                         ! total energy of electrons
      up(:, 5) = abs(2./(3.*phys%Mref)*(ua(:, 3) - 0.5*ua(:, 2)**2/ua(:, 1))) ! pressure of ions
      up(:, 6) = abs(2./(3.*phys%Mref)*ua(:, 4))                          ! pressure of electrons
      up(:, 7) = up(:, 5)/up(:, 1)                                         ! temperature of ions
      up(:, 8) = up(:, 6)/up(:, 1)                                         ! temperature of electrons
      up(:, 9) = sqrt((up(:, 7) + up(:, 8))*phys%Mref)                     ! sound speed
      up(:, 10) = up(:, 2)/up(:, 9)                                         ! Mach

   END SUBROUTINE cons2phys

   !*****************************************
   ! Jacobian matrices
   !****************************************
   SUBROUTINE jacobianMatrices(U, A)
      real*8, intent(in)  :: U(:)
      real*8, intent(out) :: A(:, :)
!                [       0,                                                1,                                    0,              0; ...
!     -2/3*U(2)**2/U(1)**2                              4/3*U(2)/U(1)                               2/3,            2/3; ...
!     -5/3*U(2)*U(3)/U(1)**2+2/3*U(2)**3/U(1)**3      5/3*U(3)/U(1)-U(2)**2/U(1)**2                5/3*U(2)/U(1),       0 ;   ...
!     -5/3*U(4)*U(2)/U(1)**2,                           5/3*U(4)/U(1),                             0,           5/3*U(2)/U(1)]

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
      end if

   END SUBROUTINE jacobianMatrices

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
      endif
      An = bn*An
   END SUBROUTINE jacobianMatricesFace

   !*****************************************
   ! Jacobian matrix for the Bohm BC
   !****************************************
   SUBROUTINE jacobianMatricesBohm(U, A)
      real*8, intent(in)  :: U(:)
      real*8, intent(out) :: A(:, :)
      real*8              :: auxi, auxe

      A = 0.
      auxi = (5.-2.*phys%Gmbohm)/3.
      auxe = (5.-2.*phys%Gmbohme)/3.
      A(1, 1) = -auxi*(U(2)**3/U(1)**3 - U(2)*U(3)/U(1)**2)
      A(1, 2) = auxi*(U(3)/U(1) - 3./2.*U(2)**2/U(1)**2)
      A(1, 3) = auxi*U(2)/U(1)

      A(2, 1) = -auxe*U(2)*U(4)/U(1)**2
      A(2, 2) = auxe*U(4)/U(1)
      A(2, 4) = auxe*U(2)/U(1)
   END SUBROUTINE jacobianMatricesBohm

   !*****************************************
   ! Set the perpendicular diffusion
   !****************************************
   SUBROUTINE setLocalDiff(xy, d_iso, d_ani, Bmod)
      real*8, intent(in)  :: xy(:, :), Bmod(:)
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
      d_iso(3, 3, :) = phys%diff_e
      d_iso(4, 4, :) = phys%diff_ee

      d_ani(1, 1, :) = phys%diff_n
      d_ani(2, 2, :) = phys%diff_u
      d_ani(3, 3, :) = phys%diff_e
      d_ani(4, 4, :) = phys%diff_ee

      !*****************************
      ! Non diagonal terms
      !*****************************
      ! No non-diagonal terms defined for this model

      if (switch%difcor .gt. 0) then
         call computeIperDiffusion(xy, iperdiff)
         d_iso(1, 1, :) = d_iso(1, 1, :)*iperdiff
         d_iso(2, 2, :) = d_iso(2, 2, :)*iperdiff
         d_iso(3, 3, :) = d_iso(3, 3, :)*iperdiff
         d_iso(4, 4, :) = d_iso(4, 4, :)*iperdiff
      endif

   END SUBROUTINE setLocalDiff

   !*******************************************
   ! Compute local diffusion in points
   !*******************************************
   SUBROUTINE computeIperDiffusion(X, ipdiff)
      real*8, intent(IN)  :: X(:, :)
      real*8, intent(OUT) :: ipdiff(:)
      real*8             :: xcorn, ycorn
      real*8             :: rad(size(X, 1))
      real*8             :: h
      real*8, parameter   :: tol = 1.e-6
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
         WRITE (6, *) "Case not valid"
         STOP
      END SELECT

                                                                !!**********************************************************
                                                                !! Gaussian around the corner
                                                                !!**********************************************************
      h = 10e-3
      rad = sqrt((X(:, 1)*phys%lscale - xcorn)**2 + (X(:, 2)*phys%lscale - ycorn)**2)
      ipdiff = 1 + numer%dc_coe*exp(-(2*rad/h)**2)

   END SUBROUTINE computeIperDiffusion

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
      if (aux < tol) aux = tol
      res = aux**phys%epn
   END FUNCTION computeAlphai

   FUNCTION computeAlphae(U) RESULT(res)
      real*8 :: U(:)
      real*8 :: res, aux
      real, parameter :: tol = 1e-5
      aux = U(4)/U(1)
      if (aux < tol) aux = tol
      res = aux**phys%epn
   END FUNCTION computeAlphae

   SUBROUTINE compute_dAlpha_dUi(U, res)
      real*8, intent(IN) :: U(:)
      real*8, intent(OUT):: res(:)
      real*8             :: aux
      real, parameter :: tol = 1e-5
      aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
      if (aux < 0) aux = tol
      res = 0.d0
      res(1) = -U(3)/U(1)**2 + U(2)**2/U(1)**3
      res(2) = -U(2)/U(1)**2
      res(3) = 1./U(1)
      res = phys%epn*aux**(phys%epn - 1)*res
   END SUBROUTINE compute_dAlpha_dUi

   SUBROUTINE compute_dAlpha_dUe(U, res)
      real*8, intent(IN) :: U(:)
      real*8, intent(OUT):: res(:)
      real*8             :: aux
      real, parameter :: tol = 1e-5
      aux = U(4)/U(1)
      if (aux < 0) aux = tol
      res = 0.d0
      res(1) = -U(4)/U(1)**2
      res(4) = 1./U(1)
      res = phys%epn*aux**(phys%epn - 1)*res
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
      s = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4 - U3 + 0.5*(U(2)**2/U1))
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
      res = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*res
   END SUBROUTINE compute_ds_dU

   !***********************************************************************
   !
   !    COMPUTATION OF THE STABILIZATION PARAMETER
   !
   !***********************************************************************

   !*******************************************
   ! Compute the stabilization tensor tau
   !*******************************************
   SUBROUTINE computeTauGaussPoints(up, uc, b, n, iel, ifa, isext, xy, tau)
      real*8, intent(in)  :: up(:), uc(:), b(:), n(:), xy(:)
      real, intent(in) :: isext
      integer, intent(in)  :: ifa, iel
      real*8, intent(out) :: tau(:, :)
      real*8              :: tau_aux(4)
      integer             :: ndim
      real*8 :: xc, yc, rad, h, aux, bn, bnorm
      real*8 :: U1, U2, U3, U4
      U1 = uc(1)
      U2 = uc(2)
      U3 = uc(3)
      U4 = uc(4)

      tau = 0.
      ndim = size(n)
      bn = dot_product(b(1:ndim), n)
      bnorm = norm2(b(1:ndim))

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
         else
#endif
            ! Toroidal face
            tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
                tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale               
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
         else
#endif
            ! Toroidal face
            tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
            tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
tau_aux(3) = tau_aux(3) + (phys%diff_e + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
                tau_aux(4) = tau_aux(4) + (phys%diff_ee + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale               
#ifdef TOR3D
         endif
#endif

      elseif (numer%stab == 4) then
         tau_aux = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(1) + sqrt(abs(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)), &
                       phys%lscale/geom%R0*abs(bn)*phys%diff_pari*up(7)**2.5, phys%lscale/geom%R0*abs(bn)*phys%diff_pare*up(8)**2.5)
      else
         write (6, *) "Wrong stabilization type: ", numer%stab
         stop
      endif
      tau(1, 1) = tau_aux(1)
      tau(2, 2) = tau_aux(2)
      tau(3, 3) = tau_aux(3)
      tau(4, 4) = tau_aux(4)
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

!
!                !***********************************************************************
!                !
!                !                           MAGNETIC FIELD
!                !
!                !***********************************************************************
!#ifdef TOR3D
!  SUBROUTINE defineMagneticField(x,y,t,b,divb,drift)
!  real*8, intent(in)      :: x(:),y(:),t(:)
!  real*8, intent(out),optional     :::: b(:,:),divb(:),drift(:,:)
!  real*8                  :: R0,q,r
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
!!  xm = -0.5
!!  ym = -0.5
!  ! Initialization
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
!                                                                                                        ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!                           Br = (yy-ym)
!                           Bz = (-xx+xm)
!                           Bt = 1.
!                                                                                                 divbp = 0.
!                                                                                                CASE(2)
!                                                                                                        IF (.not.switch%axisym) THEN
!                                                                                                                                WRITE(6,*) "This is an axisymmetric test case!"
!                                                                                                                                stop
!                                                                                                 END IF
!                                                                                                        ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!                           Br = (yy-ym)/xx
!                           Bz = (-xx+xm)/xx
!                           Bt = 1.
!                           divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!                                                                                                CASE(3)
!                                                                                                        IF (.not.switch%axisym) THEN
!                                                                                                                                WRITE(6,*) "This is an axisymmetric test case!"
!                                                                                                                                stop
!                                                                                                 END IF
!                                                                                                        ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
!                           Br = (yy-ym)/xx
!                           Bz = (-xx+xm)/xx
!                           Bt = 1.
!                           divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!
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
!

!
!
!
!
!                                                                                                        CASE DEFAULT
!                                                                                                                                        WRITE(6,*) "Error! Test case not valid"
!                                                                                                                                        STOP
!                                                                                        END SELECT
!
!                                                                                        Bp = sqrt(Br**2+Bz**2)
!                                                     BB = sqrt(Bp**2+Bt**2)
!                                                     b(ind,1) = Br/BB
!                                                     b(ind,2) = Bz/BB
!                                                     b(ind,3) = Bt/BB
!                                                     if (present(divb)) then
!                                                        divb(ind) = divbp
!                                                     endif
!                                                 END DO
!                        END DO
!  END SUBROUTINE defineMagneticField
!#else
!  SUBROUTINE defineMagneticField(x,y,b,divb,drift,Bmod)
!  real*8, intent(in)      :: x(:),y(:)
!  real*8, intent(out)     :: b(:,:)
!  real*8, intent(out),optional::divb(:),drift(:,:),Bmod(:)
!  real*8                  :: R0,q,r(size(x)),auxdivb(size(x)),auxdrift(size(b,1),size(b,2))
!  real*8                  :: xmax,xmin,ymax,ymin,xm,ym,xx,yy
!  real*8                  :: Br,Bz,Bt,BB,Bp
!  integer :: i
!  ! Initialization
!  b     = 0.
!  auxdivb(:) = 0.
!  xmax = Mesh%xmax
!  xmin = Mesh%xmin
!  ymax = Mesh%ymax
!  ymin = Mesh%ymin
!  xm = 0.5*(xmax+xmin)
!  ym = 0.5*(ymax+ymin)
!  if (present(Bmod)) then
!     Bmod  = 0.
!  endif
!!  xm = -0.5
!!  ym = -0.5
!                                SELECT CASE(switch%testcase)
!                                  CASE(1)
!      ! Cartesian case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
!!            b(:,1) = 1./30.*(x-y**2+2)
!!            b(:,2) = 1./30.*(x*y+y)
!!            auxdivb(:) = 1./15.+(1./30.)*x
!      DO i=1,size(x)
!             xx = x(i)
!             yy = y(i)
!                           Br = (yy-ym)
!                           Bz = (-xx+xm)
!                           Bt = 1.
!                                                                                                        Bp = sqrt(Br**2+Bz**2)
!                                                                                   BB = sqrt(Bp**2+Bt**2)
!                                                                                   b(i,1) = Br/BB
!                                                                                   b(i,2) = Bz/BB
!             auxdivb(i) = 0.
!      END DO

!      CASE(2)
!!      ! Axisymmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
!!            b(:,1) = 1./30.*(x-y**2+2)
!!            b(:,2) = 1./30.*(x*y+y)
!!            auxdivb(:) = 1./30.+((1./30.)*x-(1./30.)*y**2+1./15.)/x+1./30.*(x+1)
!                                        DO i=1,size(x)
!             xx = x(i)
!             yy = y(i)
!                           Br = (yy-ym)/xx
!                           Bz = (-xx+xm)/xx
!                           Bt = 1.
!                                                                                                        Bp = sqrt(Br**2+Bz**2)
!                                                                                   BB = sqrt(Bp**2+Bt**2)
!                                                                                   b(i,1) = Br/BB
!                                                                                   b(i,2) = Bz/BB
!             auxdivb(i) = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
!            END DO
!      CASE(5)
!      ! Cartesian case, square mesh, horizontal field
!            b(:,1) = 0.1
!            b(:,2) = 0.
!      CASE(6)
!      ! Cartesian case, square mesh, horizontal field
!         DO i=1,size(x)
!            IF (y(i).ge.0.5) THEN
!               b(i,1) = 0.1
!            ELSE
!               b(i,1) = -0.1
!            END IF
!         END DO
!      CASE(50:59)
!         write(6,*) "Error in defineMagneticField: you should not be here!"
!         STOP
!                                  CASE(60:69)
!
!        ! Circular case with limiter
!        R0 = geom%R0
!        q  = geom%q
!        r  = phys%lscale*sqrt((x-R0/phys%lscale)**2+y**2)
!        b(:,1) = -phys%lscale*y/sqrt(R0**2*q**2+(1-q**2)*r**2)
!        b(:,2) = phys%lscale*(x-R0/phys%lscale)/sqrt(R0**2*q**2+(1-q**2)*r**2)

!        IF (switch%axisym) THEN
!            auxdivb(:) = -y/x/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
!        ELSE
!            WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
!            STOP
!        END IF
!
!        IF (switch%driftdia) THEN
!           auxdrift(:,2) =  -1./R0*phys%lscale
!        END IF
!                                  CASE DEFAULT
!                                      WRITE(6,*) "Error! Test case not valid"
!                                      STOP
!                                END SELECT

!    IF (present(divb)) THEN
!       divb = auxdivb
!    ENDIF
!    IF (present(drift)) THEN
!       drift = auxdrift
!    ENDIF
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
!                                                                        phys%b(i,1) = interpolate( ip, yvec,jp, xvec, bx, y,x, 1e-12)
!                                                                        phys%b(i,2) = interpolate( ip, yvec,jp, xvec, by, y,x, 1e-12)
!                                                                        phys%divb(i) = interpolate( ip, yvec,jp, xvec, divb, y,x, 1e-12)
!                                                                        phys%drift(i,1) = interpolate( ip, yvec,jp, xvec,driftx, y,x, 1e-12)
!                                                                        phys%drift(i,2) = interpolate( ip, yvec,jp, xvec,drifty, y,x, 1e-12)
!                                                                        phys%flux2D(i) = interpolate( ip, yvec,jp, xvec,flux2D, y,x, 1e-12)
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
!                                        integer        :: ierr,k
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
!     write(nit, "(i10)") time%it
!     nit = trim(adjustl(nit))
!     k = INDEX(nit, " ") -1
!
!                                        IF (MPIvar%glob_size.GT.1) THEN
!                                           write(nid,*) MPIvar%glob_id+1
!                                           write(npr,*) MPIvar%glob_size
!                                           fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
!                                        ELSE
!                                           fname_complete = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
!                                        END IF
!
!     write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))
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
