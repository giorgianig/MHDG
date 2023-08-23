!**********************************
! project: MHDG
! file: adimensionalization.f90
! date: 21/04/2020
! Adimensionalization of the inputs
!**********************************
SUBROUTINE adimensionalization()
  USE globals
  USE MPI_OMP
  IMPLICIT NONE

  ! Parameters
  real*8, parameter :: mi = 3.35e-27         ! Ionic mass [kg]
  real*8, parameter :: me = 9.109e-31        ! Electronic mass [kg]
  real*8, parameter :: e = 1.60217662e-19    ! Electron charge [C]
  real*8, parameter :: kB = 1.38064852e-23   ! Boltzmann constant [m^2*kg*s^-2*K^-1]
  real*8, parameter :: eps0 = 8.85e-12       ! Permittivity of free space [F/m]

  ! Reference values
  real*8, parameter :: L0 = 1.901e-3         ! Length scale [m]
  real*8, parameter :: t0 = 1.374e-07        ! Time scale [s]
  ! *** Temperature scale [eV]. It corresponds to the background temperature
  real*8           :: Tev                    ! Temperature scale [eV].

  real*8, parameter :: n0 = 1e19             ! Reference density [m^-3]

  ! Derived reference values
  real*8           :: u0                     ! Reference speed [m*s^-1]
  real*8           :: B0                     ! Reference magnetic field [kg*C^-1*s^-1]
  real*8           :: phi0                   ! Reference electric potential [kg*m^2*C^-1*s^-2]
  real*8           :: W0                     ! Reference vorticity [kg^-1*C]
  real*8           :: J0                     ! Reference current density [C*s^-1*m^-2]
  real*8           :: D0                     ! Reference diffusion [m^2*s^-1]

  ! Values from Stengby
  real*8, parameter :: k0 = 2000             ! Stengby

  ! Other reals
  real*8           :: k_star, coef, tau_ie

  ! Set the reference temperature to the background one
#ifdef TEMPERATURE
  Tev = 50
#else
  Tev = phys%Tbg
#endif

  ! Initialize simpar to adimensional case
  simpar%refval_mass = 1.
  simpar%refval_charge = 1.
  simpar%refval_length = 1.
  simpar%refval_time = 1.
  simpar%refval_temperature = 1.
  simpar%refval_density = 1.
  simpar%refval_neutral = 1.
  simpar%refval_speed = 1.
  simpar%refval_potential = 1.
  simpar%refval_vorticity = 1.
  simpar%refval_magfield = 1.
  simpar%refval_current = 1.
  simpar%refval_diffusion = 1.
  simpar%refval_momentum = 1.
  simpar%refval_specpress = 1.
  simpar%refval_specenergy = 1.
  simpar%refval_specenergydens = 1.
  simpar%refval_length_dimensions = '-'
  simpar%refval_time_dimensions = '-'
  simpar%refval_temperature_dimensions = '-'
  simpar%refval_density_dimensions = '-'
  simpar%refval_neutral_dimensions = '-'
  simpar%refval_speed_dimensions = '-'
  simpar%refval_potential_dimensions = '-'
  simpar%refval_vorticity_dimensions = '-'
  simpar%refval_magfield_dimensions = '-'
  simpar%refval_current_dimensions = '-'
  simpar%refval_diffusion_dimensions = '-'
  simpar%refval_momentum_dimensions = '-'
  simpar%refval_specpress_dimensions = '-'
  simpar%refval_specenergy_dimensions = '-'
  simpar%refval_specenergydens_dimensions = '-'

  ! Some scaling coefficients
  phys%lscale = 1.
  phys%B0 = 1.
  phys%dfcoef = 1.
  phys%c1 = 1.
  phys%c2 = 1.
  phys%dexbcoef = 1.


  IF (switch%testcase < 10) THEN
    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) "No adimensionalization needed"
    ENDIF
    RETURN
  ELSE
    phys%Mref = 1.
  ENDIF
  IF (MPIvar%glob_id .eq. 0) THEN
    WRITE (6, *) "Adimensionalizing input values"
  ENDIF
  ! Computing derived reference values
  u0 = L0/t0
  D0 = L0**2/t0
  W0 = e/mi
  B0 = (W0*t0)**(-1)
  phi0 = u0**2/W0
  J0 = mi*W0/t0/L0**2

  ! Reference values stored in phys
  phys%lscale = L0
  phys%B0 = B0
  !
  ! Adimesional isothermal compressibility coefficient
  phys%a = 2*Tev*e/mi/(L0/t0)**2

  ! Adimensional non-isothermal compressibility coefficient
  phys%Mref = Tev*e/(mi*u0**2) !(0.5*phys%a)

  ! Constants for vorticity (Mref = C2/C1)
  phys%c1 = 1/(n0*t0*mi*W0**2)
  phys%c2 = Tev*t0/(L0**2*n0*mi*W0)

  ! Parallel temperature diffusion coefficient for non-isothermal model
  ! tau_ie = mi/me*2.4/3.*1e10*Tev**(0.5)*u0**2/n0/t0
  k_star = t0**3*Tev**(7./2.)/(n0*L0**4)*k0/mi
  coef = 3*sqrt2/e**4/12.*eps0**2*pi**1.5*mi/me*sqrt(me)*e**1.5
  !tau_ie = mi/me*2.4/3.*1e10*sqrt(Tev)*u0**2/n0/t0
  !tau_ie =(12./15.)*(kB/e)*2./3.*sqrt(Tev)*mi*u0**2/(n0*t0*kB)
  tau_ie = -(12./15.)*(kB/e)*2./3.*coef*sqrt(Tev)*mi*u0**2/(n0*t0*kB)
  !   tau_ie = 2./3.*coef*Tev**(0.5)*u0**2/n0/t0*mi/k0
  phys%diff_pari = k_star/33.333333333333336
  phys%diff_pare = k_star

  ! Temperature exchange characteristic time
  phys%tie = tau_ie

  ! Diffusion coefficients
  phys%diff_n = phys%diff_n/D0
  phys%diff_u = phys%diff_u/D0
  phys%diff_e = phys%diff_e/D0
  phys%diff_ee = phys%diff_ee/D0
  phys%diff_vort = phys%diff_vort/D0
  phys%diff_pot = phys%diff_pot/D0
  phys%diff_nn = phys%diff_nn/D0
  
  ! Pinch velocity
  phys%v_p = phys%v_p/u0

  ! Curvature drift coefficient
  phys%dfcoef = 2*Tev*t0/(L0**2*B0)
  phys%dexbcoef = phi0/B0*t0/L0**2

  ! Store reference values
  simpar%refval_mass = mi
  simpar%refval_charge = e
  simpar%refval_length = L0
  simpar%refval_time = t0
  simpar%refval_temperature = Tev
  simpar%refval_density = n0
  simpar%refval_neutral = n0
  simpar%refval_speed = u0
  simpar%refval_potential = phi0
  simpar%refval_vorticity = W0
  simpar%refval_magfield = B0
  simpar%refval_current = J0
  simpar%refval_diffusion = D0
  simpar%refval_momentum = L0*n0
  simpar%refval_specpress = Tev*n0/mi*e
  simpar%refval_specenergy = u0**2
  simpar%refval_specenergydens = u0**2*n0
  simpar%refval_length_dimensions = 'm'
  simpar%refval_time_dimensions = 's'
  simpar%refval_temperature_dimensions = 'eV'
  simpar%refval_mass_dimensions = 'kg'
  simpar%refval_density_dimensions = 'm^-3'
  simpar%refval_speed_dimensions = 'm*s^-1'
  simpar%refval_potential_dimensions = 'kg*m^2*A^-1*s^-3'
  simpar%refval_vorticity_dimensions = 'kg^-1*A*s'
  simpar%refval_magfield_dimensions = 'kg*A^-1*s^-2'
  simpar%refval_current_dimensions = 'A*m^-2'
  simpar%refval_diffusion_dimensions = 'm^2*s^-1'
  simpar%refval_momentum_dimensions = 'm^-2*s^-1'
  simpar%refval_specpress_dimensions = 'm^-1*s^-2'
  simpar%refval_specenergy_dimensions = 'm^2*s^-2'
  simpar%refval_specenergydens_dimensions = 'm^-1*s^-2'

END SUBROUTINE adimensionalization
