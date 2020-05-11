!********************************************
! project: MHDG
! file: prec_const.f90
! date: 06/09/2016
! Definition of real precision and real constants
!********************************************

module prec_const
   implicit none

   ! Precision for real
   !*****************
   integer, parameter :: float = SELECTED_REAL_KIND(13, 99)

   ! Some useful constants
   !*********************
   REAL(float), PARAMETER :: PI = 3.141592653589793238462643383279502884197_float
   REAL(float), PARAMETER :: sqrt2 = 1.414213562373095_float

end module prec_const

