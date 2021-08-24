!**********************************
! project: MHDG
! file: READ_input.f90
! date: 04/09/2016
! Subroutine for input file loading
!**********************************

!********************************
! Loads input file
!********************************
SUBROUTINE READ_input()
  USE prec_const
  USE globals
  USE MPI_OMP
  IMPLICIT NONE

  logical :: driftdia,driftexb, axisym, restart,steady,dotiming,psdtime,decoup,bxgradb
  logical :: ckeramp,saveNR,filter,saveTau,lstiming,fixdPotLim,dirivortcore,dirivortlim,convvort,logrho
  integer :: thresh, difcor, tis, stab,pertini,init
  integer :: itmax, itrace, rest, istop, sollib
  integer :: uinput, printint, testcase, nrp
  integer :: nts, tsw, freqdisp, freqsave, shockcp, limrho
  integer :: bcflags(1:10), ntor, ptor, npartor,bohmtypebc
  real*8  :: dt0, R0, diff_n, diff_u, tau(1:5), tNr, tTM, div, Tbg
  real*8  :: tfi, a, bohmth, q, diffred, diffmin
  real*8  :: sc_coe, so_coe, df_coe, thr, thrpre, minrho, dc_coe, sc_sen
  real*8  :: epn, Mref, diff_pari, diff_e, Gmbohm, Gmbohme
  real*8  :: diff_pare, diff_ee, tie, dumpnr, tmax, tol
  real*8  :: diff_vort, diff_pot, etapar, c1, c2, Potfloat,diagsource(10)
  character(100) :: msg
  character(20)  :: kmethd, ptype

  character(len=20) :: smther, smther2, prol, restr, solve, restr2, prol2, solve2, mlcycle
  character(len=20) :: aggr_prol, par_aggr_alg, aggr_ord, aggr_filter, csolve, csbsolve, cmat
  integer           :: jsweeps, novr, fill, jsweeps2, novr2, fill2, outer_sweeps, maxlevs, csize, cfill, cjswp
  real              :: thrsol, thrsol2, mncrratio, athres, cthres
  real              :: exbdump

  ! RMP and Ripple
  logical     :: RMP, Ripple
  real*8      :: amp_rmp, torElongCoils_rmp, amp_ripple, triang, ellip
  integer     :: nbCoils_rmp,  parite, nbRow, nbCoils_ripple

  ! Neutral and Ohmic heating
  real*8  :: Pohmic,diff_nn,Re,puff

  ! Defining the variables to READ from the file
  NAMELIST /SWITCH_LST/ steady, axisym, init, driftdia, driftexb, testcase, RMP, Ripple, psdtime, diffred, diffmin, &
    & shockcp, limrho, difcor, thresh, filter, decoup, ckeramp, saveNR, saveTau, fixdPotLim, dirivortcore,dirivortlim, convvort,pertini,&
    & logrho,bxgradb
  NAMELIST /NUMER_LST/ tau,nrp,tNR,tTM,div,sc_coe,sc_sen,minrho,so_coe,df_coe,dc_coe,thr,thrpre,stab,dumpnr,ntor,ptor,tmax,npartor,bohmtypebc,exbdump
  NAMELIST /GEOM_LST/ R0, q
  NAMELIST /MAGN_LST/ amp_rmp,nbCoils_rmp,torElongCoils_rmp,parite,nbRow,amp_ripple,nbCoils_ripple,triang,ellip ! RMP and Ripple
  NAMELIST /TIME_LST/ dt0, nts, tfi, tsw, tis
  NAMELIST /PHYS_LST/ diff_n, diff_u, diff_e, diff_ee, diff_vort, diff_nn, Re, puff, Pohmic, Tbg, bcflags, bohmth,&
    &Gmbohm, Gmbohme, a, Mref, tie, diff_pari, diff_pare, diff_pot, epn, etapar, Potfloat,diagsource
  NAMELIST /UTILS_LST/ PRINTint, dotiming, freqdisp, freqsave
  NAMELIST /LSSOLV_LST/ sollib, lstiming, itmax, itrace, rest, istop, tol, kmethd, ptype,&
    &smther, jsweeps,&
    &novr, restr, prol, solve, fill, thrsol, smther2, jsweeps2, novr2, restr2, prol2, solve2, fill2, thrsol2, mlcycle,&
    &outer_sweeps, maxlevs, csize, aggr_prol, par_aggr_alg, aggr_ord, aggr_filter, mncrratio, athres,&
    &csolve, csbsolve, cmat, cfill, cthres, cjswp

  ! Reading the file
  uinput = 100
  diagsource = 0.
  OPEN (uinput, file='param.txt', status='unknown')
  READ (uinput, SWITCH_LST)
  READ (uinput, NUMER_LST)
  READ (uinput, GEOM_LST)
  READ (uinput, MAGN_LST)
  READ (uinput, TIME_LST)
  READ (uinput, PHYS_LST)
  READ (uinput, UTILS_LST)
  READ (uinput, LSSOLV_LST)
  CLOSE (uinput)

  ! Storing at the right place
  switch%steady           = steady
  switch%axisym           = axisym
  switch%init             = init
  switch%driftdia         = driftdia
  switch%driftexb         = driftexb
  switch%testcase         = testcase
  switch%RMP              = RMP
  switch%Ripple           = Ripple
  switch%psdtime          = psdtime
  switch%diffred          = diffred
  switch%diffmin          = diffmin
  switch%shockcp          = shockcp
  switch%limrho           = limrho
  switch%difcor           = difcor
  switch%thresh           = thresh
  switch%filter           = filter
  switch%decoup           = decoup
  switch%ckeramp          = ckeramp
  switch%saveNR           = saveNR
  switch%saveTau          = saveTau
  switch%fixdPotLim       = fixdPotLim
  switch%dirivortcore     = dirivortcore
  switch%dirivortlim      = dirivortlim
  switch%convvort         = convvort
  switch%pertini          = pertini
  switch%logrho           = logrho
  switch%bxgradb          = bxgradb
  numer%tau               = tau
  numer%nrp               = nrp
  numer%tNR               = tNR
  numer%tTM               = tTM
  numer%div               = div
  numer%sc_coe            = sc_coe
  numer%sc_sen            = sc_sen
  numer%minrho            = minrho
  numer%so_coe            = so_coe
  numer%df_coe            = df_coe
  numer%dc_coe            = dc_coe
  numer%thr               = thr
  numer%thrpre            = thrpre
  numer%stab              = stab
  numer%dumpnr            = dumpnr
  numer%ntor              = ntor
  numer%ptor              = ptor
  numer%tmax              = tmax
  numer%npartor           = npartor
  numer%bohmtypebc        = bohmtypebc
  numer%exbdump           = exbdump
  geom%R0                 = R0
  geom%q                  = q
  magn%amp_rmp            = amp_rmp
  magn%nbCoils_rmp        = nbCoils_rmp
  magn%torElongCoils_rmp  = torElongCoils_rmp
  magn%parite             = parite
  magn%nbRow              = nbRow
  magn%amp_ripple         = amp_ripple
  magn%nbCoils_ripple     = nbCoils_ripple
  magn%triang             = triang
  magn%ellip              = ellip
  time%dt0                = dt0
  time%nts                = nts
  time%tfi                = tfi
  time%tsw                = tsw
  time%tis                = tis
  phys%diff_n             = diff_n
  phys%diff_u             = diff_u
  phys%diff_e             = diff_e
  phys%diff_ee            = diff_ee
  phys%diff_vort          = diff_vort
  phys%diff_nn            = diff_nn
  phys%Re                 = Re
  phys%puff               = puff
  phys%Pohmic             = Pohmic
  phys%Tbg                = Tbg
  phys%bcflags            = bcflags
  phys%bohmth             = bohmth
  phys%Gmbohm             = Gmbohm
  phys%Gmbohme            = Gmbohme
  phys%a                  = a
  phys%Mref               = Mref
  phys%tie                = tie
  phys%diff_pari          = diff_pari
  phys%diff_pare          = diff_pare
  phys%diff_pot           = diff_pot
  phys%epn                = epn
  phys%etapar             = etapar
  phys%Potfloat           = Potfloat
  phys%diagsource         = diagsource
  utils%PRINTint          = printint
  utils%timing            = dotiming
  utils%freqdisp          = freqdisp
  utils%freqsave          = freqsave
  lssolver%sollib         = sollib
  lssolver%timing         = lstiming
  lssolver%itmax          = itmax
  lssolver%itrace         = itrace
  lssolver%rest           = rest
  lssolver%istop          = istop
  lssolver%tol            = tol
  lssolver%kmethd         = kmethd
  lssolver%ptype          = ptype
  lssolver%smther         = smther
  lssolver%jsweeps        = jsweeps
  lssolver%novr           = novr
  lssolver%restr          = restr
  lssolver%prol           = prol
  lssolver%solve          = solve
  lssolver%fill           = fill
  lssolver%thr            = thrsol
  lssolver%smther2        = smther2
  lssolver%jsweeps2       = jsweeps2
  lssolver%novr2          = novr2
  lssolver%restr2         = restr2
  lssolver%prol2          = prol2
  lssolver%solve2         = solve2
  lssolver%fill2          = fill2
  lssolver%thr2           = thrsol2
  lssolver%mlcycle        = mlcycle
  lssolver%outer_sweeps   = outer_sweeps
  lssolver%maxlevs        = maxlevs
  lssolver%csize          = csize
  lssolver%aggr_prol      = aggr_prol
  lssolver%par_aggr_alg   = par_aggr_alg
  lssolver%aggr_ord       = aggr_ord
  lssolver%aggr_filter    = aggr_filter
  lssolver%mncrratio      = mncrratio
  lssolver%athres         = athres
  lssolver%csolve         = csolve
  lssolver%csbsolve       = csbsolve
  lssolver%cmat           = cmat
  lssolver%cfill          = cfill
  lssolver%cthres         = cthres
  lssolver%cjswp          = cjswp

  IF (switch%steady) THEN
    msg = 'Steady state simulation'
  ELSEIF (switch%psdtime) THEN
    msg = 'Pseudotime simulation for reducing diffusion'
  ELSE
    msg = 'Time advancing simulation'
  END IF

  ! Some checking of the inputs
  if (time%tis .gt. 6) then
    write (6, *) "Error: wrong time integration scheme in parameters: tis=", time%tis
    stop
  end if
  if (numer%dumpnr < 0. .or. numer%dumpnr > 1.) then
    write (6, *) "Error: wrong dumping factor for Newton-Raphson: dumpnr=", numer%dumpnr
    stop
  end if
  if (switch%testcase.ge.50 .and. .not.switch%axisym) then
    write (6, *) "Error: this should be an axisymmetric simulation"
    stop
  endif
  ! Adimensionalize units and store reference values
  CALL adimensionalization()

  ! A little message for the user...
  if (MPIvar%glob_id .eq. 0) THEN
    PRINT *, '                                                                      '
    PRINT *, '                                                                      '
#ifdef TOR3D
    PRINT *, '------------------------  MHDG SIMULATION IN 3D ----------------------'
#else
    PRINT *, '------------------------  MHDG SIMULATION IN 2D ----------------------'
#endif
    PRINT *, '                                                                      '
#ifndef TEMPERATURE
#ifdef NEUTRAL
    PRINT *, ' MODEL: N-Gamma isothermal with neutral                               '
#else
    PRINT *, ' MODEL: N-Gamma isothermal                                            '
#endif
#else
#ifdef NEUTRAL
    PRINT *, ' MODEL: N-Gamma-Ti-Te with neutral                                    '
#else
    PRINT *, ' MODEL: N-Gamma-Ti-Te                                                 '
#endif
#endif
    PRINT *, '                                                                      '
    PRINT *, '----------------------------------------------------------------------'
    PRINT *, '  Simulation type: ', Adjustl(Trim(msg))
    PRINT *, '----------------------------------------------------------------------'
    PRINT *, 'Parameter file loaded:'
    PRINT *, '        ***************** Geometry ****************************'
    PRINT *, '                - R0:                                                 ', R0
    PRINT *, '                - Security factor:                                    ', q
    PRINT *, '      	***************** Magnetic ****************************'
    PRINT *, '		            - RMP amplitude:                                      ', amp_rmp
    PRINT *, '		            - RMP coil number:                                    ', nbCoils_rmp
    PRINT *, '		            - RMP coil toroidal elongation:                       ', torElongCoils_rmp
    PRINT *, '		            - RMP parity:                                         ', parite
    PRINT *, '		            - RMP row number:                                     ', nbRow
    PRINT *, '		            - Ripple amplitude:                                   ', amp_ripple
    PRINT *, '		            - Ripple coil number:                                 ', nbCoils_ripple
    PRINT *, '		            - triangularity:                                      ', triang
    PRINT *, '		            - ellipticity:                                        ', ellip
    PRINT *, '        ***************** Time stepping ************************'
    PRINT *, '                - dt0:                                                ', time%dt0
    PRINT *, '                - dt modification:                                    ', time%tsw
    PRINT *, '                - final time:                                         ', time%tfi
    PRINT *, '                - max number of time steps:                           ', time%nts
    PRINT *, '                - time integration scheme:                            ', time%tis
    PRINT *, '        ***************** Physics *****************************'
    PRINT *, '                - perp. diffusion in the continuity equation:         ', phys%diff_n
    PRINT *, '                - perp. diffusion in the momentum equation:           ', phys%diff_u
#ifdef TEMPERATURE
    PRINT *, '                - perp. diffusion in the ions energy equation:        ', phys%diff_e
    PRINT *, '                - perp. diffusion in the electrons energy equation:   ', phys%diff_ee
    PRINT *, '                - paral. diffusion for ions temperature:              ', phys%diff_pari
    PRINT *, '                - paral. diffusion for electrons temperature:         ', phys%diff_pare
    PRINT *, '                - temperature exchange coefficient ions/electrons:    ', phys%tie
    PRINT *, '                - temperature diffusion exponential:                  ', phys%epn
    PRINT *, '                - reference Mach number:                              ', phys%Mref
    PRINT *, '                - gamma for Bohm boundary condition on ions:          ', phys%Gmbohm
    PRINT *, '                - gamma for Bohm boundary condition for electrons:    ', phys%Gmbohme
    PRINT *, '                - Ohmic heating (only for testcase=54 and core):      ', phys%Pohmic
#endif
#ifdef VORTICITY
    PRINT *, '                - perp. diffusion in the vorticity equation:          ', phys%diff_vort
    PRINT *, '                - perp. diffusion in the potential equation:          ', phys%diff_pot
#endif
#ifdef NEUTRAL
    PRINT *, '                - diffusion in the neutral equation:                  ', phys%diff_nn
    PRINT *, '                - recycling coefficient in the neutral equation:      ', phys%Re
    PRINT *, '                - puff coefficient in the neutral equation:           ', phys%puff
#endif
    PRINT *, '                - constant for the momentum equation (isoth)          ', phys%a
    PRINT *, '                - diagonal implicit sources                           ', phys%diagsource
    PRINT *, '        ***************** Switches ****************************'
    PRINT *, '                - stady state simulation:                             ', switch%steady
    PRINT *, '                - axisym:                                             ', switch%axisym
    if (switch%init.eq.1) then
      PRINT *, '  Initializing with analytical solution at nodes                    '
    else
      PRINT *, '  Initializing with L2 projection                                   '
    endif
    PRINT *, '                - driftdia:                                           ', driftdia
    PRINT *, '                - driftexb:                                           ', driftexb
    PRINT *, '                - test case:                                          ', testcase
    PRINT *, '                - RMP:                                                ', RMP
    PRINT *, '                - Ripple:                                             ', Ripple
    PRINT *, '                - shockcp:                                            ', shockcp
    PRINT *, '                - minrho:                                             ', minrho
    PRINT *, '                - logrho:                                             ', logrho
    PRINT *, '                - thresh:                                             ', thresh
    PRINT *, '                - filter:                                             ', filter
#ifdef TEMPERATURE
    PRINT *, '                - decoup:                                             ', decoup
#endif
    PRINT *, '                - ckeramp:                                            ', ckeramp
    PRINT *, '                - saveNR:                                             ', saveNR
    PRINT *, '                - saveTau:                                            ', saveTau
    PRINT *, '        ***************** Numerics ****************************'
    PRINT *, '                - stabilization type:                                 ', numer%stab
    PRINT *, '                - tau(1):                                             ', numer%tau(1)
    PRINT *, '                - tau(2):                                             ', numer%tau(2)
    PRINT *, '                - tau(3):                                             ', numer%tau(3)
    PRINT *, '                - tau(4):                                             ', numer%tau(4)
    PRINT *, '                - tau(5):                                             ', numer%tau(5)
    PRINT *, '                - max number of N-R iterations:                       ', numer%nrp
    PRINT *, '                - tolerance for the N-R scheme:                       ', numer%tNR
    PRINT *, '                - tolerance for the steady state achievement:         ', numer%tTM
    IF (switch%shockcp .gt. 0) THEN
      PRINT *, '                - shock capturing coeff:                            ', numer%sc_coe
      PRINT *, '                - shock capturing sensibility:                      ', numer%sc_sen
    END IF
    IF (switch%limrho .gt. 0) THEN
      PRINT *, '                - applying limiting of rho at value:                ', numer%minrho
    END IF
    IF (switch%limrho .eq. 1 .or. switch%limrho .eq. 3) THEN
      PRINT *, '                - coefficient of the source for limiting rho:       ', numer%so_coe
    END IF
    IF (switch%limrho .eq. 2 .or. switch%limrho .eq. 3) THEN
      PRINT *, '                - coefficient of the diffusion for limiting rho:    ', numer%df_coe
    END IF
    IF (switch%difcor .gt. 0) THEN
      PRINT *, '                - adding diffusion in corners, in position:         ', switch%difcor
      PRINT *, '                - diffusion coefficient in corners:                 ', numer%dc_coe
    END IF
    IF (switch%thresh .gt. 0) THEN
      PRINT *, '                - using a threshold at rho:                         ', numer%thr
      PRINT *, '                - using a threshold at pressure:                    ', numer%thrpre
    END IF
    PRINT *, '                - using a dumping factor for Newton-Raphson:          ', numer%dumpnr
#ifdef TOR3D
    PRINT *, '                - number of elements in the toroidal direction:       ', numer%ntor
    PRINT *, '                - polynomial degree in the toroidal direction:        ', numer%ptor
#ifdef PARALL
    PRINT *, '                - number of MPI partitions in the toroidal direction :', numer%npartor
#endif
    PRINT *, '                - max extention in the toroidal direction:            ', numer%tmax
#endif
    PRINT *, '                - dumping for ExB term:       ', numer%exbdump
    PRINT *, '        ***************** Linear solver params******************'
    IF (lssolver%sollib == 1) THEN
      PRINT *, '                - Library used for the linear system:   PASTIX      '
    ELSEIF (lssolver%sollib == 2) THEN
      PRINT *, '                - Library used for the linear system:   PSBLAS      '
      PRINT *, '                - Iterative method:                                 ', lssolver%kmethd
      PRINT *, '                - Preconditioner:                                   ', lssolver%ptype
      PRINT *, '                - Stopping criterion type                           ', lssolver%istop
      PRINT *, '                - Stopping criterion tolerance                      ', lssolver%tol
      PRINT *, '                - Restart:                                          ', lssolver%rest
    ENDIF

    PRINT *, '        '
  end if
END SUBROUTINE READ_input

