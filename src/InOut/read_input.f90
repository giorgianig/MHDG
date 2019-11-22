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
   
   logical :: driftvel, axisym, restart,steady,timing,psdtime,decoup,ckeramp,saveNR,filter,saveTau,lstiming
   integer :: thresh,difcor,tis, stab
   integer :: itmax,itrace,rest,istop,sollib
   integer :: uinput, init_cond,printint,testcase,nrp
   integer :: nts,tsw,freqdisp,freqsave,shockcp,limrho
   integer :: bcflags(1:10),ntor,ptor,npartor
   real*8  :: dt0, R0, diff_n, diff_u, tau(1:4),tNr,tTM,div
   real*8  :: tfi,a,bohmth,lscale,q,diffred,diffmin,dfcoef
   real*8  :: sc_coe,so_coe,df_coe,thr,thrpre,minrho,dc_coe,sc_sen
   real*8  :: epn,Mref,diff_pari,diff_e,Gmbohm,Gmbohme
   real*8  :: diff_pare,diff_ee,tie,dumpnr,tmax,tol
   character(100) :: msg
   character(20)  :: kmethd,ptype
   
   character(len=20) :: smther,smther2,prol,restr,solve,restr2,prol2,solve2,mlcycle
   character(len=20) :: aggr_prol,par_aggr_alg,aggr_ord,aggr_filter,csolve,csbsolve,cmat
   integer           :: jsweeps,novr,fill,jsweeps2,novr2,fill2,outer_sweeps,maxlevs,csize,cfill,cjswp
   real              :: thrsol,thrsol2,mncrratio,athres,cthres
      
   ! Defining the variables to READ from the file
   NAMELIST /SWITCH_LST/ restart,steady,axisym,init_cond,driftvel,testcase,psdtime,diffred,diffmin, &
                       & shockcp,limrho,difcor,thresh,filter,decoup,ckeramp,saveNR,saveTau
   NAMELIST /NUMER_LST/ tau,nrp,tNR,tTM,div,sc_coe,sc_sen,minrho,so_coe,df_coe,dc_coe,thr,thrpre,stab,dumpnr,ntor,ptor,tmax,npartor
   NAMELIST /GEOM_LST/ R0,q
   NAMELIST /TIME_LST/ dt0,nts, tfi,tsw,tis
   NAMELIST /PHYS_LST/ lscale,diff_n,diff_u,a,bcflags,bohmth,dfcoef,diff_e,epn,Mref,diff_pari,Gmbohm,diff_ee,diff_pare,tie,Gmbohme
   NAMELIST /UTILS_LST/ PRINTint, timing,freqdisp,freqsave
   NAMELIST /LSSOLV_LST/ sollib,lstiming,itmax,itrace,rest,istop,tol,kmethd,ptype,&
   &smther,jsweeps,&
   &novr,restr,prol,solve,fill,thrsol,smther2,jsweeps2,novr2,restr2,prol2,solve2,fill2,thrsol2,mlcycle,&
   &outer_sweeps,maxlevs,csize,aggr_prol,par_aggr_alg,aggr_ord,aggr_filter,mncrratio,athres,&
   &csolve,csbsolve,cmat,cfill,cthres,cjswp
   
    ! Reading the file
    uinput = 100
    OPEN(uinput, file='param.txt', status='unknown')
    READ(uinput,SWITCH_LST)
    READ(uinput,NUMER_LST)
    READ(uinput,GEOM_LST)
    READ(uinput,TIME_LST)
    READ(uinput,PHYS_LST)
    READ(uinput,UTILS_LST)
    READ(uinput,LSSOLV_LST)
    CLOSE(uinput)

   ! Storing at the right place
   switch%restart   = restart
   switch%steady    = steady
   switch%axisym    = axisym
   switch%driftvel  = driftvel  
   switch%init_cond = init_cond  
   switch%testcase  = testcase  
   switch%psdtime   = psdtime
   switch%diffred   = diffred
   switch%diffmin   = diffmin
   switch%shockcp   = shockcp   
			switch%limrho    = limrho   
			switch%thresh    = thresh   
			switch%filter    = filter   
			switch%difcor    = difcor
			switch%decoup    = decoup
			switch%ckeramp   = ckeramp
			switch%saveNR    = saveNR
			switch%saveTau   = saveTau
   numer%tau        = tau 
   numer%nrp        = nrp
   numer%tNR        = tNR
   numer%tTM        = tTM
   numer%div        = div
   numer%sc_coe     = sc_coe
   numer%sc_sen     = sc_sen
   numer%minrho     = minrho   
   numer%so_coe     = so_coe
   numer%df_coe     = df_coe
   numer%dc_coe     = dc_coe
   numer%thr        = thr
   numer%thrpre     = thrpre
   numer%stab       = stab
   numer%dumpnr     = dumpnr
   numer%ntor       = ntor
   numer%ptor       = ptor
   numer%tmax       = tmax
   numer%npartor    = npartor
   geom%R0          = R0
   geom%q           = q
   time%dt0         = dt0
   time%tfi         = tfi
   time%nts         = nts
   time%tsw         = tsw
   time%tis         = tis
   phys%lscale      = lscale
   phys%diff_n      = diff_n
   phys%diff_u      = diff_u
   phys%diff_e      = diff_e   
   phys%a           = a
   phys%bcflags     = bcflags
   phys%bohmth      = bohmth
   phys%dfcoef      = dfcoef
   phys%epn         = epn
   phys%Mref        = Mref
   phys%diff_pari   = diff_pari
   phys%Gmbohm      = Gmbohm
   phys%Gmbohme     = Gmbohme
   phys%diff_pare   = diff_pare
   phys%diff_ee     = diff_ee
   phys%tie         = tie
   utils%PRINTint   = printint
   utils%timing     = timing
   utils%freqdisp   = freqdisp
   utils%freqsave   = freqsave
   lssolver%sollib  = sollib
   lssolver%timing  = lstiming
   lssolver%itmax   = itmax
   lssolver%itrace  = itrace
   lssolver%rest    = rest
   lssolver%istop   = istop
   lssolver%tol     = tol
   lssolver%kmethd  = kmethd
   lssolver%ptype   = ptype
   lssolver%smther  = smther
   lssolver%jsweeps = jsweeps
   lssolver%novr    = novr
   lssolver%restr   = restr
   lssolver%prol    = prol
   lssolver%solve   = solve
   lssolver%fill    = fill
   lssolver%thr     = thrsol
   lssolver%smther2 = smther2
   lssolver%jsweeps2 = jsweeps2
   lssolver%novr2   = novr2
   lssolver%restr2  = restr2
   lssolver%prol2   = prol2
   lssolver%solve2  = solve2
   lssolver%fill2   = fill2
   lssolver%thr2    = thrsol2
   lssolver%mlcycle = mlcycle
   lssolver%outer_sweeps = outer_sweeps
   lssolver%maxlevs = maxlevs
   lssolver%csize   = csize
   lssolver%mncrratio = mncrratio
   lssolver%athres  = athres
   lssolver%aggr_prol = aggr_prol
   lssolver%par_aggr_alg = par_aggr_alg
   lssolver%aggr_ord = aggr_ord
   lssolver%aggr_filter = aggr_filter
   lssolver%csolve  = csolve
   lssolver%csbsolve = csbsolve
   lssolver%cmat    = cmat 
   lssolver%cfill   = cfill
   lssolver%cthres  = cthres
   lssolver%cjswp   = cjswp
     
  IF (switch%steady) THEN 
     msg = 'Steady state simulation'
  ELSEIF (switch%psdtime) THEN
     msg = 'Pseudotime simulation for reducing diffusion'
  ELSE
     msg = 'Time advancing simulation'
  END IF
 
	 ! Some checking of the inputs
		if (time%tis.ne.1 .and. time%tis.ne.2) then
		   write(6,*) "Error: wrong time integration scheme in parameters: tis=", time%tis
		   stop
		end if
		if (numer%dumpnr<0. .or.numer%dumpnr>1.) then
		   write(6,*) "Error: wrong dumping factor for Newton-Raphson: dumpnr=", numer%dumpnr
		   stop
		end if  
   ! A little message for the user... 
  if (MPIvar%glob_id.eq.0) THEN  
  PRINT*,'                                                            '
  PRINT*,'                                                            '    
#ifdef TOR3D  
  PRINT*,'-------------------  MHDG SIMULATION IN 3D -----------------'
#else
  PRINT*,'-------------------  MHDG SIMULATION IN 2D  -----------------'
#endif  
  PRINT*,'                                                            '
#ifndef TEMPERATURE
  PRINT*,' MODEL: N-Gamma isothermal                                  '  
#else
  PRINT*,' MODEL: N-Gamma-Ti-Te                                       '    
#endif
  PRINT*,'                                                            '
  PRINT*,'------------------------------------------------------------'
  PRINT*,'  Simulation type: ',Adjustl(Trim(msg))
  PRINT*,'------------------------------------------------------------'
		PRINT*,'Parameter file loaded:'
		PRINT*,'	***************** Geometry ****************************'
		PRINT*,'		- R0:                                               ', R0
		PRINT*,'		- Security factor:                                  ', q
		PRINT*,'	***************** Time stepping ************************'
		PRINT*,'		- dt0:                                              ', time%dt0
		PRINT*,'		- dt modification:                                  ', time%tsw
		PRINT*,'		- final time:                                       ', time%tfi
		PRINT*,'		- max number of time steps:                         ', time%nts
		PRINT*,'		- time integration scheme:                          ', time%tis
		PRINT*,'	***************** Physics *****************************'
		PRINT*,'		- length scale:                                     ', lscale
		PRINT*,'		- perp. diffusion in the continuity equation:       ', diff_n
		PRINT*,'		- perp. diffusion in the momentum equation:         ', diff_u  
#ifdef TEMPERATURE
		PRINT*,'		- perp. diffusion in the ions energy equation:      ', diff_e
		PRINT*,'		- perp. diffusion in the electrons energy equation: ', diff_ee
		PRINT*,'		- paral. diffusion for ions temperature:            ', diff_pari
		PRINT*,'		- paral. diffusion for electrons temperature:       ', diff_pare
		PRINT*,'		- temperature exchange coefficient ions/electrons:  ', tie
		PRINT*,'		- temperature diffusion exponential:                ', epn
		PRINT*,'		- reference Mach number:                            ', Mref
  PRINT*,'		- gamma for Bohm boundary condition on ions:        ', Gmbohm
  PRINT*,'		- gamma for Bohm boundary condition for electrons:  ', Gmbohme
#endif		
		PRINT*,'		- constant for the momentum equation (isoth)        ', a    
		PRINT*,'		- constant for the drift velocity:                  ', dfcoef    
		PRINT*,'	***************** Switches ****************************'
		PRINT*,'		- restart:                                          ', switch%restart
		PRINT*,'		- stady state simulation:                           ', switch%steady
		PRINT*,'		- initialization:                                   ', switch%init_cond
		PRINT*,'		- axisym:                                           ', switch%axisym
		PRINT*,'		- driftvel:                                         ', driftvel
		PRINT*,'		- test case:                                        ', testcase
  PRINT*,'		- shockcp:                                          ', shockcp
  PRINT*,'		- minrho:                                           ', minrho
  PRINT*,'		- thresh:                                           ', thresh
  PRINT*,'		- filter:                                           ', filter  		    		
#ifdef TEMPERATURE		  
  PRINT*,'		- decoup:                                           ', decoup  		  
#endif  
  PRINT*,'		- ckeramp:                                          ', ckeramp
  PRINT*,'  - saveNR:                                           ', saveNR
  PRINT*,'  - saveTau:                                          ', saveTau
		PRINT*,'	***************** Numerics ****************************'
		PRINT*,'		- stabilization type:                               ', numer%stab
		PRINT*,'		- tau(1):                                           ', numer%tau(1)
		PRINT*,'		- tau(2):                                           ', numer%tau(2)
		PRINT*,'		- tau(3):                                           ', numer%tau(3)
		PRINT*,'		- tau(4):                                           ', numer%tau(4)								
		PRINT*,'		- max number of N-R iterations:                     ', numer%nrp
		PRINT*,'		- tolerance for the N-R scheme:                     ', numer%tNR		
		PRINT*,'		- tolerance for the steady state achievement:       ', numer%tTM
  IF (switch%shockcp.gt.0) THEN
  PRINT*,'		- shock capturing coeff:                            ', numer%sc_coe
  PRINT*,'		- shock capturing sensibility:                      ', numer%sc_sen
  END IF		
  IF (switch%limrho.gt.0) THEN 
  PRINT*,'		- applying limiting of rho at value:                ', numer%minrho
  END IF
  IF (switch%limrho.eq.1 .or. switch%limrho.eq.3) THEN
  PRINT*,'		- coefficient of the source for limiting rho:       ', numer%so_coe
  END IF		
  IF (switch%limrho.eq.2 .or. switch%limrho.eq.3) THEN
  PRINT*,'		- coefficient of the diffusion for limiting rho:    ', numer%df_coe
  END IF
  IF (switch%difcor.gt.0) THEN
  PRINT*,'		- adding diffusion in corners, in position:         ', switch%difcor
  PRINT*,'		- diffusion coefficient in corners:                 ', numer%dc_coe
  END IF		  
  IF (switch%thresh.gt.0) THEN
  PRINT*,'		- using a threshold at rho:                         ', numer%thr
  PRINT*,'		- using a threshold at pressure:                    ', numer%thrpre
  END IF
  PRINT*,'		- using a dumping factor for Newton-Raphson:        ', numer%dumpnr
#ifdef TOR3D
  PRINT*,'		- number of elements in the toroidal direction:     ', numer%ntor
  PRINT*,'		- polynomial degree in the toroidal direction:      ', numer%ptor
#ifdef PARALL
  PRINT*,'		- number of MPI partitions in the toroidal direction:', numer%ptor
#endif  
  PRINT*,'		- max extention in the toroidal direction:          ', numer%tmax    
#endif
		PRINT*,'	***************** Linear solver params******************'
  IF (lssolver%sollib==1) THEN
  PRINT*,'		- Library used for the linear system:   PASTIX '
  ELSEIF		(lssolver%sollib==2) THEN
  PRINT*,'		- Library used for the linear system:   PSBLAS '
		PRINT*,'		- Iterative method:                            ', lssolver%kmethd
		PRINT*,'		- Preconditioner:                              ', lssolver%ptype
		PRINT*,'		- Stopping criterion type                      ', lssolver%istop						
		PRINT*,'		- Stopping criterion tolerance                 ', lssolver%tol	
		PRINT*,'		- Restart:                                     ', lssolver%rest	
		ENDIF
  
		PRINT*,'	'
  end if
END SUBROUTINE READ_input

