MODULE solve_psblas
  USE matrices_types
  USE types
  USE globals
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
#ifdef WITH_MLD2P4  
  use mld_prec_mod
#endif    
  IMPLICIT NONE
  
  
   TYPE PSBLAS_STRUC ! A type to store structures for the PSBLAS library
      type(psb_dspmat_type) :: mat          ! sparse matrix for PSBLAS
#ifndef WITH_MLD2P4      
      type(psb_dprec_type)  :: prec         ! preconditioner
#else      
      type(mld_dprec_type)  :: prec
#endif          
      type(psb_desc_type)   :: desc_a       ! descriptor
      type(psb_d_vect_type) :: x,b         ! dense vector
      integer(psb_ipk_)     :: ictxt,iam,np,info   ! parallel environment
   END TYPE
   
   TYPE(PSBLAS_STRUC) :: matPSBLAS
     
   CONTAINS
  
!		real*8  :: tps, tpe
!		integer :: cks,clock_rate,cke
!		if (lssolver%timing) then
!				call cpu_time(tps)					
!				call system_clock(cks,clock_rate)      
!		endif

!		if (lssolver%timing) then
!			call cpu_time(tpe)
!			call system_clock(cke,clock_rate)
!			matK%rtime1 = matK%rtime1+(cke-cks)/real(clock_rate)
!			matK%ctime1 = matK%ctime1+tpe-tps
!		end if

  
!***********************************************
! Initialization of matrix instance of the
! PSBLAS library
!***********************************************  
  SUBROUTINE init_mat_PSBLAS(matPSBLAS,matK)
   TYPE(MAT_CSR_TYP) :: matK
   TYPE(PSBLAS_STRUC) :: matPSBLAS  
   integer(psb_lpk_),allocatable   :: ind(:)
   integer(psb_lpk_),allocatable   :: rows(:),cols(:) 
   integer(psb_ipk_)               :: info 
   integer(psb_ipk_)               :: iam,np,ncoef
   integer                         :: irow,i 
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif

 !******************************************
 !Initialize parallel environment for psblas                                      
 !******************************************  
   call psb_init(matPSBLAS%ictxt)
   
  ! get information about the PSBLAS parallel environment
 ! psb_info: in{icontxt} out{iam,np}
 ! iam is the current process
 ! np number of processes  
   call psb_info(matPSBLAS%ictxt,iam,np)

   if (iam < 0) then 
     ! This should not happen, but just in case
     call psb_exit(matPSBLAS%ictxt)
     stop
   endif
 !  if(psb_get_errstatus() /= 0) goto 9999
 !  name='mld_s_pde3d'
   call psb_set_errverbosity(itwo)

 !******************************************
 !  Allocate a communication descriptor    
 !******************************************
   call psb_cdall(matPSBLAS%ictxt,matPSBLAS%desc_a,info,nl=matK%n)
   
 !******************************************
 !  Create the communication descriptor    
 !******************************************    
    do irow = 1,matK%n
       ncoef = matK%rowptr(irow+1)-matK%rowptr(irow)
       allocate(rows(ncoef),ind(ncoef),cols(ncoef))
       ind=(/(i,i=matK%rowptr(irow),matK%rowptr(irow+1)-1)/)
       rows = irow
       cols = matK%cols(ind)
       call psb_cdins(ncoef,rows,cols,matPSBLAS%desc_a,info)
       if (info .ne. psb_success_) write(6,*) "failed filling matrix!"
       deallocate(rows,ind,cols)
    end do  

 !******************************************
 !  Assembly the communication descriptor    
 !******************************************  
  call psb_cdasb(matPSBLAS%desc_a,info) 
  
 !******************************************
 !  Allocate matrix,unknowns and rhs                             
 !******************************************     
   if (info == psb_success_) call psb_spall(matPSBLAS%mat,matPSBLAS%desc_a,info,nnz=matK%nnz)
   if (info == psb_success_) call psb_geall(matPSBLAS%x,matPSBLAS%desc_a,info)
   if (info == psb_success_) call psb_geall(matPSBLAS%b,matPSBLAS%desc_a,info)
   
			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime1 = matK%rtime1+(cke-cks)/real(clock_rate)
				matK%ctime1 = matK%ctime1+tpe-tps
			end if   
  END SUBROUTINE init_mat_PSBLAS
  
  


!***********************************************
! Here the matrix is passed to the PSBLAS
! library
!***********************************************  
  SUBROUTINE build_mat_PSBLAS(matPSBLAS,matK)
   TYPE(MAT_CSR_TYP) :: matK
   TYPE(PSBLAS_STRUC) :: matPSBLAS   
   integer(psb_lpk_),allocatable   :: ind(:)
   integer(psb_lpk_),allocatable   :: rows(:),cols(:)   
   real(psb_dpk_), allocatable     :: vals(:)
   integer(psb_ipk_)               :: info,irow,i,ncoef
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif
   
!******************************************
!  Initialize matrix to zero                             
!******************************************     
   call psb_sprn(matPSBLAS%mat, matPSBLAS%desc_a, info)
   
!******************************************
!  Fill matrix                             
!******************************************     
    do irow = 1,matK%n
       ncoef = matK%rowptr(irow+1)-matK%rowptr(irow)
       allocate(rows(ncoef),ind(ncoef),cols(ncoef),vals(ncoef))
       ind=(/(i,i=matK%rowptr(irow),matK%rowptr(irow+1)-1)/)
       rows = irow
       cols = matK%cols(ind)
       vals = matK%vals(ind)
       call psb_spins(ncoef,rows,cols,vals,matPSBLAS%mat,matPSBLAS%desc_a,info)
       if (info .ne. psb_success_) write(6,*) "failed filling matrix!"
       deallocate(rows,ind,cols,vals)
    end do  

 !******************************************
 !  Assembly matrix                      
 !******************************************     
    call psb_spasb(matPSBLAS%mat,matPSBLAS%desc_a,info,dupl=psb_dupl_err_)
    if (info .ne. psb_success_) then
       write(6,*) "failed assemblying the matrix"
       stop
    endif   
    

			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime2 = matK%rtime2+(cke-cks)/real(clock_rate)
				matK%ctime2 = matK%ctime2+tpe-tps
			end if     
  END SUBROUTINE build_mat_PSBLAS
  




!***********************************************
! Build the preconditioner 
! 
!***********************************************  
#ifndef WITH_MLD2P4
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
   TYPE(PSBLAS_STRUC) :: matPSBLAS
   integer(psb_ipk_)  :: info
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif

     call matPSBLAS%prec%init(matPSBLAS%ictxt,lssolver%ptype,info)
     call matPSBLAS%prec%build(matPSBLAS%mat,matPSBLAS%desc_a,info)

    if (info .ne. psb_success_) then
       write(6,*) "Something wrong in the preconditioner construction"
       stop
    end if
			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime3 = matK%rtime3+(cke-cks)/real(clock_rate)
				matK%ctime3 = matK%ctime3+tpe-tps
			end if    
  END SUBROUTINE build_prec_PSBLAS
#else
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
   TYPE(PSBLAS_STRUC) :: matPSBLAS
   integer(psb_ipk_)  :: info
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif
			   
  call matPSBLAS%prec%init(matPSBLAS%ictxt,lssolver%ptype,info)
  select case(lssolver%ptype)
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','L1-JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
    
  case ('BJAC')
    call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
    call matPSBLAS%prec%set('sub_solve',       lssolver%solve,   info)
    call matPSBLAS%prec%set('sub_fillin',      lssolver%fill,    info)
    call matPSBLAS%prec%set('sub_iluthrs',     lssolver%thr,     info)

  case('AS')
    call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
    call matPSBLAS%prec%set('sub_ovr',         lssolver%novr,    info)
    call matPSBLAS%prec%set('sub_restr',       lssolver%restr,   info)
    call matPSBLAS%prec%set('sub_prol',        lssolver%prol,    info)
    call matPSBLAS%prec%set('sub_solve',       lssolver%solve,   info)
    call matPSBLAS%prec%set('sub_fillin',      lssolver%fill,    info)
    call matPSBLAS%prec%set('sub_iluthrs',     lssolver%thr,     info)
    
  case ('ML') 
    ! multilevel preconditioner

    call matPSBLAS%prec%set('ml_cycle',        lssolver%mlcycle,    info)
    call matPSBLAS%prec%set('outer_sweeps',    lssolver%outer_sweeps,info)
    if (lssolver%csize>0)&
         & call matPSBLAS%prec%set('min_coarse_size', lssolver%csize,      info)
    if (lssolver%mncrratio>1)&
         & call matPSBLAS%prec%set('min_cr_ratio',   lssolver%mncrratio, info)
    if (lssolver%maxlevs>0)&
         & call matPSBLAS%prec%set('max_levs',    lssolver%maxlevs,    info)
    if (lssolver%athres >= dzero) &
         & call matPSBLAS%prec%set('aggr_thresh',     lssolver%athres,  info)
!    if (lssolver%thrvsz>0) then
!      do k=1,min(lssolver%thrvsz,size(prec%precv)-1)
!        call matPSBLAS%prec%set('aggr_thresh',     lssolver%athresv(k),  info,ilev=(k+1))
!      end do
!    end if

    call matPSBLAS%prec%set('aggr_prol',       lssolver%aggr_prol,   info)
    call matPSBLAS%prec%set('par_aggr_alg',    lssolver%par_aggr_alg,   info)
    call matPSBLAS%prec%set('aggr_ord',        lssolver%aggr_ord,   info)
    call matPSBLAS%prec%set('aggr_filter',     lssolver%aggr_filter,info)


    call matPSBLAS%prec%set('smoother_type',   lssolver%smther,     info)
    call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps,    info)

    select case (psb_toupper(lssolver%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
      ! do nothing
    case default
      call matPSBLAS%prec%set('sub_ovr',         lssolver%novr,       info)
      call matPSBLAS%prec%set('sub_restr',       lssolver%restr,      info)
      call matPSBLAS%prec%set('sub_prol',        lssolver%prol,       info)
      call matPSBLAS%prec%set('sub_solve',       lssolver%solve,      info)
      call matPSBLAS%prec%set('sub_fillin',      lssolver%fill,       info)
      call matPSBLAS%prec%set('sub_iluthrs',     lssolver%thr,        info)
    end select

    if (psb_toupper(lssolver%smther2) /= 'NONE') then
      call matPSBLAS%prec%set('smoother_type',   lssolver%smther2,   info,pos='post')
      call matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps2,  info,pos='post')
      select case (psb_toupper(lssolver%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
        ! do nothing
      case default
        call matPSBLAS%prec%set('sub_ovr',         lssolver%novr2,     info,pos='post')
        call matPSBLAS%prec%set('sub_restr',       lssolver%restr2,    info,pos='post')
        call matPSBLAS%prec%set('sub_prol',        lssolver%prol2,     info,pos='post')
        call matPSBLAS%prec%set('sub_solve',       lssolver%solve2,    info,pos='post')
        call matPSBLAS%prec%set('sub_fillin',      lssolver%fill2,     info,pos='post')
        call matPSBLAS%prec%set('sub_iluthrs',     lssolver%thr2,      info,pos='post')
      end select
    end if

    call matPSBLAS%prec%set('coarse_solve',    lssolver%csolve,    info)
    if (psb_toupper(lssolver%csolve) == 'BJAC') &
         &  call matPSBLAS%prec%set('coarse_subsolve', lssolver%csbsolve,  info)
    call matPSBLAS%prec%set('coarse_mat',      lssolver%cmat,      info)
    call matPSBLAS%prec%set('coarse_fillin',   lssolver%cfill,     info)
    call matPSBLAS%prec%set('coarse_iluthrs',  lssolver%cthres,    info)
    call matPSBLAS%prec%set('coarse_sweeps',   lssolver%cjswp,     info)

  end select
  
  ! build the preconditioner
  call matPSBLAS%prec%hierarchy_build(matPSBLAS%mat,matPSBLAS%desc_a,info)
  call matPSBLAS%prec%smoothers_build(matPSBLAS%mat,matPSBLAS%desc_a,info)   
   
    if (info .ne. psb_success_) then
       write(6,*) "Something wrong in the preconditioner construction"
       stop
    end if
			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime3 = matK%rtime3+(cke-cks)/real(clock_rate)
				matK%ctime3 = matK%ctime3+tpe-tps
			end if     
  END SUBROUTINE build_prec_PSBLAS
#endif    



!***********************************************
! Fill RHS and initial guess
! 
!***********************************************  
  SUBROUTINE fill_vec_PSBLAS(matPSBLAS)
   TYPE(PSBLAS_STRUC) :: matPSBLAS
   integer(psb_lpk_),allocatable   :: ind(:)
   integer(psb_ipk_)  :: info,n,i
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif

     call matPSBLAS%b%set(rhs%vals)
     call matPSBLAS%x%set(sol%u_tilde)
     
   
!     n = matPSBLAS%mat%get_nrows() 
!     call matPSBLAS%x%zero()
!     call matPSBLAS%b%zero() 
!     allocate(ind(n))
!     ind = (/(i,i=1,n)/)
!     call psb_geins(n,ind,rhs%vals,matPSBLAS%b,matPSBLAS%desc_a,info)
!     call psb_geins(n,ind,sol%u_tilde,matPSBLAS%x,matPSBLAS%desc_a,info)
!     deallocate(ind)  

 !******************************************
 !  Assembly rhs and initial guess                            
 !******************************************   
 
     ! Init. guess
     call psb_geasb(matPSBLAS%x,matPSBLAS%desc_a,info)
		   if (info .ne. psb_success_) then
		      write(6,*) "Something wrong in the initial guess assembly"
		      stop
		   end if
		   
     ! RHS
     call psb_geasb(matPSBLAS%b,matPSBLAS%desc_a,info)
		   if (info .ne. psb_success_) then
		      write(6,*) "Something wrong in the rhs assembly"
		      stop
		   end if
		   
			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime4 = matK%rtime4+(cke-cks)/real(clock_rate)
				matK%ctime4 = matK%ctime4+tpe-tps
			end if  		   
  END SUBROUTINE fill_vec_PSBLAS     
  
  
  
  
!***********************************************
! Solve call for PSBLAS, assign output to 
! u_tilde
!***********************************************   
  SUBROUTINE solve_mat_PSBLAS(matPSBLAS)
   TYPE(PSBLAS_STRUC) :: matPSBLAS  
   integer(psb_ipk_)     :: info   ! parallel environment
   integer(psb_ipk_)     :: iter,itmax,itrace,istopc,irst
   real(psb_dpk_)        :: err
			real*8  :: tps, tpe
			integer :: cks,clock_rate,cke
			if (lssolver%timing) then
					call cpu_time(tps)					
					call system_clock(cks,clock_rate)      
			endif

 !******************************************
 !  Solve                            
 !****************************************** 
   call psb_krylov(lssolver%kmethd,matPSBLAS%mat,matPSBLAS%prec,matPSBLAS%b,matPSBLAS%x,lssolver%tol,matPSBLAS%desc_a,info,& 
        & itmax=lssolver%itmax,iter=iter,err=err,itrace=lssolver%itrace,istop=lssolver%istop,irst=lssolver%rest) 

 !******************************************
 !  Store solution                            
 !******************************************    
!   sol%u_tilde = (1.-numer%dumpnr)*sol%u_tilde + numer%dumpnr*matPSBLAS%x%get_vect()     
			if (lssolver%timing) then
				call cpu_time(tpe)
				call system_clock(cke,clock_rate)
				matK%rtime5 = matK%rtime5+(cke-cks)/real(clock_rate)
				matK%ctime5 = matK%ctime5+tpe-tps
			end if     
  END SUBROUTINE solve_mat_PSBLAS
  
  
  
!******************************************
!  Free memory                    
!******************************************    
  SUBROUTINE terminate_PSBLAS()
  integer(psb_ipk_)     :: info
  
  call psb_gefree(matPSBLAS%b,matPSBLAS%desc_a,info)
  call psb_gefree(matPSBLAS%x,matPSBLAS%desc_a,info)
  call psb_spfree(matPSBLAS%mat,matPSBLAS%desc_a,info)
  call matPSBLAS%prec%free(info)
  call psb_cdfree(matPSBLAS%desc_a,info)  
  call psb_exit(matPSBLAS%ictxt)
  END SUBROUTINE terminate_PSBLAS
  
END MODULE solve_psblas
