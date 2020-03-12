!*****************************************
! project: MHDG
! file: hdg_BC.f90
! date: 25/01/2017
! Set boundary conditions in the global 
! problem
!*****************************************

SUBROUTINE HDG_BC()
  USE globals
  USE LinearAlgebra
  USE physics
  USE analytical
  USE printUtils
  USE in_out
  USE MPI_OMP
  
  IMPLICIT NONE
#ifdef TOR3D   
!***********************************************************************
! 
!                            VERSION 3D TOROIDAL
! 
!***********************************************************************  
  integer                   :: itor,iel3,ifa,ifl,iel,fl,bc,Neq,Npel,i,j,Ndim,Fi,N2D,Ngfl
  integer*4                 :: Np2D,Np1Dpol,Np1Dtor,Npfl,ntorloc,itorg,Nfdir
  integer                   :: nod(refElPol%Nfacenodes)
  integer                   :: NGaussPol,NGaussTor,dd,delta
  integer*4                 :: ind_uf(refElTor%Nfl*phys%neq),ind_fe(refElTor%Nfl*phys%neq),ind_ff(refElTor%Nfl*phys%neq)  
  integer*4                 :: ind_fg(refElTor%Nfl*3*phys%neq)
  integer*4                 :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),ind(refElPol%Ngauss1d)
  real*8                    :: tel(refElTor%Nnodes1d)
  real*8                    :: htor,Xf(Mesh%Nnodesperface,Mesh%Ndim)
  real*8                    :: uf(refElTor%Nfl*phys%neq)
  real*8                    :: thetafg(refElTor%Ngauss1d),dsurf(refElTor%Ngl),n_g(refElTor%Ngl,3),t_g(refElPol%Ngauss1d,2)  
#ifdef TEMPERATURE
  integer*4                 :: ind_qf(refElTor%Nfl*3*phys%Neq)  
  real*8                    :: qqf(refElTor%Nfl,3*phys%Neq), qfg(refElTor%Ngl,phys%neq*3)
  real*8                    :: coefi,coefe  
#endif  
  real*8                    :: tdiv(numer%ntor+1)
  real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  real*8                    :: xyf(refElPol%Ngauss1d,2),xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
  integer                   :: ind_asf(refElTor%Nfl),ind_ash(refElTor%Nfl)

  Ndim        = 3                                               ! Number of dimensions
  Neq         = phys%neq
  N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Npel        = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Npfl        = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfdir       = Mesh%Ndir
  NGaussPol   = refElPol%Ngauss1D
  NGaussTor   = refElTor%Ngauss1D
  Ngfl        = refElTor%Ngl
  
  ! Some indices 
  ind_asf = (/(i,i = 0,Neq*(Npfl-1),Neq)/)  
  ind_ash = (/(i,i = 0,Neq*(Npfl-1)*Ndim,Neq*Ndim)/)  
    
  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i=1,numer%ntor
     tdiv(i+1) = i*htor
  END DO
    
#ifdef TEMPERATURE
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
#endif  

 
#ifdef PARALL
    IF (MPIvar%ntor.gt.1) THEN
       ntorloc = numer%ntor/MPIvar%ntor+1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
 
	   
  
  ! Initialize matrices
  elMat%Aul_dir = 0.  
  elMat%Aql_dir = 0. 
     
  !************************
  ! Loop in external faces
  !************************
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif  
     tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg+1)-tdiv(itorg))
     DO ifa = 1,Mesh%Nextfaces
        
        ! Global numbering of this face
        Fi = ifa + Mesh%Nintfaces 
     
        ! Element to which this face belongs
        iel = Mesh%extfaces(ifa,1)
        iel3 = (itor-1)*N2d+iel
        
        ! Face local numbering
        ifl = Mesh%extfaces(ifa,2)     
            
        ! Nodes in local numbering
        nod = refElPol%Face_nodes(ifl,:)

        ! Coordinates of the nodes of the face
        Xf  = Mesh%X(Mesh%T(iel,nod),:)
        
        ! Gauss points position
        xyf = matmul(refElPol%N1D,Xf)    
        thetafg = matmul(refElTor%N1D,tel)
        
        ! Exterior normal & elemental surface
        xyd_g     = matmul(refelPol%Nxi1D,Xf)
        xydNorm_g = sqrt(xyd_g(:,1)**2+xyd_g(:,2)**2)
        dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))							                
        t_g(:,1) = xyd_g(:,1)/xydNorm_g
        t_g(:,2) = xyd_g(:,2)/xydNorm_g
        n_g = 0.
        DO i=1,NGaussTor
           ind = (i-1)*NGaussPol + (/(j,j=1,NGaussPol)/) 
           n_g(ind,1) =  t_g(:,2)
           n_g(ind,2) = -t_g(:,1)
        END DO
                
						  ! Assembly indices  
        dd = 1+(itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Npfl)*Neq
        delta = dd + (N2D*Np2D+(Fi-1)*Npfl)*Neq
        ind_uf = delta+(/(i,i=0,Npfl*Neq-1)/)
        ind_ff = Np2d*Neq + (ifl-1)*Npfl*Neq+ (/(i,i=1,Npfl*Neq)/)        
        ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifl,:)-1)*Neq)) 
        ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifl,:)-1)*3*Neq))           
#ifdef TEMPERATURE
        ind_qf = (iel3-1)*Ndim*Neq*Npel + ind_fg
        qqf    = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Npfl /)))
#endif     
               
        ! Solution in this face
        uf = sol%u_tilde(ind_uf)

					   ! Type of boundary
					   fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL						
						   IF (fl.eq.0) CYCLE
#endif  					
					   bc = phys%bcflags(fl) 


					   SELECT CASE(bc)

					   CASE(bc_dirichlet)
					      CALL set_dirichletstr_bc()					   				
					   CASE(bc_dirichlet_weak_form)
					      CALL set_dirichletwf_bc()
!        CASE(bc_dirichlet_weak_form_oldvalues)
!           CALL set_dirichletwfoldval_bc()
!					   CASE(bc_NeumannH)
!					      CALL set_NeumannH()	
        CASE(bc_Bohm)           
		         CALL set_Bohm_bc()
					   CASE DEFAULT
					      WRITE(6,*) "Error: wrong boundary type"
					      STOP
					   END SELECT

     END DO
   END DO
  CONTAINS
  
  !****************************
  ! Dirichlet in weak form
  !****************************  
  SUBROUTINE set_dirichletwf_bc
       integer                   :: g,igpol,igtor
       real*8                    :: uex(refElTor%Ngl,Neq)       
							real*8                    :: f(Neq,Npfl),dsurfg
							real*8,pointer            :: Nfi(:,:)
							real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
							
       ! Set the face solution           
       CALL analytical_solution(xyf(:,1),xyf(:,2),thetafg,uex)
       							                     
       ! Face shape functions
       Nfi => refElTor%sFTF	
       
				   !*****************************
				   ! Loop in face Gauss points
				   !*****************************       
										DO igtor = 1,NGaussTor
										   DO igpol = 1,NGaussPol
										   										      
										      g = (igtor-1)*NGaussPol+igpol 
						
											    ! Calculate the integration weight
													  IF (switch%axisym) THEN										
													     dsurfg = dsurf(g)*xyf(igpol,1)
													  ELSE 
													     dsurfg = dsurf(g)
													  END IF				
               NiNi = dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
               Ni = Nfi(g,:)*dsurfg													  
               
 													 ! Assembly contribution 
	 												 CALL assembly_dirichletwf_bc(iel3,ind_asf,ind_ff,uex(g,:),NiNi,Ni)

										 END DO
				    END DO
							
    
  END SUBROUTINE set_dirichletwf_bc
  
  

  !****************************
  ! Dirichlet in strong form
  !****************************  
  SUBROUTINE set_dirichletstr_bc
       integer                   :: g,igpol,igtor
       real*8                    :: uex(refElTor%Ngl,neq),dsurfg
       real*8                    :: b_f(refElTor%Ngl,3)   
       real*8,pointer            :: Nfi(:,:)
       real*8                    :: Ni(Npfl)
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif		       
							
							! Face shape functions
       Nfi => refElTor%sFTF	
                            
       ! Set the face solution           
       CALL analytical_solution(xyf(:,1),xyf(:,2),thetafg,uex) 
		
#ifndef MOVINGEQUILIBRIUM      
       ! Compute the magnetic field
			    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			        b_f = matmul(Nfi,phys%b(Mesh%T(iel,nod),:))
			    ELSE
										CALL defineMagneticField(xyf(:,1),xyf(:,2),thetafg,b_f)  
							END IF	
#else
											b_nod = 0.
										 Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2)
										 b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod
										 b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod
										 b_nod(:,3) = phys%Bt(Mesh%T(iel,:))/Bmod
										 b_f = matmul(Nfi,b_nod(nod,:))			       
#endif  

#ifdef TEMPERATURE
         ! Gradient solution at face gauss points
         qfg    = matmul(Nfi,qqf)
#endif

     !*****************************
     ! Loop in face Gauss points
     !*****************************       
					DO igtor = 1,NGaussTor
					   DO igpol = 1,NGaussPol
					   										      
					      g = (igtor-1)*NGaussPol+igpol 
			
				      ! Calculate the integration weight
						    IF (switch%axisym) THEN										
						       dsurfg = dsurf(g)*xyf(igpol,1)
						    ELSE
						       dsurfg = dsurf(g)
						    END IF		

          ! Shaspe functions product
          Ni = Nfi(g,:)*dsurfg
#ifdef TEMPERATURE
          CALL assembly_dirichletstr_bc(iel3,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),b_f(g,:),n_g(g,:))
#else
          CALL assembly_dirichletstr_bc(iel3,ind_asf,ind_ash,ind_fe,ind_fg,Ni,uex(g,:),b_f(g,:),n_g(g,:))
#endif
							   
							END DO			
				END DO				
  END SUBROUTINE set_dirichletstr_bc    
  
  
  


  !****************************
  ! Bohm
  !****************************    
  SUBROUTINE set_Bohm_bc()
       integer                   :: i,j,k
       integer                   :: g,igpol,igtor
       real*8                    :: dsurfg   
       real*8                    :: setval,delta   
							real*8                    :: ug(Ngfl),ufg(Ngfl,Neq)
							real*8                    :: upg(Ngfl,phys%npv)							
							real*8                    :: inc,sn
							real*8,pointer            :: Nfi(:,:)
       real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
							real*8                    :: soundSpeed
							real*8                    :: b_f(Ngfl,3)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,3),Bmod(Mesh%Nnodesperelem)
#endif				
#ifdef TEMPERATURE
       real*8                    :: qfg(Ngfl,neq*3)
#endif
       real*8                    :: bn
       logical                   :: ntang
       real                      :: tau_stab(Neq,Neq)
       real*8                    :: diff_iso_fac(phys%neq,phys%neq,Ngfl),diff_ani_fac(phys%neq,phys%neq,Ngfl) 
       
							! Face shape functions
       Nfi => refElTor%sFTF		
       							
#ifdef TEMPERATURE							         
       ! Gradient solution at face gauss points
       qfg    = matmul(Nfi,qqf)
#endif
										
       ! Set the face solution           
       ufg = matmul(Nfi,transpose(reshape(uf,[neq,Npfl])))		
              
       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)

       ! Compute diffusion at faces Gauss points
					  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

       ! Compute magnetic field at Gauss points
#ifndef MOVINGEQUILIBRIUM	       
			    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,nod),:))
			    ELSE
										CALL defineMagneticField(xyf(:,1),xyf(:,2),thetafg,b_f)  
							END IF				       
#else
										 b_nod = 0.
										 Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2)
										 b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod
										 b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod
										 b_nod(:,3) = phys%Bt(Mesh%T(iel,:))/Bmod
										 b_f = matmul(refElPol%N1D,b_nod(nod,:))
#endif       

     !*****************************
     ! Loop in face Gauss points
     !*****************************  
							DO igtor = 1,NGaussTor
							   DO igpol = 1,NGaussPol
							   										      
							      g = (igtor-1)*NGaussPol+igpol 
						
										   ! Calculate the integration weight
												 IF (switch%axisym) THEN										
												    dsurfg = dsurf(g)*xyf(igpol,1)
												 ELSE
												    dsurfg = dsurf(g)
												 END IF	
												 
												 bn = dot_product(b_f(g,:),n_g(g,:))
							   
							      ! Sound speed
#ifndef TEMPERATURE
             SoundSpeed = sqrt(phys%a)
#else
						       SoundSpeed = upg(g,9)
						       if (switch%decoup) then 
						          SoundSpeed = sqrt(phys%Mref)
						       endif
#endif
													! tangency
													ntang = .TRUE.							   
													inc = bn/norm2(b_f(g,1:2)) 
													
													setval = ufg(g,2)
													
													IF (abs(inc) .le. phys%bohmth ) THEN
																	sn = 1.
																	delta = 0
#ifdef TEMPERATURE
                 ntang = .FALSE.
#endif													 
													ELSE
																	sn = sign(1.,inc)
                 delta = 1
                 IF (( abs(upg(g,2)) ).le. SoundSpeed) THEN         
                     setval = sn*SoundSpeed
   !                 delta = 1
                 ELSE IF ((    abs(upg(g,2))  ).gt. SoundSpeed) THEN  
   !                 delta = 0
                     setval = upg(g,2)
                 END IF
													END IF						   									   
          
												! Impose the velocity at the sound speed in case of subsonic 
						      ! velocity. 
													NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
													Ni   = Nfi(g,:)*dsurfg
								  
								  !****************** Stabilization part******************************
            IF (numer%stab>1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab<5) THEN
                  CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b_f(g,:),n_g(g,:),iel,ifa,1.,xyf(g,:),tau_stab)
               ELSE
                  CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),1.,iel,tau_stab)
               ENDIF
          !****************** End stabilization part******************************     
            ELSE
               tau_stab = 1.
            END IF	

           ! Assembly Bohm contribution
#ifdef TEMPERATURE
						      CALL assembly_bohm_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),ufg(g,:),&
                 &b_f(g,:),n_g(g,:),tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
#else
            CALL assembly_bohm_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,ufg(g,:),&
                 &b_f(g,:),n_g(g,:),tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
#endif
          END DO
       END DO
       
  END SUBROUTINE set_Bohm_bc
  
  

  
  
  
  
    
  
  
#else
!***********************************************************************
! 
!                            VERSION 2D 
! 
!***********************************************************************
  integer                   :: ifa,ifl,iel,fl,bc,Npfl,Neq,Npel,i,Ndim,Fi,Ng1d
  integer                   :: nod(refElPol%Nfacenodes)
  integer                   :: ind_uf(refElPol%Nfacenodes*phys%neq),ind_ff(refElPol%Nfacenodes*phys%neq)
  integer                   :: ind_fe(refElPol%Nfacenodes*phys%neq)
  real*8                    :: uf(refElPol%Nfacenodes*phys%neq)
  real*8                    :: Xf(Mesh%Nnodesperface,Mesh%Ndim)
  integer                   :: ind_fG(refElPol%Nfacenodes*phys%neq*Mesh%Ndim)
#ifdef TEMPERATURE
  integer                   :: ind_qf(refElPol%Nfacenodes*phys%neq*Mesh%Ndim)
  real*8                    :: qqf(refElPol%Nfacenodes,phys%neq*Mesh%Ndim), qfg(refElPol%Ngauss1d,phys%neq*Mesh%Ndim)
  real*8                    :: coefi,coefe  
#endif  

  real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  logical                   :: save_tau
  real*8,allocatable        :: tau_save(:,:)
  real*8,allocatable        :: xy_g_save(:,:)
  integer                   :: indtausave(refElPol%Ngauss1d)
  integer                   :: ind_asf(refElPol%Nfacenodes),ind_ash(refElPol%Nfacenodes)
  
  
  save_tau = switch%saveTau
  Ndim = 2
  Npel = refElPol%Nnodes2D
  Npfl = refElPol%Nfacenodes
  Neq  = phys%neq
  Ng1d = refElPol%Ngauss1d
  
#ifdef TEMPERATURE
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
#endif  


  ! Some indices 
  ind_asf = (/(i,i = 0,Neq*(Npfl-1),Neq)/)  
  ind_ash = (/(i,i = 0,Neq*(Npfl-1)*Ndim,Neq*Ndim)/)  
  
  if (save_tau) then
     allocate(tau_save(Mesh%Nextfaces*refElPol%Ngauss1d,phys%neq))
     allocate(xy_g_save(Mesh%Nextfaces*refElPol%Ngauss1d,2))
     tau_save = 0.
     xy_g_save = 0.
  endif 
    
    
  ! Initialize matrices
  elMat%Aul_dir = 0.  
  elMat%Aql_dir = 0.  
  !************************
  ! Loop in external faces
  !************************
  DO ifa = 1,Mesh%Nextfaces
     
     ! Global numbering of this face
     Fi = ifa + Mesh%Nintfaces 
  
     ! Element to which this face belongs
     iel = Mesh%extfaces(ifa,1)
     
     ! Face local numbering
     ifl = Mesh%extfaces(ifa,2)     
         
     ! Nodes in local numbering
     nod = refElPol%Face_nodes(Mesh%extfaces(ifa,2),:)

     ! Coordinates of the nodes of the face
     Xf = Mesh%X(Mesh%T(iel,nod),:)
          
     ! Assembly indices
     ind_uf = (Fi-1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
     ind_ff = (ifl-1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
     ind_fe = reshape(tensorSumInt((/(i,i=1,neq)/),neq*(nod-1) ) , (/neq*Npfl/))
     ind_fG = reshape(tensorSumInt((/(i,i=1,neq*Ndim)/),neq*Ndim*(refElPol%face_nodes(ifl,:)-1)),(/neq*Ndim*Npfl/))
#ifdef TEMPERATURE
     ind_qf = (iel-1)*Ndim*Neq*Npel + ind_fG
     qqf    = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Npfl /)))
#endif     
     
     ! Solution in this face
     uf = sol%u_tilde(ind_uf)

					! Type of boundary
					fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL						
						IF (fl.eq.0) CYCLE
#endif  					
					bc = phys%bcflags(fl) 

					SELECT CASE(bc)

					CASE(bc_dirichlet)
					   CALL set_dirichletstr_bc()					   				
					CASE(bc_dirichlet_weak_form)
					   CALL set_dirichletwf_bc()
!					CASE(bc_NeumannH)
!					   CALL set_NeumannH()	
     CASE(bc_Bohm)
        
					   CALL set_Bohm_bc(tau_save_el,xy_g_save_el)
						  if (save_tau) then
						     indtausave = (ifa-1)*refElPol%Ngauss1d + (/(i,i =1,refElPol%Ngauss1d)/)
						     tau_save(indtausave,:)= tau_save_el
						     xy_g_save(indtausave,:)= xy_g_save_el
						  endif

					CASE DEFAULT
					   WRITE(6,*) "Error: wrong boundary type"
					   STOP
					END SELECT
					
  END DO

   if (save_tau) then
      write(6,*) "Saving tau in the boundary faces"
      call saveMatrix(tau_save,'tau_save_bound')
      call saveMatrix(xy_g_save,'xy_g_save_bound')
      deallocate(tau_save,xy_g_save)
      write(6,*) "Done saving tau!"
   endif
     
  CONTAINS
  
  !****************************
  ! Dirichlet in weak form
  !****************************  
  SUBROUTINE set_dirichletwf_bc
       integer                   :: g
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: uex(Ng1d,neq)       
       real*8                    :: xyf(Ng1d,Ndim),xyDer(Ng1d,Ndim)
							real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)

       ! Set the face solution           
       xyf    = matmul( refElPol%N1D,Xf )
       CALL analytical_solution(xyf(:,1),xyf(:,2),uex) 
       							                     
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)
       
							! Loop in 1D Gauss points
							DO g = 1,Ng1d
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
          NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
          Ni = refElPol%N1D(g,:)*dline
          CALL assembly_dirichletwf_bc(iel,ind_asf,ind_ff,uex(g,:),NiNi,Ni)
							END DO

  END SUBROUTINE set_dirichletwf_bc




  !****************************
  ! Dirichlet in strong form
  !****************************  
  SUBROUTINE set_dirichletstr_bc
       integer                   :: g,i,j,idm
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(Ng1d,Ndim),xyDer(Ng1d,Ndim)
       real*8                    :: uex(Ng1d,neq)
							real*8                    :: t_g(Ndim),n_g(Ndim),fh(Npfl*neq)
       real*8                    :: b_f(Ng1d,Ndim),divb_f(Ng1d),drift_f(Ng1d,Ndim)
       real*8                    :: Ni(Npfl)
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,Ndim),Bmod(Mesh%Nnodesperelem)
#endif		     
       ! Set the face solution           
       xyf    = matmul( refElPol%N1D,Xf )
       CALL analytical_solution(xyf(:,1),xyf(:,2),uex) 
		
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

#ifndef MOVINGEQUILIBRIUM      
       ! Compute the magnetic field
			    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,nod),:))
			        divb_f = matmul(refElPol%N1D,phys%divb(Mesh%T(iel,nod)))
			        drift_f = matmul(refElPol%N1D,phys%drift(Mesh%T(iel,nod),:))
			    ELSE
										CALL defineMagneticField(xyf(:,1),xyf(:,2),b_f,divb_f,drift_f)  
							END IF	
#else
											b_nod = 0.
										 Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2)
										 b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod
										 b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod
										 b_f = matmul(refElPol%N1D,b_nod(nod,:))			       
#endif  
#ifdef TEMPERATURE
         ! Gradient solution at face gauss points
         qfg    = matmul(refElPol%N1D,qqf)
#endif
							! Loop in 1D Gauss points
							DO g = 1,Ng1d
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
							   ! Unit normal to the boundary
							   t_g = xyDer(g,:)/xyDerNorm_g
							   n_g = [t_g(2), -t_g(1)]

          ! Shaspe functions product
          Ni = refElPol%N1D(g,:)*dline
#ifdef TEMPERATURE
          CALL assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),b_f(g,:),n_g)
#else
          CALL assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,uex(g,:),b_f(g,:),n_g)
#endif
  
							END DO							
							
  END SUBROUTINE set_dirichletstr_bc
    


  !****************************
  ! Bohm
  !****************************    
  SUBROUTINE set_Bohm_bc(tau_save_el,xy_g_save_el)
       integer                   :: g,i,j,k,idm
       real*8                    :: dline,xyDerNorm_g
       real*8 :: setval
       real*8                    :: xyf(Ng1d,Ndim),xyDer(Ng1d,Ndim)
							real*8                    :: t_g(Ndim),n_g(Ndim)
							real*8                    :: ufg(Ng1d,Neq), upg(Ng1d,phys%npv)							
							real*8                    :: inc,sn,delta
							real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
							real*8                    :: soundSpeed
       real*8                    :: tau_ngamma(2)
							real*8                    :: b_f(Ng1d,Ndim),divb_f(Ng1d),drift_f(Ng1d,Ndim)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,Ndim),Bmod(Mesh%Nnodesperelem)
#endif				
#ifdef TEMPERATURE
       real*8                    :: qfg(Ng1d,neq*Ndim)
#endif
       real*8                    :: bn
       logical                   :: ntang
       real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)       
       real                      :: tau_stab(Neq,Neq)
       real*8                    :: diff_iso_fac(phys%neq,phys%neq,Ng1d),diff_ani_fac(phys%neq,phys%neq,Ng1d)       
       
       

#ifdef TEMPERATURE							         
       ! Gradient solution at face gauss points
       qfg    = matmul(refElPol%N1D,qqf)
#endif
										
							! Gauss points position                     
						 xyf    = matmul( refElPol%N1D,Xf )
						 
       ! Compute diffusion at faces Gauss points
					  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)
							  						 
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

       ! Set the face solution           
       ufg = matmul(refElPol%N1D,transpose(reshape(uf,[neq,Npfl])))		
              
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)

! VEDERE PERCHÈ CALCOLO DUE VOLTE LE VARIABILI FISICHE NEI PUNTI DI GAUSS!!!!!!!!!!!!!!!!!!
       ! Compute magnetic field at Gauss points
#ifndef MOVINGEQUILIBRIUM	       
			    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,nod),:))
			        divb_f = matmul(refElPol%N1D,phys%divb(Mesh%T(iel,nod)))
			        drift_f = matmul(refElPol%N1D,phys%drift(Mesh%T(iel,nod),:))
			    ELSE
										CALL defineMagneticField(xyf(:,1),xyf(:,2),b_f,divb_f,drift_f)  
							END IF				       
#else
										 b_nod = 0.
										 Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2)
										 b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod
										 b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod
										 b_f = matmul(refElPol%N1D,b_nod(nod,:))
#endif       
							! Loop in 1D Gauss points
							DO g = 1,Ng1d


							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
							   
							   ! Unit normal to the boundary
							   t_g = xyDer(g,:)/xyDerNorm_g
							   n_g = [t_g(2), -t_g(1)]
							   bn = dot_product(b_f(g,:),n_g)
							   ! Sound speed
#ifndef TEMPERATURE
          SoundSpeed = sqrt(phys%a)
#else
          SoundSpeed = upg(g,9)
          if (switch%decoup) then 
             SoundSpeed = sqrt(phys%Mref)
          endif
#endif
										! tangency
										ntang = .TRUE.							   
							   inc = bn/norm2(b_f(g,:)) 

#ifdef NGAMMA
          setval = ufg(g,2)

										IF (abs(inc) .le. phys%bohmth ) THEN
													 sn = 1.
													 delta = 0 ! TODO try to module delta in the WEST case
#ifdef TEMPERATURE
              ntang = .FALSE.
#endif													 
										ELSE
													 sn = sign(1.,inc)
              delta = 1
              IF (( abs(upg(g,2)) ).le. SoundSpeed) THEN         
                  setval = sn*SoundSpeed
!                 delta = 1
              ELSE IF ((    abs(upg(g,2))  ).gt. SoundSpeed) THEN  
!                 delta = 0
                  setval = upg(g,2)
              END IF
										END IF						   									   
#endif          
									! Impose the velocity at the sound speed in case of subsonic 
         ! velocity. 
										NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline						
										Ni   = refElPol%N1D(g,:)*dline						
								  !****************** Stabilization part******************************
            IF (numer%stab>1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab<5) THEN
                  CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b_f(g,:),n_g,iel,ifa,1.,xyf(g,:),tau_stab)
               ELSE
                  CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b_f(g,:),n_g,xyf(g,:),1.,iel,tau_stab)
               ENDIF

               if (save_tau) then
                  DO i=1,Neq
                     tau_save_el(g,:) = tau_stab(i,i)
                  END DO
                  xy_g_save_el(g,:) = xyf(g,:)
               endif                 
          !****************** End stabilization part******************************     
            ELSE
               tau_stab = 1.
            END IF			
           
           ! Assembly Bohm contribution
#ifdef TEMPERATURE
						      CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
                   &ufg(g,:),b_f(g,:),n_g,tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
#else
            CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,&
                   &ufg(g,:),b_f(g,:),n_g,tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
#endif
							END DO ! END loop in Gauss points
							
          
       
  END SUBROUTINE set_Bohm_bc
  
#endif   

  
    
  

!*******************************************
!   AUXILIARY ROUTINES FOR 2D AND 3D
!*******************************************  				 


      !*********************************
						! Assembly Dirichlet in weak form
      !*********************************
      SUBROUTINE assembly_dirichletwf_bc(iel,ind_asf,ind_ff,uexg,NiNi,Ni)
      integer*4    :: iel,ind_asf(:),ind_ff(:)
      real*8       :: NiNi(:,:),Ni(:),uexg(:)
      real*8       :: kmult(Neq*Npfl)
      integer*4    :: ind(Npfl)
      
      kmult= col(tensorProduct(uexg(:),Ni))
      DO i=1,Neq
         ind = i+ind_asf
         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
      END DO       
      END SUBROUTINE assembly_dirichletwf_bc



      !*********************************
						! Assembly Dirichlet in strong form
      !*********************************
#ifdef TEMPERATURE
						SUBROUTINE assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qg,uexg,bg,ng)
#else
      SUBROUTINE assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,uexg,bg,ng)
#endif
      integer*4        :: iel,ind_asf(:),ind_ash(:),ind_fe(:),ind_fg(:)
      real*8           :: Ni(:),uexg(:),bg(:),ng(:)
#ifdef TEMPERATURE
      real*8           :: qg(:)
#endif
      real*8           :: bn,An(Neq,Neq)
      integer          :: i,j,idm
      real*8           :: kmultstab(Neq*Npfl),kmultconv(Neq*Npfl),kmultfseq(Ndim*Neq*Npfl)           
      real*8           :: uexn(Ndim*neq)
      integer*4        :: ind(Npfl),ind_1f(Npfl)
#ifdef TEMPERATURE
						real*8           :: Qpr(Ndim,Neq)
						real*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
						real*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)							
#endif  
       bn = dot_product(bg,ng)

				   CALL jacobianMatricesFace(uexg,bn,An) 
				    
				   ! Contribution of the current integration point to the elemental matrix
       kmultstab = col(tensorProduct(uexg,Ni))
       kmultconv = col(tensorProduct(matmul(An,uexg),Ni))
       uexn      = reshape(tensorProduct(ng,uexg),(/Ndim*neq/)) ! Ndim*neq x 1
       kmultfseq = col(tensorProduct(uexn,Ni))  ! Ndim*neq*npf

       DO i=1,Neq
          ind = i+ind_asf
          elMat%Aul_dir(ind_fe(ind),iel) = elMat%Aul_dir(ind_fe(ind),iel) - numer%tau(i)*kmultstab(ind)+kmultconv(ind)
							      DO idm = 1,Ndim
							         ind_1f = ind_ash+idm+(i-1)*Ndim 
                elMat%Aql_dir(ind_fG(ind_1f),iel) = elMat%Aql_dir(ind_fG(ind_1f),iel) - kmultfseq(ind_1f)
             END DO
       END DO
#ifdef TEMPERATURE        
      ! Compute Q^T^(k-1)
      Qpr = reshape(qg,(/Ndim,Neq/))

      ! Compute V(U^(k-1)) 
      call computeVi(uexg,Vveci)
      call computeVe(uexg,Vvece)
      
      ! Compute dV_dU (k-1)
      call compute_dV_dUi(uexg,dV_dUi)
      call compute_dV_dUe(uexg,dV_dUe)
      
      ! Compute Gamma and tau
      gmi = dot_product(matmul(Qpr,Vveci),bg)   ! scalar    
      gme = dot_product(matmul(Qpr,Vvece),bg)  
      Taui = matmul(Qpr,dV_dUi)                       ! 2x3
      Taue = matmul(Qpr,dV_dUe)

			   ! Compute Alpha(U^(k-1))
			   Alphai = computeAlphai(uexg)  
			   Alphae = computeAlphae(uexg)       
			   
      ! Compute dAlpha/dU^(k-1)
      call compute_dAlpha_dUi(uexg,dAlpha_dUi)
      call compute_dAlpha_dUe(uexg,dAlpha_dUe)
      
      DO i=3,4
         DO j=1,4
            IF (i==3) THEN
               elMat%Aul_dir(ind_fe(i+ind_asf),iel) = elMat%Aul_dir(ind_fe(i+ind_asf),iel) - &
                             &coefi*Ni*bn*((gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),bg)))*uexg(j))
            ELSE
               elMat%Aul_dir(ind_fe(i+ind_asf),iel) = elMat%Aul_dir(ind_fe(i+ind_asf),iel) - &
                             &coefe*Ni*bn*((gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),bg)))*uexg(j))
            END IF
         END DO
      END DO
 
#ifdef VORTICITY 
      ! Non-linear term in the vorticity equation (\Grad// n/n b)
      ind_if = 3+ind_asf
      elMat%Aul_dir(ind_fe(ind_if),iel) = elMat%Aul_dir(ind_fe(ind_if),iel) - dot_product(Qpr(:,1),b )./uexg(1)**2*Ni*bn
#endif

      
#endif         
						END SUBROUTINE assembly_dirichletstr_bc





      !*********************************
						! Assembly Bohm
      !*********************************
#ifdef TEMPERATURE
						SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,bg,ng,tau,setval,delta,diffiso,diffani,ntang)
#else
      SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,ufg,bg,ng,tau,setval,delta,diffiso,diffani,ntang)
#endif
      integer*4        :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8           :: NiNi(:,:),Ni(:),ufg(:),bg(:),ng(:),tau(:,:),setval,delta
      real*8           :: diffiso(:,:),diffani(:,:)
      logical          :: ntang
#ifdef TEMPERATURE
      real*8           :: qfg(:)
#endif
      real*8           :: bn,Abohm(Neq,Neq)
      integer          :: i,j,k,idm
      integer*4        :: ind(Npfl),ind_1f(Npfl),ind_2f(Npfl)
#ifdef TEMPERATURE
						real*8           :: Qpr(Ndim,Neq)
						real*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
						real*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)							
#endif 


          bn = dot_product(bg,ng)
          
          ! Case diagonal matrix stabilization
          IF (numer%stab<5) THEN  			   
             DO i=1,Neq
                ind = i+ind_asf
                elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
                IF (i==2) THEN
                   elMat%Alu(ind_ff(ind),ind_fe(ind-1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind-1),iel) + tau(i,i)*(delta*setval )*NiNi
                   elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*(1-delta)*NiNi
                ELSE
                   elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
                END IF
             END DO    
          ELSE
          ! Case full matrix stabilization
             DO i=1,Neq
                ind_1f = i+ind_asf
                DO j=1,Neq
                   ind_2f = j+ind_asf
                   elMat%All(ind_ff(ind_1f),ind_ff(ind_2f),iel) = elMat%All(ind_ff(ind_1f),ind_ff(ind_2f),iel) - tau(i,j)*NiNi
                   IF (j==1) THEN
                      elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) = elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) +&
                                                                          &tau(i,j)*NiNi+tau(i,j+1)*(delta*setval )*NiNi
                   ELSE IF (j==2) THEN 
                      elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) = elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) + tau(i,j)*(1-delta)*NiNi
                   ELSE
                      elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) = elMat%Alu(ind_ff(ind_1f),ind_fe(ind_2f),iel) + tau(i,j)*NiNi
                   ENDIF 
                END DO
             END DO 
          ENDIF           			   

                        			       
#ifdef TEMPERATURE
          IF (ntang) THEN
                        
						      ! Compute Q^T^(k-1)
						      Qpr = reshape(qfg,(/Ndim,Neq/))

						      ! Compute V(U^(k-1)) 
						      call computeVi(ufg,Vveci)
            call computeVe(ufg,Vvece)
            
						      ! Compute dV_dU (k-1)
						      call compute_dV_dUi(ufg,dV_dUi)
						      call compute_dV_dUe(ufg,dV_dUe)

						      ! Compute Alpha(U^(k-1))
						      Alphai = computeAlphai(ufg)
						      Alphae = computeAlphae(ufg)

						      ! Compute dAlpha/dU^(k-1)
						      call compute_dAlpha_dUi(ufg,dAlpha_dUi)
						      call compute_dAlpha_dUe(ufg,dAlpha_dUe)
            
            ! Jacobian matrix for convection part
            CALL jacobianMatricesBohm(ufg,Abohm)

            gmi = dot_product(matmul(Qpr,Vveci),bg)  ! scalar
            gme = dot_product(matmul(Qpr,Vvece),bg)             ! scalar 
            Taui = matmul(Qpr,dV_dUi)                      ! 2x3
            Taue = matmul(Qpr,dV_dUe)      ! 2x3

             ! Parallel diffusion for temperature
             DO i=1,2
                ind_1f = ind_asf+i
                IF (i==1) THEN
                   DO j=1,Neq
                      ind_2f = ind_asf+j
                      elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) = elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) + Abohm(i,j)*NiNi*bn
                      elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) = elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) -&
                                                                        &coefi*(gmi*dAlpha_dUi(j)+Alphai*( dot_product(Taui(:,j),bg)  ))*NiNi*bn
                      DO k=1,Ndim
                         elMat%Alq(ind_ff(ind_1f+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(ind_1f+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                       &coefi*Alphai*Vveci(j)*bg(k)*NiNi*bn
                      END DO
                    END DO
                    elMat%fh(ind_ff(ind_1f+2),iel) = elMat%fh(ind_ff(ind_1f+2),iel) - coefi*Alphai*( dot_product (matmul(transpose(Taui),bg),ufg)  )*Ni*bn
                ELSE
                   DO j=1,Neq
                      ind_2f = ind_asf+j
                      elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) = elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) + Abohm(i,j)*NiNi*bn
                      elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) = elMat%All(ind_ff(ind_1f+2),ind_ff(ind_2f),iel) -&
                                                                        &coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),bg) ))*NiNi*bn
                      DO k=1,Ndim
                         elMat%Alq(ind_ff(ind_1f+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(ind_1f+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                       &coefe*Alphae*Vvece(j)*bg(k)*NiNi*bn
                      END DO
                   END DO
                   elMat%fh(ind_ff(ind_1f+2),iel) = elMat%fh(ind_ff(ind_1f+2),iel) - coefe*Alphae*( dot_product (matmul(transpose(Taue),bg),ufg)  )*Ni*bn 
                ENDIF
             END DO
									END IF ! tangency
#endif          

         ! Perpendicular diffusion 
						   DO k = 1,Neq
						      DO idm = 1,Ndim
						         ind_1f = ind_ash+idm+(k-1)*Ndim 
						         ind_2f = ind_asf+k
						         if (ntang) then
						            elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)=elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)-&
                                                               &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) ) 
						         else
						             elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)=elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)-NiNi*(ng(idm)-bn*bg(idm))
						         endif
						      END DO
									END DO

#ifdef VORTICITY
         ! Parallel diffusion 
         k = 4
			      DO idm = 1,Ndim
			         ind_1f = ind_ash+idm+(k-1)*Ndim 
			         ind_2f = ind_asf+k
			         elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)=elMat%Alq(ind_ff(ind_2f),ind_fG(ind_1f),iel)+NiNi*bn*bg(idm)
			      END DO
#endif


								END SUBROUTINE assembly_bohm_bc






END SUBROUTINE HDG_BC


 
  
