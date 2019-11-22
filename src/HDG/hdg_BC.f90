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
  integer                   :: itor,iel3,ifa,ifl,iel,fl,bc,Nf,Nfp,Neq,Np,i,j,Ndim,Fi,N2D
  integer*4                 :: Np2D,Np1Dpol,Np1Dtor,Nfl,ntorloc,itorg,Nfdir
  integer                   :: nod(refElPol%Nfacenodes)
  integer                   :: NGaussPol,NGaussTor,dd,delta
  integer*4                 :: ind_ff(refElTor%Nfl*phys%neq),ind_fe(refElTor%Nfl*phys%neq),ind_fl(refElTor%Nfl*phys%neq)  
  integer*4                 :: ind_fg(refElTor%Nfl*3*phys%neq),ind_sta(refElPol%Nfaces+2)
  integer*4                 :: ind_dim(refElPol%Nfaces+2)
  integer*4                 :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),ind(refElPol%Ngauss1d)
  real*8                    :: tel(refElTor%Nnodes1d)
  real*8                    :: htor,Xf(Mesh%Nnodesperface,Mesh%Ndim)
  real*8                    :: uf(refElTor%Nfl*phys%neq)
  real*8                    :: thetafg(refElTor%Ngauss1d),dsurf(refElTor%Ngl),n_g(refElTor%Ngl,3),t_g(refElPol%Ngauss1d,2)  
  real*8, allocatable       :: Df(:,:),Ef(:,:),fh(:),Edirf(:),Hdirf(:)
  real*8, allocatable       :: Lf(:,:),Qf(:,:),Cdirf(:)
#ifdef TEMPERATURE
  integer*4                 :: ind_qf(refElTor%Nfl*3*phys%Neq)  
  real*8                    :: qqf(refElTor%Nfl,3*phys%Neq), qfg(refElTor%Ngl,phys%neq*3)
  real*8, allocatable       :: TUhf(:,:),TQhf(:,:),Tfhf(:),Thdirf(:),Hf(:,:)
  real*8                    :: coefi,coefe  
#endif  
  real*8                    :: tdiv(numer%ntor+1)
  real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  real*8                    :: xyf(refElPol%Ngauss1d,2),xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
  logical                   :: save_tau
  real*8,allocatable        :: tau_save(:,:)
  real*8,allocatable        :: xy_g_save(:,:)
  integer                   :: indtausave(refElPol%Ngauss1d)

  Ndim        = 3                                               ! Number of dimensions
  Neq         = phys%neq
  N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np          = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Nfl         = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfp         = Np2D*2+refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
  Nf          = refElPol%Nfaces                                     ! Number of faces in the 2D mesh
  Nfdir       = Mesh%Ndir

  save_tau = switch%saveTau
  
  ind_loc = 0
  DO i = 1,Nf
     DO j = 1,Nfl
        ind_loc(i,j) = refElPol%Nnodes2D + Nfl*(i-1) + j
     END DO
  END DO
  
  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i=1,numer%ntor
     tdiv(i+1) = i*htor
  END DO
    
  ALLOCATE(Df(Neq*Nfl,Neq*Nfl))
  ALLOCATE(Ef(Neq*Nfl,Neq*Nfl))
  ALLOCATE(fh(Neq*Nfl))
  ALLOCATE(Edirf(Neq*Nfl)) 
  ALLOCATE(Hdirf(Neq*Nfl))
  ALLOCATE(Cdirf(Neq*Nfl*ndim)) 
  ALLOCATE(Lf(Neq*Nfl,Neq*Nfl*Ndim))
  ALLOCATE(Qf(Neq*Nfl,Neq*Nfl*Ndim))
#ifdef TEMPERATURE
  ALLOCATE(TUhf(Neq*Nfl,Neq*Nfl))
  ALLOCATE(TQhf(Neq*Nfl,Neq*Ndim*Nfl))
  ALLOCATE(Tfhf(Neq*Nfl))
  ALLOCATE(Thdirf(Neq*Nfl))
  ALLOCATE(Hf(Neq*Nfl,Neq*Nfl))
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
#endif  
!  
!  
!  write(6,*) "coef: ", coef
!  write(6,*) "coefe: ", coefe  
  if (save_tau) then
     allocate(tau_save(Mesh%Nextfaces*refElPol%Ngauss1d,phys%neq))
     allocate(xy_g_save(Mesh%Nextfaces*refElPol%Ngauss1d,2))
     tau_save = 0.
     xy_g_save = 0.
  endif 
 
#ifdef PARALL
    IF (MPIvar%ntor.gt.1) THEN
       ntorloc = numer%ntor/MPIvar%ntor+1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
 
  NGaussPol = refElPol%Ngauss1D
  NGaussTor = refElTor%Ngauss1D	    
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
        
        ! Initialize local matrices
        Df = 0.
        Ef = 0.
        Edirf = 0.d0
        Hdirf = 0.d0
        fh = 0.
        Cdirf = 0.
        Lf = 0.
        Qf = 0.
#ifdef TEMPERATURE
        TUhf = 0.
        TQhf = 0.
        Tfhf = 0.
        Thdirf = 0.
        Hf   = 0.
#endif      
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
        dd = 1+(itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Nfl)*Neq
        delta = dd + (N2D*Np2D+(Fi-1)*Nfl)*Neq
        ind_ff = delta+(/(i,i=0,Nfl*Neq-1)/)
        ind_fl = Np2d*Neq + (ifl-1)*Nfl*Neq+ (/(i,i=1,Nfl*Neq)/)        
        ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifl,:)-1)*Neq)) 
        ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifl,:)-1)*3*Neq))           
#ifdef TEMPERATURE
        ind_qf = (iel3-1)*Ndim*Neq*Np + ind_fg
        qqf    = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Nfl /)))
#endif     
        
!write(6,*) "itor: ",itor         
!write(6,*) "dd: ",dd
!write(6,*) "delta: ",delta
!call displayVectorInt(ind_ff)        
        ! Solution in this face
        uf = sol%u_tilde(ind_ff)

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

        CASE(bc_dirichlet_weak_form_oldvalues)
           CALL set_dirichletwfoldval_bc()
					   CASE(bc_NeumannH)
					      CALL set_NeumannH()	
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
!call displayVectorInt(ind_fl)					   
					   elMat%Df(ind_fl,ind_fe,iel3) = Df
					   elMat%Ef(ind_fl,ind_fl,iel3) = Ef 		
					   elMat%fh(ind_fl,iel3)        = fh
					   IF (Mesh%ndir.gt.0) THEN					
					      elMat%Edir(ind_fe,iel3) = Edirf
					      elMat%Hdir(ind_fe,iel3) = Hdirf
					   END IF					
       	elMat%Lf(ind_fl,ind_fG,iel3) = Lf
       	elMat%Qf(ind_fl,ind_fG,iel3) = Qf
					   IF (Mesh%ndir.gt.0) THEN					    	
       	   elMat%Cdir(ind_fG,iel3) = Cdirf
#ifdef TEMPERATURE     
           elMat%Thdir(ind_fe,iel3)= Thdirf
#endif    	   
       	END IF
#ifdef TEMPERATURE
        elMat%TQhf(ind_fl,ind_fG,iel3)  = TQhf
        elMat%Tfhf(ind_fl,iel3)         = Tfhf
        elMat%Hf(ind_fl,ind_fl,iel3)    = Hf-TUhf
#endif    	

     END DO
   END DO
   if (save_tau) then
      write(6,*) "Saving tau in the boundary faces"
      call saveMatrix(tau_save,'tau_save_bound')
      call saveMatrix(xy_g_save,'xy_g_save_bound')
      deallocate(tau_save,xy_g_save)
      write(6,*) "Done saving tau!"
   endif
     
  DEALLOCATE(Df,Ef,Edirf,Hdirf,fh)

  DEALLOCATE(Cdirf,Lf,Qf)
#ifdef TEMPERATURE
  DEALLOCATE(Thdirf,TQhf,TUhf,Tfhf,Hf) 
#endif   
  CONTAINS
  
  !****************************
  ! Dirichlet in weak form
  !****************************  
  SUBROUTINE set_dirichletwf_bc
       integer                   :: g,igpol,igtor
       real*8                    :: uex(refElTor%Ngl,neq)       
							real*8                    :: f(neq,Nfl),dsurfg
							real*8,dimension(Nfl,Nfl) :: Massf
							real*8,pointer            :: Nfi(:,:)
							
							
							Massf = 0.d0
							f = 0.

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
													  
													 ! Contribution of the current integration point to the elemental matrix
				          Massf = Massf + dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
				          f = f + dsurfg*tensorProduct(uex(g,:),Nfi(g,:))
										 END DO
				    END DO
								CALL expand_matrix_mass_stab(Massf,neq,numer%tau,Ef)
							
							! Multiply for stabilization before reshaping
							DO i=1,neq
							   f(i,:) = numer%tau(i)*f(i,:)
							END DO
       fH = -reshape(f,(/Nfl*neq/))
    
  END SUBROUTINE set_dirichletwf_bc




  !****************************
  ! Dirichlet in strong form
  !****************************  
  SUBROUTINE set_dirichletstr_bc
       integer                   :: g,igpol,igtor
       real*8                    :: uex(refElTor%Ngl,neq),dsurfg
							real*8                    :: fh(Nfl*neq)
							real*8,dimension(Nfl,neq) :: Edir,Hdir
       real*8                    :: b_f(refElTor%Ngl,3)   
       real*8                    :: Cdir(ndim*neq,Nfl)
       real*8                    :: uexn(ndim*neq),bn
       real*8,dimension(neq,neq) :: An
       real*8,pointer            :: Nfi(:,:)
#ifdef TEMPERATURE
							real*8                    :: Qpr(Ndim,Neq)
							real*8                    :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
							real*8                    :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)			
       integer                   :: ind_asf(Nfl)
#endif      
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif		       
							Edir = 0.d0
							Hdir = 0.d0
       Cdir = 0.
							
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
         ind_asf = (/ (i, i = 0, Neq*(Nfl-1),Neq) /)											
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

							   ! Jacobian matrices									   
							   bn = dot_product(b_f(g,:),n_g(g,:))
							   CALL jacobianMatricesFace(uex(g,:),bn,An) 
							    
							   ! Contribution of the current integration point to the elemental matrix
 						   Edir = Edir + dsurfg*tensorProduct(Nfi(g,:),uex(g,:))			 
          Hdir = Hdir + tensorProduct( Nfi(g,:),matmul(An,uex(g,:)) )*dsurfg  
          uexn = reshape(tensorProduct(n_g(g,:),uex(g,:)),(/ndim*neq/))
          Cdir = Cdir + dsurfg*tensorProduct(uexn,Nfi(g,:))        
#ifdef TEMPERATURE        
         
         ! Compute Q^T^(k-1)
         Qpr = reshape(qfg(g,:),(/Ndim,Neq/))

         ! Compute V(U^(k-1)) 
         call computeVi(uex(g,:),Vveci)
         call computeVe(uex(g,:),Vvece)
         
			      ! Compute dV_dU (k-1)
			      call compute_dV_dUi(uex(g,:),dV_dUi)
			      call compute_dV_dUe(uex(g,:),dV_dUe)
			      
         ! Compute Gamma and tau
		       gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))   ! scalar    
		       gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))  
		       Taui = matmul(Qpr,dV_dUi)                       ! 2x3
		       Taue = matmul(Qpr,dV_dUe)

						   ! Compute Alpha(U^(k-1))
						   Alphai = computeAlphai(uex(g,:))  
						   Alphae = computeAlphae(uex(g,:))       
						   
         ! Compute dAlpha/dU^(k-1)
         call compute_dAlpha_dUi(uex(g,:),dAlpha_dUi)
         call compute_dAlpha_dUe(uex(g,:),dAlpha_dUe)
         						   
         DO i=3,4
            DO j=1,4
               IF (i==3) THEN
                  Thdirf(i+ind_asf) = Thdirf(i+ind_asf) + coefi*Nfi(g,:)*dsurfg*bn*((gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),b_f(g,:))))*uex(g,j))
               ELSE
                  Thdirf(i+ind_asf) = Thdirf(i+ind_asf) + coefe*Nfi(g,:)*dsurfg*bn*((gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),b_f(g,:))))*uex(g,j))
               END IF
            END DO
         END DO								
         
#endif           
							END DO			
				END DO				
							
							! Elemental assembly
							! Multiply for stabilization before reshaping
							DO i=1,neq
							   Edir(:,i) = numer%tau(i)*Edir(:,i)
							END DO
       Edirf =  reshape(transpose(Edir),(/neq*Nfl/)) 
       Hdirf =  reshape(transpose(Hdir),(/neq*Nfl/)) 
       Cdirf =  reshape(Cdir,(/neq*Nfl*ndim/))
  END SUBROUTINE set_dirichletstr_bc
    



  
  !********************************************************
  ! Dirichlet in weak form using old values of the unknowns
  !********************************************************  
  SUBROUTINE set_dirichletwfoldval_bc
       integer                   :: g,igpol,igtor
       real*8                    :: ufg(refElTor%Ngl,Neq)
							real*8                    :: f(neq,Nfl),dsurfg
							real*8,dimension(Nfl,Nfl) :: Massf
							real*8,pointer            :: Nfi(:,:)
							
							
							Massf = 0.d0
							f = 0.
      							                     
       ! Face shape functions
       Nfi => refElTor%sFTF	

       ! Set the face solution           
       ufg = matmul(Nfi,transpose(reshape(uf,[neq,Nfl])))	
       
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
													  
													 ! Contribution of the current integration point to the elemental matrix
				          Massf = Massf + dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
				          f = f + dsurfg*tensorProduct(ufg(g,:),Nfi(g,:))
										 END DO
				    END DO
								CALL expand_matrix_mass_stab(Massf,neq,numer%tau,Ef)
												
							! Multiply for stabilization before reshaping
							DO i=1,neq
							   f(i,:) = numer%tau(i)*f(i,:)
							END DO
       fH = -reshape(f,(/Nfl*neq/))
    
  END SUBROUTINE set_dirichletwfoldval_bc





  !******************************************
  ! Dirichlet in weak form: velocity at 
  ! sound speed and threshold value for the
  ! rest of the variables
  !******************************************  
  SUBROUTINE set_dirichletwfbohm_bc
       integer                   :: g,igpol,igtor
       real*8                    :: dsurfg
							real*8                    :: f(neq,Nfl),soundspeed
							real*8                    :: bn,sn,delta,inc
							real*8                    :: b_f(refElTor%Ngl,3)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,3),Bmod(Mesh%Nnodesperelem)
#endif	
							real*8,dimension(Nfl,Nfl) :: Massf,zeros,BB21,BB22,NiNi
							real*8,dimension(Neq)     :: vals
							real*8                    :: ug(refElTor%Ngl),ufg(refElTor%Ngl,Neq),upgf(refElTor%Ngl,phys%npv)
							real*8                    :: upg(refElTor%Ngl,phys%npv)	
							real*8,pointer            :: Nfi(:,:)								
							
							Massf = 0.d0
							zeros = 0.
							BB21  = 0.
							BB22  = 0.
							f = 0.
							vals = 0.
							
							! Face shape functions
       Nfi => refElTor%sFTF								

       ! Set the face solution           
       ufg = matmul(Nfi,transpose(reshape(uf,[neq,Nfl])))		
              
       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)

#ifndef MOVINGEQUILIBRIUM	       
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
							   inc = bn/norm2(b_f(g,:))
										IF (abs(inc) .le. phys%bohmth ) THEN
													 sn = 1.
													 delta = 0
										ELSE
													 sn = sign(1.,inc)
               IF ((sn*upg(g,2)).lt.SoundSpeed) THEN         
                  delta = 1
               ELSE
                  delta = 0
               END IF
										END IF	               
!          vals(2) = sn*soundspeed*vals(1)               
               			      			   
									! Contribution of the current integration point to the elemental matrix
									NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
									Massf =  Massf + NiNi
					    BB21  =  BB21 + (delta*sn*SoundSpeed )*NiNi									
					    BB22  =  BB22 + (1-delta)*NiNi
         f = f + dsurfg*tensorProduct(vals,Nfi(g,:))
							END DO
				END DO
#ifdef TEMPERATURE
							CALL expand_matrix_BohmD_av(zeros,BB21,BB22,zeros,zeros,Df)
#else
							CALL expand_matrix_BohmD_av(zeros,BB21,BB22,Df)
#endif							
							CALL expand_matrix_BohmE_av(Massf,Ef)
       fH = -reshape(f,(/Nfl*neq/))
    
  END SUBROUTINE set_dirichletwfbohm_bc
  
  
      

  !****************************
  ! Neumann homogeneus
  !****************************    
  SUBROUTINE set_NeumannH
       integer                   :: g,igpol,igtor
       real*8                    :: dsurfg
							real*8,dimension(Nfl,Nfl) :: Cnx,Cny,Cnt,Me,NiNi
							real*8,dimension(neq,neq) :: W
							real*8,pointer            :: Nfi(:,:)	
							
							! Face shape functions
       Nfi => refElTor%sFTF		

							Cnx = 0.; Cny = 0.
							Me = 0.
							                     
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
																						   
						       NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
						       
												! Contribution of the current integration point to the elemental matrix
													Me = Me + NiNi
						       Cnx = Cnx + NiNi*n_g(g,1)
						       Cny = Cny + NiNi*n_g(g,2)
						       Cnt = Cnt + NiNi*n_g(g,3)
										END DO
						 END DO
							CALL expand_matrix_mass(Me,neq,Df)
							CALL expand_matrix_mass(Me,neq,Ef)							
							CALL expand_matrix_Bt(Cnx,Cny,Cnt,Lf)
							CALL multiply_for_diffusion(Lf)

  END SUBROUTINE set_NeumannH




  !****************************
  ! Bohm
  !****************************    
  SUBROUTINE set_Bohm_bc(tau_save_el,xy_g_save_el)
       integer                   :: delta,i,j,k
       integer                   :: g,igpol,igtor
       real*8                    :: dsurfg       
							real*8                    :: ug(refElTor%Ngl),ufg(refElTor%Ngl,Neq),upgf(refElTor%Ngl,phys%npv)
							real*8                    :: upg(refElTor%Ngl,phys%npv)							
							real*8                    :: inc,sn
							real*8,pointer            :: Nfi(:,:)
							real*8,dimension(Nfl,Nfl) :: NiNi
							real*8,dimension(Nfl,Nfl) :: Massf,BB12,BB21,BB22
							real*8,dimension(Nfl,Nfl) :: Cnx,Cny,Cnt
							real*8,dimension(Nfl,Nfl) :: Qbx,Qby,Qbt
							real*8                    :: Nbn(Nfl)
							real*8                    :: soundSpeed
       real*8                    :: tau_ngamma(2)
							real*8                    :: b_f(refElTor%Ngl,3)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,3),Bmod(Mesh%Nnodesperelem)
#endif				
#ifdef TEMPERATURE
       real*8                    :: qfg(refElTor%Ngl,neq*3)
							real*8                    :: NN_f(Nfl,Nfl),NNf_f(Nfl)
							real*8,dimension(2,Neq)   :: Abohm
						
       real*8,dimension(Nfl,Nfl) :: Massf_en
							real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),taui(Ndim,Neq),gmi	
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),taue(Ndim,Neq),gme
       real*8                    :: Qpr(Ndim,Neq),W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
#endif
       real*8,dimension(Nfl,Nfl) :: Df_loc(Neq*Nfl,Neq*Nfl),Ef_loc(Neq*Nfl,Neq*Nfl)
       real*8                    :: tau_temp(2),bn
       real*8,dimension(Nfl,Nfl) :: aux_loc(Neq*Nfl,Neq*Nfl)
       logical                   :: ntang
       real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)       
       real                      :: tau_stab(Neq), tau_stabmat(Neq,Neq)
       integer                   :: ind_asf(Nfl),ind_ash(Nfl),ind(Nfl)
       
       
       
       ind_asf = (/ (i, i = 0, Neq*(Nfl-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Nfl-1)*Ndim,Neq*Ndim) /)
       
							Massf = 0.d0
							BB12  = 0.
							BB21  = 0.
							BB22  = 0.
							Cnx   = 0.
							Cny   = 0.
							Cnt   = 0.
							Qbx   = 0.
							Qby   = 0.
							Qbt   = 0.
							
							! Face shape functions
       Nfi => refElTor%sFTF		
       							
#ifdef TEMPERATURE							         
        Massf_en = 0.
       ! Gradient solution at face gauss points
       qfg    = matmul(Nfi,qqf)
#endif
										
       ! Set the face solution           
       ufg = matmul(Nfi,transpose(reshape(uf,[neq,Nfl])))		
       IF (numer%stab>1) THEN
            CALL cons2phys(ufg,upgf)
       END IF
              
       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)

! VEDERE PERCHÈ CALCOLO DUE VOLTE LE VARIABILI FISICHE NEI PUNTI DI GAUSS!!!!!!!!!!!!!!!!!!
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
							
							
													Df_loc = 0.; Ef_loc = 0.
													aux_loc = 0.
						
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
													
													IF (abs(inc) .le. phys%bohmth ) THEN
																	sn = 1.
																	delta = 0
#ifdef TEMPERATURE
                 ntang = .FALSE.
#endif													 
													ELSE
																	sn = sign(1.,inc)
						            IF ((sn*upg(g,2)).lt.SoundSpeed) THEN         
						               delta = 1
						            ELSE
						               delta = 0
						            END IF
													END IF						   									   
          
												! Impose the velocity at the sound speed in case of subsonic 
						      ! velocity. 
													NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
							      ! Face mass matrix (to create the stabilization terms)
													Massf =  NiNi
											  BB21  =  (delta*sn*SoundSpeed )*NiNi
													BB22  =  (1-delta)*NiNi												
								  
								  !****************** Stabilization part******************************
            IF (numer%stab>1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab<5) THEN
                  CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),iel,ifa,1.,xyf(g,:),tau_stab)
               ELSE
                  tau_stabmat = 0.
                  CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),1.,iel,tau_stabmat)
               ENDIF
               if (save_tau) then
                  if (numer%stab<5) then
                     tau_save_el(g,:) = tau_stab
                  else
                     tau_save_el(g,1) = tau_stabmat(1,1)
                     tau_save_el(g,2) = tau_stabmat(2,2)
                     tau_save_el(g,3) = tau_stabmat(3,3)
                     tau_save_el(g,4) = tau_stabmat(4,4)
                  endif
                  xy_g_save_el(g,:) = xyf(g,:)
               endif
          !****************** End stabilization part******************************     
            ELSE
               tau_stab = 1.
            END IF			

          ! Case diagonal matrix stabilization
          IF (numer%stab<5) THEN  			   
             DO i=1,Neq
                ind = i+ind_asf
                Ef(ind,ind) = Ef(ind,ind) + tau_stab(i)*NiNi
                IF (i==2) THEN
                   Df(ind,ind-1) = Df(ind,ind-1) + tau_stab(i)*(delta*sn*SoundSpeed )*NiNi
                   Df(ind,ind) = Df(ind,ind) + tau_stab(i)*(1-delta)*NiNi
                ELSE
                   Df(ind,ind) = Df(ind,ind) + tau_stab(i)*NiNi
                END IF
             END DO    
          ELSE
          ! Case full matrix stabilization
             DO i=1,Neq
                DO j=1,Neq
                   Ef(i+ind_asf,j+ind_asf) = Ef(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi
                   IF (j==1) THEN
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi+tau_stabmat(i,j+1)*(delta*sn*SoundSpeed )*NiNi
                   ELSE IF (j==2) THEN 
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*(1-delta)*NiNi
                   ELSE
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi
                   ENDIF 
                END DO
             END DO 
          ENDIF           			   
                        			       
#ifdef TEMPERATURE
          IF (ntang) THEN
                        
						      ! Compute Q^T^(k-1)
						      Qpr = reshape(qfg(g,:),(/Ndim,Neq/))

						      ! Compute V(U^(k-1)) 
						      call computeVi(ufg(g,:),Vveci)
            call computeVe(ufg(g,:),Vvece)
            
						      ! Compute dV_dU (k-1)
						      call compute_dV_dUi(ufg(g,:),dV_dUi)
						      call compute_dV_dUe(ufg(g,:),dV_dUe)

						      ! Compute Alpha(U^(k-1))
						      Alphai = computeAlphai(ufg(g,:))
						      Alphae = computeAlphae(ufg(g,:))

						      ! Compute dAlpha/dU^(k-1)
						      call compute_dAlpha_dUi(ufg(g,:),dAlpha_dUi)
						      call compute_dAlpha_dUe(ufg(g,:),dAlpha_dUe)

											! Tensor products between shape functions    
				        NNf_f =  Nfi(g,:)*dsurfg*bn
				        NN_f  =  NiNi*bn

            CALL jacobianMatricesBohm(ufg(g,:),Abohm)

            gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))  ! scalar
            gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))             ! scalar 
            Taui = matmul(Qpr,dV_dUi)                      ! Ndim x Neq
            Taue = matmul(Qpr,dV_dUe)      ! Ndim x Neq

             DO i=1,2
                IF (i==1) THEN
                   DO j=1,Neq
                      Hf(i+2+ind_asf,j+ind_asf) = Hf(i+2+ind_asf,j+ind_asf) + Abohm(i,j)*NN_f
                      TUhf(i+2+ind_asf,j+ind_asf) = TUhf(i+2+ind_asf,j+ind_asf) + coefi*(gmi*dAlpha_dUi(j)+Alphai*( dot_product(Taui(:,j),b_f(g,:))  ))*NN_f
                      DO k=1,Ndim
                         TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) + coefi*Alphai*Vveci(j)*b_f(g,k)*NN_f
                      END DO
                    END DO
                    Tfhf(i+2+ind_asf) = Tfhf(i+2+ind_asf) + coefi*Alphai*( dot_product (matmul(transpose(Taui),b_f(g,:)),ufg(g,:))  )*NNf_f
                ELSE
                   DO j=1,Neq
                      Hf(i+2+ind_asf,j+ind_asf) = Hf(i+2+ind_asf,j+ind_asf) + Abohm(i,j)*NN_f
                      TUhf(i+2+ind_asf,j+ind_asf) = TUhf(i+2+ind_asf,j+ind_asf) + coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b_f(g,:)) ))*NN_f
                      DO k=1,Ndim
                         TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) + coefe*Alphae*Vvece(j)*b_f(g,k)*NN_f
                      END DO
                   END DO
                   Tfhf(i+2+ind_asf) = Tfhf(i+2+ind_asf) + coefe*Alphae*( dot_product (matmul(transpose(Taue),b_f(g,:)),ufg(g,:))  )*NNf_f 
                ENDIF
             END DO
           
									END IF ! tangency
#endif          

         ! This is to fill Lf: impose derivative of the variables
         ! equal to zero at boundary
          Cnx = Cnx + NiNi*n_g(g,1)
          Cny = Cny + NiNi*n_g(g,2)
          Cnt = Cnt + NiNi*n_g(g,3)
          Qbx = Qbx + b_f(g,1)*NiNi*bn
          Qby = Qby + b_f(g,2)*NiNi*bn  
          Qbt = Qbt + b_f(g,3)*NiNi*bn
							END DO
				END DO
       ! Assembly of the Neumann part of the condition
       CALL expand_matrix_Bt(Cnx,Cny,Cnt,Lf)  ! this applies to all variables
       CALL expand_matrix_Bt(Qbx,Qby,Qbt,Qf)  ! this applies to all variables

       !! ATTENTION THE FOLLOWING SHOULD BE DONE AT EACH GAUSS POINTS LEVEL
       if (ntang) then
          CALL multiply_for_diffusion(Lf)
          CALL multiply_for_diffusion(Qf)
       endif  
          
       
  END SUBROUTINE set_Bohm_bc
  
  

  
    
  

!*******************************************
!           AUXILIARY ROUTINES
!*******************************************  				 
							 
								!*******************************************
								! Expand matrix: mass type
								!*******************************************									
								SUBROUTINE expand_matrix_mass(A,n,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
										integer*4            :: i,j,k,m
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
													      B(n*(i-1)+k,n*(j-1)+k) = A(i,j)
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_mass
								
								
								
								!*******************************************
								! Expand matrix: mass type with different 
								! stabilization parameters
								!*******************************************									
								SUBROUTINE expand_matrix_mass_stab(A,n,tau,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:)
										real*8,intent(out)   :: B(1:n*size(A,1),1:n*size(A,1))
										integer*4            :: i,j,k,m
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
													      B(n*(i-1)+k,n*(j-1)+k) = tau(k)*A(i,j)
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_mass_stab								
							 
#ifdef N2EQ							 
								!*******************************************
								! Expand matrix: convection type 2eq
								!*******************************************									
								SUBROUTINE expand_matrix_conv(C11,C12,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*2,size(C11,1)*2) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*2+1,(j-1)*2+1) = C11(i,j)
             C((i-1)*2+1,(j-1)*2+2) = C12(i,j)
             C((i-1)*2+2,(j-1)*2+1) = C21(i,j)
             C((i-1)*2+2,(j-1)*2+2) = C22(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix_conv	
#endif								
								
								
								
								!*******************************************
								! Expand matrix B transpose
								!*******************************************									
								!*******************************************
								! Expand matrix B
								!*******************************************									
								SUBROUTINE expand_matrix_B(Bx,By,Bz,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:),Bz(:,:)
										real*8,intent(out)   :: B(Neq*3*size(Bx,1),Neq*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										B = 0.d0
										DO j= 1,m
													DO i = 1,n
													   DO k = 1,Neq
 																		B((i-1)*Neq*3+1+(k-1)*3,(j-1)*Neq+k) = Bx(i,j)
 																		B((i-1)*Neq*3+2+(k-1)*3,(j-1)*Neq+k) = By(i,j)
 																		B((i-1)*Neq*3+3+(k-1)*3,(j-1)*Neq+k) = Bz(i,j)
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_B								
								


								!*******************************************
								! Expand matrix B transpose
								!*******************************************									
								SUBROUTINE expand_matrix_Bt(Bx,By,Bz,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:),Bz(:,:)
										real*8,intent(out)   :: B(Neq*size(Bx,1),Neq*3*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										B = 0.d0
										DO j = 1,m
													DO i = 1,n
																DO k=1,Neq
																			B((i-1)*Neq+k,(j-1)*Neq*3+1+(k-1)*3) = Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*3+2+(k-1)*3) = By(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*3+3+(k-1)*3) = Bz(i,j)																			
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt											
        
								SUBROUTINE expand_matrix_Bt_onlyNGamma(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(inout) :: B(Neq*size(Bx,1),Neq*Ndim*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										DO j = 1,m
													DO i = 1,n
																DO k=1,2
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) + Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) + By(i,j)																																					
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt_onlyNGamma							        


								SUBROUTINE expand_matrix_Bt_onlyTiTe(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(inout) :: B(Neq*size(Bx,1),Neq*Ndim*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										DO j = 1,m
													DO i = 1,n
																DO k=3,4
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) + Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) + By(i,j)																																					
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt_onlyTiTe
								
								!*******************************************
								! Multiply matrix for diffusion
								!*******************************************						
       SUBROUTINE multiply_for_diffusion(mat)
       real*8, intent(inout) :: mat(:,:)
       integer               :: n,i
       
       n = size(mat,1)/Neq
       
       DO i = 1,n
          mat((i-1)*Neq+1,:) = mat((i-1)*Neq+1,:)*phys%diff_n
          mat((i-1)*Neq+2,:) = mat((i-1)*Neq+2,:)*phys%diff_u        
#ifdef TEMPERATURE
          mat((i-1)*Neq+3,:) = mat((i-1)*Neq+3,:)*phys%diff_e     
          mat((i-1)*Neq+4,:) = mat((i-1)*Neq+4,:)*phys%diff_ee     
#endif          
       END DO
       
       END SUBROUTINE multiply_for_diffusion			
       					



								


#ifdef TEMPERATURE

								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av(C11,C21,C22,C33,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22,C33,C44
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j) 
             C((i-1)*Neq+3,(j-1)*Neq+3) = C33(i,j)         
             C((i-1)*Neq+4,(j-1)*Neq+4) = C44(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av			
								

								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)  
             C((i-1)*Neq+3,(j-1)*Neq+3) = Mat(i,j)
             C((i-1)*Neq+4,(j-1)*Neq+4) = Mat(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab(C11,C21,C22,C33,C44,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22,C33,C44
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)*tau(2) 
             C((i-1)*Neq+3,(j-1)*Neq+3) = C33(i,j)*tau(3)         
             C((i-1)*Neq+4,(j-1)*Neq+4) = C44(i,j)*tau(4)
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab			
								




								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)*tau(2)  
             C((i-1)*Neq+3,(j-1)*Neq+3) = Mat(i,j)*tau(3)
             C((i-1)*Neq+4,(j-1)*Neq+4) = Mat(i,j)*tau(4)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab								
								
								


					



								!*****************************************************
								! Expand matrix for Bohm: Temperatures part
								!*****************************************************			
								SUBROUTINE expand_matrix_mass_stab_temperature(A,tau,C)        
								real*8,intent(in),dimension(:,:)      :: A
        real*8,intent(in),dimension(:)        :: tau
								real*8,intent(inout),dimension(size(A,1)*Neq,size(A,2)*Neq) :: C
								integer*4 :: i,j,m,k
								
								m=size(A,1)

        DO i = 1,m
           DO j = 1,m
              DO k =1,2
                 C(Neq*(i-1)+k+2,Neq*(j-1)+k+2) = C(Neq*(i-1)+k+2,Neq*(j-1)+k+2) + tau(k)*A(i,j)
              END DO
           END DO
        END DO
								END SUBROUTINE expand_matrix_mass_stab_temperature
								
        
        
								!*******************************************
								! Expand matrix for Bohm: H matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmH(C31,C32,C33,C34,C41,C42,C43,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C41,C42,C43,C44
								real*8,intent(inout),dimension(size(C31,1)*4,size(C31,2)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(C31,1)
								m=size(C31,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+3,(j-1)*4+1) = C31(i,j)
             C((i-1)*4+3,(j-1)*4+2) = C32(i,j)
             C((i-1)*4+3,(j-1)*4+3) = C33(i,j)
             C((i-1)*4+3,(j-1)*4+4) = C34(i,j)
             C((i-1)*4+4,(j-1)*4+1) = C41(i,j)
             C((i-1)*4+4,(j-1)*4+2) = C42(i,j)
             C((i-1)*4+4,(j-1)*4+3) = C43(i,j)
             C((i-1)*4+4,(j-1)*4+4) = C44(i,j)
             
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmH
								
								
								
								!*******************************************
								! Expand matrix for Bohm: TQ matrix
								!*******************************************											
								SUBROUTINE expand_matrix_BohmTq(C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48,C)        
								real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48
								real*8,intent(inout),dimension(size(C31,1)*neq,size(C31,2)*Neq*Ndim) :: C
								integer*4 :: i,j,n,m
								n=size(C31,1)
								m=size(C31,2)
								DO i=1,n
								   DO j=1,m
								     C((i-1)*4+3,(j-1)*8+1) = C31(i,j)
								     C((i-1)*4+3,(j-1)*8+2) = C32(i,j)
								     C((i-1)*4+3,(j-1)*8+3) = C33(i,j)             
								     C((i-1)*4+3,(j-1)*8+4) = C34(i,j)
								     C((i-1)*4+3,(j-1)*8+5) = C35(i,j)
								     C((i-1)*4+3,(j-1)*8+6) = C36(i,j)
								     C((i-1)*4+3,(j-1)*8+7) = C37(i,j)
								     C((i-1)*4+3,(j-1)*8+8) = C38(i,j)
								     
								     C((i-1)*4+4,(j-1)*8+1) = C41(i,j)
								     C((i-1)*4+4,(j-1)*8+2) = C42(i,j)
								     C((i-1)*4+4,(j-1)*8+3) = C43(i,j)             
								     C((i-1)*4+4,(j-1)*8+4) = C44(i,j)
								     C((i-1)*4+4,(j-1)*8+5) = C45(i,j)
								     C((i-1)*4+4,(j-1)*8+6) = C46(i,j)
								     C((i-1)*4+4,(j-1)*8+7) = C47(i,j)
								     C((i-1)*4+4,(j-1)*8+8) = C48(i,j)								     
								   END DO
								END DO
								END SUBROUTINE expand_matrix_BohmTq				
        
								!*******************************************
								! Expand matrix for Bohm (and Dirichlet) : Tf matrix
								!*******************************************										
								SUBROUTINE expandTf(Tf3,Tf4,Tf)
								real*8, intent(IN) :: Tf3(:),Tf4(:)
								real*8, intent(INOUT):: Tf(:)
								integer            :: i,n								
								n = size(Tf3)
								DO i=1,n
											Tf((i-1)*4+3) = Tf3(i)
											Tf((i-1)*4+4) = Tf4(i)
								END DO				
								END SUBROUTINE expandTf									
#else

								!*******************************************
								! Expand matrix for Bohm: D matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD(C11,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD				
								




								!*******************************************
								! Expand matrix for Bohm: E matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE		
								
								
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av(C11,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j) 
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av			
								


								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)         
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av
								
								
						
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab(C11,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)*tau(2) 
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab			
								




								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)*tau(2)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab								
								
#endif						
						
						
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with tensor stabilization 
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab_mat(Ms,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Ms,C21,C22,tau
								real*8,intent(inout),dimension(size(Ms,1)*Neq,size(Ms,2)*Neq) :: C
								integer*4 :: i,j,n,m,k
								
								n=size(Ms,1)
								m=size(Ms,2)

        DO i=1,n
           DO j=1,m
              DO k=1,Neq
						          C((i-1)*Neq+k,(j-1)*Neq+1) = Ms(i,j)*tau(k,1) + C21(i,j)*tau(k,2)
						          C((i-1)*Neq+k,(j-1)*Neq+2) = C22(i,j)*tau(k,2)
						          C((i-1)*Neq+k,(j-1)*Neq+3) = Ms(i,j)*tau(k,3)         
						          C((i-1)*Neq+k,(j-1)*Neq+4) = Ms(i,j)*tau(k,4)
              END DO
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab_mat	
								
							


								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with matrix stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab_mat(A,tau,B)        
										real*8,intent(in)    :: A(:,:),tau(:,:)
										real*8,intent(out)   :: B(1:Neq*size(A,1),1:Neq*size(A,1))
										integer*4            :: i,j,k,t,m,n
										
										n = size(tau,1)
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
                   DO t = 1,n
													         B(n*(i-1)+k,n*(j-1)+t) = tau(k,t)*A(i,j)
                   END DO
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab_mat
								
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix - 
        ! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_stab(C11,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
        real*8,intent(in),dimension(2)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C((i-1)*Neq+1,(j-1)*Neq+1) + C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C((i-1)*Neq+2,(j-1)*Neq+1) + C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C((i-1)*Neq+2,(j-1)*Neq+2) + C22(i,j)*tau(2)   
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_stab				
								




								!*******************************************
								! Expand matrix for Bohm: E matrix - 
        ! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
        real*8,intent(in),dimension(2)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C((i-1)*Neq+1,(j-1)*Neq+1) + Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C((i-1)*Neq+2,(j-1)*Neq+2) + Mat(i,j)*tau(2)   
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_stab		
        

								

END SUBROUTINE HDG_BC
  
  
  
#else
!***********************************************************************
! 
!                            VERSION 2D 
! 
!***********************************************************************
  integer                   :: ifa,ifl,iel,fl,bc,Nf,Nfp,Neq,Np,i,Ndim,Fi
  integer                   :: nod(refElPol%Nfacenodes)
  integer                   :: ind_ff(refElPol%Nfacenodes*phys%neq),ind_fl(refElPol%Nfacenodes*phys%neq)
  integer                   :: ind_fe(refElPol%Nfacenodes*phys%neq)
  real*8                    :: uf(refElPol%Nfacenodes*phys%neq)
  real*8                    :: Xf(Mesh%Nnodesperface,Mesh%Ndim)
  real*8, allocatable       :: Df(:,:),Ef(:,:),fh(:),Edirf(:),Hdirf(:)
  integer                   :: ind_fG(refElPol%Nfacenodes*phys%neq*Mesh%Ndim)
  real*8, allocatable       :: Lf(:,:),Qf(:,:),Cdirf(:)
#ifdef TEMPERATURE
  integer                   :: ind_qf(refElPol%Nfacenodes*phys%neq*Mesh%Ndim)
  real*8                    :: qqf(refElPol%Nfacenodes,phys%neq*Mesh%Ndim), qfg(refElPol%Ngauss1d,phys%neq*Mesh%Ndim)
  real*8, allocatable       :: TUhf(:,:),TQhf(:,:),Tfhf(:),Thdirf(:),Hf(:,:)
  real*8                    :: coefi,coefe  
#endif  

  real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
  logical                   :: save_tau
  real*8,allocatable        :: tau_save(:,:)
  real*8,allocatable        :: xy_g_save(:,:)
  integer                   :: indtausave(refElPol%Ngauss1d)
  
  
  save_tau = switch%saveTau
  Ndim= Mesh%Ndim
  Np  = refElPol%Nnodes2D
  Nfp = refElPol%Nfacenodes
  Neq = phys%neq
  
  ALLOCATE(Df(Neq*Nfp,Neq*Nfp))
  ALLOCATE(Ef(Neq*Nfp,Neq*Nfp))
  ALLOCATE(fh(Neq*Nfp))
  ALLOCATE(Edirf(Neq*Nfp)) 
  ALLOCATE(Hdirf(Neq*Nfp))
  ALLOCATE(Cdirf(Neq*Nfp*ndim)) 
  ALLOCATE(Lf(Neq*Nfp,Neq*Nfp*Ndim))
  ALLOCATE(Qf(Neq*Nfp,Neq*Nfp*Ndim))
#ifdef TEMPERATURE
  ALLOCATE(TUhf(Neq*Nfp,Neq*Nfp))
  ALLOCATE(TQhf(Neq*Nfp,Neq*Ndim*Nfp))
  ALLOCATE(Tfhf(Neq*Nfp))
  ALLOCATE(Thdirf(Neq*Nfp))
  ALLOCATE(Hf(Neq*Nfp,Neq*Nfp))
  coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
  coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
#endif  
!  
!  
!  write(6,*) "coef: ", coef
!  write(6,*) "coefe: ", coefe  
  if (save_tau) then
     allocate(tau_save(Mesh%Nextfaces*refElPol%Ngauss1d,phys%neq))
     allocate(xy_g_save(Mesh%Nextfaces*refElPol%Ngauss1d,2))
     tau_save = 0.
     xy_g_save = 0.
  endif 
    
  !************************
  ! Loop in external faces
  !************************
  DO ifa = 1,Mesh%Nextfaces
     
     ! Initialize local matrices
     Df = 0.
     Ef = 0.
     Edirf = 0.d0
     Hdirf = 0.d0
     fh = 0.
     Cdirf = 0.
     Lf = 0.
     Qf = 0.
#ifdef TEMPERATURE
     TUhf = 0.
     TQhf = 0.
     Tfhf = 0.
     Thdirf = 0.
     Hf   = 0.
#endif      
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
     ind_ff = (Fi-1)*Neq*Nfp + (/(i,i=1,Neq*Nfp)/)
     ind_fl = (ifl-1)*Neq*Nfp + (/(i,i=1,Neq*Nfp)/)
     ind_fe = reshape(tensorSumInt((/(i,i=1,neq)/),neq*(nod-1) ) , (/neq*Nfp/))
     ind_fG = reshape(tensorSumInt((/(i,i=1,neq*ndim)/),neq*ndim*(refElPol%face_nodes(ifl,:)-1)),(/neq*ndim*Nfp/))
#ifdef TEMPERATURE
     ind_qf = (iel-1)*Ndim*Neq*Np + ind_fG
     qqf    = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Nfp /)))
#endif     
     
     ! Solution in this face
     uf = sol%u_tilde(ind_ff)

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
					CASE(bc_NeumannH)
					   CALL set_NeumannH()	
     CASE(bc_Bohm)
        
        
!        if ( iel.eq.28725 ) then
!        CALL set_NeumannH()		
!        else 
        
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
					elMat%Df(ind_fl,ind_fe,iel) = Df
					elMat%Ef(ind_fl,ind_fl,iel) = Ef 		
					elMat%fh(ind_fl,iel)        = fh
					IF (Mesh%ndir.gt.0) THEN					
					   elMat%Edir(ind_fe,iel) = Edirf
					   elMat%Hdir(ind_fe,iel) = Hdirf
					END IF					
    	elMat%Lf(ind_fl,ind_fG,iel) = Lf
    	elMat%Qf(ind_fl,ind_fG,iel) = Qf
					IF (Mesh%ndir.gt.0) THEN					    	
    	   elMat%Cdir(ind_fG,iel) = Cdirf
#ifdef TEMPERATURE     
        elMat%Thdir(ind_fe,iel)= Thdirf
#endif    	   
    	END IF
#ifdef TEMPERATURE
     elMat%TQhf(ind_fl,ind_fG,iel)  = TQhf
     elMat%Tfhf(ind_fl,iel)         = Tfhf
     elMat%Hf(ind_fl,ind_fl,iel)    = Hf-TUhf
#endif    	
					
!if (iel==2225) then
!!  call saveMatrix(elMat%TQhf(:,:,iel),'TQhf')
!!  call saveVector(elMat%Tfhf(:,iel),'Tfhf')    
!!  call saveMatrix(elMat%Hf(:,:,iel),'Hf')    
!!  call saveMatrix(elMat%Ef(:,:,iel),'Ef')  
!!  call saveMatrix(elMat%Df(:,:,iel),'Df')
!  call saveMatrix(elMat%Lf(:,:,iel),'Lfe')
!  stop
!!  call saveMatrix(elMat%Qf(:,:,iel),'Qf')
!!  stop
!endif   

  END DO

   if (save_tau) then
      write(6,*) "Saving tau in the boundary faces"
      call saveMatrix(tau_save,'tau_save_bound')
      call saveMatrix(xy_g_save,'xy_g_save_bound')
      deallocate(tau_save,xy_g_save)
      write(6,*) "Done saving tau!"
   endif
     
  DEALLOCATE(Df,Ef,Edirf,Hdirf,fh)

  DEALLOCATE(Cdirf,Lf,Qf)
#ifdef TEMPERATURE
  DEALLOCATE(Thdirf,TQhf,TUhf,Tfhf,Hf) 
#endif   
  CONTAINS
  
  !****************************
  ! Dirichlet in weak form
  !****************************  
  SUBROUTINE set_dirichletwf_bc
       integer                   :: g,Ngauss
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: uex(refElPol%Ngauss1d,neq)       
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
							real*8                    :: f(neq,Nfp)
							real*8,dimension(Nfp,Nfp) :: Massf
							
							
							NGauss = refElPol%Ngauss1D
							Massf = 0.d0
							f = 0.

       ! Set the face solution           
       xyf    = matmul( refElPol%N1D,Xf )
       CALL analytical_solution(xyf(:,1),xyf(:,2),uex) 
       							                     
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)
       
							! Loop in 1D Gauss points
							DO g = 1,NGauss
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
									! Contribution of the current integration point to the elemental matrix
         Massf = Massf + dline*tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))
         f = f + dline*tensorProduct(uex(g,:),refElPol%N1D(g,:))
							END DO

							CALL expand_matrix_mass_stab(Massf,neq,numer%tau,Ef)
							
							! Multiply for stabilization before reshaping
							DO i=1,neq
							   f(i,:) = numer%tau(i)*f(i,:)
							END DO
       fH = -reshape(f,(/Nfp*neq/))
    
  END SUBROUTINE set_dirichletwf_bc




  !****************************
  ! Dirichlet in strong form
  !****************************  
  SUBROUTINE set_dirichletstr_bc
       integer                   :: g,Ngauss,ind_asf(Nfp),i,j
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
       real*8                    :: uex(refElPol%Ngauss1d,neq)
							real*8                    :: t_g(ndim),n_g(ndim),fh(Nfp*neq)
							real*8,dimension(neq,neq) :: Ax,Ay
							real*8,dimension(Nfp,neq) :: Edir,Hdir
       real*8                    :: b_f(refElPol%Ngauss1d,ndim),divb_f(refElPol%Ngauss1d),drift_f(refElPol%Ngauss1d,ndim)    
       real*8                    :: Cdir(ndim*neq,Nfp)
       real*8                    :: uexn(ndim*neq),bn
       real*8,dimension(neq,neq) :: An
#ifdef TEMPERATURE
							real*8                    :: Qpr(Ndim,Neq)
							real*8                    :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
							real*8                    :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)							
#endif      
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif		       
							Edir = 0.d0
							Hdir = 0.d0
       Cdir = 0.
							ind_asf = (/ (i, i = 0, Neq*(Nfp-1),Neq) /)
							NGauss = refElPol%Ngauss1D
							                     
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
							DO g = 1,NGauss
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
							   ! Unit normal to the boundary
							   t_g = xyDer(g,:)/xyDerNorm_g
							   n_g = [t_g(2), -t_g(1)]

							   ! Jacobian matrices									   
							   bn = dot_product(b_f(g,:),n_g)
							   CALL jacobianMatricesFace(uex(g,:),bn,An) 
							    
							   ! Contribution of the current integration point to the elemental matrix
 						   Edir = Edir + dline*tensorProduct(refElPol%N1D(g,:),uex(g,:))			 
		        Hdir = Hdir + tensorProduct( refElPol%N1D(g,:),matmul(An,uex(g,:)) )*dline  
		        uexn = reshape(tensorProduct(n_g,uex(g,:)),(/ndim*neq/))
		        Cdir = Cdir + dline*tensorProduct(uexn,refElPol%N1D(g,:))        
#ifdef TEMPERATURE        
         ! Compute Q^T^(k-1)
         Qpr = reshape(qfg(g,:),(/Ndim,Neq/))

         ! Compute V(U^(k-1)) 
         call computeVi(uex(g,:),Vveci)
         call computeVe(uex(g,:),Vvece)
         
			      ! Compute dV_dU (k-1)
			      call compute_dV_dUi(uex(g,:),dV_dUi)
			      call compute_dV_dUe(uex(g,:),dV_dUe)
			      
         ! Compute Gamma and tau
		       gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))   ! scalar    
		       gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))  
		       Taui = matmul(Qpr,dV_dUi)                       ! 2x3
		       Taue = matmul(Qpr,dV_dUe)

						   ! Compute Alpha(U^(k-1))
						   Alphai = computeAlphai(uex(g,:))  
						   Alphae = computeAlphae(uex(g,:))       
						   
         ! Compute dAlpha/dU^(k-1)
         call compute_dAlpha_dUi(uex(g,:),dAlpha_dUi)
         call compute_dAlpha_dUe(uex(g,:),dAlpha_dUe)
         
         DO i=3,4
            DO j=1,4
               IF (i==3) THEN
                  Thdirf(i+ind_asf) = Thdirf(i+ind_asf) + coefi*refElPol%N1D(g,:)*dline*bn*((gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),b_f(g,:))))*uex(g,j))
               ELSE
                  Thdirf(i+ind_asf) = Thdirf(i+ind_asf) + coefe*refElPol%N1D(g,:)*dline*bn*((gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),b_f(g,:))))*uex(g,j))
               END IF
            END DO
         END DO
#endif           
							END DO							
							
							! Elemental assembly
							! Multiply for stabilization before reshaping
							DO i=1,neq
							   Edir(:,i) = numer%tau(i)*Edir(:,i)
							END DO
       Edirf =  reshape(transpose(Edir),(/neq*Nfp/)) 
       Hdirf =  reshape(transpose(Hdir),(/neq*Nfp/)) 
       Cdirf =  reshape(Cdir,(/neq*Nfp*ndim/))
  END SUBROUTINE set_dirichletstr_bc
    



  !******************************************
  ! Dirichlet in weak form: velocity at 
  ! sound speed and threshold value for the
  ! rest of the variables
  !******************************************  
  SUBROUTINE set_dirichletwfbohm_bc
       integer                   :: g,Ngauss
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
							real*8                    :: f(neq,Nfp),soundspeed
							real*8                    :: t_g(ndim),n_g(ndim),bn,sn,delta,inc
							real*8                    :: b_f(refElPol%Ngauss1d,ndim),divb_f(refElPol%Ngauss1d),drift_f(refElPol%Ngauss1d,ndim)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif	
							real*8,dimension(Nfp,Nfp) :: Massf,zeros,BB21,BB22,NiNi
							real*8,dimension(Neq)     :: vals
							real*8                    :: ug(refElPol%Ngauss1d),ufg(refElPol%Ngauss1d,Neq),upgf(refElPol%Ngauss1d,phys%npv)
							real*8                    :: upg(refElPol%Ngauss1d,phys%npv)									
							
							NGauss = refElPol%Ngauss1D
							Massf = 0.d0
							zeros = 0.
							BB21  = 0.
							BB22  = 0.
							f = 0.
							vals = 0.
							
!							vals(1) = numer%thr
!							vals(3) = numer%thrpre
!							vals(4) = numer%thrpre
							

       ! Set the face solution           
       xyf    = matmul( refElPol%N1D,Xf )
       							                     
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)
       
       ! Set the face solution           
       ufg = matmul(refElPol%N1D,transpose(reshape(uf,[neq,Nfp])))		
              
       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)

#ifndef MOVINGEQUILIBRIUM	       
			    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,nod),:))
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
							DO g = 1,NGauss

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
							   inc = bn/norm2(b_f(g,:))
										IF (abs(inc) .le. phys%bohmth ) THEN
													 sn = 1.
													 delta = 0
										ELSE
													 sn = sign(1.,inc)
               IF ((sn*upg(g,2)).lt.SoundSpeed) THEN         
                  delta = 1
               ELSE
                  delta = 0
               END IF
										END IF	               
!          vals(2) = sn*soundspeed*vals(1)               
               			      			   
									! Contribution of the current integration point to the elemental matrix
									NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
									Massf =  Massf + NiNi
					    BB21  =  BB21 + (delta*sn*SoundSpeed )*NiNi									
					    BB22  =  BB22 + (1-delta)*NiNi
         f = f + dline*tensorProduct(vals,refElPol%N1D(g,:))
							END DO
#ifdef TEMPERATURE
							CALL expand_matrix_BohmD_av(zeros,BB21,BB22,zeros,zeros,Df)
#else
							CALL expand_matrix_BohmD_av(zeros,BB21,BB22,Df)
#endif							
							CALL expand_matrix_BohmE_av(Massf,Ef)
       fH = -reshape(f,(/Nfp*neq/))
    
  END SUBROUTINE set_dirichletwfbohm_bc
  
  
      















  !************************************************
  ! Bohm for the velocity,  Neumann homogeneous for 
  ! the other variables
  !************************************************
  SUBROUTINE set_BohmNM_bc
       integer                   :: g,Ngauss,delta,i
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
							real*8                    :: t_g(ndim),n_g(ndim)
							real*8                    :: ug(refElPol%Ngauss1d),ufg(refElPol%Ngauss1d,Neq)
							real*8                    :: inc,sn
							real*8,dimension(Nfp,Nfp) :: NiNi
							real*8,dimension(Nfp,Nfp) :: Massf,BB12,BB21,BB22
							real*8                    :: b_f(refElPol%Ngauss1d,ndim),divb_f(refElPol%Ngauss1d),drift_f(refElPol%Ngauss1d,ndim) 
							real*8,dimension(Nfp,Nfp) :: Cnx,Cny   
							real*8,dimension(Nfp,Nfp) :: Qbx,Qby
							real*8                    :: Nbn(Nfp)
							real*8                    :: soundSpeed
							real*8                    :: upg(refElPol%Ngauss1d,phys%npv)
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif				
#ifdef ENERGY       
       real*8                    :: aux(refElPol%Ngauss1d)
#endif
							Massf = 0.d0
							BB12  = 0.
							BB21  = 0.
							BB22  = 0.
							Cnx   = 0.
							Cny   = 0.
							
							NGauss = refElPol%Ngauss1D
							 
							! Gauss points position                     
						 xyf    = matmul( refElPol%N1D,Xf )
						 
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

       ! Set the face solution           
       ufg = matmul(refElPol%N1D,transpose(reshape(uf,[neq,Nfp])))		
       
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

!#ifndef ENERGY
!       ! Parallel velocity at Gauss points
!       ug = ufg(:,2)/ufg(:,1)
!#else  
       ! Physical variables at Gauss points
       CALL cons2phys(ufg,upg)
!#endif
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
							DO g = 1,NGauss
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN										
				         dline = dline*xyf(g,1)
				      END IF				
							   
							   ! Unit normal to the boundary
							   t_g = xyDer(g,:)/xyDerNorm_g
							   n_g = [t_g(2), -t_g(1)]
							   
							   ! Sound speed
#ifndef ENERGY
          SoundSpeed = sqrt(phys%a)
#else
#ifndef TEMPERATURE
          SoundSpeed = upg(g,6)
#else
          SoundSpeed = upg(g,9)
          if (switch%decoup) then 
             SoundSpeed = sqrt(phys%Mref)
          endif          
#endif
#endif	
							   
							   inc = dot_product(n_g,b_f(g,:))/norm2(b_f(g,:))	   
										IF (ABS(inc) .le. phys%bohmth ) THEN
													 sn = 1.
													 delta = 0
										ELSE
													 sn = sign(1.,inc)
!#ifndef ENERGY													 
!               IF ((sn*ug(g)).lt.sqrt(phys%a)) THEN
!#else
               IF ((sn*upg(g,2)).lt.SoundSpeed) THEN         
!#endif               
                  delta = 1
               ELSE
                  delta = 0
               END IF
               
										END IF						   									   


									! Contribution of the current integration point to the elemental matrix
										NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
										Massf = Massf + NiNi
!#ifndef ENERGY										
!			       BB21 = BB21 + (delta*sn*sqrt(phys%a) )*NiNi
!#else
          BB21 = BB21 + (delta*sn*SoundSpeed )*NiNi
!#endif			       
			       BB22 = BB22 + (1-delta)*NiNi
			       
									! Contribution of the current integration point to the elemental matrix
          Cnx = Cnx + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(1)*dline
          Cny = Cny + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(2)*dline	
          
							END DO

!							CALL expand_matrix_BohmD(Massf,BB21,BB22,Df)
!							CALL expand_matrix_BohmE(Massf,Ef)
#ifdef TEMPERATURE
							CALL expand_matrix_BohmD_av(Massf,BB21,BB22,Massf,Massf,Df)
							CALL expand_matrix_BohmE_av(Massf,Ef)		
#else							
							CALL expand_matrix_BohmD_av(Massf,BB21,BB22,Df)
							CALL expand_matrix_BohmE_av(Massf,Ef)									
#endif							
!							CALL expand_matrix_BohmD_av_stab(Massf,BB21,BB22,Massf,Massf,numer%tau,Df)
!							CALL expand_matrix_BohmE_av_stab(Massf,numer%tau,Ef)									
       CALL expand_matrix_Bt(Cnx,Cny,Lf)
       CALL multiply_for_diffusion(Lf)
       
  END SUBROUTINE set_BohmNM_bc



  !****************************
  ! Neumann homogeneus
  !****************************    
  SUBROUTINE set_NeumannH
       integer                   :: g,Ngauss
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
							real*8                    :: t_g(ndim),n_g(ndim)
							real*8,dimension(Nfp,Nfp) :: Cnx,Cny,Me
							real*8,dimension(neq,neq) :: W
							
							
							NGauss = refElPol%Ngauss1D
							Cnx = 0.; Cny = 0.
							Me = 0.
							                     
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)
       xyf    = matmul( refElPol%N1D,Xf )	
       
							! Loop in 1D Gauss points
							DO g = 1,NGauss
						
							   ! Calculate the integration weight
							   xyDerNorm_g = norm2(xyDer(g,:))
							   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
				      IF (switch%axisym) THEN						         							
				         dline = dline*xyf(g,1)
				      END IF				
							   
							   ! Unit normal to the boundary
							   t_g = xyDer(g,:)/xyDerNorm_g
							   n_g = [t_g(2), -t_g(1)]
							   									   
          
									! Contribution of the current integration point to the elemental matrix
										Me = Me + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
          Cnx = Cnx + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(1)*dline
          Cny = Cny + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(2)*dline			       
							END DO
							CALL expand_matrix_mass(Me,neq,Df)
							CALL expand_matrix_mass(Me,neq,Ef)							
							CALL expand_matrix_Bt(Cnx,Cny,Lf)
							CALL multiply_for_diffusion(Lf)

  END SUBROUTINE set_NeumannH




  !****************************
  ! Bohm
  !****************************    
  SUBROUTINE set_Bohm_bc(tau_save_el,xy_g_save_el)
       integer                   :: g,Ngauss,delta,i,j,k
       real*8                    :: dline,xyDerNorm_g
       real*8                    :: xyf(refElPol%Ngauss1d,ndim),xyDer(refElPol%Ngauss1d,ndim)
							real*8                    :: t_g(ndim),n_g(ndim)
							real*8                    :: ug(refElPol%Ngauss1d),ufg(refElPol%Ngauss1d,Neq),upgf(refElPol%Ngauss1d,phys%npv)
							real*8                    :: upg(refElPol%Ngauss1d,phys%npv)							
							real*8                    :: inc,sn
							real*8,dimension(Nfp,Nfp) :: NiNi
							real*8,dimension(Nfp,Nfp) :: Massf,BB12,BB21,BB22
							real*8,dimension(Nfp,Nfp) :: Cnx,Cny
							real*8,dimension(Nfp,Nfp) :: Qbx,Qby
							real*8                    :: Nbn(Nfp)
							real*8                    :: soundSpeed
       real*8                    :: tau_ngamma(2)
							real*8                    :: b_f(refElPol%Ngauss1d,ndim),divb_f(refElPol%Ngauss1d),drift_f(refElPol%Ngauss1d,ndim)        
#ifdef MOVINGEQUILIBRIUM								
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif				
#ifdef TEMPERATURE
       real*8                    :: qfg(refElPol%Ngauss1d,neq*Ndim)
							real*8                    :: NN_f(Nfp,Nfp),NNf_f(Nfp)
							real*8,dimension(2,Neq)   :: Abohm
						
       real*8,dimension(Nfp,Nfp) :: Massf_en
							real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),taui(Ndim,Neq),gmi	
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),taue(Ndim,Neq),gme
       real*8                    :: Qpr(Ndim,Neq),W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
#endif
       real*8,dimension(Nfp,Nfp) :: Df_loc(Neq*Nfp,Neq*Nfp),Ef_loc(Neq*Nfp,Neq*Nfp)
       real*8                    :: tau_temp(2),bn
       real*8,dimension(Nfp,Nfp) :: aux_loc(Neq*Nfp,Neq*Nfp)
       logical                   :: ntang
       real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)       
       real                      :: tau_stab(Neq), tau_stabmat(Neq,Neq)
       integer                   :: ind_asf(Nfp),ind_ash(Nfp),ind(Nfp)
       
       
       
       ind_asf = (/ (i, i = 0, Neq*(Nfp-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Nfp-1)*Ndim,Neq*Ndim) /)
       
							Massf = 0.d0
							BB12  = 0.
							BB21  = 0.
							BB22  = 0.
							Cnx   = 0.
							Cny   = 0.
							Qbx   = 0.
							Qby   = 0.
							
#ifdef TEMPERATURE							         
        Massf_en = 0.
       ! Gradient solution at face gauss points
       qfg    = matmul(refElPol%N1D,qqf)
#endif
										
							NGauss = refElPol%Ngauss1D
							 
							! Gauss points position                     
						 xyf    = matmul( refElPol%N1D,Xf )
						 
							! Shape function derivatives at Gauss points
							xyDer  = matmul(refElPol%Nxi1D,Xf)

       ! Set the face solution           
       ufg = matmul(refElPol%N1D,transpose(reshape(uf,[neq,Nfp])))		
       IF (numer%stab>1) THEN
            CALL cons2phys(ufg,upgf)
       END IF
              
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
							DO g = 1,NGauss
							
							
   							Df_loc = 0.; Ef_loc = 0.
   							aux_loc = 0.
						
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
										IF (abs(inc) .le. phys%bohmth ) THEN
													 sn = 1.
													 delta = 0
#ifdef TEMPERATURE
              ntang = .FALSE.
#endif													 
										ELSE
													 sn = sign(1.,inc)
               IF ((sn*upg(g,2)).lt.SoundSpeed) THEN         
                  delta = 1
               ELSE
                  delta = 0
               END IF
										END IF						   									   
          
									! Impose the velocity at the sound speed in case of subsonic 
         ! velocity. 
										NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
	         ! Face mass matrix (to create the stabilization terms)
										Massf =  NiNi
					     BB21  =  (delta*sn*SoundSpeed )*NiNi
								  BB22  =  (1-delta)*NiNi												
								  
								  !****************** Stabilization part******************************
            IF (numer%stab>1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab<5) THEN
                  CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g,iel,ifa,1.,xyf(g,:),tau_stab)
               ELSE
                  tau_stabmat = 0.
                  CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g,xyf(g,:),1.,iel,tau_stabmat)
               ENDIF
               if (save_tau) then
                  if (numer%stab<5) then
                     tau_save_el(g,:) = tau_stab
                  else
                     tau_save_el(g,1) = tau_stabmat(1,1)
                     tau_save_el(g,2) = tau_stabmat(2,2)
                     tau_save_el(g,3) = tau_stabmat(3,3)
                     tau_save_el(g,4) = tau_stabmat(4,4)
                  endif
                  xy_g_save_el(g,:) = xyf(g,:)
               endif
          !****************** End stabilization part******************************     
            ELSE
               tau_stab = 1.
            END IF			

          ! Case diagonal matrix stabilization
          IF (numer%stab<5) THEN  			   
             DO i=1,Neq
                ind = i+ind_asf
                Ef(ind,ind) = Ef(ind,ind) + tau_stab(i)*NiNi
                IF (i==2) THEN
                   Df(ind,ind-1) = Df(ind,ind-1) + tau_stab(i)*(delta*sn*SoundSpeed )*NiNi
                   Df(ind,ind) = Df(ind,ind) + tau_stab(i)*(1-delta)*NiNi
                ELSE
                   Df(ind,ind) = Df(ind,ind) + tau_stab(i)*NiNi
                END IF
             END DO    
          ELSE
          ! Case full matrix stabilization
             DO i=1,Neq
                DO j=1,Neq
                   Ef(i+ind_asf,j+ind_asf) = Ef(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi
                   IF (j==1) THEN
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi+tau_stabmat(i,j+1)*(delta*sn*SoundSpeed )*NiNi
                   ELSE IF (j==2) THEN 
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*(1-delta)*NiNi
                   ELSE
                      Df(i+ind_asf,j+ind_asf) = Df(i+ind_asf,j+ind_asf) + tau_stabmat(i,j)*NiNi
                   ENDIF 
                END DO
             END DO 
          ENDIF           			   
                        			       
#ifdef TEMPERATURE
          IF (ntang) THEN
                        
						      ! Compute Q^T^(k-1)
						      Qpr = reshape(qfg(g,:),(/Ndim,Neq/))

						      ! Compute V(U^(k-1)) 
						      call computeVi(ufg(g,:),Vveci)
            call computeVe(ufg(g,:),Vvece)
            
						      ! Compute dV_dU (k-1)
						      call compute_dV_dUi(ufg(g,:),dV_dUi)
						      call compute_dV_dUe(ufg(g,:),dV_dUe)

						      ! Compute Alpha(U^(k-1))
						      Alphai = computeAlphai(ufg(g,:))
						      Alphae = computeAlphae(ufg(g,:))

						      ! Compute dAlpha/dU^(k-1)
						      call compute_dAlpha_dUi(ufg(g,:),dAlpha_dUi)
						      call compute_dAlpha_dUe(ufg(g,:),dAlpha_dUe)

											! Tensor products between shape functions    
				        NNf_f =  refElPol%N1D(g,:)*dline*bn
				        NN_f  =  NiNi*bn

            CALL jacobianMatricesBohm(ufg(g,:),Abohm)

            gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))  ! scalar
            gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))             ! scalar 
            Taui = matmul(Qpr,dV_dUi)                      ! 2x3
            Taue = matmul(Qpr,dV_dUe)      ! 2x3

             DO i=1,2
                IF (i==1) THEN
                   DO j=1,Neq
                      Hf(i+2+ind_asf,j+ind_asf) = Hf(i+2+ind_asf,j+ind_asf) + Abohm(i,j)*NN_f
                      TUhf(i+2+ind_asf,j+ind_asf) = TUhf(i+2+ind_asf,j+ind_asf) + coefi*(gmi*dAlpha_dUi(j)+Alphai*( dot_product(Taui(:,j),b_f(g,:))  ))*NN_f
                      DO k=1,Ndim
                         TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) + coefi*Alphai*Vveci(j)*b_f(g,k)*NN_f
                      END DO
                    END DO
                    Tfhf(i+2+ind_asf) = Tfhf(i+2+ind_asf) + coefi*Alphai*( dot_product (matmul(transpose(Taui),b_f(g,:)),ufg(g,:))  )*NNf_f
                ELSE
                   DO j=1,Neq
                      Hf(i+2+ind_asf,j+ind_asf) = Hf(i+2+ind_asf,j+ind_asf) + Abohm(i,j)*NN_f
                      TUhf(i+2+ind_asf,j+ind_asf) = TUhf(i+2+ind_asf,j+ind_asf) + coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b_f(g,:)) ))*NN_f
                      DO k=1,Ndim
                         TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhf(i+2+ind_asf,k+(j-1)*Ndim+ind_ash) + coefe*Alphae*Vvece(j)*b_f(g,k)*NN_f
                      END DO
                   END DO
                   Tfhf(i+2+ind_asf) = Tfhf(i+2+ind_asf) + coefe*Alphae*( dot_product (matmul(transpose(Taue),b_f(g,:)),ufg(g,:))  )*NNf_f 
                ENDIF
             END DO
             
                  
									END IF ! tangency
#endif          

         ! This is to fill Lf: impose derivative of the variables
         ! equal to zero at boundary
          Cnx = Cnx + n_g(1)*NiNi
          Cny = Cny + n_g(2)*NiNi
          Qbx = Qbx + b_f(g,1)*NiNi*bn
          Qby = Qby + b_f(g,2)*NiNi*bn 		  
							END DO
							
       ! Assembly of the Neumann part of the condition
       CALL expand_matrix_Bt(Cnx,Cny,Lf)  ! this applies to all variables
       CALL expand_matrix_Bt(Qbx,Qby,Qf)  ! this applies to all variables

       !! ATTENTION THE FOLLOWING SHOULD BE DONE AT EACH GAUSS POINTS LEVEL
       if (ntang) then
          CALL multiply_for_diffusion(Lf)
          CALL multiply_for_diffusion(Qf)
       endif  
          
       
  END SUBROUTINE set_Bohm_bc
  
  

  
    
  

!*******************************************
!           AUXILIARY ROUTINES
!*******************************************  				 
							 
								!*******************************************
								! Expand matrix: mass type
								!*******************************************									
								SUBROUTINE expand_matrix_mass(A,n,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
										integer*4            :: i,j,k,m
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
													      B(n*(i-1)+k,n*(j-1)+k) = A(i,j)
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_mass
								
								
								
								!*******************************************
								! Expand matrix: mass type with different 
								! stabilization parameters
								!*******************************************									
								SUBROUTINE expand_matrix_mass_stab(A,n,tau,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:)
										real*8,intent(out)   :: B(1:n*size(A,1),1:n*size(A,1))
										integer*4            :: i,j,k,m
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
													      B(n*(i-1)+k,n*(j-1)+k) = tau(k)*A(i,j)
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_mass_stab								
							 
#ifdef N2EQ							 
								!*******************************************
								! Expand matrix: convection type 2eq
								!*******************************************									
								SUBROUTINE expand_matrix_conv(C11,C12,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*2,size(C11,1)*2) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*2+1,(j-1)*2+1) = C11(i,j)
             C((i-1)*2+1,(j-1)*2+2) = C12(i,j)
             C((i-1)*2+2,(j-1)*2+1) = C21(i,j)
             C((i-1)*2+2,(j-1)*2+2) = C22(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix_conv	
#endif								
								
								
								
								!*******************************************
								! Expand matrix B transpose
								!*******************************************									
								SUBROUTINE expand_matrix_Bt(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(out)   :: B(Neq*size(Bx,1),Neq*Ndim*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										B = 0.d0
										DO j = 1,m
													DO i = 1,n
																DO k=1,Neq
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) = Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) = By(i,j)																																					
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt											
        
								SUBROUTINE expand_matrix_Bt_onlyNGamma(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(inout) :: B(Neq*size(Bx,1),Neq*Ndim*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										DO j = 1,m
													DO i = 1,n
																DO k=1,2
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) + Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) + By(i,j)																																					
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt_onlyNGamma							        


								SUBROUTINE expand_matrix_Bt_onlyTiTe(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(inout) :: B(Neq*size(Bx,1),Neq*Ndim*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										DO j = 1,m
													DO i = 1,n
																DO k=3,4
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+1+(k-1)*Ndim) + Bx(i,j)
																			B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) = B((i-1)*Neq+k,(j-1)*Neq*Ndim+2+(k-1)*Ndim) + By(i,j)																																					
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt_onlyTiTe
								
								!*******************************************
								! Multiply matrix for diffusion
								!*******************************************						
       SUBROUTINE multiply_for_diffusion(mat)
       real*8, intent(inout) :: mat(:,:)
       integer               :: n,i
       
       n = size(mat,1)/Neq
       
       DO i = 1,n
          mat((i-1)*Neq+1,:) = mat((i-1)*Neq+1,:)*phys%diff_n
          mat((i-1)*Neq+2,:) = mat((i-1)*Neq+2,:)*phys%diff_u        
#ifdef TEMPERATURE
          mat((i-1)*Neq+3,:) = mat((i-1)*Neq+3,:)*phys%diff_e     
          mat((i-1)*Neq+4,:) = mat((i-1)*Neq+4,:)*phys%diff_ee     
#endif          
       END DO
       
       END SUBROUTINE multiply_for_diffusion			
       					



								


#ifdef TEMPERATURE

								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av(C11,C21,C22,C33,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22,C33,C44
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j) 
             C((i-1)*Neq+3,(j-1)*Neq+3) = C33(i,j)         
             C((i-1)*Neq+4,(j-1)*Neq+4) = C44(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av			
								

								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)  
             C((i-1)*Neq+3,(j-1)*Neq+3) = Mat(i,j)
             C((i-1)*Neq+4,(j-1)*Neq+4) = Mat(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab(C11,C21,C22,C33,C44,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22,C33,C44
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)*tau(2) 
             C((i-1)*Neq+3,(j-1)*Neq+3) = C33(i,j)*tau(3)         
             C((i-1)*Neq+4,(j-1)*Neq+4) = C44(i,j)*tau(4)
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab			
								




								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)*tau(2)  
             C((i-1)*Neq+3,(j-1)*Neq+3) = Mat(i,j)*tau(3)
             C((i-1)*Neq+4,(j-1)*Neq+4) = Mat(i,j)*tau(4)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab								
								
								


					



								!*****************************************************
								! Expand matrix for Bohm: Temperatures part
								!*****************************************************			
								SUBROUTINE expand_matrix_mass_stab_temperature(A,tau,C)        
								real*8,intent(in),dimension(:,:)      :: A
        real*8,intent(in),dimension(:)        :: tau
								real*8,intent(inout),dimension(size(A,1)*Neq,size(A,2)*Neq) :: C
								integer*4 :: i,j,m,k
								
								m=size(A,1)

        DO i = 1,m
           DO j = 1,m
              DO k =1,2
                 C(Neq*(i-1)+k+2,Neq*(j-1)+k+2) = C(Neq*(i-1)+k+2,Neq*(j-1)+k+2) + tau(k)*A(i,j)
              END DO
           END DO
        END DO
								END SUBROUTINE expand_matrix_mass_stab_temperature
								
        
        
								!*******************************************
								! Expand matrix for Bohm: H matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmH(C31,C32,C33,C34,C41,C42,C43,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C41,C42,C43,C44
								real*8,intent(inout),dimension(size(C31,1)*4,size(C31,2)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(C31,1)
								m=size(C31,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+3,(j-1)*4+1) = C31(i,j)
             C((i-1)*4+3,(j-1)*4+2) = C32(i,j)
             C((i-1)*4+3,(j-1)*4+3) = C33(i,j)
             C((i-1)*4+3,(j-1)*4+4) = C34(i,j)
             C((i-1)*4+4,(j-1)*4+1) = C41(i,j)
             C((i-1)*4+4,(j-1)*4+2) = C42(i,j)
             C((i-1)*4+4,(j-1)*4+3) = C43(i,j)
             C((i-1)*4+4,(j-1)*4+4) = C44(i,j)
             
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmH
								
								
								
								!*******************************************
								! Expand matrix for Bohm: TQ matrix
								!*******************************************											
								SUBROUTINE expand_matrix_BohmTq(C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48,C)        
								real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48
								real*8,intent(inout),dimension(size(C31,1)*neq,size(C31,2)*Neq*Ndim) :: C
								integer*4 :: i,j,n,m
								n=size(C31,1)
								m=size(C31,2)
								DO i=1,n
								   DO j=1,m
								     C((i-1)*4+3,(j-1)*8+1) = C31(i,j)
								     C((i-1)*4+3,(j-1)*8+2) = C32(i,j)
								     C((i-1)*4+3,(j-1)*8+3) = C33(i,j)             
								     C((i-1)*4+3,(j-1)*8+4) = C34(i,j)
								     C((i-1)*4+3,(j-1)*8+5) = C35(i,j)
								     C((i-1)*4+3,(j-1)*8+6) = C36(i,j)
								     C((i-1)*4+3,(j-1)*8+7) = C37(i,j)
								     C((i-1)*4+3,(j-1)*8+8) = C38(i,j)
								     
								     C((i-1)*4+4,(j-1)*8+1) = C41(i,j)
								     C((i-1)*4+4,(j-1)*8+2) = C42(i,j)
								     C((i-1)*4+4,(j-1)*8+3) = C43(i,j)             
								     C((i-1)*4+4,(j-1)*8+4) = C44(i,j)
								     C((i-1)*4+4,(j-1)*8+5) = C45(i,j)
								     C((i-1)*4+4,(j-1)*8+6) = C46(i,j)
								     C((i-1)*4+4,(j-1)*8+7) = C47(i,j)
								     C((i-1)*4+4,(j-1)*8+8) = C48(i,j)								     
								   END DO
								END DO
								END SUBROUTINE expand_matrix_BohmTq				
        
								!*******************************************
								! Expand matrix for Bohm (and Dirichlet) : Tf matrix
								!*******************************************										
								SUBROUTINE expandTf(Tf3,Tf4,Tf)
								real*8, intent(IN) :: Tf3(:),Tf4(:)
								real*8, intent(INOUT):: Tf(:)
								integer            :: i,n								
								n = size(Tf3)
								DO i=1,n
											Tf((i-1)*4+3) = Tf3(i)
											Tf((i-1)*4+4) = Tf4(i)
								END DO				
								END SUBROUTINE expandTf									
#else

								!*******************************************
								! Expand matrix for Bohm: D matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD(C11,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD				
								




								!*******************************************
								! Expand matrix for Bohm: E matrix
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE		
								
								
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av(C11,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j) 
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av			
								


								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av(Mat,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)         
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av
								
								
						
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab(C11,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C22(i,j)*tau(2) 
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab			
								




								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
								real*8,intent(in),dimension(:)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = Mat(i,j)*tau(2)          
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab								
								
#endif						
						
						
								!*******************************************
								! Expand matrix for Bohm: D matrix (all var)
								! with tensor stabilization 
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_av_stab_mat(Ms,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Ms,C21,C22,tau
								real*8,intent(inout),dimension(size(Ms,1)*Neq,size(Ms,2)*Neq) :: C
								integer*4 :: i,j,n,m,k
								
								n=size(Ms,1)
								m=size(Ms,2)

        DO i=1,n
           DO j=1,m
              DO k=1,Neq
						          C((i-1)*Neq+k,(j-1)*Neq+1) = Ms(i,j)*tau(k,1) + C21(i,j)*tau(k,2)
						          C((i-1)*Neq+k,(j-1)*Neq+2) = C22(i,j)*tau(k,2)
						          C((i-1)*Neq+k,(j-1)*Neq+3) = Ms(i,j)*tau(k,3)         
						          C((i-1)*Neq+k,(j-1)*Neq+4) = Ms(i,j)*tau(k,4)
              END DO
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_av_stab_mat	
								
							


								!*******************************************
								! Expand matrix for Bohm: E matrix (all var)
								! with matrix stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_av_stab_mat(A,tau,B)        
										real*8,intent(in)    :: A(:,:),tau(:,:)
										real*8,intent(out)   :: B(1:Neq*size(A,1),1:Neq*size(A,1))
										integer*4            :: i,j,k,t,m,n
										
										n = size(tau,1)
										
										m = size(A,1)
										B = 0.d0
										DO i = 1,m
													DO j = 1,m
													   DO k =1,n
                   DO t = 1,n
													         B(n*(i-1)+k,n*(j-1)+t) = tau(k,t)*A(i,j)
                   END DO
													   END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_BohmE_av_stab_mat
								
								
								
								!*******************************************
								! Expand matrix for Bohm: D matrix - 
        ! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmD_stab(C11,C21,C22,tau,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C21,C22
        real*8,intent(in),dimension(2)          :: tau
								real*8,intent(inout),dimension(size(C11,1)*Neq,size(C11,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C((i-1)*Neq+1,(j-1)*Neq+1) + C11(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+1) = C((i-1)*Neq+2,(j-1)*Neq+1) + C21(i,j)*tau(2)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C((i-1)*Neq+2,(j-1)*Neq+2) + C22(i,j)*tau(2)   
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmD_stab				
								




								!*******************************************
								! Expand matrix for Bohm: E matrix - 
        ! with stabilization
								!*******************************************									
								SUBROUTINE expand_matrix_BohmE_stab(Mat,tau,C)        
								real*8,intent(in),dimension(:,:)        :: Mat
        real*8,intent(in),dimension(2)          :: tau
								real*8,intent(inout),dimension(size(Mat,1)*Neq,size(Mat,2)*Neq) :: C
								integer*4 :: i,j,n,m
								
								n=size(Mat,1)
								m=size(Mat,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*Neq+1,(j-1)*Neq+1) = C((i-1)*Neq+1,(j-1)*Neq+1) + Mat(i,j)*tau(1)
             C((i-1)*Neq+2,(j-1)*Neq+2) = C((i-1)*Neq+2,(j-1)*Neq+2) + Mat(i,j)*tau(2)   
           END DO
        END DO
								END SUBROUTINE expand_matrix_BohmE_stab		
        

								

END SUBROUTINE HDG_BC


#endif  
  
