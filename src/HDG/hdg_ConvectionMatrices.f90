!*****************************************
! project: MHDG
! file: hdg_ConvectionMatrices.f90
! date: 20/12/2016
! Generate the matrices that change
! during the NR iterative process
!*****************************************


SUBROUTINE HDG_convectionmatrices()
  USE globals
  USE LinearAlgebra
  USE analytical
  USE physics
  USE printUtils
  USE MPI_OMP
  USE Debug
  
  IMPLICIT NONE
    
#ifdef TOR3D   
!***********************************************************************
! 
!                            VERSION 3D TOROIDAL
! 
!***********************************************************************

  integer*4             :: itor,iel,iel3,i,j,Ndim,Neq,N2D,Np,Nfp,Nf,Nfdir,dd
  integer*4             :: ntorloc,itorg,Np1Dpol,Np1Dtor,Np2D,Nfl
  integer*4             :: ind_ue(refElTor%Nnodes3D*phys%Neq)
  integer*4             :: ind_uf(refElTor%Nft*phys%Neq)
  real*8                :: tdiv(numer%ntor+1)
  real*8                :: htor,tel(refElTor%Nnodes1d)
  integer*4             :: F_el(refElPol%Nfaces)
  integer*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8,allocatable    :: ue(:),uf(:)
  real*8,allocatable    :: Cvel(:,:),Hel(:,:),Hfel(:,:)
  ! Matrices to be filled in case of non-constant stabilization
  real*8,allocatable    :: Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:)
  
  integer               :: delta(refElPol%Nfaces+2)
  integer               :: ind_dim(refElPol%Nfaces+2),ind_sta(refElPol%Nfaces+2),aux,ifa  

		IF (utils%printint>1) THEN
   WRITE(6,*) '*************************************************'
   WRITE(6,*) '*           CONVECTION MATRICES                 *'
   WRITE(6,*) '*************************************************' 
		END IF	   
		
  
  Ndim        = 3                                               ! Number of dimensions
  N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Neq         = phys%Neq                                        ! Number of equations
  Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np          = refElTor%Nnodes3D                               ! Number of nodes for each 3D element
  Nfl         = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfp         = Np2D*2+refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
  Nf          = refElPol%Nfaces                                 ! Number of faces in the 2D mesh
  Nfdir       = Mesh%Ndir


  ! Local indices
  ind_loc = 0
  DO i = 1,Nf
     DO j = 1,Nfl*Neq
        ind_loc(i,j) = refElPol%Nnodes2D*Neq + Neq*Nfl*(i-1) + j
     END DO
  END DO

  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i=1,numer%ntor
     tdiv(i+1) = i*htor
  END DO
  
  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol,Np1Dtor,Neq,perm)  
  
  ! Indices to recovery the faces solution
  ind_uf = 0
  ind_dim(1) =  Np2D*Neq
  ind_dim(2: refElPol%Nfaces+1) =  Nfl*Neq 
  ind_dim(refElPol%Nfaces+2) = Np2D*Neq
  ind_sta(1) = 1
  aux = 0
  DO i= 2,refElPol%Nfaces+2
     ind_sta(i) = 1+ind_dim(i-1)+aux
     aux = aux+ind_dim(i-1)
  END DO
  
#ifdef PARALL
    IF (MPIvar%ntor.gt.1) THEN
       ntorloc = numer%ntor/MPIvar%ntor+1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
  
!$OMP PARALLEL DEFAULT(SHARED) & 
!$OMP PRIVATE(itor,itorg,iel,iel3,ifa,ind_ue,ind_uf,Xe,tel,F_dir_el,F_el,i,j,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel,dd,delta) 
  ALLOCATE(Cvel(Neq*Np,Neq*Np))
  ALLOCATE(Hel(Neq*Np,Neq*Nfp))
  ALLOCATE(Hfel(Neq*Nfp,Neq*Nfp))
  ALLOCATE(Del(Neq*Np,Neq*Np))
  ALLOCATE(Eel(Neq*Np,Neq*Nfp))
  ALLOCATE(Dfel(Neq*Nfp,Neq*Np))
  ALLOCATE(Efel(Neq*Nfp,Neq*Nfp))
  ALLOCATE(ue(Neq*Np))
  ALLOCATE(uf(Nfp*Neq))
!$OMP DO SCHEDULE(STATIC) 
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif 
     tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg+1)-tdiv(itorg))
					DO iel = 1,N2D
						  Cvel    = 0.d0
						  Hel     = 0.d0
						  Hfel    = 0.d0

						  Del    = 0.d0
						  Eel    = 0.d0
						  Dfel   = 0.d0
						  Efel   = 0.d0
						  
						  ! Index of 3D element
						  iel3 = (itor-1)*N2d+iel
						  						  
						  ! Coordinates of the nodes of the element
						  Xe = Mesh%X(Mesh%T(iel,:),:)

						  ! Dirichlet boundary conditions in the current element
						  F_dir_el = Mesh%Fdir(iel,:)
						  
						  ! Indices to extract the elemental and face solution
						  F_el = Mesh%F(iel,:)     
						  ind_ue = (iel3-1)*Np*Neq + (/(i,i=1,Neq*Np)/)
        dd = 1+(itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Nfl)*Neq
        delta(1) = dd +(iel-1)*Np2D*Neq
        delta(2:refElPol%Nfaces+1) = dd + (N2D*Np2D+(F_el-1)*Nfl)*Neq
        delta(refElPol%Nfaces+2) = dd + (N2D*Np2D+(Mesh%Nfaces-Nfdir)*Nfl+(iel-1)*Np2D)*Neq
!        if (MPIvar%ntor.eq.1) then
        if (itor==numer%ntor) then
           delta(refElPol%Nfaces+2) = 1+(iel-1)*Np2D*Neq
        end if
!        endif

						  ! Elemental and face solution     
						  ue = sol%u(ind_ue)
        DO ifa = 1,Nf+2
            DO i=0,ind_dim(ifa)-1 
               ind_uf(i+ind_sta(ifa)) = delta(ifa)+i
            END DO
        END DO  	
        uf = sol%u_tilde(ind_uf)					  

						  ! Compute the matrices for the element
						  CALL elemental_matrices(itor,iel,Xe,tel,F_dir_el,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel)

						  
						  ! Flip the faces for the elements for which the face is already been counted
						  DO j = 1,Nf
						     if (Mesh%flipface(iel,j)) then
						        Hfel(ind_loc(j,:),:) = Hfel(ind_loc(j,perm(:)),:)
						        Hfel(:,ind_loc(j,:)) = Hfel(:,ind_loc(j,perm(:)))
						        Dfel(ind_loc(j,:),:) = Dfel(ind_loc(j,perm(:)),:)
						        Efel(ind_loc(j,:),:) = Efel(ind_loc(j,perm(:)),:)
						        Efel(:,ind_loc(j,:)) = Efel(:,ind_loc(j,perm(:)))           
						     END if
						  END DO

						  ! Store the matrices for all the elements
						  elMat%Cv(:,:,iel3)  = Cvel
						  elMat%H(:,:,iel3)   = Hel
						  elMat%Hf(:,:,iel3)  = Hfel   
						  IF (numer%stab > 1) THEN
						     elMat%D(:,:,iel3)  = Del
						     elMat%E(:,:,iel3)  = Eel
						     elMat%Df(:,:,iel3) = Dfel
						     elMat%Ef(:,:,iel3) = Efel     
						  END IF
					END DO
  END DO
!$OMP END DO

 DEALLOCATE(Cvel,Hel,Hfel)
 DEALLOCATE(Del,Eel,Dfel,Efel,ue,uf)  
!$OMP END PARALLEL 

   
  IF (utils%printint>1) THEN
     WRITE(6,*) "Done!"
  END IF
  
  CONTAINS
			
			
!***************************************************
! Elemental matrices computation: 3D case, 2-4 eqs
!***************************************************
					SUBROUTINE elemental_matrices(itor,iel,Xe,tel,F_dir_el,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel)
       integer,intent(IN)			     :: iel,itor 
					  real*8,intent(IN)         :: Xe(:,:),tel(:)
					  logical,intent(IN)        :: F_dir_el(:)
       real*8,intent(IN)         :: ue(:),uf(:)
					  real*8,intent(INOUT)      :: Cvel(:,:),Hel(:,:),Hfel(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:)
					  
					  ! Dimension and indices
       integer*4                  :: NGaussTor,NgaussPol
							integer*4                  :: igtor,igpol,g,ifa,i,j
							integer                    :: ind(refElPol%Ngauss1d),nodes2(refelPol%Nnodes1d),ind2(refElPol%Ngauss2d)
							integer*4,allocatable      :: ind_ass(:),ind_asf(:)
							! Volume computations
							real*8                     :: dvolu,isdir,isext,htor,dvolu1d
							real*8                     :: xy(refElPol%Ngauss2d,2),teg(refElTor%Ngauss1d)
       real*8                     :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
       real*8                     :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
       real*8                     :: detJ(refElPol%Ngauss2d)
       real*8                     :: iJ11(refElPol%Ngauss2d),iJ12(refElPol%Ngauss2d)
       real*8                     :: iJ21(refElPol%Ngauss2d),iJ22(refElPol%Ngauss2d)
       real*8                     :: Nxg(refElPol%Nnodes2D),Nyg(refElPol%Nnodes2D),Nx_ax(refElPol%Nnodes2D)
       real*8                     :: N1g(refElTor%Nnodes1D),N1xg_cart(refElTor%Nnodes1D),N1xg(refElTor%Nnodes1D),N2g(refElPol%Nnodes2D)
       real*8,allocatable         :: Ni(:),Nidvolu(:),Nr(:),Nz(:),Nt(:),auxNNxy(:,:),NiNi(:,:)
       ! Faces computations
       integer*4,allocatable      :: ind_ff(:),ind_ffd(:),ind_fe(:),ind_fg(:)
       real*8,allocatable         :: aux(:,:),xyf(:,:)
       real*8,allocatable         :: thetafg(:),dsurf(:),n_g(:,:)       
							real*8                     :: xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
							real*8                     :: dsurfg
							real*8                     :: t_g(refElPol%Ngauss1d,2)	
							real*8,pointer             :: Nfi(:,:),Nfg(:)
							integer                    :: Npifa
       ! Magnetic field, volume computations							       
       real*8                     :: b(refElPol%Ngauss2d*refElTor%Ngauss1d,3),b2d(refElPol%Ngauss2d,3)
       real*8                     :: divb(refElPol%Ngauss2d*refElTor%Ngauss1d),divb2d(refElPol%Ngauss2d)
       real*8                     :: drift(refElPol%Ngauss2d*refElTor%Ngauss1d,3),drift2d(refElPol%Ngauss2d,3)
       real*8                     :: divbg, driftg(2)
       ! Magnetic field, faces computations
       real*8,allocatable         :: b_f(:,:)
       real*8                     :: b_f2d(refElPol%Ngauss1d,3),bn
       ! Solution in volume
       real*8                     :: ueg(refElTor%NGauss3D,neq),upg(refElPol%Ngauss2d*refElTor%Ngauss1d,phys%npv)
       ! Solution in faces
       real*8,allocatable         :: ufi(:),ufg(:,:),upgf(:,:)
       ! Routine specific part
       real*8,dimension(neq,neq)  :: A
							real*8                     :: Ax(phys%Neq,phys%Neq),Ay(phys%Neq,phys%Neq),An(phys%Neq,phys%Neq)
							real*8                     :: tau(phys%Neq)
       real*8                     :: taumat(phys%Neq,phys%Neq)
							real*8,allocatable         :: NxNi(:,:),NyNi(:,:),NtNi(:,:),NNxy(:,:)
       real*8,allocatable,dimension(:,:)   :: Hloc
#ifdef TEMPERATURE
       real*8                              :: Telect
       real*8,allocatable,dimension(:,:)   :: NNi
							real*8,dimension(phys%Neq,phys%Neq) :: GG
							real*8,allocatable,dimension(:,:)   :: Dre,Drifte
#endif														
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(refElPol%Nnodes2D,3),Bmod(refElPol%Nnodes2D)
       real*8                    :: Bmod_g(refElPol%Ngauss2d),Bmod_x,Bmod_y,Nx_ax(Np)
       real*8                    :: Bt(refElPol%Nnodes2D),Bt_g(refElPol%Ngauss2d)              
#endif
       real*8,allocatable        :: Mf_loc(:,:),Massf(:,:)
 						real*8, parameter         :: tol = 1e-12      
 						
       
       ALLOCATE(ind_ass(Np))
       ind_ass = (/ (i, i = 0, Neq*(Np-1),Neq) /)
       
              
       !***********************************
       !    Volume computation
       !***********************************
       ALLOCATE(Ni(Np))
							ALLOCATE(Nr(Np))
							ALLOCATE(Nz(Np))
							ALLOCATE(Nt(Np))
							ALLOCATE(Nidvolu(Np))
							ALLOCATE(NiNi(Np,Np))							
							ALLOCATE(auxNNxy(Ndim,Np)) 
							ALLOCATE(NxNi(Np,Np))
							ALLOCATE(NyNi(Np,Np))
							ALLOCATE(NtNi(Np,Np))
							ALLOCATE(NNxy(Np,Np))
														
#ifdef TEMPERATURE
       ALLOCATE(Dre(Np,Np))
       ALLOCATE(Drifte(Np*Neq,Np*Neq))
       Drifte = 0.
       Dre = 0.
#endif	
							
							
						 ! toroidal element size
						 htor = tel(refElTor%Nnodes1d) - tel(1)
						 
							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)
							
							! Solution at Gauss points
							ueg = matmul(refElTor%N3D,transpose(reshape(ue,[neq,Np])))
							
							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)
							
#ifndef MOVINGEQUILIBRIUM					
!************ CASE OF STATIC EQUILIBRIUM ****************************
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b2d = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
           divb2d = matmul(refElPol%N2D,phys%divb(Mesh%T(iel,:)))
           drift2d = matmul(refElPol%N2D,phys%drift(Mesh%T(iel,:),:))
           DO igtor = 1,refElTor%NGauss1D
              ind2 = (igtor-1)*refElPol%NGauss2D + (/(j,j=1,refElPol%NGauss2D)/)
              b(ind2,:) = b2d
              divb(ind2)   = divb2d
              drift(ind2,:) = drift2d
           END DO           
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),teg,b,divb,drift)
							END IF
#ifdef TEMPERATURE
       IF (switch%driftvel) THEN
						    ! Physical variables at Gauss points
						    CALL cons2phys(ueg,upg)
						    
						    ! Drift at Gauss points
						    drift = phys%dfcoef*drift ! I multiply for the temperature later
       END IF
#endif							
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
       b_nod = 0. ! Nnodes
       Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2) ! Nnodes
       b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod ! Nnodes
       b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod ! Nnodes
       b_nod(:,3) = phys%Bt(Mesh%T(iel,:))/Bmod ! Nnodes
       b2d = matmul(refElPol%N2D,b_nod)         ! Ngauss
       DO igtor = 1,refElTor%NGauss1D
          ind2 = (igtor-1)*refElPol%NGauss2D + (/(j,j=1,refElPol%NGauss2D)/)
          b(ind2,:) = b2d
       END DO       
#ifdef TEMPERATURE        
       IF (switch%driftvel) THEN
          ! Physical variables at Gauss points
          CALL cons2phys(ueg,upg)  ! Ngauss
        
          ! Drift at Gauss points       
          Bmod_g = matmul(refElPol%N2D,Bmod) ! Ngauss
          Bt_g   = matmul(refElPol%N2D,phys%Bt(Mesh%T(iel,:)))  ! Ngauss 
       END IF       
#endif          
#endif		
!************ END SELECTION MOVING EQUILIBRIUM **********************

					  !***************** VOLUME COMPUTATION *******************
        J11 = matmul(refElPol%Nxi2D,Xe(:,1))                           ! ng x 1
        J12 = matmul(refElPol%Nxi2D,Xe(:,2))                           ! ng x 1
        J21 = matmul(refElPol%Neta2D,Xe(:,1))                          ! ng x 1
        J22 = matmul(refElPol%Neta2D,Xe(:,2))                          ! ng x 1
        detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
					   iJ11 =  J22/detJ
					   iJ12 = -J12/detJ
					   iJ21 = -J21/detJ
					   iJ22 =  J11/detJ
        NgaussPol = refElPol%NGauss2D
        NgaussTor = refElTor%NGauss1D
		 					DO igtor = 1,NGaussTor
			 				   N1g    = refElTor%N1D(igtor,:)         ! Toroidal shape function 
				 			   N1xg_cart   = refElTor%Nxi1D(igtor,:)*2/htor       ! Toroidal shape function derivative
					 		   dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume
					 		   
						 	   DO igpol = 1,NGaussPol
							       g = (igtor-1)*NGaussPol+igpol   
							       
									 	   ! Poloidal shape functions and derivatives    							  
    						 	  N2g = refElPol%N2D(igpol,:)										   
							       Nxg = iJ11(igpol)*refElPol%Nxi2D(igpol,:) + iJ12(igpol)*refElPol%Neta2D(igpol,:)
							       Nyg = iJ21(igpol)*refElPol%Nxi2D(igpol,:) + iJ22(igpol)*refElPol%Neta2D(igpol,:)

										    ! 3D integration weight
										    dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d
										    
              IF (switch%axisym) THEN										
                 dvolu = dvolu*xy(igpol,1)
                 N1xg = N1xg_cart/xy(igpol,1)
              ELSE
              			N1xg = N1xg_cart
              END IF

							       ! 3D shape functions 
							       Ni  = col(TensorProduct(N2g,N1g))      ! 3D shape function
							       Nr  = col(TensorProduct(Nxg,N1g))      ! 3D shape function, derivative in r for divergence
							       Nz  = col(TensorProduct(Nyg,N1g))      ! 3D shape function, derivative in z
							       Nt  = col(TensorProduct(N2g,N1xg))     ! 3D shape function, derivative in t
										    
	
							       ! Shape functions tensor product
							       Nidvolu      = Ni*dvolu			
							       
              ! Jacobian matrices
               CALL jacobianMatrices(ueg(g,:),A)   

									    ! Contribution of the current integration point to the elemental matrix
									     NxNi = tensorProduct(Nr,Nidvolu)
              NyNi = tensorProduct(Nz,Nidvolu)
              NtNi = tensorProduct(Nt,Nidvolu)
              NNxy = b(g,1)*NxNi+b(g,2)*NyNi+b(g,3)*NtNi									
           
              ! Assembly convection matrix
              DO i=1,Neq
                 DO j=1,Neq
                    Cvel(i+ind_ass,j+ind_ass) = Cvel(i+ind_ass,j+ind_ass) + A(i,j)*NNxy
                 END DO
              END DO

          
#ifdef TEMPERATURE
              ! Divergence of b at the Gauss points  
#ifndef MOVINGEQUILIBRIUM
              divbg = divb(g)
#else      
							       IF (switch%axisym) THEN
							          Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
							       ELSE
							          Nx_ax = Nxg
							       END IF  
              divbg = dot_product(Nx_ax,b_nod(:,1))+dot_product(Nyg,b_nod(:,2))                
#endif          

              ! Matrix for the curvature term
              CALL GimpMatrix(ueg(g,:),divbg,GG)   
              
              ! Assembly curvature term in the convection matrix          
              NNi = tensorProduct(Ni,Nidvolu)
              DO i=2,2
                 DO j=1,Neq
                    Cvel(i+ind_ass,j+ind_ass) = Cvel(i+ind_ass,j+ind_ass) + GG(i,j)*NNi
                 END DO
              END DO
              IF (switch%driftvel) THEN               
#ifdef MOVINGEQUILIBRIUM          
			           Bmod_x  = dot_product(Nxg,Bmod) ! 1
			           Bmod_y  = dot_product(Nyg,Bmod) ! 1
			           driftg(1) =  phys%dfcoef*Bt_g(igpol)*Bmod_y/Bmod_g(igpol)**3 ! I multiply for the temperature later
			           driftg(2) = -phys%dfcoef*Bt_g(igpol)*Bmod_x/Bmod_g(igpol)**3 ! I multiply for the temperature later          
#else
              driftg(1) = drift(g,1)
              driftg(2) = drift(g,2)
#endif			       
              Telect = upg(g,8)
              Dre = Dre + TensorProduct(Nidvolu,Telect*( Nr*driftg(1)+ Nz*driftg(2) ))
           END IF
#endif  

          END DO
							END DO
							
							
							
							
							
							
#ifdef TEMPERATURE							
						CALL expand_matrix_mass(Dre,Neq,Drifte)
      Cvel = Cvel - Drifte
      DEALLOCATE(Drifte,Dre)
#endif



      DEALLOCATE(NiNi,NxNi,NyNi,NtNi,Ni,Nidvolu,Nr,Nz,Nt,NNxy,auxNNxy)
      DEALLOCATE(ind_ass)
       
       !***********************************
  					! Faces computations
       !***********************************

												
       ! Loop in the faces of the element
							DO ifa = 1,Nf+2
							
							
							  isext = 0.
							  
          IF (ifa==1) THEN
							      Npifa = Np2D
							      NGaussPol = refElPol%Ngauss2D
							      NGaussTor = 1
							      ALLOCATE(ind_ff(Npifa*Neq))
							      ALLOCATE(ind_fe(Npifa*Neq))
							      ALLOCATE(ind_fg(Ndim*Npifa*Neq))
							      ALLOCATE(xyf(NGaussPol,2))
							      ALLOCATE(thetafg(1))
							      ALLOCATE(dsurf(NGaussPol*NGaussTor))
							      ALLOCATE(n_g(NGaussPol*NGaussTor,3))
							      ind_ff = (/(i,i=1,Npifa*Neq)/)
							      ind_fe = ind_ff
							      ind_fg = (/(i,i=1,Npifa*Ndim*Neq)/)
							      xyf    = matmul(refElPol%N2D,Xe)
							      dsurf  =  refElPol%gauss_weights2D*detJ
							      isdir  = 0.
							      thetafg = tel(1)
							      Nfi => refElPol%N2D
							      n_g = 0.; n_g(:,3) = -1
							   
							   ELSE IF (ifa==refElPol%Nfaces+2) THEN
							      Npifa = Np2D
							      NGaussPol = refElPol%Ngauss2D
							      NGaussTor = 1
								     ALLOCATE(ind_ff(Npifa*Neq)) 
							      ALLOCATE(ind_fe(Npifa*Neq))
							      ALLOCATE(ind_fg(Ndim*Npifa*Neq))	
							      ALLOCATE(xyf(NGaussPol,2))	
							      ALLOCATE(thetafg(1))
							      ALLOCATE(dsurf(NGaussPol*NGaussTor))
							      ALLOCATE(n_g(NGaussPol*NGaussTor,3))
							      ind_ff = Np2D*Neq + refElPol%Nfaces*Nfl*Neq+(/(i,i=1,Npifa*Neq)/)
							      ind_fe = Np2d*Neq*(Np1dTor-1)+(/(i,i=1,Npifa*Neq)/)
							      ind_fg = Np2d*Neq*(Np1dTor-1)*Ndim + (/(i,i=1,Npifa*Ndim*Neq)/) 
							      xyf    = matmul(refElPol%N2D,Xe)				
							      thetafg = tel(refElTor%Nnodes1d) 			   
							      dsurf  =  refElPol%gauss_weights2D*detJ
							      isdir  = 0.
							      Nfi => refElPol%N2D
							      n_g = 0.; n_g(:,3) = 1
							   ELSE
							      Npifa = Nfl 
							      nodes2  = refElPol%face_nodes(ifa-1,:)
							      NGaussPol = refElPol%Ngauss1D
							      NGaussTor = refElTor%Ngauss1D							      
							      ALLOCATE(ind_ff(Npifa*Neq))
							      ALLOCATE(ind_fe(Npifa*Neq))
							      ALLOCATE(ind_fg(Ndim*Npifa*Neq))				 
							      ALLOCATE(xyf(NGaussPol,2))
							      ALLOCATE(thetafg(NGaussTor))
							      ALLOCATE(dsurf(NGaussPol*NGaussTor))
							      ALLOCATE(n_g(NGaussPol*NGaussTor,3))
             ind_ff = Np2d*Neq + (ifa-2)*Npifa*Neq+ (/(i,i=1,Npifa*Neq)/)
             ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*Neq)) 
             ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*3*Neq))
							      Nfi => refElTor%sFTF	      
							      xyf    = matmul(refElPol%N1D,Xe(nodes2,:))											      
							      thetafg = matmul(refElTor%N1D,tel)
							      
             xyd_g     = matmul(refelPol%Nxi1D,Xe(nodes2,:))
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
             IF (F_dir_el(ifa-1)) THEN
                isdir = 1.
             ELSE
                isdir = 0.
             ENDIF
             

                       
#ifdef PARALL          
! probably useless, to be verified
						       IF ( (Mesh%F(iel,ifa-1).gt.Mesh%Nintfaces)  ) THEN
						          IF (Mesh%boundaryFlag(Mesh%F(iel,ifa-1)-Mesh%Nintfaces).ne.0) THEN
						             isext = 1.
						          END IF
						       ENDIF 
#else
			          IF (Mesh%F(iel,ifa-1) > Mesh%Nintfaces ) isext = 1.
#endif               
							   END IF
							   
							   
          ! Allocation of matrices
							   ALLOCATE(NiNi(Npifa,Npifa))
							   ALLOCATE(Hloc(Neq*Npifa,Neq*Npifa))
										IF (numer%stab>1) THEN
							      ALLOCATE(Mf_loc(Neq*Nfp,Neq*Nfp))
							      Mf_loc = 0.d0							
										END IF				  
										ALLOCATE(b_f(NGaussPol*NGaussTor,3))	 
										ALLOCATE(ufi(Npifa*Neq))
										ALLOCATE(ufg(NGaussPol*NGaussTor,Neq))
										ALLOCATE(ind_asf(Npifa))
							   ind_asf = (/ (i, i = 0, Neq*(Npifa-1),Neq) /)
							   Hloc = 0.
							   
	              
          ! Set the face solution
          IF (abs(isdir-1.)<1e-12) THEN 
             CALL analytical_solution(xyf(:,1),xyf(:,2),thetafg,ufg)
          ELSE
										   ufi = uf(ind_ff)
										   IF (ifa>1 .and. ifa<refElPol%Nfaces+2) THEN
										      IF  (Mesh%flipface(iel,ifa-1)) THEN
														     ufi = ufi(perm)
														  END IF
										   END IF             
             ufg = matmul(Nfi,transpose(reshape(ufi,[neq,Npifa])))
          END IF 	                       

          
          IF (numer%stab>1) THEN
               ALLOCATE(upgf(NGaussPol*NGaussTor,phys%npv))
               CALL cons2phys(ufg,upgf) 
          END IF

#ifndef MOVINGEQUILIBRIUM			
!************ CASE OF STATIC EQUILIBRIUM ****************************
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
						       IF (ifa==1 .or. ifa==refElPol%Nfaces+2 ) THEN
						          b_f = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
						       ELSE
						          b_f2d = matmul(refElPol%N1D,phys%b(Mesh%T(iel,refElPol%face_nodes(ifa-1,:)),:))
                DO igtor = 1,refElTor%NGauss1D
                   ind = (igtor-1)*refElPol%NGauss1D + (/(j,j=1,refElPol%NGauss1D)/)
                   b_f(ind,:)     = b_f2d
                END DO				
             END IF		        
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),thetafg,b_f)   
										END IF				
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
          b_f = matmul(refElPol%N1D,b_nod(refElPol%face_nodes(ifa,:),:))
#endif
!************ END SELECTION MOVING EQUILIBRIUM **********************
										
         !*****************************
         ! Loop in face Gauss points
         !*****************************
										DO igtor = 1,NGaussTor
										   DO igpol = 1,NGaussPol
										   										      
										      g = (igtor-1)*NGaussPol+igpol 

										      ! Face shape functions										      
										      Nfg => Nfi(g,:)
										
												    IF (switch%axisym .and. ifa.gt.1 .and. ifa.lt.refElPol%Nfaces+2 ) THEN			
!												    IF (switch%axisym ) THEN									
												       dsurfg = dsurf(g)*xyf(igpol,1)
												    ELSE
												       dsurfg = dsurf(g)
												    END IF							
										   									   
			             ! Jacobian matrices
			             bn = dot_product(b_f(g,:),n_g(g,:))
                CALL jacobianMatrices(ufg(g,:),A)  

										      ! Contribution of the current integration point to the elemental matrix
                NiNi = tensorProduct(Nfg,Nfg)*dsurfg
													   
                ! Assembly convection matrix
                DO i=1,Neq
                   DO j=1,Neq
                      Hloc(i+ind_asf,j+ind_asf) = Hloc(i+ind_asf,j+ind_asf) + bn*A(i,j)*NiNi
                   END DO
                END DO

             		!************* Non constant stabilization************************											
               IF (numer%stab>1) THEN
                  ! Compute tau in the Gauss points
                  IF (numer%stab<5) THEN
                     CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),iel,ifa,isext,xyf(g,:),tau)
                  ELSE
                     taumat = 0.
                     CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),isext,iel,taumat)
                  ENDIF
                  
										     ! Expand
               IF (numer%stab<5) THEN
                  DO i=1,Neq
                     Mf_loc(i+ind_asf,i+ind_asf) = Mf_loc(i+ind_asf,i+ind_asf) + tau(i)*NiNi
                  END DO 										        
               ELSE
                  DO i=1,Neq
                     DO j=1,Neq
                        Mf_loc(i+ind_asf,j+ind_asf) = Mf_loc(i+ind_asf,j+ind_asf) + taumat(i,j)*NiNi
                     END DO
                  END DO 
               ENDIF
               END IF					       
               !************* End non constant stabilization************************      
			
   										END DO ! Gauss points
										END DO ! Gauss points
										DEALLOCATE(b_f,n_g,dsurf,xyf,thetafg,NiNi,ufi,ufg)
										DEALLOCATE(ind_asf)
          IF (numer%stab>1) THEN
               DEALLOCATE(upgf)
          END IF										

										! Elemental assembly
										Hel(ind_fe,ind_ff)  = Hel(ind_fe,ind_ff) + Hloc*(1-isdir)															
										Hfel(ind_ff,ind_ff) = Hfel(ind_ff,ind_ff) + Hloc*(1-isext)
										DEALLOCATE(Hloc)
										
										 IF (numer%stab>1) THEN 
										    ! Elemental assembly
										    Del(ind_fe,ind_fe)  = Del(ind_fe,ind_fe) + Mf_loc										
										    Eel(ind_fe,ind_ff)  = Eel(ind_fe,ind_ff) + Mf_loc*(1-isdir)
										    Efel(ind_ff,ind_ff) = Efel(ind_ff,ind_ff) + Mf_loc*(1-isext)
										    Dfel(ind_ff,ind_fe) = Dfel(ind_ff,ind_fe) + Mf_loc*(1-isext)	
										    DEALLOCATE(Mf_loc)
										 END IF									
										 DEALLOCATE(ind_ff,ind_fe,ind_fg)
							END DO ! faces
							
							

					END SUBROUTINE elemental_matrices				
													
!*******************************************
!           AUXILIARY ROUTINES
!*******************************************				
								!*******************************************
								! Expand matrix mass type
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

#ifdef TEMPERATURE									
								!*******************************************
								! Expand matrix C: 4 equations
								!*******************************************									
								SUBROUTINE expand_matrix(C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44
								real*8,intent(inout),dimension(size(C11,1)*4,size(C11,1)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+1,(j-1)*4+1) = C((i-1)*4+1,(j-1)*4+1) + C11(i,j)
             C((i-1)*4+1,(j-1)*4+2) = C((i-1)*4+1,(j-1)*4+2) + C12(i,j)
             C((i-1)*4+1,(j-1)*4+3) = C((i-1)*4+1,(j-1)*4+3) + C13(i,j)
             C((i-1)*4+1,(j-1)*4+4) = C((i-1)*4+1,(j-1)*4+4) + C14(i,j)
             
             C((i-1)*4+2,(j-1)*4+1) = C((i-1)*4+2,(j-1)*4+1) + C21(i,j)
             C((i-1)*4+2,(j-1)*4+2) = C((i-1)*4+2,(j-1)*4+2) + C22(i,j)
             C((i-1)*4+2,(j-1)*4+3) = C((i-1)*4+2,(j-1)*4+3) + C23(i,j)
             C((i-1)*4+2,(j-1)*4+4) = C((i-1)*4+2,(j-1)*4+4) + C24(i,j)
                          
             C((i-1)*4+3,(j-1)*4+1) = C((i-1)*4+3,(j-1)*4+1) + C31(i,j)
             C((i-1)*4+3,(j-1)*4+2) = C((i-1)*4+3,(j-1)*4+2) + C32(i,j)
             C((i-1)*4+3,(j-1)*4+3) = C((i-1)*4+3,(j-1)*4+3) + C33(i,j)             
             C((i-1)*4+3,(j-1)*4+4) = C((i-1)*4+3,(j-1)*4+4) + C34(i,j)             

             C((i-1)*4+4,(j-1)*4+1) = C((i-1)*4+4,(j-1)*4+1) + C41(i,j)
             C((i-1)*4+4,(j-1)*4+2) = C((i-1)*4+4,(j-1)*4+2) + C42(i,j)
             C((i-1)*4+4,(j-1)*4+3) = C((i-1)*4+4,(j-1)*4+3) + C43(i,j)             
             C((i-1)*4+4,(j-1)*4+4) = C((i-1)*4+4,(j-1)*4+4) + C44(i,j)             
           END DO
        END DO
								END SUBROUTINE expand_matrix
#else								
								!*******************************************
								! Expand matrix C: 2 equations
								!*******************************************									
								SUBROUTINE expand_matrix(C11,C12,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*2,size(C11,1)*2) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*2+1,(j-1)*2+1) = C((i-1)*2+1,(j-1)*2+1) + C11(i,j)
             C((i-1)*2+1,(j-1)*2+2) = C((i-1)*2+1,(j-1)*2+2) + C12(i,j)
             C((i-1)*2+2,(j-1)*2+1) = C((i-1)*2+2,(j-1)*2+1) + C21(i,j)
             C((i-1)*2+2,(j-1)*2+2) = C((i-1)*2+2,(j-1)*2+2) + C22(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix
#endif

#ifdef TEMPERATURE
							
								SUBROUTINE expand_matrix_Gimp(G21,G22,G23,G24,C)        
								real*8,intent(in),dimension(:,:)        :: G21,G22,G23,G24
								real*8,intent(inout),dimension(size(G21,1)*4,size(G21,1)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(G21,1)
								m=size(G21,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+2,(j-1)*4+1) = C((i-1)*4+2,(j-1)*4+1) + G21(i,j)
             C((i-1)*4+2,(j-1)*4+2) = C((i-1)*4+2,(j-1)*4+2) + G22(i,j)
             C((i-1)*4+2,(j-1)*4+3) = C((i-1)*4+2,(j-1)*4+3) + G23(i,j) 
             C((i-1)*4+2,(j-1)*4+4) = C((i-1)*4+2,(j-1)*4+4) + G24(i,j)            
           END DO
        END DO
								END SUBROUTINE expand_matrix_Gimp
#endif
       
							!*****************************************
							! Set permutations for flipping faces
							!****************************************     
							 SUBROUTINE set_permutations(Np1Dpol,Np1Dtor,Neq,perm)
							 integer, intent(IN)  :: Np1Dpol,Np1Dtor,Neq
							 integer, intent(OUT) :: perm(:)
							 integer              :: i,j,k
							 integer              :: temp_pol(Np1Dpol),temp_tor(Np1Dtor),aux(Np1Dpol,Np1Dtor)

							 temp_pol = (/ (i, i = 1, Np1Dpol) /)
							 temp_tor = Np1Dpol*( (/ (i, i = 1, Np1Dtor) /)-1)
							 aux = TensorSumInt(temp_pol,temp_tor)
			     DO j = 1, Np1Dtor
						     DO i = 1, Np1Dpol
						        DO k = 1, Neq
   						        perm( (j-1)*Np1Dpol*Neq+(i-1)*Neq+k ) = (aux(Np1Dpol-i+1,j)-1)*Neq+k
						        END DO
						     END DO
			     END DO
							 END SUBROUTINE set_permutations       
							 
							 
								!*******************************************
								! Expand matrix with different stabilization
								! parameters
								!*******************************************									
								SUBROUTINE expand_matrix_stab(A,n,tau,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
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
								END SUBROUTINE expand_matrix_stab	
								
        
        
 							!*******************************************
								! Expand matrix to use a tensor stabilization
        ! parameter
								!*******************************************			               
								SUBROUTINE expand_matrix_stab_mat(A,n,tau,B)        
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:,:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
										integer*4            :: i,j,k,t,m
										
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
								END SUBROUTINE expand_matrix_stab_mat
								

							       
END SUBROUTINE hdg_convectionmatrices

#else
!***********************************************************************
! 
!                            VERSION 2D 
! 
!***********************************************************************

  integer*4             :: iel,i,j,Ndim,Neq,Nel,Np,Nfp,Nf
  integer*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq)
  integer*4             :: ind_uf(refElPol%Nfacenodes*phys%Neq,refElPol%Nfaces)
  integer*4             :: F_el(refElPol%Nfaces)
  integer*4             :: ind_loc(refElPol%Nfaces,refElPol%Nfacenodes*phys%Neq),perm(refElPol%Nfacenodes*phys%Neq)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8                :: ue(Mesh%Nnodesperelem*phys%Neq)
  real*8                :: uf(refElPol%Nfacenodes*phys%Neq,refElPol%Nfaces)
  real*8,allocatable    :: Cvel(:,:),Hel(:,:),Hfel(:,:)
  ! Matrices to be filled in case of non-constant stabilization
  real*8,allocatable    :: Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:)
  
  
  real*8,allocatable    :: tau_save(:,:)
  real*8,allocatable    :: xy_g_save(:,:)
  integer               :: indtausave(refElPol%Nfaces*refElPol%Ngauss1d)
  
  real*8                :: tau_save_el(refElPol%Nfaces*refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Nfaces*refElPol%Ngauss1d,2)
  logical               :: save_tau

		IF (utils%printint>1) THEN
   WRITE(6,*) '*************************************************'
   WRITE(6,*) '*           CONVECTION MATRICES                 *'
   WRITE(6,*) '*************************************************' 
		END IF	   
		
		
  
  
  save_tau = switch%saveTau
  
  
  Ndim        = Mesh%ndim
  Neq         = phys%Neq
  Nel         = Mesh%Nelems
  Np          = refElPol%Nnodes2D
  Nfp         = refElPol%Nfacenodes
  Nf          = refElPol%Nfaces



  if (save_tau) then
     allocate(tau_save(refElPol%Nfaces*Nel*refElPol%Ngauss1d,phys%neq))
     allocate(xy_g_save(refElPol%Nfaces*Nel*refElPol%Ngauss1d,2))
     tau_save = 0.
     xy_g_save = 0.
  endif

  ind_loc = 0
  DO i = 1,Nf
     DO j = 1,Neq*Nfp
        ind_loc(i,j) = Neq*Nfp*(i-1) + j
     END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Nfp,Neq,perm)


!$OMP PARALLEL DEFAULT(SHARED) & 
!$OMP PRIVATE(iel,ind_ue,ind_uf,Xe,F_dir_el,F_el,i,j,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel,tau_save_el,xy_g_save_el,indtausave) 
  ALLOCATE(Cvel(Neq*Np,Neq*Np))
  ALLOCATE(Hel(Neq*Np,Nf*Neq*Nfp))
  ALLOCATE(Hfel(Nf*Neq*Nfp,Nf*Neq*Nfp))
  ALLOCATE(Del(Neq*Np,Neq*Np))
  ALLOCATE(Eel(Neq*Np,Neq*Nf*Nfp))
  ALLOCATE(Dfel(Neq*Nf*Nfp,Neq*Np))
  ALLOCATE(Efel(Neq*Nf*Nfp,Neq*Nf*Nfp))
!$OMP DO SCHEDULE(STATIC) 
  DO iel = 1,Nel
     Cvel    = 0.d0
     Hel     = 0.d0
     Hfel    = 0.d0

     Del    = 0.d0
     Eel    = 0.d0
     Dfel   = 0.d0
     Efel   = 0.d0
     
     ! Coordinates of the nodes of the element
     Xe = Mesh%X(Mesh%T(iel,:),:)

     ! Dirichlet boundary conditions in the current element
     F_dir_el = Mesh%Fdir(iel,:)
     
     ! Indices to extract the elemental and face solution
     F_el = Mesh%F(iel,:)     
     ind_ue = (iel-1)*Np*Neq + (/ (i, i = 1, Neq*Np) /)
     ind_uf =  tensorSumInt((/(i,i=1,Neq*Nfp)/),(F_el-1)*Neq*Nfp)
    
     ! Elemental and face solution     
     DO j = 1,Nf
        uf(:,j) = sol%u_tilde(ind_uf(:,j))
     END DO
     ue = sol%u(ind_ue)

     ! Compute the matrices for the element
     CALL elemental_matrices(iel,Xe,F_dir_el,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel,tau_save_el,xy_g_save_el)

     if (save_tau) then
        indtausave = (iel-1)*refElPol%Nfaces*refElPol%Ngauss1d + (/(i,i =1,refElPol%Nfaces*refElPol%Ngauss1d)/)
        tau_save(indtausave,:)= tau_save_el
        xy_g_save(indtausave,:)= xy_g_save_el
     endif
     
     ! Flip the faces for the elements for which the face is already been counted
     DO j = 1,Nf
        if (Mesh%flipface(iel,j)) then
           Hfel(ind_loc(j,:),:) = Hfel(ind_loc(j,perm(:)),:)
           Hfel(:,ind_loc(j,:)) = Hfel(:,ind_loc(j,perm(:)))
           Dfel(ind_loc(j,:),:) = Dfel(ind_loc(j,perm(:)),:)
           Efel(ind_loc(j,:),:) = Efel(ind_loc(j,perm(:)),:)
           Efel(:,ind_loc(j,:)) = Efel(:,ind_loc(j,perm(:)))           
        END if
     END DO

     ! Store the matrices for all the elements
     elMat%Cv(:,:,iel)  = Cvel
     elMat%H(:,:,iel)   = Hel
     elMat%Hf(:,:,iel)  = Hfel   
     IF (numer%stab > 1) THEN
        elMat%D(:,:,iel)  = Del
        elMat%E(:,:,iel)  = Eel
        elMat%Df(:,:,iel) = Dfel
        elMat%Ef(:,:,iel) = Efel     
     END IF
  END DO
!$OMP END DO

 DEALLOCATE(Cvel,Hel,Hfel)
 DEALLOCATE(Del,Eel,Dfel,Efel)  
!$OMP END PARALLEL 

   if (save_tau) then
      write(6,*) "Saving tau in the faces"
      call saveMatrix(tau_save,'tau_save')
      call saveMatrix(xy_g_save,'xy_g_save')
      deallocate(tau_save,xy_g_save)
      write(6,*) "Done saving tau!"
   endif
   
  IF (utils%printint>1) THEN
     WRITE(6,*) "Done!"
  END IF
  
  CONTAINS
			
			
!***************************************************
! Elemental matrices computation: 2D case, 2-3-4 eqs
!***************************************************
					SUBROUTINE elemental_matrices(iel,Xe,F_dir_el,ue,uf,Cvel,Hel,Hfel,Del,Eel,Dfel,Efel,tau_save_el,xy_g_save_el)
       integer,intent(IN)			     :: iel 
					  real*8,intent(IN)         :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
					  logical,intent(IN)        :: F_dir_el(refElPol%Nfaces)
       real*8,intent(IN)         :: ue(Mesh%Nnodesperelem*phys%Neq)
       real*8,intent(IN)         :: uf(refElPol%Nfacenodes*phys%Neq,refElPol%Nfaces)		
					  real*8,intent(INOUT)      :: Cvel(:,:),Hel(:,:),Hfel(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:)

							integer*4                 :: g,NGauss,ifa,i,j
							real*8                    :: dvolu,dline,isdir,isext,xyDerNorm_g
							real*8                    :: xy(refElPol%Ngauss2d,ndim),ufi(Nfp*neq),ufg(refElPol%Ngauss1d,neq),ueg(refElPol%Ngauss2d,neq)
							real*8                    :: xyf(refElPol%Ngauss1d,ndim)
							real*8                    :: xyDer(refElPol%Ngauss1d,ndim)

       real*8                     :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
       real*8                     :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
       real*8                     :: detJ(refElPol%Ngauss2d)
       real*8                     :: iJ11(refElPol%Ngauss2d),iJ12(refElPol%Ngauss2d)
       real*8                     :: iJ21(refElPol%Ngauss2d),iJ22(refElPol%Ngauss2d)
       
							integer*4                 :: ind_ff(phys%Neq*refElPol%Nfacenodes),ind_fe(phys%Neq*refElPol%Nfacenodes)
							integer*4                 :: ind_ass(refElPol%Nnodes2D),ind_asf(refElPol%Nfacenodes),ind(refElPol%Nfacenodes)
							real*8                    :: t_g(ndim),n_g(ndim), bn
							real*8,dimension(Np)      :: Nxg,Nyg
							real*8,dimension(Np,Np)   :: NxNi,NyNi,NNxy
							real*8,dimension(Nfp,Nfp) :: NiNi
							real*8,dimension(neq,neq) :: A
							real*8,dimension(Neq)     :: tau
       real*8,dimension(Neq,Neq)     :: taumat
       real*8                    :: upg(refElPol%Ngauss2d,phys%npv),upgf(refElPol%Ngauss1d,phys%npv)
#ifdef TEMPERATURE
       real*8,dimension(Np,Np)   :: NNi
							real*8,dimension(neq,neq) :: GG
							
							real*8,allocatable,dimension(:,:)    ::Dre,Drifte
							real*8                    :: Telect
#endif														
							real*8                    :: Hloc(phys%Neq*refElPol%Nfacenodes,phys%Neq*refElPol%Nfacenodes)
							real*8, parameter         :: tol = 1e-12							
       real*8                    :: b(refElPol%Ngauss2d,2),divb(refElPol%Ngauss2d),drift(refElPol%Ngauss2d,2),divbg,driftg(2)
       real*8                    :: b_f(refElPol%Ngauss1d,2)				
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
       real*8                    :: Bmod_g(refElPol%Ngauss2d),Bmod_x,Bmod_y,Nx_ax(Np)
       real*8                    :: Bt(Mesh%Nnodesperelem),Bt_g(refElPol%Ngauss2d)              
#endif
       real*8,allocatable    :: Mf_loc(:,:),Massf(:,:)
       real*8,intent(out)  :: tau_save_el(:,:),xy_g_save_el(:,:)
       
       
       if (save_tau) then
          tau_save_el = 0.
          xy_g_save_el = 0.
       endif 
       
       
       ind_ass = (/ (i, i = 0, Neq*(Np-1),Neq) /)
       ind_asf = (/ (i, i = 0, Neq*(Nfp-1),Neq) /)

       !***********************************
       !    Volume computation
       !***********************************
							
#ifdef TEMPERATURE
       ALLOCATE(Dre(Np,Np))
       ALLOCATE(Drifte(Np*Neq,Np*Neq))
       Drifte = 0.
       Dre = 0.
#endif	
							
							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)

							! Solution at Gauss points
							ueg = matmul(refElPol%N2D,transpose(reshape(ue,[neq,Np])))

#ifndef MOVINGEQUILIBRIUM					
!************ CASE OF STATIC EQUILIBRIUM ****************************
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
           divb = matmul(refElPol%N2D,phys%divb(Mesh%T(iel,:)))
           drift = matmul(refElPol%N2D,phys%drift(Mesh%T(iel,:),:))
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),b,divb,drift)
							END IF
#ifdef TEMPERATURE
       IF (switch%driftvel) THEN
						    ! Physical variables at Gauss points
						    CALL cons2phys(ueg,upg)
						    
						    ! Drift at Gauss points
						    drift = phys%dfcoef*drift ! I multiply for the temperature later
       END IF
#endif							
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
       b_nod = 0. ! Nnodes
       Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2) ! Nnodes
       b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod ! Nnodes
       b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod ! Nnodes
       b = matmul(refElPol%N2D,b_nod)              ! Ngauss
#ifdef TEMPERATURE        
       IF (switch%driftvel) THEN
          ! Physical variables at Gauss points
          CALL cons2phys(ueg,upg)  ! Ngauss
        
          ! Drift at Gauss points       
          Bmod_g = matmul(refElPol%N2D,Bmod) ! Ngauss
          Bt_g   = matmul(refElPol%N2D,phys%Bt(Mesh%T(iel,:)))  ! Ngauss 
       END IF       
#endif          
#endif		
!************ END SELECTION MOVING EQUILIBRIUM **********************

							! Loop in 2D Gauss points
       Ngauss = refElPol%NGauss2D
       J11 = matmul(refElPol%Nxi2D,Xe(:,1))                           ! ng x 1
       J12 = matmul(refElPol%Nxi2D,Xe(:,2))                           ! ng x 1
       J21 = matmul(refElPol%Neta2D,Xe(:,1))                          ! ng x 1
       J22 = matmul(refElPol%Neta2D,Xe(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
				   iJ11 =  J22/detJ
				   iJ12 = -J12/detJ
				   iJ21 = -J21/detJ
				   iJ22 =  J11/detJ
				                 
							DO g = 1,NGauss
							
										! Integration weight
										dvolu = refElPol%gauss_weights2D(g)*detJ(g)
          IF (switch%axisym) THEN										
             dvolu = dvolu*xy(g,1)
          END IF

										! x and y derivatives of the shape functions
						    Nxg = iJ11(g)*refElPol%Nxi2D(g,:) + iJ12(g)*refElPol%Neta2D(g,:)
						    Nyg = iJ21(g)*refElPol%Nxi2D(g,:) + iJ22(g)*refElPol%Neta2D(g,:)     										

         ! Jacobian matrices
          CALL jacobianMatrices(ueg(g,:),A)   

									! Contribution of the current integration point to the elemental matrix
									 NxNi = tensorProduct(Nxg,refElPol%N2D(g,:))*dvolu
          NyNi = tensorProduct(Nyg,refElPol%N2D(g,:))*dvolu
          NNxy = b(g,1)*NxNi+b(g,2)*NyNi									
           
          ! Assembly convection matrix
          DO i=1,Neq
             DO j=1,Neq
                Cvel(i+ind_ass,j+ind_ass) = Cvel(i+ind_ass,j+ind_ass) + A(i,j)*NNxy
             END DO
          END DO
                                                
#ifdef TEMPERATURE
          ! Divergence of b at the Gauss points  
#ifndef MOVINGEQUILIBRIUM
          divbg = divb(g)
#else      
							   IF (switch%axisym) THEN
							      Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
							   ELSE
							      Nx_ax = Nxg
							   END IF  
          divbg = dot_product(Nx_ax,b_nod(:,1))+dot_product(Nyg,b_nod(:,2))                
#endif          
             
          ! Matrix for the curvature term
          CALL GimpMatrix(ueg(g,:),divbg,GG)   
           
          ! Assembly curvature term in the convection matrix          
          NNi = tensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu   
          DO i=2,2
             DO j=1,Neq
                Cvel(i+ind_ass,j+ind_ass) = Cvel(i+ind_ass,j+ind_ass) + GG(i,j)*NNi
             END DO
          END DO
                    
          IF (switch%driftvel) THEN               
#ifdef MOVINGEQUILIBRIUM          
			       Bmod_x  = dot_product(Nxg,Bmod) ! 1
			       Bmod_y  = dot_product(Nyg,Bmod) ! 1
			       driftg(1) =  phys%dfcoef*Bt_g(g)*Bmod_y/Bmod_g(g)**3 ! I multiply for the temperature later
			       driftg(2) = -phys%dfcoef*Bt_g(g)*Bmod_x/Bmod_g(g)**3 ! I multiply for the temperature later          
#else
          driftg(1) = drift(g,1)
          driftg(2) = drift(g,2)
#endif			       
             Telect = upg(g,8)
             Dre = Dre + TensorProduct(refElPol%N2D(g,:),Telect*( Nxg*driftg(1)+ Nyg*driftg(2) ))*dvolu
          END IF

#endif  
							END DO
						
#ifdef TEMPERATURE	
						CALL expand_matrix_mass(Dre,Neq,Drifte)
      Cvel = Cvel - Drifte
      DEALLOCATE(Drifte,Dre)
#endif

       !***********************************
  					! Faces computations
       !***********************************
							NGauss = refElPol%Ngauss1D
							ALLOCATE(Mf_loc(Neq*Nfp,Neq*Nfp))
												
       ! Loop in the faces of the element
							DO ifa = 1,Nf
							   Mf_loc = 0.d0
										Hloc = 0.d0; 
										
										! Indices
          ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Nfp/))
          ind_ff = (ifa-1)*neq*Nfp + (/ (i, i = 1, neq*Nfp) /)										
	                       
          ! Set the face solution and isdir
          isdir = 0.        
          isext = 0.  
          xyf    = matmul(refElPol%N1D,Xe(refElPol%face_nodes(ifa,:),:))
          IF (F_dir_el(ifa)) THEN
             isdir = 1.          
             CALL analytical_solution(xyf(:,1),xyf(:,2),ufg) 
          ELSE
													ufi = uf(:,ifa)
													IF (Mesh%flipface(iel,ifa)) THEN
																	ufi = ufi(perm)
													END IF             
             ufg = matmul(refElPol%N1D,transpose(reshape(ufi,[neq,Nfp])))
          END IF 
          
          IF (numer%stab>1) THEN
               CALL cons2phys(ufg,upgf)
          END IF
#ifdef PARALL    
! probably useless, to be verified      
          IF ( (Mesh%F(iel,ifa).gt.Mesh%Nintfaces)  ) THEN
             IF (Mesh%boundaryFlag(Mesh%F(iel,ifa)-Mesh%Nintfaces).ne.0) THEN
                isext = 1.
             END IF
          ENDIF 
#else
          IF (Mesh%F(iel,ifa) > Mesh%Nintfaces ) isext = 1.
#endif
#ifndef MOVINGEQUILIBRIUM			
!************ CASE OF STATIC EQUILIBRIUM ****************************
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
 					        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:))
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),b_f)  
										END IF				
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
          b_f = matmul(refElPol%N1D,b_nod(refElPol%face_nodes(ifa,:),:))
#endif
!************ END SELECTION MOVING EQUILIBRIUM **********************

										! Shape function derivatives at Gauss points
										xyDer  = matmul(refElPol%Nxi1D,Xe(refElPol%face_nodes(ifa,:),:))
										
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
             CALL jacobianMatrices(ufg(g,:),A)  

												! Contribution of the current integration point to the elemental matrix
													NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
													
            ! Assembly convection matrix
            DO i=1,Neq
               DO j=1,Neq
                  Hloc(i+ind_asf,j+ind_asf) = Hloc(i+ind_asf,j+ind_asf) + bn*A(i,j)*NiNi
               END DO
            END DO
          		
          		
          		!************* Non constant stabilization************************											
            IF (numer%stab>1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab<5) THEN
                  CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g,iel,ifa,isext,xyf(g,:),tau)
               ELSE
                  taumat = 0.
                  CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g,xyf(g,:),isext,iel,taumat)
               ENDIF
               if (save_tau) then
                  if (numer%stab<5) then
                     tau_save_el((ifa-1)*Ngauss+g,:) = tau
                  else
                     tau_save_el((ifa-1)*Ngauss+g,1) = taumat(1,1)
                     tau_save_el((ifa-1)*Ngauss+g,2) = taumat(2,2)
                     tau_save_el((ifa-1)*Ngauss+g,3) = taumat(3,3)
                     tau_save_el((ifa-1)*Ngauss+g,4) = taumat(4,4)
                  endif
                  xy_g_save_el((ifa-1)*Ngauss+g,:) = xyf(g,:)
               endif
               
										     ! Expand
               IF (numer%stab<5) THEN
                  DO i=1,Neq
                     ind = i+ind_asf
                     Mf_loc(ind,ind) = Mf_loc(ind,ind) + tau(i)*NiNi
                  END DO 										        
               ELSE
                  DO i=1,Neq
                     DO j=1,Neq
                        Mf_loc(i+ind_asf,j+ind_asf) = Mf_loc(i+ind_asf,j+ind_asf) + taumat(i,j)*NiNi
                     END DO
                  END DO 
               ENDIF
            END IF					       
            !************* End non constant stabilization************************
			
										END DO ! Gauss points

										! Elemental assembly
										Hel(ind_fe,ind_ff)  = Hel(ind_fe,ind_ff) + Hloc*(1-isdir)															
										Hfel(ind_ff,ind_ff) = Hfel(ind_ff,ind_ff) + Hloc*(1-isext)
										
										 IF (numer%stab>1) THEN 
										    ! Elemental assembly
										    Del(ind_fe,ind_fe)  = Del(ind_fe,ind_fe) + Mf_loc										
										    Eel(ind_fe,ind_ff)  = Eel(ind_fe,ind_ff) + Mf_loc*(1-isdir)
										    Efel(ind_ff,ind_ff) = Efel(ind_ff,ind_ff) + Mf_loc*(1-isext)
										    Dfel(ind_ff,ind_fe) = Dfel(ind_ff,ind_fe) + Mf_loc*(1-isext)	
										 END IF									
							END DO ! faces
							
							DEALLOCATE(Mf_loc)

					END SUBROUTINE elemental_matrices				
													
!*******************************************
!           AUXILIARY ROUTINES
!*******************************************				
								!*******************************************
								! Expand matrix mass type
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

#ifdef TEMPERATURE									
								!*******************************************
								! Expand matrix C: 4 equations
								!*******************************************									
								SUBROUTINE expand_matrix(C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C13,C14,C21,C22,C23,C24,C31,C32,C33,C34,C41,C42,C43,C44
								real*8,intent(inout),dimension(size(C11,1)*4,size(C11,1)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+1,(j-1)*4+1) = C((i-1)*4+1,(j-1)*4+1) + C11(i,j)
             C((i-1)*4+1,(j-1)*4+2) = C((i-1)*4+1,(j-1)*4+2) + C12(i,j)
             C((i-1)*4+1,(j-1)*4+3) = C((i-1)*4+1,(j-1)*4+3) + C13(i,j)
             C((i-1)*4+1,(j-1)*4+4) = C((i-1)*4+1,(j-1)*4+4) + C14(i,j)
             
             C((i-1)*4+2,(j-1)*4+1) = C((i-1)*4+2,(j-1)*4+1) + C21(i,j)
             C((i-1)*4+2,(j-1)*4+2) = C((i-1)*4+2,(j-1)*4+2) + C22(i,j)
             C((i-1)*4+2,(j-1)*4+3) = C((i-1)*4+2,(j-1)*4+3) + C23(i,j)
             C((i-1)*4+2,(j-1)*4+4) = C((i-1)*4+2,(j-1)*4+4) + C24(i,j)
                          
             C((i-1)*4+3,(j-1)*4+1) = C((i-1)*4+3,(j-1)*4+1) + C31(i,j)
             C((i-1)*4+3,(j-1)*4+2) = C((i-1)*4+3,(j-1)*4+2) + C32(i,j)
             C((i-1)*4+3,(j-1)*4+3) = C((i-1)*4+3,(j-1)*4+3) + C33(i,j)             
             C((i-1)*4+3,(j-1)*4+4) = C((i-1)*4+3,(j-1)*4+4) + C34(i,j)             

             C((i-1)*4+4,(j-1)*4+1) = C((i-1)*4+4,(j-1)*4+1) + C41(i,j)
             C((i-1)*4+4,(j-1)*4+2) = C((i-1)*4+4,(j-1)*4+2) + C42(i,j)
             C((i-1)*4+4,(j-1)*4+3) = C((i-1)*4+4,(j-1)*4+3) + C43(i,j)             
             C((i-1)*4+4,(j-1)*4+4) = C((i-1)*4+4,(j-1)*4+4) + C44(i,j)             
           END DO
        END DO
								END SUBROUTINE expand_matrix
#else								
								!*******************************************
								! Expand matrix C: 2 equations
								!*******************************************									
								SUBROUTINE expand_matrix(C11,C12,C21,C22,C)        
								real*8,intent(in),dimension(:,:)        :: C11,C12,C21,C22
								real*8,intent(inout),dimension(size(C11,1)*2,size(C11,1)*2) :: C
								integer*4 :: i,j,n,m
								
								n=size(C11,1)
								m=size(C11,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*2+1,(j-1)*2+1) = C((i-1)*2+1,(j-1)*2+1) + C11(i,j)
             C((i-1)*2+1,(j-1)*2+2) = C((i-1)*2+1,(j-1)*2+2) + C12(i,j)
             C((i-1)*2+2,(j-1)*2+1) = C((i-1)*2+2,(j-1)*2+1) + C21(i,j)
             C((i-1)*2+2,(j-1)*2+2) = C((i-1)*2+2,(j-1)*2+2) + C22(i,j)
           END DO
        END DO
								END SUBROUTINE expand_matrix
#endif

#ifdef TEMPERATURE
							
								SUBROUTINE expand_matrix_Gimp(G21,G22,G23,G24,C)        
								real*8,intent(in),dimension(:,:)        :: G21,G22,G23,G24
								real*8,intent(inout),dimension(size(G21,1)*4,size(G21,1)*4) :: C
								integer*4 :: i,j,n,m
								
								n=size(G21,1)
								m=size(G21,2)

        DO i=1,n
           DO j=1,m
             C((i-1)*4+2,(j-1)*4+1) = C((i-1)*4+2,(j-1)*4+1) + G21(i,j)
             C((i-1)*4+2,(j-1)*4+2) = C((i-1)*4+2,(j-1)*4+2) + G22(i,j)
             C((i-1)*4+2,(j-1)*4+3) = C((i-1)*4+2,(j-1)*4+3) + G23(i,j) 
             C((i-1)*4+2,(j-1)*4+4) = C((i-1)*4+2,(j-1)*4+4) + G24(i,j)            
           END DO
        END DO
								END SUBROUTINE expand_matrix_Gimp
#endif
     
							!*****************************************
							! Set permutations for flipping faces
							!****************************************     
							 SUBROUTINE set_permutations(n,m,perm)
							 integer, intent(IN)  :: n,m
							 integer, intent(OUT) :: perm(:)
							 integer              :: i
							 integer              :: temp(m,n/m),templr(m,n/m)

							 IF (mod(n,m) .ne. 0) then
							    WRITE(6,*) 'Error! n must be a multiple of m'
							    STOP
							 END IF

							 templr = 0
							 temp = reshape((/(i,i=1,n)/),(/m,n/m/))
							 DO i = 1,n/m
							    templr(:,i) = temp(:,n/m-i+1)
							 END DO
							 perm = reshape(templr,(/n/))
							 END SUBROUTINE set_permutations     
							 
							 
								!*******************************************
								! Expand matrix with different stabilization
								! parameters
								!*******************************************									
								SUBROUTINE expand_matrix_stab(A,n,tau,B)
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
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
								END SUBROUTINE expand_matrix_stab	
								
        
        
 							!*******************************************
								! Expand matrix to use a tensor stabilization
        ! parameter
								!*******************************************			               
								SUBROUTINE expand_matrix_stab_mat(A,n,tau,B)        
										integer*4,intent(in) :: n
										real*8,intent(in)    :: A(:,:),tau(:,:)
										real*8,intent(out)   :: B(n*size(A,1),n*size(A,1))
										integer*4            :: i,j,k,t,m
										
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
								END SUBROUTINE expand_matrix_stab_mat
								

							       
END SUBROUTINE hdg_convectionmatrices

#endif



