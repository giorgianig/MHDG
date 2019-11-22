!*****************************************
! project: MHDG
! file: hdg_ParallDiffusionMatrices.f90
! date: 21/09/2017
! Generate the matrices related to the 
! parallel diffusion term
!*****************************************

SUBROUTINE HDG_paralldiffmatrices()
  USE globals
  USE LinearAlgebra
  USE analytical
  USE physics
  USE printUtils
  USE MPI_OMP
  
  IMPLICIT NONE
#ifdef TOR3D   
!***********************************************************************
! 
!                            VERSION 3D TOROIDAL
! 
!***********************************************************************
  integer*4             :: itor,iel,iel3,i,j,k,Ndim,Neq,N2D,Np,Nfp,Nf,Nfdir,dd
  integer*4             :: ntorloc,itorg,Np1Dpol,Np1Dtor,Np2D,Nfl
  integer*4             :: ind_ue(refElTor%Nnodes3D*phys%Neq), ind_qe(refElTor%Nnodes3D*phys%Neq*3)
  integer*4             :: ind_uf(refElTor%Nft*phys%Neq)
  integer               :: ind_dim(refElPol%Nfaces+2),ind_sta(refElPol%Nfaces+2),aux,ifa   
  real*8                :: tdiv(numer%ntor+1)
  real*8                :: htor,tel(refElTor%Nnodes1d)
  integer*4             :: F_el(refElPol%Nfaces)
  integer*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8,allocatable    :: qe(:),ue(:),uf(:)
  real*8,allocatable    :: TU(:,:),TUh(:,:),TUhf(:,:),TQ(:,:),TQh(:,:),TQhf(:,:),Tf(:),Tfh(:),Tfhf(:)
  real*8                :: coefi,coefe
  integer               :: delta(refElPol%Nfaces+2)

		IF (utils%printint>1) THEN
   WRITE(6,*) '*************************************************'
   WRITE(6,*) '*      PARALLEL DIFFUSION MATRICES              *'
   WRITE(6,*) '*************************************************' 
		END IF	   
		
		coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
		coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
		
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

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(itor,itorg,tel,iel,iel3,ifa,ind_qe,ind_ue,ind_uf,Xe,F_dir_el,&
!$OMP&                                 F_el,i,j,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf,dd,delta) 

  ALLOCATE(TU(Neq*Np,Neq*Np))
  ALLOCATE(TUh(Neq*Np,Neq*Nfp))
  ALLOCATE(TQ(Neq*Np,Neq*Np*Ndim))
  ALLOCATE(TQh(Neq*Np,Neq*Np*Ndim))
  ALLOCATE(TUhf(Neq*Nfp,Neq*Nfp))
  ALLOCATE(TQhf(Neq*Nfp,Neq*Ndim*Np))
  ALLOCATE(Tf(Neq*Np))
  ALLOCATE(Tfh(Neq*Np))
  ALLOCATE(Tfhf(Neq*Nfp))  
  ALLOCATE(qe(Ndim*Neq*Np))
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

        TU    = 0.d0
        TUh   = 0.d0
        TQ    = 0.d0
        TQh   = 0.d0
        TUhf  = 0.d0
        TQhf  = 0.d0
        Tf    = 0.d0
        Tfh   = 0.d0
        Tfhf  = 0.d0

        ! Coordinates of the nodes of the element
        Xe = Mesh%X(Mesh%T(iel,:),:)

        ! Dirichlet boundary conditions in the current element
        F_dir_el = Mesh%Fdir(iel,:)
        
						  ! Index of 3D element
						  iel3 = (itor-1)*N2d+iel
						          
        ! Indices to extract the elemental and face solution
        F_el = Mesh%F(iel,:)    
        ind_qe = (iel3-1)*Neq*Np*Ndim + (/(i,i=1,Neq*Np*Ndim) /)
        ind_ue = (iel3-1)*Np*Neq + (/(i,i=1,Neq*Np) /)

        ! Elemental and face solution     
        dd = 1+(itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Nfl)*Neq
        delta(1) = dd +(iel-1)*Np2D*Neq
        delta(2:refElPol%Nfaces+1) = dd + (N2D*Np2D+(F_el-1)*Nfl)*Neq
        delta(refElPol%Nfaces+2) = dd + (N2D*Np2D+(Mesh%Nfaces-Nfdir)*Nfl+(iel-1)*Np2D)*Neq
!        if (MPIvar%ntor.eq.1) then
        if (itor==numer%ntor) then
           delta(refElPol%Nfaces+2) = 1+(iel-1)*Np2D*Neq
        end if
!        endif
        DO ifa = 1,Nf+2
            DO i=0,ind_dim(ifa)-1 
               ind_uf(i+ind_sta(ifa)) = delta(ifa)+i
            END DO
        END DO  	
  
                        	
        uf = sol%u_tilde(ind_uf)		
        qe = sol%q(ind_qe)
        ue = sol%u(ind_ue)

        ! Compute the matrices for the element
        CALL elemental_matrices(itor,iel,Xe,tel,F_dir_el,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)

        ! Flip the faces for the elements for which the face is already been counted
        DO j = 1,Nf
           if (Mesh%flipface(iel,j)) then          
               TUhf(ind_loc(j,:),:) = TUhf(ind_loc(j,perm),:)
               TUhf(:,ind_loc(j,:)) = TUhf(:,ind_loc(j,perm))
               TQhf(ind_loc(j,:),:) = TQhf(ind_loc(j,perm),:)
               Tfhf(ind_loc(j,:))   = Tfhf(ind_loc(j,perm))            
           END if
        END DO
        
        ! Store the matrices for all the elements
        elMat%Cv(:,:,iel3)  = elMat%Cv(:,:,iel3) - TU
        elMat%H(:,:,iel3)   = elMat%H(:,:,iel3)  - TUh
        elMat%Hf(:,:,iel3)  = elMat%Hf(:,:,iel3) - TUhf     
        elMat%TQ(:,:,iel3)  = TQ-TQh
        elMat%Tf(:,iel3)    = Tf-Tfh
        elMat%TQhf(:,:,iel3)= TQhf
        elMat%Tfhf(:,iel3)  = Tfhf
        
   !call saveMatrix(elMat%Cv(:,:,iel),'Cv2')     
   !stop     
     END DO
  END DO
!$OMP END DO

  DEALLOCATE(qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)
  
!$OMP END PARALLEL 
  
  
  
  CONTAINS
!*******************************************************
! 3D elemental matrices computation
!*******************************************************
					SUBROUTINE elemental_matrices(itor,iel,Xe,tel,F_dir_el,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)

       integer,intent(IN)			     :: iel,itor 
					  real*8,intent(IN)         :: Xe(:,:),tel(:)
					  logical,intent(IN)        :: F_dir_el(:)
					  real*8,intent(IN)         :: qe(:)
       real*8,intent(IN)         :: ue(:)
       real*8,intent(IN)         :: uf(:)		
					  real*8,intent(INOUT)      :: TU(:,:),TUh(:,:),TQ(:,:),TQh(:,:),TUhf(:,:),TQhf(:,:),Tf(:),Tfh(:),Tfhf(:)

					  ! Dimension and indices
       integer*4                  :: NGaussTor,NgaussPol
							integer*4                  :: igtor,igpol,g,ifa,i,j,k
							integer                    :: ind(refElPol%Ngauss1d),nodes2(refelPol%Nnodes1d),ind2(refElPol%Ngauss2d)
							integer*4,allocatable      :: ind_ff(:),ind_fe(:),ind_fg(:)
							integer*4,allocatable      :: ind_ass(:),ind_asg(:),ind_asf(:),ind_ash(:)
							! Volume computations
							real*8                     :: dvolu,htor,dvolu1d
							real*8                     :: xy(refElPol%Ngauss2d,2),teg(refElTor%Ngauss1d)
       real*8                     :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
       real*8                     :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
       real*8                     :: detJ(refElPol%Ngauss2d)
       real*8                     :: iJ11(refElPol%Ngauss2d),iJ12(refElPol%Ngauss2d)
       real*8                     :: iJ21(refElPol%Ngauss2d),iJ22(refElPol%Ngauss2d)
       real*8                     :: Nxg(Np2d),Nyg(Np2d),Nx_ax(Np2d)
       real*8                     :: N1g(refElTor%Nnodes1D),N1xg(refElTor%Nnodes1D),N1xg_cart(refElTor%Nnodes1D),N2g(Np2d)   
       real*8,allocatable         :: Ni(:),Nidvolu(:),Nr(:),Nrr(:),Nz(:),Nt(:)    
       ! Faces computations
       real*8,allocatable         :: aux(:,:),xyf(:,:)
       real*8,allocatable         :: thetafg(:),dsurf(:),n_g(:,:)       
							real*8                     :: xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
							real*8                     :: isdir,dsurfg,isext
							real*8                     :: t_g(refElPol%Ngauss1d,2)	
							real*8,pointer             :: Nfi(:,:),Nfg(:)
							integer                    :: Npifa
       ! Magnetic field, volume computations							       
       real*8                     :: b(refElPol%Ngauss2d*refElTor%Ngauss1d,Ndim),b2d(refElPol%Ngauss2d,Ndim)
       ! Magnetic field, faces computations
       real*8                     :: bn
       real*8,allocatable         :: b_f(:,:)
       real*8                     :: b_f2d(refElPol%Ngauss1d,Ndim)
       ! Solution in volume
       real*8                     :: ueg(refElPol%Ngauss2d*refElTor%Ngauss1d,Neq),upg(refElPol%Ngauss2d*refElTor%Ngauss1d,phys%npv)
       real*8                     :: qeg(refElPol%Ngauss2d*refElTor%Ngauss1d,Neq*Ndim)
       real*8                     :: qer(Np,Neq*Ndim)
       ! Solution in faces
       real*8,allocatable         :: ufi(:),ufg(:,:),upfg(:,:),qfg(:,:)
       ! Routine specific part
							real*8                    :: Qpr(3,Neq),Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq),W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(Ndim,Neq)       							
							real*8,allocatable        :: NNf(:),NN(:,:),NN_f(:,:),NNf_f(:),NNi(:,:)
							real*8,allocatable        :: TUhloc(:,:),TQhloc(:,:),Tfhloc(:)

       



							real*8, parameter         :: tol = 1e-12			
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif       											

						 ! toroidal element size
						 htor = tel(refElTor%Nnodes1d) - tel(1)


       !***********************************
       !    Volume computation
       !***********************************
       ALLOCATE(ind_ass(Np))
       ALLOCATE(ind_asg(Np))
       ind_ass = (/ (i, i = 0, Neq*(Np-1),Neq) /)
       ind_asg = (/ (i, i = 0, Neq*(Np-1)*Ndim,Neq*Ndim) /)


							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)

							! Solution at Gauss points
							qer = transpose(reshape(qe,[Neq*Ndim,Np]))                      ! np x neq*ndim
							ueg = matmul(refElTor%N3D,transpose(reshape(ue,[Neq,Np])))      ! ng x neq 
							qeg = matmul(refElTor%N3D,qer)                                  ! ng x neq*ndim

							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)
!**********************************************************************										
!                 Magnetic field at elements Gauss points
!**********************************************************************
#ifndef MOVINGEQUILIBRIUM					
!************ CASE OF STATIC EQUILIBRIUM ****************************
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b2d = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
           DO igtor = 1,refElTor%NGauss1D
              ind2 = (igtor-1)*refElPol%NGauss2D + (/(j,j=1,refElPol%NGauss2D)/)
              b(ind2,:) = b2d
           END DO  
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),teg,b)
							END IF					
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
#endif		
!************ END SELECTION MOVING EQUILIBRIUM **********************

							! Loop in 2D Gauss points
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

														IF (switch%axisym) THEN
															  Nx_ax = Nxg+1./xy(igpol,1)*N2g
															  N1xg = N1xg_cart/xy(igpol,1) 
														ELSE
															  Nx_ax = Nxg
															  N1xg = N1xg_cart
														END IF
														
							       ! 3D shape functions 
							       Ni  = col(TensorProduct(N2g,N1g))      ! 3D shape function
							       Nrr = col(TensorProduct(Nx_ax,N1g))   ! 3D shape function, derivative in r for gradient
							       Nr  = col(TensorProduct(Nxg,N1g))      ! 3D shape function, derivative in r for divergence
							       Nz  = col(TensorProduct(Nyg,N1g))      ! 3D shape function, derivative in z
							       Nt  = col(TensorProduct(N2g,N1xg))     ! 3D shape function, derivative in t
										    
										    ! 3D integration weight
										    dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d
              IF (switch%axisym) THEN										
                 dvolu = dvolu*xy(igpol,1)
              END IF
              							       
              ! Compute Q^T^(k-1)
              Qpr = reshape(qeg(g,:),(/Ndim,Neq/))

              ! Compute V(U^(k-1)) 
              call computeVi(ueg(g,:),Vveci)
              call computeVe(ueg(g,:),Vvece)

              ! Compute dV_dU (k-1)
              call compute_dV_dUi(ueg(g,:),dV_dUi)
              call compute_dV_dUe(ueg(g,:),dV_dUe)

              ! Compute Alpha(U^(k-1))
              Alphai = computeAlphai(ueg(g,:))
              Alphae = computeAlphae(ueg(g,:))

              ! Compute dAlpha/dU^(k-1)
              call compute_dAlpha_dUi(ueg(g,:),dAlpha_dUi)
              call compute_dAlpha_dUe(ueg(g,:),dAlpha_dUe)
              
               gmi = dot_product(matmul(Qpr,Vveci),b(g,:))             ! scalar     
               gme = dot_product(matmul(Qpr,Vvece),b(g,:))             ! scalar         
               Taui = matmul(Qpr,dV_dUi)      ! Ndim x Neq
               Taue = matmul(Qpr,dV_dUe)      ! Ndim x Neq
				           
				           
				           ! Compute W(U^(k-1)) 
				           call compute_W(ueg(g,:),W)						      
				           
				           ! Compute dW_dU(U^(k-1)) 
				           call compute_dW_dU(ueg(g,:),dW_dU)	
				           						      
				           ! Compute s(U^(k-1)) 
				           call compute_S(ueg(g,:),s)		

				           ! Compute s(U^(k-1)) 
				           call compute_dS_dU(ueg(g,:),ds_dU)	
         						      						      
               Zet  = matmul(Qpr,dW_dU)       ! Ndim x Neq
          
							       ! Shape functions tensor product
							       Nidvolu  = Ni*dvolu			
							       NNi      = tensorProduct(Ni,Nidvolu)		
              NNf      = (b(g,1)*Nr+b(g,2)*Nz+b(g,3)*Nt)*dvolu
              NN       = tensorProduct(NNf,Ni) 
              
              							                 
									! Tensor products between shape functions    
									
!          NNf =  (Nxg*b(g,1)+Nyg*b(g,2))*dvolu
!          NN  =  tensorProduct((Nxg*b(g,1)+Nyg*b(g,2)),refElPol%N2D(g,:))*dvolu          
!          NNi = tensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu
          
          
        
          DO i=3,4
             IF (i==3) THEN
                DO j=1,4
                   TU(i+ind_ass,j+ind_ass) = TU(i+ind_ass,j+ind_ass) + coefi*(gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),b(g,:))))*NN + ( dot_product(Zet(:,j),b(g,:)) +ds_dU(j))*NNi
                   DO k=1,Ndim
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + coefi*Alphai*Vveci(j)*b(g,k)*NN
                      IF (j==4) THEN
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + W*NNi*b(g,k)
                      END IF
                   END DO
                END DO 
                Tf(i+ind_ass) = Tf(i+ind_ass) + coefi*Alphai*( dot_product (matmul(transpose(Taui),b(g,:)),ueg(g,:)) )*NNf + s*Nidvolu
               
             ELSE
                DO j=1,4
                   TU(i+ind_ass,j+ind_ass) = TU(i+ind_ass,j+ind_ass) + coefe*(gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),b(g,:))))*NN - (dot_product(Zet(:,j),b(g,:))+ds_dU(j))*NNi
                   DO k=1,Ndim
                      TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + coefe*Alphae*Vvece(j)*b(g,k)*NN
                      IF (j==4) THEN
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) - W*NNi*b(g,k)
                      END IF     
                   END DO                              
                END DO  
                Tf(i+ind_ass) = Tf(i+ind_ass) + coefe*Alphae*( dot_product (matmul(transpose(Taue),b(g,:)),ueg(g,:))  )*NNf - s*Nidvolu
             END IF           
          END DO
           
							END DO
							END DO
       DEALLOCATE(ind_ass,ind_asg)
       
 
 
 
 
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
							   
							   
							   ALLOCATE(NN_f(Npifa,Npifa))			
							   ALLOCATE(NNf_f(Npifa))
							   ALLOCATE(TUhloc(Neq*Npifa,Neq*Npifa))
							   ALLOCATE(TQhloc(Neq*Npifa,Neq*Npifa*Ndim))
							   ALLOCATE(Tfhloc(Neq*Npifa))
	         ALLOCATE(ind_asf(Npifa))
          ALLOCATE(ind_ash(Npifa))  
          ALLOCATE(b_f(NGaussPol*NGaussTor,3))	   
          ALLOCATE(ufi(Npifa*Neq))   
          ALLOCATE(ufg(NGaussPol*NGaussTor,Neq))
          ALLOCATE(qfg(NGaussPol*NGaussTor,Neq*Ndim))
          
          ind_asf = (/ (i, i = 0, Neq*(Npifa-1),Neq) /)
          ind_ash = (/ (i, i = 0, Neq*(Npifa-1)*Ndim,Neq*Ndim) /)
          
          TUhloc = 0.
          TQhloc = 0.
          Tfhloc = 0.           
                     
                     
          ! Set the face solution and isdir
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
           
         ! Gradient solution at face gauss points
         qfg    = matmul(Nfi, transpose(reshape(qe(ind_fg),([Neq*Ndim,Npifa]) ))) ! Ng_face x Neq*Ndim

!**********************************************************************										
!                 Magnetic field at faces Gauss points
!**********************************************************************
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

		             gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))             ! scalar    
		             gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))  
		             Taui = matmul(Qpr,dV_dUi)      ! 2x3
		             Taue = matmul(Qpr,dV_dUe)      ! 2x3

											    ! Tensor products between shape functions    
				            NNf_f =  bn*Nfg*dsurfg
				            NN_f  =  bn*tensorProduct(Nfg,Nfg)*dsurfg

                ! Contribution of the current integration point to the elemental matrix
                DO i=3,4
                   IF (i==3) THEN
                      DO j=1,4
                         TUhloc(i+ind_asf,j+ind_asf)         = TUhloc(i+ind_asf,j+ind_asf) + coefi*(gmi*dAlpha_dUi(j)+Alphai*(  dot_product(Taui(:,j),b_f(g,:))  ))*NN_f
                         DO k=1,Ndim
                            TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) + coefi*Alphai*Vveci(j)*b_f(g,k)*NN_f
                         END DO
                      END DO 
                      Tfhloc(i+ind_asf) = Tfhloc(i+ind_asf) + coefi*Alphai*(  dot_product (matmul(transpose(Taui),b_f(g,:)),ufg(g,:))   )*NNf_f
                   ELSE
                      DO j=1,4
                         TUhloc(i+ind_asf,j+ind_asf)         = TUhloc(i+ind_asf,j+ind_asf) + coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b_f(g,:)) ))*NN_f
                         DO k=1,Ndim
                            TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) + coefe*Alphae*Vvece(j)*b_f(g,k)*NN_f
                         END DO
                      END DO  
                      Tfhloc(i+ind_asf) = Tfhloc(i+ind_asf) + coefe*Alphae*(dot_product (matmul(transpose(Taue),b_f(g,:)),ufg(g,:))    )*NNf_f
                   END IF           
                END DO
                     
											
   										END DO ! Gauss points
										END DO ! Gauss points
										
										DEALLOCATE(b_f,n_g,dsurf,xyf,thetafg,NN_f,NNf_f,ind_asf,ind_ash,ufg,ufi,qfg)

						 ! elemental assembly
						 TUh(ind_fe,ind_ff)  = TUh(ind_fe,ind_ff)  + (1-isdir)*TUhloc
						 TQh(ind_fe,ind_fg)  = TQh(ind_fe,ind_fg)  + TQhloc
						 TUhf(ind_ff,ind_ff) = TUhf(ind_ff,ind_ff) + (1-isext)*TUhloc
						 Tfh(ind_fe)         = Tfh(ind_fe)         + Tfhloc
						 TQhf(ind_ff,ind_fg) = TQhf(ind_ff,ind_fg) + (1-isext)*TQhloc
						 Tfhf(ind_ff)        = Tfhf(ind_ff)        + (1-isext)*Tfhloc                                
       
       DEALLOCATE(TUhloc,TQhloc,Tfhloc)
       DEALLOCATE(ind_ff,ind_fe,ind_fg)
       END DO ! End loop in faces

					END SUBROUTINE elemental_matrices




				!*******************************************
				! Expand matrix T: 4 equations
				!*******************************************									
				SUBROUTINE expand_matrixT(C31,C32,C33,C34,C41,C42,C43,C44,C)        
				real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C41,C42,C43,C44
				real*8,intent(inout),dimension(size(C31,1)*4,size(C31,1)*4) :: C
				integer*4 :: i,j,n,m
				
				n=size(C31,1)
				m=size(C31,2)

    DO i=1,n
       DO j=1,m
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
				END SUBROUTINE expand_matrixT
								
								
								
				!*******************************************
				! Expand matrix TQ: 8 equations
				!*******************************************									
				SUBROUTINE expand_matrixTQ(C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48,C)        
				real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48
				real*8,intent(inout),dimension(size(C31,1)*4,size(C31,1)*8) :: C
				integer*4 :: i,j,n,m
				
				n=size(C31,1)
				m=size(C31,2)

    DO i=1,n
       DO j=1,m
         C((i-1)*4+3,(j-1)*8+1) = C((i-1)*4+3,(j-1)*8+1) + C31(i,j)
         C((i-1)*4+3,(j-1)*8+2) = C((i-1)*4+3,(j-1)*8+2) + C32(i,j)
         C((i-1)*4+3,(j-1)*8+3) = C((i-1)*4+3,(j-1)*8+3) + C33(i,j)             
         C((i-1)*4+3,(j-1)*8+4) = C((i-1)*4+3,(j-1)*8+4) + C34(i,j)
         C((i-1)*4+3,(j-1)*8+5) = C((i-1)*4+3,(j-1)*8+5) + C35(i,j)
         C((i-1)*4+3,(j-1)*8+6) = C((i-1)*4+3,(j-1)*8+6) + C36(i,j)                     
         C((i-1)*4+3,(j-1)*8+7) = C((i-1)*4+3,(j-1)*8+7) + C37(i,j)
         C((i-1)*4+3,(j-1)*8+8) = C((i-1)*4+3,(j-1)*8+8) + C38(i,j)

         C((i-1)*4+4,(j-1)*8+1) = C((i-1)*4+4,(j-1)*8+1) + C41(i,j)
         C((i-1)*4+4,(j-1)*8+2) = C((i-1)*4+4,(j-1)*8+2) + C42(i,j)
         C((i-1)*4+4,(j-1)*8+3) = C((i-1)*4+4,(j-1)*8+3) + C43(i,j)             
         C((i-1)*4+4,(j-1)*8+4) = C((i-1)*4+4,(j-1)*8+4) + C44(i,j)
         C((i-1)*4+4,(j-1)*8+5) = C((i-1)*4+4,(j-1)*8+5) + C45(i,j)
         C((i-1)*4+4,(j-1)*8+6) = C((i-1)*4+4,(j-1)*8+6) + C46(i,j)                     
         C((i-1)*4+4,(j-1)*8+7) = C((i-1)*4+4,(j-1)*8+7) + C47(i,j)
         C((i-1)*4+4,(j-1)*8+8) = C((i-1)*4+4,(j-1)*8+8) + C48(i,j)         
       END DO
    END DO
				END SUBROUTINE expand_matrixTQ					
				
				
				
				!*******************************************
				! Expand matrix Tf
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

				
				
							!*****************************************
							! Set permutations for flipping faces
							!****************************************     
							 SUBROUTINE set_permutations(Np1Dpol,Np1Dtor,Neq,perm)
							 integer, intent(IN)  :: Np1Dpol,Np1Dtor,Neq
							 integer, intent(OUT) :: perm(:)
							 integer              :: i,j
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
							
END SUBROUTINE hdg_paralldiffmatrices

#else
!***********************************************************************
! 
!                            VERSION 2D 
! 
!***********************************************************************

  integer*4             :: iel,i,j,Ndim,Neq,Nel,Np,Nfp,Nf
  integer*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq), ind_qe(Mesh%Nnodesperelem*phys%Neq*Mesh%Ndim)
  integer*4             :: ind_uf(refElPol%Nfacenodes*phys%Neq,refElPol%Nfaces)
  integer*4             :: F_el(refElPol%Nfaces)
  integer*4             :: ind_loc(refElPol%Nfaces,refElPol%Nfacenodes*phys%Neq),perm(refElPol%Nfacenodes*phys%Neq)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8                :: ue(Mesh%Nnodesperelem*phys%Neq),qe(Mesh%Nnodesperelem*phys%Neq*Mesh%Ndim)
  real*8                :: uf(refElPol%Nfacenodes*phys%Neq,refElPol%Nfaces)
  real*8,allocatable    :: TU(:,:),TUh(:,:),TUhf(:,:),TQ(:,:),TQh(:,:),TQhf(:,:),Tf(:),Tfh(:),Tfhf(:)
  real*8                :: coefi,coefe

		IF (utils%printint>1) THEN
   WRITE(6,*) '*************************************************'
   WRITE(6,*) '*      PARALLEL DIFFUSION MATRICES              *'
   WRITE(6,*) '*************************************************' 
		END IF	   
		
		coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
		coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
  Ndim        = Mesh%ndim
  Neq         = phys%Neq
  Nel         = Mesh%Nelems
  Np          = refElPol%Nnodes2D
  Nfp         = refElPol%Nfacenodes
  Nf          = refElPol%Nfaces

  ind_loc = 0
  DO i = 1,Nf
     DO j = 1,Neq*Nfp
        ind_loc(i,j) = Neq*Nfp*(i-1) + j
     END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Nfp,Neq,perm)



!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iel,ind_qe,ind_ue,ind_uf,Xe,F_dir_el,F_el,i,j,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf) 
  ALLOCATE(TU(Neq*Np,Neq*Np))
  ALLOCATE(TUh(Neq*Np,Nf*Neq*Nfp))
  ALLOCATE(TQ(Neq*Np,Neq*Np*Ndim))
  ALLOCATE(TQh(Neq*Np,Neq*Np*Ndim))
  ALLOCATE(TUhf(Neq*Nf*Nfp,Neq*Nf*Nfp))
  ALLOCATE(TQhf(Neq*Nf*Nfp,Neq*Ndim*Np))
  ALLOCATE(Tf(Neq*Np))
  ALLOCATE(Tfh(Neq*Np))
  ALLOCATE(Tfhf(Neq*Nf*Nfp))  


!$OMP DO SCHEDULE(STATIC) 
  DO iel = 1,Nel
     TU    = 0.d0
     TUh   = 0.d0
     TQ    = 0.d0
     TQh   = 0.d0
     TUhf  = 0.d0
     TQhf  = 0.d0
     Tf    = 0.d0
     Tfh   = 0.d0
     Tfhf  = 0.d0

     ! Coordinates of the nodes of the element
     Xe = Mesh%X(Mesh%T(iel,:),:)

     ! Dirichlet boundary conditions in the current element
     F_dir_el = Mesh%Fdir(iel,:)
     
     ! Indices to extract the elemental and face solution
     F_el = Mesh%F(iel,:)    
     ind_qe = (iel-1)*Neq*Np*Ndim + (/ (i, i = 1,Neq*Np*Ndim) /)
     ind_ue = (iel-1)*Np*Neq + (/ (i, i = 1, Neq*Np) /)
     ind_uf =  tensorSumInt((/(i,i=1,Neq*Nfp)/),(F_el-1)*Neq*Nfp)
    
     ! Elemental and face solution     
     DO j = 1,Nf
        uf(:,j) = sol%u_tilde(ind_uf(:,j))
     END DO
     qe = sol%q(ind_qe)
     ue = sol%u(ind_ue)

     ! Compute the matrices for the element
     CALL elemental_matrices(iel,Xe,F_dir_el,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)

     ! Flip the faces for the elements for which the face is already been counted
     DO j = 1,Nf
        if (Mesh%flipface(iel,j)) then          
            TUhf(ind_loc(j,:),:) = TUhf(ind_loc(j,perm),:)
            TUhf(:,ind_loc(j,:)) = TUhf(:,ind_loc(j,perm))
            TQhf(ind_loc(j,:),:) = TQhf(ind_loc(j,perm),:)
            Tfhf(ind_loc(j,:))   = Tfhf(ind_loc(j,perm))            
        END if
     END DO
     
     ! Store the matrices for all the elements
     elMat%Cv(:,:,iel)  = elMat%Cv(:,:,iel) - TU
     elMat%H(:,:,iel)   = elMat%H(:,:,iel)  - TUh
     elMat%Hf(:,:,iel)  = elMat%Hf(:,:,iel) - TUhf     
     elMat%TQ(:,:,iel)  = TQ-TQh
     elMat%Tf(:,iel)    = Tf-Tfh
     elMat%TQhf(:,:,iel)= TQhf
     elMat%Tfhf(:,iel)  = Tfhf
     
!call saveMatrix(elMat%Cv(:,:,iel),'Cv2')     
!stop     
  END DO
!$OMP END DO

  DEALLOCATE(TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)
  
 !$OMP END PARALLEL 
  
  
  
  CONTAINS
!*******************************************************
! Elemental matrices computation
!*******************************************************
					SUBROUTINE elemental_matrices(iel,Xe,F_dir_el,qe,ue,uf,TU,TUh,TQ,TQh,TUhf,TQhf,Tf,Tfh,Tfhf)

       integer,intent(IN)			     :: iel 
					  real*8,intent(IN)         :: Xe(:,:)
					  logical,intent(IN)        :: F_dir_el(:)
					  real*8,intent(IN)         :: qe(:)
       real*8,intent(IN)         :: ue(:)
       real*8,intent(IN)         :: uf(:,:)		
					  real*8,intent(INOUT)      :: TU(:,:),TUh(:,:),TQ(:,:),TQh(:,:),TUhf(:,:),TQhf(:,:),Tf(:),Tfh(:),Tfhf(:)

							integer*4                 :: g,NGauss,ifa,ind_ass(Np),ind_asg(Np),ind_asf(Nfp),ind_ash(Nfp),i,j,k
							real*8                    :: dvolu,dline,isdir,isext,xyDerNorm_g
							real*8                    :: xy(refElPol%Ngauss2d,ndim),ufi(Nfp*neq),ufg(refElPol%Ngauss1d,neq),ueg(refElPol%Ngauss2d,neq)
							real*8                    :: qeg(refElPol%Ngauss2d,neq*Ndim), qfg(refElPol%Ngauss1d,neq*Ndim),qer(Np,neq*Ndim)
							real*8                    :: xyf(refElPol%Ngauss1d,ndim)
							real*8                    :: xyDer(refElPol%Ngauss1d,ndim)
       
       real*8                     :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
       real*8                     :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
       real*8                     :: detJ(refElPol%Ngauss2d)
       real*8                     :: iJ11(refElPol%Ngauss2d),iJ12(refElPol%Ngauss2d)
       real*8                     :: iJ21(refElPol%Ngauss2d),iJ22(refElPol%Ngauss2d)
       							
							integer*4                 :: ind_ff(neq*Nfp),ind_fe(neq*Nfp),ind_fg(neq*Nfp*Ndim)
							real*8                    :: t_g(ndim),n_g(ndim), bn
							real*8,dimension(Np)      :: Nxg,Nyg,NNf
				   real*8                    :: NN(Np,Np)
							real*8                    :: NN_f(Nfp,Nfp),NNf_f(Nfp)
				   real*8                    :: NNi(Np,Np)
							real*8,dimension(Neq*Nfp,Neq*Nfp) :: TUhloc
							real*8,dimension(Neq*Nfp,Neq*Nfp*Ndim) :: TQhloc
							real*8,dimension(Neq*Nfp) :: Tfhloc
							real*8                    :: Qpr(Ndim,Neq),Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq),W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)

							real*8, parameter         :: tol = 1e-12			
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
#endif       											
       real*8                    :: b(refElPol%Ngauss2d,ndim)
       real*8                    :: b_f(1:refElPol%Ngauss1d,1:ndim)


       !***********************************
       !    Volume computation
       !***********************************

       ind_ass = (/ (i, i = 0, Neq*(Np-1),Neq) /)
       ind_asg = (/ (i, i = 0, Neq*(Np-1)*Ndim,Neq*Ndim) /)
       ind_asf = (/ (i, i = 0, Neq*(Nfp-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Nfp-1)*Ndim,Neq*Ndim) /)
							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)

							! Solution at Gauss points
							qer = transpose(reshape(qe,[Neq*Ndim,Np]))                   ! np x neq*ndim
							ueg = matmul(refElPol%N2D,transpose(reshape(ue,[Neq,Np])))      ! ng x nvar 
							qeg = matmul(refElPol%N2D,qer)                                  ! ng x nvar*ndim

!**********************************************************************										
!                 Magnetic field at elements Gauss points
!**********************************************************************
#ifndef MOVINGEQUILIBRIUM					
!************ CASE OF STATIC EQUILIBRIUM ****************************
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),b)
							END IF					
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
       b_nod = 0.
       Bmod   = sqrt(phys%Br(Mesh%T(iel,:))**2+phys%Bz(Mesh%T(iel,:))**2+phys%Bt(Mesh%T(iel,:))**2)
       b_nod(:,1) = phys%Br(Mesh%T(iel,:))/Bmod
       b_nod(:,2) = phys%Bz(Mesh%T(iel,:))/Bmod
       b = matmul(refElPol%N2D,b_nod)        
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
     										
         ! Compute Q^T^(k-1)
         Qpr = reshape(qeg(g,:),(/Ndim,Neq/))

         ! Compute V(U^(k-1)) 
         call computeVi(ueg(g,:),Vveci)
         call computeVe(ueg(g,:),Vvece)

         ! Compute dV_dU (k-1)
         call compute_dV_dUi(ueg(g,:),dV_dUi)
         call compute_dV_dUe(ueg(g,:),dV_dUe)

         ! Compute Alpha(U^(k-1))
         Alphai = computeAlphai(ueg(g,:))
         Alphae = computeAlphae(ueg(g,:))

         ! Compute dAlpha/dU^(k-1)
         call compute_dAlpha_dUi(ueg(g,:),dAlpha_dUi)
         call compute_dAlpha_dUe(ueg(g,:),dAlpha_dUe)
         
          gmi = dot_product(matmul(Qpr,Vveci),b(g,:))             ! scalar     
          gme = dot_product(matmul(Qpr,Vvece),b(g,:))             ! scalar         
          Taui = matmul(Qpr,dV_dUi)      ! Ndim x Neq
          Taue = matmul(Qpr,dV_dUe)      ! Ndim x Neq
				      
				      
				      ! Compute W(U^(k-1)) 
				      call compute_W(ueg(g,:),W)						      
				      
				      ! Compute dW_dU(U^(k-1)) 
				      call compute_dW_dU(ueg(g,:),dW_dU)	
				      						      
				      ! Compute s(U^(k-1)) 
				      call compute_S(ueg(g,:),s)		

				      ! Compute s(U^(k-1)) 
				      call compute_dS_dU(ueg(g,:),ds_dU)	
    						      						      
          Zet  = matmul(Qpr,dW_dU)       ! Ndim x Neq
          
          
									! Tensor products between shape functions    
          NNf =  (Nxg*b(g,1)+Nyg*b(g,2))*dvolu
          NN  =  tensorProduct(NNf,refElPol%N2D(g,:))
          NNi = tensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu
          
          DO i=3,4
             IF (i==3) THEN
                DO j=1,4
                   TU(i+ind_ass,j+ind_ass) = TU(i+ind_ass,j+ind_ass) + coefi*(gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),b(g,:))))*NN + ( dot_product(Zet(:,j),b(g,:)) +ds_dU(j))*NNi
                   DO k=1,Ndim
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + coefi*Alphai*Vveci(j)*b(g,k)*NN
                      IF (j==4) THEN
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + W*NNi*b(g,k)
                      END IF
                   END DO
                END DO 
                Tf(i+ind_ass) = Tf(i+ind_ass) + coefi*Alphai*( dot_product (matmul(transpose(Taui),b(g,:)),ueg(g,:)) )*NNf + s*refElPol%N2D(g,:)*dvolu
             ELSE
                DO j=1,4
                   TU(i+ind_ass,j+ind_ass) = TU(i+ind_ass,j+ind_ass) + coefe*(gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),b(g,:))))*NN - (dot_product(Zet(:,j),b(g,:))+ds_dU(j))*NNi
                   DO k=1,Ndim
                      TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) + coefe*Alphae*Vvece(j)*b(g,k)*NN
                      IF (j==4) THEN
                         TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) = TQ(i+ind_ass,k+(j-1)*Ndim+ind_asg) - W*NNi*b(g,k)
                      END IF     
                   END DO                              
                END DO  
                Tf(i+ind_ass) = Tf(i+ind_ass) + coefe*Alphae*( dot_product (matmul(transpose(Taue),b(g,:)),ueg(g,:))  )*NNf - s*refElPol%N2D(g,:)*dvolu
             END IF           
          END DO
          
          

							END DO

       !***********************************
  					! Faces computations
       !***********************************
							NGauss = refElPol%Ngauss1D
							
       ! Loop in the faces of the element
							DO ifa = 1,Nf
							   TUhloc = 0.
							   TQhloc = 0.
							   Tfhloc = 0.							
										! Indices
          ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Nfp/))
          ind_fg = reshape(tensorSumInt( (/ (i, i = 1, neq*ndim) /),neq*ndim*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Nfp*ndim/))
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
           
         ! Gradient solution at face gauss points
         qfg    = matmul(refElPol%N1D,qer(refElPol%face_nodes(ifa,:),:))
										
										
										
!**********************************************************************										
!                 Magnetic field at faces Gauss points
!**********************************************************************
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

		          gmi = dot_product(matmul(Qpr,Vveci),b_f(g,:))             ! scalar    
		          gme = dot_product(matmul(Qpr,Vvece),b_f(g,:))  
		          Taui = matmul(Qpr,dV_dUi)      ! 2x3
		          Taue = matmul(Qpr,dV_dUe)      ! 2x3

											! Tensor products between shape functions    
				        NNf_f =  bn*refElPol%N1D(g,:)*dline
				        NN_f  =  bn*tensorProduct( refElPol%N1D(g,:),refElPol%N1D(g,:))*dline    

           ! Contribution of the current integration point to the elemental matrix
          DO i=3,4
             IF (i==3) THEN
                DO j=1,4
                   TUhloc(i+ind_asf,j+ind_asf)         = TUhloc(i+ind_asf,j+ind_asf) + coefi*(gmi*dAlpha_dUi(j)+Alphai*(  dot_product(Taui(:,j),b_f(g,:))  ))*NN_f
                   DO k=1,Ndim
                      TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) + coefi*Alphai*Vveci(j)*b_f(g,k)*NN_f
                   END DO
                END DO 
                Tfhloc(i+ind_asf) = Tfhloc(i+ind_asf) + coefi*Alphai*(  dot_product (matmul(transpose(Taui),b_f(g,:)),ufg(g,:))   )*NNf_f
             ELSE
                DO j=1,4
                   TUhloc(i+ind_asf,j+ind_asf)         = TUhloc(i+ind_asf,j+ind_asf) + coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b_f(g,:)) ))*NN_f
                   DO k=1,Ndim
                      TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) = TQhloc(i+ind_asf,k+(j-1)*Ndim+ind_ash) + coefe*Alphae*Vvece(j)*b_f(g,k)*NN_f
                   END DO
                END DO  
                Tfhloc(i+ind_asf) = Tfhloc(i+ind_asf) + coefe*Alphae*(dot_product (matmul(transpose(Taue),b_f(g,:)),ufg(g,:))    )*NNf_f
             END IF           
          END DO
                     
											
          END DO  ! End loop in Gauss points

						 ! elemental assembly
						 TUh(ind_fe,ind_ff)  = TUh(ind_fe,ind_ff)  + (1-isdir)*TUhloc
						 TQh(ind_fe,ind_fg)  = TQh(ind_fe,ind_fg)  + TQhloc
						 TUhf(ind_ff,ind_ff) = TUhf(ind_ff,ind_ff) + (1-isext)*TUhloc
						 Tfh(ind_fe)         = Tfh(ind_fe)         + Tfhloc
						 TQhf(ind_ff,ind_fg) = TQhf(ind_ff,ind_fg) + (1-isext)*TQhloc
						 Tfhf(ind_ff)        = Tfhf(ind_ff)        + (1-isext)*Tfhloc                                
    
       END DO ! End loop in faces

					END SUBROUTINE elemental_matrices




				!*******************************************
				! Expand matrix T: 4 equations
				!*******************************************									
				SUBROUTINE expand_matrixT(C31,C32,C33,C34,C41,C42,C43,C44,C)        
				real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C41,C42,C43,C44
				real*8,intent(inout),dimension(size(C31,1)*4,size(C31,1)*4) :: C
				integer*4 :: i,j,n,m
				
				n=size(C31,1)
				m=size(C31,2)

    DO i=1,n
       DO j=1,m
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
				END SUBROUTINE expand_matrixT
								
								
								
				!*******************************************
				! Expand matrix TQ: 8 equations
				!*******************************************									
				SUBROUTINE expand_matrixTQ(C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48,C)        
				real*8,intent(in),dimension(:,:)        :: C31,C32,C33,C34,C35,C36,C37,C38,C41,C42,C43,C44,C45,C46,C47,C48
				real*8,intent(inout),dimension(size(C31,1)*4,size(C31,1)*8) :: C
				integer*4 :: i,j,n,m
				
				n=size(C31,1)
				m=size(C31,2)

    DO i=1,n
       DO j=1,m
         C((i-1)*4+3,(j-1)*8+1) = C((i-1)*4+3,(j-1)*8+1) + C31(i,j)
         C((i-1)*4+3,(j-1)*8+2) = C((i-1)*4+3,(j-1)*8+2) + C32(i,j)
         C((i-1)*4+3,(j-1)*8+3) = C((i-1)*4+3,(j-1)*8+3) + C33(i,j)             
         C((i-1)*4+3,(j-1)*8+4) = C((i-1)*4+3,(j-1)*8+4) + C34(i,j)
         C((i-1)*4+3,(j-1)*8+5) = C((i-1)*4+3,(j-1)*8+5) + C35(i,j)
         C((i-1)*4+3,(j-1)*8+6) = C((i-1)*4+3,(j-1)*8+6) + C36(i,j)                     
         C((i-1)*4+3,(j-1)*8+7) = C((i-1)*4+3,(j-1)*8+7) + C37(i,j)
         C((i-1)*4+3,(j-1)*8+8) = C((i-1)*4+3,(j-1)*8+8) + C38(i,j)

         C((i-1)*4+4,(j-1)*8+1) = C((i-1)*4+4,(j-1)*8+1) + C41(i,j)
         C((i-1)*4+4,(j-1)*8+2) = C((i-1)*4+4,(j-1)*8+2) + C42(i,j)
         C((i-1)*4+4,(j-1)*8+3) = C((i-1)*4+4,(j-1)*8+3) + C43(i,j)             
         C((i-1)*4+4,(j-1)*8+4) = C((i-1)*4+4,(j-1)*8+4) + C44(i,j)
         C((i-1)*4+4,(j-1)*8+5) = C((i-1)*4+4,(j-1)*8+5) + C45(i,j)
         C((i-1)*4+4,(j-1)*8+6) = C((i-1)*4+4,(j-1)*8+6) + C46(i,j)                     
         C((i-1)*4+4,(j-1)*8+7) = C((i-1)*4+4,(j-1)*8+7) + C47(i,j)
         C((i-1)*4+4,(j-1)*8+8) = C((i-1)*4+4,(j-1)*8+8) + C48(i,j)         
       END DO
    END DO
				END SUBROUTINE expand_matrixTQ					
				
				
				
				!*******************************************
				! Expand matrix Tf
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
							
END SUBROUTINE hdg_paralldiffmatrices

#endif
