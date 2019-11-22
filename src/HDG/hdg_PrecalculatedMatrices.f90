!*****************************************
! project: MHDG
! file: hdg_precalculatedMatrices.f90
! date: 20/12/2016
! Generate the matrices that remain constant
! during the whole computation
!*****************************************


SUBROUTINE HDG_precalculatedmatrices()
  USE globals
  USE analytical
  USE LinearAlgebra
  USE printUtils
  USE MPI_OMP
  
  IMPLICIT NONE

#ifdef TOR3D  
!***********************************************************************
! 
!                            VERSION 3D TOROIDAL
! 
!***********************************************************************
  integer*4             :: iel,i,j,Ndim,Neq,Np,Nfp,Nf,Nfl,ierr,ith,ntorloc,itorg
  integer*4             :: N2d,itor,iel3,Np1Dpol,Np1Dtor,Np2D
  integer*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8                :: tdiv(numer%ntor+1)
  real*8                :: htor,tel(refElTor%Nnodes1d)
  real*8,allocatable    :: Ael(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:), fe(:)
  real*8,allocatable    :: Bel(:,:),Cel(:,:),Gel(:,:),iLel(:,:),Pel(:,:),Qel(:,:), Lfel(:,:),Qfel(:,:)

  
  integer :: OMP_GET_THREAD_NUM
  ith=0
  
  IF (MPIvar%glob_id.eq.0) THEN
					IF (utils%printint>0) THEN
						WRITE(6,*) '*************************************************'
						WRITE(6,*) '*           PRECALCULATED MATRICES              *'
						WRITE(6,*) '*************************************************' 
					END IF	   
		END IF
		
  
  Ndim        = 3                                               ! Number of dimensions
  N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Neq         = phys%Neq                                        ! Number of equations
  Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np          = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Nfl         = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfp         = Np2D*2+refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
  Nf          = refElPol%Nfaces                                 ! Number of faces in the 2D mesh

 
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
 
!#ifdef PARALL
!  itor_start = (MPIvar%itor-1)*ntor/MPIvar%npartor+1
!  itor_end   = MPIvar%itor*ntor/MPIvar%npartor
!#else
!  itor_start = 1
!  itor_end   = ntor
!#endif
#ifdef PARALL
    IF (MPIvar%ntor.gt.1) THEN
       ntorloc = numer%ntor/MPIvar%ntor+1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif

!*****************
! Loop in elements
!*****************  

!$OMP PARALLEL DEFAULT(SHARED) &
#ifndef MOVINGEQUILIBRIUM
!$OMP PRIVATE(itor,itorg,iel3,iel,Xe,tel,F_dir_el,j,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Lfel,Qfel,Gel)
#else
!$OMP PRIVATE(itor,itorg,iel3,iel,Xe,tel,F_dir_el,j,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif

  ALLOCATE(Ael(Neq*Np,Neq*Np))
  ALLOCATE(Del(Neq*Np,Neq*Np))
  ALLOCATE(Eel(Neq*Np,Neq*Nfp))
  ALLOCATE(Dfel(Neq*Nfp,Neq*Np))
  ALLOCATE(Efel(Neq*Nfp,Neq*Nfp))
  ALLOCATE(fe(Neq*Np))
  ALLOCATE(Bel(Neq*Ndim*Np,Neq*Np))
  ALLOCATE(Cel(Neq*Ndim*Np,Neq*Nfp))
  ALLOCATE(iLel(Neq*Ndim*Np,Neq*Ndim*Np))
  ALLOCATE(Lfel(Neq*Nfp,Neq*Ndim*Np))  
#ifndef MOVINGEQUILIBRIUM
  ALLOCATE(Gel(Neq*Np,Neq*Np))
  ALLOCATE(Pel(Neq*Np,Ndim*Neq*Np))
  ALLOCATE(Qel(Neq*Np,Ndim*Neq*Np))
  ALLOCATE(Qfel(Neq*Nfp,Neq*Ndim*Np))
#endif

!$OMP DO SCHEDULE(STATIC)
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif     
     tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg+1)-tdiv(itorg))
     DO iel = 1,N2d
     
     Ael    = 0.d0     
     Del    = 0.d0
     Eel    = 0.d0
     Dfel   = 0.d0
     Efel   = 0.d0
     fe     = 0.d0     
     Bel  = 0.
     Cel  = 0.
     iLel = 0.
     Lfel = 0.
#ifndef MOVINGEQUILIBRIUM
     Gel  = 0.    
     Pel  = 0.
     Qel  = 0.     
     Qfel = 0.     
#endif
     
     ! Index of 3D element
     iel3 = (itor-1)*N2d+iel
     
     ! Coordinates of the nodes of the element
     Xe = Mesh%X(Mesh%T(iel,:),:)

     ! Dirichlet boundary conditions in the current element
     F_dir_el = Mesh%Fdir(iel,:)
 
     ! Compute the matrices for the element
#ifndef MOVINGEQUILIBRIUM     
     CALL elemental_matrices(itor,iel,Xe,tel,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Gel,Lfel,Qfel)
#else
     CALL elemental_matrices(itor,iel,Xe,tel,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif     
     ! Flip the faces for the elements for which the face is already been counted
     DO j = 1,Nf
        IF (Mesh%flipface(iel,j)) THEN
           Dfel(ind_loc(j,:),:) = Dfel(ind_loc(j,perm(:)),:)
           Efel(ind_loc(j,:),:) = Efel(ind_loc(j,perm(:)),:)
           Efel(:,ind_loc(j,:)) = Efel(:,ind_loc(j,perm(:)))
           Lfel(ind_loc(j,:),:) = Lfel(ind_loc(j,perm),:)
#ifndef MOVINGEQUILIBRIUM
           Qfel(ind_loc(j,:),:) = Qfel(ind_loc(j,perm),:)           
#endif           
        END if
     END DO

     ! Store the matrices for all the elements
     IF (.not.switch%steady) THEN
        elMat%M(:,:,iel3)  = Ael
     END IF
     elMat%D(:,:,iel3)  = Del
     elMat%E(:,:,iel3)  = Eel
     elMat%Df(:,:,iel3) = Dfel
     elMat%Ef(:,:,iel3) = Efel
     elMat%S(:,iel3)    = fe
     elMat%B(:,:,iel3)  = Bel
     elMat%C(:,:,iel3)  = Cel
     elMat%iL(:,:,iel3) = iLel
     elMat%Lf(:,:,iel3)   = Lfel
#ifndef MOVINGEQUILIBRIUM
     elMat%G(:,:,iel3)  = Gel 
     elMat%P(:,:,iel3)  = Pel-Qel
     elMat%Qf(:,:,iel3)   = Qfel
#endif     
     END DO
  END DO
!$OMP END DO

  DEALLOCATE(Ael,Del,Eel,Dfel,Efel,fe)
  DEALLOCATE(Bel,Cel,iLel,Lfel)
#ifndef MOVINGEQUILIBRIUM
  DEALLOCATE(Pel,Qel,Qfel,Gel)
#endif  

!$OMP END PARALLEL


 
  IF (MPIvar%glob_id.eq.0) THEN
					IF (utils%printint>0) THEN
						  WRITE(6,*) "Done!"
					END IF
  END IF
  
  
  CONTAINS
  

!*****************************************
! Elemental matrices computation: 3D case
!*****************************************
#ifndef MOVINGEQUILIBRIUM
					SUBROUTINE elemental_matrices(itor,iel,Xe,tel,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Gel,Lfel,Qfel )
#else
     SUBROUTINE elemental_matrices(itor,iel,Xe,tel,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif					
							integer*4,intent(IN)  :: itor,iel
       real*8,intent(IN)     :: Xe(:,:),tel(:)							
       logical,intent(IN)    :: F_dir_el(:)
       real*8,intent(INOUT)  :: Ael(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:), fe(:)     
#ifndef MOVINGEQUILIBRIUM
       real*8,intent(INOUT)  :: Bel(:,:),Cel(:,:),iLel(:,:),Pel(:,:),Qel(:,:),Lfel(:,:),Qfel(:,:),Gel(:,:) 
#else
       real*8,intent(INOUT)  :: Bel(:,:),Cel(:,:),iLel(:,:),Lfel(:,:) 
#endif       
					  ! Dimension and indices
       integer*4                  :: Ndim,Np2d,Np1Dpol,Np1Dtor,Np,Nfl,Nfp,Nf,NGaussTor,NgaussPol,Npifa
							integer*4                  :: igtor,igpol,g,ifa,i
							integer                    :: ind(refElPol%Ngauss1d),ind2(refElPol%Ngauss2d),nodes2(refelPol%Nnodes1d)
							! Volume computations
							real*8                     :: dvolu,isdir,isext,htor,dvolu1d
							real*8                     :: xy(refElPol%Ngauss2d,2),teg(refElTor%Ngauss1d)
       real*8                     :: J11(refElPol%Ngauss2d),J12(refElPol%Ngauss2d)
       real*8                     :: J21(refElPol%Ngauss2d),J22(refElPol%Ngauss2d)
       real*8                     :: detJ(refElPol%Ngauss2d)
       real*8                     :: iJ11(refElPol%Ngauss2d),iJ12(refElPol%Ngauss2d)
       real*8                     :: iJ21(refElPol%Ngauss2d),iJ22(refElPol%Ngauss2d)
       real*8                     :: Nxg(refElPol%Nnodes2D),Nyg(refElPol%Nnodes2D),Nx_ax(refElPol%Nnodes2D)
       real*8                     :: N1g(refElTor%Nnodes1D),N1xg(refElTor%Nnodes1D),N1xg_cart(refElTor%Nnodes1D),N2g(refElPol%Nnodes2D)
       real*8,allocatable         :: Ni(:),Nidvolu(:),Nr(:),Nrr(:),Nz(:),Nt(:),auxNNxy(:,:),NiNi(:,:),NNbb(:),NNxy(:)
       ! Faces computations
       integer*4,allocatable      :: ind_ff(:),ind_fe(:),ind_fg(:)       
       real*8,allocatable         :: aux(:,:),xyf(:,:)
       real*8,allocatable         :: thetafg(:),dsurf(:),n_g(:,:)       
							real*8                     :: xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
							real*8                     :: dsurfg,bn
							real*8                     :: t_g(refElPol%Ngauss1d,2)	
							real*8,pointer             :: Nfi(:,:),Nfg(:)
       ! Magnetic field, volume computations							       
       real*8                     :: b(refElPol%Ngauss2d*refElTor%Ngauss1d,3),b2d(refElPol%Ngauss2d,3)
       real*8                     :: divb(refElPol%Ngauss2d*refElTor%Ngauss1d),divb2d(refElPol%Ngauss2d)
       real*8                     :: drift(refElPol%Ngauss2d*refElTor%Ngauss1d,3),drift2d(refElPol%Ngauss2d,3)
       real*8                     :: fluxg(refElPol%Ngauss2d*refElTor%Ngauss1d)
       ! Magnetic field, faces computations
       real*8,allocatable         :: b_f(:,:)
       real*8                     :: b_f2d(refElPol%Ngauss1d,3)
       ! Routine specific part
							real*8                     :: force(refElPol%Ngauss2d*refElTor%Ngauss1d,phys%neq)
							real*8,allocatable         :: Mf_loc(:,:),Mass(:,:),Massf(:,:),floc(:,:),NNf(:,:)
							real*8, parameter          :: tol = 1e-12
       real*8,allocatable,dimension(:,:)    :: Bx,By,Bt,Px,Py,Pt,Qb_loc,Qbx,Qby,Qbt,Cn_loc,Cnx,Cny,Cnt,Ge,Lel,Dre,Drifte
       
       Ndim        = 3                                               ! Number of dimensions
       Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
       Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
       Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
       Np          = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
       Nfl         = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
       Nfp         = Np2D*2+refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
       Nf          = Mesh%Nfaces                                     ! Number of faces in the 2D mesh


       !***********************************
       !    Volume computation
       !***********************************
							ALLOCATE(Mass(Np,Np))
							ALLOCATE(floc(Np,Neq))
       ALLOCATE(Lel(Neq*Ndim*Np,Neq*Ndim*Np)) 
       ALLOCATE(Bx(Np,Np))
       ALLOCATE(By(Np,Np))
       ALLOCATE(Bt(Np,Np))
							Mass = 0.d0
							floc = 0.d0       
       Lel = 0.
       Bx = 0.; By = 0.; Bt = 0.
#ifndef MOVINGEQUILIBRIUM       
       ALLOCATE(Px(Np,Np))
       ALLOCATE(Py(Np,Np))
       ALLOCATE(Pt(Np,Np))
       Px = 0.; Py = 0.; Pt = 0.
#ifndef TEMPERATURE       
       ALLOCATE(Ge(Np,Np))
       ALLOCATE(Dre(Np,Np))
       ALLOCATE(Drifte(Np*Neq,Np*Neq))
       Ge = 0.
       Drifte = 0.
       Dre = 0.        
#endif       
#endif   
    
       ALLOCATE(Ni(Np))
							ALLOCATE(Nr(Np))
							ALLOCATE(Nrr(Np))
							ALLOCATE(Nz(Np))
							ALLOCATE(Nt(Np))
							ALLOCATE(Nidvolu(Np))
							ALLOCATE(NiNi(Np,Np))							
							ALLOCATE(NNbb(Np))
							ALLOCATE(NNxy(Ndim*Np))
							ALLOCATE(auxNNxy(Ndim,Np))
							      
						 ! toroidal element size
						 htor = tel(refElTor%Nnodes1d) - tel(1)
						 
							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)
							
							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)
														
							! Body force at the integration points
							CALL body_force(xy(:,1),xy(:,2),teg,force)							

!write(6,*) "xy: "
!call displayMatrix(xy)
!write(6,*) "teg: "
!call displayVector(teg)
!write(6,*) "force: "
!call displayMatrix(force)
!stop
!******************************************************
!          MOVING EQUILIBRIUM SECTION
!****************************************************** 
#ifndef MOVINGEQUILIBRIUM
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b2d = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
           divb2d = matmul(refElPol%N2D,phys%divb(Mesh%T(iel,:)))
#ifndef TEMPERATURE           
           drift2d = matmul(refElPol%N2D,phys%drift(Mesh%T(iel,:),:))
#endif           
          DO igtor = 1,refElTor%NGauss1D
             ind2 = (igtor-1)*refElPol%NGauss2D + (/(j,j=1,refElPol%NGauss2D)/)
             b(ind2,:) = b2d
             divb(ind2)   = divb2d
             drift(ind2,:) = drift2d
          END DO
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),teg,b,divb,drift)
							END IF
#ifndef TEMPERATURE							
							drift = phys%dfcoef*drift
#endif						
							!! Some sources for West cases
							IF (switch%testcase	.ge.51 .and. switch%testcase	.le.54) THEN						
							   fluxg = matmul(refElPol%N2D,phys%flux2D(Mesh%T(iel,:)))
							   DO g=1,refElPol%NGauss2D
													IF (switch%testcase==51) THEN
													   IF (fluxg(g).le.-0.88 .and. fluxg(g).ge.-0.90) THEN
													      !force(g,1) = 3.20119388718018e-05
                    force(g,1) = 4.782676673609557e-05
													   END IF
													ELSE IF (switch%testcase==52) THEN
													   IF (fluxg(g).le.-0.90 .and. fluxg(g).ge.-1.) THEN
													      force(g,1) = 9.45155008295538e-06
													   END IF													
													ELSE IF (switch%testcase==53) THEN
													   IF (fluxg(g).le.-0.90) THEN
													      force(g,1) = 7.24032211339971e-06
													   END IF													
													ELSE IF (switch%testcase==54) THEN
													   IF (fluxg(g).le.-1.03) THEN
													      force(g,1) = 0.000115575293741846
													   END IF
												 ELSE	IF (switch%testcase==55) THEN
													   IF (fluxg(g).le.-0.88 .and. fluxg(g).ge.-0.90) THEN
													      force(g,1) = 10
													   END IF													
													END IF
#ifdef TEMPERATURE
             force(g,3) = 18*force(g,1)
             force(g,4) = force(g,3)
#endif
							   END DO
							END IF
							!! end sources
#endif
!******************************************************
!          END MOVING EQUILIBRIUM SECTION
!******************************************************

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

														IF (switch%axisym) THEN
															  Nx_ax = Nxg+1./xy(igpol,1)*N2g
															  N1xg = N1xg_cart/xy(igpol,1) 
														ELSE
															  Nx_ax = Nxg
															  N1xg = N1xg_cart
														END IF
														
							       ! 3D shape functions 
							       Ni  = col(TensorProduct(N2g,N1g))      ! 3D shape function
							       Nrr = col(TensorProduct(Nx_ax,N1g))    ! 3D shape function, derivative in r for divergence 
							       Nr  = col(TensorProduct(Nxg,N1g))      ! 3D shape function, derivative in r for gradient 
							       Nz  = col(TensorProduct(Nyg,N1g))      ! 3D shape function, derivative in z
							       Nt  = col(TensorProduct(N2g,N1xg))     ! 3D shape function, derivative in t
							       
!							       auxNNxy(1,:) = Nr
!							       auxNNxy(2,:) = Nz
!							       auxNNxy(3,:) = Nt
!							       NNxy         = col(auxNNxy)
							       NNbb         = b(g,1)*Nr+b(g,2)*Nz+b(g,3)*Nt
							       							      
										    ! 3D integration weight
										    dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d
              IF (switch%axisym) THEN										
                 dvolu = dvolu*xy(igpol,1)
              END IF
              							       
							       ! Shape functions tensor product
							       Nidvolu      = Ni*dvolu			
							       NiNi = tensorProduct(Ni,Nidvolu)				      							       
								      
    										! Contribution of the current integration point to the elemental matrix
										    CALL TensorProductCumul(Bx,Nrr,Nidvolu) 
										    CALL TensorProductCumul(By,Nz,Nidvolu) 
										    CALL TensorProductCumul(Bt,Nt,Nidvolu) 

#ifndef MOVINGEQUILIBRIUM      
										    CALL TensorProductCumul(Px,Nr-b(g,1)*NNbb,Nidvolu) !Nb1  = Nb1  + b(g,1)*tensorProduct(NNbb,Nidvolu)
										    CALL TensorProductCumul(Py,Nz-b(g,2)*NNbb,Nidvolu) !Nb2  = Nb2  + b(g,2)*tensorProduct(NNbb,Nidvolu)
										    CALL TensorProductCumul(Pt,Nt-b(g,3)*NNbb,Nidvolu) !Nb3  = Nb3  + b(g,3)*tensorProduct(NNbb,Nidvolu)
#ifndef TEMPERATURE 
              Ge = Ge- phys%a*divb(igpol)*NiNi
              Dre = Dre + TensorProduct(Nidvolu,( Nr*drift(igpol,1)+ Nz*drift(g,2) +Nt*drift(g,3) ))
#endif          
#endif										
												  Mass = Mass + NiNi
												  floc = floc + tensorProduct(Nidvolu,force(g,:))
							    END DO
							END DO

							CALL expand_matrix(Mass,Neq,Ael)					
							fe = reshape(transpose(floc),(/Neq*Np/))
       CALL expand_matrix(Mass,Neq*Ndim,Lel)
       CALL invert_matrix(Lel,iLel)
       CALL expand_matrix_B(Bx,By,Bt,Bel)
#ifndef MOVINGEQUILIBRIUM        
       CALL expand_matrix_Bt(Px,Py,Pt,Pel)
#ifndef TEMPERATURE       
       CALL expand_matrix(Dre,Neq,Drifte)
#endif       
#endif


       DEALLOCATE(Ni,Nr,Nrr,Nz,Nt,Nidvolu,NiNi,NNbb,NNxy,auxNNxy)
	 

       !***********************************
  					! Faces computations
       !***********************************
       ! Loop in the faces of the element
							DO ifa = 1,refElPol%Nfaces+2

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
							      ind_fe = Np2d*(Np1dTor-1)*Neq+(/(i,i=1,Npifa*Neq)/)
							      ind_fg = Np2d*(Np1dTor-1)*Ndim*Neq + (/(i,i=1,Npifa*Ndim*Neq)/) 
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
							ALLOCATE(NNf(Npifa,Npifa))
							ALLOCATE(Mf_loc(Neq*Npifa,Neq*Npifa))
							ALLOCATE(Massf(Npifa,Npifa))
       ALLOCATE(Cn_loc(Neq*Ndim*Npifa,Neq*Npifa))
       ALLOCATE(Cnx(Npifa,Npifa))
       ALLOCATE(Cny(Npifa,Npifa))
       ALLOCATE(Cnt(Npifa,Npifa))
#ifndef MOVINGEQUILIBRIUM        
       ALLOCATE(Qb_loc(Neq*Npifa,Ndim*Neq*Npifa))
       ALLOCATE(Qbx(Npifa,Npifa))
       ALLOCATE(Qby(Npifa,Npifa))
       ALLOCATE(Qbt(Npifa,Npifa))
       ALLOCATE(b_f(NGaussPol*NGaussTor,3))		
#endif   
        
          
#ifndef MOVINGEQUILIBRIUM 
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
										Qbx = 0.; Qby = 0.; Qbt = 0.  
#endif										
          Cnx = 0.; Cny = 0.; Cnt = 0.
										Massf = 0.d0
          Mf_loc = 0.d0
#ifndef MOVINGEQUILIBRIUM          
          Qb_loc = 0.						
#endif 
          Cn_loc = 0.  

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
										
										      ! Contribution of the current integration point to the elemental matrix
                NNf = tensorProduct(Nfg,Nfg)*dsurfg

						          IF (numer%stab==1) THEN
																   ! Face mass matrix (to create the stabilization terms)
																   Massf = Massf + NNf
																END IF
																
                Cnx = Cnx + n_g(g,1)*NNf
                Cny = Cny + n_g(g,2)*NNf
                Cnt = Cnt + n_g(g,3)*NNf
#ifndef MOVINGEQUILIBRIUM
                bn = ( b_f(g,1)*n_g(g,1)+b_f(g,2)*n_g(g,2)+b_f(g,3)*n_g(g,3) )
                Qbx = Qbx + b_f(g,1)*NNf*bn
                Qby = Qby + b_f(g,2)*NNf*bn
                Qbt = Qbt + b_f(g,3)*NNf*bn
#endif										   
										   END DO ! Gauss points
										END DO
										
										! Assembly of the stabilization terms
          IF (numer%stab==1) THEN
										   CALL expand_matrix_stab(Massf,Neq,numer%tau,Mf_loc)
										   ! Elemental assembly
										   Del(ind_fe,ind_fe)  = Del(ind_fe,ind_fe) + Mf_loc										
										   Eel(ind_fe,ind_ff)  = Eel(ind_fe,ind_ff) + Mf_loc*(1-isdir)
										   Efel(ind_ff,ind_ff) = Efel(ind_ff,ind_ff) + Mf_loc*(1-isext)
										   Dfel(ind_ff,ind_fe) = Dfel(ind_ff,ind_fe) + Mf_loc*(1-isext)
										END IF


          CALL expand_matrix_B(Cnx,Cny,Cnt,Cn_loc)
          Cel(ind_fG,ind_ff) = Cel(ind_fG,ind_ff) + Cn_loc*(1-isdir)
          Lfel(ind_ff,ind_fG) = Lfel(ind_ff,ind_fG) + transpose(Cn_loc)*(1-isext)
#ifndef MOVINGEQUILIBRIUM
          CALL expand_matrix_Bt(Qbx,Qby,Qbt,Qb_loc)
          Qel(ind_fe,ind_fG) = Qel(ind_fe,ind_fG) + transpose(Cn_loc) - Qb_loc
          Qfel(ind_ff,ind_fG)  = Qfel(ind_ff,ind_fG) + Qb_loc*(1-isext)
#endif         
							DEALLOCATE(Mf_loc,Massf,Cn_loc,Cnx,Cny,Cnt)
#ifndef MOVINGEQUILIBRIUM        
       DEALLOCATE(Qb_loc,Qbx,Qby,Qbt,b_f)		
#endif  
       DEALLOCATE(n_g,dsurf,xyf,thetafg,ind_ff,ind_fe,ind_fg,NNf)
							END DO ! Faces

!call saveMatrix(Dfel,'Df_int')
!stop
      ! Multiply by the diffusion coefficient
      ! I am considering constant diffusion: 
      ! this should be done at Gauss point level 
      ! for variable diffusion
      CALL multiply_for_diffusion(Lfel)
#ifndef MOVINGEQUILIBRIUM      
      CALL multiply_for_diffusion(Qfel)
      CALL multiply_for_diffusion(Pel)
      CALL multiply_for_diffusion(Qel)
      
#ifndef TEMPERATURE 
      ! Ge and Drifte are filled only in the isothermal case.      
      CALL expand_matrix_G(Ge,Gel)    
      Gel = Gel + Drifte
      
      ! Deallocate stuff
      DEALLOCATE(Ge,Drifte,Dre)
#endif       
      DEALLOCATE(Px,Py,Pt)
#endif      
      DEALLOCATE(Bx,By,Bt,Lel)
      DEALLOCATE(Mass,floc)
					END SUBROUTINE elemental_matrices				
					
													
!*******************************************
!           AUXILIARY ROUTINES
!*******************************************

									
								!*******************************************
								! Expand matrix
								!*******************************************									
								SUBROUTINE expand_matrix(A,n,B)
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
								END SUBROUTINE expand_matrix
								
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
																
						
								
              
								!*******************************************
								! Expand matrix G 
								!*******************************************									
								SUBROUTINE expand_matrix_G(Gl,G)
										real*8,intent(in)    :: Gl(:,:)
										real*8,intent(inout) :: G(2*size(Gl,1),2*size(Gl,2))
										integer*4            :: i,j,n,m
										
										n = size(Gl,1)
										m = size(Gl,2)
				      DO i=1,n
				         DO j=1,m
				           G((i-1)*2+2,(j-1)*2+1) = G((i-1)*2+2,(j-1)*2+1) + Gl(i,j)
				         END DO
				      END DO
								END SUBROUTINE expand_matrix_G	
															
								
     
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



END SUBROUTINE HDG_precalculatedmatrices

#else


!***********************************************************************
! 
!                            VERSION 2D 
! 
!***********************************************************************
  integer*4             :: iel,i,j,Ndim,Neq,Nel,Np,Nfp,Nf,ierr,ith
  integer*4,allocatable :: ind_loc(:,:),perm(:)
  logical               :: F_dir_el(refElPol%Nfaces)
  real*8                :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
  real*8,allocatable    :: Ael(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:), fe(:)
  real*8,allocatable    :: Bel(:,:),Cel(:,:),Gel(:,:),iLel(:,:),Pel(:,:),Qel(:,:), Lfel(:,:),Qfel(:,:)

  
  integer :: OMP_GET_THREAD_NUM
  ith=0
  
  IF (MPIvar%glob_id.eq.0) THEN
					IF (utils%printint>0) THEN
						WRITE(6,*) '*************************************************'
						WRITE(6,*) '*           PRECALCULATED MATRICES              *'
						WRITE(6,*) '*************************************************' 
					END IF	   
		END IF
		
  Ndim        = Mesh%ndim
  Neq         = phys%Neq
  Nel         = Mesh%Nelems
  Np          = refElPol%Nnodes2D
  Nfp         = refElPol%Nfacenodes
  Nf          = refElPol%Nfaces

 
  ALLOCATE(ind_loc(Nf,Neq*Nfp))
  ind_loc = 0
  DO i = 1,Nf
     DO j = 1,Neq*Nfp
        ind_loc(i,j) = Neq*Nfp*(i-1) + j
     END DO
  END DO

  ! Set perm for flipping faces
  ALLOCATE(perm(Neq*Nfp))
  perm = 0
  
  CALL set_permutations(Neq*Nfp,Neq,perm)
 

!*****************
! Loop in elements
!*****************  

!$OMP PARALLEL DEFAULT(SHARED) &
#ifndef MOVINGEQUILIBRIUM
!$OMP PRIVATE(iel,Xe,F_dir_el,j,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Lfel,Qfel,Gel)
#else
!$OMP PRIVATE(iel,Xe,F_dir_el,j,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif

  ALLOCATE(Ael(Neq*Np,Neq*Np))
  ALLOCATE(Del(Neq*Np,Neq*Np))
  ALLOCATE(Eel(Neq*Np,Neq*Nf*Nfp))
  ALLOCATE(Dfel(Neq*Nf*Nfp,Neq*Np))
  ALLOCATE(Efel(Neq*Nf*Nfp,Neq*Nf*Nfp))
  ALLOCATE(fe(Neq*Np))
  ALLOCATE(Bel(Neq*Ndim*Np,Neq*Np))
  ALLOCATE(Cel(Neq*Ndim*Np,Nf*Neq*Nfp))
  ALLOCATE(iLel(Neq*Ndim*Np,Neq*Ndim*Np))
  ALLOCATE(Lfel(Neq*Nf*Nfp,Neq*Ndim*Np))  
#ifndef MOVINGEQUILIBRIUM
  ALLOCATE(Gel(Neq*Np,Neq*Np))
  ALLOCATE(Pel(Neq*Np,Ndim*Neq*Np))
  ALLOCATE(Qel(Neq*Np,Ndim*Neq*Np))
  ALLOCATE(Qfel(Neq*Nf*Nfp,Neq*Ndim*Np))
#endif

!$OMP DO SCHEDULE(STATIC)
  DO iel = 1,Nel 
     Ael    = 0.d0     
     Del    = 0.d0
     Eel    = 0.d0
     Dfel   = 0.d0
     Efel   = 0.d0
     fe     = 0.d0     
     Bel  = 0.
     Cel  = 0.
     iLel = 0.
     Lfel = 0.
#ifndef MOVINGEQUILIBRIUM
     Gel  = 0.    
     Pel  = 0.
     Qel  = 0.     
     Qfel = 0.     
#endif
     ! Coordinates of the nodes of the element
     Xe = Mesh%X(Mesh%T(iel,:),:)

!ith =  omp_get_thread_num()
!write(6,*) "TEST 1 - iel:" , iel ," thread: ", ith
     ! Dirichlet boundary conditions in the current element
     F_dir_el = Mesh%Fdir(iel,:)
 
     ! Compute the matrices for the element
#ifndef MOVINGEQUILIBRIUM     
     CALL elemental_matrices(iel,Xe,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Gel,Lfel,Qfel)
#else
     CALL elemental_matrices(iel,Xe,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif     
     ! Flip the faces for the elements for which the face is already been counted
     DO j = 1,Nf
        IF (Mesh%flipface(iel,j)) THEN
           Dfel(ind_loc(j,:),:) = Dfel(ind_loc(j,perm(:)),:)
           Efel(ind_loc(j,:),:) = Efel(ind_loc(j,perm(:)),:)
           Efel(:,ind_loc(j,:)) = Efel(:,ind_loc(j,perm(:)))
           Lfel(ind_loc(j,:),:) = Lfel(ind_loc(j,perm),:)
#ifndef MOVINGEQUILIBRIUM
           Qfel(ind_loc(j,:),:) = Qfel(ind_loc(j,perm),:)           
#endif           
        END if
     END DO

     ! Store the matrices for all the elements
     IF (.not.switch%steady) THEN
        elMat%M(:,:,iel)  = Ael
     END IF
     elMat%D(:,:,iel)  = Del
     elMat%E(:,:,iel)  = Eel
     elMat%Df(:,:,iel) = Dfel
     elMat%Ef(:,:,iel) = Efel
     elMat%S(:,iel)    = fe
     elMat%B(:,:,iel)  = Bel
     elMat%C(:,:,iel)  = Cel
     elMat%iL(:,:,iel) = iLel
     elMat%Lf(:,:,iel)   = Lfel
#ifndef MOVINGEQUILIBRIUM
     elMat%G(:,:,iel)  = Gel 
     elMat%P(:,:,iel)  = Pel-Qel
     elMat%Qf(:,:,iel)   = Qfel
#endif     
  END DO
!$OMP END DO

  DEALLOCATE(Ael,Del,Eel,Dfel,Efel,fe)
  DEALLOCATE(Bel,Cel,iLel,Lfel)
#ifndef MOVINGEQUILIBRIUM
  DEALLOCATE(Pel,Qel,Qfel,Gel)
#endif  

!$OMP END PARALLEL

  DEALLOCATE(ind_loc,perm)

 
  IF (MPIvar%glob_id.eq.0) THEN
					IF (utils%printint>0) THEN
						  WRITE(6,*) "Done!"
					END IF
  END IF
  
  
  CONTAINS
  

!*****************************************
! Elemental matrices computation: 2D case
!*****************************************
#ifndef MOVINGEQUILIBRIUM
					SUBROUTINE elemental_matrices(iel,Xe,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Pel,Qel,Gel,Lfel,Qfel )
#else
     SUBROUTINE elemental_matrices(iel,Xe,F_dir_el,Ael,Del,Eel,Dfel,Efel,fe,Bel,Cel,iLel,Lfel)
#endif					
							integer*4,intent(IN)  :: iel
       real*8,intent(IN)     :: Xe(:,:)							
       logical,intent(IN)    :: F_dir_el(:)
       real*8,intent(INOUT)  :: Ael(:,:),Del(:,:),Eel(:,:),Dfel(:,:),Efel(:,:), fe(:)     
#ifndef MOVINGEQUILIBRIUM
       real*8,intent(INOUT)  :: Bel(:,:),Cel(:,:),iLel(:,:),Pel(:,:),Qel(:,:),Lfel(:,:),Qfel(:,:),Gel(:,:) 
#else
       real*8,intent(INOUT)  :: Bel(:,:),Cel(:,:),iLel(:,:),Lfel(:,:) 
#endif       
							integer*4             :: g,NGauss,ifa,Np,Neq,Ndim,Nfp,Nf,i,j
							real*8                :: detJ,dvolu,dline,isdir,isext,xyDerNorm_g
							real*8                :: xy(refElPol%Ngauss2d,2)
							real*8                :: xyf(refElPol%Ngauss1d,2)
							real*8                :: xyDer(refElPol%Ngauss1d,2)
							real*8                :: force(refElPol%Ngauss2d,phys%neq)
							real*8                :: Jacob(2,2)
							integer*4             :: ind_ff(phys%neq*refElPol%Nfacenodes),ind_fe(phys%neq*refElPol%Nfacenodes)
							real*8,allocatable    :: Mf_loc(:,:),Mass(:,:),Massf(:,:),floc(:,:)
							real*8, parameter     :: tol = 1e-12
       integer*4             :: ind_fG(phys%neq*2*refElPol%Nfacenodes)
       real*8,dimension(refElPol%Nnodes2D)  :: Nxg,Nyg,Nx_ax
       real*8                :: invJ(2,2),t_g(2),n_g(2)
       real*8                :: b(refElPol%Ngauss2d,2),divb(refElPol%Ngauss2d),drift(refElPol%Ngauss2d,2)
       real*8                :: b_f(refElPol%Ngauss1d,2)
       real*8                :: Nbn(refElPol%Nfacenodes)
       real*8                :: fluxg(refElPol%Ngauss2d)
       real*8,allocatable,dimension(:,:)    :: Bx,By,Px,Py,Qb_loc,Qbx,Qby,Cn_loc,Cnx,Cny,Ge,Lel,Dre,Drifte

       integer*4             :: ind_ass(refElPol%Nnodes2D)
       !***********************************
       !    Volume computation
       !***********************************
							Ndim        = Mesh%ndim
							Neq         = phys%Neq
							Np          = refElPol%Nnodes2D
							Nfp         = refElPol%Nfacenodes
							Nf          = refElPol%Nfaces
							       
       ind_ass = (/ (i, i = 0, Neq*(Np-1),Neq) /)							       
       
							ALLOCATE(Mass(Np,Np))
							ALLOCATE(floc(Np,Neq))
							Mass = 0.d0
							floc = 0.d0
       ALLOCATE(Lel(Neq*Ndim*Np,Neq*Ndim*Np)) 
       ALLOCATE(Bx(Np,Np))
       ALLOCATE(By(Np,Np))
#ifndef MOVINGEQUILIBRIUM       
       ALLOCATE(Px(Np,Np))
       ALLOCATE(Py(Np,Np))
#ifndef TEMPERATURE       
       ALLOCATE(Ge(Np,Np))
       ALLOCATE(Dre(Np,Np))
       ALLOCATE(Drifte(Np*Neq,Np*Neq))
       Ge = 0.
       Drifte = 0.
       Dre = 0.        
#endif       
       Px = 0.; Py = 0.      
#endif       
       Lel = 0.
       Bx = 0.; By = 0.

							! Gauss points position
							xy = matmul(refElPol%N2D,Xe)
							
							! Body force at the integration points
							CALL body_force(xy(:,1),xy(:,2),force)							

#ifndef MOVINGEQUILIBRIUM
       IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
           b = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
           divb = matmul(refElPol%N2D,phys%divb(Mesh%T(iel,:)))
#ifndef TEMPERATURE           
           drift = matmul(refElPol%N2D,phys%drift(Mesh%T(iel,:),:))
#endif           
       ELSE
							   CALL defineMagneticField(xy(:,1),xy(:,2),b,divb,drift)
							END IF
#ifndef TEMPERATURE							
							drift = phys%dfcoef*drift
#endif						
							!! Some sources for West cases
							IF (switch%testcase	.ge.51 .and. switch%testcase	.le.54) THEN						
							   fluxg = matmul(refElPol%N2D,phys%flux2D(Mesh%T(iel,:)))
							   DO g=1,refElPol%NGauss2D
													IF (switch%testcase==51) THEN
													   IF (fluxg(g).le.-0.88 .and. fluxg(g).ge.-0.90) THEN
													      !force(g,1) = 3.20119388718018e-05
                    force(g,1) = 4.782676673609557e-05
													   END IF
													ELSE IF (switch%testcase==52) THEN
													   IF (fluxg(g).le.-0.90 .and. fluxg(g).ge.-1.) THEN
													      force(g,1) = 9.45155008295538e-06
													   END IF													
													ELSE IF (switch%testcase==53) THEN
													   IF (fluxg(g).le.-0.90) THEN
													      force(g,1) = 7.24032211339971e-06
													   END IF													
													ELSE IF (switch%testcase==54) THEN
													   IF (fluxg(g).le.-1.03) THEN
													      force(g,1) = 0.000115575293741846
													   END IF
												 ELSE	IF (switch%testcase==55) THEN
													   IF (fluxg(g).le.-0.88 .and. fluxg(g).ge.-0.90) THEN
													      force(g,1) = 10
													   END IF													
													END IF
#ifdef TEMPERATURE
             force(g,3) = 18*force(g,1)
             force(g,4) = force(g,3)
#endif
							   END DO
							END IF
							!! end sources
#endif
							! Loop in 2D Gauss points
       Ngauss = refElPol%NGauss2D
							DO g = 1,NGauss
 							
										! Jacobian										
										Jacob = 0.d0
										Jacob(1,1) = dot_product(refElPol%Nxi2D(g,:),Xe(:,1))
										Jacob(1,2) = dot_product(refElPol%Nxi2D(g,:),Xe(:,2))
										Jacob(2,1) = dot_product(refElPol%Neta2D(g,:),Xe(:,1))
										Jacob(2,2) = dot_product(refElPol%Neta2D(g,:),Xe(:,2))
										detJ = Jacob(1,1)*Jacob(2,2) - Jacob(1,2)*Jacob(2,1)

										IF (detJ < tol) THEN
										   error stop "Negative jacobian"
										END if

										! x and y derivatives of the shape functions
										call invert_matrix(Jacob,invJ)
							   Nxg = invJ(1,1)*refElPol%Nxi2D(g,:) + invJ(1,2)*refElPol%Neta2D(g,:)
							   Nyg = invJ(2,1)*refElPol%Nxi2D(g,:) + invJ(2,2)*refElPol%Neta2D(g,:)
							   IF (switch%axisym) THEN
							      Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
							   ELSE
							      Nx_ax = Nxg
							   END IF
										! Integration weight
										dvolu = refElPol%gauss_weights2D(g)*detJ
          IF (switch%axisym) THEN										
             dvolu = dvolu*xy(g,1)
          END IF
										! Contribution of the current integration point to the elemental matrix
!										DO i=1,Ndim
!										   if (i==1) then
!										   CALL TensorProductCumul(Bel(ind_ass+i,:),Nx_ax*dvolu,refElPol%N2D(g,:))
!										   else
!										   CALL TensorProductCumul(Bel(ind_ass+i,:),Nyg*dvolu,refElPol%N2D(g,:))
!										   endif
!										END DO
!										CALL TensorProductCumul(Bx,Nx_ax*dvolu,refElPol%N2D(g,:))
!										CALL TensorProductCumul(By,Nyg*dvolu,refElPol%N2D(g,:))  
          Bx = Bx + TensorProduct(Nx_ax,refElPol%N2D(g,:))*dvolu       
          By = By + TensorProduct(Nyg,refElPol%N2D(g,:))*dvolu
#ifndef MOVINGEQUILIBRIUM           
          Px = Px + TensorProduct(Nxg-b(g,1)*(b(g,1)*Nxg+b(g,2)*Nyg),refElPol%N2D(g,:))*dvolu
          Py = Py + TensorProduct(Nyg-b(g,2)*(b(g,1)*Nxg+b(g,2)*Nyg),refElPol%N2D(g,:))*dvolu
#ifndef TEMPERATURE          
          Ge = Ge - phys%a*divb(g)*TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))*dvolu       
          Dre = Dre + TensorProduct(refElPol%N2D(g,:),( Nxg*drift(g,1)+ Nyg*drift(g,2) ))*dvolu
#endif          
#endif										
				      Mass = Mass + dvolu*TensorProduct(refElPol%N2D(g,:),refElPol%N2D(g,:))
				      floc = floc + tensorProduct(refElPol%N2D(g,:),force(g,:))*dvolu
							END DO

							CALL expand_matrix(Mass,Neq,Ael)					
							fe = reshape(transpose(floc),(/Neq*Np/))
       CALL expand_matrix(Mass,Neq*Ndim,Lel)
       CALL invert_matrix(Lel,iLel)
       CALL expand_matrix_B(Bx,By,Bel)
#ifndef MOVINGEQUILIBRIUM        
       CALL expand_matrix_Bt(Px,Py,Pel)
#ifndef TEMPERATURE       
       CALL expand_matrix(Dre,Neq,Drifte)
#endif       
#endif

       !***********************************
  					! Faces computations
       !***********************************
							NGauss = refElPol%Ngauss1D
							ALLOCATE(Mf_loc(Neq*Nfp,Neq*Nfp))
							ALLOCATE(Massf(Nfp,Nfp))
       ALLOCATE(Cn_loc(Neq*Ndim*Nfp,Neq*Nfp))
       ALLOCATE(Cnx(Nfp,Nfp))
       ALLOCATE(Cny(Nfp,Nfp))
#ifndef MOVINGEQUILIBRIUM        
       ALLOCATE(Qb_loc(Neq*Nfp,Ndim*Neq*Nfp))
       ALLOCATE(Qbx(Nfp,Nfp))
       ALLOCATE(Qby(Nfp,Nfp))	
#endif   

       ! Loop in the faces of the element
							DO ifa = 1,Nf
          Mf_loc = 0.d0
#ifndef MOVINGEQUILIBRIUM          
       Qb_loc = 0.						
#endif 
       Cn_loc = 0.
          ! Set isdir for Dirichlet faces and 
          ! isext for exterior faces
          isdir = 0.
          isext = 0.
          IF (F_dir_el(ifa)) isdir = 1.
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
          ind_fG = reshape(tensorSumInt( (/ (i, i = 1, neq*ndim) /),neq*ndim*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*ndim*Nfp/))
          ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Nfp/))
          ind_ff = (ifa-1)*neq*Nfp + (/ (i, i = 1, neq*Nfp) /)
										xyf    = matmul(refElPol%N1D,Xe(refElPol%face_nodes(ifa,:),:))
          xyDer  = matmul(refElPol%Nxi1D,Xe(refElPol%face_nodes(ifa,:),:))
          
#ifndef MOVINGEQUILIBRIUM 
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
						        b_f = matmul(refElPol%N1D,phys%b(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:))
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),b_f)  
										END IF				
										Qbx = 0.; Qby = 0.  
#endif										
          Cnx = 0.; Cny = 0.
										Massf = 0.d0
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
             IF (numer%stab==1) THEN
										      ! Face mass matrix (to create the stabilization terms)
										      Massf = Massf + dline*tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))
										   END IF
             Cnx = Cnx + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(1)*dline
             Cny = Cny + tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(2)*dline
#ifndef MOVINGEQUILIBRIUM
             Nbn = refElPol%N1D(g,:)*( b_f(g,1)*n_g(1)+b_f(g,2)*n_g(2) )
             Qbx = Qbx + b_f(g,1)*tensorProduct(refElPol%N1D(g,:),Nbn)*dline
             Qby = Qby + b_f(g,2)*tensorProduct(refElPol%N1D(g,:),Nbn)*dline     
#endif										   
										END DO ! Gauss points
										
										! Assembly of the stabilization terms
          IF (numer%stab==1) THEN
										   CALL expand_matrix_stab(Massf,Neq,numer%tau,Mf_loc)
										   ! Elemental assembly
										   Del(ind_fe,ind_fe)  = Del(ind_fe,ind_fe) + Mf_loc										
										   Eel(ind_fe,ind_ff)  = Eel(ind_fe,ind_ff) + Mf_loc*(1-isdir)
										   Efel(ind_ff,ind_ff) = Efel(ind_ff,ind_ff) + Mf_loc*(1-isext)
										   Dfel(ind_ff,ind_fe) = Dfel(ind_ff,ind_fe) + Mf_loc*(1-isext)
										END IF

          CALL expand_matrix_B(Cnx,Cny,Cn_loc)
          Cel(ind_fG,ind_ff) = Cel(ind_fG,ind_ff) + Cn_loc*(1-isdir)
          Lfel(ind_ff,ind_fG) = Lfel(ind_ff,ind_fG) + transpose(Cn_loc)*(1-isext)
#ifndef MOVINGEQUILIBRIUM
          CALL expand_matrix_Bt(Qbx,Qby,Qb_loc)
          Qel(ind_fe,ind_fG) = Qel(ind_fe,ind_fG) + transpose(Cn_loc) - Qb_loc
          Qfel(ind_ff,ind_fG)  = Qfel(ind_ff,ind_fG) + Qb_loc*(1-isext)
#endif          
							END DO ! Faces


     DEALLOCATE(Mf_loc,Mass,Massf,floc)		
      CALL multiply_for_diffusion(Lfel)
#ifndef MOVINGEQUILIBRIUM      
      CALL multiply_for_diffusion(Qfel)
      CALL multiply_for_diffusion(Pel)
      CALL multiply_for_diffusion(Qel)
#ifndef TEMPERATURE 
      ! Ge and Drifte are filled only in the isothermal case.      
      CALL expand_matrix_G(Ge,Gel)    
      Gel = Gel + Drifte
      DEALLOCATE(Ge,Drifte,Dre)
#endif       
      DEALLOCATE(Px,Py,Qb_loc,Qbx,Qby)
#endif      
      DEALLOCATE(Bx,By,Lel)
      DEALLOCATE(Cn_loc,Cnx,Cny)
					END SUBROUTINE elemental_matrices				
					
													
!*******************************************
!           AUXILIARY ROUTINES
!*******************************************

									
								!*******************************************
								! Expand matrix
								!*******************************************									
								SUBROUTINE expand_matrix(A,n,B)
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
								END SUBROUTINE expand_matrix
								
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
								! Expand matrix B
								!*******************************************									
								SUBROUTINE expand_matrix_B(Bx,By,B)
										real*8,intent(in)    :: Bx(:,:),By(:,:)
										real*8,intent(out)   :: B(Neq*Ndim*size(Bx,1),Neq*size(Bx,2))
										integer*4            :: i,j,k,n,m
										
										n = size(Bx,1)
										m = size(Bx,2)
										B = 0.d0
										DO j= 1,m
													DO i = 1,n
													   DO k = 1,Neq
 																		B((i-1)*Neq*Ndim+1+(k-1)*Ndim,(j-1)*Neq+k) = Bx(i,j)
 																		B((i-1)*Neq*Ndim+2+(k-1)*Ndim,(j-1)*Neq+k) = By(i,j)
!																B((i-1)*4+1,(j-1)*2+1) = Bx(i,j)
!																B((i-1)*4+2,(j-1)*2+1) = By(i,j)
!																B((i-1)*4+3,(j-1)*2+2) = Bx(i,j)
!																B((i-1)*4+4,(j-1)*2+2) = By(i,j)
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_B								
								


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
!																B((i-1)*2+1,(j-1)*4+1) = Bx(i,j)
!																B((i-1)*2+1,(j-1)*4+2) = By(i,j)
!																B((i-1)*2+2,(j-1)*4+3) = Bx(i,j)
!																B((i-1)*2+2,(j-1)*4+4) = By(i,j)
																END DO
													END DO
										END DO
								END SUBROUTINE expand_matrix_Bt			
																
								
!! OLD VERSION OF EXPAND MATRIX B AND Bt 								
!								!*******************************************
!								! Expand matrix B
!								!*******************************************									
!								SUBROUTINE expand_matrix_B(Bx,By,B)
!										real*8,intent(in)    :: Bx(:,:),By(:,:)
!										real*8,intent(out)   :: B(4*size(Bx,1),2*size(Bx,2))
!										integer*4            :: i,j,n,m
!										
!										n = size(Bx,1)
!										m = size(Bx,2)
!										B = 0.d0
!										DO j= 1,m
!													DO i = 1,n
!																B((i-1)*4+1,(j-1)*2+1) = Bx(i,j)
!																B((i-1)*4+2,(j-1)*2+1) = By(i,j)
!																B((i-1)*4+3,(j-1)*2+2) = Bx(i,j)
!																B((i-1)*4+4,(j-1)*2+2) = By(i,j)
!													END DO
!										END DO
!								END SUBROUTINE expand_matrix_B								
!								


!								!*******************************************
!								! Expand matrix B transpose
!								!*******************************************									
!								SUBROUTINE expand_matrix_Bt(Bx,By,B)
!										real*8,intent(in)    :: Bx(:,:),By(:,:)
!										real*8,intent(out)   :: B(2*size(Bx,1),4*size(Bx,2))
!										integer*4            :: i,j,n,m
!										
!										n = size(Bx,1)
!										m = size(Bx,2)
!										B = 0.d0
!										DO j = 1,m
!													DO i = 1,n
!																B((i-1)*2+1,(j-1)*4+1) = Bx(i,j)
!																B((i-1)*2+1,(j-1)*4+2) = By(i,j)
!																B((i-1)*2+2,(j-1)*4+3) = Bx(i,j)
!																B((i-1)*2+2,(j-1)*4+4) = By(i,j)
!													END DO
!										END DO
!								END SUBROUTINE expand_matrix_Bt					
								
								

								!*******************************************
								! Expand matrix G
								!*******************************************									
								SUBROUTINE expand_matrix_G(Gl,G)
										real*8,intent(in)    :: Gl(:,:)
										real*8,intent(inout) :: G(2*size(Gl,1),2*size(Gl,2))
										integer*4            :: i,j,n,m
										
										n = size(Gl,1)
										m = size(Gl,2)
				      DO i=1,n
				         DO j=1,m
				           G((i-1)*2+2,(j-1)*2+1) = G((i-1)*2+2,(j-1)*2+1) + Gl(i,j)
				         END DO
				      END DO
								END SUBROUTINE expand_matrix_G	
															
								
     
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
							 temp = reshape( (/ (i, i = 1, n) /), (/ m, n/m /) )
							 DO i = 1,n/m
							    templr(:,i) = temp(:,n/m-i+1)
							 END DO
							 perm = reshape(templr,(/n/))
							 END SUBROUTINE set_permutations     
     



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



END SUBROUTINE HDG_precalculatedmatrices

#endif


