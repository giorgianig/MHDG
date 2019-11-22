SUBROUTINE HDG_mapping(um1,um2)
  USE globals
  USE LinearAlgebra
  USE analytical
  USE printUtils
  USE MPI_OMP
  
  IMPLICIT NONE
  real*8,intent(in)  :: um1(size(sol%u))
  real*8,intent(in),optional  :: um2(size(sol%u))
  integer*4          :: Ndim,N2d,neq,Np,Nfp,Nfe,i,j,ifa,iel,ntorloc,iel3,Nfg
  integer*4          :: Np2D,Np1Dpol,Np1Dtor,Nfl,itor
#ifdef TOR3D  
  integer*4          :: ind(refElPol%Nfaces,phys%neq*refElTor%Nfl) 
  integer*4          :: perm(refElTor%Nfl*phys%Neq)
#else
  integer*4          :: ind(refElPol%Nfaces,phys%neq*Mesh%Nnodesperface)
  integer*4          :: perm(refElPol%Nfacenodes*phys%Neq)
#endif  
  integer*4,allocatable :: ind_u(:)
  real*8,allocatable :: u(:)
  real*8,allocatable :: u_2(:)
  real*8,pointer     :: M(:,:),Cv(:,:),D(:,:),E(:,:)
  real*8,pointer     :: H(:,:),Hd(:)
  real*8,allocatable :: S(:),Ed(:)
  real*8,allocatable :: UU(:,:),U0(:)
  real*8,allocatable :: Mu(:,:),Mt(:,:),M0(:)
  real*8,pointer     :: B(:,:),C(:,:),iL(:,:),Cd(:),G(:,:)
  real*8,allocatable :: PinvL(:,:),P(:,:)
  real*8,allocatable :: LL(:,:),L0(:)
#ifdef TEMPERATURE
  real*8,pointer     ::TQ(:,:),Tf(:),Td(:)
#endif  


#ifdef TOR3D		
  Ndim        = 3                                               ! Number of dimensions
  neq         = phys%neq                                        ! Number of equations
  N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np          = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
  Nfl         = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
  Nfg         = Np2D*2+refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
  Nfe         = refElPol%Nfaces                                 ! Number of faces in the 2D element
#else  
  Ndim        = 2
  N2D         = Mesh%Nelems
  Np          = Mesh%Nnodesperelem
  neq         = phys%neq
  Nfe         = refElPol%Nfaces
  Nfg         = Mesh%Nnodesperface*Nfe
  Nfp         = Mesh%Nnodesperface
#endif


#ifdef TOR3D
  ind = 0
  DO i = 1,Nfe
     DO j = 1,Nfl*Neq
        ind(i,j) = refElPol%Nnodes2D*Neq + Neq*Nfl*(i-1) + j
     END DO
  END DO
    
  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol,Np1Dtor,Neq,perm)   
#else
		! local assembly indexes
		ind = 0
		DO ifa = 1,Nfe
				  ind(ifa,:) = neq*Nfp*(ifa-1)+ ((/(i,i=1,neq*Nfp)/))
		END DO

  ! Set perm for flipping faces
  perm = 0 
  CALL set_permutations(Neq*Nfp,Neq,perm)
#endif  

#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%ntor.gt.1) THEN
       ntorloc = numer%ntor/MPIvar%ntor+1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
#endif


!$OMP PARALLEL DEFAULT(SHARED) &
#ifndef TEMPERATURE
!$OMP PRIVATE(iel,itor,iel3,ind_u,i,u,u_2,M,Cv,D,E,H,S,Ed,Hd,iL,B,C,P,G,Cd,Mu,Mt,M0,PinvL,LL,L0,UU,U0,ifa)   
#else
!$OMP  PRIVATE(iel,itor,iel3,ind_u,i,u,u_2,M,Cv,D,E,H,S,Ed,Hd,iL,B,C,P,G,Cd,Mu,Mt,M0,PinvL,LL,L0,UU,U0,ifa,TQ,Tf,Td)   
#endif
  ALLOCATE(u(Neq*Np))
  ALLOCATE(ind_u(Neq*Np))
  ALLOCATE(UU(Neq*Np,Neq*Nfg))
  ALLOCATE(U0(Neq*Np))
  ALLOCATE(Mu(Neq*Np,Neq*Np))
  ALLOCATE(Mt(Neq*Np,Neq*Nfg)) 
  ALLOCATE(M0(Neq*Np))
  ALLOCATE(S(Neq*Np))
  ALLOCATE(Ed(Neq*Np))
 IF (time%tis.eq.2 .and. time%ik.gt.1) THEN
    ALLOCATE(u_2(phys%neq*Mesh%Nnodesperelem))
 END IF
  ALLOCATE(P(Neq*Np,Ndim*Neq*Np))
  ALLOCATE(PinvL(Neq*Np,Neq*Np*Ndim))
  ALLOCATE(LL(Neq*Np*Ndim,Neq*Nfg))
  ALLOCATE(L0(Neq*Np*Ndim))


!$OMP DO SCHEDULE(STATIC)
#ifdef TOR3D
  DO itor = 1,ntorloc
#endif     
     DO iel = 1,N2d
     
#ifdef TOR3D     
        iel3 = (itor-1)*N2d+iel
#else
        iel3 = iel
#endif

     ! Solution index
     ind_u = (iel3-1)*neq*Np + ((/ (i,i=1,neq*Np) /))

     ! Solution in the element
     u = um1(ind_u)
     
					IF (time%tis.eq.2 .and. time%ik.gt.1) THEN
								u_2 = um2(ind_u)
					END IF          

     ! Elemental matrices
     M  => elMat%M(:,:,iel3)
     Cv => elMat%Cv(:,:,iel3)
     D  => elMat%D(:,:,iel3)
     E  => elMat%E(:,:,iel3)
     H  => elMat%H(:,:,iel3)     
     Hd => elMat%Hdir(:,iel3)
     S  =  elMat%S(:,iel3)
IF (Mesh%ndir.gt.0) THEN     
     Ed =  elMat%Edir(:,iel3)
ELSE
     Ed = 0.
ENDIF     
     iL => elMat%iL(:,:,iel3)
     B  => elMat%B(:,:,iel3)
     C  => elMat%C(:,:,iel3)
     Cd => elMat%Cdir(:,iel3)
     G  => elMat%G(:,:,iel3)
     P  =  elMat%P(:,:,iel3)   
#ifdef TEMPERATURE
     TQ  => elMat%TQ(:,:,iel3)
     Tf  => elMat%Tf(:,iel3)
     Td  => elMat%Thdir(:,iel3)
     S = S+Tf
     P = P+TQ
IF (Mesh%ndir.gt.0) THEN     
     Ed = Ed+Td      
ENDIF             
#endif

     ! Use limiting 
     IF (switch%limrho.gt.0) THEN
						  IF (switch%limrho.eq.1) THEN
						     IF (Mesh%flag_elems_rho(iel3).ne.0) THEN
						        S = S+elMat%S_lrho(:,Mesh%flag_elems_rho(iel3))
						     END IF
						  ELSE IF (switch%limrho.eq.2) THEN
						     IF (Mesh%flag_elems_rho(iel3).ne.0) THEN
						        P = P+elMat%P_lrho(:,:,Mesh%flag_elems_rho(iel3))
						     END IF
						  ELSE 
						     IF (Mesh%flag_elems_rho(iel3).ne.0) THEN
						        S = S+elMat%S_lrho(:,Mesh%flag_elems_rho(iel3))
						        P = P+elMat%P_lrho(:,:,Mesh%flag_elems_rho(iel3))
						     END IF     
						  END IF
     END IF      
     ! Use shock capturing
     IF (switch%shockcp.gt.0) THEN
			     IF (Mesh%flag_elems_sc(iel3).ne.0) THEN
			        P = P+elMat%P_sc(:,:,Mesh%flag_elems_sc(iel3))
			     END IF     
     END IF     
     Mu = M/time%dt-Cv+D
     Mu = Mu+G
     IF (time%tis.eq.2 .and. time%ik.gt.1) THEN
        Mu = Mu+0.5*M/time%dt
     END IF
     Mt = E-H
     IF (Mesh%ndir.gt.0) THEN
        M0 = matmul(M/time%dt,u)+S+Ed-Hd
     ELSE
        M0 = matmul(M/time%dt,u)+S
     END IF
     IF (time%tis.eq.2 .and. time%ik.gt.1) THEN
        M0 = M0 + matmul(M/time%dt,u)-0.5*matmul(M/time%dt,u_2)
     END IF
     PinvL =  matmul(P,iL)       
     Mu = Mu - matmul(PinvL,B)
     Mt = Mt - matmul(PinvL,C) 
     IF (Mesh%ndir.gt.0) THEN
        M0 = M0 - matmul(PinvL,Cd)
     END IF


     ! Compute elemental mapping U and U0
     CALL solve_linear_system(Mu,Mt,UU)
     CALL solve_linear_system_sing(Mu,M0,U0)

     LL = matmul(iL,( -matmul(B,UU)+C))
     IF (Mesh%ndir.gt.0) THEN
        L0 = matmul(iL,( -matmul(B,U0)+Cd))
     ELSE
        L0 = matmul(iL,( -matmul(B,U0)))
     ENDIF
     ! Flip faces
		   DO ifa = 1,refElPol%Nfaces
		       IF (Mesh%flipFace(iel,ifa)) THEN
		           UU(:,ind(ifa,:)) = UU(:,ind(ifa,perm))
             LL(:,ind(ifa,:)) = LL(:,ind(ifa,perm))
		       END IF
		   END DO

     ! Store the mapping
     elMat%UU(:,:,iel3) = UU
     elMat%U0(:,iel3)   = U0
     elMat%LL(:,:,iel3) = LL
     elMat%L0(:,iel3)   = L0     
  END DO  
#ifdef TOR3D     
  END DO    
#endif


  
!$OMP END DO  
  DEALLOCATE(S,UU,U0,Mu,Mt,M0,Ed,ind_u,u)
  NULLIFY(M,Cv,D,E,H,Hd)
  IF (time%tis.eq.2 .and. time%ik.gt.1) THEN
     DEALLOCATE(u_2)
  END IF  
  DEALLOCATE(LL,L0,PinvL,P) 
  NULLIFY(iL,B,C,Cd,G)   
!$OMP END PARALLEL
 


  CONTAINS

#ifdef TOR3D  
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
#else		  
					!*****************************************
					! Set permutations for flipping faces
					!****************************************     
					 SUBROUTINE set_permutations(n,m,perm)
					 integer*4, intent(IN)  :: n,m
					 integer*4, intent(OUT) :: perm(:)
					 integer*4              :: i
					 integer*4              :: temp(m,n/m),templr(m,n/m)

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
#endif							 

END SUBROUTINE HDG_mapping
