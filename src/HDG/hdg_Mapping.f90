SUBROUTINE HDG_mapping()
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE MPI_OMP
  USE in_out
  
  IMPLICIT NONE
  integer*4             :: Ndim,N2d,neq,Np,Nfp,Nfe,i,j,ifa,iel,ntorloc,iel3,Nfg,ierr
  integer*4             :: Np2D,Np1Dpol,Np1Dtor,Nfl,itor
#ifdef TOR3D  
  integer*4             :: ind(refElPol%Nfaces,phys%neq*refElTor%Nfl) 
  integer*4             :: perm(refElTor%Nfl*phys%Neq)
#else
  integer*4             :: ind(refElPol%Nfaces,phys%neq*Mesh%Nnodesperface)
  integer*4             :: perm(refElPol%Nfacenodes*phys%Neq)
#endif  
  real*8,allocatable    :: LL(:,:),UU(:,:),L0(:),U0(:)
  real*8,allocatable    :: Auq_iAqq(:,:),Mu(:,:),Ml(:,:),Md(:),aux(:,:),auxAqq(:,:),auxAquUU(:,:)
  real*8,pointer        :: iAqq(:,:),Aqu(:,:),Aql(:,:)
  real*8,pointer        :: Auq(:,:),Auu(:,:),Aul(:,:)
  real*8,pointer        :: Aql_dir(:),Aul_dir(:),f(:)
  
  IF (MPIvar%glob_id.eq.0) THEN
					IF (utils%printint>1) THEN
						WRITE(6,*) '*************************************************'
						WRITE(6,*) '*            MAPPING                            *'
						WRITE(6,*) '*************************************************' 
					END IF	   
		END IF
		Neq = phys%neq
#ifdef TOR3D		
  Ndim        = 3                                               ! Number of dimensions
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
!$OMP PRIVATE(iel,itor,iel3,LL,UU,L0,U0,Auq_iAqq,Auu,Mu,Ml,Md,f,aux,auxAqq,auxAquUU,iAqq,Aqu,Aql,Auq,Aul,Aql_dir,Aul_dir,ifa)      
  ALLOCATE(LL(Neq*Np*Ndim,Neq*Nfg))
  ALLOCATE(L0(Neq*Np*Ndim))
  ALLOCATE(UU(Neq*Np,Neq*Nfg))
  ALLOCATE(U0(Neq*Np))
  ALLOCATE(Auq_iAqq(Neq*Np,Neq*Ndim*Np))
  ALLOCATE(Ml(Neq*Np,Neq*Nfp))
  ALLOCATE(Mu(Neq*Np,Neq*Np))
  ALLOCATE(Md(Neq*Np))
		ALLOCATE(aux(Neq*Ndim*Np,Neq*Np))	   
		ALLOCATE(auxAqq(Neq*Ndim*Np,Neq*Ndim*Np))
		ALLOCATE(auxAquUU(Neq*Ndim*Np,Neq*Nfg))
!*****************
! Loop in elements
!*****************
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
        iAqq =>elmat%iAqq(:,:,iel3)
        Aqu =>elmat%Aqu(:,:,iel3)
        Aql =>elmat%Aql(:,:,iel3)
        Auq =>elmat%Auq(:,:,iel3)        
        Auu =>elmat%Auu(:,:,iel3)
        Aul =>elmat%Aul(:,:,iel3)
        Aql_dir =>elmat%Aql_dir(:,iel3)
        Aul_dir =>elmat%Aul_dir(:,iel3)
        f   =>elmat%S(:,iel3)
						  Auq_iAqq = matmul(Auq,iAqq)
						  Mu = Auu-matmul(Auq_iAqq,Aqu)
						  Ml = Aul-matmul(Auq_iAqq,Aql)
        Md = f-(Aul_dir-matmul(Auq_iAqq,Aql_dir))
 							call solve_linear_system(Mu,-Ml,UU)
	 						call solve_linear_system_sing(Mu,Md,U0)


 						 auxAquUU=matmul(Aqu,UU) ! Avoid stack overflow
 						 auxAquUU = auxAquUU+Aql
		 					LL = -matmul(iAqq,auxAquUU)
			 				L0 = -matmul(iAqq, Aql_dir+(matmul(Aqu,U0)))
        ! Flip the faces
        DO ifa = 1,refElPol%Nfaces
           IF (Mesh%flipface(iel,ifa)) THEN
               LL(:,ind(ifa,:))     = LL(:,ind(ifa,perm(:)))
               UU(:,ind(ifa,:))     = UU(:,ind(ifa,perm(:)))           
           END if
        END DO
        
        ! Store the mapping
        elMat%LL(:,:,iel3)   = LL
        elMat%UU(:,:,iel3)   = UU
        elMat%L0(:,iel3)     = L0
        elMat%U0(:,iel3)     = U0
        NULLIFY(iAqq,Aqu,Aql,Auq,Aul,Aql_dir,Aul_dir,f)
     END DO
#ifdef TOR3D     
  END DO    
#endif

      
!$OMP END DO
  DEALLOCATE(LL,UU,L0,U0,Auq_iAqq,Mu,Ml,Md,aux,auxAqq,auxAquUU)
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
