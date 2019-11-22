!*****************************************
! project: MHDG
! file: compute_solution.f90
! date: 24/01/2017
! Solve the global system, compute face
! and element solution
!*****************************************

SUBROUTINE compute_solution
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE solve_pastix
  USE MPI_OMP
#ifdef PARALL
  USE Communications
#endif
  IMPLICIT NONE
  
  integer*4             :: i,j  
  real,allocatable      :: rhspert(:)
  real                  :: pertamp,errsol
#ifdef PARALL
  integer*4             :: ct,indg(refElPol%Nfacenodes*phys%Neq),indl(refElPol%Nfacenodes*phys%Neq)
#endif  
  
  !******************************** 
  ! Save rhs for accuracy check 
  !********************************
  IF (switch%ckeramp) THEN
     errsol = 0.
     pertamp = 1.e-6
     ALLOCATE(rhspert(matK%n))
     CALL RANDOM_NUMBER(rhspert)               ! rhspert is between 0-1 
     rhspert = (1. + pertamp*2.*(0.5-rhspert)) ! now it is between 1-eps and 1+eps 
     rhspert = rhspert*rhs%vals
  END IF
      
  !******************************** 
  ! Compute the face solution
  !********************************
  ! Solver dependent part
  IF (matPASTIX%start) THEN
					call init_mat_PASTIX(matK,matPASTIX)
					call check_mat_PASTIX(matPASTIX)
					call anal_mat_PASTIX(matPASTIX)
					matPASTIX%start = .false.
		ELSE
					call build_mat_PASTIX(matK, matPASTIX)
		END IF
  call LU_mat_pastix(matPASTIX)
  call solve_mat_PASTIX(matPASTIX)


  !******************************** 
  ! Store face solution
  !********************************
#ifdef PARALL
  ct = 0
  DO i=1,Mesh%Nfaces
     IF (Mesh%ghostfaces(i).eq.1) CYCLE
     ct = ct+1
     indg = (i-1) *Neq*Nfp+(/(j, j=1, Neq*Nfp)/)
     indl = (ct-1)*Neq*Nfp+(/(j, j=1, Neq*Nfp)/)
     sol%u_tilde(indg) = matPASTIX%rhs(indl)
     sol%u_tilde(indg) = (1.-numer%dumpnr)*sol%u_tilde(indg) + numer%dumpnr*matPASTIX%rhs(indl)
  END DO
#else  
  DO i = 1, matPASTIX%n
      sol%u_tilde(i) = (1.-numer%dumpnr)*sol%u_tilde(i) + numer%dumpnr*matPASTIX%rhs(i)
  ENDDO
#endif
 
 
  !******************************** 
  ! Check accuracy of the solution
  !********************************
  IF (switch%ckeramp) THEN
     ! use the perturbed rhs
     matPASTIX%rhs = rhspert
     
     ! solve again..
     CALL solve_mat_PASTIX(matPASTIX)
#ifdef PARALL
					ct = 0
					DO i=1,Mesh%Nfaces
						  IF (Mesh%ghostfaces(i).eq.1) CYCLE
						  ct = ct+1
						  indg = (i-1) *Neq*Nfp+(/(j, j=1, Neq*Nfp)/)
						  indl = (ct-1)*Neq*Nfp+(/(j, j=1, Neq*Nfp)/)
						  errsol = max(errsol,maxval(abs(sol%u_tilde(indg)-matPASTIX%rhs(indl))))
					END DO
#else  
					DO i = 1, matPASTIX%n
						   errsol = max(errsol,abs(sol%u_tilde(i)-matPASTIX%rhs(i)))
					ENDDO
#endif     
     errsol = errsol/pertamp
#ifdef PARALL
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,errsol,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif        
     IF (MPIvar%glob_id.eq.0) THEN
        WRITE(6,*) "ERROR AMPLIFICATION: ", errsol
     ENDIF
     DEALLOCATE(rhspert)  
  END IF
  
  
  !**********************************************
  ! Deallocate the Pastix instance of the matrix
  !**********************************************
!  CALL free_mat_PASTIX(matPASTIX)
  CALL free_mat()

  !**********************************************
  ! MPI communications
  !**********************************************  
#ifdef PARALL
  CALL exchangeSol()
!  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif  
  
  CALL retrieveElementSolution()

  CONTAINS
  
  !**********************************
  ! Retrieve the element solution
  !**********************************   
#ifdef TOR3D  
!***********************************************************************
!                            VERSION 3D TOROIDAL
!***********************************************************************

  SUBROUTINE retrieveElementSolution
  integer               :: itor,iel,iel3d,ifa
  integer               :: Ndim,Neq,N2D,Np2D,Np,Nfl,Ntor,Np1Dpol,Np1Dtor,Nf,Nfdir
  integer               :: ind_ue(refElTor%Nnodes3D*phys%neq),ind_uf(refElTor%Nft*phys%neq)
  integer               :: dd,delta(refElPol%Nfaces+2)
  integer               :: ind_dim(refElPol%Nfaces+2),ind_sta(refElPol%Nfaces+2),aux
  integer*4             :: Fe(refElPol%Nfaces)

   Ndim        = 3                                               ! Number of dimensions
   Neq         = phys%neq                                        ! Number of equations
   N2D         = Mesh%Nelems                                     ! Number elements in the 2D mesh
   Np2D        = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
   Np1Dpol     = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
   Np1Dtor     = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
   Np          = refElTor%Nnodes3D                               ! Number of nodes for each 3D element
   Nfl         = refElTor%Nfl                                    ! Number of nodes in the lateral faces
   Ntor        = numer%ntor
   Nf          = refElPol%Nfaces
   Nfdir = Mesh%Ndir
   
   ind_uf = 0
   ind_dim(1)                   =  Np2D*Neq
   ind_dim(2:refElPol%Nfaces+1) =  Nfl*Neq 
   ind_dim(refElPol%Nfaces+2)   = Np2D*Neq
   ind_sta(1) = 1
   aux = 0
   DO i= 2,refElPol%Nfaces+2
      ind_sta(i) = 1+ind_dim(i-1)+aux
      aux = aux+ind_dim(i-1)
   END DO
   
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(itor,iel,iel3d,ind_ue,ind_uf,Fe,dd,delta,ifa)
!$OMP DO    
   DO itor =1,numer%ntor
      DO iel = 1,N2d
         iel3d = (itor-1)*N2d+iel  !3d numbering of the element
         ind_ue = (iel3d-1)*Np*Neq + (/(i,i=1,Np*Neq)/)
         
         ! Index for the face solution
         Fe     = Mesh%F(iel,:)
         dd = 1+(itor-1)*Neq*(N2d*Np2D+(Mesh%Nfaces-Nfdir)*Nfl)
         delta(1) = dd +(iel-1)*Np2D*Neq
         delta(2:refElPol%Nfaces+1) = dd + N2d*Np2D*Neq+(Fe-1)*Nfl*Neq
         delta(refElPol%Nfaces+2) = dd + N2d*Np2D*Neq+(Mesh%Nfaces-Nfdir)*Nfl*Neq+(iel-1)*Np2D*Neq
         if (itor==numer%ntor) then
            delta(refElPol%Nfaces+2) = 1+(iel-1)*Np2D*Neq
         end if
          
         do ifa = 1,Nf+2
             do i=0,ind_dim(ifa)-1 
                ind_uf(i+ind_sta(ifa)) = delta(ifa)+i
             end do
         end do        
        
        ! elemental solutions
        sol%u(ind_ue) = matmul(elMat%UU(:,:,iel3d),sol%u_tilde(ind_uf)) + elMat%U0(:,iel3d)

      END DO
   END DO
!$OMP END DO
!$OMP END PARALLEL   

END SUBROUTINE retrieveElementSolution
#else 
!***********************************************************************
!                            VERSION 2D
!***********************************************************************
SUBROUTINE retrieveElementSolution   
  integer*4             :: i,j,iel,N2d,Np,Nfp,Neq,Nf,Ndim,ierr,Nfdir
  integer*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq)
  integer*4             :: ind_uf(refElPol%Nfacenodes*phys%Neq*refElPol%Nfaces)  
#ifdef TEMPERATURE
  integer*4             :: ind_ug(Mesh%Nnodesperelem*phys%Neq*Mesh%Ndim)
#endif     
  integer*4             :: Fe(refElPol%Nfaces)
  
   Ndim        = 2
   N2D         = Mesh%Nelems
   Np          = Mesh%Nnodesperelem
   neq         = phys%neq
   Nfp         = Mesh%Nnodesperface
   Nf          = refElPol%Nfaces
     
   DO iel = 1,N2d
     ! Index for the face solution
     Fe     = Mesh%F(iel,:)
     ind_uf = reshape(tensorSumInt((/(i, i=1, Neq*Nfp)/),(Fe-1)*Neq*Nfp),(/Neq*Nfp*Nf/))

     ! Index for the element solution
     ind_ue = (iel-1)*neq*Np + (/(i,i=1,Neq*Np)/)
    
     ! elemental solutions
     sol%u(ind_ue) = matmul(elMat%UU(:,:,iel),sol%u_tilde(ind_uf)) + elMat%U0(:,iel)

#ifdef TEMPERATURE
     ! Index for the element gradient
     ind_ug = (iel-1)*neq*Np*Ndim + (/(i,i=1,Neq*Np*Ndim)/)

     ! elemental solutions
     sol%q(ind_ug) = matmul(elMat%LL(:,:,iel),sol%u_tilde(ind_uf)) + elMat%L0(:,iel)
#endif     
   END DO   
END SUBROUTINE retrieveElementSolution   
   
#endif

   
END SUBROUTINE compute_solution



