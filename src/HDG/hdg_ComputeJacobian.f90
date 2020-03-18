!*****************************************
! project: MHDG
! file: hdg_ConvectionMatrices.f90
! date: 20/12/2016
! Generate the matrices that change
! during the NR iterative process
!*****************************************


SUBROUTINE HDG_computeJacobian()
  USE globals
  USE LinearAlgebra
  USE analytical
  USE physics
  USE printUtils
  USE MPI_OMP
  USE Debug
  
  IMPLICIT NONE


!***********************************************************************
! 
!              COMPUTATION OF THE JACOBIANS
! 
!***********************************************************************
  integer*4             :: Ndim,Neq,N2D,Npel,Npfl,Nfre,Ng1d,Ng2d
  integer*4             :: iel,ifa,iface,i,j,els(2),fas(2)
  integer*4             :: sizeu,sizel
  real*8,allocatable    :: ures(:,:),lres(:,:),u0res(:,:,:)
  real*8                :: Xel(Mesh%Nnodesperelem,2),Xfl(refElPol%Nfacenodes,2)
  real*8,allocatable    :: tau_save(:,:)
  real*8,allocatable    :: xy_g_save(:,:)  
  logical               :: isdir
#ifdef TEMPERATURE
  real*8                :: coefi,coefe
#endif  
#ifdef TOR3D
 ! Definitions in 3D
  integer*4             :: itor,iel3,ntorloc,itorg,Np1Dpol,Np1Dtor,Np2D,Npfp,Ng1dpol,Ng1dtor,Ngfp,Ngfl,Nfdir,Ngvo
  integer*4             :: inde(refElTor%Nnodes3D)
  integer*4             :: indft(refElTor%Nfl),indfp(Mesh%Nnodesperelem),indl(Mesh%Nnodesperelem)  
  integer*4             :: ind_loc(refElPol%Nfaces,refElTor%Nfl*phys%Neq),perm(refElTor%Nfl*phys%Neq)
  integer               :: dd
  integer               :: ind_dim(refElPol%Nfaces+2),ind_sta(refElPol%Nfaces+2),aux
  real*8                :: tdiv(numer%ntor+1),tel(refElTor%Nnodes1d),tg(1),htor
  real*8                :: ue(refElTor%Nnodes3D,phys%Neq),u0e(refElTor%Nnodes3D,phys%Neq,time%tis)
  real*8                :: ufp(Mesh%Nnodesperelem,phys%Neq),uft(refElTor%Nfl,phys%Neq)
  real*8                :: uefp(Mesh%Nnodesperelem,phys%Neq,2),ueft(refElTor%Nfl,phys%Neq,2)
  real*8                :: qe(refElTor%Nnodes3D,phys%Neq*3)
  real*8                :: qefp(Mesh%Nnodesperelem,phys%Neq*3,2),qeft(refElTor%Nfl,phys%Neq*3,2)
  real*8,allocatable    :: qres(:,:)
#else    
 ! Definitions in 2D
  logical               :: save_tau
  integer*4             :: inde(Mesh%Nnodesperelem)
  integer*4             :: indf(refElPol%Nfacenodes)
  integer*4             :: ind_loc(refElPol%Nfaces,refElPol%Nfacenodes*phys%Neq),perm(refElPol%Nfacenodes*phys%Neq)
  real*8                :: ue(Mesh%Nnodesperelem,phys%Neq),u0e(Mesh%Nnodesperelem,phys%Neq,time%tis)
  real*8                :: uf(refElPol%Nfacenodes,phys%Neq),uef(refElPol%Nfacenodes,phys%Neq,2)
  integer               :: indtausave(refElPol%Nfaces*refElPol%Ngauss1d)
  real*8                :: tau_save_el(refElPol%Nfaces*refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Nfaces*refElPol%Ngauss1d,2)
  real*8                :: qe(Mesh%Nnodesperelem,phys%Neq*2),qef(refElPol%Nfacenodes,phys%Neq*2,2)
  real*8,allocatable    :: qres(:,:)
#endif

		IF (utils%printint>1) THEN
   WRITE(6,*) '*************************************************'
   WRITE(6,*) '*          COMPUTING JACOBIAN                   *'
   WRITE(6,*) '*************************************************' 
		END IF	   

 ! Reset matrices
 elMat%Auq = 0.
 elMat%Auu = 0.
 elMat%Aul = 0.
 elMat%Alq = 0.
 elMat%Alu = 0.
 elMat%All = 0.
 elMat%S   = 0.
 elMat%fh  = 0.


#ifdef TEMPERATURE
		coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1+phys%epn)
		coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1+phys%epn)
#endif  

#ifdef TOR3D		
  !*************************************************************
  !               3D stuff
  !*************************************************************
  Ndim        = 3                         ! Number of dimensions
  Neq         = phys%Neq                  ! Number of equations
  N2D         = Mesh%Nelems               ! Number of 2D elements
  Np2D        = refElPol%Nnodes2D         ! Number of nodes in the 2D elements
  Np1Dpol     = refElPol%Nnodes1D         ! Number of nodes in the 1D poloidal segments
  Np1Dtor     = refElTor%Nnodes1D         ! Number of nodes in the 1D toroidal segments  
  Npel        = Np2D*refElTor%Nnodes1D    ! Number of nodes of each element
  Npfl        = Np1Dpol*Np1Dtor           ! Number of nodes of each lateral face
  Npfp        = Np2D                      ! Number of nodes of each poloidal face
  Ng2D        = refElPol%Ngauss2d         ! Number of Gauss points in the 2D element
  Ng1Dpol     = refElPol%Ngauss1d         ! Number of Gauss points in the 1D poloidal segments
  Ng1Dtor     = refEltor%Ngauss1d         ! Number of Gauss points in the 1D toroidal segments
  Ngvo        = Ng2D*Ng1dtor              ! Number of Gauss points for volume computations
  Ngfl        = Ng1dtor*Ng1dpol           ! Number of Gauss points for toroidal faces computations 
  Ngfp        = Ng2D                      ! Number of Gauss points for poloidal faces computations     
		Nfre        = refElPol%Nfaces           ! Number of faces in the reference element		
		Nfdir       = Mesh%Ndir
  ! Toroidal discretization
  htor = numer%tmax/numer%ntor
  tdiv = 0.
  DO i=1,numer%ntor
     tdiv(i+1) = i*htor
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

  ! Local indices
  ind_loc = 0
  DO i = 1,Nfre
     DO j = 1,Npfl*Neq
        ind_loc(i,j) = refElPol%Nnodes2D*Neq + Neq*Npfl*(i-1) + j
     END DO
  END DO
  
  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Np1Dpol,Np1Dtor,Neq,perm)    

#else
  !*************************************************************
  !               2D stuff
  !*************************************************************
  Ndim        = 2
  Neq         = phys%Neq
  N2D         = Mesh%Nelems
  Npel        = refElPol%Nnodes2D
  Npfl        = refElPol%Nfacenodes
  Nfre        = refElPol%Nfaces
  Ng1d        = refElPol%Ngauss1d
  Ng2d        = refElPol%Ngauss2d


  ind_loc = 0
  DO i = 1,Nfre
     DO j = 1,Neq*Npfl
        ind_loc(i,j) = Neq*Npfl*(i-1) + j
     END DO
  END DO

  ! Set perm for flipping faces
  perm = 0
  CALL set_permutations(Neq*Npfl,Neq,perm)

  save_tau = switch%saveTau
#endif  
  
  
  

  




  ! reshape
  sizeu = size(sol%u)
  sizel = size(sol%u_tilde)
  allocate(ures(sizeu/Neq,Neq))
  allocate(lres(sizel/Neq,Neq))
  allocate(u0res(sizeu/Neq,Neq,time%tis))
  ures = transpose(reshape(sol%u,[Neq,sizeu/Neq]))
  lres = transpose(reshape(sol%u_tilde,[Neq,sizel/Neq]))
  do i=1,time%tis
     u0res(:,:,i) =  transpose(reshape(sol%u0(:,i),[Neq,sizeu/Neq]))
  end do
  allocate(qres(sizeu/Neq,Neq*Ndim))
  qres = transpose(reshape(sol%q,[Neq*Ndim,sizeu/Neq]))









#ifdef TOR3D
!********************************************
! 
!                 3D routines
! 
!********************************************


!************************************
!   Loop in elements
!************************************  
!$OMP PARALLEL DEFAULT(SHARED) &
#ifdef TEMPERATURE
!$OMP PRIVATE(iel,iel3,itor,itorg,tel,Xel,i,j,inde,qe,ue)
#else
!$OMP PRIVATE(iel,iel3,itor,itorg,tel,Xel,i,j,inde,ue)
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
     
  DO iel = 1,N2D   

     ! Index of 3D element
     iel3 = (itor-1)*N2d+iel
     
     ! Coordinates of the nodes of the element
     Xel = Mesh%X(Mesh%T(iel,:),:)

     ! Indices to extract the elemental and face solution
     inde = (iel3-1)*Npel + (/(i,i=1,Npel)/)
    
     qe = qres(inde,:)
     ue = ures(inde,:)
     u0e = u0res(inde,:,:)

     ! Compute the matrices for the element
     CALL elemental_matrices_volume(iel3,Xel,tel,qe,ue,u0e)
  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL 






!****************************
! Loop in poloidal faces
!****************************  
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(iface,els,fas,itor,itorg,tel,i,j,Xfl,indfp,inde,indl,ufp,qefp,uefp)
!$OMP DO SCHEDULE(STATIC)
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif     
     tg = tdiv(itorg)
  DO iface = 1,N2D

     els(1) = (itor-1)*N2D+iface
     fas(1) = 1
     els(2) = (itor-2)*N2D+iface
     if (itor==1 .and. MPIvar%ntor==1) then
        els(2) = N2D*(ntorloc-1)+iface
     endif 
     fas(2) = refElPol%Nfaces+2
     
     ! Coordinates of the nodes of the face
     Xel = Mesh%X(Mesh%T(iface,:),:)

     ! Face solution
     dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Npfl)
     indfp = dd+(iface-1)*Np2D+(/(i,i=1,Np2d)/)
     ufp = lres(indfp,:)

     ! Elements solution
     uefp = 0.
     qefp = 0.
     DO j=1,2
        if (els(j)<1) CYCLE
        inde = (els(j)-1)*Npel + (/(i,i=1,Npel)/)
        if (j==1) then
           indl = (/(i,i=1,Np2d)/)
        else
           indl = (/(i,i=Npel-Np2d+1,Npel)/)
        endif
        uefp(:,:,j) = ures( inde( indl ),:)
        qefp(:,:,j) = qres( inde( indl ),:)
     END DO
     
     ! Compute the matrices for the element
     CALL elemental_matrices_faces_pol(els,fas,iface,Xel,tg,qefp,uefp,ufp)
  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL











!************************************
!   Loop in toroidal interior faces
!************************************ 
!!$OMP PARALLEL DEFAULT(SHARED) & 
!!$OMP PRIVATE(iface,itor,itorg,tel,els,fas,inde,indft,Xfl,i,j,dd,qeft,ueft,uft) 
!!$OMP DO SCHEDULE(STATIC) 
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif     
     tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg+1)-tdiv(itorg))
   
  DO iface = 1,Mesh%Nintfaces

     els = Mesh%intfaces(iface,(/1,3/))
     fas = Mesh%intfaces(iface,(/2,4/))
            
     ! Coordinates of the nodes of the face
     Xfl = Mesh%X(Mesh%T(els(1),refElPol%face_nodes(fas(1),:)),:)
     
     ! Face solution
     dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Npfl)+N2D*Np2D
     indft = dd+(iface-1)*Npfl+(/(i,i=1,Npfl)/)
     uft = lres(indft,:)
     
     ! Elements solution
     DO j=1,2
        inde = ((itor-1)*N2d+els(j)-1)*Npel + (/(i,i=1,Npel)/)
        ueft(:,:,j) = ures( inde( refElTor%faceNodes3(fas(j),:) ),:)
        qeft(:,:,j) = qres( inde( refElTor%faceNodes3(fas(j),:) ),:)
     END DO
     
     ! Compute the matrices for the element
     CALL elemental_matrices_faces_int((itor-1)*N2d+els,fas+1,els,Xfl,tel,qeft,ueft,uft)
         
     ! Flip faces 
     DO j=1,2
        IF (Mesh%flipface(els(j),fas(j)) ) then
           iel = (itor-1)*N2d+els(j)
           elMat%Alq(ind_loc(fas(j),:),:,iel)   = elMat%Alq(ind_loc(fas(j),perm),:,iel)
           elMat%Alu(ind_loc(fas(j),:),:,iel)   = elMat%Alu(ind_loc(fas(j),perm),:,iel)
           elMat%All(ind_loc(fas(j),:),:,iel)   = elMat%All(ind_loc(fas(j),perm),:,iel)
           elMat%All(:,ind_loc(fas(j),:),iel)   = elMat%All(:,ind_loc(fas(j),perm),iel)
           elMat%fh(ind_loc(fas(j),:),iel)      = elMat%fh(ind_loc(fas(j),perm),iel)
        END if
     END DO
  END DO
  END DO
!!$OMP END DO
!!$OMP END PARALLEL 







!*********************************
! Loop in toroidal exterior faces
!*********************************   
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(iface,iel,iel3,ifa,itor,itorg,tel,Xfl,indft,dd,inde,uft,ueft,qeft)
!$OMP DO SCHEDULE(STATIC)
  DO itor = 1,ntorloc
#ifdef PARALL  
     itorg = itor+(MPIvar%itor-1)*numer%ntor/MPIvar%ntor
     if (itorg==numer%ntor+1) itorg=1
#else
     itorg = itor
#endif     
     tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg+1)-tdiv(itorg))
   
  DO iface = 1,Mesh%Nextfaces

     iel = Mesh%extfaces(iface,1)
     ifa = Mesh%extfaces(iface,2)
     
     ! Index of 3D element
     iel3 = (itor-1)*N2d+iel  
          
     ! Coordinates of the nodes of the element
     Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)

     ! Face solution
     isdir = Mesh%Fdir(iel,ifa)
     IF (isdir) THEN
        uft = 0.
     ELSE
        dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Nfdir)*Npfl)+N2D*Np2D
        indft = dd+(iface+Mesh%Nintfaces-1)*Npfl+(/(i,i=1,Npfl)/)
        uft = lres(indft,:)
     ENDIF
     
     ! Elements solution
     inde = (iel3-1)*Npel + (/(i,i=1,Npel)/)
     ueft(:,:,1) = ures( inde( refElTor%faceNodes3(ifa,:) ),:)
     qeft(:,:,1) = qres( inde( refElTor%faceNodes3(ifa,:) ),:)
      
     ! Compute the matrices for the element
     CALL elemental_matrices_faces_ext(iel3,ifa+1,iel,Xfl,tel,qeft(:,:,1),ueft(:,:,1),uft,isdir)
#ifdef PARALL
     ! Flip faces (only for parallel)
     IF (Mesh%flipface(iel,ifa) ) then
        elMat%Alq(ind_loc(ifa,:),:,iel3)   = elMat%Alq(ind_loc(ifa,perm),:,iel3)
        elMat%Alu(ind_loc(ifa,:),:,iel3)   = elMat%Alu(ind_loc(ifa,perm),:,iel3)
        elMat%All(ind_loc(ifa,:),:,iel3)   = elMat%All(ind_loc(ifa,perm),:,iel3)
        elMat%All(:,ind_loc(ifa,:),iel3)   = elMat%All(:,ind_loc(ifa,perm),iel3)
        elMat%fh(ind_loc(ifa,:),iel3)      = elMat%fh(ind_loc(ifa,perm),iel3)
     END if
#endif

  END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  deallocate(ures,lres,u0res)
  deallocate(qres)

 CONTAINS




			
!***************************************************
! Volume computation in 3D
!***************************************************
					SUBROUTINE elemental_matrices_volume(iel,Xel,tel,qe,ue,u0e)
       integer,intent(IN)			     :: iel 
					  real*8,intent(IN)         :: Xel(:,:),tel(:)
					  real*8,intent(IN)         :: qe(:,:)
       real*8,intent(IN)         :: ue(:,:)
       real*8,intent(IN)         :: u0e(:,:,:)
							integer*4                 :: g,NGaussPol,NGaussTor,igtor,igpol,i,j,k,iord
							real*8                    :: dvolu,dvolu1d,htor
       real*8                    :: J11(Ng2d),J12(Ng2d)
       real*8                    :: J21(Ng2d),J22(Ng2d)
       real*8                    :: detJ(Ng2d)
       real*8                    :: iJ11(Ng2d),iJ12(Ng2d)
       real*8                    :: iJ21(Ng2d),iJ22(Ng2d)
       real*8                    :: fluxg(Ng2d)
							real*8                    :: xy(Ng2d,2),teg(Ng1dtor)
       real*8                    :: ueg(Ngvo,neq),upg(Ngvo,phys%npv),ueg(Ngvo,neq,time%tis)
       real*8                    :: qeg(Ngvo,neq*Ndim)
       real*8                      :: force(Ngvo,Neq)
							integer*4,dimension(Npel)   :: ind_ass,ind_asg,ind_asq
							real*8                      :: ktis(time%tis+1)
       real*8                      :: Nxg(Np2D),Nyg(Np2D),Nx_ax(Np2D)
       real*8,pointer              :: N1g(:),N2g(:),N3g(:)
       real*8                      :: N1xg(Np1Dtor),N1xg_cart(Np1Dtor)
       real*8,dimension(Npel)      :: Nr,Nz,Nt,Ni,NNbb
       real*8,dimension(Npel,Npel) :: NNi,NxNi,NyNi,NtNi,NNxy
       real*8                      :: NxNy(Npel,Npel,Ndim)

       real*8                    :: kmult(Npfl,Npfl)
							real*8, parameter         :: tol = 1e-12					
       integer*4                 :: ind2(Ng2d)		
       real*8                    :: b(Ngvo,3),divb(Ngvo),drift(Ngvo,2),divbg,driftg(2)
       real*8                    :: b2d(Ng2d,3),divb2d(Ng2d),drift2d(Ng2d,2)
       real*8                    :: diff_iso_vol(Neq,Neq,Ngvo),diff_ani_vol(Neq,Neq,Ngvo)
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
       real*8                    :: Bmod_g(Ng2d),Bmod_x,Bmod_y,Nx_ax(Npel)
       real*8                    :: Bt(Mesh%Nnodesperelem),Bt_g(Ng2d)              
#endif
	       
       ind_ass = (/ (i, i = 0, Neq*(Npel-1),Neq) /)
       ind_asg = (/ (i, i = 0, Neq*(Npel-1)*Ndim,Neq*Ndim) /)
       ind_asq = (/ (i, i = 0, Neq*(Npel-1)*Ndim,Neq*Ndim) /)

       !***********************************
       !    Volume computation
       !***********************************
							
							! Gauss points position
							xy = matmul(refElPol%N2D,Xel)

							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)
							
						 ! toroidal element size
						 htor = tel(Np1Dtor) - tel(1)
						 
							! Compute diffusion at Gauss points
							CALL setLocalDiff(xy,diff_iso_vol,diff_ani_vol) ! TODO: change this for 3D

							! Solution at Gauss points
							ueg = matmul(refElTor%N3D,ue)
							qeg = matmul(refElTor%N3D,qe) 
							
       ! Solution at previous time steps, at Gauss points
       do i=1,time%tis
          u0eg(:,:,i) = matmul(refElTor%N3D,u0e(:,:,i))
       end do

			    ! Physical variables at Gauss points
			    CALL cons2phys(ueg,upg)
       
       ! Constant sources
							! Body force at the integration points
							CALL body_force(xy(:,1),xy(:,2),teg,force)			

							!! Some sources for West cases
							IF (switch%testcase	.ge.51 .and. switch%testcase	.le.54) THEN						
							   fluxg = matmul(refElPol%N2D,phys%flux2D(Mesh%T(iel,:)))
							   DO g=1,Ngvo
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
       b = matmul(refElPol%N2D,b_nod)              ! Ngauss
#ifdef TEMPERATURE        
       IF (switch%driftvel) THEN
        
          ! Drift at Gauss points       
          Bmod_g = matmul(refElPol%N2D,Bmod) ! Ngauss
          Bt_g   = matmul(refElPol%N2D,phys%Bt(Mesh%T(iel,:)))  ! Ngauss 
       END IF       
#endif          
#endif		
!************ END SELECTION MOVING EQUILIBRIUM **********************

							! Loop in 2D Gauss points
       J11 = matmul(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
       J12 = matmul(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
       J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
       J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
				   iJ11 =  J22/detJ
				   iJ12 = -J12/detJ
				   iJ21 = -J21/detJ
				   iJ22 =  J11/detJ
				                 
       ! Coefficient time integration scheme
       call setTimeIntegrationCoefficients(ktis)

       NgaussPol = refElPol%NGauss2D
       NgaussTor = refElTor%NGauss1D
	 					DO igtor = 1,NGaussTor
		 				   N1g    => refElTor%N1D(igtor,:)         ! Toroidal shape function 
			 			   N1xg_cart   = refElTor%Nxi1D(igtor,:)*2/htor       ! Toroidal shape function derivative
				 		   dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume
				 		   
					 	   DO igpol = 1,NGaussPol
						       g = (igtor-1)*NGaussPol+igpol   



								 	   ! Poloidal shape functions and derivatives    							  
   						 	  N2g => refElPol%N2D(igpol,:)										   
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
						       N3g => refElTor%N3D(g,:)               ! 3D shape function
						       Nr  = col(TensorProduct(Nxg,N1g))      ! 3D shape function, derivative in r for divergence
						       Nz  = col(TensorProduct(Nyg,N1g))      ! 3D shape function, derivative in z
						       Nt  = col(TensorProduct(N2g,N1xg))     ! 3D shape function, derivative in t    										
						       
						       
						       ! Shape functions products
						       Ni   = N3g*dvolu                                                ! Npel x 1
						       NNi  = tensorProduct(N3g,Ni)                                    ! Npel x Npel 
						       NxNi = tensorProduct(Nr,Ni)                                     ! Npel x Npel 
             NyNi = tensorProduct(Nz,Ni)                                     ! Npel x Npel 
             NtNi = tensorProduct(Nt,Ni)                                     ! Npel x Npel
             NNxy = b(g,1)*NxNi+b(g,2)*NyNi+b(g,3)*NtNi                      ! Npel x Npel
             NxNy(:,:,1)    = NxNi 
             NxNy(:,:,2)    = NyNi                                            
             NxNy(:,:,3)    = NtNi                                           ! Npel x Npel x 3
             NNbb =  (Nr*b(g,1)+Nz*b(g,2)+Nt*b(g,3))*dvolu                   ! Npel x 1


            ! Divergence of b at the Gauss points  
#ifndef MOVINGEQUILIBRIUM
            divbg = divb(g)
						      driftg(1) = drift(g,1)
						      driftg(2) = drift(g,2)
#else      
						      IF (switch%axisym) THEN
						         Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
						      ELSE
						         Nx_ax = Nxg
						      END IF  
            divbg = dot_product(Nx_ax,b_nod(:,1))+dot_product(Nyg,b_nod(:,2))    
								    Bmod_x  = dot_product(Nxg,Bmod) ! 1
								    Bmod_y  = dot_product(Nyg,Bmod) ! 1
								    driftg(1) =  phys%dfcoef*Bt_g(g)*Bmod_y/Bmod_g(g)**3 ! I multiply for the temperature later
								    driftg(2) = -phys%dfcoef*Bt_g(g)*Bmod_x/Bmod_g(g)**3 ! I multiply for the temperature later             
#endif 

          CALL assemblyVolumeContribution(iel,ind_ass,ind_asg,ind_asq,b(g,:),divbg,driftg,force(g,:),&
          &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxg,Nyg,NNxy,NxNy,NNbb,upg(g,:),ueg(g,:),qeg(g,:),u0eg(g,:,:))
							   END DO ! END loop in volume Gauss points
							END DO
							
					END SUBROUTINE elemental_matrices_volume






!*****************************************
! Poloidal faces computations in 3D
!*****************************************
     SUBROUTINE elemental_matrices_faces_pol(els,fas,iel2,Xfp,tg,qef,uef,uf)
							integer*4,intent(IN)  :: els(:),fas(:),iel2
       real*8,intent(IN)     :: Xfp(:,:),tg(:)
       real*8,intent(IN)     :: qef(:,:,:)	
       real*8,intent(IN)     :: uef(:,:,:),uf(:,:)
       real*8                :: ufg(Ng2d,Neq),uefg(Ng2d,Neq),upgf(Ng2d,phys%npv)       
       real*8                :: qfg(Ng2d,Neq*Ndim)
       integer*4             :: ind_asf(Np2D),ind_ash(Np2D)
							integer*4             :: g,NGauss,i,j,lel
							real*8                :: dsurf(Ng2d),bn
							real*8                :: xyf(Ng2d,2)
       real*8                :: J11(Ng2d),J12(Ng2d)
       real*8                :: J21(Ng2d),J22(Ng2d)
       real*8                :: detJ(Ng2d)
       real*8                :: iJ11(Ng2d),iJ12(Ng2d)
       real*8                :: iJ21(Ng2d),iJ22(Ng2d)						
       real*8,pointer        :: Nfg(:)
       real*8                :: NNif(Np2D,Np2D),Nfbn(Np2D)
       real*8                :: n_g(Ng2d,3)
							integer*4             :: ind_ff(Neq*Np2D),ind_fe(Neq*Np2D),ind_fg(Neq*Ndim*Np2D)       
       real*8, parameter     :: tol = 1e-12							       
       real*8                :: b_f(Ng2d,3),tau(Neq,Neq)			
       real*8                :: diff_iso_fac(Neq,Neq,Ng2D),diff_ani_fac(Neq,Neq,Ng2D)

       
       
       ind_asf = (/ (i, i = 0, Neq*(Np2D-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Np2D-1)*Ndim,Neq*Ndim) /)
       							       
							! Gauss points position
							xyf = matmul(refElPol%N2D,Xfp)
							
       ! Compute diffusion at faces Gauss points
					  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

							! Loop in 2D Gauss points
       Ngauss = Ng2d
       
       ! Trace solution at face Gauss points
       ufg = matmul(refElPol%N2D,uf)       


       ! Physical variables related to the trace solution
       CALL cons2phys(ufg,upgf) 
       
       ! Magnetic field in the poloidal face
#ifndef MOVINGEQUILIBRIUM			
!************ CASE OF STATIC EQUILIBRIUM ****************************
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
						       b_f = matmul(refElPol%N2D,phys%b(Mesh%T(iel,:),:))
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),tg,b_f)   
										END IF				
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
          b_f = matmul(refElPol%N1D,b_nod(refElPol%face_nodes(ifa,:),:))
#endif
!************ END SELECTION MOVING EQUILIBRIUM **********************


       J11 = matmul(refElPol%Nxi2D,Xfp(:,1))                           ! ng x 1
       J12 = matmul(refElPol%Nxi2D,Xfp(:,2))                           ! ng x 1
       J21 = matmul(refElPol%Neta2D,Xfp(:,1))                          ! ng x 1
       J22 = matmul(refElPol%Neta2D,Xfp(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
				   iJ11 =  J22/detJ
				   iJ12 = -J12/detJ
				   iJ21 = -J21/detJ
				   iJ22 =  J11/detJ 
				   dsurf  =  refElPol%gauss_weights2D*detJ      
       DO  lel = 1,2
          iel = els(lel)
          if (iel<1) CYCLE
          ifa = fas(lel)
          IF (ifa==1) THEN
	            ind_ff = (/(i,i=1,Np2D*Neq)/)
             ind_fe = (/(i,i=1,Np2D*Neq)/)
	            ind_fg = (/(i,i=1,Np2D*Ndim*Neq)/)	
	            ! Exterior normal
							      n_g = 0.; n_g(:,3) = -1			   
	         ELSE
	            ind_ff = Np2D*Neq + refElPol%Nfaces*Npfl*Neq+(/(i,i=1,Np2D*Neq)/)
	            ind_fe = Np2d*(Np1dTor-1)*Neq + (/(i,i=1,Np2D*Neq)/) 	
	            ind_fg = Np2d*(Np1dTor-1)*Ndim*Neq + (/(i,i=1,Np2D*Ndim*Neq)/) 	     
	            ! Exterior normal
							      n_g = 0.; n_g(:,3) = 1	          
	         ENDIF
	      
             ! Element solution at face Gauss points										
             uefg = matmul(refElPol%N2D,uef(:,:,lel))
             ! Gradient solution at face gauss points
             qfg    = matmul(refElPol%N2D,qef(:,:,lel))


							   DO g = 1,NGauss

							      ! Shape functions product
             Nfg => refElPol%N2D(g,:)
							      bn = dot_product(b_f(g,:),n_g(g,:))
						       NNif = tensorProduct(Nfg,Nfg)*dsurf(g)
						       Nfbn = bn*Nfg*dsurf(g)  
						    
						       ! Compute the stabilization term			
						       tau = 0.			   			
						       IF (numer%stab==1) THEN
						     		! Constant stabilization
						          DO i=1,Neq
						             tau(i,i) = numer%tau(i)
						          END DO
						       ELSE
						       ! Non constant stabilization
						          ! Compute tau in the Gauss points
						          IF (numer%stab<5) THEN
						             CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),iel2,ifa,0.,xyf(g,:),tau)
						          ELSE
						             CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),0.,iel2,tau)
						          ENDIF
						       END IF
				       
						       ! Assembly local contributions              				        										
						       CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
						       n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau,ifa) 
							   END DO
							END DO

      
					END SUBROUTINE elemental_matrices_faces_pol
					
					
					



!***************************************************
! Interior faces computation in 3D
!***************************************************
					SUBROUTINE elemental_matrices_faces_int(els,fas,els2,Xfl,tel,qef,uef,uf)
       integer,intent(IN)			     :: els(:),fas(:),els2(:)
					  real*8,intent(IN)         :: Xfl(:,:),tel(:)
					  real*8,intent(IN)         :: qef(:,:,:)
       real*8,intent(IN)         :: uef(:,:,:),uf(:,:)	

							integer*4                 :: g,ifa,i,j,k,lel,iel,ieln,igtor,igpol
							real*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
							real*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)						
							real*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
							real*8                    :: dsurf(Ngfl),dsurfg
       real*8                    :: qfg(Ngfl,neq*Ndim)
							integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
       integer*4,dimension(Npfl) :: ind_asf,ind_ash
       integer*4,dimension(Npfl) :: indf,ind_if,ind_jf,ind_kf
       integer*4                 :: ind(Ng1Dpol)
       integer*4                 :: permsing(Npfl)
							real*8                    :: t_g(Ng1dpol,2),n_g(Ngfl,3), bn
							real*8                    :: NNif(Npfl,Npfl),Nfbn(Npfl)
							real*8,pointer            :: Nfg(:)
       real*8                    :: tau(Neq,Neq)
							real*8, parameter         :: tol = 1e-12							
       real*8                    :: b_f(Ngfl,3),b_f2d(Ng1Dpol,3)
       real*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

      
       ind_asf = (/ (i, i = 0, Neq*(Npfl-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Npfl-1)*Ndim,Neq*Ndim) /)

						 ! toroidal element size
						 htor = tel(Np1Dtor) - tel(1)
						 
							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)				 
       
       ! Loop in the faces of the element
							DO  lel = 1,2
										
          iel = els(lel)
          ifa = fas(lel)
          if (lel==1) then
             ieln = els(2)
          else
             ieln = els(1)
          endif

										! Indices
          ind_ff = Np2d*Neq + (ifa-2)*Npfl*Neq+ (/(i,i=1,Npfl*Neq)/)
          ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*Neq)) 
          ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*3*Neq))
                 
          ! Coordinates, derivatives and trace solution at face Gauss points
										IF (Mesh%flipFace(iel,ifa-1)) THEN
										    call set_permutations(Np1Dpol,Np1Dtor,1,permsing)
										    ufg = matmul(refElTor%sFTF,uf(permsing,:))
              xyf = matmul(refElPol%N1D,Xfl((/(i,i=Np1Dpol,1,-1)/),:))
														! Shape function derivatives at Gauss points
              xyDer  = matmul(refElPol%Nxi1D,Xfl((/(i,i=Np1Dpol,1,-1)/),:))
														
          ELSE
              ufg = matmul(refElTor%sFTF,uf)
              xyf = matmul(refElPol%N1D,Xfl)
														! Shape function derivatives at Gauss points
														xyDer  = matmul(refElPol%Nxi1D,Xfl)
										END IF 
										
          ! Element solution at face Gauss points										
          uefg = matmul(refElTor%sFTF,uef(:,:,lel))
          ! Gradient solution at face gauss points
          qfg    = matmul(refElTor%sFTF,qef(:,:,lel))
          ! Compute diffusion at faces Gauss points
							   CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

          ! Physical variables at face Gauss points
          CALL cons2phys(ufg,upgf)

										! Compute dsurf
          xydNorm_g = sqrt(xyDer(:,1)**2+xyDer(:,2)**2)										
          dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))
          
          ! Compute exterior normal
          t_g(:,1) = xyDer(:,1)/xydNorm_g
          t_g(:,2) = xyDer(:,2)/xydNorm_g          
          n_g = 0.
          DO i=1,Ng1dTor
             ind = (i-1)*Ng1dPol + (/(j,j=1,Ng1dPol)/) 
             n_g(ind,1) =  t_g(:,2)
             n_g(ind,2) = -t_g(:,1)
          END DO


#ifndef MOVINGEQUILIBRIUM			
!************ CASE OF STATIC EQUILIBRIUM ****************************
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
						          b_f2d = matmul(refElPol%N1D,phys%b(Mesh%T(iel,refElPol%face_nodes(ifa-1,:)),:))
                DO igtor = 1,refElTor%NGauss1D
                   ind = (igtor-1)*refElPol%NGauss1D + (/(j,j=1,refElPol%NGauss1D)/)
                   b_f(ind,:)     = b_f2d
                END DO				
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),teg,b_f)   
										END IF				
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
          b_f = matmul(refElPol%N1D,b_nod(refElPol%face_nodes(ifa,:),:))
#endif
!************ END SELECTION MOVING EQUILIBRIUM **********************
										
										
										
          !*****************************
          ! Loop in face Gauss points
          !*****************************
										DO igtor = 1,Ng1dTor
										   DO igpol = 1,Ng1dPol
										
										      g = (igtor-1)*Ng1dPol+igpol
										
										      ! Face shape functions										      
										      Nfg => refElTor%sFTF(g,:)
										
												    IF (switch%axisym) THEN			
												       dsurfg = dsurf(g)*xyf(igpol,1)
												    ELSE
												       dsurfg = dsurf(g)
												    END IF									

										      ! Shape functions product
										      bn = dot_product(b_f(g,:),n_g(g,:))
                NNif = tensorProduct(Nfg,Nfg)*dsurfg
                Nfbn = bn*Nfg*dsurfg   
 				        
    				        ! Compute the stabilization term			
    				        tau = 0.			   			
                IF (numer%stab==1) THEN
              		! Constant stabilization
						             DO i=1,Neq
						                tau(i,i) = numer%tau(i)
						             END DO
                ELSE
                ! Non constant stabilization
                   ! Compute tau in the Gauss points
                   IF (numer%stab<5) THEN
                      CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),els2(lel),ifa,0.,xyf(g,:),tau)
                   ELSE
                      CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),0.,els2(lel),tau)
                   ENDIF
                END IF
                
                ! Assembly local contributions              				        										  
                CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
                n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau,ifa) 
             
										   END DO ! Gauss points
										END DO
							END DO ! 2 elements

					END SUBROUTINE elemental_matrices_faces_int	
					
					
					
					
	
	
	
	



			
!***************************************************
! Exterior faces computation in 3D
!***************************************************
					SUBROUTINE elemental_matrices_faces_ext(iel,ifa,iel2,Xfl,tel,qef,uef,uf,isdir)
       integer,intent(IN)			     :: iel ,ifa,iel2
					  real*8,intent(IN)         :: Xfl(:,:),tel(:)
					  logical,intent(IN)        :: isdir
					  real*8,intent(IN)         :: qef(:,:)
       real*8,intent(INOUT)      :: uef(:,:),uf(:,:)
	
							integer*4                 :: i,j,k,g,igtor,igpol
							real*8                    :: xyf(Ng1Dpol,2),teg(Ng1dtor)
							real*8                    :: xyDer(Ng1Dpol,2),xydNorm_g(Ng1Dpol)								
							real*8                    :: ufg(Ngfl,neq),uefg(Ngfl,neq),upgf(Ngfl,phys%npv)
							real*8                    :: dsurf(Ngfl),dsurfg
       real*8                    :: qfg(Ngfl,neq*Ndim)
							integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
       integer*4,dimension(Npfl) :: ind_asf,ind_ash
       integer*4,dimension(Npfl) :: indf,ind_if,ind_jf,ind_kf
       integer*4                 :: ind(Ng1Dpol)
       integer*4                 :: permsing(Npfl)
       real                      :: isext
							real*8                    :: t_g(Ng1Dpol,2),n_g(Ngfl,3), bn
							real*8                    :: NNif(Npfl,Npfl),Nfbn(Npfl)
							real*8,pointer            :: Nfg(:)
							real*8                    :: tau(Neq,Neq) 
													
							real*8, parameter         :: tol = 1e-12							
       real*8                    :: b_f(Ngfl,3),b_f2d(Ng1Dpol,3)			
       real*8                    :: diff_iso_fac(Neq,Neq,Ngfl),diff_ani_fac(Neq,Neq,Ngfl)

       
       ind_asf = (/ (i, i = 0, Neq*(Npfl-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Npfl-1)*Ndim,Neq*Ndim) /)

						 ! toroidal element size
						 htor = tel(Np1Dtor) - tel(1)
						 
							! Gauss points position in the toroidal direction
							teg = matmul(refElTor%N1d,tel)	
									
							! Indices
       ind_ff = Np2d*Neq + (ifa-2)*Npfl*Neq+ (/(i,i=1,Npfl*Neq)/)
       ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*Neq)) 
       ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifa-1,:)-1)*3*Neq))
              
										
       ! Trace solution at face Gauss points
       xyf = matmul(refElPol%N1D,Xfl)
       xyDer  = matmul(refElPol%Nxi1D,Xfl)
       IF (isdir) THEN
          CALL analytical_solution(xyf(:,1),xyf(:,2),tel,ufg) 
       ELSE
#ifdef PARALL
						    IF  (Mesh%flipface(iel2,ifa-1)) THEN
						       call set_permutations(Np1Dpol,Np1Dtor,1,permsing)
						       uf = uf(permsing,:)
						    END IF
#endif 
							   ufg = matmul(refElTor%sFTF,uf)
       END IF
        
       ! Element solution at face Gauss points										
       uefg = matmul(refElTor%sFTF,uef)
       ! Gradient solution at face gauss points
       qfg    = matmul(refElTor%sFTF,qef)
          
      ! Compute diffusion at faces Gauss points
				  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

      ! Physical variables at face Gauss points
      CALL cons2phys(ufg,upgf)

#ifndef MOVINGEQUILIBRIUM			
!************ CASE OF STATIC EQUILIBRIUM ****************************
						    IF (switch%testcase	.ge.50 .and. switch%testcase	.le.59) THEN
			          b_f2d = matmul(refElPol%N1D,phys%b(Mesh%T(iel,refElPol%face_nodes(ifa-1,:)),:))
             DO igtor = 1,refElTor%NGauss1D
                ind = (igtor-1)*refElPol%NGauss1D + (/(j,j=1,refElPol%NGauss1D)/)
                b_f(ind,:)     = b_f2d
             END DO				
						    ELSE
													CALL defineMagneticField(xyf(:,1),xyf(:,2),teg,b_f)   
										END IF				
#else
!************ CASE OF MOVING EQUILIBRIUM ****************************
          b_f = matmul(refElPol%N1D,b_nod(refElPol%face_nodes(ifa,:),:))
#endif
!************ END SELECTION MOVING EQUILIBRIUM **********************


							! Compute dsurf
       xydNorm_g = sqrt(xyDer(:,1)**2+xyDer(:,2)**2)										
       dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))
       
       ! Compute exterior normal
       t_g(:,1) = xyDer(:,1)/xydNorm_g
       t_g(:,2) = xyDer(:,2)/xydNorm_g     
       n_g = 0.     
       DO i=1,Ng1dTor
          ind = (i-1)*Ng1dPol + (/(j,j=1,Ng1dPol)/) 
          n_g(ind,1) =  t_g(:,2)
          n_g(ind,2) = -t_g(:,1)
       END DO
				  
				  
				  
      !*****************************
      ! Loop in face Gauss points
      !*****************************
					  DO igtor = 1,Ng1dTor
					     DO igpol = 1,Ng1dPol
					  
					        g = (igtor-1)*Ng1dPol+igpol
					  
					        ! Face shape functions										      
					        Nfg => refElTor%sFTF(g,:)
					  
							      IF (switch%axisym) THEN			
							         dsurfg = dsurf(g)*xyf(igpol,1)
							      ELSE
							         dsurfg = dsurf(g)
							      END IF									
					        
					        ! Shape functions product
					        bn = dot_product(b_f(g,:),n_g(g,:))
             NNif = tensorProduct(Nfg,Nfg)*dsurfg
             Nfbn = bn*Nfg*dsurfg   
          
			          ! Compute the stabilization term		
			          isext = 1.
#ifdef PARALL    
             IF (Mesh%boundaryFlag(Mesh%F(iel,ifa-1)-Mesh%Nintfaces).eq.0) THEN
                isext = 0.
             END IF
#endif			          	
			          tau = 0.			   			
             IF (numer%stab==1) THEN
           		! Constant stabilization
	               DO i=1,Neq
	                  tau(i,i) = numer%tau(i)
	               END DO
             ELSE
             ! Non constant stabilization
                ! Compute tau in the Gauss points
                IF (numer%stab<5) THEN
                   CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),iel2,ifa,isext,xyf(g,:),tau)
                ELSE
                   CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g(g,:),xyf(g,:),isext,iel2,tau)
                ENDIF
             END IF

            ! Assembly local contributions
#ifdef PARALL
             IF (Mesh%boundaryFlag(Mesh%F(iel2,ifa-1)-Mesh%Nintfaces).eq.0) THEN
						          ! Ghost face: assembly it as interior
						          CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
						          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau,ifa) 
             ELSE
						          CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
						          n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau,ifa) 
             ENDIF
#else
             CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
             n_g(g,:),diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau,ifa) 
#endif
		        END DO ! Gauss points
		     END DO 
										

					END SUBROUTINE elemental_matrices_faces_ext		
										
					
					
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




!********************************************
! 
!                 2D routines
! 
!********************************************


!************************************
!   Loop in elements in 2D
!************************************  
!$OMP PARALLEL DEFAULT(SHARED) & 
!$OMP PRIVATE(iel,inde,Xel,i,qe,ue) 
!$OMP DO SCHEDULE(STATIC) 
  DO iel = 1,N2D   

     ! Coordinates of the nodes of the element
     Xel = Mesh%X(Mesh%T(iel,:),:)

     ! Indices to extract the elemental and face solution
     inde = (iel-1)*Npel + (/(i,i=1,Npel)/)
    
     qe = qres(inde,:)
     ue = ures(inde,:)
     u0e = u0res(inde,:,:)

     ! Compute the matrices for the element
     CALL elemental_matrices_volume(iel,Xel,qe,ue,u0e)
  END DO
!$OMP END DO
!$OMP END PARALLEL 


 uef = 0.
 qef = 0.
!************************************
!   Loop in interior faces in 2D
!************************************ 
!!$OMP PARALLEL DEFAULT(SHARED) & 
!#ifdef TEMPERATURE
!!$OMP PRIVATE(iface,els,fas,inde,indf,Xfl,i,j,qef,uef,uf,tau_save_el,xy_g_save_el,indtausave) 
!#else
!!$OMP PRIVATE(iface,els,fas,inde,indf,Xfl,i,j,uef,uf,tau_save_el,xy_g_save_el,indtausave) 
!#endif
!!$OMP DO SCHEDULE(STATIC) 
  DO iface = 1,Mesh%Nintfaces

     els = Mesh%intfaces(iface,(/1,3/))    
     fas = Mesh%intfaces(iface,(/2,4/))
            
     ! Coordinates of the nodes of the face
     Xfl = Mesh%X(Mesh%T(els(1),refElPol%face_nodes(fas(1),:)),:)
     
     ! Face solution
     indf = (iface-1)*Npfl + (/(i,i=1,Npfl)/)
     uf = lres(indf,:)
     
     ! Elements solution
     DO j=1,2
        inde = (els(j)-1)*Npel + (/(i,i=1,Npel)/)
        uef(:,:,j) = ures( inde( refElPol%face_nodes(fas(j),:) ),:)
        qef(:,:,j) = qres( inde( refElPol%face_nodes(fas(j),:) ),:)
     END DO
     
     ! Compute the matrices for the element
     CALL elemental_matrices_faces_int(els,fas,Xfl,qef,uef,uf,tau_save_el,xy_g_save_el)

     if (save_tau) then
        indtausave = (iel-1)*Nfre*Ng1d + (/(i,i =1,Nfre*Ng1d)/)
        tau_save(indtausave,:)= tau_save_el
        xy_g_save(indtausave,:)= xy_g_save_el
     endif
     
     ! Flip faces 
     DO j=1,2
        IF (Mesh%flipface(els(j),fas(j)) ) then
           elMat%Alq(ind_loc(fas(j),:),:,els(j))   = elMat%Alq(ind_loc(fas(j),perm),:,els(j))
           elMat%Alu(ind_loc(fas(j),:),:,els(j))   = elMat%Alu(ind_loc(fas(j),perm),:,els(j))
           elMat%All(ind_loc(fas(j),:),:,els(j))   = elMat%All(ind_loc(fas(j),perm),:,els(j))
           elMat%All(:,ind_loc(fas(j),:),els(j))   = elMat%All(:,ind_loc(fas(j),perm),els(j))
           elMat%fh(ind_loc(fas(j),:),els(j))      = elMat%fh(ind_loc(fas(j),perm),els(j))
        END if
     END DO
  END DO
!!$OMP END DO
!!$OMP END PARALLEL 



 uef = 0.
 qef = 0.

!************************************
!   Loop in exterior faces in 2D
!************************************ 
!$OMP PARALLEL DEFAULT(SHARED) & 
!$OMP PRIVATE(iface,iel,ifa,isdir,inde,indf,Xfl,i,j,qef,uef,uf,tau_save_el,xy_g_save_el,indtausave) 
!$OMP DO SCHEDULE(STATIC) 
  DO iface = 1,Mesh%Nextfaces

     iel = Mesh%extfaces(iface,1)
     ifa = Mesh%extfaces(iface,2)
     isdir =  Mesh%Fdir(iel,ifa)
     
     ! Coordinates of the nodes of the face
     Xfl = Mesh%X(Mesh%T(iel,refElPol%face_nodes(ifa,:)),:)
     
     ! Face solution
     indf = (iface+Mesh%Nintfaces-1)*Npfl + (/(i,i=1,Npfl)/)
     uf = lres(indf,:)
     
     ! Elements solution
     inde = (iel-1)*Npel + (/(i,i=1,Npel)/)
     uef(:,:,1) = ures( inde( refElPol%face_nodes(ifa,:) ),:)
     qef(:,:,1) = qres( inde( refElPol%face_nodes(ifa,:) ),:)
     
     ! Compute the matrices for the element
     CALL elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,qef(:,:,1),uef(:,:,1),uf,tau_save_el,xy_g_save_el)

#ifdef PARALL
     ! Flip faces 
      IF (Mesh%flipface(iel,ifa) ) then
         elMat%Alq(ind_loc(ifa,:),:,iel)   = elMat%Alq(ind_loc(ifa,perm),:,iel)
         elMat%Alu(ind_loc(ifa,:),:,iel)   = elMat%Alu(ind_loc(ifa,perm),:,iel)
         elMat%All(ind_loc(ifa,:),:,iel)   = elMat%All(ind_loc(ifa,perm),:,iel)
         elMat%All(:,ind_loc(ifa,:),iel)   = elMat%All(:,ind_loc(ifa,perm),iel)
         elMat%fh(ind_loc(ifa,:),iel)      = elMat%fh(ind_loc(ifa,perm),iel)
      END if
#endif
     if (save_tau) then
        indtausave = (iel-1)*Nfre*Ng1d + (/(i,i =1,Nfre*Ng1d)/)
        tau_save(indtausave,:)= tau_save_el
        xy_g_save(indtausave,:)= xy_g_save_el
     endif
  END DO
!$OMP END DO
!$OMP END PARALLEL 

  deallocate(ures,lres,u0res)
#ifdef TEMPERATURE
  deallocate(qres)
#endif

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
! Volume computation in 2D 
!***************************************************
					SUBROUTINE elemental_matrices_volume(iel,Xel,qe,ue,u0e)
       integer,intent(IN)			     :: iel 
					  real*8,intent(IN)         :: Xel(:,:)
					  real*8,intent(IN)         :: qe(:,:)
       real*8,intent(IN)         :: ue(:,:),u0e(:,:,:)
							integer*4                 :: g,NGauss,i,j,k
							real*8                    :: dvolu
							real*8                    :: xy(Ng2d,ndim),ueg(Ng2d,neq),u0eg(Ng2d,neq,time%tis)
       real*8                    :: force(Ng2d,Neq)
       real*8                    :: qeg(Ng2d,neq*Ndim)
       real*8                     :: J11(Ng2d),J12(Ng2d)
       real*8                     :: J21(Ng2d),J22(Ng2d)
       real*8                     :: detJ(Ng2d)
       real*8                     :: iJ11(Ng2d),iJ12(Ng2d)
       real*8                     :: iJ21(Ng2d),iJ22(Ng2d)
       real*8                     :: fluxg(Ng2d)
							integer*4,dimension(Npel)   :: ind_ass,ind_asg,ind_asq
							real*8                    :: ktis(time%tis+1)
							real*8,dimension(Npel)      :: Ni,Nxg,Nyg,NNbb
							real*8,dimension(Npel,Npel)   :: NxNi,NyNi,NNxy,NNi
       real*8                    :: NxNy(Npel,Npel,Ndim)
							
       real*8                    :: upg(Ng2d,phys%npv)
       real*8                    :: b(Ng2d,2),divb(Ng2d),drift(Ng2d,2),divbg,driftg(2)
       real*8                    :: diff_iso_vol(Neq,Neq,Ng2d),diff_ani_vol(Neq,Neq,Ng2d)
#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
       real*8                    :: Bmod_g(Ng2d),Bmod_x,Bmod_y,Nx_ax(Npel)
       real*8                    :: Bt(Mesh%Nnodesperelem),Bt_g(Ng2d)              
#endif
	       
       ind_ass = (/ (i, i = 0, Neq*(Npel-1),Neq) /)
       ind_asg = (/ (i, i = 0, Neq*(Npel-1)*Ndim,Neq*Ndim) /)
       ind_asq = (/ (i, i = 0, Neq*(Npel-1)*Ndim,Neq*Ndim) /)

       !***********************************
       !    Volume computation
       !***********************************
							
							! Gauss points position
							xy = matmul(refElPol%N2D,Xel)

							! Compute diffusion at Gauss points
							CALL setLocalDiff(xy,diff_iso_vol,diff_ani_vol)

							! Solution at Gauss points
							ueg = matmul(refElPol%N2D,ue)
							qeg = matmul(refElPol%N2D,qe) 

       ! Solution at previous time steps, at Gauss points
       do i=1,time%tis
          u0eg(:,:,i) = matmul(refElPol%N2D,u0e(:,:,i))
       end do

			    ! Physical variables at Gauss points
			    CALL cons2phys(ueg,upg)
       
       ! Constant sources
							! Body force at the integration points
							CALL body_force(xy(:,1),xy(:,2),force)			

							!! Some sources for West cases
							IF (switch%testcase	.ge.51 .and. switch%testcase	.le.54) THEN						
							   fluxg = matmul(refElPol%N2D,phys%flux2D(Mesh%T(iel,:)))
							   DO g=1,Ng2d
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
        
          ! Drift at Gauss points       
          Bmod_g = matmul(refElPol%N2D,Bmod) ! Ngauss
          Bt_g   = matmul(refElPol%N2D,phys%Bt(Mesh%T(iel,:)))  ! Ngauss 
       END IF       
#endif          
#endif		
!************ END SELECTION MOVING EQUILIBRIUM **********************

							! Loop in 2D Gauss points
       Ngauss = Ng2d
       J11 = matmul(refElPol%Nxi2D,Xel(:,1))                           ! ng x 1
       J12 = matmul(refElPol%Nxi2D,Xel(:,2))                           ! ng x 1
       J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
       J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
				   iJ11 =  J22/detJ
				   iJ12 = -J12/detJ
				   iJ21 = -J21/detJ
				   iJ22 =  J11/detJ
				                 
       ! Coefficient time integration scheme
       call setTimeIntegrationCoefficients(ktis)

      ! Loop in 2D Gauss points
							DO g = 1,NGauss
							
										! Integration weight
										dvolu = refElPol%gauss_weights2D(g)*detJ(g)
          IF (switch%axisym) THEN										
             dvolu = dvolu*xy(g,1)
          END IF

										! x and y derivatives of the shape functions
						    Nxg = iJ11(g)*refElPol%Nxi2D(g,:) + iJ12(g)*refElPol%Neta2D(g,:)
						    Nyg = iJ21(g)*refElPol%Nxi2D(g,:) + iJ22(g)*refElPol%Neta2D(g,:)     										

									 ! Shape functions products
									 Ni  = refElPol%N2D(g,:)*dvolu
          NNi = tensorProduct(Ni,refElPol%N2D(g,:))                        ! Npel x Npel 
									 NxNi = tensorProduct(Nxg,Ni)                                     ! Npel x Npel 
          NyNi = tensorProduct(Nyg,Ni)                                     ! Npel x Npel 
          NNxy = b(g,1)*NxNi+b(g,2)*NyNi			                                ! Npel x Npel 
          NxNy(:,:,1)    = NxNi 
          NxNy(:,:,2)    = NyNi                                            ! Npel x Npel x 2
          NNbb =  (Nxg*b(g,1)+Nyg*b(g,2))*dvolu                             ! Npel x 1

          ! Divergence of b at the Gauss points  
#ifndef MOVINGEQUILIBRIUM
          divbg = divb(g)
						    driftg(1) = drift(g,1)
						    driftg(2) = drift(g,2)
#else      
							   IF (switch%axisym) THEN
							      Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
							   ELSE
							      Nx_ax = Nxg
							   END IF  
          divbg = dot_product(Nx_ax,b_nod(:,1))+dot_product(Nyg,b_nod(:,2))    
									 Bmod_x  = dot_product(Nxg,Bmod) ! 1
									 Bmod_y  = dot_product(Nyg,Bmod) ! 1
									 driftg(1) =  phys%dfcoef*Bt_g(g)*Bmod_y/Bmod_g(g)**3 ! I multiply for the temperature later
									 driftg(2) = -phys%dfcoef*Bt_g(g)*Bmod_x/Bmod_g(g)**3 ! I multiply for the temperature later             
#endif 


          CALL assemblyVolumeContribution(iel,ind_ass,ind_asg,ind_asq,b(g,:),divbg,driftg,force(g,:),&
          &ktis,diff_iso_vol(:,:,g),diff_ani_vol(:,:,g),Ni,NNi,Nxg,Nyg,NNxy,NxNy,NNbb,upg(g,:),ueg(g,:),qeg(g,:),u0eg(g,:,:))
							END DO ! END loop in volume Gauss points
					END SUBROUTINE elemental_matrices_volume








!***************************************************
! Interior faces computation in 2D
!***************************************************
					SUBROUTINE elemental_matrices_faces_int(els,fas,Xfl,qef,uef,uf,tau_save_el,xy_g_save_el)
       integer,intent(IN)			     :: els(:),fas(:) 
					  real*8,intent(IN)         :: Xfl(:,:)
					  real*8,intent(IN)         :: qef(:,:,:)
       real*8,intent(IN)         :: uef(:,:,:),uf(:,:)	
       real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)

							integer*4                 :: g,NGauss,ifa,i,j,k,lel,iel,ieln
							real*8                    :: dline,xyDerNorm_g
							real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
							real*8                    :: xyf(Ng1d,ndim)
							real*8                    :: xyDer(Ng1d,ndim)
       real*8                    :: qfg(Ng1d,neq*Ndim)
							integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
       integer*4,dimension(Npfl)  :: ind_asf,ind_ash
							real*8                    :: t_g(ndim),n_g(ndim), bn
							real*8                    :: NNif(Npfl,Npfl),Nfbn(Npfl)
       real*8                    :: upgf(Ng1d,phys%npv)
       real*8                    :: tau(Neq,Neq)
													
       real*8                    :: b_f(Ng1d,2)				
       real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)

#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
       real*8                    :: Bmod_g(Ng2d),Bmod_x,Bmod_y,Nx_ax(Npel)
       real*8                    :: Bt(Mesh%Nnodesperelem),Bt_g(Ng2d)              
#endif

       
       if (save_tau) then
          tau_save_el = 0.
          xy_g_save_el = 0.
       endif 
       ind_asf = (/(i,i=0,Neq*(Npfl-1),Neq)/)
       ind_ash = (/(i,i=0,Neq*(Npfl-1)*Ndim,Neq*Ndim)/)

       !***********************************
  					! Faces computations
       !***********************************
							NGauss = Ng1d
												
       ! Loop in the faces of the element
							DO  lel = 1,2
										
          iel = els(lel)
          ifa = fas(lel)
          if (lel==1) then
             ieln = els(2)
          else
             ieln = els(1)
          endif

										! Indices
          ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Npfl/))
          ind_ff = (ifa-1)*neq*Npfl + (/ (i, i = 1, neq*Npfl) /)										
	         ind_fg = reshape(tensorSumInt( (/ (i, i = 1, neq*ndim) /),neq*ndim*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Npfl*ndim/))
       
          ! Trace solution at face Gauss points
										IF (Mesh%flipFace(iel,ifa)) THEN
														ufg = matmul(refElPol%N1D,uf((/(i,i=Npfl,1,-1)/),:))
              xyf = matmul(refElPol%N1D,Xfl((/(i,i=Npfl,1,-1)/),:))
														! Shape function derivatives at Gauss points
              xyDer  = matmul(refElPol%Nxi1D,Xfl((/(i,i=Npfl,1,-1)/),:))
														
          ELSE
 													ufg = matmul(refElPol%N1D,uf)
              xyf = matmul(refElPol%N1D,Xfl)
														! Shape function derivatives at Gauss points
														xyDer  = matmul(refElPol%Nxi1D,Xfl)
										END IF       
										
          ! Element solution at face Gauss points										
          uefg = matmul(refElPol%N1D,uef(:,:,lel))
          ! Gradient solution at face gauss points
          qfg    = matmul(refElPol%N1D,qef(:,:,lel))
          
         ! Compute diffusion at faces Gauss points
							  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

         ! Physical variables at face Gauss points
         CALL cons2phys(ufg,upgf)

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
										   									   
												! Shape functions products
												 bn = dot_product(b_f(g,:),n_g)
													NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
 				        Nfbn =  bn*refElPol%N1D(g,:)*dline
 				        
 				        
 				        ! Compute the stabilization term			
 				        tau = 0.			   			
             IF (numer%stab==1) THEN
           		! Constant stabilization
						          DO i=1,Neq
						             tau(i,i) = numer%tau(i)
						          END DO
             ELSE
             ! Non constant stabilization
                ! Compute tau in the Gauss points
                IF (numer%stab<5) THEN
                   CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g,iel,ifa,0.,xyf(g,:),tau)
                ELSE
                   CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g,xyf(g,:),0.,iel,tau)
                ENDIF
             END IF

             ! Assembly local contributions              				        										   									   
             CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
             n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau) 

             if (save_tau) then
                DO i=1,Neq
                   tau_save_el((ifa-1)*Ngauss+g,i) = tau(i,i)
                END DO
                xy_g_save_el((ifa-1)*Ngauss+g,:) = xyf(g,:)
             endif 
                        
										END DO ! Gauss points
							END DO ! 2 elements

					END SUBROUTINE elemental_matrices_faces_int	
										





			
!***************************************************
! Exterior faces computation in 2D
!***************************************************
					SUBROUTINE elemental_matrices_faces_ext(iel,ifa,isdir,Xfl,qef,uef,uf,tau_save_el,xy_g_save_el)
       integer,intent(IN)			     :: iel ,ifa
					  real*8,intent(IN)         :: Xfl(:,:)
					  logical,intent(IN)        :: isdir
					  real*8,intent(IN)         :: qef(:,:)
       real*8,intent(INOUT)         :: uef(:,:),uf(:,:)
       real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)
	
							integer*4                 :: g,NGauss,i,j,k
							real*8                    :: dline,xyDerNorm_g
							real*8                    :: ufg(Ng1d,neq),uefg(Ng1d,neq)
							real*8                    :: xyf(Ng1d,ndim)
							real*8                    :: xyDer(Ng1d,ndim)
       real*8                    :: qfg(Ng1d,neq*Ndim)
							integer*4                 :: ind_ff(Neq*Npfl),ind_fe(Neq*Npfl),ind_fg(Neq*Ndim*Npfl)
       integer*4,dimension(Npfl)  :: ind_asf,ind_ash
       real                      :: isext
							real*8                    :: t_g(ndim),n_g(ndim), bn
							real*8                    :: NNif(Npfl,Npfl),Nfbn(Npfl)
							real*8                    :: tau(Neq,Neq) 
       real*8                    :: upgf(Ng1d,phys%npv)
													
       real*8                    :: b_f(Ng1d,2)				
       real*8                    :: diff_iso_fac(Neq,Neq,Ng1d),diff_ani_fac(Neq,Neq,Ng1d)

#ifdef MOVINGEQUILIBRIUM
       real*8                    :: b_nod(Mesh%Nnodesperelem,ndim),Bmod(Mesh%Nnodesperelem)
       real*8                    :: Bmod_g(Ng2d),Bmod_x,Bmod_y,Nx_ax(Npel)
       real*8                    :: Bt(Mesh%Nnodesperelem),Bt_g(Ng2d)              
#endif
       


       if (save_tau) then
          tau_save_el = 0.
          xy_g_save_el = 0.
       endif 
       ind_asf = (/ (i, i = 0, Neq*(Npfl-1),Neq) /)
       ind_ash = (/ (i, i = 0, Neq*(Npfl-1)*Ndim,Neq*Ndim) /)

       !***********************************
  					! Faces computations
       !***********************************
							NGauss = Ng1d
									

							! Indices
       ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Npfl/))
       ind_ff = (ifa-1)*neq*Npfl + (/ (i, i = 1, neq*Npfl) /)										
       ind_fg = reshape(tensorSumInt( (/ (i, i = 1, neq*ndim) /),neq*ndim*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Npfl*ndim/))
    
       ! Trace solution at face Gauss points
       xyf = matmul(refElPol%N1D,Xfl)
       IF (isdir) THEN
          CALL analytical_solution(xyf(:,1),xyf(:,2),ufg) 
       ELSE
#ifdef PARALL
						    IF (Mesh%flipFace(iel,ifa)) THEN
						       uf = uf((/(i,i=Npfl,1,-1)/),:)
						    ENDIF
#endif
							   ufg = matmul(refElPol%N1D,uf)
       END IF
       
							
       ! Element solution at face Gauss points										
       uefg = matmul(refElPol%N1D,uef)
       ! Gradient solution at face gauss points
       qfg    = matmul(refElPol%N1D,qef)
          
      ! Compute diffusion at faces Gauss points
				  CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac)

      ! Physical variables at face Gauss points
      CALL cons2phys(ufg,upgf)

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
							xyDer  = matmul(refElPol%Nxi1D,Xfl)
							
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
							   									   
							   								   
		        ! Compute the stabilization term		
		        isext = 1.
#ifdef PARALL    
          IF (Mesh%boundaryFlag(Mesh%F(iel,ifa)-Mesh%Nintfaces).eq.0) THEN
             isext = 0.
          END IF
#endif		        
		        tau = 0.				   			
          IF (numer%stab==1) THEN
        		! Constant stabilization
			          DO i=1,Neq
			             tau(i,i) = numer%tau(i)
			          END DO
          ELSE
          ! Non constant stabilization
             ! Compute tau in the Gauss points
						       IF (numer%stab<5) THEN
						          CALL computeTauGaussPoints(upgf(g,:),ufg(g,:),b_f(g,:),n_g,iel,ifa,isext,xyf(g,:),tau)
						       ELSE
						          CALL computeTauGaussPoints_matrix(upgf(g,:),ufg(g,:),b_f(g,:),n_g,xyf(g,:),isext,iel,tau)
						       ENDIF
          END IF
          							  
									! Shape functions products
									 bn = dot_product(b_f(g,:),n_g)
										NNif = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
		        Nfbn =  bn*refElPol%N1D(g,:)*dline

         ! Assembly local contributions
#ifdef PARALL
         IF (Mesh%boundaryFlag(Mesh%F(iel,ifa)-Mesh%Nintfaces).eq.0) THEN
            ! Ghost face: assembly it as interior
            CALL assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
            n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau)          
         ELSE
            CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
            n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau)         
         ENDIF
#else
          CALL assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b_f(g,:),&
          n_g,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),NNif,Nfbn,ufg(g,:),qfg(g,:),tau)
#endif         
          if (save_tau) then
             DO i=1,Neq
                tau_save_el((ifa-1)*Ngauss+g,i) = tau(i,i)
             END DO
             xy_g_save_el((ifa-1)*Ngauss+g,:) = xyf(g,:)
          endif  
                       
							END DO ! Gauss points
										

					END SUBROUTINE elemental_matrices_faces_ext		




			
!*******************************************
!           AUXILIARY ROUTINES
!*******************************************				
     
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
							 
							 
#endif							 
							 
							 
							 
							 
							 
							 
							 
							 
!********************************************************************************************							 
! 
!							 
!                             ROUTINES FOR 2D AND 3D COMPUTATIONS 
! 
! 
!********************************************************************************************							 
        SUBROUTINE setTimeIntegrationCoefficients(ktis)			
        real*8,intent(out) :: ktis(:)
        integer :: it

        ktis = 0.


        if (time%it.lt.time%tis) then
            it = time%it
        else
            it=time%tis
        end if

        SELECT CASE(it)
            case(1)
                ktis(1) = 1.
                ktis(2) = 1.
            case(2)
                ktis(1) = 1.5
                ktis(2) = 2
                ktis(3) = -0.5
            case(3)
                ktis(1) = 11./6.
                ktis(2) = 3.
                ktis(3) = -1.5
                ktis(4) = 1./3.
            case(4)
                ktis(1) = 25./12.
                ktis(2) = 4.
                ktis(3) = -3.
                ktis(4) = 4./3.
                ktis(4) = -0.25
            case(5)
                ktis(1) = 137./60.
                ktis(2) = 5.
                ktis(3) = -5.
                ktis(4) = 10./3.
                ktis(5) = -1.25 
                ktis(6) = 0.2               
            case(6)
                ktis(1) = 147./60.
                ktis(2) = 6.
                ktis(3) = -7.5
                ktis(4) = 20./3.
                ktis(5) = -3.75 
                ktis(6) = 1.2
                ktis(7) = -1.6
            case default
                write(6,*) 'Formula not available'
                stop
        END SELECT
        END SUBROUTINE setTimeIntegrationCoefficients
				       
				       
!********************************************************************				   
! 
!         ASSEMBLY VOLUME CONTRIBUTION
! 
!********************************************************************		
				   SUBROUTINE assemblyVolumeContribution(iel,ind_ass,ind_asg,ind_asq,b,divb,drift,f,ktis,diffiso,diffani,Ni,NNi,Nxg,Nyg,NNxy,NxNy,NNbb,upe,ue,qe,u0e) 
				   integer*4,intent(IN)      :: iel,ind_ass(:),ind_asg(:),ind_asq(:)    
       real*8,intent(IN)         :: b(:),divb,drift(:),f(:),ktis(:)
       real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
       real*8,intent(IN)         :: Ni(:),NNi(:,:),Nxg(:),Nyg(:),NNxy(:,:),NxNy(:,:,:),NNbb(:)
       real*8,intent(IN)         :: upe(:),ue(:),u0e(:,:)
       real*8,intent(IN)         :: qe(:)			
       integer*4                 :: i,j,k,iord,ii
       integer*4,dimension(Npel) :: ind_i,ind_j,ind_k
       real*8,dimension(neq,neq) :: A
       real*8                    :: nlcorr
       real*8                    :: Qpr(Ndim,Neq)
#ifdef TEMPERATURE
							real*8,dimension(neq,neq) :: GG
							real*8                    :: Telect
							real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
       real*8                    :: W,dW_dU(Neq,Neq),s,ds_dU(Neq),Zet(ndim,Neq)
#endif       
         ! Jacobian for convection term
         CALL jacobianMatrices(ue,A)   

         ! Compute Q^T^(k-1)
         Qpr = reshape(qe,(/Ndim,Neq/))
  

#ifdef TEMPERATURE
         ! Jacobian for the curvature term
         CALL GimpMatrix(ue,divb,GG)

         ! Compute V(U^(k-1)) 
         call computeVi(ue,Vveci)
         call computeVe(ue,Vvece)

         ! Compute dV_dU (k-1)
         call compute_dV_dUi(ue,dV_dUi)
         call compute_dV_dUe(ue,dV_dUe)

         ! Compute Alpha(U^(k-1))
         Alphai = computeAlphai(ue)
         Alphae = computeAlphae(ue)

         ! Compute dAlpha/dU^(k-1)
         call compute_dAlpha_dUi(ue,dAlpha_dUi)
         call compute_dAlpha_dUe(ue,dAlpha_dUe)
         
         gmi = dot_product(matmul(Qpr,Vveci),b)             ! scalar     
         gme = dot_product(matmul(Qpr,Vvece),b)             ! scalar         
         Taui = matmul(Qpr,dV_dUi)      ! Ndim x Neq
         Taue = matmul(Qpr,dV_dUe)      ! Ndim x Neq
			      
			      
         ! Parallel current term
			      ! Compute W(U^(k-1)) 
			      call compute_W(ue,W)						      
			      ! Compute dW_dU(U^(k-1)) 
			      call compute_dW_dU(ue,dW_dU)	
			      		

         ! Temperature exchange terms				      
			      ! s(U^(k-1)) 
			      call compute_S(ue,s)		
			      ! ds_du(U^(k-1)) 
			      call compute_dS_dU(ue,ds_dU)	
   						      						      
         Zet  = matmul(Qpr,dW_dU)       ! Ndim x Neq

#endif

          ! Assembly local matrix
          ! Loop in equations
          DO i=1,Neq
             ind_i = i+ind_ass
             IF (.not.switch%steady) THEN             
                ! Time derivative contribution
#ifdef VORTICITY
             ! In the vorticity model I don't assemble the mass matrix for the potential equation
             IF (i .ne. 4) THEN 
#endif                
                elMat%Auu(ind_i,ind_i,iel) = elMat%Auu(ind_i,ind_i,iel)+ktis(1)*NNi/time%dt
#ifdef VORTICITY
             END IF
#endif                
             END IF
#ifndef TEMPERATURE
             IF (i==2) THEN
                ! Curvature contribution (isothermal)
                elMat%Auu(ind_i,ind_i-1,iel) = elMat%Auu(ind_i,ind_i-1,iel) -phys%a*divb*NNi
             END IF
             IF (switch%driftvel) THEN
                ! B x GradB drift (isothermal)
                elMat%Auu(ind_i,ind_i,iel) = elMat%Auu(ind_i,ind_i,iel) + TensorProduct(Ni,( Nxg*drift(1)+ Nyg*drift(2) ))
             END IF
#else
						       IF (i==2) THEN
						          DO j=1,Neq
                   ind_j = j+ind_ass
                   ! Curvature contribution (non-isothermal)
						             elMat%Auu(ind_i,ind_j,iel) = elMat%Auu(ind_i,ind_j,iel) - GG(i,j)*NNi 
						          END DO
						       END IF
                    
             IF (switch%driftvel) THEN               
									       Telect = upe(8)
                ! B x GradB drift (non-isothermal)
												    elMat%Auu(ind_i,ind_i,iel) = elMat%Auu(ind_i,ind_i,iel) + TensorProduct(Ni,Telect*( Nxg*drift(1)+ Nyg*drift(2) ))
             END IF

             ! Parallel diffusion for the temperature
             IF (i==3) THEN
                DO j=1,4
                   ind_j = j+ind_ass
                   elMat%Auu(ind_i,ind_j,iel) = elMat%Auu(ind_i,ind_j,iel) + &
                   &coefi*(gmi*dAlpha_dUi(j)+Alphai*(dot_product(Taui(:,j),b)))*NNxy + ( dot_product(Zet(:,j),b) +ds_dU(j))*NNi   
                   DO k=1,Ndim
                      ind_k = k+(j-1)*Ndim+ind_asg
                      elMat%Auq(ind_i,ind_k,iel) = elMat%Auq(ind_i,ind_k,iel) + coefi*Alphai*Vveci(j)*b(k)*NNxy   
                      IF (j==4) THEN
                         elMat%Auq(ind_i,ind_k,iel) = elMat%Auq(ind_i,ind_k,iel) + W*NNi*b(k)                      
                      END IF
                   END DO
                END DO 
                elMat%S(ind_i,iel) = elMat%S(ind_i,iel) + coefi*Alphai*( dot_product (matmul(transpose(Taui),b),ue) )*NNbb + s*Ni
             ELSEIF (i==4) THEN
                DO j=1,4
                   ind_j = j+ind_ass
                   elMat%Auu(ind_i,ind_j,iel) = elMat%Auu(ind_i,ind_j,iel) + &
                   &coefe*(gme*dAlpha_dUe(j)+Alphae*(dot_product(Taue(:,j),b)))*NNxy - (dot_product(Zet(:,j),b)+ds_dU(j))*NNi   
                   DO k=1,Ndim
                      ind_k = k+(j-1)*Ndim+ind_asg
                      elMat%Auq(ind_i,ind_k,iel) = elMat%Auq(ind_i,ind_k,iel) + coefe*Alphae*Vvece(j)*b(k)*NNxy      
                      IF (j==4) THEN
                         elMat%Auq(ind_i,ind_k,iel) = elMat%Auq(ind_i,ind_k,iel) - W*NNi*b(k)                    
                      END IF     
                   END DO                              
                END DO  
                elMat%S(ind_i,iel) =  elMat%S(ind_i,iel) + coefe*Alphae*( dot_product (matmul(transpose(Taue),b),ue)  )*NNbb - s*Ni
             END IF           
#endif
             ! Convection contribution
             DO j=1,Neq
                ind_j = j+ind_ass
                elMat%Auu(ind_i,ind_j,iel) = elMat%Auu(ind_i,ind_j,iel) - A(i,j)*NNxy
             END DO

             ! Perpendicular diffusion contribution
							      DO k = 1,Ndim
							         ind_k = ind_asq+k+(i-1)*Ndim 
							         ind_j = ind_ass+i
                ! Diagonal terms for perpendicular diffusion
                elMat%Auq(ind_j,ind_k,iel) = elMat%Auq(ind_j,ind_k,iel) + diffiso(i,i)*NxNy(:,:,k) - diffani(i,i)*NNxy*b(k) 

#ifdef VORTICITY
                ! Non-diagonal terms for perpendicular diffusion
                DO ii=1,Neq
                   IF (ii==i) CYCLE ! diagonal alredy assembled
                   IF (abs(diffiso(i,ii))<1e-12 .and. abs(diffani(i,ii))<1e-12) CYCLE
                   nlcorr = 1.
                   ! Non-linear correction for non-linear diffusive terms.
                   ! TODO: find a smarter way to include it, avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)
                   IF (i==3 .and. ii==1) then
                      nlcorr = 1./ue(1)
                   ENDIF
                   ind_k = ind_asq+k+(ii-1)*Ndim 
                   elMat%Auq(ind_j,ind_k,iel) = elMat%Auq(ind_j,ind_k,iel) + nlcorr*diffiso(i,ii)*NxNy(:,:,k) - nlcorr*diffani(i,ii)*NNxy*b(k)
                END DO
#endif
										   END DO
										   
#ifdef VORTICITY
                ! Non-linear term in the vorticity equation (\Grad// n/n b)
                IF (i==3) THEN
                   j=1
                   ind_j = ind_ass+j
                   elMat%Auu(ind_i,ind_j,iel) = elMat%Auu(ind_i,ind_j,iel) - (dot_product(Qpr(:,j),b))/(ue(1)**2)*NNxy
                   elMat%S(ind_i,iel) = elMat%S(ind_i,iel) - dot_product(Qpr(:,j),b )/(ue(1))*NNbb
                ENDIF
                
                ! The vorticity is the source term in the potential equation
                IF (i==4) THEN
                   elMat%S(ind_i,iel) = elMat%S(ind_i,iel) + ue(3)*Ni
                ENDIF
#endif										   
          END DO ! Loop in equations
           
          ! Assembly RHS
          IF (.not.switch%steady) THEN
             DO iord =1,time%tis
                ! Time derivative contribution
                elMat%S(:,iel) = elMat%S(:,iel)+ktis(iord+1)*col(tensorProduct(u0e(:,iord),Ni))/time%dt
             END DO
          END IF
          ! Linear body force contribution
          elMat%S(:,iel) = elMat%S(:,iel)+col(tensorProduct(f,Ni))
          END SUBROUTINE assemblyVolumeContribution
          
          
          
          
          
          
          
          
!********************************************************************				   
! 
!         ASSEMBLY INTERIOR FACES CONTRIBUTION
! 
!********************************************************************	        
       SUBROUTINE assemblyIntFacesContribution(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b,n,diffiso,diffani,NNif,Nfbn,uf,qf,tau,ifa)
				   integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
       real*8,intent(IN)         :: b(:),n(:)
       real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
       real*8,intent(IN)         :: NNif(:,:),Nfbn(:)
       real*8,intent(IN)         :: uf(:)
       real*8,intent(IN)         :: qf(:)		
       real*8,optional,intent(IN) :: tau(:,:)
       real*8                     :: nlcorr
       integer*4,optional         :: ifa       	
       integer*4                  :: i,j,k,ii
       integer*4,dimension(size(ind_asf))  :: ind_if,ind_jf,ind_kf
       real*8,dimension(neq,neq) :: A
       real*8                    :: bn,kmult(size(ind_asf),size(ind_asf)), kmultf(size(ind_asf))
       real*8                    :: Qpr(Ndim,Neq)
#ifdef TEMPERATURE
							real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
#endif            



			          ! Jacobian matrices
			          bn = dot_product(b,n)
             CALL jacobianMatrices(uf,A)  

								     ! Compute Q^T^(k-1)
								     Qpr = reshape(qf,(/Ndim,Neq/))

#ifdef TEMPERATURE


								     ! Compute V(U^(k-1)) 
								     call computeVi(uf,Vveci)
								     call computeVe(uf,Vvece)

								     ! Compute dV_dU (k-1)
								     call compute_dV_dUi(uf,dV_dUi)
								     call compute_dV_dUe(uf,dV_dUe)

								     ! Compute Alpha(U^(k-1))
								     Alphai = computeAlphai(uf)
								     Alphae = computeAlphae(uf)

								     ! Compute dAlpha/dU^(k-1)
								     call compute_dAlpha_dUi(uf,dAlpha_dUi)
								     call compute_dAlpha_dUe(uf,dAlpha_dUe)

				         gmi = dot_product(matmul(Qpr,Vveci),b)             ! scalar    
				         gme = dot_product(matmul(Qpr,Vvece),b)  
				         Taui = matmul(Qpr,dV_dUi)      ! 2x3
				         Taue = matmul(Qpr,dV_dUe)      ! 2x3
				         
#endif

             ! Assembly local matrix
             DO i=1,Neq
                ind_if = ind_asf+i
						 						   DO k = 1,Ndim
						 						      ind_kf = ind_ash+k+(i-1)*Ndim 
                   kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))

 						            ! Diagonal terms for perpendicular diffusion
 						            elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult 		 						            
 						            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult 						 
#ifdef VORTICITY
                   DO ii=1,Neq
                      IF (ii==i) CYCLE ! diagonal alredy assembled
                      IF (abs(diffiso(i,ii))<1e-12 .and. abs(diffani(i,ii))<1e-12) CYCLE                
				 		               ind_kf = ind_ash+k+(ii-1)*Ndim
				 		               nlcorr = 1.
                      ! Non-linear correction for non-linear diffusive terms.
                      ! TODO: find a smarter way to include it, avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)
                      IF (i==3 .and. ii==1) then
                         nlcorr = 1./uf(1)
                      ENDIF				 		               
                      kmult = NNif*nlcorr*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii) )
 						               elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fG(ind_kf),iel) - kmult 					 		               
					 	               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
						             END DO
#endif
								 				   END DO ! k-loop
						       
                ! Convection contribution
                DO j=1,Neq
                   ind_jf = ind_asf+j
                   kmult = bn*A(i,j)*NNif
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult	  
                   elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult    
                END DO ! j-loop
                             
                
#ifdef TEMPERATURE
										      ! Parallel diffusion for the temperature
							         IF (i==3) THEN
							            DO j=1,4
                      ind_jf = ind_asf+j
                      kmult = coefi*(gmi*dAlpha_dUi(j)+Alphai*(  dot_product(Taui(:,j),b)  ))*NNif*bn                      
							               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
							               elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)  = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
							               DO k=1,Ndim 
                         ind_kf = k+(j-1)*Ndim+ind_ash
                         kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
					 		                 elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel)  - kmult
						 	                 elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel)  - kmult
							               END DO
							            END DO 
                   kmultf = coefi*Alphai*(  dot_product (matmul(transpose(Taui),b),uf)   )*Nfbn
							            elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) - kmultf
                   elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
							         ELSEIF (i==4) THEN
							            DO j=1,4
                      ind_jf = ind_asf+j
                      kmult = coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b) ))*NNif*bn
							               elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)   = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
							               elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel)   = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
							               DO k=1,Ndim
                         ind_kf = k+(j-1)*Ndim+ind_ash
                         kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
							                  elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) -kmult
							                  elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) = elMat%Alq(ind_ff(ind_if),ind_fg(ind_kf),iel) -kmult
							               END DO
							            END DO  
                   kmultf = coefe*Alphae*(dot_product (matmul(transpose(Taue),b),uf))*Nfbn
							            elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) - kmultf
                   elMat%fh(ind_ff(ind_if),iel) = elMat%fh(ind_ff(ind_if),iel) - kmultf
							         END IF     
#endif
#ifdef VORTICITY
                ! Non-linear term in the vorticity equation (\Grad// n/n b)
                IF (i==3) THEN
                   j=1
                   ind_jf = ind_asf+j
                   kmult = dot_product(Qpr(:,1),b )/uf(1)**2*NNif*bn
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult 
                   elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) + kmult
                   elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) + dot_product(Qpr(:,1),b )/uf(1)*Nfbn
                ENDIF 
#endif			

             END DO  ! i-Loop
          		

            
            
							      ! Assembly stabilization terms
             IF (numer%stab<5) THEN
                DO i=1,Neq
                   ind_if = i+ind_asf
                   kmult = tau(i,i)*NNif
                   elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
                   elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_if),iel) - kmult
                   elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_if),iel) + kmult
                END DO 										        
             ELSE
                DO i=1,Neq
                   ind_if = i+ind_asf
                   DO j=1,Neq
                      ind_jf = j+ind_asf
                      kmult = tau(i,j)*NNif
                      elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
                      elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                      elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) = elMat%All(ind_ff(ind_if),ind_ff(ind_jf),iel) - kmult
                      elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) = elMat%Alu(ind_ff(ind_if),ind_fe(ind_jf),iel) + kmult
                   END DO
                END DO 
             ENDIF
            !************* End stabilization terms************************          
          
       END SUBROUTINE assemblyIntFacesContribution
          
          
                    
!********************************************************************				   
! 
!         ASSEMBLY EXTERIOR FACES CONTRIBUTION
! 
!********************************************************************	      
       SUBROUTINE assemblyExtFacesContribution(iel,isdir,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,b,n,diffiso,diffani,NNif,Nfbn,uf,qf,tau,ifa)
				   integer*4,intent(IN)      :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
				   logical                   :: isdir
       real*8,intent(IN)         :: b(:),n(:)
       real*8,intent(IN)         :: diffiso(:,:),diffani(:,:)
       real*8,intent(IN)         :: NNif(:,:),Nfbn(:)
       real*8,intent(IN)         :: uf(:)
       real*8,intent(IN)         :: qf(:)		
       real*8,optional,intent(IN) :: tau(:,:)
       integer*4,optional         :: ifa
       real*8                    :: nlcorr
       integer*4                 :: i,j,k,ii
       integer*4,dimension(Npfl)  :: ind_if,ind_jf,ind_kf
       real*8,dimension(neq,neq) :: A
       real*8                    :: bn,kmult(Npfl,Npfl), kmultf(Npfl)
       real*8                    :: Qpr(Ndim,Neq)
#ifdef TEMPERATURE
							real*8                    :: Vveci(Neq),dV_dUi(Neq,Neq),Alphai,dAlpha_dUi(Neq),gmi,taui(Ndim,Neq)
       real*8                    :: Vvece(Neq),dV_dUe(Neq,Neq),Alphae,dAlpha_dUe(Neq),gme,taue(Ndim,Neq)
#endif          
          
          ! Jacobian matrices
          bn = dot_product(b,n)
          CALL jacobianMatrices(uf,A)  

					     ! Compute Q^T^(k-1)
					     Qpr = reshape(qf,(/Ndim,Neq/))
					     
#ifdef TEMPERATURE

					     ! Compute V(U^(k-1)) 
					     call computeVi(uf,Vveci)
					     call computeVe(uf,Vvece)

					     ! Compute dV_dU (k-1)
					     call compute_dV_dUi(uf,dV_dUi)
					     call compute_dV_dUe(uf,dV_dUe)

					     ! Compute Alpha(U^(k-1))
					     Alphai = computeAlphai(uf)
					     Alphae = computeAlphae(uf)

					     ! Compute dAlpha/dU^(k-1)
					     call compute_dAlpha_dUi(uf,dAlpha_dUi)
					     call compute_dAlpha_dUe(uf,dAlpha_dUe)

	         gmi = dot_product(matmul(Qpr,Vveci),b)             ! scalar    
	         gme = dot_product(matmul(Qpr,Vvece),b)  
	         Taui = matmul(Qpr,dV_dUi)      ! 2x3
	         Taue = matmul(Qpr,dV_dUe)      ! 2x3
#endif
          ! Assembly local matrix
          DO i=1,Neq
             ind_if = ind_asf+i
			 						   DO k = 1,Ndim
			 						      ind_kf = ind_ash+k+(i-1)*Ndim 
                kmult = NNif*(n(k)*diffiso(i,i) - bn*b(k)*diffani(i,i))
				            ! Diagonal terms for perpendicular diffusion
				            elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult 						 
#ifdef VORTICITY
                   DO ii=1,Neq
                      IF (ii==i) CYCLE ! diagonal alredy assembled
                      IF (abs(diffiso(i,ii))<1e-12 .and. abs(diffani(i,ii))<1e-12) CYCLE                
				 		               ind_kf = ind_ash+k+(ii-1)*Ndim
				 		               nlcorr = 1.
                      ! Non-linear correction for non-linear diffusive terms.
                      ! TODO: find a smarter way to include it, avoiding if statements and model dependencies (i==3,ii==1 only holds for Isothermal+Vorticity model)
                      IF (i==3 .and. ii==1) then
                         nlcorr = 1./uf(1)
                      ENDIF					 		               
                      kmult = NNif*nlcorr*(n(k)*diffiso(i,ii) - bn*b(k)*diffani(i,ii) )
					 	               elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fG(ind_kf),iel) - kmult
						             END DO
#endif
					 				   END DO ! k-loop
									   
             ! Convection contribution
             IF (.not.isdir) THEN
                DO j=1,Neq
                   ind_jf = ind_asf+j
                   kmult = bn*A(i,j)*NNif
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult	  
                END DO ! j-loop
             ENDIF
#ifdef TEMPERATURE
							      ! Parallel diffusion for the temperature
				         IF (i==3) THEN
				            DO j=1,4
                   ind_jf = ind_asf+j
                   IF (.not.isdir) THEN
                      kmult = coefi*(gmi*dAlpha_dUi(j)+Alphai*(  dot_product(Taui(:,j),b)  ))*NNif*bn                        
				                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)  = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
				               END IF
				               DO k=1,Ndim 
                      ind_kf = k+(j-1)*Ndim+ind_ash                      
                      kmult = coefi*Alphai*Vveci(j)*b(k)*NNif*bn
		 		                 elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel)  - kmult
				               END DO
				            END DO 
                kmultf = coefi*Alphai*(  dot_product (matmul(transpose(Taui),b),uf)   )*Nfbn
				            elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) - kmultf
				         ELSEIF (i==4) THEN
				            DO j=1,4
                   ind_jf = ind_asf+j
                   IF (.not.isdir) THEN
                      kmult = coefe*(gme*dAlpha_dUe(j)+Alphae*( dot_product(Taue(:,j),b) ))*NNif*bn
				                  elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel)   = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
				               END IF
				               DO k=1,Ndim
                      ind_kf = k+(j-1)*Ndim+ind_ash
                      kmult = coefe*Alphae*Vvece(j)*b(k)*NNif*bn
				                  elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) = elMat%Auq(ind_fe(ind_if),ind_fg(ind_kf),iel) -kmult
				               END DO
				            END DO  
                kmultf = coefe*Alphae*(dot_product (matmul(transpose(Taue),b),uf))*Nfbn
				            elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) - kmultf
				         END IF           
#endif
#ifdef VORTICITY
                ! Non-linear term in the vorticity equation (\Grad// n/n b)
                IF (i==3) THEN
                   j=1
                   ind_jf = ind_asf+j
                   kmult = dot_product(Qpr(:,1),b )/uf(1)**2*NNif*bn
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) + kmult
                   elMat%S(ind_fe(ind_if),iel)  = elMat%S(ind_fe(ind_if),iel) + dot_product(Qpr(:,1),b )/uf(1)*Nfbn
                ENDIF 
#endif		
          END DO  ! i-Loop
       		
        
         
				      ! Assembly stabilization terms
          IF (numer%stab<5) THEN
             DO i=1,Neq
                ind_if = i+ind_asf
                kmult = tau(i,i)*NNif
                elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_if),iel) + kmult
                IF (.not.isdir) THEN
                   elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_if),iel) - kmult
                ENDIF
             END DO 										        
          ELSE
             DO i=1,Neq
                ind_if = i+ind_asf
                DO j=1,Neq
                   ind_jf = j+ind_asf
                   kmult = tau(i,j)*NNif
                   elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) = elMat%Auu(ind_fe(ind_if),ind_fe(ind_jf),iel) + kmult
                   IF (.not.isdir) THEN
                      elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) = elMat%Aul(ind_fe(ind_if),ind_ff(ind_jf),iel) - kmult
                   ENDIF
                END DO
             END DO 
          ENDIF
         !************* End stabilization terms************************             
       END SUBROUTINE assemblyExtFacesContribution          
          
          		       
END SUBROUTINE HDG_computeJacobian





