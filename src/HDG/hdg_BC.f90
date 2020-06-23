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
   integer*4                 :: ind_qf(refElTor%Nfl*3*phys%Neq)
   real*8                    :: qqf(refElTor%Nfl,3*phys%Neq),qfg(refElTor%Ngl,phys%neq*3)
   real*8                    :: coefi,coefe
   real*8                    :: tdiv(numer%ntor + 1)
   real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
   real*8                    :: xyf(refElPol%Ngauss1d,2),xyd_g(refElPol%Ngauss1d,2),xydNorm_g(refElPol%Ngauss1d)
   integer                   :: ind_asf(refElTor%Nfl),ind_ash(refElTor%Nfl)
   real*8                    :: Bfl(refElTor%Nfl,3)
   integer                   :: indbt(refElTor%Nfl)
   real*8                    :: Bmod_nod(refElTor%Nfl),b_nod(refElTor%Nfl,3)

   if (utils%timing) then
      call cpu_time(timing%tps1)
      call system_clock(timing%cks1,timing%clock_rate1)
   end if

   Ndim = 3                                               ! Number of dimensions
   Neq = phys%neq
   N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
   Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
   Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
   Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
   Npel = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
   Npfl = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
   Nfdir = Mesh%Ndir
   NGaussPol = refElPol%Ngauss1D
   NGaussTor = refElTor%Ngauss1D
   Ngfl = refElTor%Ngl

   ! Some indices
   ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
   ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

   ! Toroidal discretization
   htor = numer%tmax/numer%ntor
   tdiv = 0.
   DO i = 1,numer%ntor
      tdiv(i + 1) = i*htor
   END DO

   coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
   coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)

#ifdef PARALL
   IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
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
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
      DO ifa = 1,Mesh%Nextfaces

         ! Global numbering of this face
         Fi = ifa + Mesh%Nintfaces

         ! Element to which this face belongs
         iel = Mesh%extfaces(ifa,1)
         iel3 = (itor - 1)*N2d+iel

         ! Face local numbering
         ifl = Mesh%extfaces(ifa,2)

         ! Nodes in local numbering
         nod = refElPol%Face_nodes(ifl,:)

         ! Coordinates of the nodes of the face
         Xf = Mesh%X(Mesh%T(iel,nod),:)

         !****************************************************
         !                      Magnetic field
         !****************************************************
         ! Magnetic field of the nodes of the face
         indbt = colint(tensorSumInt(Mesh%T(iel,nod),(itor - 1)*(Np1Dtor - 1)*Mesh%Nnodes +&
                &Mesh%Nnodes*((/(i,i=1,Np1Dtor)/) - 1)))
         Bfl = phys%B(indbt,:)
         ! Magnetic field norm and direction at element nodes
         Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
         b_nod(:,1) = Bfl(:,1)/Bmod_nod
         b_nod(:,2) = Bfl(:,2)/Bmod_nod
         b_nod(:,3) = Bfl(:,3)/Bmod_nod

         ! Gauss points position
         xyf = matmul(refElPol%N1D,Xf)
         thetafg = matmul(refElTor%N1D,tel)

         ! Exterior normal & elemental surface
         xyd_g = matmul(refelPol%Nxi1D,Xf)
         xydNorm_g = sqrt(xyd_g(:,1)**2 + xyd_g(:,2)**2)
         dsurf = col(tensorProduct(refElPol%gauss_weights1D*xydNorm_g,refElTor%gauss_weights1D*0.5*htor))
         t_g(:,1) = xyd_g(:,1)/xydNorm_g
         t_g(:,2) = xyd_g(:,2)/xydNorm_g
         n_g = 0.
         DO i = 1,NGaussTor
            ind = (i - 1)*NGaussPol + (/(j,j=1,NGaussPol)/)
            n_g(ind,1) = t_g(:,2)
            n_g(ind,2) = -t_g(:,1)
         END DO

         ! Assembly indices
         dd = 1 + (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Nfdir)*Npfl)*Neq
         delta = dd + (N2D*Np2D+(Fi - 1)*Npfl)*Neq
         ind_uf = delta + (/(i,i=0,Npfl*Neq - 1)/)
         ind_ff = Np2d*Neq + (ifl - 1)*Npfl*Neq + (/(i,i=1,Npfl*Neq)/)
         ind_fe = colint(tensorSumInt((/(i,i=1,Neq)/),(refElTor%faceNodes3(ifl,:) - 1)*Neq))
         ind_fg = colint(tensorSumInt((/(i,i=1,3*Neq)/),(refElTor%faceNodes3(ifl,:) - 1)*3*Neq))
         ind_qf = (iel3 - 1)*Ndim*Neq*Npel + ind_fg
         qqf = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Npfl/)))

         ! Solution in this face
         uf = sol%u_tilde(ind_uf)

         ! Type of boundary
         fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
         IF (fl .eq. 0) CYCLE
#endif
         bc = phys%bcflags(fl)

         SELECT CASE (bc)

         CASE (bc_dirichlet)
            CALL set_dirichletstr_bc()
         CASE (bc_dirichlet_weak_form)
            CALL set_dirichletwf_bc()
         CASE (bc_dirichlet_and_Neumann)
            CALL set_diric_neum_bc()
!        CASE(bc_dirichlet_weak_form_oldvalues)
!           CALL set_dirichletwfoldval_bc()
!                                           CASE(bc_NeumannH)
!                                              CALL set_NeumannH()
         CASE (bc_Bohm)
            CALL set_Bohm_bc()
         CASE DEFAULT
            WRITE (6,*) "Error: wrong boundary type"
            STOP
         END SELECT

      END DO
   END DO


   if (utils%timing) then
      call cpu_time(timing%tpe1)
      call system_clock(timing%cke1,timing%clock_rate1)
      timing%runtbcd = timing%runtbcd + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
      timing%cputbcd = timing%cputbcd + timing%tpe1 - timing%tps1
   end if


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

            g = (igtor - 1)*NGaussPol + igpol

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

   !***************************************************************
   ! Dirichlet in weak form in some variables and Neumann in others
   !***************************************************************
   SUBROUTINE set_diric_neum_bc
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

            g = (igtor - 1)*NGaussPol + igpol

            ! Calculate the integration weight
            IF (switch%axisym) THEN
               dsurfg = dsurf(g)*xyf(igpol,1)
            ELSE
               dsurfg = dsurf(g)
            END IF
            NiNi = dsurfg*tensorProduct(Nfi(g,:),Nfi(g,:))
            Ni = Nfi(g,:)*dsurfg

            ! Assembly contribution
            CALL assembly_diric_neum_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fg,uex(g,:),NiNi,Ni,n_g(g,:))

         END DO
      END DO

   END SUBROUTINE set_diric_neum_bc

   !****************************
   ! Dirichlet in strong form
   !****************************
   SUBROUTINE set_dirichletstr_bc
      integer                   :: g,igpol,igtor
      real*8                    :: uex(refElTor%Ngl,neq),dsurfg
      real*8,pointer            :: Nfi(:,:)
      real*8                    :: Ni(Npfl)
      real*8                    :: b(Ngfl,3)

      ! Magnetic field norm and direction at Gauss points
      b = matmul(refElTor%sFTF,b_nod)

      ! Face shape functions
      Nfi => refElTor%sFTF

      ! Set the face solution
      CALL analytical_solution(xyf(:,1),xyf(:,2),thetafg,uex)

      ! Gradient solution at face gauss points
      qfg = matmul(Nfi,qqf)

      !*****************************
      ! Loop in face Gauss points
      !*****************************
      DO igtor = 1,NGaussTor
         DO igpol = 1,NGaussPol

            g = (igtor - 1)*NGaussPol + igpol

            ! Calculate the integration weight
            IF (switch%axisym) THEN
               dsurfg = dsurf(g)*xyf(igpol,1)
            ELSE
               dsurfg = dsurf(g)
            END IF

            ! Shaspe functions product
            Ni = Nfi(g,:)*dsurfg
            CALL assembly_dirichletstr_bc(iel3,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),b(g,:),n_g(g,:))

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
      real*8                    :: qfg(Ngfl,neq*3)
      real*8                    :: bn
      logical                   :: ntang
      real                      :: tau_stab(Neq,Neq)
      real*8                    :: diff_iso_fac(phys%neq,phys%neq,Ngfl),diff_ani_fac(phys%neq,phys%neq,Ngfl)
      real*8                    :: Bmod(Ngfl),b(Ngfl,3)

      ! Magnetic field norm and direction at Gauss points
      Bmod = matmul(refElTor%sFTF,Bmod_nod)
      b = matmul(refElTor%sFTF,b_nod)

      ! Face shape functions
      Nfi => refElTor%sFTF

      ! Gradient solution at face gauss points
      qfg = matmul(Nfi,qqf)

      ! Set the face solution
      ufg = matmul(Nfi,transpose(reshape(uf,[neq,Npfl])))

      ! Physical variables at Gauss points
      CALL cons2phys(ufg,upg)

      ! Compute diffusion at faces Gauss points
      CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac,Bmod)

      ! Compute magnetic field at Gauss points

      !*****************************
      ! Loop in face Gauss points
      !*****************************
      DO igtor = 1,NGaussTor
         DO igpol = 1,NGaussPol

            g = (igtor - 1)*NGaussPol + igpol

            ! Calculate the integration weight
            IF (switch%axisym) THEN
               dsurfg = dsurf(g)*xyf(igpol,1)
            ELSE
               dsurfg = dsurf(g)
            END IF

            bn = dot_product(b(g,:),n_g(g,:))

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
            inc = bn/norm2(b(g,1:2))

            setval = ufg(g,2)

            IF (abs(inc) .le. phys%bohmth) THEN
               sn = 1.
               delta = 0
               ntang = .FALSE.
            ELSE
               sn = sign(1.,inc)
               delta = 1
               IF ((abs(upg(g,2))) .le. SoundSpeed) THEN
                  setval = sn*SoundSpeed
                  delta = 1
               ELSE IF ((abs(upg(g,2))) .gt. SoundSpeed) THEN
                  delta = 0
                  setval = upg(g,2)
               END IF
            END IF

            ! Impose the velocity at the sound speed in case of subsonic
            ! velocity.
            NiNi = tensorProduct(Nfi(g,:),Nfi(g,:))*dsurfg
            Ni = Nfi(g,:)*dsurfg

            !****************** Stabilization part******************************
            IF (numer%stab > 1) THEN
               ! Compute tau in the Gauss points
               IF (numer%stab < 5) THEN
                  CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),iel,ifa,1.,xyf(g,:),tau_stab)
               ELSE
                  CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,:),n_g(g,:),xyf(g,:),1.,iel,tau_stab)
               ENDIF
               !****************** End stabilization part******************************
            ELSE
               tau_stab=0.
               DO i = 1,Neq
                  tau_stab(i,i) = numer%tau(i)
               END DO
            END IF

            ! Assembly Bohm contribution
            CALL assembly_bohm_bc(iel3,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),ufg(g,:),&
                 &b(g,:),n_g(g,:),tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
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
   integer                   :: ind_qf(refElPol%Nfacenodes*phys%neq*Mesh%Ndim)
   real*8                    :: qqf(refElPol%Nfacenodes,phys%neq*Mesh%Ndim),qfg(refElPol%Ngauss1d,phys%neq*Mesh%Ndim)
   real*8                    :: coefi,coefe
   real*8                    :: tau_save_el(refElPol%Ngauss1d,phys%neq),xy_g_save_el(refElPol%Ngauss1d,2)
   logical                   :: save_tau
   real*8,allocatable        :: tau_save(:,:)
   real*8,allocatable        :: xy_g_save(:,:)
   integer                   :: indtausave(refElPol%Ngauss1d)
   integer                   :: ind_asf(refElPol%Nfacenodes),ind_ash(refElPol%Nfacenodes)
   real*8                    :: Bfl(refElPol%Nfacenodes,3)
   real*8                    :: Bmod_nod(refElPol%Nfacenodes),b_nod(refElPol%Nfacenodes,3)

   save_tau = switch%saveTau
   Ndim = 2
   Npel = refElPol%Nnodes2D
   Npfl = refElPol%Nfacenodes
   Neq = phys%neq
   Ng1d = refElPol%Ngauss1d

   coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
   coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)

   ! Some indices
   ind_asf = (/(i,i=0,Neq*(Npfl - 1),Neq)/)
   ind_ash = (/(i,i=0,Neq*(Npfl - 1)*Ndim,Neq*Ndim)/)

   if (save_tau) then
      allocate (tau_save(Mesh%Nextfaces*refElPol%Ngauss1d,phys%neq))
      allocate (xy_g_save(Mesh%Nextfaces*refElPol%Ngauss1d,2))
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

      !****************************************************
      !                      Magnetic field
      !****************************************************
      ! Magnetic field of the nodes of the face
      Bfl = phys%B(Mesh%T(iel,nod),:)
      ! Magnetic field norm and direction at element nodes
      Bmod_nod = sqrt(Bfl(:,1)**2 + Bfl(:,2)**2 + Bfl(:,3)**2)
      b_nod(:,1) = Bfl(:,1)/Bmod_nod
      b_nod(:,2) = Bfl(:,2)/Bmod_nod
      b_nod(:,3) = Bfl(:,3)/Bmod_nod

      ! Assembly indices
      ind_uf = (Fi - 1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
      ind_ff = (ifl - 1)*Neq*Npfl + (/(i,i=1,Neq*Npfl)/)
      ind_fe = reshape(tensorSumInt((/(i,i=1,neq)/),neq*(nod - 1)),(/neq*Npfl/))
      ind_fG = reshape(tensorSumInt((/(i,i=1,neq*Ndim)/),neq*Ndim*(refElPol%face_nodes(ifl,:) - 1)),(/neq*Ndim*Npfl/))
      ind_qf = (iel - 1)*Ndim*Neq*Npel + ind_fG
      qqf = transpose(reshape(sol%q(ind_qf),(/Ndim*Neq,Npfl/)))

      ! Solution in this face
      uf = sol%u_tilde(ind_uf)

      ! Type of boundary
      fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
      IF (fl .eq. 0) CYCLE
#endif
      bc = phys%bcflags(fl)

      SELECT CASE (bc)

      CASE (bc_dirichlet)
         CALL set_dirichletstr_bc()
      CASE (bc_dirichlet_weak_form)
         CALL set_dirichletwf_bc()
      CASE (bc_dirichlet_and_Neumann)
         CALL set_diric_neum_bc()
!                                        CASE(bc_NeumannH)
!                                           CALL set_NeumannH()
      CASE (bc_Bohm)

         CALL set_Bohm_bc(tau_save_el,xy_g_save_el)
         if (save_tau) then
            indtausave = (ifa - 1)*refElPol%Ngauss1d+(/(i,i=1,refElPol%Ngauss1d)/)
            tau_save(indtausave,:) = tau_save_el
            xy_g_save(indtausave,:) = xy_g_save_el
         endif

      CASE DEFAULT
         WRITE (6,*) "Error: wrong boundary type"
         STOP
      END SELECT

   END DO

   if (save_tau) then
      write (6,*) "Saving tau in the boundary faces"
      call saveMatrix(tau_save,'tau_save_bound')
      call saveMatrix(xy_g_save,'xy_g_save_bound')
      deallocate (tau_save,xy_g_save)
      write (6,*) "Done saving tau!"
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
      xyf = matmul(refElPol%N1D,Xf)
      CALL analytical_solution(xyf(:,1),xyf(:,2),uex)

      ! Shape function derivatives at Gauss points
      xyDer = matmul(refElPol%Nxi1D,Xf)

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

   !***************************************************************
   ! Dirichlet in weak form in some variables and Neumann in others
   !***************************************************************
   SUBROUTINE set_diric_neum_bc
      integer                   :: g
      real*8                    :: dline,xyDerNorm_g
      real*8                    :: t_g(Ndim),n_g(Ndim)
      real*8                    :: uex(Ng1d,neq)
      real*8                    :: xyf(Ng1d,Ndim),xyDer(Ng1d,Ndim)
      real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)

      ! Set the face solution
      xyf = matmul(refElPol%N1D,Xf)
      CALL analytical_solution(xyf(:,1),xyf(:,2),uex)

      ! Shape function derivatives at Gauss points
      xyDer = matmul(refElPol%Nxi1D,Xf)

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
         n_g = [t_g(2),-t_g(1)]

         NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
         Ni = refElPol%N1D(g,:)*dline
         CALL assembly_diric_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,uex(g,:),NiNi,Ni,n_g)
      END DO

   END SUBROUTINE set_diric_neum_bc

   !****************************
   ! Dirichlet in strong form
   !****************************
   SUBROUTINE set_dirichletstr_bc
      integer                   :: g,i,j,idm
      real*8                    :: dline,xyDerNorm_g
      real*8                    :: xyf(Ng1d,Ndim),xyDer(Ng1d,Ndim)
      real*8                    :: uex(Ng1d,neq)
      real*8                    :: t_g(Ndim),n_g(Ndim),fh(Npfl*neq)
      real*8                    :: Ni(Npfl)
      real*8                    :: b(Ng1d,3)

      ! Magnetic field norm and direction at Gauss points
      b = matmul(refElPol%N1d,b_nod)

      ! Set the face solution
      xyf = matmul(refElPol%N1D,Xf)
      CALL analytical_solution(xyf(:,1),xyf(:,2),uex)

      ! Shape function derivatives at Gauss points
      xyDer = matmul(refElPol%Nxi1D,Xf)

      ! Gradient solution at face gauss points
      qfg = matmul(refElPol%N1D,qqf)

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
         n_g = [t_g(2),-t_g(1)]

         ! Shaspe functions product
         Ni = refElPol%N1D(g,:)*dline
         CALL assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qfg(g,:),uex(g,:),b(g,1:2),n_g)

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
      real*8                    :: ufg(Ng1d,Neq),upg(Ng1d,phys%npv)
      real*8                    :: inc,sn,delta
      real*8                    :: NiNi(Npfl,Npfl),Ni(Npfl)
      real*8                    :: soundSpeed
      real*8                    :: tau_ngamma(2)
!                                                        real*8                    :: b_f(Ng1d,Ndim),divb_f(Ng1d),drift_f(Ng1d,Ndim)
!#ifdef MOVINGEQUILIBRIUM
!       real*8                    :: b_nod(Mesh%Nnodesperelem,Ndim),Bmod(Mesh%Nnodesperelem)
!#endif
      real*8                    :: qfg(Ng1d,neq*Ndim)
      real*8                    :: bn
!      logical                   :: ntang
      real*8,intent(out)        :: tau_save_el(:,:),xy_g_save_el(:,:)
      real                      :: tau_stab(Neq,Neq)
      real*8                    :: diff_iso_fac(phys%neq,phys%neq,Ng1d),diff_ani_fac(phys%neq,phys%neq,Ng1d)
      real*8                    :: Bmod(Ng1d),b(Ng1d,3)

      ! Magnetic field norm and direction at Gauss points
      Bmod = matmul(refElPol%N1d,Bmod_nod)
      b = matmul(refElPol%N1d,b_nod)

      ! Gradient solution at face gauss points
      qfg = matmul(refElPol%N1D,qqf)

      ! Gauss points position
      xyf = matmul(refElPol%N1D,Xf)

      ! Compute diffusion at faces Gauss points
      CALL setLocalDiff(xyf,diff_iso_fac,diff_ani_fac,Bmod)

      ! Shape function derivatives at Gauss points
      xyDer = matmul(refElPol%Nxi1D,Xf)

      ! Set the face solution
      ufg = matmul(refElPol%N1D,transpose(reshape(uf,[neq,Npfl])))

      ! Shape function derivatives at Gauss points
      xyDer = matmul(refElPol%Nxi1D,Xf)

      ! Physical variables at Gauss points
      CALL cons2phys(ufg,upg)

! VEDERE PERCHÈ CALCOLO DUE VOLTE LE VARIABILI FISICHE NEI PUNTI DI GAUSS!!!!!!!!!!!!!!!!!!
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
         n_g = [t_g(2),-t_g(1)]
         bn = dot_product(b(g,1:2),n_g)
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
         inc = bn/norm2(b(g,1:2))

#ifdef NGAMMA
         setval = ufg(g,2)

         IF (abs(inc) .le. phys%bohmth) THEN
            sn = 1.
            delta = 0 ! TODO try to module delta in the WEST case
            ntang = .FALSE.
         ELSE
            sn = sign(1.,bn)
            delta = 1
            IF ((abs(upg(g,2))) .le. SoundSpeed) THEN
               setval = sn*SoundSpeed
                 delta = 1
            ELSE IF ((abs(upg(g,2))) .gt. SoundSpeed) THEN
                 delta = 0
               setval = upg(g,2)
            END IF
         END IF
#endif
         ! Shape functions
         NiNi = tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*dline
         Ni = refElPol%N1D(g,:)*dline
         !****************** Stabilization part******************************
         IF (numer%stab > 1) THEN
            ! Compute tau in the Gauss points
            IF (numer%stab < 5) THEN
               CALL computeTauGaussPoints(upg(g,:),ufg(g,:),b(g,1:2),n_g,iel,ifa,1.,xyf(g,:),tau_stab)
            ELSE
               CALL computeTauGaussPoints_matrix(upg(g,:),ufg(g,:),b(g,1:2),n_g,xyf(g,:),1.,iel,tau_stab)
            ENDIF

            if (save_tau) then
               DO i = 1,Neq
                  tau_save_el(g,:) = tau_stab(i,i)
               END DO
               xy_g_save_el(g,:) = xyf(g,:)
            endif
            !****************** End stabilization part******************************
         ELSE 
            tau_stab=0.
            DO i = 1,Neq
               tau_stab(i,i) = numer%tau(i)
            END DO
         END IF

         ! Assembly Bohm contribution
!         CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
!              &ufg(g,:),b(g,1:2),n_g,tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g))

         CALL assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg(g,:),&
              &ufg(g,:),b(g,1:2),n_g,tau_stab,setval,delta,diff_iso_fac(:,:,g),diff_ani_fac(:,:,g),ntang)
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
      integer*4    :: ind(Npfl),i

      kmult = col(tensorProduct(uexg(:),Ni))
      DO i = 1,Neq
         ind = i + ind_asf
         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
      END DO

   END SUBROUTINE assembly_dirichletwf_bc

   !*********************************
   ! Assembly Dirichlet and Neumann
   !*********************************
   SUBROUTINE assembly_diric_neum_bc(iel,ind_asf,ind_ash,ind_ff,ind_fg,uexg,NiNi,Ni,n)
      integer*4    :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fg(:)
      real*8       :: NiNi(:,:),Ni(:),uexg(:),n(:)
      real*8       :: kmult(Neq*Npfl)
      integer*4    :: i,j,k,idm,ndir
      integer*4    :: ind(Npfl),indi(Npfl),indj(Npfl)

!      ndir = 2
!      ! Dirichlet part
!      kmult= col(tensorProduct(uexg(:),Ni))
!      DO i=1,Neq
!         ind = i+ind_asf
!         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
!         if (i.le.ndir) then
!            elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
!         else
!            elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + numer%tau(i)*NiNi
!         endif
!      END DO
! 
!      ! Neumann part
!                           DO k = ndir+1,Neq
!                              DO idm = 1,Ndim
!                                 indi = ind_asf+k
!            indj = ind_ash+idm+(k-1)*Ndim
!            elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)+NiNi*n(idm)
!                              END DO
!                                                END DO

      ! Dirichlet values
      kmult = col(tensorProduct(uexg(:),Ni))

      ! Dirichlet weak form for the first equation
      i = 1
      ind = i + ind_asf
      elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)

      ! Dirichlet weak form for the second equation
      i = 2
      ind = i + ind_asf
      elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)

      if (switch%dirivortcore) then
         ! Dirichlet weak form for the vorticity equation
         i = 3
         ind = i + ind_asf
         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
         elMat%fh(ind_ff(ind),iel) = elMat%fh(ind_ff(ind),iel) - numer%tau(i)*kmult(ind)
      else
         ! Neumann homogeneous for the vorticity
         i = 3
         ind = i + ind_asf
         elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
         elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + numer%tau(i)*NiNi
         DO idm = 1,Ndim
            indi = ind_asf + i
            indj = ind_ash + idm + (i - 1)*Ndim
            elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*n(idm)
         END DO
      endif

      ! Neumann homogeneous for the potential
      i = 4
      ind = i + ind_asf
      elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - numer%tau(i)*NiNi
      elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + numer%tau(i)*NiNi
      DO idm = 1,Ndim
         indi = ind_asf + i
         indj = ind_ash + idm + (i - 1)*Ndim
         elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*n(idm)
      END DO
      j = 1
      DO idm = 1,Ndim
         indi = ind_asf + i
         indj = ind_ash + idm + (j - 1)*Ndim
         elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*n(idm)/uexg(1)*phys%Mref
      END DO

   END SUBROUTINE assembly_diric_neum_bc

   !*********************************
   ! Assembly Dirichlet in strong form
   !*********************************
   SUBROUTINE assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Ni,qg,uexg,bg,ng)
      integer*4        :: iel,ind_asf(:),ind_ash(:),ind_fe(:),ind_fg(:)
      real*8           :: Ni(:),uexg(:),bg(:),ng(:)
      real*8           :: qg(:)
      real*8           :: bn,An(Neq,Neq)
      integer          :: i,j,idm
      real*8           :: kmultstab(Neq*Npfl),kmultconv(Neq*Npfl),kmultfseq(Ndim*Neq*Npfl)
      real*8           :: uexn(Ndim*neq)
      integer*4        :: ind(Npfl),ind_1f(Npfl)
      real*8           :: Qpr(Ndim,Neq)
#ifdef TEMPERATURE
      real*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
      real*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
#endif
      bn = dot_product(bg,ng)

      CALL jacobianMatricesFace(uexg,bn,An)

      ! Contribution of the current integration point to the elemental matrix
      kmultstab = col(tensorProduct(uexg,Ni))
      kmultconv = col(tensorProduct(matmul(An,uexg),Ni))
      uexn = reshape(tensorProduct(ng,uexg),(/Ndim*neq/)) ! Ndim*neq x 1
      kmultfseq = col(tensorProduct(uexn,Ni))  ! Ndim*neq*npf

      DO i = 1,Neq
         ind = i + ind_asf
         elMat%Aul_dir(ind_fe(ind),iel) = elMat%Aul_dir(ind_fe(ind),iel) - numer%tau(i)*kmultstab(ind) + kmultconv(ind)
         DO idm = 1,Ndim
            ind_1f = ind_ash + idm + (i - 1)*Ndim
            elMat%Aql_dir(ind_fG(ind_1f),iel) = elMat%Aql_dir(ind_fG(ind_1f),iel) - kmultfseq(ind_1f)
         END DO
      END DO

      ! Compute Q^T^(k-1)
      Qpr = reshape(qg,(/Ndim,Neq/))

#ifdef TEMPERATURE
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

      DO i = 3,4
         DO j = 1,4
            IF (i == 3) THEN
               elMat%Aul_dir(ind_fe(i + ind_asf),iel) = elMat%Aul_dir(ind_fe(i + ind_asf),iel) - &
                             &coefi*Ni*bn*((gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),bg)))*uexg(j))
            ELSE
               elMat%Aul_dir(ind_fe(i + ind_asf),iel) = elMat%Aul_dir(ind_fe(i + ind_asf),iel) - &
                             &coefe*Ni*bn*((gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),bg)))*uexg(j))
            END IF
         END DO
      END DO
#endif
   END SUBROUTINE assembly_dirichletstr_bc

   !*********************************
   ! Assembly Bohm
   !*********************************
   SUBROUTINE assembly_bohm_bc(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,ufg,bg,ng,tau,setval,delta,diffiso,diffani,ntang)
      integer*4        :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8           :: NiNi(:,:),Ni(:),ufg(:),bg(:),ng(:),tau(:,:),setval,delta
      real*8           :: diffiso(:,:),diffani(:,:)
      logical          :: ntang
      real*8           :: qfg(:)
      real*8           :: bn,Abohm(Neq,Neq)
      integer          :: i,j,k,idm,Neqstab,Neqgrad
      integer*4        :: ind(Npfl),indi(Npfl),indj(Npfl),indk(Npfl)
      real*8           :: Qpr(Ndim,Neq),kcoeff
#ifdef TEMPERATURE
      real*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
      real*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
#endif

      Neqstab = Neq
      Neqgrad = Neq
#ifdef VORTICITY
      Neqgrad = 2
      IF (switch%fixdPotLim) THEN
         Neqstab = 3
      ENDIF
      IF (.not. (ntang)) THEN
         tau = 1.
      ENDIF
#endif

      bn = dot_product(bg,ng)

      ! Compute Q^T^(k-1)
      Qpr = reshape(qfg,(/Ndim,Neq/))

      ! Case diagonal matrix stabilization
      IF (numer%stab < 5) THEN
         DO i = 1,Neqstab
            ind = i + ind_asf
            elMat%All(ind_ff(ind),ind_ff(ind),iel) = elMat%All(ind_ff(ind),ind_ff(ind),iel) - tau(i,i)*NiNi
            IF (i == 2) THEN
               elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) + tau(i,i)*(delta*setval)*NiNi
               elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*(1 - delta)*NiNi
            ELSE
               elMat%Alu(ind_ff(ind),ind_fe(ind),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind),iel) + tau(i,i)*NiNi
            END IF
         END DO
      ELSE
         ! Case full matrix stabilization
         DO i = 1,Neqstab
            indi = i + ind_asf
            DO j = 1,Neqstab
               indj = j + ind_asf
               elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - tau(i,j)*NiNi
               IF (j == 1) THEN
                  elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) +&
                                                                      &tau(i,j)*NiNi + tau(i,j + 1)*(delta*setval)*NiNi
               ELSE IF (j == 2) THEN
                elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + tau(i,j)*(1 - delta)*NiNi
               ELSE
                  elMat%Alu(ind_ff(indi),ind_fe(indj),iel) = elMat%Alu(ind_ff(indi),ind_fe(indj),iel) + tau(i,j)*NiNi
               ENDIF
            END DO
         END DO
      ENDIF

#ifdef TEMPERATURE
      IF (ntang) THEN

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
         DO i = 1,2
            indi = ind_asf + i
            IF (i == 1) THEN
               DO j = 1,Neq
                  indj = ind_asf + j
               elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                  elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) -&
                                                           &coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),bg)))*NiNi*bn
                  DO k = 1,Ndim
                         elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                   &coefi*Alphai*Vveci(j)*bg(k)*NiNi*bn
                  END DO
               END DO
                    elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefi*Alphai*( dot_product (matmul(transpose(Taui),bg),ufg)  )*Ni*bn
            ELSE
               DO j = 1,Neq
                  indj = ind_asf + j
               elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                  elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) -&
                                                           &coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),bg)))*NiNi*bn
                  DO k = 1,Ndim
                         elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                   &coefe*Alphae*Vvece(j)*bg(k)*NiNi*bn
                  END DO
               END DO
                   elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefe*Alphae*( dot_product (matmul(transpose(Taue),bg),ufg)  )*Ni*bn 
            ENDIF
         END DO
      END IF ! tangency
#endif

!         ! Perpendicular diffusion
!                                                   DO k = 1,Neq
!                                                      DO idm = 1,Ndim
!                                                         indi = ind_asf+k
!               indj = ind_ash+idm+(k-1)*Ndim
!                                                         if (ntang) then
!                  ! Non-tangent case
!                                                            elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
!                                                               &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) )
!                                                         else
!                  ! Tangent case
!!                                                             elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*(ng(idm)-bn*bg(idm))
!                   elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-NiNi*ng(idm)
!                                                         endif
!                                                      END DO
!                                                                        END DO

      ! Perpendicular diffusion
      IF (ntang) THEN
         DO k = 1,Neqgrad
            DO idm = 1,Ndim
               indi = ind_asf + k
               indj = ind_ash + idm + (k - 1)*Ndim
               ! Non-tangent case
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                      &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) ) 
            END DO
         END DO
      ELSE
         DO k = 1,Neq
            DO idm = 1,Ndim
               indi = ind_asf + k
               indj = ind_ash + idm + (k - 1)*Ndim
               ! Tangent case
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
            END DO
         END DO
      ENDIF

#ifdef VORTICITY
      IF (ntang) THEN
         k = 3
         ! Vorticity equation         
         if (switch%dirivortlim) then
            ! Dirichlet weak form for the vorticity equation: set vorticity to 0!!!
            indi = k + ind_asf
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - numer%tau(k)*NiNi
            !elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(k)*kmult(indi)
         else
            DO idm = 1,Ndim
               indi = ind_asf + k
               indj = ind_ash + idm + (k - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
            END DO
         end if
      END IF

      ! Potential equation
      IF (ntang) THEN
         IF (switch%fixdPotLim) THEN
            !*****************************************
            ! Fixed potential on the limiter
            !*****************************************
            i = 4
            j = 4
            indi = ind_asf + i
            indj = ind_asf + j
            elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - numer%tau(4)*NiNi
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(4)*Ni*phys%Mref*phys%Potfloat
         ELSE
            !*****************************************
            ! Bohm condition for the potential (Robin)
            !*****************************************
            ! Diagonal part (\Gradpar Phi)
            i = 4
            indi = ind_asf + i
            DO idm = 1,Ndim
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*bg(idm)
            END DO

            ! Non diagonal part --> non linear part (\Gradpar n / n)
            i = 4
            j = 1
            indi = ind_asf + i
            indj = ind_asf + j
            DO idm = 1,Ndim
               indk = ind_ash + idm + (j - 1)*Ndim
               kcoeff = phys%Mref*bg(idm)
               ! Linear part
               elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kcoeff*NiNi/ufg(1)
               ! Non linear correction
               elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) + &
                                                          &kcoeff*Qpr(idm,1)*NiNi/ufg(1)**2
               elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) + kcoeff*Qpr(idm,1)*Ni/ufg(1)
            END DO

            ! Non diagonal part --> linear part (\Gammma * \Lambda)
            i = 4
            j = 2
            indi = ind_asf + i
            indj = ind_asf + j
            elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) + &
                                                       &phys%etapar/phys%c2*NiNi*phys%Mref*phys%Potfloat

            ! Diagonal part --> non linear part (\Gammma * \phi)
            ! \phi^(k+1)*\Gamma^k
            i = 4
            j = 4
            indi = ind_asf + i
            indj = ind_asf + j
            elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - &
                                                       &phys%etapar/phys%c2*NiNi*ufg(2)
            ! \phi^k*\Gamma^(k+1)
            i = 4
            j = 2
            indi = ind_asf + i
            indj = ind_asf + j
           elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - &
                                                       &phys%etapar/phys%c2*NiNi*ufg(4)
            ! \phi^k*\Gamma^k
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - phys%etapar/phys%c2*Ni*ufg(4)*ufg(2)

         END IF
      END IF

#endif

   END SUBROUTINE assembly_bohm_bc












   !*********************************
   ! Assembly Bohm
   !*********************************
   SUBROUTINE assembly_bohm_bc_new(iel,ind_asf,ind_ash,ind_ff,ind_fe,ind_fg,NiNi,Ni,qfg,&
                                   &ufg,bg,ng,tau,setval,delta,diffiso,diffani,ntang)
      integer*4        :: iel,ind_asf(:),ind_ash(:),ind_ff(:),ind_fe(:),ind_fg(:)
      real*8           :: NiNi(:,:),Ni(:),ufg(:),bg(:),ng(:),tau(:,:),setval,delta
      real*8           :: diffiso(:,:),diffani(:,:)
      logical          :: ntang
      real*8           :: qfg(:)
      real*8           :: bn,Abohm(Neq,Neq)
      integer          :: i,j,k,idm,Neqstab,Neqgrad
      integer*4        :: ind(Npfl),indi(Npfl),indj(Npfl),indk(Npfl)
      real*8           :: Qpr(Ndim,Neq),kcoeff
#ifdef TEMPERATURE
      real*8           :: Vveci(Neq),Alphai,taui(Ndim,Neq),dV_dUi(Neq,Neq),gmi,dAlpha_dUi(Neq)
      real*8           :: Vvece(Neq),Alphae,taue(Ndim,Neq),dV_dUe(Neq,Neq),gme,dAlpha_dUe(Neq)
#endif

      
      bn = dot_product(bg,ng)
      inc = bn/norm2(bg(1:2))

      ! Compute Q^T^(k-1)
      Qpr = reshape(qfg,(/Ndim,Neq/))


      !*****************************************************
      ! First equation
      !*****************************************************
      i = 1
      indi = i + ind_asf
      ! Stabilization part
      elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
      elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
      ! Gradient part
      DO idm = 1,Ndim
         indj = ind_ash + idm + (i - 1)*Ndim
         elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) ) 
      END DO



      !*****************************************************
      ! Second equation
      !*****************************************************
      i = 2
      indi = i + ind_asf
      if (numer%bohmtypebc.eq.1) then
         ! Dirichlet type always if non-tangent 
         IF (ntang) THEN
            ! >>>>>>>> Non tangent case (delta=1) <<<<<<<<<<
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*Ni*setval ! TODO: verify this
!            elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) + numer%tau(i)*setval*NiNi
         ELSE
            ! >>>>>>>> Tangent case (delta=0) <<<<<<<<<<
            ! Stabilization part
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
            ! Gradient part
            DO idm = 1,Ndim
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                      &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) ) 
            END DO
         END IF   
      elseif (numer%bohmtypebc.eq.2) then
         ! Dirichlet if u<=soundspeed, Neumann if u>soundspeed if non-tangent 
         IF (ntang) THEN
            ! >>>>>>>> Non tangent case  <<<<<<<<<<
            if (delta==1) then
            ! Velocity is subsonic--> Dirichlet
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*Ni*setval ! TODO: verify this
!            elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) = elMat%Alu(ind_ff(ind),ind_fe(ind - 1),iel) + numer%tau(i)*setval*NiNi
            else
            ! Velocity is supersonic--> Neumann (delta=0)
            ! Stabilization part
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
            ! Gradient part
            DO idm = 1,Ndim
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                      &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) ) 
            END DO
         ELSE
            ! >>>>>>>> Tangent case (delta=0) <<<<<<<<<<
            ! Stabilization part
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
            ! Gradient part
            DO idm = 1,Ndim
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                      &NiNi*(ng(idm)*diffiso(i,i)-bn*bg(idm)*diffani(i,i) ) 
            END DO
         END IF  
      else 
         write(6,*) "Wrong type in Bohm bc"
         stop
      endif

#ifdef TEMPERATURE
      !*****************************************************
      ! Third and fourth equation - TEMPERATURE MODEL
      !*****************************************************
      IF (ntang) THEN

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
         DO i = 1,2
            indi = ind_asf + i
            IF (i == 1) THEN
               DO j = 1,Neq
                  indj = ind_asf + j
               elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                  elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) -&
                                                           &coefi*(gmi*dAlpha_dUi(j) + Alphai*(dot_product(Taui(:,j),bg)))*NiNi*bn
                  DO k = 1,Ndim
                         elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                   &coefi*Alphai*Vveci(j)*bg(k)*NiNi*bn
                  END DO
               END DO
                    elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefi*Alphai*( dot_product (matmul(transpose(Taui),bg),ufg)  )*Ni*bn
            ELSE
               DO j = 1,Neq
                  indj = ind_asf + j
               elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) + Abohm(i,j)*NiNi*bn
                  elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) = elMat%All(ind_ff(indi + 2),ind_ff(indj),iel) -&
                                                           &coefe*(gme*dAlpha_dUe(j) + Alphae*(dot_product(Taue(:,j),bg)))*NiNi*bn
                  DO k = 1,Ndim
                         elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) = elMat%Alq(ind_ff(indi+2),ind_fG(k+(j-1)*Ndim+ind_ash),iel) -&
                                                                                   &coefe*Alphae*Vvece(j)*bg(k)*NiNi*bn
                  END DO
               END DO
                   elMat%fh(ind_ff(indi+2),iel) = elMat%fh(ind_ff(indi+2),iel) - coefe*Alphae*( dot_product (matmul(transpose(Taue),bg),ufg)  )*Ni*bn 
            ENDIF
         END DO

      
         DO i=3,4
            ! Stabilization part--> use tau for parallel diffusion
            indi = i + ind_asf
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
         END DO
      ELSE
         DO i=3,4
            ! Stabilization part--> use tau for perpendicular diffusion
            indi = i + ind_asf
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(1,1)*NiNi ! TODO: I assume the diffusion is the same in all eqs.
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(1,1)*NiNi ! TODO: I assume the diffusion is the same in all eqs.
         END DO

      END IF ! tangency

      ! Perpendicular diffusion for energies
      DO k = 3,4
         DO idm = 1,Ndim
            indi = ind_asf + k
            indj = ind_ash + idm + (k - 1)*Ndim
            ! Non-tangent case
            elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
                                                   &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) ) 
         END DO
      END DO
#endif


!      ! Perpendicular diffusion
!      IF (ntang) THEN
!         DO k = 1,Neqgrad
!            DO idm = 1,Ndim
!               indi = ind_asf + k
!               indj = ind_ash + idm + (k - 1)*Ndim
!               ! Non-tangent case
!               elMat%Alq(ind_ff(indi),ind_fG(indj),iel)=elMat%Alq(ind_ff(indi),ind_fG(indj),iel)-&
!                                                      &NiNi*(ng(idm)*diffiso(k,k)-bn*bg(idm)*diffani(k,k) ) 
!            END DO
!         END DO
!      ELSE
!         DO k = 1,Neq
!            DO idm = 1,Ndim
!               indi = ind_asf + k
!               indj = ind_ash + idm + (k - 1)*Ndim
!               ! Tangent case
!               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
!            END DO
!         END DO
!      ENDIF

#ifdef VORTICITY

      !*****************************************************
      ! Third equation - VORTICITY MODEL
      !*****************************************************
      IF (ntang) THEN
         i = 3
         ! Vorticity equation         
         if (switch%dirivortlim) then
            ! Dirichlet weak form for the vorticity equation: set vorticity to 0!!!
            indi = i + ind_asf
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
            !elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(k)*kmult(indi)
         else
            ! Neumann for the vorticity equation
            ! Stabilization part
            indi = i + ind_asf
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - tau(i,i)*NiNi
            elMat%Alu(ind_ff(indi),ind_fe(indi),iel) = elMat%Alu(ind_ff(indi),ind_fe(indi),iel) + tau(i,i)*NiNi
            DO idm = 1,Ndim
               indi = ind_asf + i
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) - NiNi*ng(idm)
            END DO
         end if
      END IF

      ! Potential equation
      IF (ntang) THEN
         IF (switch%fixdPotLim) THEN
            !*****************************************
            ! Fixed potential on the limiter
            !*****************************************
            i = 4
            indi = ind_asf + i
            elMat%All(ind_ff(indi),ind_ff(indi),iel) = elMat%All(ind_ff(indi),ind_ff(indi),iel) - numer%tau(i)*NiNi
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - numer%tau(i)*Ni*phys%Mref*phys%Potfloat
         ELSE
            !*****************************************
            ! Bohm condition for the potential (Robin)
            !*****************************************
            ! Diagonal part (\Gradpar Phi)
            i = 4
            indi = ind_asf + i
            DO idm = 1,Ndim
               indj = ind_ash + idm + (i - 1)*Ndim
               elMat%Alq(ind_ff(indi),ind_fG(indj),iel) = elMat%Alq(ind_ff(indi),ind_fG(indj),iel) + NiNi*bg(idm)
            END DO

            ! Non diagonal part --> non linear part (\Gradpar n / n)
            i = 4
            j = 1
            indi = ind_asf + i
            indj = ind_asf + j
            DO idm = 1,Ndim
               indk = ind_ash + idm + (j - 1)*Ndim
               kcoeff = phys%Mref*bg(idm)
               ! Linear part
               elMat%Alq(ind_ff(indi),ind_fG(indk),iel) = elMat%Alq(ind_ff(indi),ind_fG(indk),iel) - kcoeff*NiNi/ufg(1)
               ! Non linear correction
               elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) + &
                                                          &kcoeff*Qpr(idm,1)*NiNi/ufg(1)**2
               elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) + kcoeff*Qpr(idm,1)*Ni/ufg(1)
            END DO

            ! Non diagonal part --> linear part (\Gammma * \Lambda)
            i = 4
            j = 2
            indi = ind_asf + i
            indj = ind_asf + j
            elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) + &
                                                       &phys%etapar/phys%c2*NiNi*phys%Mref*phys%Potfloat

            ! Diagonal part --> non linear part (\Gammma * \phi)
            ! \phi^(k+1)*\Gamma^k
            i = 4
            j = 4
            indi = ind_asf + i
            indj = ind_asf + j
            elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - &
                                                       &phys%etapar/phys%c2*NiNi*ufg(2)
            ! \phi^k*\Gamma^(k+1)
            i = 4
            j = 2
            indi = ind_asf + i
            indj = ind_asf + j
           elMat%All(ind_ff(indi),ind_ff(indj),iel) = elMat%All(ind_ff(indi),ind_ff(indj),iel) - &
                                                       &phys%etapar/phys%c2*NiNi*ufg(4)
            ! \phi^k*\Gamma^k
            elMat%fh(ind_ff(indi),iel) = elMat%fh(ind_ff(indi),iel) - phys%etapar/phys%c2*Ni*ufg(4)*ufg(2)

         END IF
      END IF

#endif

   END SUBROUTINE assembly_bohm_bc_new
END SUBROUTINE HDG_BC

