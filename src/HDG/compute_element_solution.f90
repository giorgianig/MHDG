!*****************************************
! project: MHDG
! file: compute_solution.f90
! date: 24/01/2017
! Solve the global system, compute face
! and element solution
!*****************************************

SUBROUTINE compute_element_solution
  USE globals
  USE LinearAlgebra
  USE printUtils
  USE solve_pastix
  USE MPI_OMP
  IMPLICIT NONE

  !**********************************
  ! Retrieve the element solution
  !**********************************
#ifdef TOR3D
  !***********************************************************************
  !                            VERSION 3D TOROIDAL
  !***********************************************************************

  integer               :: itor, iel, iel3d, ifa, i
  integer               :: Ndim, Neq, N2D, Np2D, Np, Nfl, Ntorloc, Np1Dpol, Np1Dtor, Nf, Nfdir
  integer               :: ind_ue(refElTor%Nnodes3D*phys%neq), ind_uf(refElTor%Nft*phys%neq)
  integer*4             :: ind_ug(refElTor%Nnodes3D*phys%Neq*3)
  integer               :: dd, delta(refElPol%Nfaces + 2)
  integer               :: ind_dim(refElPol%Nfaces + 2), ind_sta(refElPol%Nfaces + 2), aux
  integer*4             :: Fe(refElPol%Nfaces)

  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1, timing%clock_rate1)
  end if


  Ndim = 3                                               ! Number of dimensions
  Neq = phys%neq                                        ! Number of equations
  N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
  Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
  Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
  Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
  Np = refElTor%Nnodes3D                               ! Number of nodes for each 3D element
  Nfl = refElTor%Nfl                                    ! Number of nodes in the lateral faces
  Nf = refElPol%Nfaces
  Nfdir = Mesh%Ndir

  ind_uf = 0
  ind_dim(1) = Np2D*Neq
  ind_dim(2:refElPol%Nfaces + 1) = Nfl*Neq
  ind_dim(refElPol%Nfaces + 2) = Np2D*Neq
  ind_sta(1) = 1
  aux = 0
  DO i = 2, refElPol%Nfaces + 2
    ind_sta(i) = 1 + ind_dim(i - 1) + aux
    aux = aux + ind_dim(i - 1)
  END DO

#ifdef PARALL
  IF (MPIvar%ntor .gt. 1) THEN
    ntorloc = numer%ntor/MPIvar%ntor + 1
  ELSE
    ntorloc = numer%ntor
  ENDIF
#else
  ntorloc = numer%ntor
#endif

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(itor,iel,iel3d,ind_ue,ind_uf,Fe,dd,delta,ifa)
  !$OMP DO
  DO itor = 1, ntorloc
    DO iel = 1, N2d
      iel3d = (itor - 1)*N2d+iel  !3d numbering of the element
      ind_ue = (iel3d-1)*Np*Neq + (/(i, i=1, Np*Neq)/)

      ! Index for the face solution
      Fe = Mesh%F(iel, :)
      dd = 1 + (itor - 1)*Neq*(N2d*Np2D+(Mesh%Nfaces - Nfdir)*Nfl)
      delta(1) = dd + (iel - 1)*Np2D*Neq
      delta(2:refElPol%Nfaces + 1) = dd + N2d*Np2D*Neq + (Fe - 1)*Nfl*Neq
      delta(refElPol%Nfaces + 2) = dd + N2d*Np2D*Neq + (Mesh%Nfaces - Nfdir)*Nfl*Neq + (iel - 1)*Np2D*Neq
      if (MPIvar%ntor .eq. 1) then
        if (itor == numer%ntor) then
          delta(refElPol%Nfaces + 2) = 1 + (iel - 1)*Np2D*Neq
        end if
      endif

      do ifa = 1, Nf + 2
        do i = 0, ind_dim(ifa) - 1
          ind_uf(i + ind_sta(ifa)) = delta(ifa) + i
        end do
      end do

      ! elemental solutions
      sol%u(ind_ue) = matmul(elMat%UU(:, :, iel3d), sol%u_tilde(ind_uf)) + elMat%U0(:, iel3d)

      ! Index for the element gradient
      ind_ug = (iel3d-1)*neq*Np*Ndim + (/(i, i=1, Neq*Np*Ndim)/)

      ! elemental solutions
      sol%q(ind_ug) = matmul(elMat%LL(:, :, iel), sol%u_tilde(ind_uf)) + elMat%L0(:, iel)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

#else
  !***********************************************************************
  !                            VERSION 2D
  !***********************************************************************
  integer*4             :: i, j, iel, N2d, Np, Nfp, Neq, Nf, Ndim, ierr, Nfdir
  integer*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq)
  integer*4             :: ind_uf(refElPol%Nfacenodes*phys%Neq*refElPol%Nfaces)
  integer*4             :: ind_ug(Mesh%Nnodesperelem*phys%Neq*Mesh%Ndim)
  integer*4             :: Fe(refElPol%Nfaces)

  if (utils%timing) then
    call cpu_time(timing%tps1)
    call system_clock(timing%cks1, timing%clock_rate1)
  end if

  Ndim = 2
  N2D = Mesh%Nelems
  Np = Mesh%Nnodesperelem
  neq = phys%neq
  Nfp = Mesh%Nnodesperface
  Nf = refElPol%Nfaces

  DO iel = 1, N2d
    ! Index for the face solution
    Fe = Mesh%F(iel, :)
    ind_uf = reshape(tensorSumInt((/(i, i=1, Neq*Nfp)/), (Fe - 1)*Neq*Nfp), (/Neq*Nfp*Nf/))

    ! Index for the element solution
    ind_ue = (iel - 1)*neq*Np + (/(i, i=1, Neq*Np)/)

    ! elemental solutions
    sol%u(ind_ue) = matmul(elMat%UU(:, :, iel), sol%u_tilde(ind_uf)) + elMat%U0(:, iel)

    ! Index for the element gradient
    ind_ug = (iel - 1)*neq*Np*Ndim + (/(i, i=1, Neq*Np*Ndim)/)

    ! elemental solutions
    sol%q(ind_ug) = matmul(elMat%LL(:, :, iel), sol%u_tilde(ind_uf)) + elMat%L0(:, iel)
  END DO

#endif


  if (utils%timing) then
    call cpu_time(timing%tpe1)
    call system_clock(timing%cke1, timing%clock_rate1)
    timing%runtsol = timing%runtsol + (timing%cke1 - timing%cks1)/real(timing%clock_rate1)
    timing%cputsol = timing%cputsol + timing%tpe1 - timing%tps1
  end if

END SUBROUTINE compute_element_solution

