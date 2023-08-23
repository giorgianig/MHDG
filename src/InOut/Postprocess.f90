!*****************************************
! project: MHDG
! file: Postprocess.f90
! date: 13/02/2017
! Set of postprocess routines
!*****************************************
MODULE Postprocess

  USE globals
  USE analytical
  USE Printutils

CONTAINS

#ifdef TOR3D
  !********************************************
  !           3D toroidal version
  !********************************************
  SUBROUTINE computeL2ErrorAnalyticSol(L2err)
    real*8, intent(out)    :: L2err(:)
    integer*4             :: i, iel, Ndim, Nel, Np2D, iel3, itor, Neq
    integer*4             :: Np1Dpol, Np1Dtor, Np, Nfl, Nfp, Nf, ntor, N2D
    integer*4             :: ind_ue(refElPol%Nnodes2D*refElTor%Nnodes1D*phys%Neq)
    real*8                :: Xe(Mesh%Nnodesperelem, Mesh%Ndim)
    real*8                :: ue(refElPol%Nnodes2D*refElTor%Nnodes1D*phys%Neq)
    real*8                :: htor, tel(refElTor%Nnodes1d)
    real*8                :: tdiv(numer%ntor + 1)
    real*8                :: err_sol(phys%Neq), int_sol(phys%Neq)
    real*8                :: err_sole(phys%Neq), int_sole(phys%Neq)

    Ndim = 3                                               ! Number of dimensions
    Neq = phys%Neq
    N2D = Mesh%Nelems                                     ! Number elements in the 2D mesh
    Np2D = refElPol%Nnodes2D                               ! Number of nodes in the 2D elements
    Np1Dpol = refElPol%Nnodes1D                               ! Number of nodes in the 1D poloidal segments
    Np1Dtor = refElTor%Nnodes1D                               ! Number of nodes in the 1D toroidal segments
    Nel = N2D*numer%ntor                                  ! Number of 3D elements
    Np = Np2D*refElTor%Nnodes1D                          ! Number of nodes for each 3D element
    Nfl = Np1Dpol*Np1Dtor                                 ! Number of nodes in the lateral faces
    Nfp = Np2D*2 + refElPol%Nfaces*Nfl                      ! Number of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
    ntor = numer%ntor

    ! Initialization
    L2err = 0.
    err_sol = 0.
    int_sol = 0.

    ! Toroidal discretization
    htor = numer%tmax/ntor
    tdiv = 0.
    DO i = 1, ntor
      tdiv(i + 1) = i*htor
    END DO

    ! Element loop
    DO itor = 1, ntor
      tel = tdiv(itor) + 0.5*(refElTor%coord1d+1)*(tdiv(itor + 1) - tdiv(itor))
      DO iel = 1, N2D

        iel3 = (itor - 1)*N2d+iel

        ! Coordinates of the nodes of the element
        Xe = Mesh%X(Mesh%T(iel, :), :)

        ! Indices to extract the elemental and face solution
        ind_ue = (iel3 - 1)*Np*Neq + (/(i, i=1, Np*Neq)/)

        ! Elemental solution
        ue = sol%u(ind_ue)

        ! Compute the matrices for the element
        CALL elemental_error()

        ! Add elemental contribution
        err_sol = err_sol + err_sole
        int_sol = int_sol + int_sole
      END DO
    END DO

    !  L2err = sqrt(err_sol)
    L2err = sqrt(err_sol/int_sol)

  CONTAINS

    SUBROUTINE elemental_Error()
      USE LinearAlgebra
      USE analytical
      integer            :: g, NgaussPol, NgaussTor, igtor, igpol
      real*8             :: dvolu1d, dvolu
      real*8             :: Ni(Np), N1g(refElTor%Nnodes1D), N2g(refElPol%Nnodes2D)
      real*8             :: uan(refElPol%Ngauss2d*refElTor%Ngauss1d, Neq), ueg(Neq)
      real*8             :: xy(refElPol%Ngauss2d, Mesh%Ndim)
      real*8             :: J11(refElPol%Ngauss2d), J12(refElPol%Ngauss2d), J21(refElPol%Ngauss2d), J22(refElPol%Ngauss2d)
      real*8             :: iJ11(refElPol%Ngauss2d), iJ12(refElPol%Ngauss2d), iJ21(refElPol%Ngauss2d), iJ22(refElPol%Ngauss2d)
      real*8             :: detJ(refElPol%Ngauss2d), teg(refElTor%Ngauss1d)
      real*8, parameter  :: tol = 1e-12

      ! Initialization
      err_sole = 0.
      int_sole = 0.

      ! toroidal element size
      htor = tel(refElTor%Nnodes1d) - tel(1)

      ! Gauss points position in the poloidal plane
      xy = matmul(refElPol%N2D, Xe)

      ! Gauss points position in the toroidal direction
      teg = matmul(refElTor%N1d, tel)

      !                                                        ! Solution at Gauss points
      !                                                        ueg = matmul(refElPol%N2D,ue)

      ! Analytical solution at Gauss points
      CALL analytical_solution(iel,xy(:, 1), xy(:, 2), teg, uan)

      ! Error in the Gauss points
      !                                                        err_g = ueg-uan

      J11 = matmul(refElPol%Nxi2D, Xe(:, 1))                           ! ng x 1
      J12 = matmul(refElPol%Nxi2D, Xe(:, 2))                           ! ng x 1
      J21 = matmul(refElPol%Neta2D, Xe(:, 1))                          ! ng x 1
      J22 = matmul(refElPol%Neta2D, Xe(:, 2))                          ! ng x 1
      detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
      iJ11 = J22/detJ
      iJ12 = -J12/detJ
      iJ21 = -J21/detJ
      iJ22 = J11/detJ
      NgaussPol = refElPol%NGauss2D
      NgaussTor = refElTor%NGauss1D
      DO igtor = 1, NGaussTor
        N1g = refElTor%N1D(igtor, :)         ! Toroidal shape function
        dvolu1d = 0.5*refElTor%gauss_weights1D(igtor)*htor ! Toroidal contribution to the elemental volume
        DO igpol = 1, NGaussPol
          g = (igtor - 1)*NGaussPol + igpol
          ! Poloidal shape functions and derivatives
          N2g = refElPol%N2D(igpol, :)

          ! 3D shape functions
          Ni = col(TensorProduct(N2g, N1g))   ! 3D shape function

          ! Solution at courrent Gauss point
          ueg = matmul(Ni, transpose(reshape(ue, [neq, Np])))

          ! 3D integration weight
          dvolu = refElPol%gauss_weights2D(igpol)*detJ(igpol)*dvolu1d

          IF (switch%axisym) THEN
            dvolu = dvolu*xy(igpol, 1)
          ENDIF

          ! Contribution of the current integration point to the elemental error
          err_sole = err_sole + (uan(g, :) - ueg)**2*dvolu
          int_sole = int_sole + uan(g, :)**2*dvolu
        END DO
      END DO

    END SUBROUTINE elemental_Error

  END SUBROUTINE computeL2ErrorAnalyticSol
#else
  !********************************************
  !           2D version
  !********************************************
  SUBROUTINE computeL2ErrorAnalyticSol(L2err)
    real*8, intent(out)    :: L2err(phys%Neq)
    integer*4             :: i, iel, Ndim, Neq, Nel, Np
    integer*4             :: ind_ue(Mesh%Nnodesperelem*phys%Neq)
    real*8                :: Xe(Mesh%Nnodesperelem, Mesh%Ndim)
    real*8                :: ue(Mesh%Nnodesperelem*phys%Neq)
    real*8                :: err_sol(phys%Neq), int_sol(phys%Neq)
    real*8                :: err_sole(phys%Neq), int_sole(phys%Neq)

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D

    ! Initialization
    L2err = 0.
    err_sol = 0.
    int_sol = 0.

    ! Element loop
    DO iel = 1, Nel

      ! Coordinates of the nodes of the element
      Xe = Mesh%X(Mesh%T(iel, :), :)

      ! Indices to extract the elemental and face solution
      ind_ue = (iel - 1)*Np*Neq + (/(i, i=1, Neq*Np)/)

      ! Elemental solution
      ue = sol%u(ind_ue)

      ! Compute the matrices for the element
      CALL elemental_error()

      ! Add elemental contribution
      err_sol = err_sol + err_sole
      int_sol = int_sol + int_sole
    END DO

    L2err = sqrt(err_sol/int_sol)

  CONTAINS

    SUBROUTINE elemental_Error()
      integer            :: g, Ngauss
      real*8             :: uan(refElPol%NGauss2D, Neq), ueg(refElPol%NGauss2D, Neq)
      real*8             :: xy(refElPol%Ngauss2d, ndim)
      real*8             :: Jacob(ndim, ndim), detJ, dvolu
      real*8             :: err_g(refElPol%Ngauss2d, Neq)
      real*8, parameter  :: tol = 1e-12

      ! Initialization
      err_sole = 0.
      int_sole = 0.

      ! Gauss points position
      xy = matmul(refElPol%N2D, Xe)

      ! Solution at Gauss points
      ueg = matmul(refElPol%N2D, transpose(reshape(ue, [neq, Np])))

      ! Analytical solution at Gauss points
      CALL analytical_solution(iel,xy(:, 1), xy(:, 2), uan)

      ! Error in the Gauss points
      err_g = ueg - uan

      ! Loop in 2D Gauss points
      Ngauss = refElPol%NGauss2D
      DO g = 1, NGauss

        ! Jacobian
        Jacob(1, 1) = dot_product(refElPol%Nxi2D(g, :), Xe(:, 1))
        Jacob(1, 2) = dot_product(refElPol%Nxi2D(g, :), Xe(:, 2))
        Jacob(2, 1) = dot_product(refElPol%Neta2D(g, :), Xe(:, 1))
        Jacob(2, 2) = dot_product(refElPol%Neta2D(g, :), Xe(:, 2))
        detJ = Jacob(1, 1)*Jacob(2, 2) - Jacob(1, 2)*Jacob(2, 1)
        IF (detJ < tol) THEN
          error stop "Negative jacobian"
        END if

        ! Integration weight
        dvolu = refElPol%gauss_weights2D(g)*detJ
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(g, 1)
        ENDIF

        ! Contribution of the current integration point to the elemental error
        err_sole = err_sole + err_g(g, :)**2*dvolu
        int_sole = int_sole + uan(g, :)**2*dvolu
      END DO

    END SUBROUTINE elemental_Error

  END SUBROUTINE computeL2ErrorAnalyticSol

#endif
END MODULE Postprocess
