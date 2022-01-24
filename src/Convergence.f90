PROGRAM Convergence
  USE in_out
  USE reference_element
  USE preprocess
  USE MPI_OMP
  USE printutils
  USE debug
  USE initialization
#ifdef WITH_PASTIX
  USE solve_pastix
#endif
#ifdef WITH_PSBLAS
  USE solve_psblas
#endif
  USE Postprocess
#ifdef PARALL
  USE Communications
#endif
  USE HDG_LimitingTechniques
  IMPLICIT NONE

  integer             :: it, ir, nts, nu, nb_args, n_points, n_polyns
  integer             :: ieq
  integer             :: ipol, ipts
  integer             :: Neq, Np, Nel, Ndim, Nfp, Nf, ntorloc
  integer             :: Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor
  integer             :: ifa, nut, k, i, j, g, is
  integer, allocatable :: faceNodes1d(:)
  real*8              :: dt, errNR, errTime
  real*8, allocatable :: uiter(:), u0(:)
  character(LEN=1024) :: cname, mesh_name, n_points_c, n_polyns_c
  integer, allocatable:: Nelv(:, :), Nunk(:, :)
  real*8, allocatable  :: L2err(:, :, :), L2err_pt(:), elSize(:, :), slope(:, :, :)

  ! Check the number of input arguments
  nb_args = iargc()
  IF (nb_args .lt. 1) THEN
    PRINT *, " Error: please provide the case name, the number of points and the polynomial degrees"
    stop
  ELSEIF (nb_args .lt. 2) THEN
    PRINT *, " Error: please provide the number of points and the polynomial degrees"
    stop
  ELSEIF (nb_args .lt. 3) THEN
    PRINT *, " Error: please provide the polynomial degrees"
    stop
  ELSEIF (nb_args .gt. 3) THEN
    PRINT *, " Error: too many input"
    stop
  END IF

  CALL getarg(1, cname)
  CALL getarg(2, n_points_c)
  CALL getarg(3, n_polyns_c)
  READ (n_points_c, *) n_points
  READ (n_polyns_c, *) n_polyns

  ! Initialize MPI
  CALL init_MPI_OMP()

  ! Read input file param.txt
  CALL read_input()

  ! Check input
  IF (switch%testcase > 10) THEN
    write (6, *) 'Testcase should be <10 for convergence check'
    stop
  ENDIF
#ifdef VORTICITY
  IF (.not. switch%convvort) THEN
    write (6, *) 'Consider convective term in vorticity equation for convergence check'
    stop
  ENDIF
#endif

  ! Initialization of the simulation parameters !TODO check if I need to pass dt to init_sim
  CALL init_sim(nts, dt)

  ! Allocate result vectors
  ALLOCATE (L2err_pt(phys%Neq))
  ALLOCATE (L2err(n_points, n_polyns, phys%Neq))
  ALLOCATE (elSize(n_points, n_polyns))
  ALLOCATE (Nelv(n_points, n_polyns))
  ALLOCATE (Nunk(n_points, n_polyns))
  ALLOCATE (Slope(n_points, n_polyns, phys%Neq))

  ! Initialization of the results
  L2err = 0.
  elSize = 0.
  Nelv = 0
  Nunk = 0
  slope = 0.

  WRITE(6,*) "Computing a convergence on the mesh ", trim(cname), "with ", trim(n_points_c),  " points and ", trim(n_polyns_c), " polynomial degrees"

  DO ipol = 1, n_polyns
    DO ipts = 1, n_points

      WRITE (mesh_name, '(A,A,I0,A,I0,A)') Trim(cname), '_', ipts, '_P', ipol
      WRITE (6, *) "Computing ", trim(mesh_name)

      ! Load the mesh file
      CALL load_mesh(Trim(mesh_name))

      ! Pastix: set the start to true
      matK%start = .true.

      ! Create the reference element based on the mesh type
      CALL create_reference_element(refElPol, 2)
#ifdef TOR3D
      !*************************************************
      !                REFERENCE ELEMENT 3D
      !*************************************************
      ! Create the reference element for the toroidal interpolation
      CALL create_reference_element(refElTor, 1, numer%ptor)

      !***********************************************
      ! Generate 3D basis functions
      !***********************************************
      ALLOCATE (refElTor%N3D(refElPol%NGauss2D*refElTor%NGauss1D, refElPol%Nnodes2D*refElTor%Nnodes1D))
      DO i = 1, refElTor%NGauss1D
        DO j = 1, refElPol%NGauss2D
          g = (i - 1)*refElPol%NGauss2D+j
          refElTor%N3D(g, :) = col(TensorProduct(refElPol%N2D(j, :), refElTor%N1D(i, :)))
        END DO
      END DO

      !***********************************************
      ! Generate face nodes for toroidal faces
      !***********************************************
      ALLOCATE (refElTor%faceNodes3(refElPol%Nfaces, size(refElPol%face_nodes, 2)*size(refElTor%face_nodes, 2)))
      ALLOCATE (faceNodes1d(refElTor%Nnodes1D))
      DO ifa = 1, refElPol%Nfaces
        faceNodes1d = (/(i, i=1, refElTor%Nnodes1D)/)
        refElTor%faceNodes3(ifa, :) = colint(TensorSumInt(refElPol%face_nodes(ifa, :), (faceNodes1d-1)*refElPol%Nnodes2D))
      END DO
      DEALLOCATE (faceNodes1d)

      !***********************************************
      ! Generate shape functions for toroidal faces
      !***********************************************
      Np1dPol = size(refElPol%N1d, 2)
      Np1dTor = size(refElTor%N1d, 2)
      Ng1dPol = size(refElPol%N1d, 1)
      Ng1dTor = size(refElTor%N1d, 1)
      refElTor%Nfl = Np1dPol*Np1dTor ! Number of nodes for lateral faces
      refElTor%Ngl = Ng1dPol*Ng1dTor ! Number of Gauss points for lateral faces
      refElTor%Nnodes3D = refElPol%Nnodes2D*refElTor%Nnodes1D ! Number of nodes in the 3D element
      refElTor%NGauss3D = refElPol%NGauss2D*refElTor%NGauss1D ! Number of Gauss points in the 3D element
      refElTor%Nft = refElTor%Nfl*refElPol%Nfaces + 2*refElPol%Nnodes2d         ! Number of nodes in all the faces of the 3D element
      ALLOCATE (refElTor%sFTF(Ng1dPol*Ng1dTor, Np1dPol*Np1dTor))
      DO i = 1, Np1dTor
        DO j = 1, Np1dPol
          k = (i - 1)*Np1dPol + j
          refElTor%sFTF(:, k) = col(tensorProduct(refElPol%N1d(:, j), refElTor%N1d(:, i)))
        END DO
      END DO
#endif

      ! Mesh preprocess: create the mesh related structures
      ! used in the HDG scheme
      CALL mesh_preprocess()

#ifdef TOR3D
      ! Define toroidal discretization
      CALL define_toroidal_discretization
#endif

      ! Initialize magnetic field (the Mesh is needed)
      CALL initialize_magnetic_field()

      ! Load magnetic field
      CALL load_magnetic_field()

#ifdef TOR3D
      Ndim = 3                                               ! Number of dimensions
      ntorloc = numer%ntor
      Nel = Mesh%Nelems*ntorloc                             ! Number of 3D elements
      Np = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
      Nfl = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
      Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
      Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
      nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
      Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
      nu = phys%Neq*Nel*Np
#else
      Ndim = Mesh%ndim
      Nel = Mesh%Nelems
      Np = refElPol%Nnodes2D
      Nf = refElPol%Nfaces
      Nfp = refElPol%Nfacenodes*Nf
      nut = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
      nu = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif

      ! Allocation and initialization of the elemental matrices
      CALL init_elmat()

      ! Compute the elemental matrices that depend only on the
      ! mesh, i.e. fixed in time
      CALL HDG_precalculatedfirstequation()

      ! Initialize the solution
      CALL init_sol()

      ! Allocate and initialize uiter and u0
      ALLOCATE(uiter(nu),u0(nu))
      ALLOCATE (sol%u0(nu, time%tis))
      sol%u0 = 0.
      sol%u0(:, 1) = sol%u
      uiter = 0.
      !u0 = 0.
      !*******************************************************
      !                  TIME LOOP
      !*******************************************************
      DO it = 1, nts

        ! Compute time step

        ! Actualization of time
        time%t = time%t + dt
        time%it = time%it + 1 ! is it of any use??
        time%ik = time%ik + 1
        WRITE (6, '(" *", 60("*"), "**")')
        WRITE (6, '(" *", 20X,    "Time iteration   = ", I5, 16X, " *")') time%it
        WRITE (6, '(" *", 60("*"), "**")')
        !*******************************************************
        !             Newton-Raphson iterations
        !*******************************************************
        uiter = sol%u0(:, 1)
        DO ir = 1, numer%nrp

          ! Compute the elemental matrices that depends on the N-R iteration
          CALL HDG_computeJacobian()

          ! Set boundary conditions
          CALL hdg_BC()

          ! Compute elemental mapping
          CAlL hdg_Mapping()

          ! Assembly the global matrix
          CALL hdg_Assembly()

          ! Solve linear system
          CALL solve_global_system()

          CALL compute_element_solution()

          ! Check convergence of Newton-Raphson
          errNR = computeResidual(sol%u, uiter, 1.)

          WRITE (6, *) "NR iteration: ", ir, ' - error: ', errNR
          IF (errNR .lt. numer%tNR) THEN
            EXIT
          ELSEIF (errNR .gt. numer%div) THEN
            WRITE (6, *) 'Problem in the N-R procedure'
            STOP
          ELSE
            uiter = sol%u
          END IF

        END DO ! End of NR loop
        WRITE (6, *) " "
        WRITE (6, *) " "

        ! Check convergence in time advancing and update
        errTime = computeResidual(sol%u, sol%u0(:, 1), time%dt)

        IF (errTime .lt. numer%tTM) THEN
          WRITE (6, *) "**********************"
          WRITE (6, *) "Time scheme converged!"
          WRITE (6, *) "**********************"
          EXIT
        ELSEIF (errTime .gt. numer%div) THEN
          WRITE (6, *) 'Problem in the time advancing scheme'
          STOP
        ELSE
          u0 = sol%u
        END IF
        WRITE (6, *) " "
      END DO ! End time loop
      DEALLOCATE (uiter, u0)

      ! Store results
      CALL computeL2ErrorAnalyticSol(L2err_pt)
      L2err(ipts, ipol, :) = L2err_pt
      elSize(ipts, ipol) = maxval(Mesh%elemSize)
      Nelv(ipts, ipol) = Mesh%Nelems
      Nunk(ipts, ipol) = MatK%n
      IF (ipts .gt. 1) THEN
        slope(ipts, ipol, :) = (log10(L2err(ipts, ipol, :)) - log10(L2err(ipts - 1, ipol, :)))/ &
          (log10(elSize(ipts, ipol)) - log10(elSize(ipts - 1, ipol)))
        !                                                                            (log10(sqrt(dble(Nunk(ipts-1,ipol))))-log10(sqrt(dble(Nunk(ipts,ipol)))))
      END IF
      DO ieq = 1, phys%Neq
        WRITE (6, *) "Error for ", Trim(phys%conVarNam(ieq)), ": ", L2err_pt(ieq)
      END DO

      ! Free allocated variables

      CALL free_all()
      time%it = 0

    END DO ! Loop in number of convergence points
  END DO ! Loop in number of polynomial degrees

  !*******************************************
  ! Display results
  !*******************************************
  ! I reinitialize to allocate some deallocates
  ! variables
  CALL init_sim(nts, dt)
  DO ieq = 1, phys%Neq

    WRITE (6, '(" * ", 55("-"), "*")')
    WRITE (6, '(" * Convergence for ", A7, 31X, " *")') Trim(phys%conVarNam(ieq))
    DO ipol = 1, n_polyns
      WRITE (6, '(" * ", 55("-"), "*")')
      WRITE (6, '(" * Interpolation: p", I0, 38X, "*")'), ipol
      WRITE (6, '(" * ",5X, " h ", 7X, " Nel ", 5X, " Nunk ", 2X, " L2err", 6X, "Slope", 5X, "*")')
      DO ipts = 1, n_points
        WRITE (6, '(" * ",2X,ES10.2E3, 2X, I5, 2X, I8, 2X, ES10.2E3, 2X, ES10.2E3, 2X, "*")'), &
          elSize(ipts, ipol), Nelv(ipts, ipol), Nunk(ipts, ipol), L2err(ipts, ipol, ieq), Slope(ipts, ipol, ieq)
      END DO
    END DO
    WRITE (6, '(" * ", 55("-"), "*")')
    WRITE (6, *) " "
  END DO

  DEALLOCATE (L2err, Nunk, elSize, Slope, Nelv)

contains

  !************************************************
  ! Compute the residual
  !************************************************
  FUNCTION computeResidual(u, uref, coeff) RESULT(res)
    real*8   :: u(nu), uref(nu)
    integer  :: nglo, ierr
    real     :: res, sum2, coeff

    sum2 = sum((u - uref)**2)
    nglo = nu

#ifdef PARALL
    call mpi_allreduce(MPI_IN_PLACE, sum2, 1, mpi_double_precision, mpi_sum, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(nu, nglo, 1, mpi_integer, mpi_sum, MPI_COMM_WORLD, ierr)
#endif
    res = sqrt(sum2)/sqrt(dble(nglo))/coeff
  END FUNCTION computeResidual

  !**********************************************
  ! Definition of the toroidal discretization
  !**********************************************
  SUBROUTINE define_toroidal_discretization
    integer :: i, ntorloc, itor, itorg, nnodes_toroidal
    integer :: ind(numer%ptor + 1)
    real*8  :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

    ! Toroidal discretization
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1, numer%ntor
      tdiv(i + 1) = i*htor
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

    nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
    Mesh%Nnodes_toroidal = nnodes_toroidal
    ALLOCATE (Mesh%toroidal(nnodes_toroidal))

    DO itor = 1, ntorloc
#ifdef PARALL
      itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
      if (itorg == numer%ntor + 1) itorg = 1
#else
      itorg = itor
#endif
      tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
      ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
      Mesh%toroidal(ind) = tel
    END DO
  END SUBROUTINE define_toroidal_discretization

END PROGRAM
