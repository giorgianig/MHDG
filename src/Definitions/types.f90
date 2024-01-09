!************************************************************
! project: MHDG
! file: types.f90
! date: 06/09/2016
! Declaration of all the types used
! in the code
!************************************************************

MODULE types
  USE prec_const
  IMPLICIT NONE

  !*******************************************************
  ! Boundary condition flags
  !*******************************************************
  INTEGER, PARAMETER, PUBLIC :: max_num_diff_bc = 100
  ! 1-10 Dirichlet type
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet = 1
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_weak_form = 2
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_weak_form_oldvalues = 3
  INTEGER, PARAMETER, PUBLIC :: bc_Transmission = 4
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_and_Neumann = 5
   INTEGER, PARAMETER, PUBLIC :: bc_iter_core = 6
  ! 20-.. In out type
  INTEGER, PARAMETER, PUBLIC :: bc_inout = 20
  ! 30-.. Neumann type
  INTEGER, PARAMETER, PUBLIC :: bc_NeumannH = 30
  ! 50-.. Bohm type
  INTEGER, PARAMETER, PUBLIC :: bc_Bohm = 50
  INTEGER, PARAMETER, PUBLIC :: bc_BohmNM = 51
  INTEGER, PARAMETER, PUBLIC :: bc_BohmSimp = 52
  INTEGER, PARAMETER, PUBLIC :: bc_BohmPump = 55 ! bohm for plasma, pump for neutrals
  INTEGER, PARAMETER, PUBLIC :: bc_BohmPuff = 56 ! bohm for plasma, puff for neutrals

  ! 60-.. Slip wall type
  INTEGER, PARAMETER, PUBLIC :: bc_slip_wall = 60
  ! 70-.. Periodic
  INTEGER, PARAMETER, PUBLIC :: bc_periodic = 70


  ! Boundary conditions flags and names
  CHARACTER(100) ::  bc_flag_type(10), bc_flag_name(max_num_diff_bc)

  !*******************************************************
  ! Reference element
  !*******************************************************
  TYPE :: Reference_element_type
    integer*4 :: ndim          ! Number of dimensions of the reference element
    integer*4 :: elemType      ! type of mesh elements (0 = triangles, 1 = quadrilaterals, 2 = thetrahedra, 3 = hexahedra)
    integer*4 :: Nvertices     ! Number of vertices of the element: 3 for triangles, 4 for quads and thetraedra, 8 for hexahedra
    integer*4 :: Ndeg          ! degree of interpolation
    integer*4 :: Nnodes3D      ! number of 3D nodes in the element
    integer*4 :: Nnodes2D      ! number of 2D nodes in the element
    integer*4 :: Nnodes1D      ! number of 1D nodes in the element
    integer*4 :: Nfaces        ! number of faces in the element: a face in a 2D mesh is a segment, in a 3D mesh is a 2D entity
    integer*4 :: Nbords        ! number of bords (only in 3D), is the number of 1D entities of the element
    integer*4 :: Nextnodes     ! number of exterior nodes of the element
    integer*4 :: Nfacenodes    ! number of nodes in each face of the element
    integer*4 :: Nfacenodeslin ! number of nodes in each linear face of the element
    integer*4 :: Nbordnodes    ! number of nodes in each bord of the element (only in 3D)
    integer*4 :: Ninnernodes   ! number of inner nodes of the element
    integer*4 :: Ninnernodesface  ! number of inner nodes of the 2D faces (only for 3D meshes)
    integer*4 :: NGauss1D      ! number of Gauss points of the 1D quadrature
    integer*4 :: NGauss2D      ! number of Gauss points of the 2D quadrature
    integer*4 :: NGauss3D      ! number of Gauss points of the 3D quadrature
    integer*4, allocatable :: Face_nodes(:, :)  ! numbering of the face nodes
    integer*4, allocatable :: Bord_nodes(:)    ! numbering of the 1D edge nodes (for 3D meshes)
    integer*4, allocatable :: inner_nodes(:)   ! numbering of the inner nodes
    integer*4, allocatable :: inner_nodes_face(:)   ! numbering of the inner nodes for 2D faces (for 3D meshes)
    real*8, pointer :: coord3d(:, :) => null()       ! spatial coordinates of the 3D nodes
    real*8, pointer :: coord2d(:, :) => null()       ! spatial coordinates of the 2D nodes
    real*8, pointer :: coord1d(:) => null()       ! spatial coordinates of the 1D nodes
    real*8, allocatable :: gauss_points3D(:, :) ! Gauss integration points for the 3D quadrature
    real*8, allocatable :: gauss_points2D(:, :) ! Gauss integration points for the 2D quadrature
    real*8, allocatable :: gauss_points1D(:)   ! Gauss integration points for the 1D quadrature
    real*8, allocatable :: gauss_weights3D(:)  ! Weights of the integration points for the 3D quadrature
    real*8, allocatable :: gauss_weights2D(:)  ! Weights of the integration points for the 2D quadrature
    real*8, allocatable :: gauss_weights1D(:)  ! Weights of the integration points for the 1D quadrature
    real*8, allocatable :: N3D(:, :)      ! Shape functions 3D at the Gauss points
    real*8, allocatable :: Nxi3D(:, :)    ! Derivative with respect to xi of the shape functions 3D at the Gauss points
    real*8, allocatable :: Neta3D(:, :)   ! Derivative with respect to eta of the shape functions 3D at the Gauss points
    real*8, allocatable :: Nzeta3D(:, :)  ! Derivative with respect to zeta of the shape functions 3D at the Gauss points
    real*8, allocatable :: N2D(:, :)      ! Shape functions 2D at the Gauss points
    real*8, allocatable :: Nxi2D(:, :)    ! Derivative with respect to xi of the shape functions 2D at the Gauss points
    real*8, allocatable :: Neta2D(:, :)   ! Derivative with respect to eta of the shape functions 2D at the Gauss points
    real*8, allocatable :: N1D(:, :)      ! Shape functions 1D at the Gauss points
    real*8, allocatable :: Nxi1D(:, :)    ! Derivative of the shape functions 1D at the Gauss points
    real*8, allocatable :: Nlin(:, :)     ! Linear shape functions 2D at the nodes of the element (only used in shock capturing so far...)
    real*8, allocatable :: sFTF(:, :)     ! Shape functions for the toroidal faces
    integer*4, allocatable :: faceNodes3(:, :) ! Face numbering for toroidal faces
    integer*4          :: Nfl           ! Number of nodes in lateral faces for toroidal 3d computations
    integer*4          :: Ngl           ! Number of Gauss points in lateral faces for toroidal 3d computations
    integer*4          :: Nft           ! Number of nodes in all faces for toroidal 3d computations
  END TYPE Reference_element_type

  !*******************************************************
  ! Mesh
  !*******************************************************
  TYPE :: Mesh_type
    integer*4 :: Ndim                ! Number of dimensions of the mesh
    integer*4 :: Nnodes              ! Number of nodes in the mesh
    integer*4 :: Nnodesperelem       ! Number of nodes per element in the mesh (dimension 2 of the connectvity matrix)
    integer*4 :: Nnodesperface       ! Number of nodes per face in the mesh (dimension 2 of boundary connectvity matrix )
    integer*4 :: Nelems              ! Number of elements in the mesh
    integer*4 :: Nfaces              ! Number of faces in the mesh
    integer*4 :: Nextfaces           ! Number of exterior faces
    integer*4 :: Nintfaces           ! Number of interior faces
    integer*4 :: elemType            ! 0 for quads - 1 for triangles - 2 for thetra - 3 for hexa
    integer*4 :: ndir                ! number of Dirichlet faces
    integer*4 :: ukf                 ! number of faces in the mesh minus the number of Dirichlet faces
    integer*4, pointer :: T(:, :) => null()    ! Elements connectivity matrix
    integer*4, pointer :: Tlin(:, :) => null()    ! Elements linear connectivity matrix
    integer*4, pointer :: Tb(:, :) => null()    ! Outer faces connectivity matrix
    integer*4, pointer :: boundaryFlag(:)     ! Flag for the boundary condition for each external face (set in the mesh generator)
    integer*4, allocatable :: F(:, :)          ! Faces connectivity matrix
    integer*4, allocatable :: N(:, :)          ! Nodes connectivity matrix
    integer*4, allocatable :: faces(:, :, :)    ! for each triangle i, stores info k on each face j: faces(i,j,1) = # of neighbouring triangle (0 if external
    ! boundary), faces(i,j,2) = type of boundary (), faces(i,j,3) = type of boundary condition
    integer*4, allocatable :: extfaces(:, :)   ! for each exterior face, stores the number of the triangle, the number of the face, and the type of BC
    integer*4, allocatable :: intfaces(:, :)   ! for each interior face, stores the number of the triangle, the number of the face, the number of the
    ! neighboring triangle, the number of its face and the number of the node of the neighboring triangle that
    ! matches the first knot of the triangle
    logical, allocatable :: flipface(:, :)    ! for each triangle, and for each face, 0 if the order of the numbering in the face is to be kept, 1 if the
    ! order is to be reversed
    logical, allocatable    :: Fdir(:, :)         ! for each element states if each local face is of Dirichlet type
    integer*4, allocatable  :: periodic_faces(:)  ! Mapping for periodic faces
    integer*4, allocatable  :: Diric(:)
    integer*4, allocatable  :: numberbcs(:)
    real*8, allocatable     :: elemSize(:)   ! element size (area in 2D, volume in 3D) [Number of elements]
    real*8, pointer         :: X(:, :) => null()    ! nodes coordinates
    integer*4              :: Nnodes_toroidal     ! Number of nodes in the toroidal direction
    real*8, pointer         :: toroidal(:) => null()    ! nodes coordinates in the toroidal direction
    ! Limiting & shock capturing stuff
    integer, allocatable    :: flag_elems_rho(:)      ! Flagged elements for limiting rho [Number of elements]
    integer, allocatable    :: flag_elems_sc(:)     !  Flagged elements for shock-capturing [Number of elements]
    real*8, allocatable     :: minrho_elems(:)        ! Minimum value of density in the flagged elements[Number of elements]
    real*8, allocatable     :: sour_elems(:)          ! Source to limit rho in the flagged elements[Number of elements]
    real*8, allocatable     :: diff_elems(:)          ! Diffusion to limit rho in the flagged elements[Number of elements]
    real*8, allocatable     :: scdiff_nodes(:, :)      ! Shock capturing diffusion in each node [Number of elements,Number of nodes per element]
    real*8                 :: xmax, xmin, ymax, ymin    ! Limit of the GLOBAL matrix, across mpi partitions
    real*8                  :: puff_area         ! area of the puff bounday condition
    real*8                  :: core_area         ! area of the puff bounday condition
    real*8,allocatable      :: Xg(:,:)               ! 2D Gauss point coordinates
    real*8,allocatable      :: Xgf(:,:)              ! 1D Gauss point coordinates at interior faces
    real*8,allocatable      :: Xgb(:,:)          ! 1D Gauss point  coordinates at boundary faces
#ifdef PARALL
    integer, pointer       :: loc2glob_fa(:)     ! mapping number of the faces for creating the global matrix [number of faces in the mesh]
    integer, pointer       :: loc2glob_el(:)     ! mapping number of the elements from local to global [number of elements in the mesh]
    integer, pointer       :: ghostfaces(:)      ! integer that states if a face is to be assembled locally or not [number of faces in the mesh]
    integer, pointer       :: ghostelems(:)      ! integer that states if an element is to be assembled locally or not (for 3D) [number of elements in the mesh]
    integer               :: nghostfaces        ! number of ghost faces
    integer               :: nghostelems        ! number of ghost elements (used only for 3D)
    integer               :: Nel_glob           ! number of elements of global mesh
    integer               :: Nfa_glob           ! number of faces of global mesh
    integer               :: Nno_glob           ! number of nodes of global mesh
    integer               :: Ndir_glob          ! number of Dirichlet faces in the global mesh
    integer               :: Ngho_glob          ! number of ghost faces in the global mesh
    ! readed from input
    integer, pointer       :: ghostflp(:)        ! flipFaces for the ghost faces [number of ghost faces]
    integer, pointer       :: ghostloc(:)        ! local numbering of the ghost face in the process that assemble it [number of ghost faces]
    integer, pointer       :: ghostpro(:)        ! the process that assemble the ghost face [number of ghost faces]
    integer, pointer       :: ghelsloc(:)        ! local numbering of the ghost element in the process that assemble it (only 3D) [number of ghost elements]
    integer, pointer       :: ghelspro(:)        ! the process that assemble the ghost element (only 3D) [number of ghost elements]
    ! built after reading from input
    integer, allocatable   :: fc2sd(:)          ! face 2 send: faces computed locally that the local process has to send (vector)
    integer, allocatable   :: pr2sd(:)          ! process 2 send: to which process the faces computed locally need to be sent (vector)
    integer, allocatable   :: fc2rv(:)          ! face 2 receive: ghost faces computed by other processes that the local process need to receive[number of ghost faces]
    integer, allocatable   :: pr2rv(:)          ! process 2 receive: from which process the faces computed externally need to be received [number of ghost faces] (it is the same as ghostpro)

    integer, allocatable   :: el2sd(:)          ! element 2 send: elements computed locally that the local process has to send (only 3D) (vector)
    integer, allocatable   :: pe2sd(:)          ! process 2 send: to which process the elements computed locally need to be sent (only 3D) (vector)
    integer, allocatable   :: el2rv(:)          ! element 2 receive: ghost elements computed by other processes that the local process need to receive (only 3D) [number of ghost elements]
    integer, allocatable   :: pe2rv(:)          ! process 2 receive: from which process the elements computed externally need to be received (only 3D) [number of ghost elements] (it is the same as ghelspro)

    !     integer,allocatable   :: connpro(:)        ! processes connected to the local process
#endif
  END TYPE Mesh_type

  !*******************************************************
  ! Physics: type for physical model info
  !*******************************************************
  TYPE Physics_type
    integer         :: neq            ! Number of equations
    integer         :: npv            ! Number of physical variables
    real*8          :: diff_n, diff_u ! Perpendicular diffusion in the continuity and momentum equation
    real*8          :: nu0 = 0.       ! Max collision frequency for rho_n < 0.95: for selconsistent core pinch
    real*8          :: v_p            ! Pinch velocity in the continuity equation
    real*8          :: v_pmax = 0.    ! Max pinch velocity
    real*8          :: a              ! Proportionality constant between pressure and density for isothermal model (p = a*rho)
    real*8          :: dfcoef         ! Constant related to the diamagnetic drift velocity
    real*8          :: dexbcoef         ! Constant related to the ExB drift velocity
    integer*4       :: bcflags(1:10)          ! Set the correspondence between mesh boundary flag (Mesh%boundaryFlag) and the boundary condition
    real*8          :: diagsource(1:10)       ! Diagonal implicit sources
    character(LEN=20), pointer:: phyVarNam(:) => Null() ! Names of the physical variables (set in initPhys)
    character(LEN=20), pointer:: conVarNam(:) => Null() ! Names of the conservative variables (set in initPhys)
    real*8          :: lscale          ! Length scale for the non-dimensionalization of the equations
    ! Magnetic field defined for each node of the mesh.
    real*8, pointer :: B(:, :)            ! Magnetic field, Br,Bz,Bphi [n of nodes  x 3]
    real*8          :: B0                 ! Reference value for the magnetic field [Tesla]
    real*8, pointer :: magnetic_flux(:)   ! Magnetic flux   [n of nodes]
    real*8, pointer :: magnetic_psi(:)    ! Magnetic flux normalized to separatrix magnetic flux   [n of nodes]
    real*8, pointer :: safety_factor(:)   ! Safety factor profile as a function of psi   [n of nodes]
    real*8          :: Flux2Dmin          ! Minimum of the magnetic flux, across the MPI partitions
    real*8          :: Flux2Dmax          ! Maximum of the magnetic flux, across the MPI partitions
    real*8, pointer :: Bperturb(:, :)     ! Magnetic perturbation, Br,Bz,Bphi [n of nodes  x 3]
    real*8          :: Tbg                ! Background temperature in the isothermal model
    real*8, pointer :: Jtor(:)            ! Toroidal Current
    real*8          :: I_p                ! Total plasma current
    real*8          :: bohmth             ! Threshold for imposing the Bohm boundary condition
    ! Energy equation coefficients
    real*8          :: diff_e             ! Perpendicular diffusion in the energy equation
    real*8          :: epn                ! Exponential of the parallel diffusion (usually 5/2)
    real*8          :: Mref               ! Reference Mach number
    real*8          :: diff_pari          ! Parallel diffusion for the temperature (usually 1e7)
    real*8          :: Gmbohm             ! gamma for Bohm boundary condition on energy:
    ! Temperature equations coefficients (the ions coefficients are the ones defined previously)
    real*8          :: diff_ee            ! Perpendicular diffusion in the elcetron energy equation
    real*8          :: diff_pare          ! Parallel diffusion for the electron temperature
    real*8          :: tie                ! Temperature exchange coefficient between ions and electrons
    real*8          :: Gmbohme            ! gamma for Bohm boundary condition on electron energy:
    real*8          :: Pohmic             ! Ohmic heating power
    ! Coefficients for the vorticity equations
    real*8          :: diff_vort
    real*8          :: diff_pot
    real*8          :: etapar
    real*8          :: c1, c2             ! coefficients coming from the adimensionalization
    real*8          :: Potfloat
    ! Coefficients for the neutral equations
    real*8          :: diff_nn            ! Diffusion in the neutral equation
    real*8,allocatable:: diff_nn_Vol(:)   ! Diffusion in the neutral equation at 2D Gauss points
    real*8,allocatable:: diff_nn_Fac(:)   ! Diffusion in the neutral equation at 1D Gauss points on interior faces
    real*8,allocatable:: diff_nn_Bou(:)   ! Diffusion in the neutral equation at 1D Gauss points on boundary faces 
    real*8,allocatable:: v_nn_Vol(:,:)    ! Convective velocity in the neutral equation at 2D Gauss points
    real*8,allocatable:: v_nn_Fac(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on interior faces
    real*8,allocatable:: v_nn_Bou(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on boundary faces 
    real*8          :: Re                 ! Recycling for the neutral equation
    real*8          :: puff               ! Puff coefficient
    real*8          :: puff_slope         ! Puff increment coefficient (only for moving equilibrium for ITER)
    real*8,pointer  :: puff_exp(:)        ! Puff experimental coefficient (only for moving equilibriums)
    real*8,pointer  :: n_li(:)            ! Central line integrated density (only for moving equilibriums)
    real*8,pointer  :: n_lit(:)           ! Central line integrated density target (only for moving equilibriums)
    real*8          :: part_source        ! Particle source for ITER
    real*8          :: ener_source        ! Particle source for ITER
    real*8          :: density_source     ! Density source for WEST (2D, case 52)
    real*8          :: ener_source_e      ! Ion energy source for WEST and ITER (2D, case 52 and 81)
    real*8          :: ener_source_ee     ! Electron source for WEST and ITER (2D, case 52 and 81)
    real*8          :: sigma_source       ! Sigma for the gaussian sources for WEST and ITER (2D, case 52 and 81)
    real*8          :: fluxg_trunc        ! Value of the NORMALISED magnetic flux at which to truncate the gaussian sources for WEST (2D, case 52), refer to Source_shape.m file
    ! Diffusion coefficients ITER evolving equilibrium
    real*8          :: ME_diff_n
    real*8          :: ME_diff_u
    real*8          :: ME_diff_e
    real*8          :: ME_diff_ee
  END TYPE Physics_type

  !*******************************************************
  ! Geometry: type for geometry info
  !*******************************************************
  TYPE Geometry_type
    integer     :: ndim     ! Number of dimensions of the problem (can be different from the ndim of the mesh)
    real*8      :: R0       ! Major radius at the magnetic axis
    real*8      :: q        ! Safety factor
  END TYPE Geometry_type

  !*******************************************************
  ! Magnetic: type for 3D magnetic perturbation
  !*******************************************************
  TYPE Magnetic_type
    real*8          :: amp_rmp            ! amplitude RMP
    integer         :: nbCoils_rmp        ! number coils RMP (full torus)
    real*8          :: torElongCoils_rmp  ! Toroidal elongation of RMP coils (rectangular coils)
    integer         :: parite             ! parite RMP (for 2 row, -1 even parite, 1 odd)
    integer         :: nbRow              ! number of rows (1,2 or 3) of RMP coils, default is 2
    real*8          :: amp_ripple         ! amplitude Ripple
    integer         :: nbCoils_ripple     ! number coils Ripple (full torus)
    real*8          :: triang             ! triangularity (0: None)
    real*8          :: ellip              ! ellipticity (1: None)
    real*8, pointer :: coils_rmp(:, :, :) ! Coil coordinates for RMP (nbCoils*4*Discr,start-stop*(xyz)=6,rowNb) (4 for square coils)
    real*8, pointer :: coils_ripple(:, :) ! Coil coordinates for Ripple (nbCoils*Discr,start-stop*(xyz)=6)
  END TYPE Magnetic_type

  !*******************************************************
  ! Switches: type for main code switches
  !*******************************************************
  TYPE Switches_type
    logical :: axisym ! Is it an axisymmetric simulation?
    ! true =
    ! false =
    logical :: rmp      ! To activate resonant magnetic perturbation
    logical :: ripple   ! To activate ripple
    logical :: ohmicsrc ! Set to TRUE to consider ohmic source of energy
    logical :: ME       ! Set to TRUE to allow magnetic equilibrium evolution in time
    logical :: driftdia ! Set to TRUE to consider diamagnetic drift
    logical :: driftexb ! Set to TRUE to consider ExB drift
    logical :: steady
    logical :: time_init ! true if it is a time initialization simulation. The time counter "it" does not increment  (i.e. when the analitical initialisation is not good enough). Used for moving equilibrium (case 59)
    integer :: init     ! 1-init. analy. solution at nodes; 2-L2 projection
    ! Set to TRUE for a steady state computation
    ! Set to FALSE for a transient computation
    integer :: testcase  ! Define the testcase ( which analytical solution, body force, magnetic field...)
    logical :: psdtime ! Reduce the diffusion every time we reach the steady state
    ! condition (only works if steady=.false.)
    real*8  :: diffred ! Reduction factor of the diffusion for psdtime simulation
    real*8  :: diffmin ! Minimum value of the diffusion for a psdtime simulation
    integer :: shockcp ! Shock capturing option
    integer :: limrho  ! Add a source for limiting the min value of rho
    integer :: difcor  ! Add diffusion in corners
    integer :: thresh  ! Use a threshold for limiting the min value of rho
    ! (rho-rho*E in case of N-Gamma-Energy model, rho-rho*Ei-rho*Ee in case of N-Gamma-Ti-Te model)
    logical :: filter  ! Filter solution to avoid oscillation in empty zones
    logical :: decoup  ! Decouple N-Gamma from Te-Ti (only used for N-Gamma-Ti-Te model)
    logical :: ckeramp ! Chech the error amplification in the linear system solution (for very ill-conditioned matrices)
    logical :: saveNR  ! Save solution at each NR iteration
    logical :: saveTau ! Save tau on faces
    logical :: fixdPotLim
    logical :: dirivortcore
    logical :: dirivortlim
    logical :: convvort ! consider the convective term in the vorticity equation
    logical :: bxgradb  ! consider the term in BxGradB in the vorticity equation
    integer :: pertini  ! add perturbation in the initial solution
    ! 1 -add sinusoidal perturbation
    ! 2 -add density blob
    logical :: logrho   ! solve for the density logarithm instead of density
  END TYPE Switches_type

  !*******************************************************
  ! Time: type for the time stepping information
  !*******************************************************
  TYPE Time_type
    real*8      :: dt0  ! initial time step
    real*8      :: dt   ! current time step
    real*8      :: tfi  ! final time of the simulation
    integer     :: it   ! the number of the current time step
    integer     :: ik   ! same as it but always incrementing (also in case of pseudotime..)
    integer     :: ndt  ! max number of time steps to do in the current session
    integer     :: tsw  ! switch to modify the time step
    integer     :: nts  ! max number of time iterations to do in the current session (only for transient simulations)
    integer     :: tis  ! time integration scheme
    ! 1 - first order
    ! 2 - second order
    real*8      :: t    ! time of the simulation (initialized to finish time of previous simulation if restart, to 0 if new simulation)
  END TYPE Time_type

  !*******************************************************
  ! Numerics: type for numeric scheme parameters
  !*******************************************************
  TYPE Numeric_type
    integer        :: nrp      ! Max number of Newton-Raphson iterations
    real*8         :: tNR      ! Tolerance of the Newton-Raphson scheme
    real*8         :: tTM      ! Tolerance for the steady state achievement
    real*8         :: div      ! Divergence detector
    real*8         :: tau(1:6) ! Stabilization parameter for each equation (4 values max for now...)
    real*8         :: sc_coe   ! Shock capturing coefficient
    real*8         :: sc_sen   ! Shock capturing sensibility
    real*8         :: minrho   ! Value of rho to start applying limiting
    real*8         :: so_coe   ! Source coefficient for limiting rho
    real*8         :: df_coe   ! Diffusion coefficient for limiting rho
    real*8         :: dc_coe   ! Diffusion coefficient in corners
    real*8         :: thr      ! Threshold to limit rho
    real*8         :: thrpre   ! Threshold to limit pressure
    integer        :: stab     ! Stabilization type
    ! 1 - constant tau (one for each equation) in the whole domain
    ! 2 -
    ! 3 -
    real*8         :: dumpnr   ! dumping factor for Newton-Raphson. 0<dumpnr<1
    integer        :: ntor     ! Number of elements in the toroidal direction
    integer        :: ptor     ! Polynomial degree in the toroidal direction
    real*8         :: tmax     ! Max extention in the toroidal direction
    integer        :: npartor  ! Number of MPI divisions in the toroidal direction
    integer        :: bohmtypebc ! Implementation of the Bohm bc for Gamma
    real*8         :: exbdump ! Dumping for ExB drifts
  END TYPE Numeric_type

  !*******************************************************
  ! Utilities: type for printing/debugging/saving...
  !*******************************************************
  TYPE Utils_type
    integer :: printint       ! Integer for printing
    logical :: timing         ! Timing of the code
    integer :: freqdisp       ! Frequency of results display
    integer :: freqsave       ! Frequency of solution save
  END TYPE Utils_type

  !*******************************************************
  ! Linear system solver parameters
  !*******************************************************
  TYPE Lssolver_type
    integer           :: sollib    ! Solver library to be used
    ! 1-Pastix
    ! 2-PSBLAS
    logical           :: timing    ! timing of the linear solver
    ! Parameters relative to the library PSBLAS
    character(len=20) :: kmethd    ! Krylov method (see list on the library manual)
    integer           :: istop     ! Stopping criterion (see spec on the library manual)
    integer           :: itmax     ! Max number of iterations
    integer           :: itrace    ! Display convergence at each iteration
    integer           :: rest      ! Restart
    real*8            :: tol       ! Stopping tolerance
    character(len=20) :: ptype     ! Preconditioner type
    ! Parameters relative to the library MLD2P4
    ! First smoother / 1-lev preconditioner
    character(len=20) :: smther       ! smoother type
    integer           :: jsweeps      ! (pre-)smoother / 1-lev prec sweeps
    integer           :: novr         ! Number of olverlap layers, for Additive Schwartz only
    character(len=20) :: restr        ! restriction  over application for Additive Schwartz only
    character(len=20) :: prol         ! Type of prolongation operator for Additive Schwartz only
    character(len=20) :: solve        ! local subsolver
    integer           :: fill         ! fill-in level p of the ILU factorizations
    real              :: thr          ! threshold for ILUT
    ! Second smoother/ AMG post-smoother (if NONE ignored in main)
    character(len=20) :: smther2      ! smoother type
    integer           :: jsweeps2     ! (post-)smoother sweeps
    integer           :: novr2        ! number of overlap layers
    character(len=20) :: restr2       ! restriction  over application of AS
    character(len=20) :: prol2        ! prolongation over application of AS
    character(len=20) :: solve2       ! local subsolver
    integer           :: fill2        ! fill-in for incomplete LU
    real              :: thr2         ! threshold for ILUT
    ! general AMG data
    character(len=20) :: mlcycle      ! multi-level cycle
    integer           :: outer_sweeps ! number of multilevel cycles
    integer           :: maxlevs      ! Maximum number of levels
    integer           :: csize        ! Coarse size threshold
    ! aggregation
    real              :: mncrratio    ! Minimum coarsening ratio
    real              :: athres       ! Aggregation threshold
    character(len=20) :: aggr_prol    ! Prolongator used by the aggregation algorithm
    character(len=20) :: par_aggr_alg ! Parallel aggregation algorithm
    character(len=20) :: aggr_ord     ! Initial ordering of indices for the aggregation algorithm
    character(len=20) :: aggr_filter  ! Matrix used in computing the smoothed prolongator
    ! coasest-level solver
    character(len=20) :: csolve       ! coarsest-lev solver
    character(len=20) :: csbsolve     ! coarsest-lev solver
    character(len=20) :: cmat         ! coarsest mat layout
    integer           :: cfill        ! fill-in for incompl LU
    real              :: cthres       ! Threshold for ILUT
    integer           :: cjswp        ! sweeps for GS/JAC subsolver
  END TYPE Lssolver_type

  !**********************************************************
  ! Solution: contains the solution at the current time step
  !**********************************************************
  TYPE Sol_type
    real*8, pointer :: u(:) => null()  ! Elemental solution
    real*8, pointer :: u_tilde(:) => null()  ! Face solution
    real*8, pointer :: q(:) => null()  ! Elemental solution for the gradient
    real*8, allocatable :: u0(:, :)           ! Elemental solution at previous time steps
    real*8, allocatable :: tres(:)           ! Time residual
    real*8, allocatable :: time(:)           ! Time evolution
    integer             :: Nt                ! Number of time-steps
  END TYPE Sol_type

  !**********************************************************
  ! Simulation parameters: for saving purpose
  !**********************************************************
  TYPE Simulationparams_type
    character(len=50) :: model
    integer   :: Ndim
    integer   :: Neq
    real, allocatable :: consvar_refval(:)
    real, allocatable :: physvar_refval(:)
    ! Reference values
    real*8    :: refval_length
    real*8    :: refval_mass
    real*8    :: refval_charge
    real*8    :: refval_time
    real*8    :: refval_temperature
    real*8    :: refval_density
    real*8    :: refval_neutral
    real*8    :: refval_speed
    real*8    :: refval_potential
    real*8    :: refval_vorticity
    real*8    :: refval_magfield
    real*8    :: refval_current
    real*8    :: refval_diffusion
    real*8    :: refval_momentum
    real*8    :: refval_specpress
    real*8    :: refval_specenergy
    real*8    :: refval_specenergydens
    ! Adimesional isothermal compressibility coefficient
    real*8    :: compress_coeff

    ! Dimensions used in the simulation
    character(len=20)    :: refval_length_dimensions
    character(len=20)    :: refval_mass_dimensions
    character(len=20)    :: refval_charge_dimensions
    character(len=20)    :: refval_time_dimensions
    character(len=20)    :: refval_temperature_dimensions
    character(len=20)    :: refval_density_dimensions
    character(len=20)    :: refval_neutral_dimensions
    character(len=20)    :: refval_speed_dimensions
    character(len=20)    :: refval_potential_dimensions
    character(len=20)    :: refval_vorticity_dimensions
    character(len=20)    :: refval_magfield_dimensions
    character(len=20)    :: refval_current_dimensions
    character(len=20)    :: refval_diffusion_dimensions
    character(len=20)    :: refval_momentum_dimensions
    character(len=20)    :: refval_specpress_dimensions
    character(len=20)    :: refval_specenergy_dimensions
    character(len=20)    :: refval_specenergydens_dimensions

    ! Physical parameters used in the computation
    real*8    :: a
    real*8    :: Mref
    real*8    :: c1
    real*8    :: c2
    real*8    :: diff_pari
    real*8    :: diff_pare
    real*8    :: diff_n
    real*8    :: diff_u
    real*8    :: diff_e
    real*8    :: diff_ee
  END TYPE Simulationparams_type

  !**********************************************************
  ! Elemental matrices: type to store the elemental matrices
  ! used during the computation
  !**********************************************************
  TYPE :: elmat_type
    real*8, allocatable :: iAqq(:, :, :)
    real*8, allocatable :: Aqu(:, :, :)
    real*8, allocatable :: Aql(:, :, :)
    real*8, allocatable :: Auq(:, :, :)
    real*8, allocatable :: Auu(:, :, :)
    real*8, allocatable :: Aul(:, :, :)
    real*8, allocatable :: Alq(:, :, :)
    real*8, allocatable :: Alu(:, :, :)
    real*8, allocatable :: All(:, :, :)
    real*8, allocatable :: Aql_dir(:, :)
    real*8, allocatable :: Aul_dir(:, :)

    real*8, allocatable :: M(:, :, :)
    real*8, allocatable :: Cv(:, :, :)
    real*8, allocatable :: H(:, :, :)
    real*8, allocatable :: Hdir(:, :)
    real*8, allocatable :: D(:, :, :)
    real*8, allocatable :: E(:, :, :)
    real*8, allocatable :: Edir(:, :)
    real*8, allocatable :: S(:, :)
    real*8, allocatable :: UU(:, :, :)
    real*8, allocatable :: U0(:, :)
    real*8, allocatable :: Hf(:, :, :)
    real*8, allocatable :: Df(:, :, :)
    real*8, allocatable :: Ef(:, :, :)
    real*8, allocatable :: fH(:, :)
    real*8, allocatable :: B(:, :, :)
    real*8, allocatable :: C(:, :, :)
    real*8, allocatable :: Cdir(:, :)
    real*8, allocatable :: P(:, :, :)
    real*8, allocatable :: G(:, :, :)
    real*8, allocatable :: IL(:, :, :)
    real*8, allocatable :: Lf(:, :, :)
    real*8, allocatable :: Qf(:, :, :)
    real*8, allocatable :: LL(:, :, :)
    real*8, allocatable :: L0(:, :)
    ! Limiting rho
    real*8, allocatable :: S_lrho(:, :)
    real*8, allocatable :: P_lrho(:, :, :)
    ! Shock capturing
    real*8, allocatable :: P_sc(:, :, :)
    real*8, allocatable :: Lf_sc(:, :, :)
    real*8, allocatable :: TQ(:, :, :)
    real*8, allocatable :: TQhf(:, :, :)
    real*8, allocatable :: Tfhf(:, :)
    real*8, allocatable :: Tf(:, :)
    real*8, allocatable :: Thdir(:, :)

  END TYPE elmat_type

  !**********************************************************
  ! Timing: structure used to store execution times
  !**********************************************************
  TYPE Timing_type
    integer :: cks1, clock_rate1, cke1, clock_start1, clock_end1
    integer :: cks2, clock_rate2, cke2, clock_start2, clock_end2
    real*8  :: time_start1, time_finish1, tps1, tpe1
    real*8  :: time_start2, time_finish2, tps2, tpe2
    real*8  :: cputpre, cputmap, cputass, cputbcd, cputsol, cputjac, cputglb
    real*8  :: runtpre, runtmap, runtass, runtbcd, runtsol, runtjac, runtglb
    real(8) :: clstime1, clstime2, clstime3, clstime4, clstime5, clstime6
    real(8) :: rlstime1, rlstime2, rlstime3, rlstime4, rlstime5, rlstime6
    real(8) :: cputcom
    real(8) :: runtcom
  END TYPE Timing_type


CONTAINS

  SUBROUTINE set_boundary_flag_names()
    bc_flag_type(1) = 'Tb_Dirichlet'
    bc_flag_type(2) = 'Tb_LEFT'
    bc_flag_type(3) = 'Tb_RIGHT'
    bc_flag_type(4) = 'Tb_UP'
    bc_flag_type(5) = 'Tb_DOWN'
    bc_flag_type(6) = 'Tb_WALL'
    bc_flag_type(7) = 'Tb_LIM'
    bc_flag_type(8) = 'Tb_IN'
    bc_flag_type(9) = 'Tb_OUT'
    bc_flag_type(10) = 'Tb_ULIM'

    bc_flag_name(1) = 'Dirichlet strong form'
    bc_flag_name(2) = 'Dirichlet weak form'
    bc_flag_name(3) = 'Dirichlet weak form using old values of the variables'
    bc_flag_name(4) = 'Transmission boundary conditions'
    bc_flag_name(5) = 'Dirichlet in some equations and Neumann in others '
    bc_flag_name(6) = 'ITER core BC: plasma flux = plasma flux, energy flux imposed '
    bc_flag_name(20) = 'Inlet-Outlet'
    bc_flag_name(30) = 'Neumann homogeneus'
    bc_flag_name(50) = 'Bohm'
    bc_flag_name(51) = 'Bohm for the velocity, Neumann homogeneous for the other variables'
    bc_flag_name(52) = 'Simplified Bohm bc: normal Bohm for the velocity, Grad//T=0 for Ti, Te'
    bc_flag_name(55) = 'Bohm for plasma, pump for neutrals'
    bc_flag_name(56) = 'Bohm for plasma, puff for neutrals'
    bc_flag_name(60) = 'Slip wall'
    bc_flag_name(70) = 'Periodic'

  END SUBROUTINE set_boundary_flag_names

END MODULE types
