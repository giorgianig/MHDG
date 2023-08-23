!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module to load/save files
! in the code
!************************************************************

MODULE in_out
  USE HDF5
  USE HDF5_io_module
  USE globals
  USE printutils
  USE MPI_OMP
  IMPLICIT NONE

CONTAINS

  !********************************
  ! Loads mesh from an hdf5 file
  ! external file
  !********************************
  SUBROUTINE load_mesh(fname)
    USE MPI_OMP
    character(LEN=*) :: fname
    character(len=1000) :: fname_complete
    character(10)  :: str
    character(70)  :: npr, nid
    real*8, parameter::tol = 1e-6
    real*8 :: xmin
    integer :: elemType, ndim, Nnodes, Nelems, Nnodesperelem, Nfaces
    integer :: Nextfaces, Nnodesperface, IERR
    integer(HID_T) :: file_id
#ifdef PARALL
    integer :: ghfa, ghel, i, Nel_glob, Nfa_glob, Ndir_glob, Ngho_glob
#endif
#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%npol .GT. 1) THEN
      write (nid, *) MPIvar%ipol
      write (npr, *) MPIvar%npol
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
#else
    fname_complete = trim(adjustl(fname))//'.h5'
#endif
#else
    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
#endif
    IF (utils%printint > 0) THEN
      print *, 'Loading mesh.'
      print *, '        '
    ENDIF

    CALL HDF5_open(fname_complete, file_id, IERR)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error opening mesh file: ", fname_complete
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, elemType, 'elemType', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: elemType"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, ndim, 'Ndim', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Ndim"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodes, 'Nnodes', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nnodes"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nelems, 'Nelems', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nelems"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperelem, 'Nnodesperelem', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nnodesperelem"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperface, 'Nnodesperface', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nnodesperface"
      STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nextfaces, 'Nextfaces', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nextfaces"
      STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_integer_reading(file_id, Nfaces, 'Nfaces', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading integer: Nfaces"
      STOP
    ENDIF
#endif
    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    ALLOCATE (Mesh%X(Nnodes, ndim))
    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))
#ifdef PARALL
    ALLOCATE (Mesh%ghostFaces(Nfaces))
    ALLOCATE (Mesh%loc2glob_fa(Nfaces))
    ALLOCATE (Mesh%loc2glob_el(Nelems))
    ALLOCATE (Mesh%ghostElems(Nelems))
#endif
    CALL HDF5_array2D_reading_int(file_id, Mesh%T, 'T', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading mesh connectivity T"
      STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%Tb, 'Tb', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading boundary connectivity Tb"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%boundaryFlag, 'boundaryFlag', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading boundaryFlag"
      STOP
    ENDIF
    CALL HDF5_array2D_reading(file_id, Mesh%X, 'X', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading coordinate matrix X"
      STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_fa, 'loc2glob_fa', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading loc2glob_fa"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_el, 'loc2glob_el', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading loc2glob_el"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostFaces, 'ghostFaces', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghostFaces"
      STOP
    ENDIF

    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostElems, 'ghostElems', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghostElems"
      STOP
    ENDIF

    ! Find the number of ghost faces
    ghfa = sum(Mesh%ghostFaces)
    Mesh%nghostfaces = ghfa

    ! Find the number of ghost elements
    ghel = sum(Mesh%ghostElems)
    Mesh%nghostElems = ghel

    ALLOCATE (Mesh%ghostflp(ghfa))
    ALLOCATE (Mesh%ghostpro(ghfa))
    ALLOCATE (Mesh%ghostloc(ghfa))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostflp, 'ghostFlp', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghostFlp"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostLoc, 'ghostLoc', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghostLoc"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostPro, 'ghostPro', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghostPro"
      STOP
    ENDIF

#ifdef TOR3D
    IF (MPIvar%ntor > 1) THEN
      DO i = 1, size(Mesh%ghostPro)
        IF (Mesh%ghostPro(i) .gt. -1) THEN
          Mesh%ghostPro(i) = Mesh%ghostPro(i) + (MPIvar%itor - 1)*MPIvar%npol
        ENDIF
      END DO
    ENDIF
    ALLOCATE (Mesh%ghelspro(ghel))
    ALLOCATE (Mesh%ghelsloc(ghel))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsLoc, 'ghelsLoc', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghelsLoc"
      STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsPro, 'ghelsPro', ierr)
    IF (IERR .ne. 0) THEN
      WRITE (6, *) "Error reading ghelsPro"
      STOP
    ENDIF
    IF (MPIvar%ntor .gt. 1) then
      DO i = 1, size(Mesh%ghelspro)
        IF (Mesh%ghelspro(i) .gt. -1) THEN
          Mesh%ghelsPro(i) = Mesh%ghelsPro(i) + (MPIvar%itor - 1)*MPIvar%npol
        END IF
      END DO
    END IF
#endif
#endif
    CALL HDF5_close(file_id)

    !************************************************************************
    !   CONFIRMATION MESSAGE FOR THE USER
    !************************************************************************
#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    WRITE (6, *) "Process: ", MPIvar%glob_id, "-- readed mesh file: ", trim(adjustl(fname_complete))
#else
    WRITE (6, *) "Readed mesh file: ", trim(adjustl(fname_complete))
#endif

#ifdef PARALL
    CALL MPI_ALLREDUCE(maxval(Mesh%loc2glob_el), Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(maxval(Mesh%loc2glob_fa), Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%ndir, Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%nghostfaces, Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    Mesh%Nel_glob = Nel_glob
    Mesh%Nfa_glob = Nfa_glob
    Mesh%Ndir_glob = Ndir_glob
    Mesh%Ngho_glob = Ngho_glob
#endif
    Mesh%Ndim = ndim
    Mesh%Nnodes = Nnodes
    Mesh%Nelems = Nelems
    Mesh%Nnodesperelem = Nnodesperelem
    Mesh%Nnodesperface = Nnodesperface
    Mesh%elemType = elemType
    Mesh%Nextfaces = Nextfaces

    xmin = minval(Mesh%X(:,1))
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
    ! Apply shift if axisymmetric case
    IF ((switch%axisym .and. switch%testcase .ge. 60 .and. switch%testcase .lt. 80) .or. (switch%axisym .and. xmin < tol)) THEN
      IF (MPIvar%glob_id .eq. 0) THEN
        WRITE (6, *) "*** Applying translation in axisymmetric case!"
      ENDIF
      Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
    END IF

    ! Apply length scale
    Mesh%X = Mesh%X/phys%lscale

    Mesh%xmax = maxval(Mesh%X(:, 1))
    Mesh%xmin = minval(Mesh%X(:, 1))
    Mesh%ymax = maxval(Mesh%X(:, 2))
    Mesh%ymin = minval(Mesh%X(:, 2))

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    IF (utils%printint > 0) then
      IF (MPIvar%glob_id .eq. 0) THEN
        IF (elemType == 0) then
          WRITE (str, '(A)') 'triangles'
        ELSEIF (elemType == 1) then
          WRITE (str, '(A)') 'quads'
        ELSEIF (elemType == 2) then
          WRITE (str, '(A)') 'thetra'
        ELSEIF (elemType == 3) then
          WRITE (str, '(A)') 'hexa'
        END IF
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*                    MESH                       *'
        WRITE (6, *) '*************************************************'
        WRITE (6, '(A,I18)') ' Number of dimensions:         ', ndim
        WRITE (6, '(A,A34)') ' Element type: ', trim(str)
        WRITE (6, '(A,I18)') ' Number of elements:           ', Nelems
        WRITE (6, '(A,I18)') ' Number of nodes:              ', Nnodes
        WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Nnodesperelem
        WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Nnodesperface
        WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Nextfaces
        WRITE (6, *) ' '
        WRITE (6, *) ' '
        IF (utils%printint > 1) THEN
          WRITE (6, *) "Connectivity matrix T:"
          CALL displayMatrixInt(Mesh%T)
          WRITE (6, *) "Boundary connectivity matrix Tb:"
          CALL displayMatrixInt(Mesh%Tb)
        END IF
      ENDIF
    ENDIF

  END SUBROUTINE load_mesh

  !**********************************************************************
  ! Save solution in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_solution(fname)
    USE globals
    implicit none

    character(LEN=*) :: fname
    character(70)  :: npr, nid
#ifdef TOR3D
    character(70)  :: nip, nit, ngd
#endif
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

#ifdef TOR3D
    IF (MPIvar%glob_size .GT. 1) THEN
      write (nip, *) MPIvar%ipol
      write (nit, *) MPIvar%itor
      write (ngd, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_ip'//trim(adjustl(nip))//'_it'//trim(adjustl(nit))//'_np'//trim(adjustl(ngd))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
#else
    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
#endif
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_array1D_saving(file_id, sol%u, size(sol%u), 'u')
    call HDF5_array1D_saving(file_id, sol%u_tilde, size(sol%u_tilde), 'u_tilde')
    if (phys%Neq .ge. 5) then
       if (switch%saveTau) then
          call HDF5_array1D_saving(file_id, phys%diff_nn_Vol, size(phys%diff_nn_Vol), 'DnnVol')
          call HDF5_array1D_saving(file_id, phys%diff_nn_Fac, size(phys%diff_nn_Fac), 'DnnFac')
          call HDF5_array1D_saving(file_id, phys%diff_nn_Bou, size(phys%diff_nn_Bou), 'DnnBou')
          call HDF5_array2D_saving(file_id, phys%v_nn_Vol, size(phys%v_nn_Vol, 1), size(phys%v_nn_Vol, 2), 'VnnVol')
          call HDF5_array2D_saving(file_id, phys%v_nn_Fac, size(phys%v_nn_Fac, 1), size(phys%v_nn_Fac, 2), 'VnnFac')
          call HDF5_array2D_saving(file_id, phys%v_nn_Bou, size(phys%v_nn_Bou, 1), size(phys%v_nn_Bou, 2), 'VnnBou')
          call HDF5_array2D_saving(file_id, Mesh%Xg, size(Mesh%Xg, 1), size(Mesh%Xg, 2), 'Xg')
          call HDF5_array2D_saving(file_id, Mesh%Xgf, size(Mesh%Xgf, 1), size(Mesh%Xgf, 2), 'Xgf')
          call HDF5_array2D_saving(file_id, Mesh%Xgb, size(Mesh%Xgb, 1), size(Mesh%Xgb, 2), 'Xgb')
       endif
    endif   
    !      call HDF5_array1D_saving(file_id,sol%tres,sol%Nt,'tres')
    !      call HDF5_array1D_saving(file_id,sol%time,sol%Nt,'time')
    if (switch%steady .or. switch%psdtime) then
      call HDF5_integer_saving(file_id,0,'it')
    else
      call HDF5_integer_saving(file_id,time%it,'it')
    endif
    !!!      call HDF5_integer_saving(file_id,sol%Nt,'Nt')
    call HDF5_array1D_saving(file_id, sol%q, size(sol%q), 'q')
    ! Save magnetic field
    call HDF5_array2D_saving(file_id, phys%B, size(phys%B, 1), size(phys%B, 2), 'magnetic_field')
    ! Save normalized psi
    call HDF5_array1D_saving(file_id, phys%magnetic_psi, size(phys%magnetic_psi), 'magnetic_psi')
    ! Save toroidal current
    call HDF5_array1D_saving(file_id, phys%Jtor, size(phys%Jtor), 'Jtor')
    ! Save magnetic perturbation and related fields
    if ((switch%rmp).or.(switch%ripple)) then
      call HDF5_array2D_saving(file_id, phys%Bperturb, size(phys%Bperturb, 1), size(phys%Bperturb, 2), 'magnetic_perturbation')
    endif
    if (switch%rmp) then
      call HDF5_array3D_saving(file_id, magn%coils_rmp, size(magn%coils_rmp, 1), size(magn%coils_rmp, 2), size(magn%coils_rmp, 3), 'coils_rmp')
    endif
    if (switch%ripple) then
      call HDF5_array2D_saving(file_id, magn%coils_ripple, size(magn%coils_ripple, 1), size(magn%coils_ripple, 2), 'coils_ripple')
    endif
    ! Save boundary structure
    call HDF5_array2D_saving_int(file_id, Mesh%extfaces, size(Mesh%extfaces, 1), size(Mesh%extfaces, 2), 'exterior_faces')
    call HDF5_array1D_saving_int(file_id, Mesh%boundaryFlag, size(Mesh%boundaryFlag, 1), 'boundary_flags')
    ! Save simulation parameters
    call save_simulation_parameters()
    IF (switch%shockcp .eq. 3) THEN
      call HDF5_array2D_saving(file_id, Mesh%scdiff_nodes, size(Mesh%scdiff_nodes, 1), size(Mesh%scdiff_nodes, 2), 'scdiff_nodes')
    END IF
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    IF (MPIvar%glob_id .eq. 0) THEN
      print *, 'Output written to file ', trim(adjustl(fname_complete))
      print *, '        '
    END IF

  contains

    !**********************************************************************
    ! Save simulation parameters
    !**********************************************************************
    subroutine save_simulation_parameters()
      integer(HID_T) :: group_id1, group_id2, group_id3

      ! Create simulation parameters group
      CALL HDF5_group_create('simulation_parameters', file_id, group_id1, ierr)
      ! Save model definition
      call HDF5_string_saving(group_id1, simpar%model, 'model')
      ! Save Ndim and Neq
      call HDF5_integer_saving(group_id1, simpar%Ndim, 'Ndim')
      call HDF5_integer_saving(group_id1, simpar%Neq, 'Neq')
      ! Create adimensionalization subgroup
      CALL HDF5_group_create('adimensionalization', group_id1, group_id2, ierr)
      ! Save adimensionalization
      call HDF5_string_saving(group_id2, simpar%refval_time_dimensions, 'time_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_mass_dimensions, 'mass_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_length_dimensions, 'length_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_temperature_dimensions, 'temperature_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_density_dimensions, 'density_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_neutral_dimensions, 'density_neutral_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_speed_dimensions, 'speed_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_potential_dimensions, 'potential_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_vorticity_dimensions, 'vorticity_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_magfield_dimensions, 'magfield_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_current_dimensions, 'current_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_diffusion_dimensions, 'diffusion_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_momentum_dimensions, 'momentum_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_specpress_dimensions, 'specific_pressure_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_specenergy_dimensions, 'specific_energy_scale_dimensions')
      call HDF5_string_saving(group_id2, simpar%refval_specenergydens_dimensions, 'specific_energy_density_scale_dimensions')
      call HDF5_real_saving(group_id2, simpar%refval_length, 'length_scale')
      call HDF5_real_saving(group_id2, simpar%refval_time, 'time_scale')
      call HDF5_real_saving(group_id2, simpar%refval_mass, 'mass_scale')
      call HDF5_real_saving(group_id2, simpar%refval_temperature, 'temperature_scale')
      call HDF5_real_saving(group_id2, simpar%refval_density, 'density_scale')
      call HDF5_real_saving(group_id2, simpar%refval_neutral, 'neutral_scale')
      call HDF5_real_saving(group_id2, simpar%refval_speed, 'speed_scale')
      call HDF5_real_saving(group_id2, simpar%refval_potential, 'potential_scale')
      call HDF5_real_saving(group_id2, simpar%refval_vorticity, 'vorticity_scale')
      call HDF5_real_saving(group_id2, simpar%refval_magfield, 'magfield_scale')
      call HDF5_real_saving(group_id2, simpar%refval_current, 'current_scale')
      call HDF5_real_saving(group_id2, simpar%refval_diffusion, 'diffusion_scale')
      call HDF5_real_saving(group_id2, simpar%refval_momentum, 'momentum_scale')
      call HDF5_real_saving(group_id2, simpar%refval_specpress, 'specific_pressure_scale')
      call HDF5_real_saving(group_id2, simpar%refval_specenergy, 'specific_energy_scale')
      call HDF5_real_saving(group_id2, simpar%refval_specenergydens, 'specific_energy_density_scale')
      call HDF5_array1d_saving(group_id2, simpar%physvar_refval, phys%npv, 'reference_values_physical_variables')
      call HDF5_array1d_saving(group_id2, simpar%consvar_refval, phys%Neq, 'reference_values_conservative_variables')
      call HDF5_group_close(group_id2, ierr)
      ! Close group adimensionalization

      ! Create physics parameters group
      CALL HDF5_group_create('physics', group_id1, group_id2, ierr)
      call HDF5_string_array1D_saving(group_id2, phys%phyVarNam, 'physical_variable_names')
      call HDF5_string_array1D_saving(group_id2, phys%conVarNam, 'conservative_variable_names')
      call HDF5_real_saving(group_id2, phys%a, 'a')
      call HDF5_real_saving(group_id2, phys%Mref, 'Mref')
      call HDF5_real_saving(group_id2, phys%c1, 'c1')
      call HDF5_real_saving(group_id2, phys%c2, 'c2')
      call HDF5_real_saving(group_id2, phys%diff_pari, 'diff_pari')
      call HDF5_real_saving(group_id2, phys%diff_pare, 'diff_pare')
      call HDF5_real_saving(group_id2, phys%etapar, 'eta_parallel')
      call HDF5_real_saving(group_id2, phys%diff_n, 'diff_n')
      call HDF5_real_saving(group_id2, phys%diff_u, 'diff_u')
      call HDF5_real_saving(group_id2, phys%diff_e, 'diff_e')
      call HDF5_real_saving(group_id2, phys%diff_ee, 'diff_ee')
      call HDF5_real_saving(group_id2, phys%diff_vort, 'diff_vort')
      call HDF5_real_saving(group_id2, phys%diff_pot, 'diff_pot')
      call HDF5_real_saving(group_id2, phys%diff_nn, 'diff_nn')
	     call HDF5_real_saving(group_id2, phys%Re, 'recycling')
	     call HDF5_real_saving(group_id2, phys%puff, 'puff')
	     IF (switch%ME) THEN
	        call HDF5_array1d_saving(group_id2, phys%puff_exp, time%nts, 'puff_exp')
	     END IF
      call HDF5_real_saving(group_id2, phys%tie, 'tau_ie')
      call HDF5_real_saving(group_id2, phys%dfcoef, 'dfcoef')
      call HDF5_real_saving(group_id2, phys%dexbcoef, 'dexbcoef')
      call HDF5_real_saving(group_id2, phys%bohmth, 'bohmth')
      call HDF5_real_saving(group_id2, phys%epn, 'epn')
      call HDF5_real_saving(group_id2, phys%Gmbohm, 'Gmbohm')
      call HDF5_real_saving(group_id2, phys%Gmbohme, 'Gmbohme')
      call HDF5_real_saving(group_id2, phys%Potfloat, 'Potfloat')
      call HDF5_array1d_saving_int(group_id2, phys%bcflags, 10, 'boundary_flags')
      call HDF5_group_close(group_id2, ierr)

      ! Create switches parameters group
      CALL HDF5_group_create('switches', group_id1, group_id2, ierr)
      call HDF5_logical_saving(group_id2, switch%driftdia, 'diamagnetic_drift')
      call HDF5_logical_saving(group_id2, switch%driftexb, 'ExB_drift')
      call HDF5_logical_saving(group_id2, switch%steady, 'steady')
      call HDF5_integer_saving(group_id2, switch%testcase, 'testcase')
      call HDF5_logical_saving(group_id2, switch%ohmicsrc, 'ohmicsrc')
      call HDF5_logical_saving(group_id2, switch%ohmicsrc, 'ME')
      call HDF5_logical_saving(group_id2, switch%rmp, 'RMP')
      call HDF5_logical_saving(group_id2, switch%ripple, 'Ripple')
      call HDF5_logical_saving(group_id2, switch%psdtime, 'psdtime')
      call HDF5_real_saving(group_id2, switch%diffred, 'diffred')
      call HDF5_real_saving(group_id2, switch%diffmin, 'diffmin')
      call HDF5_integer_saving(group_id2, switch%shockcp, 'shockcp')
      call HDF5_integer_saving(group_id2, switch%limrho, 'limrho')
      call HDF5_integer_saving(group_id2, switch%difcor, 'difcor')
      call HDF5_integer_saving(group_id2, switch%thresh, 'thresh')
      call HDF5_logical_saving(group_id2, switch%filter, 'filter')
      call HDF5_logical_saving(group_id2, switch%decoup, 'decoup')
      call HDF5_logical_saving(group_id2, switch%ckeramp, 'ckeramp')
      call HDF5_logical_saving(group_id2, switch%saveNR, 'saveNR')
      call HDF5_logical_saving(group_id2, switch%saveTau, 'saveTau')
      call HDF5_logical_saving(group_id2, switch%fixdPotLim, 'fixdPotLim')
      call HDF5_logical_saving(group_id2, switch%dirivortcore, 'dirivortcore')
      call HDF5_logical_saving(group_id2, switch%dirivortlim, 'dirivortlim')
      call HDF5_logical_saving(group_id2, switch%convvort, 'convvort')
      call HDF5_logical_saving(group_id2, switch%logrho, 'logrho')
      call HDF5_group_close(group_id2, ierr)

      ! Create numerics parameters group
      CALL HDF5_group_create('numerics', group_id1, group_id2, ierr)
      call HDF5_integer_saving(group_id2, numer%nrp, 'Max_number_of_NR_iterations')
      call HDF5_real_saving(group_id2, numer%tnr, 'NR_convergence_criterium')
      call HDF5_real_saving(group_id2, numer%ttm, 'Time_convergence_criterium')
      call HDF5_real_saving(group_id2, numer%div, 'Divergence_criterium')
      call HDF5_array1d_saving(group_id2, numer%tau, 4, 'Stabilization_parameter')
      call HDF5_real_saving(group_id2, numer%sc_coe, 'Shock_capturing_parameter')
      call HDF5_real_saving(group_id2, numer%sc_sen, 'Shock_capturing_sensibility')
      call HDF5_real_saving(group_id2, numer%minrho, 'Value_of_rho_to_start_applying_limiting')
      call HDF5_real_saving(group_id2, numer%so_coe, 'Source_coefficient_for_limiting_rho')
      call HDF5_real_saving(group_id2, numer%df_coe, 'Diffusion_coefficient_for_limiting_rho')
      call HDF5_real_saving(group_id2, numer%dc_coe, 'Diffusion_coefficient_in_corners')
      call HDF5_real_saving(group_id2, numer%thr, 'Threshold_to_limit_rho')
      call HDF5_real_saving(group_id2, numer%thrpre, 'Threshold_to_limit_pressure')
      call HDF5_integer_saving(group_id2, numer%stab, 'Stabilization_type')
      call HDF5_real_saving(group_id2, numer%dumpnr, 'dumping_factor_for_Newton_Raphson')
      call HDF5_integer_saving(group_id2, numer%ntor, 'Number_of_elements_in_the_toroidal_direction')
      call HDF5_integer_saving(group_id2, numer%ptor, 'Polynomial_degree_in_the_toroidal_direction')
      call HDF5_real_saving(group_id2, numer%tmax, 'Max_extention_in_the_toroidal_direction')
      call HDF5_integer_saving(group_id2, numer%npartor, 'Number_of_MPI_divisions_in_the_toroidal_direction')
      call HDF5_real_saving(group_id2, numer%exbdump, 'Dumping_for_ExB_drift')
      call HDF5_group_close(group_id2, ierr)

      ! Create time parameters group
      CALL HDF5_group_create('time', group_id1, group_id2, ierr)
      call HDF5_real_saving(group_id2, time%dt0, 'Initial_time_step')
      call HDF5_real_saving(group_id2, time%dt, 'Current_time_step')
      call HDF5_real_saving(group_id2, time%tfi, 'Final_time')
      call HDF5_integer_saving(group_id2, time%it, 'Current_time_step_number')
      call HDF5_integer_saving(group_id2, time%ik, 'Current_pseudo_time_step_number')
      call HDF5_integer_saving(group_id2, time%nts, 'Number_of_time_steps')
      call HDF5_integer_saving(group_id2, time%tis, 'Time_integration_scheme')
      call HDF5_real_saving(group_id2, time%t, 'Current_time')
      call HDF5_group_close(group_id2, ierr)

      ! Create geometry parameters group
      CALL HDF5_group_create('geometry', group_id1, group_id2, ierr)
      call HDF5_real_saving(group_id2, geom%R0, 'Major_radius')
      call HDF5_real_saving(group_id2, geom%q, 'Safety_factor')
      call HDF5_group_close(group_id2, ierr)

      ! Create time parameters group
      CALL HDF5_group_create('magnetic', group_id1, group_id2, ierr)
      call HDF5_real_saving(group_id2, magn%amp_rmp, 'Amplitude_RMP')
      call HDF5_integer_saving(group_id2, magn%nbCoils_rmp, 'Number_coils_RMP')
      call HDF5_integer_saving(group_id2, magn%parite, 'Parity_RMP')
      call HDF5_integer_saving(group_id2, magn%nbRow, 'number_rows_RMP')
      call HDF5_real_saving(group_id2, magn%amp_ripple, 'Amplitude_Ripple')
      call HDF5_integer_saving(group_id2, magn%nbCoils_ripple, 'Number_coils_Ripple')
      call HDF5_real_saving(group_id2, magn%triang, 'Triangularity')
      call HDF5_real_saving(group_id2, magn%ellip, 'Ellipticity')
      call HDF5_group_close(group_id2, ierr)

      call HDF5_group_close(group_id1, ierr)


    end subroutine save_simulation_parameters

  end subroutine HDF5_save_solution

  !**********************************************************************
  ! Load solution in HDF5 file format
  !**********************************************************************
  subroutine HDF5_load_solution(fname)
    USE globals
    USE LinearAlgebra, ONLY: tensorsumint, colint,col
    implicit none

    character(LEN=1000) :: fname
    character(len=100), pointer :: mod_ptr
    character(len=100), target :: model_string

    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id, group_id,group_id2
    real, pointer :: oldtres(:), oldtime(:)
    real, pointer:: u_aux(:), u_tilde_aux(:), q_aux(:)
#ifdef TOR3D
    character(70)  :: nip, nit, ngd
    integer :: ntorloc
#endif
    integer :: Neq, Ndim, N2d, Nel, Np1d, Np2d, Np, Nfl, Nfp, Nfg, Nf, sizeutilde, sizeu
    integer :: iel, ifa, iface, Fi, itor, dd, delta, i, j, it, iel3
    real, allocatable    :: u2D(:), q2D(:), ut2D(:), q2D_add(:, :),uaux(:,:),utaux(:,:),qaux(:,:)
    integer, allocatable :: indu2D(:), indq2D(:), indu3D(:), indq3D(:), indufp(:), induf2D(:), induf3D(:), indul(:), indql(:), indutl(:)
    integer, allocatable :: ind_q_add(:)
    integer              :: logrho_ptr
    real*8               :: t
    Neq = phys%Neq
    mod_ptr => model_string
#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%ntor .gt. 1) THEN
      ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
      ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    Ndim = 3                             ! N. of dimensions
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                      ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
    sizeu = Neq*Nel*Np                    ! Size of u
#ifdef PARALL
    if (MPIvar%ntor .gt. 1) then
      sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d) + Neq*Np2d*N2d! Size of utilde
    else
      sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
    endif
#else
    sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
#endif
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
    sizeu = Neq*Nel*Np
    sizeutilde = Neq*Mesh%Nfaces*Mesh%Nnodesperface
#endif

    ALLOCATE (sol%u(sizeu))
    ALLOCATE (sol%u_tilde(sizeutilde))
    ALLOCATE (sol%q(sizeu*Ndim))

#ifdef TOR3D
    !*************************************
    !              3D case
    !*************************************

    IF (fname(1:5) == 'Sol3D') THEN
      ! Initialization with a 3D solution
      WRITE (6, *) "3D initial solution"

      IF (MPIvar%glob_size .GT. 1) THEN
        write (nip, *) MPIvar%ipol
        write (nit, *) MPIvar%itor
        write (ngd, *) MPIvar%glob_size
        fname_complete = trim(adjustl(fname))//'_ip'//trim(adjustl(nip))//'_it'//trim(adjustl(nit))//'_np'//trim(adjustl(ngd))//'.h5'
      ELSE
        fname_complete = trim(adjustl(fname))//'.h5'
      END IF

      CALL HDF5_open(fname_complete, file_id, IERR)
      CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
      CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
      CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
      CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
      CALL HDF5_group_close(group_id2, ierr)
#else
      logrho_ptr = switch%logrho
#endif
      CALL HDF5_group_close(group_id, ierr)
      ! Check if the readed solution corresponds to the right model
      IF (simpar%model .ne. model_string) THEN
        WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
        STOP
      ENDIF

      CALL HDF5_array1D_reading(file_id, sol%u, 'u')
      CALL HDF5_array1D_reading(file_id, sol%u_tilde, 'u_tilde')
      CALL HDF5_integer_reading(file_id,time%it,'it')
      !      if (.not.switch%steady) then
      !         CALL HDF5_integer_reading(file_id,sol%Nt,'Nt')
      !         allocate(oldtres(sol%Nt))
      !         allocate(oldtime(sol%Nt))
      !         call HDF5_array1D_reading(file_id,oldtres,'tres')
      !         call HDF5_array1D_reading(file_id,oldtime,'time')
      !         if (allocated(sol%tres)) then
      !            deallocate(sol%tres)
      !         endif
      !         if (allocated(sol%time)) then
      !            deallocate(sol%time)
      !         endif
      !         allocate(sol%tres(time%nts+sol%Nt))
      !         allocate(sol%time(time%nts+sol%Nt))
      !         sol%tres = 0.
      !         sol%time = 0.
      !         sol%tres(1:sol%Nt) = oldtres
      !         sol%time(1:sol%Nt) = oldtime
      !         time%t = sol%time(sol%Nt)
      !         deallocate(oldtres,oldtime)
      !      endif
      CALL HDF5_array1D_reading(file_id, sol%q, 'q')
      CALL HDF5_close(file_id)
    ELSEIF (fname(1:5) == 'Sol2D') THEN
      ! Initialization with a 2D solution
      WRITE (6, *) "2D initial solution: propagating in the torus..."

      IF (MPIvar%glob_size .GT. 1) THEN
        write (nid, *) MPIvar%ipol
        write (npr, *) MPIvar%npol
        fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
      ELSE
        fname_complete = trim(adjustl(fname))//'.h5'
      END IF

      CALL HDF5_open(fname_complete, file_id, IERR)
      CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
      CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
      CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
      CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
      CALL HDF5_group_close(group_id2, ierr)
#else
      logrho_ptr = switch%logrho
#endif
      CALL HDF5_group_close(group_id, ierr)
      ! Check if the readed solution corresponds to the right model
      IF (simpar%model .ne. model_string) THEN
        WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
        STOP
      ENDIF

      ALLOCATE (u_aux(Neq*N2d*Np))
      ALLOCATE (u_tilde_aux(Neq*Mesh%Nfaces*Mesh%Nnodesperface))
      ALLOCATE (u2D(Np2D*Neq))
      ALLOCATE (ut2D(refElPol%Nnodes1D*Neq))
      ALLOCATE (indu2D(Np2D*Neq))
      ALLOCATE (indu3D(Np*Neq))
      ALLOCATE (indufp(Np2D*Neq))
      ALLOCATE (induf2D(refElPol%Nnodes1D*Neq))
      ALLOCATE (induf3D(Nfl*Neq))
      ALLOCATE (indul(Np2D*Neq))
      ALLOCATE (indutl(refElPol%Nnodes1D*Neq))
      ALLOCATE (q_aux(Neq*N2d*Np*(Ndim - 1)))
      ALLOCATE (q2D(Np2D*Neq*(Ndim - 1)))
      ALLOCATE (q2D_add(Np2D, Ndim*Neq))
      ALLOCATE (indq2D(Np2D*Neq*(Ndim - 1)))
      ALLOCATE (indq3D(Np*Neq*Ndim))
      ALLOCATE (indql(Np2D*Neq*Ndim))
      ALLOCATE (ind_q_add(Neq*(Ndim - 1)))
      q2D_add = 0.
      ind_q_add = colint(tensorsumint((/(j, j=1, (ndim - 1))/), ndim*(/(i, i=0, (Neq - 1))/)))

      CALL HDF5_open(fname_complete, file_id, IERR)
      CALL HDF5_array1D_reading(file_id, u_aux, 'u')
      CALL HDF5_array1D_reading(file_id, u_tilde_aux, 'u_tilde')
      CALL HDF5_array1D_reading(file_id, q_aux, 'q')
      CALL HDF5_close(file_id)

#ifdef PARALL
      WRITE (6, *) "Process: ", MPIvar%glob_id, "-- readed solution file: ", trim(adjustl(fname_complete))
#else
      WRITE (6, *) "Readed solution file: ", trim(adjustl(fname_complete))
#endif

      DO iel = 1, N2D
        indu2D = (iel - 1)*Np2d*Neq + (/(i, i=1, Np2d*Neq)/)
        u2D = u_aux(indu2D)
        indq2D = (iel - 1)*Np2d*Neq*(Ndim - 1) + (/(i, i=1, Np2d*Neq*(Ndim - 1))/)
        q2D = q_aux(indq2D)
        q2d_add(:, ind_q_add) = transpose(reshape(q2D, [(Ndim - 1)*Neq, Np2d]))

        DO itor = 1, ntorloc
          iel3 = (itor - 1)*N2d+iel
          indu3D = (iel3 - 1)*Np*Neq + (/(i, i=1, Np*Neq)/)
          indq3D = (iel3 - 1)*Np*Neq*Ndim + (/(i, i=1, Np*Neq*Ndim)/)
          dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (iel - 1)*Np2D*Neq
          indufp = dd + (/(i, i=1, Np2D*Neq)/)
          sol%u_tilde(indufp) = u2d
          DO it = 1, Np1d
            indul = (it - 1)*Np2D*Neq + (/(i, i=1, Np2D*Neq)/)
            sol%u(indu3D(indul)) = u2D
            indql = (it - 1)*Np2D*Neq*Ndim + (/(i, i=1, Np2D*Neq*Ndim)/)
            sol%q(indq3D(indql)) = reshape(transpose(q2d_add), (/Np2d*Neq*Ndim/))
          END DO
        END DO
      END DO

      DO iface = 1, Mesh%Nintfaces
        Fi = iface
        induf2D = (Fi - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
        ut2d = u_tilde_aux(induf2D)
        DO itor = 1, ntorloc
          dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          induf3D = dd + (/(i, i=1, Nfl*Neq)/)
          DO it = 1, Np1d
            indutl = (it - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
            sol%u_tilde(induf3D(indutl)) = ut2d
          END DO
        END DO
      END DO

      DO iface = 1, Mesh%Nextfaces
        iel = Mesh%extfaces(iface, 1)
        ifa = Mesh%extfaces(iface, 2)
        IF (Mesh%Fdir(iel, ifa)) CYCLE
        Fi = iface + Mesh%Nintfaces
        induf2D = (Fi - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
        ut2d = u_tilde_aux(induf2D)
        DO itor = 1, ntorloc
          dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
          induf3D = dd + (/(i, i=1, Nfl*Neq)/)
          DO it = 1, Np1d
            indutl = (it - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
            sol%u_tilde(induf3D(indutl)) = ut2d
          END DO
        END DO
      END DO

#ifdef PARALL
      ! Add solution on toroidal ghost faces
      IF (MPIvar%ntor .gt. 1) THEN

        DO iel = 1, N2D
          indu2D = (iel - 1)*Np2d*Neq + (/(i, i=1, Np2d*Neq)/)
          u2D = u_aux(indu2D)
          indq2D = (iel - 1)*Np2d*Neq*Ndim + (/(i, i=1, Np2d*Neq*Ndim)/)
          q2D = q_aux(indq2D)
          DO itor = 1, ntorloc
            dd = ntorloc*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (iel - 1)*Np2D*Neq
            indufp = dd + (/(i, i=1, Np2D*Neq)/)
            sol%u_tilde(indufp) = u2d
          END DO
        END DO
      ENDIF
#endif

      WRITE (6, *) "Done!"

      DEALLOCATE (u_aux, u_tilde_aux, u2D, ut2D, indu2D, indu3D, indufp, induf2D, induf3D, indul, indutl)
      DEALLOCATE (q_aux, q2D, q2D_add, indq2D, indq3D, indql, ind_q_add)
    END IF

#else
    !*************************************
    !              2D case
    !*************************************

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    CALL HDF5_open(fname_complete, file_id, IERR)
    CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
    CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
    CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
    CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
    CALL HDF5_group_close(group_id2, ierr)
#else  
    logrho_ptr = switch%logrho
#endif
    IF (switch%ME) THEN
       CALL HDF5_group_open(group_id, 'time', group_id2, ierr)
       CALL HDF5_integer_reading(group_id2, it, 'Current_time_step_number')
       CALL HDF5_real_reading(group_id2, t, 'Current_time')
       IF (it .gt. 1) THEN
          time%it = it
          time%ik = it
          time%t = t
          sol%Nt = it
          sol%time(it) = t
       ELSE
       END IF
       CALL HDF5_group_close(group_id2, ierr)
       IF (time%it .ne. 0) THEN
          CALL HDF5_group_open(group_id, 'physics', group_id2, ierr)
          CALL HDF5_array1D_reading(group_id2, phys%puff_exp, 'puff_exp')
          CALL HDF5_group_close(group_id2, ierr)
       END IF
    END IF
    CALL HDF5_group_close(group_id, ierr)


    ! Check if the readed solution corresponds to the right model
    IF (simpar%model .ne. model_string) THEN
      WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
      STOP
    ENDIF
    CALL HDF5_array1D_reading(file_id, sol%u, 'u')
    CALL HDF5_array1D_reading(file_id, sol%u_tilde, 'u_tilde')
    CALL HDF5_array1D_reading(file_id, sol%q, 'q')
    CALL HDF5_close(file_id)
#endif

    if (switch%logrho .and. logrho_ptr.eq.0 ) then
      write(6,*) "Readed solution without logrho but switch logrho set to true: "
      write(6,*) "Computing log of density"
      allocate(uaux(size(sol%u)/phys%neq,phys%neq))
      allocate(utaux(size(sol%u_tilde)/phys%neq,phys%neq))
      allocate(qaux(size(sol%q)/phys%neq/Ndim,phys%neq*Ndim))
      uaux = transpose(reshape(sol%u,[phys%neq,size(sol%u)/phys%neq]))
      utaux = transpose(reshape(sol%u_tilde,[phys%neq,size(sol%u_tilde)/phys%neq]))
      qaux = transpose(reshape(sol%q,[phys%neq*Ndim,size(sol%q)/phys%neq/Ndim]))
      qaux(:,1) = qaux(:,1)/uaux(:,1)
      qaux(:,2) = qaux(:,2)/uaux(:,1)
      uaux(:,1) = log(uaux(:,1))
      utaux(:,1) = log(utaux(:,1))
      sol%u = col(transpose(uaux))
      sol%u_tilde = col(transpose(utaux))
      sol%q = col(transpose(qaux))
      deallocate(uaux,utaux,qaux)
      elseif (.not.switch%logrho .and. logrho_ptr.eq.1 ) then
      write(6,*) "Readed solution with logrho but switch logrho set to false: "
      write(6,*) "Computing exp of density"
      allocate(uaux(size(sol%u)/phys%neq,phys%neq))
      allocate(utaux(size(sol%u_tilde)/phys%neq,phys%neq))
      allocate(qaux(size(sol%q)/phys%neq/Ndim,phys%neq*Ndim))
      uaux = transpose(reshape(sol%u,[phys%neq,size(sol%u)/phys%neq]))
      utaux = transpose(reshape(sol%u_tilde,[phys%neq,size(sol%u_tilde)/phys%neq]))
      qaux = transpose(reshape(sol%q,[phys%neq*Ndim,size(sol%q)/phys%neq/Ndim]))
      uaux(:,1) = exp(uaux(:,1))
      utaux(:,1) = exp(utaux(:,1))
      qaux(:,1) = qaux(:,1)*uaux(:,1)
      qaux(:,2) = qaux(:,2)*uaux(:,1)
      sol%u = col(transpose(uaux))
      sol%u_tilde = col(transpose(utaux))
      sol%q = col(transpose(qaux))
      deallocate(uaux,utaux,qaux)
    endif


    ! Message to confirm succesful reading of file
    IF (MPIvar%glob_id .eq. 0) THEN
      print *, 'Solution read from file ', trim(adjustl(fname_complete))
      print *, '        '
    END IF

  end subroutine HDF5_load_solution

  !**********************************************************************
  ! Save HDG matrix (CSR) in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_CSR_matrix(fname)
    USE globals
    implicit none

    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_integer_saving(file_id, MatK%n, 'n')
    call HDF5_integer_saving(file_id, MatK%nnz, 'nnz')
    call HDF5_array1D_saving_int(file_id, MatK%cols, MatK%nnz, 'cols')
    call HDF5_array1D_saving_int(file_id, MatK%rowptr, MatK%n + 1, 'rowptr')
    call HDF5_array1D_saving_int(file_id, MatK%loc2glob, MatK%n, 'loc2glob')
    call HDF5_array1D_saving(file_id, MatK%vals, MatK%nnz, 'vals')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_CSR_matrix

  !**********************************************************************
  ! Save HDG vector (CSR) in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_CSR_vector(fname)
    USE globals
    implicit none

    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_integer_saving(file_id, rhs%n, 'n')
    call HDF5_array1D_saving_int(file_id, rhs%loc2glob, rhs%n, 'loc2glob')
    call HDF5_array1D_saving(file_id, rhs%vals, rhs%n, 'vals')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_CSR_vector

  !**********************************************************************
  ! Save 3D array in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_array(Arr, fname)
    USE globals
    implicit none

    REAL, DIMENSION(:, :, :), INTENT(IN) :: Arr
    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_array3D_saving(file_id, Arr, size(Arr, 1), size(Arr, 2), size(Arr, 3), 'array')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_array

  !**********************************************************************
  ! Save 2D array in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_matrix(Mat, fname)
    USE globals
    implicit none

    REAL, DIMENSION(:, :), INTENT(IN) :: Mat
    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_array2D_saving(file_id, Mat, size(Mat, 1), size(Mat, 2), 'mat')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_matrix

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_vector(Vec, fname)
    USE globals
    implicit none

    REAL, DIMENSION(:), INTENT(IN) :: Vec
    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_array1D_saving(file_id, Vec, size(Vec, 1), 'vec')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_vector

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  subroutine HDF5_save_vector_int(Vec, fname)
    USE globals
    implicit none

    INTEGER, DIMENSION(:), INTENT(IN) :: Vec
    character(LEN=*) :: fname
    character(70)  :: npr, nid
    integer :: ierr
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
      write (nid, *) MPIvar%glob_id + 1
      write (npr, *) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
    END IF
    call HDF5_create(fname_complete, file_id, ierr)
    call HDF5_array1D_saving_int(file_id, Vec, size(Vec, 1), 'vec')
    call HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  end subroutine HDF5_save_vector_int

END MODULE in_out
