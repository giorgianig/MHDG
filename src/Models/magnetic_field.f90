!*****************************************
! project: MHDG
! file: magnetic_field.f90
! date: 06/04/2020
! Definition of the magnetic field
! The magnetic field is defined at the mesh
! nodes
!*****************************************

MODULE magnetic_field
   USE prec_const
   USE globals
CONTAINS

   !**********************************************
   ! Allocate storing space for the magnetic field
   !**********************************************
   SUBROUTINE initialize_magnetic_field
      INTEGER  :: nnodes
#ifdef TOR3D
      nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
      nnodes = Mesh%Nnodes
#endif
      ! Allocate storing space in phys
      ALLOCATE (phys%B(nnodes, 3))
      ALLOCATE (phys%magnetic_flux(nnodes))
   END SUBROUTINE initialize_magnetic_field

   !**********************************************
   ! Load magnetic field in phys
   !**********************************************
   SUBROUTINE load_magnetic_field

      phys%B = 0.
      phys%magnetic_flux = 0.

      SELECT CASE (switch%testcase)
      CASE (1:49)
         ! Analytic definition of the magnetic field
         CALL load_magnetic_field_analytical

      CASE (50:59)
         ! Magnetic field loaded from file in a cartesian grid
         ! Interpolation is needed
         CALL load_magnetic_field_grid

      CASE (60:69)
         ! Analytic definition of the magnetic field
         CALL load_magnetic_field_analytical

      CASE (70:79)
         ! Magnetic field loaded from file in the mesh nodes
         CALL load_magnetic_field_nodes
      END SELECT

      ! Adimensionalization of the magnetic field
      phys%B = phys%B/phys%B0
   END SUBROUTINE load_magnetic_field

   !***********************************************************************
   ! Magnetic field defined analytically
   !***********************************************************************
   SUBROUTINE load_magnetic_field_analytical
      real*8   :: x(Mesh%Nnodes), y(Mesh%Nnodes)  ! Coordinates in the plane
#ifdef TOR3D
      real*8   :: t(Mesh%Nnodes_toroidal)       ! Toroidal coordinates
#endif
      real*8                  :: xc, yc, R0, q, r
      real*8                  :: xmax, xmin, ymax, ymin, xm, ym, p, divbp
      real*8                  :: xx, yy, tt, Bp, Bt, Br, Bz, BB, dmax, B0, xr, yr
      integer*4               :: i, j, ind, N2D, N1D

      x = Mesh%X(:, 1)
      y = Mesh%X(:, 2)
      N2d = size(X, 1)
#ifdef TOR3D
      t = Mesh%toroidal
      N1d = size(t, 1)
#endif
      xmax = Mesh%xmax
      xmin = Mesh%xmin
      ymax = Mesh%ymax
      ymin = Mesh%ymin
      xm = 0.5*(xmax + xmin)
      ym = 0.5*(ymax + ymin)

      xc = 0.
      yc = 0.
      DO i = 1, N2d
         xx = x(i)
         yy = y(i)
         ind = i
#ifdef TOR3D
         DO j = 1, N1d
            tt = t(j)
            ind = (j - 1)*N2d+i
#endif
            SELECT CASE (switch%testcase)
            CASE (1)
               IF (switch%axisym) THEN
                  WRITE (6, *) "This is NOT an axisymmetric test case!"
                  stop
               END IF
               ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
               phys%B(ind, 1) = (yy - yc)
               phys%B(ind, 2) = (-xx + xc)
               phys%B(ind, 3) = 1.

            CASE (2)
               IF (.not. switch%axisym) THEN
                  WRITE (6, *) "This is an axisymmetric test case!"
                  stop
               END IF
               ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
               phys%B(ind, 1) = (yy - ym)/xx
               phys%B(ind, 2) = (-xx + xm)/xx
               phys%B(ind, 3) = 1.
            CASE (50:59)
               write (6, *) "Error in defineMagneticField: you should not be here!"
               STOP
            CASE (60:69)

               ! Circular case with limiter
               R0 = geom%R0
               q = geom%q
               B0 = 2*0.1522 ! not influential
               xr = xx*phys%lscale
               yr = yy*phys%lscale

               r = sqrt((xr - R0)**2 + yr**2)
               phys%B(ind, 1) = -B0*yr/(xr*q*sqrt(1 - (r/R0)**2))
               phys%B(ind, 2) = B0*(xr - R0)/(xr*q*sqrt(1 - (r/R0)**2))
               phys%B(ind, 3) = B0*R0/xr

            CASE DEFAULT
               WRITE (6, *) "Error! Test case not valid"
               STOP
            END SELECT
#ifdef TOR3D
         END DO
#endif
      END DO

   END SUBROUTINE load_magnetic_field_analytical

   !***********************************************************************
   ! Magnetic field loaded by a hdf5 file in a cartesian grid
   !***********************************************************************
   SUBROUTINE load_magnetic_field_grid()
      USE interpolation
      USE HDF5
      USE HDF5_io_module
      integer        :: i, j, ierr, ip, jp, ind
      integer(HID_T) :: file_id
      real*8, pointer, dimension(:, :) :: r2D, z2D, flux2D, Br2D, Bz2D, Bphi2D
      real*8, allocatable, dimension(:)   :: xvec, yvec
      real*8                            :: x, y, t
      real*8                            :: Br, Bz, Bt, flux

      WRITE (6, *) "******* Loading magnetic field *******"

      ! Dimensions of the file storing the magnetic field for West
      ip = 541
      jp = 391
      ALLOCATE (r2D(ip, jp))
      ALLOCATE (z2D(ip, jp))
      ALLOCATE (flux2D(ip, jp))
      ALLOCATE (Br2D(ip, jp))
      ALLOCATE (Bz2D(ip, jp))
      ALLOCATE (Bphi2D(ip, jp))

      ! Read file
      CALL HDF5_open('WEST_far_465.h5', file_id, IERR)
      CALL HDF5_array2D_reading(file_id, r2D, 'r2D')
      CALL HDF5_array2D_reading(file_id, z2D, 'z2D')
      CALL HDF5_array2D_reading(file_id, flux2D, 'flux2D')
      CALL HDF5_array2D_reading(file_id, Br2D, 'Br2D')
      CALL HDF5_array2D_reading(file_id, Bz2D, 'Bz2D')
      CALL HDF5_array2D_reading(file_id, Bphi2D, 'Bphi2D')
      CALL HDF5_close(file_id)

      ! Apply length scale
      r2D = r2D/phys%lscale
      z2D = z2D/phys%lscale

      ! Interpolate
      ALLOCATE (xvec(jp))
      ALLOCATE (yvec(ip))
      xvec = r2D(1, :)
      yvec = z2D(:, 1)
      DO i = 1, Mesh%Nnodes
         x = Mesh%X(i, 1)
         y = Mesh%X(i, 2)
         Br = interpolate(ip, yvec, jp, xvec, Br2D, y, x, 1e-12)
         Bz = interpolate(ip, yvec, jp, xvec, Bz2D, y, x, 1e-12)
         Bt = interpolate(ip, yvec, jp, xvec, Bphi2D, y, x, 1e-12)
         flux = interpolate(ip, yvec, jp, xvec, flux2D, y, x, 1e-12)
         ind = i
#ifdef TOR3D
         DO j = 1, Mesh%Nnodes_toroidal
            ind = (j - 1)*Mesh%Nnodes + i
#endif
            phys%B(ind, 1) = Br
            phys%B(ind, 2) = Bz
            phys%B(ind, 3) = Bt
            phys%magnetic_flux(ind) = flux
#ifdef TOR3D
         END DO
#endif
      END DO

      ! Free memory
      DEALLOCATE (Br2D, Bz2D, Bphi2D, xvec, yvec)
      DEALLOCATE (r2D, z2D, flux2D)

   END SUBROUTINE load_magnetic_field_grid

   !***********************************************************************
   ! Magnetic field loaded by a hdf5 file in the nodes !TODO modify for 3D
   !***********************************************************************
   SUBROUTINE load_magnetic_field_nodes()
      USE HDF5
      USE HDF5_io_module
      USE MPI_OMP
      integer        :: i, ierr, k
      character(LEN=20) :: fname = 'Evolving_equilibrium'
      character(10)  :: npr, nid, nit
      character(len=1000) :: fname_complete
      integer(HID_T) :: file_id
      real*8, pointer, dimension(:) :: Br, Bz, Bt, flux
      INTEGER  :: nnodes
#ifdef TOR3D
      nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
      nnodes = Mesh%Nnodes
#endif
      WRITE (6, *) "******* Loading magnetic field *******"

      ALLOCATE (flux(nnodes))
      ALLOCATE (Br(nnodes))
      ALLOCATE (Bz(nnodes))
      ALLOCATE (Bt(nnodes))

      ! File name
      write (nit, "(i10)") time%it
      nit = trim(adjustl(nit))
      k = INDEX(nit, " ") - 1

      IF (MPIvar%glob_size .GT. 1) THEN
         write (nid, *) MPIvar%glob_id + 1
         write (npr, *) MPIvar%glob_size
                                           fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
      ELSE
         fname_complete = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
      END IF

      write (6, *) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))

      ! Read file
      CALL HDF5_open(fname_complete, file_id, IERR)
      CALL HDF5_array1D_reading(file_id, Br, 'Br')
      CALL HDF5_array1D_reading(file_id, Bz, 'Bz')
      CALL HDF5_array1D_reading(file_id, Bt, 'Bt')
      CALL HDF5_array1D_reading(file_id, flux, 'flux')
      CALL HDF5_close(file_id)

      phys%B(:, 1) = Br
      phys%B(:, 2) = Bz
      phys%B(:, 3) = Bt
      phys%magnetic_flux = flux
      DEALLOCATE (Br, Bz, Bt, flux)
   END SUBROUTINE load_magnetic_field_nodes

END MODULE Magnetic_field
