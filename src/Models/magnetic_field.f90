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
  USE in_out
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
    ALLOCATE (phys%magnetic_psi(nnodes))
    IF (switch%ohmicsrc) THEN
      ALLOCATE (phys%Jtor(nnodes))
    END IF
    IF ((switch%RMP).or.(switch%Ripple)) THEN
      ALLOCATE (phys%Bperturb(nnodes, 3))
    END IF
  END SUBROUTINE initialize_magnetic_field

  !**********************************************
  ! Load magnetic field in phys
  !**********************************************
  SUBROUTINE load_magnetic_field

    phys%B = 0.
    phys%magnetic_flux = 0.
    phys%magnetic_psi = 0.

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

    CASE (80:89)
      ! Magnetic field loaded from file in a cartesian grid
      ! Interpolation is needed
      CALL load_magnetic_field_grid


    END SELECT
#ifdef TOR3D
    ! Magnetic perturbations
    ! RMP part testcase 60-69.
    ! We need to be here to have refElTor%Nodes1d, refElTor%coord1d and numer%ntor (for tdiv)
    IF ((switch%RMP).or.(switch%Ripple)) THEN
      phys%Bperturb = 0.
      CALL addMagneticPerturbation()
    END IF
#endif
    ! Adimensionalization of the magnetic field
    !phys%B = phys%B/phys%B0
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

        CASE (5)
          IF (switch%axisym) THEN
            WRITE (6, *) "This is NOT an axisymmetric test case!"
            stop
          END IF
          !
          phys%B(ind, 1) = 0.
          phys%B(ind, 2) = 0.
          phys%B(ind, 3) = 1.
        CASE (6:7)
          IF (.not.switch%axisym) THEN
            WRITE (6, *) "This is an axisymmetric test case!"
            stop
          END IF
          ! Axysimmetric case
          phys%B(ind, 1) = 0.
          phys%B(ind, 2) = 0.
          phys%B(ind, 3) = 1.+xx
        CASE (50:59)
          write (6, *) "Error in defineMagneticField: you should not be here!"
          STOP
        CASE (60:68)

          ! Circular case with limiter
          R0 = geom%R0
          q = geom%q
          B0 = 2!*0.1522
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

    integer                           :: i, j, ierr, ip, jp, ind, k
    integer(HID_T)                    :: file_id
    real*8, pointer, dimension(:, :)  :: r2D, z2D, flux2D, Br2D, Bz2D, Bphi2D
    real*8, allocatable, dimension(:) :: xvec, yvec
    real*8                            :: x, y, t
    real*8                            :: Br, Bz, Bt, flux, psiSep

    character(LEN=1000) :: fname
    character(50)  :: npr,nid,nit

    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE (6, *) "******* Loading magnetic field *******"
    ENDIF
    ! Read file
    if (switch%testcase>=50 .and. switch%testcase<60) then
       ! WEST case
				   ! Dimensions of the file storing the magnetic field for West
				   ip = 457
       jp = 457
       !ip = 541
       !jp = 391 
       if(switch%testcase .ne. 59) then !if not a moving equilibrium simulation
         fname = 'WEST_54487_0201.h5'
       else
        fname = 'B_field_exp/WEST_54487'
        write(nit, "(i10)") int(time%it + 1)
        nit = trim(adjustl(nit))
        k = INDEX(nit, " ") -1
          ! I don't think this is actually needed for parallel read of the magnetic files
          !IF (MPIvar%glob_size.GT.1) THEN
          !   write(nid,*) MPIvar%glob_id+1
          !   write(npr,*) MPIvar%glob_size
          !   fname = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
          !ELSE
        fname = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
          !END IF
        endif
        IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname))
        ENDIF
        CALL HDF5_open(fname, file_id, IERR)
    elseif (switch%testcase>=80 .and. switch%testcase<90) then
        ! ITER case
    		  if (switch%ME .eq. .FALSE.) then !if not a moving equilibrium simulation
			        !ip = 513
           !jp = 257 
    		     !fname = 'ITER_2008_MagField.h5'
    		     fname = 'ITER_135011_0000.h5'
			     else
			       fname = 'B_field_exp/ITER_135011'
          write(nit, "(i10)") int(time%it + 1)
          nit = trim(adjustl(nit))
          k = INDEX(nit, " ") -1
          fname = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
        endif
        IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname))
        ENDIF
        CALL HDF5_open(fname, file_id, IERR)
        CALL HDF5_integer_reading(file_id, ip, 'ip')
        CALL HDF5_integer_reading(file_id, jp, 'jp') 
    endif
    
    ALLOCATE (r2D(ip, jp))
    ALLOCATE (z2D(ip, jp))
    ALLOCATE (flux2D(ip, jp))
    ALLOCATE (Br2D(ip, jp))
    ALLOCATE (Bz2D(ip, jp))
    ALLOCATE (Bphi2D(ip, jp))
    CALL HDF5_array2D_reading(file_id, r2D, 'r2D')
    CALL HDF5_array2D_reading(file_id, z2D, 'z2D')
    CALL HDF5_real_reading(file_id, psiSep, 'psiSep')
    CALL HDF5_array2D_reading(file_id, flux2D, 'flux2D')
    CALL HDF5_array2D_reading(file_id, Br2D, 'Br2D')
    CALL HDF5_array2D_reading(file_id, Bz2D, 'Bz2D')
    CALL HDF5_array2D_reading(file_id, Bphi2D, 'Bphi2D')
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Min and Max flux for inizialization
    !phys%Flux2Dmin = minval(flux2D)
    !phys%Flux2Dmax = maxval(flux2D)

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
     
    ! Min and Max flux for inizialization
    phys%Flux2Dmin = minval(phys%magnetic_flux)      
    phys%Flux2Dmax = maxval(phys%magnetic_flux) 

#IFDEF PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
#ENDIF 
  
    ! Magnetic flux normalized to separatrix: PSI
    phys%magnetic_psi = (phys%magnetic_flux - phys%Flux2Dmin)/(psiSep - phys%Flux2Dmin) 

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
    character(len=500) :: fname_complete
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
    IF (MPIvar%glob_id .eq. 0) THEN
      write (6, *) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))
    ENDIF
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

#ifdef TOR3D
  !********************************
  ! A subroutine adding RMP or Ripple to an equilibrium not normalized B field
  ! Input:
  !        x, y, t: coordinates x, y and toroidal one
  !        b: (R,Z,Phi) in: computation of B at equilibrium (NOT normalized). out: normalized b with perturbation
  !********************************
  SUBROUTINE addMagneticPerturbation()
    real*8                                   :: xc,yc
    real*8                                   :: xmax,xmin,ymax,ymin,xm,ym,divbp
    real*8                                   :: xx,yy,tt,Bp,Bt,Br,Bz,BB
    integer*4                                :: i,j,ind,N2D,N1D
    real*8,  dimension(Mesh%Nnodes)          :: x, y
    real*8,  dimension(Mesh%Nnodes_toroidal) :: t
    real*8, allocatable                      :: brmp(:,:), bripple(:,:)
    real*8, dimension(2,2)                   :: coilCoord
    integer                                  :: elDiscr, rowNb


    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    N2d = size(x,1)
    N1d = size(t,1)
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    !xmax = 2314.57127827459!Mesh%xmax
    !xmin = 1262.49342451341!Mesh%xmin
    !ymax = 526.038926880589!Mesh%ymax
    !ymin = -526.038926880589!Mesh%ymin
    xm = 0.5*(xmax+xmin)
    ym = 0.5*(ymax+ymin)

    xc = -0.5
    yc = -0.5

    ! RMP part testcase 61-69
    allocate(brmp(N2d*N1d,3))
    allocate(bripple(N2d*N1d,3))
    brmp = 0.
    bripple = 0.
    if (switch%RMP) then ! RMP
      elDiscr = 16
      allocate(magn%coils_rmp(magn%nbCoils_rmp*4*elDiscr,6,magn%nbRow)) !4 because square, 6: 2 positions (start/stop) and 3 coordinates
      if (magn%nbRow.eq.2) then
        ! Upper row
        rowNb = 1
        coilCoord(1,1) = 0.95*xmax ! R coordinate top, upper row
        coilCoord(2,1) = ym + 3./4.*(ymax - ym) ! Z coordinate top, upper row
        coilCoord(1,2) = 1.05*xmax ! R coordinate bottom, upper row
        coilCoord(2,2) = ym + 1./4.*(ymax - ym) ! Z coordinate bottom, upper row
        call calcRMPField(brmp, coilCoord, rowNb, elDiscr)
        ! Lower row
        rowNb = 2
        coilCoord(1,1) = 1.05*xmax ! R coordinate top, lower row
        coilCoord(2,1) = ym - 1./4.*(ymax - ym) ! Z coordinate top, lower row
        coilCoord(1,2) = 0.95*xmax ! R coordinate bottom, lower row
        coilCoord(2,2) = ym - 3./4.*(ymax - ym) ! Z coordinate bottom, lower row
        call calcRMPField(brmp, coilCoord, rowNb, elDiscr)
      else
        WRITE(6, *) 'TODO: not implemented yet'
        stop
      endif
    endif

    if (switch%Ripple) then ! Ripple
      elDiscr = 32
      allocate(magn%coils_ripple(magn%nbCoils_ripple*elDiscr,6))
      call calcRippleField(bripple, elDiscr)
    endif

    DO i=1,N2d
      DO j=1,N1d
        xx = x(i)
        yy = y(i)
        tt = t(j)
        ind = (j-1)*N2d+i

        phys%B(ind,1) = phys%B(ind,1) + brmp(ind,1) + bripple(ind,1)
        phys%B(ind,2) = phys%B(ind,2) + brmp(ind,2) + bripple(ind,2)
        phys%B(ind,3) = phys%B(ind,3) + brmp(ind,3) + bripple(ind,3)

        BB = sqrt(phys%B(ind,1)**2 + phys%B(ind,2)**2 + phys%B(ind,3)**2)
        !phys%B(ind,1) = phys%B(ind,1)/BB
        !phys%B(ind,2) = phys%B(ind,2)/BB
        !phys%B(ind,3) = phys%B(ind,3)/BB
        phys%bperturb(ind,1) = (brmp(ind,1) + bripple(ind,1))!/BB
        phys%bperturb(ind,2) = (brmp(ind,2) + bripple(ind,2))!/BB
        phys%bperturb(ind,3) = (brmp(ind,3) + bripple(ind,3))!/BB
      END DO
    END DO
    deallocate(brmp,bripple)
  END SUBROUTINE addMagneticPerturbation


  !********************************
  ! A subroutine drawing the (rectangular) coils for RMP and computing B RMP field
  ! Input:
  !        brmp: (R,Z,Phi) computation of B (not normalized) generated by coils.
  !        coilCoord: (2,2) with in each colum the (R,Z) coordinates of top (bottom) of the coils in a row
  !           ex: (1,1): R top, (2,1): Z top, (1,2): R bottom, (2,2): Z bottom
  !        rowNb: current number of row
  !********************************
  SUBROUTINE calcRMPField(brmp, coilCoord, rowNb, elDiscr)
    real*8,  dimension(:,:), intent(inout)       :: brmp
    real*8,  dimension(2,2), intent(in)          :: coilCoord
    integer, intent(in)                          :: rowNb,elDiscr

    real*8,  dimension(Mesh%Nnodes)              :: x, y
    real*8,  dimension(Mesh%Nnodes_toroidal)     :: t
    real*8                                       :: tmax, aLcoil, phiCoil
    real*8                                       :: spaceBetwCoil, dlxx, dlyy, dlzz
    real*8                                       :: xxc0, xxc1, yyc0, yyc1, zzc0, zzc1
    real*8                                       :: x0, y0, z0, xxc, yyc, zzc
    real*8                                       :: rrx, rry, rrz, rr2, rr3
    real*8                                       :: Bx, By, Bz, dBx, dBy, dBz
    integer                                      :: i, j, k, l, N2d, N1d, i2d, i1d, ind, par
    real*8, dimension(magn%nbCoils_RMP, 4)       :: xx, yy, zz
    character(len=18)                            :: filename

    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    ! Size of the domain
    tmax = 2*pi !numer%tmax
    N2d = size(x,1)
    N1d = size(t,1)

    aLcoil = (coilCoord(1,1) + coilCoord(1,2))/2.0*magn%torElongCoils_rmp
    spaceBetwCoil = (tmax - magn%nbCoils_rmp*magn%torElongCoils_rmp)/magn%nbCoils_rmp
    if (spaceBetwCoil.le.-0.01) then
      write(6,*) "Error in RMP coils: widths of coils in a row to much for the chosen toroidal expansion"
      write(6,*) "Negative space between coils: ", spaceBetwCoil
      STOP
    endif
    ! Be careful to the direction of the coil (4 corners from bottom right then trigo)
    ! Creation of the n_coils_row by rotation (in row/ toroidal direction)
    do i=1,magn%nbCoils_RMP
      phiCoil = tmax*(i-1)/magn%nbCoils_rmp + spaceBetwCoil

      xx(i,1) = coilCoord(1,2)*cos(phiCoil) - (-aLcoil/2.0)*sin(phiCoil)
      yy(i,1) = coilCoord(1,2)*sin(phiCoil) + (-alcoil/2.0)*cos(phiCoil)
      zz(i,1) = coilCoord(2,2)

      xx(i,2) = coilCoord(1,2)*cos(phiCoil) - (aLcoil/2.0)*sin(phiCoil)
      yy(i,2) = coilCoord(1,2)*sin(phiCoil) + (alcoil/2.0)*cos(phiCoil)
      zz(i,2) = coilCoord(2,2)

      xx(i,3) = coilCoord(1,1)*cos(phiCoil) - (aLcoil/2.0)*sin(phiCoil)
      yy(i,3) = coilCoord(1,1)*sin(phiCoil) + (alcoil/2.0)*cos(phiCoil)
      zz(i,3) = coilCoord(2,1)

      xx(i,4) = coilCoord(1,1)*cos(phiCoil) - (-aLcoil/2.0)*sin(phiCoil)
      yy(i,4) = coilCoord(1,1)*sin(phiCoil) + (-alcoil/2.0)*cos(phiCoil)
      zz(i,4) = coilCoord(2,1)
    enddo

    !Saving the coils coordinates for drawing
    !Loop on coils
    ind = 0
    do i=1,magn%nbCoils_rmp
      !Loop on the 4 parts of a coil: 1->2, 2->3, 3->4, 4->1
      do j=1,4
        if (j.le.3) then
          k = j + 1
        else
          k = 1
        endif
        dlxx = (xx(i,k)-xx(i,j))/elDiscr
        dlyy = (yy(i,k)-yy(i,j))/elDiscr
        dlzz = (zz(i,k)-zz(i,j))/elDiscr
        ! Loop on elements of coils for writing coils coordinates
        do l=1,elDiscr
          ind = ind + 1
          magn%coils_rmp(ind,1,rowNb) = (l-1)*dlxx + xx(i,j)
          magn%coils_rmp(ind,3,rowNb) = (l-1)*dlyy + yy(i,j)
          magn%coils_rmp(ind,5,rowNb) = (l-1)*dlzz + zz(i,j)

          magn%coils_rmp(ind,2,rowNb) = l*dlxx + xx(i,j)
          magn%coils_rmp(ind,4,rowNb) = l*dlyy + yy(i,j)
          magn%coils_rmp(ind,6,rowNb) = l*dlzz + zz(i,j)
        enddo
      enddo
    enddo

    ! RMP hard coded with Biot and Savard law (see phd E. Nardon and ERGOS)
    do i2d=1,N2d
      do i1d=1,N1d
        ! In cartesian coordinates (for Biot and Savard)
        Bx = 0.0
        By = 0.0
        Bz = 0.0

        x0 = sqrt(x(i2d)**2 + y(i2d)**2)*cos(t(i1d))
        y0 = sqrt(x(i2d)**2 + y(i2d)**2)*sin(t(i1d))
        z0 = y(i2d)
        ind = (i1d-1)*N2d+i2d

        !Loop on coils
        do i=1,magn%nbCoils_rmp
          if ((magn%parite.eq.1).or.(magn%parite.eq.-1)) then
            par = magn%parite*(-1)**i
          else
            par = (-1)**i
          endif
          !Loop on the 4 parts of a coil: 1->2, 2->3, 3->4, 4->1
          do j=1,4
            if (j.le.3) then
              k = j + 1
            else
              k = 1
            endif

            dlxx = (xx(i,k)-xx(i,j))/elDiscr
            dlyy = (yy(i,k)-yy(i,j))/elDiscr
            dlzz = (zz(i,k)-zz(i,j))/elDiscr
            ! Loop on elements of coils for writing coils coordinates
            ! rr is the r vector in Biot and Savard
            do l=1,elDiscr
              xxc = (l-1)*dlxx + xx(i,j)
              yyc = (l-1)*dlyy + yy(i,j)
              zzc = (l-1)*dlzz + zz(i,j)

              rrx = x0 - xxc
              rry = y0 - yyc
              rrz = z0 - zzc

              rr2 = rrx**2 + rry**2 + rrz**2
              rr3 = rr2**(3./2)

              dBx = par*(dlyy*rrz-dlzz*rry)/rr3/phys%lscale
              dBy = par*(dlzz*rrx-dlxx*rrz)/rr3/phys%lscale
              dBz = par*(dlxx*rry-dlyy*rrx)/rr3/phys%lscale

              Bx = Bx + dBx
              By = By + dBy
              Bz = Bz + dBz
            enddo ! end elements of one coil
          enddo ! end 4 parts of one coil
        enddo ! end loop on coils
        ! Br, Bz, Bt
        brmp(ind,1) = brmp(ind,1) + magn%amp_rmp*(-Bx*sin(t(i1d)) + By*cos(t(i1d)))
        brmp(ind,2) = brmp(ind,2) + magn%amp_rmp*Bz
        brmp(ind,3) = brmp(ind,3) + magn%amp_rmp*(Bx*cos(t(i1d)) + By*sin(t(i1d)))
      enddo
    enddo
  END SUBROUTINE calcRMPField

  !********************************
  ! A subroutine drawing the toroidal coils for ripple and computing B ripple field
  ! Input:
  !        bripple: (R,Z,Phi) computation of B (not normalized) generated by coils minus averaged value of B.
  !        amp: amplitude of the magnetic field defined as the ratio of the coil current and the plasma current (in TK3X)
  !        triang: triangularity
  !        ellip: ellipticity
  !********************************
  SUBROUTINE calcRippleField(bripple, elDiscr)
    real*8,  dimension(:,:), intent(inout)       :: bripple
    integer, intent(in)                          :: elDiscr

    real*8,  dimension(Mesh%Nnodes)              :: x, y
    real*8,  dimension(Mesh%Nnodes_toroidal)     :: t
    real*8,  dimension(Mesh%Nnodes,3)            :: bripple_av
    real*8                                       :: minRadius, tmax, thetaCoil, phiCoil, csteR
    real*8                                       :: dlxx, dlyy, dlzz
    real*8                                       :: x0, y0, z0, xxc, yyc, zzc
    real*8                                       :: rrx, rry, rrz, rr2, rr3
    real*8                                       :: Bx, By, Bz, dBx, dBy, dBz
    integer                                      :: i, j, k, N2d, N1d, i2d, i1d, ind, loc, par
    real*8, dimension(:,:), allocatable          :: xx, yy, zz
    character(len=15)                            :: filename

    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    ! Size of the domain
    tmax = 2*pi !numer%tmax
    N2d = size(x,1)
    N1d = size(t,1)
    csteR = 1.4
    minRadius = (Mesh%xmax - Mesh%xmin)/2.0

    allocate(xx(magn%nbCoils_ripple, elDiscr))
    allocate(yy(magn%nbCoils_ripple, elDiscr))
    allocate(zz(magn%nbCoils_ripple, elDiscr))

    ! Creation of the N_coils by rotation (in row/ toroidal direction)
    ! Need to create coils a0round full torus to avoid B-field inconsistency
    do i=1,magn%nbCoils_ripple
      ! Shift from 0 for first phi to avoid non-axisymmetry
      phiCoil = tmax*(i-1)/magn%nbCoils_ripple + tmax/(2*magn%nbCoils_ripple)
      ! Theta discretization on Ndiscr points and minor radius of 2*a for toroidal coils
      do j = 1,elDiscr
        thetaCoil = 2*PI*j/elDiscr
        xx(i, j) = (geom%R0/phys%lscale + csteR*minRadius*cos(thetaCoil + magn%triang*sin(thetaCoil)))*cos(phiCoil)
        yy(i, j) = (geom%R0/phys%lscale + csteR*minRadius*cos(thetaCoil + magn%triang*sin(thetaCoil)))*sin(phiCoil)
        zz(i, j) = magn%ellip*csteR*minRadius*sin(thetaCoil)
      enddo
    enddo

    !Saving the coils coordinates for drawing
    !Loop on coils
    ind = 0
    do i=1,magn%nbCoils_ripple
      do j=1,elDiscr
        if (j.lt.elDiscr) then
          k = j + 1
        else
          k = 1
        endif
        ind = ind + 1
        magn%coils_ripple(ind,1) = xx(i,j)
        magn%coils_ripple(ind,3) = yy(i,j)
        magn%coils_ripple(ind,5) = zz(i,j)

        magn%coils_ripple(ind,2) = xx(i,k)
        magn%coils_ripple(ind,4) = yy(i,k)
        magn%coils_ripple(ind,6) = zz(i,k)
      enddo
    enddo

    ! Ripple hard coded with Biot and Savard law (see phd E. Nardon and ERGOS)
    do i2d=1,N2d
      do i1d=1,N1d
        ! In cartesian coordinates (for Biot and Savard)
        Bx = 0.0
        By = 0.0
        Bz = 0.0

        x0 = sqrt(x(i2d)**2 + y(i2d)**2)*cos(t(i1d))
        y0 = sqrt(x(i2d)**2 + y(i2d)**2)*sin(t(i1d))
        z0 = y(i2d)
        ind = (i1d-1)*N2d+i2d

        !Loop on coils
        do i=1,magn%nbCoils_ripple
          par = 1 !(-1)**i
          ! Loop on elements of coils for writing coils coordinates
          ! rr is the r vector in Biot and Savard
          do j=1,elDiscr
            if (j.lt.elDiscr) then
              k = j + 1
            else
              k = 1
            endif

            rrx = x0 - xx(i,j)
            rry = y0 - yy(i,j)
            rrz = z0 - zz(i,j)

            dlxx = xx(i,k)-xx(i,j)
            dlyy = yy(i,k)-yy(i,j)
            dlzz = zz(i,k)-zz(i,j)

            rr2 = rrx**2 + rry**2 + rrz**2
            rr3 = rr2**(3./2)

            dBx = par*(dlyy*rrz-dlzz*rry)/rr3/phys%lscale
            dBy = par*(dlzz*rrx-dlxx*rrz)/rr3/phys%lscale
            dBz = par*(dlxx*rry-dlyy*rrx)/rr3/phys%lscale

            Bx = Bx + dBx
            By = By + dBy
            Bz = Bz + dBz
          enddo ! end elements of one coil
        enddo ! end loop on coils
        ! Br, Bz, Bt
        bripple(ind,1) = bripple(ind,1) + magn%amp_ripple*(-Bx*sin(t(i1d)) + By*cos(t(i1d)))
        bripple(ind,2) = bripple(ind,2) + magn%amp_ripple*Bz
        bripple(ind,3) = bripple(ind,3) + magn%amp_ripple*(Bx*cos(t(i1d)) + By*sin(t(i1d)))
      enddo
    enddo
    ! Compute the ripple average
    bripple_av = 0.
    do i2d=1,N2d
      do i1d=1,N1d
        ind = (i1d-1)*N2d+i2d
        bripple_av(i2d,1) = bripple_av(i2d,1) + bripple(ind,1)
        bripple_av(i2d,2) = bripple_av(i2d,2) + bripple(ind,2)
        bripple_av(i2d,3) = bripple_av(i2d,3) + bripple(ind,3)
      enddo
    enddo
    bripple_av = bripple_av/size(t,1)
    ! Substract the mean along phi to the magnetic field to only keep the ripple
    ! minloc can be necessary to find the correct index in the global matrix bripple_av to substract to the local matrix bripple
    ! When bripple and phys%bripple are the same size, loc must be equal to ind (which unfortunately is NOT the case for a few
    ! index).
    do i2d=1,N2d
      !loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1)
      ! Fortran 2008 only !!!
      !if (any(mesh%T(1:size(T,1)/2,:).eq.loc)) then
      !   loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1, BACK=.FALSE.)
      !else
      !   loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1, BACK=.TRUE.)
      !endif
      do i1d=1,N1d
        ind = (i1d-1)*N2d+i2d
        !bripple(ind,1) = bripple(ind,1) - phys%bripple(loc, 1)
        !bripple(ind,2) = bripple(ind,2) - phys%bripple(loc, 2)
        !bripple(ind,3) = bripple(ind,3) - phys%bripple(loc, 3)
        bripple(ind,1) = bripple(ind,1) - bripple_av(i2d, 1)
        bripple(ind,2) = bripple(ind,2) - bripple_av(i2d, 2)
        bripple(ind,3) = bripple(ind,3) - bripple_av(i2d, 3)
      enddo
    enddo
    deallocate(xx,yy,zz)
  END SUBROUTINE calcRippleField
#endif

  ! Below are routines from Manuel MHDG v2.1. Copy as it without any check: TODO adapt it to global magnetic field

  SUBROUTINE loadJtorMap()
    USE interpolation
    USE HDF5
    USE HDF5_io_module
    !USE MPI_OMP
    integer        :: i,ierr,ip,jp,ind,j,k
    integer(HID_T) :: file_id

    character(LEN=1000)    :: fname
    character(70)        :: npr,nid,nit

    real*8,pointer,dimension(:,:) :: r2D,z2D,Jtor
    real*8,allocatable,dimension(:)   :: xvec,yvec
    real*8                            :: x,y

    IF (MPIvar%glob_id .eq. 0) THEN
      WRITE(6,*) "******* Loading Toroidal Current *******"
    ENDIF
    ! Allocate storing space in phys
    ALLOCATE(phys%Jtor(Mesh%Nnodes))
    phys%Jtor = 0.

    ! Read file
    if (switch%testcase>=50 .and. switch%testcase<60) then
       ! WEST case
			 ! Dimensions of the file storing the magnetic field for West
			 ip = 457
       jp = 457
       !ip = 541
       !jp = 391 
       if (switch%ME .eq. .FALSE.) then !if not a moving equilibrium simulation
          fname = 'WEST_54487_Jtor_0201.h5'
       else
          fname = 'B_field_exp/WEST_54487_Jtor'
          write(nit, "(i10)") int(time%it + 1)
          nit = trim(adjustl(nit))
          k = INDEX(nit, " ") -1
          fname = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
        endif
        IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'Toroidal current loaded from file: ', trim(adjustl(fname))
        ENDIF
        CALL HDF5_open(fname, file_id, IERR)
    elseif (switch%testcase>=80 .and. switch%testcase<90) then
        ! ITER case
    		if(switch%ME .eq. .FALSE.) then !if not a moving equilibrium simulation
    	    !fname = 'ITER_2008_MagField.h5'
    	    fname = 'ITER_135011_Jtor_0000.h5'
			  else
			    fname = 'B_field_exp/ITER_135011_Jtor'
			    write(nit, "(i10)") int(time%it + 1)
          nit = trim(adjustl(nit))
          k = INDEX(nit, " ") -1
          fname = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
        endif
        IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'Toroidal current loaded from file: ', trim(adjustl(fname))
        ENDIF
        CALL HDF5_open(fname, file_id, IERR)
        CALL HDF5_integer_reading(file_id, ip, 'ip')
        CALL HDF5_integer_reading(file_id, jp, 'jp')
    endif
    
    ALLOCATE(r2D(ip,jp))
    ALLOCATE(z2D(ip,jp))
    ALLOCATE(Jtor(ip,jp))

    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
    CALL HDF5_array2D_reading(file_id,Jtor,'Jtor')
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Interpolate
    ALLOCATE(xvec(jp))
    ALLOCATE(yvec(ip))
    xvec = r2D(1,:)
    yvec = z2D(:,1)
    DO i = 1,Mesh%Nnodes
      x = Mesh%X(i,1)
      y = Mesh%X(i,2)
      ind = i
#ifdef TOR3D
      DO j = 1, Mesh%Nnodes_toroidal
        ind = (j - 1)*Mesh%Nnodes + i
#endif
        phys%Jtor(ind) = interpolate(ip, yvec,jp, xvec,Jtor, y,x, 1e-12)
#ifdef TOR3D
      END DO
#endif
    END DO

    !Compute Ip
    CALL computeIplasma()
    
    IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'I_p =  ', phys%I_p, '[MA]'
    ENDIF    

    ! Free memory
    DEALLOCATE(r2D,z2D,Jtor,xvec,yvec)

  END SUBROUTINE loadJtorMap

  SUBROUTINE loadMagneticFieldFromExperimentalData()
    USE interpolation
    USE HDF5
    USE HDF5_io_module
    USE MPI_OMP
    integer        :: ierr,k,ip,jp,i,j,ind
    real*8,pointer,dimension(:,:) :: r2D,z2D,flux2D,Br2D,Bz2D,Bphi2D
    real*8,allocatable,dimension(:)   :: xvec,yvec,Bmod
    real*8                            :: x,y, Br, Bz, Bt, flux
    character(LEN=25) :: fname = 'WEST_54487'
    character(10)  :: npr,nid,nit
    character(len=1000) :: fname_complete
    integer(HID_T) :: file_id

    WRITE(6,*) "******* Loading magnetic field *******"

    ! Dimensions of the file storing the magnetic field for West
    ip = 457
    jp = 457
    ALLOCATE(r2D(ip,jp))
    ALLOCATE(z2D(ip,jp))
    ALLOCATE(flux2D(ip,jp))
    ALLOCATE(Br2D(ip,jp))
    ALLOCATE(Bz2D(ip,jp))
    ALLOCATE(Bphi2D(ip,jp))
    ALLOCATE(Bmod(Mesh%Nnodes))

    ! File name
    write(nit, "(i10)") time%it+1
    nit = trim(adjustl(nit))
    k = INDEX(nit, " ") -1

    IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
    ELSE
      fname_complete = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
    END IF

    write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))

    ! Read file
    CALL HDF5_open(fname_complete,file_id,IERR)
    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
    CALL HDF5_array2D_reading(file_id,flux2D,'flux2D')
    CALL HDF5_array2D_reading(file_id,Br2D,'Br2D')
    CALL HDF5_array2D_reading(file_id,Bz2D,'Bz2D')
    CALL HDF5_array2D_reading(file_id,Bphi2D,'Bphi2D')
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
      x = Mesh%X(i,1)
      y = Mesh%X(i,2)
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
    DEALLOCATE(Br2D,Bz2D,Bphi2D,xvec,yvec)
    DEALLOCATE(r2D,z2D,flux2D)

  END SUBROUTINE loadMagneticFieldFromExperimentalData

  SUBROUTINE computeIplasma()
    real*8  	:: Xel(refElPol%Nnodes2D,2),xyg(refElPol%NGauss2D,2),xyg_d(refElPol%NGauss2D,2),dvolu
    real*8		:: Jtorel(refElPol%Nnodes2D),Jtorg(refElPol%NGauss2D)
    real*8		:: J11(refElPol%NGauss2D),J12(refElPol%NGauss2D)
    real*8		:: J21(refElPol%NGauss2D),J22(refElPol%NGauss2D)
    real*8		:: detJ(refElPol%NGauss2D)
    integer	:: iel,g
    real*8  	:: xyDerNorm_g
    
    phys%I_p = 0.
    
    DO iel = 1, Mesh%Nelems
    	! Coordinates of the nodes of the element
    	Xel = Mesh%X(Mesh%T(iel,:),:)
        
   	! Toroidal current of the nodes of the element
        Jtorel = phys%Jtor(Mesh%T(iel,:))
   
       ! Gauss points position
       xyg = matmul(refElPol%N2D,Xel)
  
       ! Toroidal current at Gauss points
      Jtorg = matmul(refElPol%N2D,Jtorel)
       
       ! Jacobian
       J11 = matmul(refElPol%Nxi2D,Xel(:,1))                             ! ng x 1
       J12 = matmul(refElPol%Nxi2D,Xel(:,2))                             ! ng x 1
       J21 = matmul(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
       J22 = matmul(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                 			           ! determinant of the Jacobian

#ifdef PARALL
       IF (Mesh%ghostElems(iel) .eq. 0) THEN
#endif
       
          ! Loop in 2D Gauss points
          DO g = 1, refElPol%NGauss2D
             ! Integration weight
             dvolu = refElPol%gauss_weights2D(g)*detJ(g)

             ! Compute I plasma
             phys%I_p = phys%I_p + Jtorg(g)*dvolu*phys%lscale**2 	 
       
          END DO
      		
#ifdef PARALL
       ENDIF
#endif
     
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%I_p, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    
    ! Toroidal current in [MA]
    phys%I_p = phys%I_p/1.e6                                                   
    
  END SUBROUTINE computeIplasma


  SUBROUTINE SetPuff()
    USE HDF5
    USE HDF5_io_module
    USE interpolation
    integer           :: ierr
    character(LEN=25) :: fname = 'Puff_54487.h5'
    integer(HID_T)    :: file_id
    integer           :: qp, Nn2D
    integer           :: T(Mesh%Nelems,refElPol%Nnodes2D)  
    real*8            :: lower, upper, nli, n_Gw, n_la, pi = 3.1416, a = 2.
    real*8            :: X(Mesh%Nnodes,2), u(Mesh%Nelems*refElPol%Nnodes2D,phys%Neq)
    real*8            :: linex(1000), liney(1000), n_i(Mesh%Nelems*refElPol%Nnodes2D)

    ! Allocate storing space in phys (puff for WEST, 403 entries)
    IF (switch%testcase .ge. 50 .and. switch%testcase .le. 59) THEN
       !ALLOCATE(phys%puff_exp(403))
       IF (MPIvar%glob_id .eq. 0) THEN
          WRITE (6, *) "******* Loading puff *******"
       ENDIF

       ! Read file
       CALL HDF5_open(fname,file_id,IERR)
       CALL HDF5_array1D_reading(file_id,phys%puff_exp,'puff')
       IF (MPIvar%glob_id .eq. 0) THEN
          write(6,*) 'Puff loaded from file: ', trim(adjustl(fname))
       ENDIF
       CALL HDF5_close(file_id)
    END IF
    
    ! ITER puff: linear increase up to nli = 4.00E+19
    IF (switch%testcase .ge. 80 .and. switch%testcase .le. 89) THEN
        ! Puff feedback: check if central line integrated has reached the target value (4e19 for ITER ohmic)
        qp = size(linex)
        Nn2D = refElPol%Nnodes2D
        X = mesh%X
        T = mesh%T
        nli = 0. 
        lower = minval(Mesh%X(:,1))
        upper = maxval(Mesh%X(:,1))
        linex = (/(lower + (upper - lower)/1000.*(i-1), i=1, 1000)/)
        liney = 0.5/phys%lscale
        u = transpose(reshape(sol%u,[phys%Neq,size(sol%u)/phys%Neq]))
        n_i = u(:,1)
#ifdef PARALL
        ! Check in the case of horizontal partition to not waste time
        IF (maxval(Mesh%X(:, 2)) .gt. liney(1) .and. minval(Mesh%X(:, 2)) .lt. liney(1) ) THEN
           CALL lineintegration(qp, linex, liney, n_i, X, T, Nn2D, nli)
        END IF
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, nli, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        n_la = nli/(8.3659 - 4.04)*simpar%refval_density
        IF (MPIvar%glob_id .eq. 0) THEN
           WRITE(6,*) 'nli = ', nli, 'E+19 [m^-2]'
           WRITE(6,*) 'n_la = ', n_la, '[m^-3]'
        END IF
#else
        CALL lineintegration(qp, linex, liney, n_i, X, T, Nn2D, nli)
        n_la = nli/(8.3659 - 4.04)*simpar%refval_density
        WRITE(6,*) 'nli = ', nli, 'E+19 [m^-2]'
        WRITE(6,*) n_la, '[m^-3]'
#endif
        ! Upgrade puff
        n_Gw = phys%I_p/(pi*a**2)*10.*simpar%refval_density
        !IF (nli .lt. 4) THEN
        !    IF (time%it .eq. 0) THEN
        !       phys%puff_exp(time%it+1) = phys%puff_slope*time%dt*simpar%refval_time
        !       phys%puff = phys%puff_exp(time%it+1)
        !    ELSE
        !        phys%puff_exp(time%it+1) = phys%puff_exp(time%it) + phys%puff_slope*time%dt*simpar%refval_time
        !       phys%puff = phys%puff_exp(time%it+1)
        !    END IF
        !    IF (MPIvar%glob_id .eq. 0) THEN
        !       WRITE (6, '(" * Puff = ", E10.3, 27X, " *")')  phys%puff
        !    END IF
        !ELSE
        !    phys%puff_exp(time%it+1) = phys%puff_exp(time%it) - exp(nli - 4)*simpar%refval_density*time%dt*simpar%refval_time
        !    phys%puff = max(phys%puff_exp(time%it+1),0.)
        !    IF (MPIvar%glob_id .eq. 0) THEN
        !        WRITE (6, '(" * Puff = ", E10.3, 27X, " *")')  phys%puff
        !    END IF
        !END IF
        IF (time%it .eq. 0) THEN
           phys%puff_exp(time%it+1) = max(10.*(phys%puff_slope*n_Gw - n_la), 0.)
                phys%puff = phys%puff_exp(time%it+1)
        ELSE
           phys%puff_exp(time%it+1) = max(phys%puff_exp(time%it) + 50.*(2 - sign(1.,phys%puff_slope*n_Gw - n_la))*(phys%puff_slope*n_Gw - n_la), 0.)
           !phys%puff_exp(time%it+1) = max(phys%puff_exp(time%it) + 50.*(phys%puff_slope*n_Gw - n_la), 0.)   
           phys%puff = phys%puff_exp(time%it+1)
       END IF
       IF (MPIvar%glob_id .eq. 0) THEN
	       WRITE (6, '(" * Puff = ", E10.3, 27X, " *")')  phys%puff
       END IF
   
    END IF
          
  END SUBROUTINE SetPuff

END MODULE Magnetic_field
