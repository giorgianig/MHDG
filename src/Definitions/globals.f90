!*****************************************
! project: MHDG
! file: globals.f90
! date: 06/09/2016
! Definition of the global variable structure
!*****************************************
MODULE globals
  USE prec_const
  USE types
  USE matrices_types
  IMPLICIT NONE

  !*******************************************
  ! Declaretion of the GLOBAL variables
  !*******************************************

  TYPE(Reference_element_type), target :: refElPol
  TYPE(Reference_element_type), target :: refElTor
  TYPE(Mesh_type), target :: mesh
  TYPE(Physics_type), target :: phys
  TYPE(Geometry_type), target :: geom
  TYPE(Magnetic_type), target :: magn
  TYPE(Switches_type), target :: switch
  TYPE(Time_type), target :: time
  TYPE(Numeric_type), target :: numer
  TYPE(Utils_type), target :: utils
  TYPE(Lssolver_type), target :: lssolver
  TYPE(Elmat_type), target :: elmat
  TYPE(Sol_type), target :: sol
  TYPE(MAT_CSR_TYP), target :: MatK
  TYPE(RHS_TYP), target :: rhs
  TYPE(Simulationparams_type), target :: simpar
  TYPE(Timing_type), target :: timing

CONTAINS

  !***********************************************
  ! Routine to deallocate all previously allocated
  ! structures
  !***********************************************
  SUBROUTINE free_all()

    ! Reference element
    IF (ALLOCATED(refElPol%Face_nodes)) THEN
      DEALLOCATE (refElPol%Face_nodes)
    END IF
    IF (ALLOCATED(refElPol%Bord_nodes)) THEN
      DEALLOCATE (refElPol%Bord_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes)) THEN
      DEALLOCATE (refElPol%inner_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes_face)) THEN
      DEALLOCATE (refElPol%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElPol%coord3D)) THEN
      DEALLOCATE (refElPol%coord3D)
    END IF
    IF (ASSOCIATED(refElPol%coord2D)) THEN
      DEALLOCATE (refElPol%coord2D)
    END IF
    IF (ASSOCIATED(refElPol%coord1D)) THEN
      DEALLOCATE (refElPol%coord1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points3D)) THEN
      DEALLOCATE (refElPol%gauss_points3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points2D)) THEN
      DEALLOCATE (refElPol%gauss_points2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points1D)) THEN
      DEALLOCATE (refElPol%gauss_points1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights3D)) THEN
      DEALLOCATE (refElPol%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights2D)) THEN
      DEALLOCATE (refElPol%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights1D)) THEN
      DEALLOCATE (refElPol%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElPol%N3D)) THEN
      DEALLOCATE (refElPol%N3D)
    END IF
    IF (ALLOCATED(refElPol%Nxi3D)) THEN
      DEALLOCATE (refElPol%Nxi3D)
    END IF
    IF (ALLOCATED(refElPol%Neta3D)) THEN
      DEALLOCATE (refElPol%Neta3D)
    END IF
    IF (ALLOCATED(refElPol%Nzeta3D)) THEN
      DEALLOCATE (refElPol%Nzeta3D)
    END IF
    IF (ALLOCATED(refElPol%N2D)) THEN
      DEALLOCATE (refElPol%N2D)
    END IF
    IF (ALLOCATED(refElPol%Nxi2D)) THEN
      DEALLOCATE (refElPol%Nxi2D)
    END IF
    IF (ALLOCATED(refElPol%Neta2D)) THEN
      DEALLOCATE (refElPol%Neta2D)
    END IF
    IF (ALLOCATED(refElPol%N1D)) THEN
      DEALLOCATE (refElPol%N1D)
    END IF
    IF (ALLOCATED(refElPol%Nxi1D)) THEN
      DEALLOCATE (refElPol%Nxi1D)
    END IF

    ! Reference element
    IF (ALLOCATED(refElTor%Face_nodes)) THEN
      DEALLOCATE (refElTor%Face_nodes)
    END IF
    IF (ALLOCATED(refElTor%Bord_nodes)) THEN
      DEALLOCATE (refElTor%Bord_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes)) THEN
      DEALLOCATE (refElTor%inner_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes_face)) THEN
      DEALLOCATE (refElTor%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElTor%coord3D)) THEN
      DEALLOCATE (refElTor%coord3D)
    END IF
    IF (ASSOCIATED(refElTor%coord2D)) THEN
      DEALLOCATE (refElTor%coord2D)
    END IF
    IF (ASSOCIATED(refElTor%coord1D)) THEN
      DEALLOCATE (refElTor%coord1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points3D)) THEN
      DEALLOCATE (refElTor%gauss_points3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points2D)) THEN
      DEALLOCATE (refElTor%gauss_points2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points1D)) THEN
      DEALLOCATE (refElTor%gauss_points1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights3D)) THEN
      DEALLOCATE (refElTor%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights2D)) THEN
      DEALLOCATE (refElTor%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights1D)) THEN
      DEALLOCATE (refElTor%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElTor%N3D)) THEN
      DEALLOCATE (refElTor%N3D)
    END IF
    IF (ALLOCATED(refElTor%Nxi3D)) THEN
      DEALLOCATE (refElTor%Nxi3D)
    END IF
    IF (ALLOCATED(refElTor%Neta3D)) THEN
      DEALLOCATE (refElTor%Neta3D)
    END IF
    IF (ALLOCATED(refElTor%Nzeta3D)) THEN
      DEALLOCATE (refElTor%Nzeta3D)
    END IF
    IF (ALLOCATED(refElTor%N2D)) THEN
      DEALLOCATE (refElTor%N2D)
    END IF
    IF (ALLOCATED(refElTor%Nxi2D)) THEN
      DEALLOCATE (refElTor%Nxi2D)
    END IF
    IF (ALLOCATED(refElTor%Neta2D)) THEN
      DEALLOCATE (refElTor%Neta2D)
    END IF
    IF (ALLOCATED(refElTor%N1D)) THEN
      DEALLOCATE (refElTor%N1D)
    END IF
    IF (ALLOCATED(refElTor%Nxi1D)) THEN
      DEALLOCATE (refElTor%Nxi1D)
    ENDIF

    ! Mesh
    IF (ASSOCIATED(Mesh%T)) THEN
      DEALLOCATE (Mesh%T)
    END IF
    IF (ASSOCIATED(Mesh%Tlin)) THEN
      DEALLOCATE (Mesh%Tlin)
    END IF
    IF (ASSOCIATED(Mesh%Tb)) THEN
      DEALLOCATE (Mesh%Tb)
    END IF
    IF (ASSOCIATED(Mesh%boundaryFlag)) THEN
      DEALLOCATE (Mesh%boundaryFlag)
    END IF
    IF (ALLOCATED(Mesh%F)) THEN
      DEALLOCATE (Mesh%F)
    END IF
    IF (ALLOCATED(Mesh%N)) THEN
      DEALLOCATE (Mesh%N)
    END IF
    IF (ALLOCATED(Mesh%faces)) THEN
      DEALLOCATE (Mesh%faces)
    END IF
    IF (ALLOCATED(Mesh%extfaces)) THEN
      DEALLOCATE (Mesh%extfaces)
    END IF
    IF (ALLOCATED(Mesh%intfaces)) THEN
      DEALLOCATE (Mesh%intfaces)
    END IF
    IF (ALLOCATED(Mesh%flipFace)) THEN
      DEALLOCATE (Mesh%flipface)
    END IF
    IF (ALLOCATED(Mesh%Fdir)) THEN
      DEALLOCATE (Mesh%Fdir)
    END IF
    IF (ALLOCATED(Mesh%Diric)) THEN
      DEALLOCATE (Mesh%Diric)
    END IF
    IF (ALLOCATED(Mesh%numberbcs)) THEN
      DEALLOCATE (Mesh%numberbcs)
    END IF
    IF (ASSOCIATED(Mesh%X)) THEN
      DEALLOCATE (Mesh%X)
    END IF
    IF (ALLOCATED(Mesh%elemSize)) THEN
      DEALLOCATE (Mesh%elemSize)
    END IF
    IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
      DEALLOCATE (Mesh%scdiff_nodes)
    END IF
    IF (ASSOCIATED(Mesh%toroidal)) THEN
      DEALLOCATE (Mesh%toroidal)
    END IF
    IF (ALLOCATED(Mesh%periodic_faces)) THEN
      DEALLOCATE (Mesh%periodic_faces)
    END IF
    ! sol type
    IF (ASSOCIATED(sol%u)) THEN
      DEALLOCATE (sol%u)
    END IF
    IF (ASSOCIATED(sol%u_tilde)) THEN
      DEALLOCATE (sol%u_tilde)
    END IF
    IF (ALLOCATED(sol%u0)) THEN
      DEALLOCATE (sol%u0)
    END IF
    IF (ALLOCATED(sol%tres)) THEN
      DEALLOCATE (sol%tres)
    END IF
    IF (ALLOCATED(sol%time)) THEN
      DEALLOCATE (sol%time)
    END IF

    IF (ALLOCATED(elMat%iAqq)) THEN
      DEALLOCATE (elMat%iAqq)
    END IF
    IF (ALLOCATED(elMat%Aqu)) THEN
      DEALLOCATE (elMat%Aqu)
    END IF
    IF (ALLOCATED(elMat%Aql)) THEN
      DEALLOCATE (elMat%Aql)
    END IF
    IF (ALLOCATED(elMat%Auq)) THEN
      DEALLOCATE (elMat%Auq)
    END IF
    IF (ALLOCATED(elMat%Auu)) THEN
      DEALLOCATE (elMat%Auu)
    END IF
    IF (ALLOCATED(elMat%Aul)) THEN
      DEALLOCATE (elMat%Aul)
    END IF
    IF (ALLOCATED(elMat%Alq)) THEN
      DEALLOCATE (elMat%Alq)
    END IF
    IF (ALLOCATED(elMat%Alu)) THEN
      DEALLOCATE (elMat%Alu)
    END IF
    IF (ALLOCATED(elMat%All)) THEN
      DEALLOCATE (elMat%All)
    END IF
    IF (ALLOCATED(elMat%Aql_dir)) THEN
      DEALLOCATE (elMat%Aql_dir)
    END IF
    IF (ALLOCATED(elMat%Aul_dir)) THEN
      DEALLOCATE (elMat%Aul_dir)
    END IF
    IF (ALLOCATED(elMat%fH)) THEN
      DEALLOCATE (elMat%fH)
    END IF
    IF (ALLOCATED(elMat%LL)) THEN
      DEALLOCATE (elMat%LL)
    END IF
    IF (ALLOCATED(elMat%L0)) THEN
      DEALLOCATE (elMat%L0)
    ENDIF
    IF (ALLOCATED(elMat%S)) THEN
      DEALLOCATE (elMat%S)
    END IF
    IF (ALLOCATED(elMat%UU)) THEN
      DEALLOCATE (elMat%UU)
    END IF
    IF (ALLOCATED(elMat%U0)) THEN
      DEALLOCATE (elMat%U0)
    END IF

    ! MatK
    IF (ASSOCIATED(MatK%rowptr)) THEN
      DEALLOCATE (MatK%rowptr)
    END IF
    IF (ASSOCIATED(MatK%cols)) THEN
      DEALLOCATE (MatK%cols)
    END IF
    IF (ASSOCIATED(MatK%vals)) THEN
      DEALLOCATE (MatK%vals)
    END IF
    IF (ASSOCIATED(MatK%loc2glob)) THEN
      DEALLOCATE (MatK%loc2glob)
    END IF

    ! RHS
    IF (ASSOCIATED(RHS%vals)) THEN
      DEALLOCATE (RHS%vals)
    END IF
    IF (ASSOCIATED(RHS%loc2glob)) THEN
      DEALLOCATE (RHS%loc2glob)
    END IF

    ! Physics
    IF (ASSOCIATED(phys%B)) THEN
      DEALLOCATE (phys%B)
    END IF
    IF (ASSOCIATED(phys%magnetic_flux)) THEN
      DEALLOCATE (phys%magnetic_flux)
    END IF
    IF (ASSOCIATED(phys%magnetic_psi)) THEN
      DEALLOCATE (phys%magnetic_psi)
    END IF
    IF (ASSOCIATED(phys%Bperturb)) THEN
      DEALLOCATE (phys%Bperturb)
    END IF

    ! magnetic
    IF (ASSOCIATED(magn%coils_rmp)) THEN
      DEALLOCATE (magn%coils_rmp)
    END IF
    IF (ASSOCIATED(magn%coils_ripple)) THEN
      DEALLOCATE (magn%coils_ripple)
    END IF

    IF (ALLOCATED(simpar%physvar_refval)) THEN
      DEALLOCATE (simpar%physvar_refval)
    END IF
    IF (ALLOCATED(simpar%consvar_refval)) THEN
      DEALLOCATE (simpar%consvar_refval)
    END IF

  END SUBROUTINE free_all
  
  SUBROUTINE free_mesh
    ! Mesh
    IF (ASSOCIATED(Mesh%T)) THEN
      DEALLOCATE (Mesh%T)
    END IF
    IF (ASSOCIATED(Mesh%Tlin)) THEN
      DEALLOCATE (Mesh%Tlin)
    END IF
    IF (ASSOCIATED(Mesh%Tb)) THEN
      DEALLOCATE (Mesh%Tb)
    END IF
    IF (ASSOCIATED(Mesh%boundaryFlag)) THEN
      DEALLOCATE (Mesh%boundaryFlag)
    END IF
    IF (ALLOCATED(Mesh%F)) THEN
      DEALLOCATE (Mesh%F)
    END IF
    IF (ALLOCATED(Mesh%N)) THEN
      DEALLOCATE (Mesh%N)
    END IF
    IF (ALLOCATED(Mesh%faces)) THEN
      DEALLOCATE (Mesh%faces)
    END IF
    IF (ALLOCATED(Mesh%extfaces)) THEN
      DEALLOCATE (Mesh%extfaces)
    END IF
    IF (ALLOCATED(Mesh%intfaces)) THEN
      DEALLOCATE (Mesh%intfaces)
    END IF
    IF (ALLOCATED(Mesh%flipFace)) THEN
      DEALLOCATE (Mesh%flipface)
    END IF
    IF (ALLOCATED(Mesh%Fdir)) THEN
      DEALLOCATE (Mesh%Fdir)
    END IF
    IF (ALLOCATED(Mesh%Diric)) THEN
      DEALLOCATE (Mesh%Diric)
    END IF
    IF (ALLOCATED(Mesh%numberbcs)) THEN
      DEALLOCATE (Mesh%numberbcs)
    END IF
    IF (ASSOCIATED(Mesh%X)) THEN
      DEALLOCATE (Mesh%X)
    END IF
    IF (ALLOCATED(Mesh%elemSize)) THEN
      DEALLOCATE (Mesh%elemSize)
    END IF
    IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
      DEALLOCATE (Mesh%scdiff_nodes)
    END IF
    IF (ASSOCIATED(Mesh%toroidal)) THEN
      DEALLOCATE (Mesh%toroidal)
    END IF
    IF (ALLOCATED(Mesh%periodic_faces)) THEN
      DEALLOCATE (Mesh%periodic_faces)
    END IF  
  END SUBROUTINE free_mesh
  
  SUBROUTINE free_reference_element
     ! Reference element
    IF (ALLOCATED(refElPol%Face_nodes)) THEN
      DEALLOCATE (refElPol%Face_nodes)
    END IF
    IF (ALLOCATED(refElPol%Bord_nodes)) THEN
      DEALLOCATE (refElPol%Bord_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes)) THEN
      DEALLOCATE (refElPol%inner_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes_face)) THEN
      DEALLOCATE (refElPol%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElPol%coord3D)) THEN
      DEALLOCATE (refElPol%coord3D)
    END IF
    IF (ASSOCIATED(refElPol%coord2D)) THEN
      DEALLOCATE (refElPol%coord2D)
    END IF
    IF (ASSOCIATED(refElPol%coord1D)) THEN
      DEALLOCATE (refElPol%coord1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points3D)) THEN
      DEALLOCATE (refElPol%gauss_points3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points2D)) THEN
      DEALLOCATE (refElPol%gauss_points2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points1D)) THEN
      DEALLOCATE (refElPol%gauss_points1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights3D)) THEN
      DEALLOCATE (refElPol%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights2D)) THEN
      DEALLOCATE (refElPol%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights1D)) THEN
      DEALLOCATE (refElPol%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElPol%N3D)) THEN
      DEALLOCATE (refElPol%N3D)
    END IF
    IF (ALLOCATED(refElPol%Nxi3D)) THEN
      DEALLOCATE (refElPol%Nxi3D)
    END IF
    IF (ALLOCATED(refElPol%Neta3D)) THEN
      DEALLOCATE (refElPol%Neta3D)
    END IF
    IF (ALLOCATED(refElPol%Nzeta3D)) THEN
      DEALLOCATE (refElPol%Nzeta3D)
    END IF
    IF (ALLOCATED(refElPol%N2D)) THEN
      DEALLOCATE (refElPol%N2D)
    END IF
    IF (ALLOCATED(refElPol%Nxi2D)) THEN
      DEALLOCATE (refElPol%Nxi2D)
    END IF
    IF (ALLOCATED(refElPol%Neta2D)) THEN
      DEALLOCATE (refElPol%Neta2D)
    END IF
    IF (ALLOCATED(refElPol%N1D)) THEN
      DEALLOCATE (refElPol%N1D)
    END IF
    IF (ALLOCATED(refElPol%Nxi1D)) THEN
      DEALLOCATE (refElPol%Nxi1D)
    END IF

    ! Reference element
    IF (ALLOCATED(refElTor%Face_nodes)) THEN
      DEALLOCATE (refElTor%Face_nodes)
    END IF
    IF (ALLOCATED(refElTor%Bord_nodes)) THEN
      DEALLOCATE (refElTor%Bord_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes)) THEN
      DEALLOCATE (refElTor%inner_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes_face)) THEN
      DEALLOCATE (refElTor%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElTor%coord3D)) THEN
      DEALLOCATE (refElTor%coord3D)
    END IF
    IF (ASSOCIATED(refElTor%coord2D)) THEN
      DEALLOCATE (refElTor%coord2D)
    END IF
    IF (ASSOCIATED(refElTor%coord1D)) THEN
      DEALLOCATE (refElTor%coord1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points3D)) THEN
      DEALLOCATE (refElTor%gauss_points3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points2D)) THEN
      DEALLOCATE (refElTor%gauss_points2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points1D)) THEN
      DEALLOCATE (refElTor%gauss_points1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights3D)) THEN
      DEALLOCATE (refElTor%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights2D)) THEN
      DEALLOCATE (refElTor%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights1D)) THEN
      DEALLOCATE (refElTor%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElTor%N3D)) THEN
      DEALLOCATE (refElTor%N3D)
    END IF
    IF (ALLOCATED(refElTor%Nxi3D)) THEN
      DEALLOCATE (refElTor%Nxi3D)
    END IF
    IF (ALLOCATED(refElTor%Neta3D)) THEN
      DEALLOCATE (refElTor%Neta3D)
    END IF
    IF (ALLOCATED(refElTor%Nzeta3D)) THEN
      DEALLOCATE (refElTor%Nzeta3D)
    END IF
    IF (ALLOCATED(refElTor%N2D)) THEN
      DEALLOCATE (refElTor%N2D)
    END IF
    IF (ALLOCATED(refElTor%Nxi2D)) THEN
      DEALLOCATE (refElTor%Nxi2D)
    END IF
    IF (ALLOCATED(refElTor%Neta2D)) THEN
      DEALLOCATE (refElTor%Neta2D)
    END IF
    IF (ALLOCATED(refElTor%N1D)) THEN
      DEALLOCATE (refElTor%N1D)
    END IF
    IF (ALLOCATED(refElTor%Nxi1D)) THEN
      DEALLOCATE (refElTor%Nxi1D)
    ENDIF 
  END SUBROUTINE free_reference_element  
  
  SUBROUTINE free_mat
    DEALLOCATE (MatK%cols)
    DEALLOCATE (MatK%rowptr)
    DEALLOCATE (MatK%vals)
    DEALLOCATE (MatK%loc2glob)
    DEALLOCATE (rhs%loc2glob)
    DEALLOCATE (rhs%vals)
  END SUBROUTINE

END MODULE globals
