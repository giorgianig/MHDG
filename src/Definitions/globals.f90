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
   TYPE(Mesh_type)             , target :: mesh
   TYPE(Physics_type)          , target :: phys
   TYPE(Geometry_type)         , target :: geom
   TYPE(Switches_type)         , target :: switch      
   TYPE(Time_type)             , target :: time
   TYPE(Numeric_type)          , target :: numer
   TYPE(Utils_type)            , target :: utils      
   TYPE(Lssolver_type)         , target :: lssolver
   TYPE(Elmat_type)            , target :: elmat   
   TYPE(Sol_type)              , target :: sol
   TYPE(MAT_CSR_TYP)           , target :: MatK
   TYPE(RHS_TYP)               , target :: rhs
  


   CONTAINS
   
  
   
   !***********************************************
   ! Routine to deallocate all previously allocated
   ! structures
   !***********************************************
   SUBROUTINE free_all()
   
						! Reference element
						IF (ALLOCATED(refElPol%Face_nodes)) THEN
						  DEALLOCATE(refElPol%Face_nodes)
						END IF
						IF (ALLOCATED(refElPol%Bord_nodes)) THEN
						  DEALLOCATE(refElPol%Bord_nodes)
						END IF   
						IF (ALLOCATED(refElPol%inner_nodes)) THEN
						  DEALLOCATE(refElPol%inner_nodes)
						END IF 
						IF (ALLOCATED(refElPol%inner_nodes_face)) THEN
						  DEALLOCATE(refElPol%inner_nodes_face)
						END IF     
						IF (ASSOCIATED(refElPol%coord3D)) THEN
						  DEALLOCATE(refElPol%coord3D)
						END IF
						IF (ASSOCIATED(refElPol%coord2D)) THEN
						  DEALLOCATE(refElPol%coord2D)
						END IF      
						IF (ASSOCIATED(refElPol%coord1D)) THEN
						  DEALLOCATE(refElPol%coord1D)   
      END IF
						IF (ALLOCATED(refElPol%gauss_points3D)) THEN
						  DEALLOCATE(refElPol%gauss_points3D)
						END IF
						IF (ALLOCATED(refElPol%gauss_points2D)) THEN
						  DEALLOCATE(refElPol%gauss_points2D)
						END IF      
						IF (ALLOCATED(refElPol%gauss_points1D)) THEN
						  DEALLOCATE(refElPol%gauss_points1D)
						END IF      
						IF (ALLOCATED(refElPol%gauss_weights3D)) THEN
						  DEALLOCATE(refElPol%gauss_weights3D)
						END IF  
						IF (ALLOCATED(refElPol%gauss_weights2D)) THEN
						  DEALLOCATE(refElPol%gauss_weights2D)
						END IF  
						IF (ALLOCATED(refElPol%gauss_weights1D)) THEN
						  DEALLOCATE(refElPol%gauss_weights1D)
						END IF           
						IF (ALLOCATED(refElPol%N3D)) THEN
						  DEALLOCATE(refElPol%N3D)
						END IF 
						IF (ALLOCATED(refElPol%Nxi3D)) THEN
						  DEALLOCATE(refElPol%Nxi3D)
						END IF 
						IF (ALLOCATED(refElPol%Neta3D)) THEN
						  DEALLOCATE(refElPol%Neta3D)
						END IF   
						IF (ALLOCATED(refElPol%Nzeta3D)) THEN
						  DEALLOCATE(refElPol%Nzeta3D)
						END IF        
						IF (ALLOCATED(refElPol%N2D)) THEN
						  DEALLOCATE(refElPol%N2D)
						END IF 
						IF (ALLOCATED(refElPol%Nxi2D)) THEN
						  DEALLOCATE(refElPol%Nxi2D)
						END IF  
						IF (ALLOCATED(refElPol%Neta2D)) THEN
						  DEALLOCATE(refElPol%Neta2D)
						END IF 
						IF (ALLOCATED(refElPol%N1D)) THEN
						  DEALLOCATE(refElPol%N1D)
						END IF     
						IF (ALLOCATED(refElPol%Nxi1D)) THEN
						  DEALLOCATE(refElPol%Nxi1D)
						END IF                    
						
						
						
					! Reference element
						IF (ALLOCATED(refElTor%Face_nodes)) THEN
						  DEALLOCATE(refElTor%Face_nodes)
						END IF
						IF (ALLOCATED(refElTor%Bord_nodes)) THEN
						  DEALLOCATE(refElTor%Bord_nodes)
						END IF   
						IF (ALLOCATED(refElTor%inner_nodes)) THEN
						  DEALLOCATE(refElTor%inner_nodes)
						END IF 
						IF (ALLOCATED(refElTor%inner_nodes_face)) THEN
						  DEALLOCATE(refElTor%inner_nodes_face)
						END IF     
						IF (ASSOCIATED(refElTor%coord3D)) THEN
						  DEALLOCATE(refElTor%coord3D)
						END IF
						IF (ASSOCIATED(refElTor%coord2D)) THEN
						  DEALLOCATE(refElTor%coord2D)
						END IF      
						IF (ASSOCIATED(refElTor%coord1D)) THEN
						  DEALLOCATE(refElTor%coord1D)   
      END IF
						IF (ALLOCATED(refElTor%gauss_points3D)) THEN
						  DEALLOCATE(refElTor%gauss_points3D)
						END IF
						IF (ALLOCATED(refElTor%gauss_points2D)) THEN
						  DEALLOCATE(refElTor%gauss_points2D)
						END IF      
						IF (ALLOCATED(refElTor%gauss_points1D)) THEN
						  DEALLOCATE(refElTor%gauss_points1D)
						END IF      
						IF (ALLOCATED(refElTor%gauss_weights3D)) THEN
						  DEALLOCATE(refElTor%gauss_weights3D)
						END IF  
						IF (ALLOCATED(refElTor%gauss_weights2D)) THEN
						  DEALLOCATE(refElTor%gauss_weights2D)
						END IF  
						IF (ALLOCATED(refElTor%gauss_weights1D)) THEN
						  DEALLOCATE(refElTor%gauss_weights1D)
						END IF           
						IF (ALLOCATED(refElTor%N3D)) THEN
						  DEALLOCATE(refElTor%N3D)
						END IF 
						IF (ALLOCATED(refElTor%Nxi3D)) THEN
						  DEALLOCATE(refElTor%Nxi3D)
						END IF 
						IF (ALLOCATED(refElTor%Neta3D)) THEN
						  DEALLOCATE(refElTor%Neta3D)
						END IF   
						IF (ALLOCATED(refElTor%Nzeta3D)) THEN
						  DEALLOCATE(refElTor%Nzeta3D)
						END IF        
						IF (ALLOCATED(refElTor%N2D)) THEN
						  DEALLOCATE(refElTor%N2D)
						END IF 
						IF (ALLOCATED(refElTor%Nxi2D)) THEN
						  DEALLOCATE(refElTor%Nxi2D)
						END IF  
						IF (ALLOCATED(refElTor%Neta2D)) THEN
						  DEALLOCATE(refElTor%Neta2D)
						END IF 
						IF (ALLOCATED(refElTor%N1D)) THEN
						  DEALLOCATE(refElTor%N1D)
						END IF     
						IF (ALLOCATED(refElTor%Nxi1D)) THEN
						  DEALLOCATE(refElTor%Nxi1D)
					 ENDIF
						
						! Mesh
						IF (ASSOCIATED(Mesh%T)) THEN
						   DEALLOCATE(Mesh%T) 
						END IF
						IF (ASSOCIATED(Mesh%Tlin)) THEN
						   DEALLOCATE(Mesh%Tlin) 
						END IF
						IF (ASSOCIATED(Mesh%Tb)) THEN
						   DEALLOCATE(Mesh%Tb) 
						END IF      
						IF (ASSOCIATED(Mesh%boundaryFlag)) THEN
						   DEALLOCATE(Mesh%boundaryFlag) 
						END IF
						IF (ALLOCATED(Mesh%F)) THEN
						   DEALLOCATE(Mesh%F) 
						END IF   
						IF (ALLOCATED(Mesh%N)) THEN
						   DEALLOCATE(Mesh%N) 
						END IF   
						IF (ALLOCATED(Mesh%faces)) THEN
						   DEALLOCATE(Mesh%faces) 
						END IF   
						IF (ALLOCATED(Mesh%extfaces)) THEN
						   DEALLOCATE(Mesh%extfaces) 
						END IF   
						IF (ALLOCATED(Mesh%intfaces)) THEN
						   DEALLOCATE(Mesh%intfaces) 
						END IF   
						IF (ALLOCATED(Mesh%flipFace)) THEN
						   DEALLOCATE(Mesh%flipface) 
						END IF   
						IF (ALLOCATED(Mesh%Fdir)) THEN
						   DEALLOCATE(Mesh%Fdir) 
						END IF   
						IF (ALLOCATED(Mesh%Diric)) THEN
						   DEALLOCATE(Mesh%Diric) 
						END IF   
						IF (ALLOCATED(Mesh%numberbcs)) THEN
						   DEALLOCATE(Mesh%numberbcs) 
						END IF   
						IF (ASSOCIATED(Mesh%X)) THEN
						   DEALLOCATE(Mesh%X) 
						END IF   
						IF (ALLOCATED(Mesh%elemSize)) THEN
						   DEALLOCATE(Mesh%elemSize) 
						END IF
						IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
						   DEALLOCATE(Mesh%scdiff_nodes) 
						END IF												
						! physics
!						IF (ASSOCIATED(phys%phyVarNam)) THEN
!						   DEALLOCATE(phys%phyVarNam) 
!						END IF  
!						IF (ASSOCIATED(phys%conVarNam)) THEN
!						   DEALLOCATE(phys%conVarNam) 
!						END IF   
						     

						! sol type
						IF (ASSOCIATED(sol%u)) THEN
						   DEALLOCATE(sol%u) 
						END IF 
						IF (ASSOCIATED(sol%u_tilde)) THEN
						   DEALLOCATE(sol%u_tilde) 
						END IF    
						IF (ALLOCATED(sol%tres)) THEN
						   DEALLOCATE(sol%tres) 
						END IF       
						IF (ALLOCATED(sol%time)) THEN
						   DEALLOCATE(sol%time) 
						END IF   						
						
						! elemental matrices
						IF (ALLOCATED(elMat%M)) THEN
						   DEALLOCATE(elMat%M) 
						END IF       
						 IF (ALLOCATED(elMat%Cv)) THEN
						   DEALLOCATE(elMat%Cv) 
						END IF       
						 IF (ALLOCATED(elMat%H)) THEN
						   DEALLOCATE(elMat%H) 
						END IF       
						 IF (ALLOCATED(elMat%Hdir)) THEN
						   DEALLOCATE(elMat%Hdir) 
						END IF       
						 IF (ALLOCATED(elMat%D)) THEN
						   DEALLOCATE(elMat%D) 
						END IF       
						 IF (ALLOCATED(elMat%E)) THEN
						   DEALLOCATE(elMat%E) 
						END IF       
						 IF (ALLOCATED(elMat%Edir)) THEN
						   DEALLOCATE(elMat%Edir) 
						END IF       
						 IF (ALLOCATED(elMat%S)) THEN
						   DEALLOCATE(elMat%S) 
						END IF       
						IF (ALLOCATED(elMat%UU)) THEN
						   DEALLOCATE(elMat%UU) 
						END IF       
						 IF (ALLOCATED(elMat%U0)) THEN
						   DEALLOCATE(elMat%U0) 
						END IF       
						 IF (ALLOCATED(elMat%Hf)) THEN
						   DEALLOCATE(elMat%Hf) 
						END IF       
						IF (ALLOCATED(elMat%Df)) THEN
						   DEALLOCATE(elMat%Df) 
						END IF       
						IF (ALLOCATED(elMat%Ef)) THEN
						   DEALLOCATE(elMat%Ef) 
						END IF       
						IF (ALLOCATED(elMat%fH)) THEN
						   DEALLOCATE(elMat%fH) 
						END IF       
						IF (ALLOCATED(elMat%LL)) THEN
						   DEALLOCATE(elMat%LL) 
						END IF       
						 IF (ALLOCATED(elMat%L0)) THEN
						   DEALLOCATE(elMat%L0) 
						END IF       
						 IF (ALLOCATED(elMat%Lf)) THEN
						   DEALLOCATE(elMat%Lf) 
						END IF       
						 IF (ALLOCATED(elMat%Qf)) THEN
						   DEALLOCATE(elMat%Qf) 
						END IF  
						 IF (ALLOCATED(elMat%B)) THEN
						   DEALLOCATE(elMat%B) 
						END IF 						
						 IF (ALLOCATED(elMat%C)) THEN
						   DEALLOCATE(elMat%C) 
						END IF 						
						 IF (ALLOCATED(elMat%Cdir)) THEN
						   DEALLOCATE(elMat%Cdir) 
						END IF 						
						 IF (ALLOCATED(elMat%P)) THEN
						   DEALLOCATE(elMat%P) 
						END IF 	
						 IF (ALLOCATED(elMat%G)) THEN
						   DEALLOCATE(elMat%G) 
						END IF 												
						 IF (ALLOCATED(elMat%iL)) THEN
						   DEALLOCATE(elMat%iL) 
						END IF 												
						 IF (ALLOCATED(elMat%Lf)) THEN
						   DEALLOCATE(elMat%Lf) 
						END IF
						 IF (ALLOCATED(elMat%S_lrho)) THEN
						   DEALLOCATE(elMat%S_lrho) 
						END IF						 		
						 IF (ALLOCATED(elMat%P_lrho)) THEN
						   DEALLOCATE(elMat%P_lrho) 
						END IF							
						 IF (ALLOCATED(elMat%P_sc)) THEN
						   DEALLOCATE(elMat%P_sc) 
						END IF	
						 IF (ALLOCATED(elMat%Lf_sc)) THEN
						   DEALLOCATE(elMat%Lf_sc) 
						END IF								
#ifdef TEMPERATURE
						 IF (ALLOCATED(elMat%TQhf)) THEN
						   DEALLOCATE(elMat%TQhf) 
						END IF		
						 IF (ALLOCATED(elMat%TQ)) THEN
						   DEALLOCATE(elMat%TQ) 
						END IF		
						 IF (ALLOCATED(elMat%Tf)) THEN
						   DEALLOCATE(elMat%Tf) 
						END IF		
						 IF (ALLOCATED(elMat%Tfhf)) THEN
						   DEALLOCATE(elMat%Tfhf) 
						END IF																		
						 IF (ALLOCATED(elMat%Thdir)) THEN
						   DEALLOCATE(elMat%Thdir) 
						END IF																										
#endif											
						! MatK
						IF (ASSOCIATED(MatK%rowptr)) THEN
						   DEALLOCATE(MatK%rowptr) 
						END IF 
						IF (ASSOCIATED(MatK%cols)) THEN
						   DEALLOCATE(MatK%cols) 
						END IF 
						IF (ASSOCIATED(MatK%vals)) THEN
						   DEALLOCATE(MatK%vals) 
						END IF 
						IF (ASSOCIATED(MatK%loc2glob)) THEN
						   DEALLOCATE(MatK%loc2glob) 
						END IF 
						
						! RHS
						IF (ASSOCIATED(RHS%vals)) THEN
						   DEALLOCATE(RHS%vals) 
						END IF 
						IF (ASSOCIATED(RHS%loc2glob)) THEN
						   DEALLOCATE(RHS%loc2glob) 
						END IF 					  						
 					
    
   END SUBROUTINE free_all
   
   SUBROUTINE free_mat
			 DEALLOCATE(MatK%cols)
				DEALLOCATE(MatK%rowptr)
				DEALLOCATE(MatK%vals)
				DEALLOCATE(MatK%loc2glob)
				DEALLOCATE(rhs%loc2glob)
				DEALLOCATE(rhs%vals)
   END SUBROUTINE
    
END MODULE globals      
