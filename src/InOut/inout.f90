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
   character(70)  :: npr,nid
   real*8,parameter::tol=1e-6

   integer :: elemType, ndim, Nnodes, Nelems, Nnodesperelem, Nfaces
   integer :: Nextfaces, Nnodesperface, IERR
   integer(HID_T) :: file_id
#ifdef PARALL
   integer :: ghfa,ghel,i,Nel_glob,Nfa_glob,Ndir_glob,Ngho_glob
#endif
#ifdef TOR3D
#ifdef PARALL
   IF (MPIvar%npol.GT.1) THEN
      write(nid,*) MPIvar%ipol
      write(npr,*) MPIvar%npol
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF
#else
   fname_complete = trim(adjustl(fname))//'.h5' 
#endif
#else
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF
#endif
   IF (utils%printint>0) THEN
      print*,'Loading mesh.'
      print*,'	'
   ENDIF

      CALL HDF5_open(fname_complete,file_id,IERR)     
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error opening mesh file: ", fname_complete
         STOP
      ENDIF
      CALL HDF5_integer_reading(file_id,elemType,'elemType',ierr) 
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: elemType" 
         STOP
      ENDIF      
      CALL HDF5_integer_reading(file_id,ndim,'Ndim',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Ndim" 
         STOP
      ENDIF       
      CALL HDF5_integer_reading(file_id,Nnodes,'Nnodes',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nnodes" 
         STOP
      ENDIF       
      CALL HDF5_integer_reading(file_id,Nelems,'Nelems',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nelems" 
         STOP
      ENDIF       
      CALL HDF5_integer_reading(file_id,Nnodesperelem,'Nnodesperelem',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nnodesperelem" 
         STOP
      ENDIF       
      CALL HDF5_integer_reading(file_id,Nnodesperface,'Nnodesperface',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nnodesperface" 
         STOP
      ENDIF       
      CALL HDF5_integer_reading(file_id,Nextfaces,'Nextfaces',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nextfaces" 
         STOP
      ENDIF       
#ifdef PARALL 
      CALL HDF5_integer_reading(file_id,Nfaces,'Nfaces',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading integer: Nfaces" 
         STOP
      ENDIF      
#endif      
      ALLOCATE(Mesh%T(Nelems,Nnodesperelem))
      ALLOCATE(Mesh%X(Nnodes,ndim))
      ALLOCATE(Mesh%Tb(Nextfaces,Nnodesperface))
      ALLOCATE(Mesh%boundaryFlag(Nextfaces))   
#ifdef PARALL 
      ALLOCATE(Mesh%ghostFaces(Nfaces))
      ALLOCATE(Mesh%loc2glob_fa(Nfaces))
      ALLOCATE(Mesh%loc2glob_el(Nelems))
#ifdef TOR3D
      ALLOCATE(Mesh%ghostElems(Nelems))
#endif      
#endif      
      CALL HDF5_array2D_reading_int(file_id,Mesh%T,'T',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading mesh connectivity T" 
         STOP
      ENDIF      
      CALL HDF5_array2D_reading_int(file_id,Mesh%Tb,'Tb',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading boundary connectivity Tb" 
         STOP
      ENDIF       
      CALL HDF5_array1D_reading_int(file_id,Mesh%boundaryFlag,'boundaryFlag',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading boundaryFlag" 
         STOP
      ENDIF               
      CALL HDF5_array2D_reading(file_id,Mesh%X,'X',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading coordinate matrix X" 
         STOP
      ENDIF       
#ifdef PARALL 
      CALL HDF5_array1D_reading_int(file_id,Mesh%loc2glob_fa,'loc2glob_fa',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading loc2glob_fa" 
         STOP
      ENDIF       
      CALL HDF5_array1D_reading_int(file_id,Mesh%loc2glob_el,'loc2glob_el',ierr)   
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading loc2glob_el" 
         STOP
      ENDIF                
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghostFaces,'ghostFaces',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghostFaces" 
         STOP
      ENDIF         
#ifdef TOR3D
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghostElems,'ghostElems',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghostElems" 
         STOP
      ENDIF         
#endif      
      ! Find the number of ghost faces
      ghfa = sum(Mesh%ghostFaces)
      Mesh%nghostfaces = ghfa      
#ifdef TOR3D      
      ! Find the number of ghost elements
      ghel = sum(Mesh%ghostElems)    
      Mesh%nghostElems = ghel
#endif
      ALLOCATE(Mesh%ghostflp(ghfa))
      ALLOCATE(Mesh%ghostpro(ghfa))
      ALLOCATE(Mesh%ghostloc(ghfa))
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghostflp,'ghostFlp',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghostFlp" 
         STOP
      ENDIF       
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghostLoc,'ghostLoc',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghostLoc" 
         STOP
      ENDIF       
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghostPro,'ghostPro',ierr) 
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghostPro" 
         STOP
      ENDIF       
#ifdef TOR3D
      ALLOCATE(Mesh%ghelspro(ghfa))
      ALLOCATE(Mesh%ghelsloc(ghfa))
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghelsLoc,'ghelsLoc',ierr)
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghelsLoc" 
         STOP
      ENDIF       
      CALL HDF5_array1D_reading_int(file_id,Mesh%ghelsPro,'ghelsPro',ierr) 
      IF (IERR .ne. 0) THEN
         WRITE(6,*) "Error reading ghelsPro" 
         STOP
      ENDIF       
#endif 
#endif      
      CALL HDF5_close(file_id)  
!************************************************************************      
!   CONFIRMATION MESSAGE FOR THE USER
!************************************************************************
#ifdef PARALL      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      WRITE(6,*) "Process: ", MPIvar%glob_id, "-- readed mesh file: ",trim(adjustl(fname_complete))
#else
      WRITE(6,*) "Readed mesh file: ",trim(adjustl(fname_complete))
#endif      

#ifdef PARALL
      CALL MPI_ALLREDUCE(maxval(Mesh%loc2glob_el),Nel_glob,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(maxval(Mesh%loc2glob_fa),Nfa_glob,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(Mesh%ndir,Ndir_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(Mesh%nghostfaces,Ngho_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      Mesh%Nel_glob = Nel_glob
      Mesh%Nfa_glob = Nfa_glob
      Mesh%Ndir_glob = Ndir_glob
      Mesh%Ngho_glob = Ngho_glob
#endif      
      Mesh%Ndim   = ndim
      Mesh%Nnodes = Nnodes
      Mesh%Nelems = Nelems
      Mesh%Nnodesperelem = Nnodesperelem
      Mesh%Nnodesperface = Nnodesperface
      Mesh%elemType = elemType
      Mesh%Nextfaces = Nextfaces
   
   ! Apply shift if axisymmetric case
   IF ( (switch%axisym .and. switch%testcase.ge.60) .or. (switch%axisym .and. minval(Mesh%X(:,1))<tol )) THEN
      IF (MPIvar%glob_id.eq.0) THEN
         WRITE(6,*) "*** Applying translation in axisymmetric case!"
      ENDIF
      Mesh%X(:,1) = Mesh%X(:,1) + geom%R0
   END IF
   
   ! Apply length scale
   Mesh%X = Mesh%X/phys%lscale
   
   IF (utils%printint>0) then
     IF (MPIvar%glob_id.eq.0) THEN
        IF (elemType==0) then
         WRITE(str,'(A)') 'triangles'
         ELSEIF (elemType==1) then
         WRITE(str,'(A)') 'quads'
         ELSEIF (elemType==2) then
         WRITE(str,'(A)')'thetra'
         ELSEIF (elemType==3) then
         WRITE(str,'(A)') 'hexa'
         end IF         
         WRITE(6,*) '*************************************************'
         WRITE(6,*) '*                    MESH                       *'
         WRITE(6,*) '*************************************************'
         WRITE(6,'(A,I18)')  ' Number of dimensions:         ', ndim
         WRITE(6,'(A,A34)')  ' Element type: ', trim(str)
         WRITE(6,'(A,I18)')  ' Number of elements:           ', Nelems
         WRITE(6,'(A,I18)')  ' Number of nodes:              ', Nnodes
         WRITE(6,'(A,I18)')  ' Number of nodes per element:  ', Nnodesperelem
         WRITE(6,'(A,I18)')  ' Number of nodes per face:     ', Nnodesperface
         WRITE(6,'(A,I18)')  ' Number of exterior faces:     ', Nextfaces
         WRITE(6,*)  ' '
         WRITE(6,*)  ' '
      IF (utils%printint>1) THEN
         WRITE(6,*) "Connectivity matrix T:"
         CALL displayMatrixInt(Mesh%T)
         WRITE(6,*) "Boundary connectivity matrix Tb:"
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
   character(70)  :: npr,nid
#ifdef TOR3D
    character(70)  :: nip,nit,ngd
#endif   
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   
   
#ifdef TOR3D
   IF (MPIvar%glob_size.GT.1) THEN
      write(nip,*) MPIvar%ipol
      write(nit,*) MPIvar%itor
      write(ngd,*) MPIvar%glob_size   
      fname_complete = trim(adjustl(fname))//'_ip'//trim(adjustl(nip))//'_it'//trim(adjustl(nit))//'_np'//trim(adjustl(ngd))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF  
#else   
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF      
#endif   
      call HDF5_create(fname_complete,file_id,ierr)
      call HDF5_array1D_saving(file_id,sol%u,size(sol%u),'u')
      call HDF5_array1D_saving(file_id,sol%u_tilde,size(sol%u_tilde),'u_tilde')
!      call HDF5_array1D_saving(file_id,sol%tres(1:sol%Nt),sol%Nt,'tres')
!      call HDF5_array1D_saving(file_id,sol%time(1:sol%Nt),sol%Nt,'time')
!      call HDF5_integer_saving(file_id,sol%Nt,'Nt')
#ifdef TEMPERATURE      
      call HDF5_array1D_saving(file_id,sol%q,size(sol%q),'q')
#endif      
      IF (switch%shockcp .eq.3) THEN
         call HDF5_array2D_saving(file_id,Mesh%scdiff_nodes,size(Mesh%scdiff_nodes,1),size(Mesh%scdiff_nodes,2),'scdiff_nodes')
      END IF
      call HDF5_close(file_id)
      ! Message to confirm succesful creation and filling of file
      IF (MPIvar%glob_id.eq.0) THEN
						   print*,'Output written to file ', trim(adjustl(fname_complete))
						   print*,'	'
      END IF
   
end subroutine HDF5_save_solution


!**********************************************************************
! Load solution in HDF5 file format
!**********************************************************************
subroutine HDF5_load_solution(fname)
   USE globals
   implicit none
   
   character(LEN=*) :: fname
   character(70)  :: npr,nid
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   real,pointer :: oldtres(:),oldtime(:)
   real,pointer:: u_aux(:),u_tilde_aux(:),q_aux(:)
   
   integer :: Neq,Ndim,Ntor,N2d,Nel,Np1d,Np2d,Np,Nfl,Nfp,Nfg,Nf,sizeutilde,sizeu
   integer :: iel,ifa,iface,Fi,itor,dd,delta,i,it,iel3
   real,allocatable    :: u2D(:),q2D(:),ut2D(:)
   integer,allocatable :: indu2D(:),indq2D(:),indu3D(:),indq3D(:),indufp(:),induf2D(:),induf3D(:),indul(:),indql(:),indutl(:)
    Neq         = phys%Neq   
#ifdef TOR3D     
    Ndim        = 3                             ! N. of dimensions
    Ntor        = numer%ntor                    ! N. of toroidal divisions
    N2d         = Mesh%Nelems                   ! N. of 2D elements
    Nel         = N2d*ntor                      ! N. of 3D elements
    Np1d        = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element  
    Np2d        = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np          = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl         = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg         = Np2d*2+refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf          = Mesh%Nfaces                   ! N. of faces in the 2D mesh
    sizeu       = Neq*Nel*Np                    ! Size of u    
    sizeutilde  = Neq*numer%ntor*(Nfl*Nf + Np2d*N2d)! Size of utilde    
#else				 
    Ndim        = 2
    Nel         = Mesh%Nelems
    Np          = refElPol%Nnodes2D
    Nf          = refElPol%Nfaces 
    Nfg         = refElPol%Nfacenodes*Nf 
    sizeu       = Neq*Nel*Np
    sizeutilde  = Neq*Mesh%Nfaces*Mesh%Nnodesperface
#endif 
  
  		  ALLOCATE(sol%u(sizeu))
		    ALLOCATE(sol%u_tilde(sizeutilde))
#ifdef TEMPERATURE  		    
      ALLOCATE(sol%q(sizeu*Ndim))		    
#endif
						IF (MPIvar%glob_size.GT.1) THEN
						   write(nid,*) MPIvar%glob_id+1
						   write(npr,*) MPIvar%glob_size
						   fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
						ELSE
						   fname_complete = trim(adjustl(fname))//'.h5'
						END IF      		    
      
#ifdef TOR3D     
!*************************************
!              3D case
!*************************************
 
      IF (fname_complete(1:5)=='Sol3D') THEN
! Initialization with a 3D solution      
         WRITE(6,*) "3D initial solution"
         CALL HDF5_open(fname_complete,file_id,IERR)
         CALL HDF5_array1D_reading(file_id,sol%u,'u')
         CALL HDF5_array1D_reading(file_id,sol%u_tilde,'u_tilde')
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
#ifdef TEMPERATURE      
         CALL HDF5_array1D_reading(file_id,sol%q,'q')
#endif       
         CALL HDF5_close(file_id)
      ELSEIF (fname_complete(1:5)=='Sol2D') THEN
! Initialization with a 2D solution      
         WRITE(6,*) "2D initial solution: propagating in the torus..."
  		     ALLOCATE(u_aux(Neq*N2d*Np))
		       ALLOCATE(u_tilde_aux(Neq*Mesh%Nfaces*Mesh%Nnodesperface))
         ALLOCATE(u2D(Np2D*Neq))
         ALLOCATE(ut2D(refElPol%Nnodes1D*Neq))
         ALLOCATE(indu2D(Np2D*Neq))
         ALLOCATE(indu3D(Np*Neq))
         ALLOCATE(indufp(Np2D*Neq))
         ALLOCATE(induf2D(refElPol%Nnodes1D*Neq))        
         ALLOCATE(induf3D(Nfl*Neq)) 
         ALLOCATE(indul(Np2D*Neq))
         ALLOCATE(indutl(refElPol%Nnodes1D*Neq))        
        



#ifdef TEMPERATURE  		    
         ALLOCATE(q_aux(Neq*N2d*Np*Ndim))		  
         ALLOCATE(q2D(Np2D*Neq*Ndim))  
         ALLOCATE(indq2D(Np2D*Neq*Ndim))
         ALLOCATE(indq3D(Np*Neq*Ndim))
         ALLOCATE(indql(Np2D*Neq*Ndim))
#endif




         CALL HDF5_open(fname_complete,file_id,IERR)
         CALL HDF5_array1D_reading(file_id,u_aux,'u')
         CALL HDF5_array1D_reading(file_id,u_tilde_aux,'u_tilde')
#ifdef TEMPERATURE      
         CALL HDF5_array1D_reading(file_id,q_aux,'q')
#endif          
         CALL HDF5_close(file_id)         
         
         DO iel = 1,N2D
            indu2D = (iel-1)*Np2d*Neq + (/(i,i=1,Np2d*Neq)/)
            u2D    = u_aux(indu2D)
#ifdef TEMPERATURE
            indq2D = (iel-1)*Np2d*Neq*Ndim + (/(i,i=1,Np2d*Neq*Ndim)/)
            q2D    = q_aux(indq2D)
#endif            
            DO itor = 1,Ntor
               iel3   = (itor-1)*N2d+iel
               indu3D = (iel3-1)*Np*Neq + (/(i,i=1,Np*Neq)/)
#ifdef TEMPERATURE            
               indq3D = (iel3-1)*Np*Neq*Ndim + (/(i,i=1,Np*Neq*Ndim)/)
#endif            
               dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Mesh%Ndir)*Nfl)*Neq+(iel-1)*Np2D*Neq
               indufp = dd+(/(i,i=1,Np2D*Neq)/)
               sol%u_tilde(indufp) = u2d
               DO it = 1,Np1d
                  indul = (it-1)*Np2D*Neq + (/(i,i=1,Np2D*Neq)/)
                  sol%u(indu3D(indul)) = u_aux(indu2D)
#ifdef TEMPERATURE   
                  indql = (it-1)*Np2D*Neq*Ndim + (/(i,i=1,Np2D*Neq*Ndim)/)
                  sol%q(indq3D(indql)) = q_aux(indq2D)         
#endif               
               END DO
            END DO
         END DO

         DO iface = 1,Mesh%Nintfaces
            Fi = iface
            induf2D = (Fi-1)*refElPol%Nnodes1D*Neq+(/(i,i=1,refElPol%Nnodes1D*Neq)/)
            ut2d = u_tilde_aux(induf2D)
            DO itor = 1,Ntor
               dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Mesh%Ndir)*Nfl)*Neq+(N2D*Np2D+(Fi-1)*Nfl)*Neq
               induf3D = dd+(/(i,i=1,Nfl*Neq)/)
               DO it = 1,Np1d
                  indutl = (it-1)*refElPol%Nnodes1D*Neq + (/(i,i=1,refElPol%Nnodes1D*Neq)/)
                  sol%u_tilde(induf3D(indutl)) = ut2d
               END DO
            END DO
         END DO

         DO iface = 1,Mesh%Nextfaces
            iel = Mesh%extfaces(iface,1)
            ifa = Mesh%extfaces(iface,2)
            IF ( Mesh%Fdir(iel,ifa) ) CYCLE
            Fi = iface+Mesh%Nintfaces
            induf2D = (Fi-1)*refElPol%Nnodes1D*Neq+(/(i,i=1,refElPol%Nnodes1D*Neq)/)
            ut2d = u_tilde_aux(induf2D)
            DO itor = 1,Ntor
               dd = (itor-1)*(N2D*Np2D+(Mesh%Nfaces-Mesh%Ndir)*Nfl)*Neq+(N2D*Np2D+(Fi-1)*Nfl)*Neq
               induf3D = dd+(/(i,i=1,Nfl*Neq)/)
               DO it = 1,Np1d
                  indutl = (it-1)*refElPol%Nnodes1D*Neq + (/(i,i=1,refElPol%Nnodes1D*Neq)/)
                  sol%u_tilde(induf3D(indutl)) = ut2d
               END DO
            END DO
         END DO

         
          WRITE(6,*) "Done!"

         DEALLOCATE(u_aux,u_tilde_aux,u2D,ut2D,indu2D,indu3D,indufp,induf2D,induf3D,indul,indutl)      
#ifdef TEMPERATURE
         DEALLOCATE(q_aux,q2D,indq2D,indq3D,indql)
#endif
      END IF


#else
!*************************************
!              2D case
!*************************************
         CALL HDF5_open(fname_complete,file_id,IERR)
         CALL HDF5_array1D_reading(file_id,sol%u,'u')
         CALL HDF5_array1D_reading(file_id,sol%u_tilde,'u_tilde')
#ifdef TEMPERATURE      
         CALL HDF5_array1D_reading(file_id,sol%q,'q')
#endif
         CALL HDF5_close(file_id)
#endif         

      ! Message to confirm succesful reading of file
      IF (MPIvar%glob_id.eq.0) THEN
						   print*,'Solution read from file ', trim(adjustl(fname_complete))
						   print*,'	'
      END IF
   
end subroutine HDF5_load_solution


!**********************************************************************
! Save HDG matrix (CSR) in HDF5 file format
!**********************************************************************
subroutine HDF5_save_CSR_matrix(fname)
   USE globals
   implicit none
   
   character(LEN=*) :: fname
   character(70)  :: npr,nid
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF      
      call HDF5_create(fname_complete,file_id,ierr)
      call HDF5_integer_saving(file_id,MatK%n,'n')
      call HDF5_integer_saving(file_id,MatK%nnz,'nnz')
      call HDF5_array1D_saving_int(file_id,MatK%cols,MatK%nnz,'cols')
      call HDF5_array1D_saving_int(file_id,MatK%rowptr,MatK%n+1,'rowptr')
      call HDF5_array1D_saving_int(file_id,MatK%loc2glob,MatK%n,'loc2glob')
      call HDF5_array1D_saving(file_id,MatK%vals,MatK%nnz,'vals')
      call HDF5_close(file_id)
      ! Message to confirm succesful creation and filling of file
!      IF (MPIvar%glob_id.eq.0) THEN
!						   print*,'Output written to file ', trim(adjustl(fname_complete))
!						   print*,'	'
!      END IF
   
end subroutine HDF5_save_CSR_matrix



!**********************************************************************
! Save 3D array in HDF5 file format
!**********************************************************************
subroutine HDF5_save_array(Arr,fname)
   USE globals
   implicit none
   
   REAL, DIMENSION(:,:,:), INTENT(IN) :: Arr
   character(LEN=*) :: fname
   character(70)  :: npr,nid
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF      
      call HDF5_create(fname_complete,file_id,ierr)
      call HDF5_array3D_saving(file_id,Arr,size(Arr,1),size(Arr,2),size(Arr,3),'array')
      call HDF5_close(file_id)
      ! Message to confirm succesful creation and filling of file
!      IF (MPIvar%glob_id.eq.0) THEN
!						   print*,'Output written to file ', trim(adjustl(fname_complete))
!						   print*,'	'
!      END IF
   
end subroutine HDF5_save_array





!**********************************************************************
! Save 2D array in HDF5 file format
!**********************************************************************
subroutine HDF5_save_matrix(Mat,fname)
   USE globals
   implicit none
   
   REAL, DIMENSION(:,:), INTENT(IN) :: Mat
   character(LEN=*) :: fname
   character(70)  :: npr,nid
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF      
      call HDF5_create(fname_complete,file_id,ierr)
      call HDF5_array2D_saving(file_id,Mat,size(Mat,1),size(Mat,2),'mat')
      call HDF5_close(file_id)
      ! Message to confirm succesful creation and filling of file
!      IF (MPIvar%glob_id.eq.0) THEN
!						   print*,'Output written to file ', trim(adjustl(fname_complete))
!						   print*,'	'
!      END IF
   
end subroutine HDF5_save_matrix


!**********************************************************************
! Save 2D array in HDF5 file format
!**********************************************************************
subroutine HDF5_save_vector(Vec,fname)
   USE globals
   implicit none
   
   REAL, DIMENSION(:), INTENT(IN) :: Vec
   character(LEN=*) :: fname
   character(70)  :: npr,nid
   integer :: ierr
   character(len=1000) :: fname_complete
   integer(HID_T) :: file_id
   
   IF (MPIvar%glob_size.GT.1) THEN
      write(nid,*) MPIvar%glob_id+1
      write(npr,*) MPIvar%glob_size
      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'.h5'
   ELSE
      fname_complete = trim(adjustl(fname))//'.h5'
   END IF      
      call HDF5_create(fname_complete,file_id,ierr)
      call HDF5_array1D_saving(file_id,Vec,size(Vec,1),'vec')
      call HDF5_close(file_id)
      ! Message to confirm succesful creation and filling of file
!      IF (MPIvar%glob_id.eq.0) THEN
!						   print*,'Output written to file ', trim(adjustl(fname_complete))
!						   print*,'	'
!      END IF
   
end subroutine HDF5_save_vector


END MODULE in_out
