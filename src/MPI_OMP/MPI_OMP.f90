!*************************************************
! project: TOKAM3X
! file: MPI_OMP.f90
! date: 02/07/2012
! Initialisation of MPI and OpenMP parallelization
!*************************************************




MODULE MPI_OMP
   USE mpi
   IMPLICIT NONE

   ! Type for MPI parallelization info
   TYPE MPI_typ
      integer :: glob_id, glob_size ! global communicator size and id of the local process in the communicator
#ifdef TOR3D      
      integer :: ntor, npol, itor, ipol
#endif      
      ! timing stuffs
      real(8) :: ctime1,ctime2,ctime3,ctime4,ctime5,ctime6
      real(8) :: rtime1,rtime2,rtime3,rtime4,rtime5,rtime6      

   END TYPE

   ! Type for OMP parallelization info
   TYPE OMP_typ
      integer :: Nthreads ! number of openMP threads
   END TYPE

   TYPE(MPI_typ) :: MPIvar
   TYPE(OMP_typ) :: OMPvar

   PUBLIC :: MPIvar, OMPvar

CONTAINS

!********************************************
! Initialization of MPI/OMP
!********************************************
subroutine init_MPI_OMP()

   integer :: IERR, OMP_GET_MAX_THREADS, MPI_THREAD_provided, MPI_THREAD_required

   ! Initialization of the MPI communicator
!  MPI_THREAD_required = MPI_THREAD_SINGLE
   MPI_THREAD_required = MPI_THREAD_MULTIPLE
!#ifdef THREAD_FUNNELED
!   MPI_THREAD_required = MPI_THREAD_FUNNELED
!#else
!   MPI_THREAD_required = MPI_THREAD_MULTIPLE
!#endif
   call MPI_init_thread(MPI_THREAD_required,MPI_THREAD_provided,IERR)
   call MPI_Comm_size(MPI_COMM_WORLD,MPIvar%glob_size,IERR)
#ifndef PARALL
   if (MPIvar%glob_size>1) then
      write(6,*) 'Error: using serial version in a parallel environment'
      stop
   endif
#endif   
   call MPI_Comm_rank(MPI_COMM_WORLD,MPIvar%glob_id,IERR)
   if ((MPI_THREAD_provided.NE.MPI_THREAD_required).AND.(MPIvar%glob_id.EQ.0)) then
      print*,'Problem with initialization of multi-threaded MPI.'
      print*,'Required: ', MPI_THREAD_required
      print*,'Provided: ', MPI_THREAD_provided
      print*,'(MPI_THREAD_SINGLE = ',MPI_THREAD_SINGLE,', MPI_THREAD_FUNNELED = ',MPI_THREAD_FUNNELED, &
		', MPI_THREAD_SERIALIZED = ',MPI_THREAD_SERIALIZED,', MPI_THREAD_MULTIPLE = ', MPI_THREAD_MULTIPLE,')'
      print*,'Exiting...'
      stop
   endif
   
   ! Makes sure that all the processes in the communicator have the same number of openMP threads
   if (MPIvar%glob_id.EQ.0) then
      OMPvar%Nthreads = OMP_GET_MAX_THREADS()
   endif
   call MPI_BCast(OMPvar%Nthreads,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   call OMP_SET_NUM_THREADS(OMPvar%Nthreads)
#ifdef TOR3D
   MPIvar%ntor =1
   MPIvar%npol =1
   MPIvar%itor =0
   MPIvar%ipol =0
#endif  
   MPIvar%ctime1 = 0.
   MPIvar%ctime2 = 0.
   MPIvar%ctime3 = 0.
   MPIvar%ctime4 = 0.
   MPIvar%ctime5 = 0.
   MPIvar%ctime6 = 0.
   MPIvar%rtime1 = 0.
   MPIvar%rtime2 = 0.
   MPIvar%rtime3 = 0.
   MPIvar%rtime4 = 0.
   MPIvar%rtime5 = 0.
   MPIvar%rtime6 = 0.

end subroutine init_MPI_OMP


END MODULE MPI_OMP

