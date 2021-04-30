#include "cppdefs.h"
      MODULE mod_parallel
!
!svn $Id: mod_parallel.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none

#ifdef MPI
      include 'mpif.h'
#endif
!
!  Switch to identify master processor. In serial and shared-memory
!  applications it is always true.
!
      logical :: Master
!
!  Switch to identify which thread is processing input/output files.
!  In distributed-memory applications, this thread can be the master
!  thread or all threads in case of parallel output. In serial and
!  shared-memory applications it is always true.
!
      logical :: InpThread
      logical :: OutThread

!$OMP THREADPRIVATE (InpThread, OutThread)
!
!  Number of shared-memory parallel threads or distributed-memory
!  parallel nodes.
!
      integer :: numthreads
#ifndef DISJOINTED
      integer :: ForkSize
#endif
!
!  First and last tile to process in a tiled application.
!
      integer, allocatable :: first_tile(:)
      integer, allocatable :: last_tile(:)

!$OMP THREADPRIVATE (first_tile, last_tile)

#ifdef AIR_OCEAN
!
!  Parallel nodes assined to the atmosphere model.
!
      integer :: peATM_frst            ! first atmosphere parallel node
      integer :: peATM_last            ! last  atmosphere parallel node
#endif
#ifdef WAVES_OCEAN
!
!  Parallel nodes assined to the wave model.
!
      integer :: peWAV_frst            ! first atmosphere parallel node
      integer :: peWAV_last            ! last  atmosphere parallel node
#endif
!
!  Parallel nodes assined to the ocean model.
!
      integer :: peOCN_frst            ! first ocean parallel node
      integer :: peOCN_last            ! last  ocean parallel node
!
!  Parallel threads/nodes counters used in critical parallel regions.
!
      integer :: tile_count = 0
      integer :: block_count = 0
      integer :: thread_count = 0
!
!  Profiling variables as function of parallel thread:
!
!    proc          Parallel process ID.
!    Cstr          Starting time for program region.
!    Cend          Ending time for program region.
!    Csum          Accumulated time for progam region.
!    Ctotal        Total time profiled
!    total_cpu     Total elapsed CPU
!
      integer, allocatable :: proc(:,:,:)

      real(r8) :: Ctotal, total_cpu, total_model(4)

      real(r8), allocatable :: Cstr(:,:,:)
      real(r8), allocatable :: Cend(:,:,:)
      real(r8), allocatable :: Csum(:,:,:)

!$OMP THREADPRIVATE (proc)
!$OMP THREADPRIVATE (Cstr, Cend)

#if defined DISTRIBUTE && defined PROFILE
!
!  Switch manage time clock in "mp_bcasts". During initialization is
!  set to .FALSE. because the profiling variables cannot be allocated
!  and initialized before the "Ngrids" parameter is known.
!
      logical :: Lwclock = .FALSE.
#endif
!
!  Distributed-memory master process.
!
      integer :: MyMaster = 0

#if defined DISTRIBUTE && defined DISJOINTED
!
!  Split communicator parameters.
!
      logical :: FullMaster     ! full communicator master process

      integer :: ForkColor = -1 ! Fork split criteria: FullRank/ForkSize
      integer :: ForkKey        ! Fork split rank ordering: FullRank
      integer :: ForkSize       ! number of processes in FORK_COMM_WORLD
      integer :: FullSize       ! number of precesses in FULL_COMM_WORLD
      integer :: NsubGroups     ! number of disjointed subgroups
      integer :: TaskColor = -1 ! Task split criteria: FullRank/TaskSize
      integer :: TaskKey        ! Task split rank ordering: FullRank
      integer :: TaskSize       ! number of processes in TASK_COMM_WORLD
#endif
!
!  Rank of the parallel local process.
!
      integer :: FullRank = 0
      integer :: TaskRank = 0
      integer :: MyRank = 0
      integer :: MyThread = 0

!$OMP THREADPRIVATE (MyThread)

#ifdef DISTRIBUTE
# ifdef MPI
!
!  Distributed-memory group communicator handles.
!
#  ifdef DISJOINTED
      integer :: FULL_COMM_WORLD            ! full communicator
      integer :: FORK_COMM_WORLD            ! fork communicator
      integer :: TASK_COMM_WORLD            ! task communicator
#  endif
      integer :: OCN_COMM_WORLD             ! internal ROMS communicator
!
!  Set mpi_info opaque object handle.
!
      integer :: MP_INFO = MPI_INFO_NULL
# endif
!
!  Type of message-passage floating point bindings.
!
# ifdef DOUBLE_PRECISION
#  ifdef MPI
      integer, parameter :: MP_FLOAT = MPI_DOUBLE_PRECISION
!!    integer, parameter :: MP_FLOAT = MPI_REAL8
#  endif
# else
#  ifdef MPI
      integer, parameter :: MP_FLOAT = MPI_REAL
!!    integer, parameter :: MP_FLOAT = MPI_REAL4
#  endif
# endif
# ifdef MPI
      integer, parameter :: MP_DOUBLE = MPI_DOUBLE_PRECISION
# endif
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_parallel
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_strings, ONLY: Nregion
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
!  First and last tile to process in a tiled application.
!
      IF (.not.allocated(first_tile)) THEN
        allocate ( first_tile(Ngrids) )
      END IF
      IF (.not.allocated(last_tile)) THEN
        allocate ( last_tile(Ngrids) )
      END IF
!
!  Time profiling variables.
!
      IF (.not.allocated(proc)) THEN
        allocate ( proc(0:1,4,Ngrids) )
        proc(0:1,1:4,1:Ngrids)=0
      END IF

      IF (.not.allocated(Cstr)) THEN
        allocate ( Cstr(0:Nregion,4,Ngrids) )
        Cstr(0:Nregion,1:4,1:Ngrids)=0.0_r8
      END IF

      IF (.not.allocated(Cend)) THEN
        allocate ( Cend(0:Nregion,4,Ngrids) )
        Cend(0:Nregion,1:4,1:Ngrids)=0.0_r8
      END IF
!$OMP END PARALLEL

      IF (.not.allocated(Csum)) THEN
        allocate ( Csum(0:Nregion,4,Ngrids) )
        Csum(0:Nregion,1:4,1:Ngrids)=0.0_r8
      END IF
!
! Initialize other profiling variables.
!
      Ctotal=0.0_r8
      total_cpu=0.0_r8
      total_model=0.0_r8

#if defined DISTRIBUTE && defined PROFILE
!
! Activate wall clock switch used only in "mp_bcasts". This switch
! is set to .FALSE. during initialization before calling "inp_par.F"
! because the above profiling variables are allocated and initialized
! after the value of "Ngrids" is known.
!
      Lwclock=.TRUE.
#endif

      RETURN
      END SUBROUTINE allocate_parallel
!
      SUBROUTINE initialize_parallel
!
!=======================================================================
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: i
#ifdef DISTRIBUTE
      integer :: MyError
#else
      integer :: my_numthreads, my_threadnum
#endif

#if defined _OPENMP
!
!-----------------------------------------------------------------------
!  Initialize shared-memory (OpenMP) configuration.
!-----------------------------------------------------------------------
!
!  Disable dynamic adjustment, by the run-time environment, of the
!  number of threads available to execute parallel regions.
!
      CALL omp_set_dynamic (.FALSE.)
!
!  Inquire number of threads in parallel region.  Set master and I/O
!  switches.
!
      numthreads=my_numthreads()
      IF (my_threadnum().eq.0) THEN
        InpThread=.TRUE.
        OutThread=.TRUE.
      ELSE
        InpThread=.FALSE.
        OutThread=.FALSE.
      END IF
      Master=.TRUE.

#elif defined DISTRIBUTE
# ifdef MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration.
!-----------------------------------------------------------------------
!
!  Get the number of processes in the group associated with the world
!  communicator.  Here FullRank and MyRank are the same. It is computed
!  for consistency wity concurrent applications.
!
      CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' ROMS/TOMS - Unable to inquire number of',           &
     &            ' processors in the group.')
        exit_flag=6
        RETURN
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, FullRank, MyError)
      ForkSize=numthreads
!
!  Identify master, input and output threads.
!
#  ifdef PARALLEL_IO
      Master=.FALSE.
      InpThread=.TRUE.
      OutThread=.TRUE.
      IF (MyRank.eq.MyMaster) THEN
        Master=.TRUE.
      END IF
#  else
      Master=.FALSE.
      InpThread=.FALSE.
      OutThread=.FALSE.
      IF (MyRank.eq.MyMaster) THEN
        Master=.TRUE.
        InpThread=.TRUE.
        OutThread=.TRUE.
      END IF
#  endif
# endif
#else
!
!-----------------------------------------------------------------------
!  Initialize serial configuration.
!-----------------------------------------------------------------------
!
      numthreads=my_numthreads()
      Master=.TRUE.
      InpThread=.TRUE.
      OutThread=.TRUE.
#endif
!
      RETURN
      END SUBROUTINE initialize_parallel

#if defined DISTRIBUTE && defined DISJOINTED
!
      SUBROUTINE split_communicator (Nsplit, Ntasks)
!
!=======================================================================
!                                                                      !
!  This routine splits the main distributed-memory communicator object !
!  (OCN_COMM_WORLD saved as FULL_COMM_WORLD) into several disjointed   !
!  subgroups processes with new communicator handle (FORK_COMM_WORLD). !
!  Notice that at entry, the OCN_COMM_WORLD handle is the same object  !
!  as MPI_COMM_WORLD, if no multi-model coupling.                      !
!                                                                      !
!  If Ntask=2, a second split is done. Both tasks (1 and 2) have equal !
!  number processors. Task 1 uses the first half, whereas Task 2 uses  !
!  the second half of all the processes assigned to ROMS. The extra    !
!  communicator handle (TASK_COMM_WORLD) can be used for collective    !
!  operations within each task section.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nsplit     Number of disjointed by related subgroups (integer)   !
!     Ntasks     Number of different code sections tasks (integer)     !
!                  Usually, Ntasks=1 or Ntasks=2                       !
!                                                                      !
!  Therefore,                                                          !
!                                                                      !
!     ForkSize        Number of processes assigned for each subgroup:  !
!                       ForkSize = FullSize/(Nplit*Ntasks)             !
!                                                                      !
!     FullSize        Total number of processes assigned to ROMS       !
!                                                                      !
!     TaskSize        Number of processes assigned for each task       !
!                       TaskSize=FullSize/Ntasks                       !
!                                                                      !
!     NsubGroups      Total number of disjointed subgroups             !
!                       NsubGroups=Nsplit*Ntasks                       !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations
!
      integer, intent(in) :: Nsplit, Ntasks
!
!  Local variable declarations.
!
      integer :: MyError, MySize
      integer :: i, Lstr, Serror
!
      integer, allocatable :: GroupRanks(:)
!
      character (len=MPI_MAX_ERROR_STRING) :: string
!
!-----------------------------------------------------------------------
!  Split ROMS communicator into disjointed subgroups, FORK_COMM_WORLD.
!  It is used, for example, for the concurrent time-stepping of the
!  ROMS kernel.
!-----------------------------------------------------------------------
!
!  Save full communicator (starting handle) and inquire about its size
!  and rank.
!
      FULL_COMM_WORLD=OCN_COMM_WORLD
      CALL mpi_comm_rank (FULL_COMM_WORLD, FullRank, MyError)
      CALL mpi_comm_size (FULL_COMM_WORLD, FullSize, MyError)
      FullMaster=FullRank.eq.MyMaster
!
!  Split the full communicator into sub-communicators based on color and
!  key. It is a collective operation.
!
!  For example, if FullSize=8 and Nsplit=2, we get ForkSize=4 and:
!
!  FullRank:    0    1    2    3    4    5    6    7    FULL_COMM_WORLD
!  ForkColor:   0    0    0    0    1    1    1    1    split criteria
!  ForkKey:     0    1    2    3    4    5    6    7    same as FullRank
!  ForkRank:    0    1    2    3    0    1    2    3    FORK_COMM_WORLD
!
!
!  The color determines in which of the sub-communicators the current
!  process (spawn by FULL_COMM_WORLD) will fall.
!
!  The key argument determines the rank ordering within the split
!  communicator. The process that passes in the smallest value for
!  key will be rank 0, the next smallest will be rank 1, and so on.
!  Here, the key is set to the rank of the full communicator since we
!  want all of the processes in the split communicator to be in the
!  same order that they were in the full communicator.
!
!  The resulting communicator, FORK_COMM_WORLD, has the same handle
!  value (context ID) in all sub-communicators. However, each process
!  only have reference to the sub-communicator that it belongs. It is
!  an opaque object.
!
      NsubGroups=Nsplit*Ntasks
      IF (MOD(FullSize, NsubGroups).ne.0) THEN
        IF (FullMaster) THEN
          WRITE (stdout,10) NsubGroups, FullSize
        END IF
        exit_flag=2
        RETURN
      END IF
      ForkSize=FullSize/NsubGroups               ! integer operation
      ForkColor=FullRank/ForkSize                ! integer operation
      ForkKey=FullRank
      CALL mpi_comm_split (FULL_COMM_WORLD, ForkColor, ForkKey,         &
     &                     FORK_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        WRITE (stdout,20) 'MPI_COMM_SPLIT', 'TASK_COMM_WORLD',          &
     &                     FullRank, MyError, TRIM(string)
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  If Ntasks=2, split the FULL_COMM_WORLD into disjointed tasks
!  communicators.
!-----------------------------------------------------------------------
!
      IF (Ntasks.gt.1) THEN
        TaskSize=FullSize/Ntasks
        TaskColor=FullRank/TaskSize
        TaskKey=FullRank
        CALL mpi_comm_split (FULL_COMM_WORLD, TaskColor, TaskKey,       &
     &                       TASK_COMM_WORLD, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          WRITE (stdout,20) 'MPI_COMM_SPLIT', 'TASK_COMM_WORLD',        &
     &                      FullRank, MyError, TRIM(string)
          exit_flag=2
          RETURN
        END IF
        CALL mpi_comm_rank (TASK_COMM_WORLD, TaskRank, MyError)
      ELSE
        TASK_COMM_WORLD=FULL_COMM_WORLD
        TaskSize=FullSize
        TaskRanK=FullRank
        TaskKey=0
        TaskColor=0
      END IF
!
 10   FORMAT (/,' SPLIT_COMMUNICATOR - illegal configuration, ',        &
     &        'Nsplit * Ntasks = ',i0,/,22x,                            &
     &        'needs to be less, equal, or a multiple of the ',         &
     &        'communicator size, FullSize = ',i0,/,22x,                &
     &        'That is, MOD(FullSise, Nsplit*Ntask)=0')
 20   FORMAT (/,' SPLIT_COMMUNICATOR - error during ',a,                &
     &        ' call, Comm = ',a,', Rank = ',i0,' Error = ',i0,/,22x,a)
!
      RETURN
      END SUBROUTINE split_communicator
!
      SUBROUTINE assign_communicator (choice)
!
!=======================================================================
!                                                                      !
!  This routine selects request distributed memory communicator and    !
!  it associated parameters and control switches.  It assigns the      !
!  requested communicator to the OCN_COMM_WORLD handle.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     choice       Requested communicator keyword (string)             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: choice
!
!  Local variable declarations.
!
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Assing internal ROMS communicator, OCN_COMM_WORLD.
!-----------------------------------------------------------------------
!
!  Set communicator handle and get its size and rank.
!
      SELECT CASE (TRIM(choice))
        CASE ('FORK', 'SPLIT', 'CONCURRENT')
          OCN_COMM_WORLD=FORK_COMM_WORLD
          CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
          CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
          Master=MyRank.eq.MyMaster
        CASE ('FULL')
          OCN_COMM_WORLD=FULL_COMM_WORLD
          CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
          CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
          Master=MyRank.eq.MyMaster
      END SELECT
      CALL mpi_barrier (OCN_COMM_WORLD, MyError)
!
!  Initialize parallel control switches.
!
# ifdef PARALLEL_IO
      InpThread=.TRUE.
      OutThread=.TRUE.
# else
      IF (MyRank.eq.0) THEN
        InpThread=.TRUE.
        OutThread=.TRUE.
      ELSE
        InpThread=.FALSE.
        OutThread=.FALSE.
      END IF
# endif
!
      RETURN
      END SUBROUTINE assign_communicator
#endif

      END MODULE mod_parallel
