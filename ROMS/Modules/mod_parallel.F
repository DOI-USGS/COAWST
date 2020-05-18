#include "cppdefs.h"
      MODULE mod_parallel
!
!svn $Id: mod_parallel.F 995 2020-01-10 04:01:28Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
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
!
!  Rank of the parallel local process.
!
      integer :: MyRank = 0
      integer :: MyThread = 0

!$OMP THREADPRIVATE (MyThread)

#ifdef DISTRIBUTE
# ifdef MPI
!
!  Ocean model MPI group communicator handle.
!
      integer :: OCN_COMM_WORLD
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
!  communicator.
!
      CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
      IF (MyError.ne.0) THEN
        WRITE (stdout,10)
  10    FORMAT (/,' ROMS/TOMS - Unable to inquire number of',           &
     &            ' processors in the group.')
        exit_flag=6
        RETURN
      END IF
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

      RETURN
      END SUBROUTINE initialize_parallel

      END MODULE mod_parallel
