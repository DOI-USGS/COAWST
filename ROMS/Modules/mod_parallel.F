#include "cppdefs.h"
      MODULE mod_parallel
!
!svn $Id: mod_parallel.F 1607 2008-07-27 03:03:34Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains all variables used for parallelization         !
!                                                                      !
!=======================================================================
!
        USE mod_param
        USE mod_strings, ONLY: Nregion
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
!
!  Number of shared-memory parallel threads.  In distributed memory
!  configurations, its value must be equal to one.
!
        integer :: numthreads = 1
!
!  Number distributed memory nodes.
!
        integer :: numnodes = 0

#ifdef AIR_OCEAN
!
!  Parallel nodes assined to the atmosphere model.
!
        integer :: peATM_frst          ! first atmosphere parallel node
        integer :: peATM_last          ! last  atmosphere parallel node 
#endif
#ifdef WAVES_OCEAN
!
!  Parallel nodes assined to the wave model.
!
        integer :: peWAV_frst          ! first atmosphere parallel node
        integer :: peWAV_last          ! last  atmosphere parallel node 
#endif
!
!  Parallel nodes assined to the ocean model.
!
        integer :: peOCN_frst          ! first ocean parallel node
        integer :: peOCN_last          ! last  ocean parallel node
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
!
        integer  :: proc(0:1,4,Ngrids)

        real(r8) :: Cstr(0:Nregion,4,Ngrids)
        real(r8) :: Cend(0:Nregion,4,Ngrids)
        real(r8) :: Csum(0:Nregion,4,Ngrids)

#ifdef _OPENMP
!$OMP THREADPRIVATE (proc)
!$OMP THREADPRIVATE (Cstr, Cend)
#endif
!
!  Distributed-memory master process and rank of the local process.
!
        integer, parameter :: MyMaster = 0
        integer :: MyRank = 0

#ifdef DISTRIBUTE
# ifdef MPI
!
!  Ocean model MPI group communicator handle.
!
        integer :: OCN_COMM_WORLD
# endif
!
!  Type of message-passage floating point bindings.
!
# ifdef DOUBLE_PRECISION
#  ifdef MPI
        integer, parameter :: MP_FLOAT = MPI_DOUBLE_PRECISION
!!      integer, parameter :: MP_FLOAT = MPI_REAL8
#  endif
# else
#  ifdef MPI
        integer, parameter :: MP_FLOAT = MPI_REAL
!!      integer, parameter :: MP_FLOAT = MPI_REAL4
#  endif
# endif
#endif

        CONTAINS

        SUBROUTINE initialize_parallel
!
!=======================================================================
!                                                                      !
!  This routine initializes and spawn distribute-memory nodes.         !
!                                                                      !
!=======================================================================
!
          USE mod_param
          USE mod_iounits
          USE mod_scalars
          USE mod_strings, ONLY: Nregion
!
!  Local variable declarations.
!
          integer :: i
#ifdef DISTRIBUTE
          integer :: MyError
#endif
#ifndef DISTRIBUTE
          integer :: my_numthreads
!
!-----------------------------------------------------------------------
!  Initialize shared-memory (OpenMP) or serial configuration.
!-----------------------------------------------------------------------
!
!  Inquire number of threads in parallel region.
!
          numthreads=my_numthreads()
          Master=.TRUE.
          InpThread=.TRUE.
          OutThread=.TRUE.
#endif
#ifdef DISTRIBUTE
# ifdef MPI
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration.
!-----------------------------------------------------------------------
!
!  Get the number of processes in the group associated with the world
!  communicator.
!
          numthreads=1
          CALL mpi_comm_size (OCN_COMM_WORLD, numnodes, MyError)
          IF (MyError.ne.0) THEN
            WRITE (stdout,10)
  10        FORMAT (/,' ROMS/TOMS - Unable to inquire number of',       &
     &              ' processors in the group.')
            exit_flag=6
            RETURN
          END IF
!
!  Identify master thread and output thread.
!
          Master=.FALSE.
          InpThread=.FALSE.
          OutThread=.FALSE.
          IF (MyRank.eq.MyMaster) THEN
            Master=.TRUE.
            InpThread=.TRUE.
            OutThread=.TRUE.
          END IF
# endif
#endif
!
!-----------------------------------------------------------------------
!  Initialize profiling variables.
!-----------------------------------------------------------------------
!
          proc(0:1,1:4,1:Ngrids)=0
          Cstr(0:Nregion,1:4,1:Ngrids)=0.0_r8
          Cend(0:Nregion,1:4,1:Ngrids)=0.0_r8
          Csum(0:Nregion,1:4,1:Ngrids)=0.0_r8

          RETURN
        END SUBROUTINE initialize_parallel
      END MODULE mod_parallel
