#include "cppdefs.h"
!
!svn $Id: mp_routines.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package contains multi-processing routines used during         !
!  parallel applications:                                              !
!                                                                      !
!     my_flush         Flushes the contents of a unit buffer.          !
!     my_getarg        Returns the argument from command-line.         !
!     my_getpid        Returns process ID of the calling process.      !
!     my_numthreads    Returns number of threads that would            !
!                        execute in parallel regions.                  !
!     my_threadnum     Returns which thread number is working          !
!                        in a parallel region.                         !
!     my_wtime         Returns an elapsed wall time in seconds since   !
!                        an arbitrary time in the past.                !
!                                                                      !
!=======================================================================
!
!
!-----------------------------------------------------------------------
      SUBROUTINE my_flush (unit)
!-----------------------------------------------------------------------
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: unit
!
#if defined CYGWIN
#elif defined AIX
      CALL flush_ (unit)
#else
      CALL flush (unit)
#endif

      RETURN
      END SUBROUTINE my_flush

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
      SUBROUTINE my_getarg (Iarg, Carg)
!-----------------------------------------------------------------------
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Iarg
      character (len=*), intent(inout) :: Carg
!
!  Local variable declarations.
!
      integer :: Lstr, ierror

# ifdef CRAY
      CALL pxfgetarg (Iarg, Carg, Lstr, ierror)
# else
      CALL getarg (Iarg, Carg)
# endif

      RETURN
      END SUBROUTINE my_getarg
#endif
!
!-----------------------------------------------------------------------
      FUNCTION my_getpid ()
!-----------------------------------------------------------------------
!
      USE mod_kinds

      implicit none

      integer :: getpid, my_getpid

#if defined CRAYX1
      my_getpid=0
#elif defined CYGWIN
      my_getpid=0
#else
      my_getpid=getpid()
#endif

      RETURN
      END FUNCTION my_getpid
!
!-----------------------------------------------------------------------
      FUNCTION my_numthreads ()
!-----------------------------------------------------------------------
!
      USE mod_kinds
!
      implicit none
!
      integer :: my_numthreads

#ifdef DISTRIBUTE
      my_numthreads=1
#elif defined _OPENMP
      integer :: omp_get_max_threads
!!    integer :: omp_get_num_threads

      my_numthreads=omp_get_max_threads()
!!    my_numthreads=omp_get_num_threads()
#elif defined SGI
      integer :: mp_numthreads

      my_numthreads=mp_numthreads()
#else
      my_numthreads=1
#endif

      RETURN
      END FUNCTION my_numthreads
!
!-----------------------------------------------------------------------
      FUNCTION my_threadnum ()
!-----------------------------------------------------------------------
!
      USE mod_kinds
!
      implicit none
!
      integer :: my_threadnum

#if defined _OPENMP
      integer :: omp_get_thread_num

      my_threadnum=omp_get_thread_num()
#elif defined SGI
      integer :: mp_my_threadnum

      my_threadnum=mp_my_threadnum()
#else
      my_threadnum=0
#endif

      RETURN
      END FUNCTION my_threadnum
!
!-----------------------------------------------------------------------
      FUNCTION my_wtime (wtime)
!-----------------------------------------------------------------------
!
      USE mod_kinds
!
      implicit none
!
      real(r8) :: wtime(2)
      real(r8) :: my_wtime
#ifdef _OPENMP
      real(r8) :: omp_get_wtime
#endif
#ifdef CRAYX1
      integer :: count, count_rate, count_max
#endif

#if defined _OPENMP
      my_wtime=omp_get_wtime()
#elif defined CRAYX1
      CALL system_clock (count, count_rate, count_max)
      my_wtime=REAL(count,r8)/REAL(count_rate,r8)
#else
      CALL cpu_time (wtime(1))
      my_wtime=wtime(1)
#endif

      RETURN
      END FUNCTION my_wtime
