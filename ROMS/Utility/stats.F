#include "cppdefs.h"
      MODULE stats_mod
!
!git $Id$
!svn $Id: stats.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to compute the statistics of  !
!  provided field array, F.  The following information is computed:    !
!                                                                      !
!    S % count                   processed values count                !
!    S % min                     minimum value                         !
!    S % max                     maximum value                         !
!    S % avg                     arithmetic mean                       !
!    S % rms                     root mean square                      !
!                                                                      !
!  Notice that to compute high order moments (standard deviation,      !
!  variance, skewness, and kurosis) cannot be computed in parallel     !
!  with a single call because we need to know the mean first.          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    stats_2dfld      Statistic information for 2D state field         !
!    stats_3dfld      Statistic information for 3D state field         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PUBLIC :: stats_2dfld
      PUBLIC :: stats_3dfld
!
      CONTAINS
!
      SUBROUTINE stats_2dfld (ng, tile, model, gtype, S,                &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        F, Fmask, debug)
!
!=======================================================================
!                                                                      !
!  This routine computes requested 2D-field statistics.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     tile         Domain partition (integer)                          !
!     model        Calling model identifier (integer)                  !
!     gtype        Grid type (integer)                                 !
!     LBi          I-dimension Lower bound (integer)                   !
!     UBi          I-dimension Upper bound (integer)                   !
!     LBj          J-dimension Lower bound (integer)                   !
!     UBj          J-dimension Upper bound (integer)                   !
!     F            2D state field (2D tiled array)                     !
!     Fmask        Land/sea mask (OPTIONAL, 2D tiled array)            !
!     debug        Switch to debug (OPTIONAL, logical)                 !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     S            Field statistics (structure)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: debug

      integer, intent(in) :: ng, tile, model, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: F(LBi:,LBj:)
      real(r8), intent(in), optional :: Fmask(LBi:,LBj:)
#else
      real(r8), intent(in) :: F(LBi:UBi,LBj:UBj)
      real(r8), intent(in), optional :: Fmask(LBi:UBi,LBj:UBj)
#endif
      TYPE(T_STATS), intent(inout) :: S
!
!  Local variable declarations.
!
      logical :: Lprint

      integer :: Imin, Imax, Jmin, Jmax, NSUB
      integer :: i, j

      integer :: my_threadnum

      real(r8) :: fac
      real(r8) :: my_count, my_max, my_min
      real(r8) :: my_avg, my_rms

#ifdef DISTRIBUTE
      real(r8), dimension(5) :: rbuffer
      character (len=3), dimension(5) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  Set tile indices.
!-----------------------------------------------------------------------
!
!  Determine I- and J-indices according to staggered C-grid type.
!
      SELECT CASE (ABS(gtype))
        CASE (p2dvar)
          Imin=BOUNDS(ng)%IstrP(tile)
          Imax=BOUNDS(ng)%IendP(tile)
          Jmin=BOUNDS(ng)%JstrP(tile)
          Jmax=BOUNDS(ng)%JendP(tile)
        CASE (r2dvar)
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE (u2dvar)
          Imin=BOUNDS(ng)%IstrP(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE (v2dvar)
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrP(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE DEFAULT
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
      END SELECT
!
!-----------------------------------------------------------------------
!  Compute field statistics.
!-----------------------------------------------------------------------
!
!  Initialize local values.
!
      my_count=0.0_r8
      my_min= 1.0E+37_r8
      my_max=-1.0E+37_r8
      my_avg=0.0_r8
      my_rms=0.0_r8
!
      IF (PRESENT(debug)) THEN
        Lprint=debug
      ELSE
        Lprint=.FALSE.
      END IF
!
!  Compute field mean, range, and processed values count for each tile.
!
      IF (PRESENT(Fmask)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            IF (Fmask(i,j).gt.0.0_r8) THEN
              my_count=my_count+1.0_r8
              my_min=MIN(my_min, F(i,j))
              my_max=MAX(my_max, F(i,j))
              my_avg=my_avg+F(i,j)
              my_rms=my_rms+F(i,j)*F(i,j)
            END IF
          END DO
        END DO
      ELSE
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            my_count=my_count+1.0_r8
            my_min=MIN(my_min, F(i,j))
            my_max=MAX(my_max, F(i,j))
            my_avg=my_avg+F(i,j)
            my_rms=my_rms+F(i,j)*F(i,j)
          END DO
        END DO
      END IF
!
!  Compute global field mean, range, and processed values count.
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      NSUB=NtileX(ng)*NtileE(ng)         ! tiled application
#endif
!$OMP CRITICAL (F_STATS)
      S%count=S%count+my_count
      S%min=MIN(S%min,my_min)
      S%max=MAX(S%max,my_max)
      S%avg=S%avg+my_avg
      S%rms=S%rms+my_rms
      IF (Lprint) THEN
        WRITE (stdout,10)                                               &
     &            'thread = ', my_threadnum(),                          &
     &            'tile = ', tile,                                      &
     &            'tile_count = ', tile_count,                          &
     &            'my_count = ', my_count,                              &
     &            'S%count  = ', S%count,                               &
     &            'MINVAL   = ', MINVAL(F(Imin:Imax,Jmin:Jmax)),        &
     &            'MAXVAL   = ', MAXVAL(F(Imin:Imax,Jmin:Jmax)),        &
     &            'SUM      = ', SUM(F(Imin:Imax,Jmin:Jmax)),           &
     &            'MEAN     = ', SUM(F(Imin:Imax,Jmin:Jmax))/my_count,  &
     &            'my_min   = ', my_min,                                &
     &            'my_max   = ', my_max,                                &
     &            'S%min    = ', S%min,                                 &
     &            'S%max    = ', S%max,                                 &
     &            'my_avg   = ', my_avg,                                &
     &            'S%avg    = ', S%avg,                                 &
     &            'my_rms   = ', my_rms,                                &
     &            'S%rms    = ', S%rms
  10    FORMAT (10x,3(5x,a,i4),/,                                       &
     &          15x,a,f15.0,5x,a,f15.0,/,                               &
     &          6(15x,a,1p,e15.8,0p,5x,a,1p,e15.8,0p,/))
        CALL my_flush (stdout)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        rbuffer(1)=S%count
        rbuffer(2)=S%min
        rbuffer(3)=S%max
        rbuffer(4)=S%avg
        rbuffer(5)=S%rms
        op_handle(1)='SUM'
        op_handle(2)='MIN'
        op_handle(3)='MAX'
        op_handle(4)='SUM'
        op_handle(5)='SUM'
        CALL mp_reduce (ng, model, 5, rbuffer, op_handle)
        S%count=rbuffer(1)
        S%min=rbuffer(2)
        S%max=rbuffer(3)
        S%avg=rbuffer(4)
        S%rms=rbuffer(5)
#endif
!
!  Finalize computation of the mean and root mean square.
!
        IF (S%count.gt.0) THEN
          fac=1.0_r8/S%count
          S%avg=S%avg*fac
          S%rms=SQRT(S%rms*fac)
        ELSE
          S%avg=1.0E+37_r8
          S%rms=1.0E+37_r8
        END IF
      END IF
!$OMP END CRITICAL (F_STATS)
!$OMP BARRIER

      RETURN
      END SUBROUTINE stats_2dfld
!
      SUBROUTINE stats_3dfld (ng, tile, model, gtype, S,                &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        F, Fmask, debug)
!
!=======================================================================
!                                                                      !
!  This routine computes requested 3D-field statistics.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     tile         Domain partition (integer)                          !
!     model        Calling model identifier (integer)                  !
!     gtype        Grid type (integer)                                 !
!     LBi          I-dimension Lower bound (integer)                   !
!     UBi          I-dimension Upper bound (integer)                   !
!     LBj          J-dimension Lower bound (integer)                   !
!     UBj          J-dimension Upper bound (integer)                   !
!     LBk          K-dimension Lower bound (integer)                   !
!     UBk          K-dimension Upper bound (integer)                   !
!     F            3D state field (3D tiled array)                     !
!     Fmask        Land/sea mask (OPTIONAL, 2D tiled array)            !
!     debug        Switch to debug (OPTIONAL, logical)                 !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     S            Field statistics (structure)                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#endif
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: debug

      integer, intent(in) :: ng, tile, model, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: F(LBi:,LBj:,LBk:)
      real(r8), intent(in), optional :: Fmask(LBi:,LBj:)
#else
      real(r8), intent(in) :: F(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(in), optional :: Fmask(LBi:UBi,LBj:UBj)
#endif
      TYPE(T_STATS), intent(inout) :: S
!
!  Local variable declarations.
!
      logical :: Lprint

      integer :: Imin, Imax, Jmin, Jmax, NSUB
      integer :: i, j, k

      integer :: my_threadnum

      real(r8) :: fac
      real(r8) :: my_count, my_max, my_min
      real(r8) :: my_avg, my_rms

#ifdef DISTRIBUTE
      real(r8), dimension(5) :: rbuffer
      character (len=3), dimension(5) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  Set tile indices.
!-----------------------------------------------------------------------
!
!  Determine I- and J-indices according to staggered C-grid type.
!
      SELECT CASE (ABS(gtype))
        CASE (p3dvar)
          Imin=BOUNDS(ng)%IstrP(tile)
          Imax=BOUNDS(ng)%IendP(tile)
          Jmin=BOUNDS(ng)%JstrP(tile)
          Jmax=BOUNDS(ng)%JendP(tile)
        CASE (r3dvar)
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE (u3dvar)
          Imin=BOUNDS(ng)%IstrP(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE (v3dvar)
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrP(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
        CASE DEFAULT
          Imin=BOUNDS(ng)%IstrT(tile)
          Imax=BOUNDS(ng)%IendT(tile)
          Jmin=BOUNDS(ng)%JstrT(tile)
          Jmax=BOUNDS(ng)%JendT(tile)
      END SELECT
!
!-----------------------------------------------------------------------
!  Compute field statistics.
!-----------------------------------------------------------------------
!
!  Initialize local values.
!
      my_count=0.0_r8
      my_min= 1.0E+37_r8
      my_max=-1.0E+37_r8
      my_avg=0.0_r8
      my_rms=0.0_r8
!
      IF (PRESENT(debug)) THEN
        Lprint=debug
      ELSE
        Lprint=.FALSE.
      END IF
!
!  Compute field mean, range, and processed values count for each tile.
!
      IF (PRESENT(Fmask)) THEN
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              IF (Fmask(i,j).gt.0.0_r8) THEN
                my_count=my_count+1.0_r8
                my_min=MIN(my_min, F(i,j,k))
                my_max=MAX(my_max, F(i,j,k))
                my_avg=my_avg+F(i,j,k)
                my_rms=my_rms+F(i,j,k)*F(i,j,k)
              END IF
            END DO
          END DO
        END DO
      ELSE
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              my_count=my_count+1.0_r8
              my_min=MIN(my_min, F(i,j,k))
              my_max=MAX(my_max, F(i,j,k))
              my_avg=my_avg+F(i,j,k)
              my_rms=my_rms+F(i,j,k)*F(i,j,k)
            END DO
          END DO
        END DO
      END IF
!
!  Compute global field mean, range, and processed values count.
!  Notice that the critical region F_STATS is the same as in
!  "stats_2dfld" but we are usinf different counter "thread_count"
!  to avoid interference and race conditions in shared-memory.
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      NSUB=NtileX(ng)*NtileE(ng)         ! tiled application
#endif
!$OMP CRITICAL (F_STATS)
      S%count=S%count+my_count
      S%min=MIN(S%min,my_min)
      S%max=MAX(S%max,my_max)
      S%avg=S%avg+my_avg
      S%rms=S%rms+my_rms
      IF (Lprint) THEN
        WRITE (stdout,10)                                               &
     &            'thread = ', my_threadnum(),                          &
     &            'tile = ', tile,                                      &
     &            'tile_count = ', thread_count,                        &
     &            'my_count = ', my_count,                              &
     &            'S%count  = ', S%count,                               &
     &            'MINVAL   = ', MINVAL(F(Imin:Imax,Jmin:Jmax,:)),      &
     &            'MAXVAL   = ', MAXVAL(F(Imin:Imax,Jmin:Jmax,:)),      &
     &            'SUM      = ', SUM(F(Imin:Imax,Jmin:Jmax,:)),         &
     &            'MEAN     = ', SUM(F(Imin:Imax,Jmin:Jmax,:))/my_count,&
     &            'my_min   = ', my_min,                                &
     &            'my_max   = ', my_max,                                &
     &            'S%min    = ', S%min,                                 &
     &            'S%max    = ', S%max,                                 &
     &            'my_avg   = ', my_avg,                                &
     &            'S%avg    = ', S%avg,                                 &
     &            'my_rms   = ', my_rms,                                &
     &            'S%rms    = ', S%rms
  10    FORMAT (10x,3(5x,a,i4),/,                                       &
     &          15x,a,f15.0,5x,a,f15.0,/,                               &
     &          6(15x,a,1p,e15.8,0p,5x,a,1p,e15.8,0p,/))
        CALL my_flush (stdout)
      END IF
      thread_count=thread_count+1
      IF (thread_count.eq.NSUB) THEN
        thread_count=0
#ifdef DISTRIBUTE
        rbuffer(1)=S%count
        rbuffer(2)=S%min
        rbuffer(3)=S%max
        rbuffer(4)=S%avg
        rbuffer(5)=S%rms
        op_handle(1)='SUM'
        op_handle(2)='MIN'
        op_handle(3)='MAX'
        op_handle(4)='SUM'
        op_handle(5)='SUM'
        CALL mp_reduce (ng, model, 5, rbuffer, op_handle)
        S%count=rbuffer(1)
        S%min=rbuffer(2)
        S%max=rbuffer(3)
        S%avg=rbuffer(4)
        S%rms=rbuffer(5)
#endif
!
!  Finalize computation of the mean and root mean square.
!
        IF (S%count.gt.0) THEN
          fac=1.0_r8/S%count
          S%avg=S%avg*fac
          S%rms=SQRT(S%rms*fac)
        ELSE
          S%avg=1.0E+37_r8
          S%rms=1.0E+37_r8
        END IF
      END IF
!$OMP END CRITICAL (F_STATS)
!$OMP BARRIER

      RETURN
      END SUBROUTINE stats_3dfld

      END MODULE stats_mod
