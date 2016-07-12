#include "cppdefs.h"
#if (defined FOUR_DVAR    || defined VERIFICATION) && \
     defined OBSERVATIONS && \
     defined DISTRIBUTE   && defined SOLVE3D
      SUBROUTINE obs_depth (ng, tile, model)
!
!svn $Id: obs_depth.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine proccess the observation vertical position, Zobs. If   !
!  the input position is in meters (Zobs < 0) instead of  fractional   !
!  level (Zobs > 0; 1 <= Zobs <= N),  zero out its Zobs value in all   !
!  nodes in group but itself to faciliate exchange between all tiles   !
!  with "mp_collect", after model extraction.                          !
!                                                                      !
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      logical :: r_bound, u_bound, v_bound

      integer :: Mstr, Mend, iobs, itrc

      real(r8) :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  If distributed-memory and Zobs < 0, zero-out Zobs values in all
!  nodes but itself to facilitate exchages between tiles latter.
!-----------------------------------------------------------------------
!
!  Set starting index of obervation vectors for reading.  In weak
!  constraint, the entire observation data is loaded. Otherwise,
!  only the observation for the current time window are loaded
!  and started from vector index one.
!
# ifdef WEAK_CONSTRAINT
      Mstr=NstrObs(ng)
      Mend=NendObs(ng)
# else
      Mstr=1
      Mend=Nobs(ng)
# endif
!
!  Process vertical positions.
!
      DO iobs=Mstr,Mend
        IF (Zobs(iobs).lt.0.0_r8) THEN
          r_bound=((rXmin(ng) .le.Xobs(iobs)).and.                      &
     &             (Xobs(iobs).lt.rXmax(ng))).and.                      &
     &            ((rYmin(ng) .le.Yobs(iobs)).and.                      &
     &             (Yobs(iobs).lt.rYmax(ng)))
          u_bound=((uXmin(ng) .le.Xobs(iobs)).and.                      &
     &             (Xobs(iobs).lt.uXmax(ng))).and.                      &
     &            ((uYmin(ng) .le.Yobs(iobs)).and.                      &
     &             (Yobs(iobs).lt.uYmax(ng)))
          v_bound=((vXmin(ng) .le.Xobs(iobs)).and.                      &
     &             (Xobs(iobs).lt.vXmax(ng))).and.                      &
     &            ((vYmin(ng) .le.Yobs(iobs)).and.                      &
     &             (Yobs(iobs).lt.vYmax(ng)))
          IF (ObsType(iobs).eq.isUvel) THEN
            IF (.not.u_bound) THEN
              Zobs(iobs)=IniVal
            END IF
          ELSE IF (ObsType(iobs).eq.isVvel) THEN
            IF (.not.v_bound) THEN
              Zobs(iobs)=IniVal
            END IF
          ELSE
            DO itrc=1,NT(ng)
              IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                IF (.not.r_bound) THEN
                  Zobs(iobs)=IniVal
                END IF
              END IF
            END DO
          END IF
        END IF
      END DO

      RETURN
      END SUBROUTINE obs_depth
#else
      SUBROUTINE obs_depth
      RETURN
      END SUBROUTINE obs_depth
#endif
