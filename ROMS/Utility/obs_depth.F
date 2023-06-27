#include "cppdefs.h"
#if (defined FOUR_DVAR    || defined VERIFICATION) && \
     defined OBSERVATIONS && \
     defined DISTRIBUTE   && defined SOLVE3D
      SUBROUTINE obs_depth (ng, tile, model)
!
!git $Id$
!svn $Id: obs_depth.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
!  Set observations depth (Zobs) to zero if not bounded in the current
!  parallel tile to facilitate global reduction in "mp_collect".
!
      DO iobs=Mstr,Mend
        r_bound=((rXmin(ng) .le.Xobs(iobs)).and.                        &
     &           (Xobs(iobs).lt.rXmax(ng))).and.                        &
     &          ((rYmin(ng) .le.Yobs(iobs)).and.                        &
     &           (Yobs(iobs).lt.rYmax(ng)))
        u_bound=((uXmin(ng) .le.Xobs(iobs)).and.                        &
     &           (Xobs(iobs).lt.uXmax(ng))).and.                        &
     &          ((uYmin(ng) .le.Yobs(iobs)).and.                        &
     &           (Yobs(iobs).lt.uYmax(ng)))
        v_bound=((vXmin(ng) .le.Xobs(iobs)).and.                        &
     &           (Xobs(iobs).lt.vXmax(ng))).and.                        &
     &          ((vYmin(ng) .le.Yobs(iobs)).and.                        &
     &           (Yobs(iobs).lt.vYmax(ng)))
        IF (.not.(r_bound.or.u_bound.or.v_bound)) THEN
          Zobs(iobs)=IniVal                              ! not bouded
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_depth
#else
      SUBROUTINE obs_depth
      RETURN
      END SUBROUTINE obs_depth
#endif
