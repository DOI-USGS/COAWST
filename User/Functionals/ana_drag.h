      SUBROUTINE ana_drag (ng, tile, model)
!
!! svn $Id: ana_drag.h 795 2016-05-11 01:42:43Z arango $
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical, spatially varying bottom roughness    !
!  length (m), or linear drag coefficients (m/s), or quadratic drag    !
!  coefficients (nondimensional) at RHO-points. It  depends on  the    !
!  activated bottom stress formulation.                                !
!                                                                      !
!  There are many ways to compute spatially varying drag parameters:   !
!                                                                      !
!    * Partition the grid into different provinces with different      !
!      with different values (regimes).                                !
!    * A piecewise value that depends on the water depth.              !
!    * Empirical formulas in terms of water depth (Chezy formula)      !
!    * Inverse techniques using adjoint parameter estimation, but      !
!      it is beyond the scope of this routine.                         !
!                                                                      !
!  The User should experiment to get the appropriate distribution      !
!  for their application.                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_drag_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
#if defined SEDIMENT || defined BBL_MODEL
     &                    SEDBED(ng) % bottom,                          &
#endif
#if defined UV_LOGDRAG
     &                    GRID(ng) % ZoBot)
#elif defined UV_LDRAG
     &                    GRID(ng) % rdrag)
#elif defined UV_QDRAG
     &                    GRID(ng) % rdrag2)
#endif
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 2)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_drag
!
!***********************************************************************
      SUBROUTINE ana_drag_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#if defined SEDIMENT || defined BBL_MODEL
     &                          bottom,                                 &
#endif
#if defined UV_LOGDRAG
     &                          ZoBot)
#elif defined UV_LDRAG
     &                          rdrag)
#elif defined UV_QDRAG
     &                          rdrag2)
#endif
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
# if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(out) :: bottom(LBi:,LBj:,:)
# endif
# if defined UV_LOGDRAG
      real(r8), intent(out) :: ZoBot(LBi:,LBj:)
# elif defined UV_LDRAG
      real(r8), intent(out) :: rdrag(LBi:,LBj:)
# elif defined UV_QDRAG
      real(r8), intent(out) :: rdrag2(LBi:,LBj:)
# endif

#else

# if defined SEDIMENT || defined BBL_MODEL
      real(r8), intent(out) :: bottom(LBi:UBi,LBj:UBj,MBOTP)
# endif
# if defined UV_LOGDRAG
      real(r8), intent(out) :: ZoBot(LBi:UBi,LBj:UBj)
# elif defined UV_LDRAG
      real(r8), intent(out) :: rdrag(LBi:UBi,LBj:UBj)
# elif defined UV_QDRAG
      real(r8), intent(out) :: rdrag2(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
#if defined UV_LOGDRAG
!  Set spatially varying bottom roughness length (m).
#elif defined UV_LDRAG
!  Set spatially varying linear drag coefficient (m/s).
#elif defined UV_QDRAG
!  Set spatially varying quadratic drag coefficient (nondimensional)
#endif
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ZoBot(i,j)=???
        END DO
      END DO
# elif defined UV_LDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rdrag(i,j)=???
        END DO
      END DO
# elif defined UV_QDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rdrag2(i,j)=???
        END DO
      END DO
# endif
#else
      ana_drag.h: no values provided for either ZoBot, rdrag, or rdrag2
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
#if defined UV_LOGDRAG
     &                          ZoBot)
#elif defined UV_LDRAG
     &                          rdrag)
#elif defined UV_QDRAG
     &                          rdrag2)
#endif
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
# if defined UV_LOGDRAG
     &                    ZoBot)
# elif defined UV_LDRAG
     &                    rdrag)
# elif defined UV_QDRAG
     &                    rdrag2)
# endif
#endif

#if defined UV_LOGDRAG && (defined SEDIMENT || defined BBL_MODEL)
!
!-----------------------------------------------------------------------
!  Load bottom roughness length into bottom properties array.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          bottom(i,j,izdef)=ZoBot(i,j)
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE ana_drag_tile
