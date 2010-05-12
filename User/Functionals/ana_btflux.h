      SUBROUTINE ana_btflux (ng, tile, model, itrc)
!
!! svn $Id: ana_btflux.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets kinematic bottom flux of tracer type variables    !
!  (tracer units m/s).                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc

#include "tile.h"
!
      CALL ana_btflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % btflx)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 3)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_btflux
!
!***********************************************************************
      SUBROUTINE ana_btflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
#else
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic bottom heat flux (degC m/s) at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
#if defined MY_APPLICATION
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=???
          END DO
        END DO
#else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s) at horizontal RHO-points,
!  scaling by bottom salinity is done elsewhere.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
#if defined MY_APPLICATION
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=???
          END DO
        END DO
#else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set kinematic bottom flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
#if defined MY_APPLICATION
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=???
          END DO
        END DO
      END IF
#else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
#endif

      RETURN
      END SUBROUTINE ana_btflux_tile
