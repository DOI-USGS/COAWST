      SUBROUTINE ana_btflux (ng, tile, model, itrc)
!
!! git $Id$
!! svn $Id: ana_btflux.h 1151 2023-02-09 03:08:53Z arango $
!!======================================================================
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Sets bottom flux of tracer type variables btflux(:,:,itrc) using    !
!  analytical expressions (TracerUnits m/s).  The surface fluxes are   !
!  processed and loaded to state variable "btflx" in "set_vbc".        !
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
!
! Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
      CALL ana_btflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % btflux)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 3)=MyFile
      END IF
!
      RETURN
      END SUBROUTINE ana_btflux
!
!***********************************************************************
      SUBROUTINE ana_btflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            btflux)
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
      real(r8), intent(inout) :: btflux(LBi:,LBj:,:)
#else
      real(r8), intent(inout) :: btflux(LBi:UBi,LBj:UBj,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set bottom heat flux (degC m/s) at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
#if defined MY_APPLICATION
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=???
          END DO
        END DO
#else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=0.0_r8
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set bottom salt flux (m/s) at horizontal RHO-points. The scaling
!  by bottom salinity is done in "set_vbc".
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
#if defined MY_APPLICATION
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=???
          END DO
        END DO
#else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=0.0_r8
          END DO
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Set bottom flux (Tunits m/s) of passive tracers at RHO-point,
!  if any.
!-----------------------------------------------------------------------
!
      ELSE
#if defined MY_APPLICATION
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=???
          END DO
        END DO
      END IF
#else
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflux(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
#endif
!
      RETURN
      END SUBROUTINE ana_btflux_tile
