      SUBROUTINE ana_m3clima (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical 3D momentum climatology fields.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m3clima_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(13)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m3clima
!
!***********************************************************************
      SUBROUTINE ana_m3clima_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_scalars
!
      USE exchange_3d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: i, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology.
!-----------------------------------------------------------------------
!
      IF (Lm3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrP,IendT
              CLIMA(ng)%uclm(i,j,k)=???
            END DO
          END DO
          DO j=JstrP,JendT
            DO i=IstrT,IendT
              CLIMA(ng)%vclm(i,j,k)=???
            END DO
          END DO
        END DO
!
!  Exchange boundary data.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng) % uclm)
          CALL exchange_v3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng) % vclm)
        END IF

#ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 2,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % uclm,                           &
     &                      CLIMA(ng) % vclm)
#endif
      END IF

      RETURN
      END SUBROUTINE ana_m3clima_tile
