      SUBROUTINE ana_m3clima (ng, tile, model)
!
!! svn $Id: ana_m3clima.h 737 2008-09-07 02:06:44Z jcwarner $
!!======================================================================
!! Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical 3D momentum climatology fields.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_m3clima_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng) % uclm,                          &
     &                       CLIMA(ng) % vclm)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(13)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_m3clima
!
!***********************************************************************
      SUBROUTINE ana_m3clima_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             uclm, vclm)
!***********************************************************************
!
      USE mod_param
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_3d_mod
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: uclm(LBi:,LBj:,:)
      real(r8), intent(out) :: vclm(LBi:,LBj:,:)
#else
      real(r8), intent(out) :: uclm(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: vclm(LBi:UBi,LBj:UBj,N(ng))
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
# ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
# else
      logical :: EWperiodic=.FALSE.
# endif
# ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
# else
      logical :: NSperiodic=.FALSE.
# endif
#endif
      integer :: i, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set 3D momentum climatology.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            uclm(i,j,k)=???
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vclm(i,j,k)=???
          END DO
        END DO
      END DO
#if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_u3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        uclm)
      CALL exchange_v3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        vclm)
#endif
#ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    uclm, vclm)
#endif
      RETURN
      END SUBROUTINE ana_m3clima_tile
