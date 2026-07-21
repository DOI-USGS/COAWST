#include "cppdefs.h"
      MODULE ice_enthalpy_mod

#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the ice enthalpy (ti*hi) prior to advection.  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ice
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
      implicit none
!
      PRIVATE
      PUBLIC  ice_enthalpy
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_enthalpy (ng, tile, model)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_enthalpy_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold(ng), linew(ng),                     &
     &                        ICE(ng)%Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_enthalpy
!
!***********************************************************************
      SUBROUTINE ice_enthalpy_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              liold, linew,                       &
     &                              Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
#else
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
#endif
!
!  Local variables
!
      integer :: i, j

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Computes the ice enthalpy.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Si(i,j,liold,isEnth)=Si(i,j,liold,isHice)*                    &
     &                         Si(i,j,liold,isTice)
          Si(i,j,linew,isEnth)=Si(i,j,linew,isHice)*                    &
     &                         Si(i,j,linew,isTice)
          IF(Si(i,j,liold,isHice).le.min_hi(ng)) THEN
            Si(i,j,liold,isEnth)=0.0_r8
          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,liold,isEnth))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,linew,isEnth))
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,liold,isEnth),                         &
     &                    Si(:,:,linew,isEnth))
# endif
!
      RETURN
      END SUBROUTINE ice_enthalpy_tile
#endif
      END MODULE ice_enthalpy_mod
