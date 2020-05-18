      SUBROUTINE ana_ncep (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine computes values for NCEP-type surface fluxes           !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
      implicit none

      integer, intent(in) :: ng, tile, model

#include "tile.h"

      CALL ana_ncep_tile (ng, tile, model,                              &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       FORCES(ng) % nustr,                        &
     &                       FORCES(ng) % nvstr,                        &
     &                       FORCES(ng) % cloud,                        &
     &                       FORCES(ng) % srflx,                        &
     &                       FORCES(ng) % lrflx,                        &
     &                       FORCES(ng) % shflx,                        &
     &                       FORCES(ng) % lhflx,                        &
     &                       FORCES(ng) % Pair,                         &
     &                       FORCES(ng) % runoff,                       &
     &                       FORCES(ng) % rain,                         &
     &                       FORCES(ng) % skt,                          &
     &                       FORCES(ng) % icec                          &
     &                   )
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(47)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_ncep
!
!***********************************************************************
       SUBROUTINE ana_ncep_tile (ng, tile, model, LBi, UBi, LBj, UBj,   &
     &                             nustr, nvstr, cloud, srflx, lrflx,   &
     &                             shflx, lhflx, Pair, runoff, rain,    &
     &                             skt, icec                            &
     &                             )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#if defined DISTRIBUTE
      USE USE mp_exchange_mod
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: nustr(LBi:,LBj:)
      real(r8), intent(inout) :: nvstr(LBi:,LBj:)
      real(r8), intent(inout) :: cloud(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: Pair(LBi:,LBj:)
      real(r8), intent(inout) :: runoff(LBi:,LBj:)
      real(r8), intent(inout) :: rain(LBi:,LBj:)
      real(r8), intent(inout) :: skt(LBi:,LBj:)
      real(r8), intent(inout) :: icec(LBi:,LBj:)
#else
      real(r8), intent(inout) :: nustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: nvstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: cloud(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: runoff(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rain(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: skt(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: icec(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      logical :: EWperiodic=.FALSE.
      logical :: NSperiodic=.FALSE. 

#include "set_bounds.h"

      DO j=JstrT,JendT
        DO i=IstrT,IendT
           nustr(i,j) = 0.2_r8
           nvstr(i,j) = 0.2_r8
           cloud(i,j) = 0.3_r8
           srflx(i,j) = 0._r8
           lrflx(i,j) = 40._r8
           shflx(i,j) = -6._r8
           lhflx(i,j) = 0.3_r8
           Pair(i,j) = 1.02E+5_r8
           runoff(i,j) = 0.0_r8
           rain(i,j) = 1.0E-7_r8
           skt(i,j) = 235._r8
           icec(i,j) = 1.0_r8
        ENDDO
      ENDDO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, nustr)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, nvstr)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, cloud)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, lrflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, Pair)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, runoff)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, rain)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, skt)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, icec)
      END IF
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    nustr, nvstr, cloud, srflx)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    lrflx, shflx, lhflx, Pair)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    runoff, rain, skt, icec)
#endif
      RETURN
      END SUBROUTINE ana_ncep_tile
