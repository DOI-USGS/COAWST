      SUBROUTINE ana_albedo (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets surface albedo for the ocean and optionally,      !
!  for the ice.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_albedo_tile (ng, tile, model,                            &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                     GRID(ng) % latr,                             &   
# endif
# ifdef ICE_MODEL
     &                     FORCES(ng) % albedo_ice,                     &
# endif
     &                     FORCES(ng) % albedo)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME( 9)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_albedo
!
!***********************************************************************
      SUBROUTINE ana_albedo_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                            latr,                                 &
# endif
# ifdef ICE_MODEL
     &                            albedo_ice,                           &
# endif
     &                            albedo)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
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
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(out) :: albedo_ice(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: albedo(LBi:,LBj:)
#else
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(out) :: albedo_ice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: albedo(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: iday, month, year
      real(r8) :: hour, yday
      real(r8), parameter :: alb(12) =                                &
     &           (/ .85, .85, .83, .81, .82, .78,                     &
     &              .64, .69, .84, .85, .85, .85 /)
      real(r8), parameter :: alb_w=0.06_r8

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical surface albedo.
!-----------------------------------------------------------------------
!
      CALL caldate(r_date, tdays(ng), year, yday, month, iday, hour)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
#ifdef ICE_MODEL
          albedo_ice(i,j)=alb(month)
#endif
#ifdef ALBEDO_CURVE
# ifdef BIO_1D
!using lat for M2 for whole domain
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*56.877))
# else
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*latr(i,j)))
# endif
#else
          albedo(i,j)=alb_w
#endif
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#ifdef ICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ice)
#endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo)
      END IF

#ifdef DISTRIBUTE
# ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo_ice)
# endif
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo)
#endif

      RETURN
      END SUBROUTINE ana_albedo_tile
