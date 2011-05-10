      SUBROUTINE ana_wwave (ng, tile, model)
!
!! svn $Id: ana_wwave.h 429 2009-12-20 17:30:26Z arango $
!!======================================================================
!! Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets wind induced wave amplitude, direction and     !
!  period to be used in the bottom boundary layer formulation.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_wwave_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#if defined BBL_MODEL || defined WEC
     &                     FORCES(ng) % Dwave,                          &
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
     &                     FORCES(ng) % Dissip_fric,                    &
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS
     &                     FORCES(ng) % Dissip_break,                   &
     &                     FORCES(ng) % Dissip_wcap,                    &
#endif
#ifdef WAVES_HEIGHT
     &                     FORCES(ng) % Hwave,                          &
#endif
#ifdef WAVES_LENGTH
     &                     FORCES(ng) % Lwave,                          &
#endif
#ifdef WAVES_LENGTHP
     &                     FORCES(ng) % Lwavep,                         &
#endif
#ifdef WAVES_TOP_PERIOD
     &                     FORCES(ng) % Pwave_top,                      &
#endif
#ifdef WAVES_BOT_PERIOD
     &                     FORCES(ng) % Pwave_bot,                      &
#endif
#ifdef WAVES_UB
     &                     FORCES(ng) % Uwave_rms,                      &
#endif
#ifdef TKE_WAVEDISS
     &                     FORCES(ng) % wave_dissip,                    &
#endif
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(37)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_wwave
!
!***********************************************************************
      SUBROUTINE ana_wwave_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
#if defined BBL_MODEL || defined WEC
     &                           Dwave,                                 &
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
     &                           Dissip_fric,                           &
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WAVES_DISS
     &                           Dissip_break, Dissip_wcap,             &
#endif
#ifdef WAVES_HEIGHT
     &                           Hwave,                                 &
#endif
#ifdef WAVES_LENGTH
     &                           Lwave,                                 &
#endif
#ifdef WAVES_LENGTHP
     &                           Lwavep,                                &
#endif
#ifdef WAVES_TOP_PERIOD
     &                           Pwave_top,                             &
#endif
#ifdef WAVES_BOT_PERIOD
     &                           Pwave_bot,                             &
#endif
#ifdef WAVES_UB
     &                           Uwave_rms,                             &
#endif
#ifdef TKE_WAVEDISS
     &                           Dissip_break,                          &
#endif
     &                           angler, h)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#endif
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
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
# if defined BBL_MODEL || defined WEC
      real(r8), intent(inout) :: Dwave(LBi:,LBj:)
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      real(r8), intent(inout) :: Dissip_fric(LBi:,LBj:)
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
     defined WAVES_DISS
      real(r8), intent(inout) :: Dissip_break(LBi:,LBj:)
      real(r8), intent(inout) :: Dissip_wcap(LBi:,LBj:)
# endif
# ifdef WAVES_HEIGHT
      real(r8), intent(inout) :: Hwave(LBi:,LBj:)
# endif
# ifdef WAVES_LENGTH
      real(r8), intent(inout) :: Lwave(LBi:,LBj:)
# endif
# ifdef WAVES_LENGTHP
      real(r8), intent(inout) :: Lwavep(LBi:,LBj:)
# endif
# ifdef WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:,LBj:)
# endif
# ifdef WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:,LBj:)
# endif
# ifdef WAVES_UB
      real(r8), intent(inout) :: Uwave_rms(LBi:,LBj:)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(inout) :: wave_dissip(LBi:,LBj:)
# endif

#else

      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# if defined BBL_MODEL || defined WEC
      real(r8), intent(inout) :: Dwave(LBi:UBi,LBj:UBj)
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      real(r8), intent(inout) :: Dissip_fric(LBi:UBi,LBj:UBj)
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
     defined WAVES_DISS
      real(r8), intent(inout) :: Dissip_break(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Dissip_wcap(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_HEIGHT
      real(r8), intent(inout) :: Hwave(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_LENGTH
      real(r8), intent(inout) :: Lwave(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_LENGTHP
      real(r8), intent(inout) :: Lwavep(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_UB
      real(r8), intent(inout) :: Uwave_rms(LBi:UBi,LBj:UBj)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(inout) :: wave_dissip(LBi:UBi,LBj:UBj)
# endif
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
      integer :: i, j
      real(r8) :: cff, wdir
#if defined LAKE_SIGNELL
      real(r8) :: cff1, mxst, ramp_u, ramp_time, ramp_d
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set wind induced wave amplitude (m), direction (radians) and
!  period (s) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined LENTZ_TEST
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          Dwave(i,j)=180.0_r8*deg2rad
!         Dissip_fric(i,j)=0.0_r8
          Dissip_break(i,j)=0.0_r8
          Dissip_wcap(i,j)=0.0_r8
          Hwave(i,j)=2.0_r8
          Lwave(i,j)=20.0_r8
          Pwave_top(i,j)=7.0_r8
        END DO
      END DO
#else
      ana_wwave: No values provided for Hwave, Dwave, Pwave, Lwave.
#endif
#if defined EW_PERIODIC || defined NS_PERIODIC
# if defined WAVES_DIR
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dwave)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Dwave)
#  endif
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dissip_fric)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Dissip_fric)
#  endif
# endif
# if defined TKE_WAVEDISS || defined WAVES_OCEAN || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
     defined WAVES_DISS
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dissip_break)
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Dissip_wcap)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Dissip_break, Dissip_wcap)
#  endif
# endif
# ifdef WAVES_HEIGHT
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Hwave)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Hwave)
#  endif
# endif
# ifdef WAVES_LENGTH
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lwave)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Lwave)
#  endif
# endif
# ifdef WAVES_LENGTHP
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Lwavep)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Lwavep)
#  endif
# endif
# ifdef WAVES_TOP_PERIOD
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pwave_top)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Pwave_top)
#  endif
# endif
# ifdef WAVES_BOT_PERIOD
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Pwave_bot)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Pwave_bot)
#  endif
# endif
# ifdef WAVES_UB
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Uwave_rms)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    Uwave_rms)
#  endif
# endif
# ifdef TKE_WAVEDISS
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        wave_dissip)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    wave_dissip)
#  endif
# endif
#endif
      RETURN
      END SUBROUTINE ana_wwave_tile
