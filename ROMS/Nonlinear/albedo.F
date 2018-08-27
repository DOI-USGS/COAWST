#include "cppdefs.h"
      MODULE albedo_mod
#if defined ALBEDO && !defined ALBEDO_FILE && !defined ANA_ALBEDO
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the albedo                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: albedo_eval

contains
!
!***********************************************************************
      SUBROUTINE albedo_eval (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17)
# endif
      CALL albedo_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef ICE_MODEL
     &                  liold(ng),                                      &
     &                  linew(ng),                                      &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                  GRID(ng) % latr,                                &
# endif
# ifdef ICE_MODEL
     &                  ICE(ng) % ai,                                   &
     &                  ICE(ng) % hi,                                   &
     &                  FORCES(ng) % albedo_ice,                        &
#  ifdef ICE_THERMO
     &                  ICE(ng) % hsn,                                  &
     &                  ICE(ng) % tis,                                  &
#   ifdef MELT_PONDS
     &                  ICE(ng) % apond,                                &
     &                  ICE(ng) % hpond,                                &
#   endif
#  endif
# endif
     &                  FORCES(ng) % albedo                             &
     &                  )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17)
# endif
      RETURN
      END SUBROUTINE albedo_eval
!
!***********************************************************************
      SUBROUTINE albedo_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef ICE_MODEL
     &                        liold, linew,                             &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                        latr,                                     &
# endif
# ifdef ICE_MODEL
     &                        ai, hi, albedo_ice,                       &
#  ifdef ICE_THERMO
     &                        hsn, tis,                                 &
#   ifdef MELT_PONDS
     &                        apond, hpond,                             &
#   endif
#  endif
# endif
     &                        albedo                                    &
     &                        )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined BEST_NPZ
      USE mod_biology
# endif
!
# if defined BEST_NPZ && defined CLIM_ICE_1D
      USE mod_clima
# endif

      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ICE_MODEL
      integer, intent(in) :: liold
      integer, intent(in) :: linew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(out) :: albedo_ice(LBi:,LBj:)
#   ifdef ICE_THERMO
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
#    ifdef MELT_PONDS
      real(r8), intent(in) :: apond(LBi:,LBj:,:)
      real(r8), intent(in) :: hpond(LBi:,LBj:,:)
#    endif
#   endif
#  endif
      real(r8), intent(out) :: albedo(LBi:,LBj:)

# else
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: albedo_ice(LBi:UBi,LBj:UBj)
#   ifdef ICE_THERMO
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
#    ifdef MELT_PONDS
      real(r8), intent(in) :: apond(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hpond(LBi:UBi,LBj:UBj,2)
#    endif
#   endif
#  endif
      real(r8), intent(out) :: albedo(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, listp

      real(r8) :: cff1, cff2
      real(r8) :: cff
# ifdef ICE_BULK_FLUXES
#  ifdef ALBEDO_CSIM
      real(r8), parameter :: alb_i_thick=0.54_r8
      real(r8) :: alb_i_dry
      real(r8) :: alb_ice, alb_snow, fh, fsn
      real(r8), parameter :: alb_s_dry=0.83_r8
      real(r8), parameter :: alb_s_wet=0.70_r8
#  else
#ifdef ICE_BOX
      real(r8), parameter :: alb_i_dry=0.75_r8
      real(r8), parameter :: alb_i_wet=0.64_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.82_r8
#else
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
#endif
#  endif
      real(r8) :: albs, qlwi, qlh_i, qsh_i
      real(r8) :: le_i, dq_i,fqlat1, slp, Qsati
      real(r8) :: vap_p_i
#  ifdef ICE_ALB_EC92
      real(r8) :: albi, albsn, thki_n, thksn_n
#  endif
# endif
      real(r8), parameter :: alb_w=0.06_r8
# if defined ICE_MODEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
# endif

# include "set_bounds.h"
!-------------------------------------------------------------------------------
! PURPOSE:
!   computes albedo over snow/ice/water
!-------------------------------------------------------------------------------

# if defined ICE_MODEL
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        listp = liold
      ELSE
        listp = linew
      END IF
# endif

! Note that this loop needs to be cleaned of all global arrays for
! OpenMP.
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# if defined ICE_MODEL
! Calculate the ice/snow albedo
          ice_thick(i,j) = hi(i,j,listp)/(ai(i,j,listp)+0.001)
          snow_thick(i,j) = hsn(i,j,listp)/(ai(i,j,listp)+0.001)
!
#  ifdef ICE_ALB_EC92
! Ice and snow albedo is calculated from Ebert and Curry,1992b
          albi=0._r8
          albsn=0._r8
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            thki_n = ice_thick(i,j)
            thki_n = MAX(thki_n,0.00001_r8)
            thksn_n = snow_thick(i,j)
            albi=0.082409_r8*LOG(thki_n)+0.485472_r8
            IF (thki_n.GE.1._r8) albi=0.07616_r8*thki_n+0.414492_r8
            IF (thki_n.GE.2._r8) albi=0.561632_r8
!approximated values for albsn(depends on COSZ, but small variation)
            albsn=0.83_r8
#   ifdef MELT_PONDS
            IF (apond(i,j,listp).GT.0._r8) THEN
              IF (thksn_n.GE.0.1_r8) THEN
                albsn=0.701009_r8
              ELSE
                albsn=albi+(0.701009_r8-albi)*thksn_n/0.1
              ENDIF
            ENDIF
#   endif
            albedo_ice(i,j)=albi
            IF (hsn(i,j,listp).GT.0._r8) albedo_ice(i,j)=albsn
!            IF (sfwat(i,j,listp).GT.0._r8) albs=0.10737_r8              &
!     &         +0.518_r8*EXP(-8.1_r8 *sfwat(i,j,listp)-0.47_r8)         &
!     &         +0.341_r8*EXP(-31.8_r8*sfwat(i,j,listp)-0.94_r8)         &
!     &         +0.131_r8*EXP(-2.6_r8 *sfwat(i,j,listp)-3.82_r8)
!            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  elif defined ALBEDO_CSIM
          fh = min(atan(4.0*ice_thick(i,j))/atan(2.0), 1.0)
          fsn = snow_thick(i,j) / (snow_thick(i,j) + 0.02)
          alb_i_dry = alb_w*(1-fh) + alb_i_thick*fh
          cff1 = alb_s_wet - alb_s_dry
          cff2 = -0.075      ! alb_i_wet - alb_i_dry
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (tis(i,j) .gt. -1.0_r8) THEN
              alb_snow = cff1*(tis(i,j)+1.0_r8)+alb_s_dry
            ELSE
              alb_snow = alb_s_dry
            ENDIF
            IF (tis(i,j) .gt. -1.0_r8) THEN
              alb_ice = cff2*(tis(i,j)+1.0_r8)+alb_i_dry
            ELSE
              alb_ice = alb_i_dry
            ENDIF
            albedo_ice(i,j) = fsn*alb_snow + (1-fsn)*alb_ice
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  elif defined ICE_BOX
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (hsn(i,j,listp).gt.0._r8) THEN
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_s_wet
              ELSE
                albedo_ice(i,j) = alb_s_dry
              ENDIF
            ELSE
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_i_wet
              ELSE
                albedo_ice(i,j) = alb_i_dry
              ENDIF
            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (hsn(i,j,listp).gt.0._r8) THEN
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = cff1*(tis(i,j)+1.0_r8)+alb_s_dry
              ELSE
                albedo_ice(i,j) = alb_s_dry
              ENDIF
            ELSE
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = cff2*(tis(i,j)+1.0_r8)+alb_i_dry
              ELSE
                albedo_ice(i,j) = alb_i_dry
              ENDIF
            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  endif
# endif
! Compute ocean albedo
# ifdef ALBEDO_CURVE
#  ifdef BIO_1D
!using lat for M2 for whole domain
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*56.877))
#  else
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*latr(i,j)))
#  endif
# else
          albedo(i,j)=alb_w
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo)
# ifdef ICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ice)
# endif
      END IF
# ifdef DISTRIBUTE
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo, albedo_ice)
#  else
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo)
#  endif
# endif

      RETURN
      END SUBROUTINE albedo_tile
#endif

      END module albedo_mod
