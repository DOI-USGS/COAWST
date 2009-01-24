#include "cppdefs.h"
      MODULE tl_set_depth_mod
#ifdef TANGENT
!
!svn $Id: tl_set_depth.F 694 2008-08-08 18:33:05Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  BASIC STATE variables needed: NONE                                  !
!  Independent Variables: tl_Hz, tl_z_r, tl_z_w                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef SOLVE3D
      PUBLIC  :: tl_set_depth, tl_set_depth_tile
# endif
      PUBLIC  :: tl_bath, tl_bath_tile

      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE tl_set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
      CALL tl_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % tl_h,                          &
#  ifdef ICESHELF
     &                        GRID(ng) % zice,                          &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        GRID(ng) % tl_bed_thick,                  &
#  endif
     &                        COUPLING(ng) % Zt_avg1,                   &
     &                        COUPLING(ng) % tl_Zt_avg1,                &
     &                        GRID(ng) % tl_Hz,                         &
     &                        GRID(ng) % tl_z_r,                        &
     &                        GRID(ng) % tl_z_w)
      RETURN
      END SUBROUTINE tl_set_depth
!
!***********************************************************************
      SUBROUTINE tl_set_depth_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              h, tl_h,                            &
#  ifdef ICESHELF
     &                              zice,                               &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                              tl_bed_thick,                       &
#  endif
     &                              Zt_avg1, tl_Zt_avg1,                &
     &                              tl_Hz, tl_z_r, tl_z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
      USE exchange_3d_mod
#  endif
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: tl_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout):: bed_thick0(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(in) :: tl_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: tl_h(LBi:,LBj:)
      real(r8), intent(out) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: tl_z_w(LBi:,LBj:,0:)
#  else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: tl_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout):: bed_thick0(LBi:UBi,LBj:UBi)
#   endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_h(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  endif
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
#   ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#   else
      logical :: EWperiodic=.FALSE.
#   endif
#   ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#   else
      logical :: NSperiodic=.FALSE.
#   endif
#  endif
      integer :: i, j, k

      real(r8) :: cff, cff_r, cff1_r, cff_w, cff1_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: tl_hinv, tl_hwater, tl_z_r0, tl_z_w0

#  ifdef WET_DRY
      real(r8), parameter :: eps = 1.0E-14_r8
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time evolving depths and vertical thicknesses.
!-----------------------------------------------------------------------
!
!  Compute vertical depths (meters, negative) at RHO- and W-points,
!  and vertical grid thicknesses.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
# if defined SEDIMENT && defined SED_MORPH
!>        h(i,j)=h(i,j)-bed_thick(i,j,nstp)+bed_thick(i,j,nnew)
!>
          tl_h(i,j)=tl_h(i,j)-                                          &
     &              tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
# endif
# if defined WET_DRY
          IF (h(i,j).eq.0.0_r8) THEN
!>          h(i,j)=eps
!>
            tl_h(i,j)=0.0_r8           
          END IF
# endif
!>        z_w(i,j,0)=-h(i,j)
!>
          tl_z_w(i,j,0)=-tl_h(i,j)
        END DO
        DO k=1,N(ng)
          cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
          cff1_w=SCALARS(ng)%Cs_w(k)
          cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
          cff1_r=SCALARS(ng)%Cs_r(k)
          DO i=IstrR,IendR
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/hwater
            tl_hinv=-hinv*hinv*tl_hwater

            z_w0=cff_w+cff1_w*hwater
            tl_z_w0=cff1_w*tl_hwater

!>          z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
!>
            tl_z_w(i,j,k)=tl_z_w0+                                      &
     &                    tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+           &
     &                    Zt_avg1(i,j)*(tl_z_w0*hinv+z_w0*tl_hinv)

            z_r0=cff_r+cff1_r*hwater
            tl_z_r0=cff1_r*tl_hwater

!>          z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
!>
            tl_z_r(i,j,k)=tl_z_r0+                                      &
     &                    tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+           &
     &                    Zt_avg1(i,j)*(tl_z_r0*hinv+z_r0*tl_hinv)
#  ifdef ICESHELF
!>          z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
!>          z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
#  endif
!>          Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
!>
            tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
          END DO
        END DO
      END DO
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#   if defined EW_PERIODIC || defined NS_PERIODIC
!>    CALL exchange_r2d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj,                       &
!>   &                        h)
!>
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_h)
!>    CALL exchange_w3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
!>   &                        z_w)
!>
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        tl_z_w)
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        z_r)
!>
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_z_r)
!>    CALL exchange_r3d_tile (ng, tile,                                 &
!>   &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
!>   &                        Hz)
!>
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        tl_Hz)
#   endif
#   ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    h)
!>
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_h)
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    z_w)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_z_w)
!>    CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints, EWperiodic, NSperiodic,         &
!>   &                    z_r, Hz)
!>
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_z_r, tl_Hz)
#   endif
#  endif
      RETURN
      END SUBROUTINE tl_set_depth_tile
# endif
!
!***********************************************************************
      SUBROUTINE tl_bath (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL tl_bath_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   GRID(ng) % tl_h)
      RETURN
      END SUBROUTINE tl_bath
!
!***********************************************************************
      SUBROUTINE tl_bath_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         tl_h)
!***********************************************************************
!
      USE mod_param
!
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tl_h(LBi:,LBj:)
# else
      real(r8), intent(out) :: tl_h(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry to zero.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          tl_h(i,j)=0.0_r8
        END DO
      END DO

# if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        tl_h)
# endif
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    tl_h)
# endif
      RETURN
      END SUBROUTINE tl_bath_tile
#endif
      END MODULE tl_set_depth_mod
