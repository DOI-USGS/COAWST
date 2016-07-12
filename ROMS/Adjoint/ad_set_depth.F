#include "cppdefs.h"
      MODULE ad_set_depth_mod
#if defined ADJOINT && defined SOLVE3D
!
!svn $Id: ad_set_depth.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  Currently, two vertical coordinate transformations are available    !
!  with various possible vertical stretching, C(s), functions, (see    !
!  routine "set_scoord.F" for details).                                !
!                                                                      !
!  BASIC STATE variables needed: NONE                                  !
!  Independent Variables: ad_Hz, ad_z_r, ad_z_w                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_set_depth, ad_set_depth_tile
# ifdef ADJUST_BOUNDARY
      PUBLIC  :: ad_set_depth_bry
# endif

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sediment
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
      CALL wclock_on (ng, iADM, 12)
# endif
      CALL ad_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % ad_h,                          &
# ifdef ICESHELF
     &                        GRID(ng) % zice,                          &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        SEDIMENTS(ng) % ad_bed,                   &
     &                        SEDIMEMTS(ng) % ad_bed_thick0,            &
# endif
     &                        COUPLING(ng) % Zt_avg1,                   &
     &                        COUPLING(ng) % ad_Zt_avg1,                &
     &                        GRID(ng) % ad_Hz,                         &
     &                        GRID(ng) % ad_z_r,                        &
     &                        GRID(ng) % ad_z_w)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
# endif
      RETURN
      END SUBROUTINE ad_set_depth

!
!***********************************************************************
      SUBROUTINE ad_set_depth_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              h, ad_h,                            &
# ifdef ICESHELF
     &                              zice,                               &
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                              ad_bed, bed_thick0,                 &
# endif
     &                              Zt_avg1, ad_Zt_avg1,                &
     &                              ad_Hz, ad_z_r, ad_z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      USE mod_sediment
# endif
!
      USE ad_exchange_2d_mod
      USE ad_exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d, ad_mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(in) :: ad_bed(LBi:,LBj:,:,:)
      real(r8), intent(inout):: ad_bed_thick0(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
      real(r8), intent(inout) :: ad_bed(LBi:UBi,LBj:UBj,Nbed,MBEDP)
      real(r8), intent(inout) :: ad_bed_thick0(LBi:UBi,LBj:UBi)
#  endif
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff, cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: adfac, ad_hinv, ad_hwater, ad_z_r0, ad_z_w0
      real(r8) :: ad_cff2_r, ad_cff2_w

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff2_r=0.0_r8
      ad_cff2_w=0.0_r8
      ad_z_r0=0.0_r8
      ad_z_w0=0.0_r8
      ad_hinv=0.0_r8
      ad_hwater=0.0_r8
!
!-----------------------------------------------------------------------
!  Compute time evolving adjoint depths and vertical thicknesses.
!-----------------------------------------------------------------------
!
!  Exchange boundary information.
!
# ifdef DISTRIBUTE
!>    CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_z_r, tl_Hz)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_z_r, ad_Hz)
!>    CALL mp_exchange3d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_z_w)
!>
      CALL ad_mp_exchange3d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_z_w)
!>    CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    tl_h)
!>
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_h)
!
# endif

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_Hz)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_Hz)
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          tl_z_r)
!>
        CALL ad_exchange_r3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 1, N(ng),        &
     &                             ad_z_r)
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          tl_z_w)
!>
        CALL ad_exchange_w3d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, 0, N(ng),        &
     &                             ad_z_w)
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          tl_h)
!>
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             ad_h)
      END IF
!
!-----------------------------------------------------------------------
!  Original formulation: Compute vertical depths (meters, negative) at
!                        RHO- and W-points, and vertical grid
!  thicknesses. Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!-----------------------------------------------------------------------
!
      IF (Vtransform(ng).eq.1) THEN
        DO j=JstrT,JendT
          DO k=N(ng),1,-1
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_r0=cff_r+cff1_r*hwater

!>            tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
!>
              ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+ad_Hz(i,j,k)
              ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
              ad_Hz(i,j,k)=0.0_r8
!>            tl_z_r(i,j,k)=tl_z_r0+                                    &
!>   &                      tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+         &
!>   &                      Zt_avg1(i,j)*(tl_z_r0*hinv+z_r0*tl_hinv)
!>
              adfac=Zt_avg1(i,j)*ad_z_r(i,j,k)
              ad_z_r0=ad_z_r0+hinv*adfac+ad_z_r(i,j,k)
              ad_hinv=ad_hinv+z_r0*adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                          &
     &                        (1.0_r8+z_r0*hinv)*ad_z_r(i,j,k)
              ad_z_r(i,j,k)=0.0_r8
!>            tl_z_r0=cff1_r*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_r*ad_z_r0
              ad_z_r0=0.0_r8
!>            tl_z_w(i,j,k)=tl_z_w0+                                    &
!>   &                      tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+         &
!>   &                      Zt_avg1(i,j)*(tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=Zt_avg1(i,j)*ad_z_w(i,j,k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_z_w(i,j,k)
              ad_hinv=ad_hinv+z_w0*adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+                          &
     &                        (1.0_r8+z_w0*hinv)*ad_z_w(i,j,k)
              ad_z_w(i,j,k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
!>            tl_hinv=-hinv*hinv*tl_hwater
!>
              ad_hwater=ad_hwater-hinv*hinv*ad_hinv
              ad_hinv=0.0_r8
!>            tl_hwater=tl_h(i,j)
!>
              ad_h(i,j)=ad_h(i,j)+ad_hwater
              ad_hwater=0.0_r8
            END DO
          END DO

          DO i=IstrT,IendT
!>          tl_z_w(i,j,0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_z_w(i,j,0)
            ad_z_w(i,j,0)=0.0_r8
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
!>            tl_h(i,j)=0.0_r8
!>
              ad_h(i,j)=0.0_r8
            END IF
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>          tl_h(i,j)=tl_h(i,j)-                                        &
!>   &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
!>
            ad_bed_thick(i,j,nnew)=ad_bed_thick(i,j,nnew)+              &
     &                             ad_h(i,j)

            ad_bed_thick(i,j,nstp)=ad_bed_thick(i,j,nstp)-              &
     &                             ad_h(i,j)
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform(ng).eq.2) THEN
        DO j=JstrT,JendT
          DO k=N(ng),1,-1
            cff_r=hc(ng)*SCALARS(ng)%sc_r(k)
            cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/(hc(ng)+hwater)
              cff2_r=(cff_r+cff1_r*hwater)*hinv
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
!>
              ad_z_w(i,j,k  )=ad_z_w(i,j,k  )+ad_Hz(i,j,k)
              ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
              ad_Hz(i,j,k)=0.0_r8
!>            tl_z_r(i,j,k)=tl_Zt_avg1(i,j)+                            &
!>   &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_r+         &
!>   &                      (Zt_avg1(i,j)+hwater)*tl_cff2_r
!>
              adfac=cff2_r*ad_z_r(i,j,k)
              ad_cff2_r=ad_cff2_r+                                      &
     &                  (Zt_avg1(i,j)+hwater)*ad_z_r(i,j,k)
              ad_hwater=ad_hwater+adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_z_r(i,j,k)+adfac
              ad_z_r(i,j,k)=0.0_r8
!>            tl_cff2_r=cff1_r*tl_hwater*hinv+                          &
!>   &                  (cff_r+cff1_r*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_r+cff1_r*hwater)*ad_cff2_r
              ad_hwater=ad_hwater+                                      &
     &                  cff1_r*hinv*ad_cff2_r
              ad_cff2_r=0.0_r8
!>            tl_z_w(i,j,k)=tl_Zt_avg1(i,j)+                            &
!>   &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_w+         &
!>   &                      (Zt_avg1(i,j)+hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_z_w(i,j,k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (Zt_avg1(i,j)+hwater)*ad_z_w(i,j,k)
              ad_hwater=ad_hwater+adfac
              ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_z_w(i,j,k)+adfac
              ad_z_w(i,j,k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
!>            tl_hinv=-hinv*hinv*tl_hwater
!>
              ad_hwater=ad_hwater-hinv*hinv*ad_hinv
              ad_hinv=0.0_r8
!>            tl_hwater=tl_h(i,j)
!>
              ad_h(i,j)=ad_h(i,j)+ad_hwater
              ad_hwater=0.0_r8
            END DO
          END DO

          DO i=IstrT,IendT
!>          tl_z_w(i,j,0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_z_w(i,j,0)
            ad_z_w(i,j,0)=0.0_r8
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
!>            tl_h(i,j)=0.0_r8
!>
              ad_h(i,j)=0.0_r8
            END IF
# endif
# if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
!>          tl_h(i,j)=tl_h(i,j)-                                        &
!>   &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
!>
            ad_bed_thick(i,j,nnew)=ad_bed_thick(i,j,nnew)+              &
     &                             ad_h(i,j)

            ad_bed_thick(i,j,nstp)=ad_bed_thick(i,j,nstp)-              &
     &                             ad_h(i,j)
# endif
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE ad_set_depth_tile

# ifdef ADJUST_BOUNDARY
!
!***********************************************************************
      SUBROUTINE ad_set_depth_bry (ng, tile)
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
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 12)
#  endif
      CALL ad_set_depth_bry_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % ad_h,                      &
#  ifdef ICESHELF
     &                            GRID(ng) % zice,                      &
#  endif
     &                            GRID(ng) % ad_Hz_bry)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 12)
#  endif
      RETURN
      END SUBROUTINE ad_set_depth_bry
!
!***********************************************************************
      SUBROUTINE ad_set_depth_bry_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj, LBij, UBij, &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  h, ad_h,                        &
#  ifdef ICESHELF
     &                                  zice,                           &
#  endif
     &                                  ad_Hz_bry)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d_bry
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: ad_h(LBi:,LBj:)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: ad_Hz_bry(LBij:,:,:)
#  else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_h(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: ad_Hz_bry(LBij:UBij,N(ng),4)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ibry, j, k

      real(r8) :: cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_w0
      real(r8) :: ad_cff2_w, ad_hinv, ad_hwater, ad_z_w0
      real(r8) :: adfac

      real(r8), dimension(0:N(ng)) :: ad_Zw

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
      ad_cff2_w=0.0_r8
      ad_z_w0=0.0_r8
      ad_hinv=0.0_r8
      ad_hwater=0.0_r8
      ad_Zw(0:N(ng))=0.0_r8

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      DO ibry=1,4
!>      CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ibry,                &
!>   &                          LBij, UBij, 1, N(ng),                   &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          tl_Hz_bry(:,:,ibry))
!>
        CALL ad_mp_exchange3d_bry (ng, tile, iTLM, 1, ibry,             &
     &                             LBij, UBij, 1, N(ng),                &
     &                             NghostPoints,                        &
     &                             EWperiodic(ng), NSperiodic(ng),      &
     &                             ad_Hz_bry(:,:,ibry))
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Original formulation: Compute vertical depths (meters, negative) at
!                        RHO- and W-points, and vertical grid
!  thicknesses. Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,t) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!-----------------------------------------------------------------------
!
      IF (Vtransform(ng).eq.1) THEN

        IF (ad_LBC(inorth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
!>            tl_Hz_bry(i,k,inorth)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,inorth)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,inorth)
              ad_Hz_bry(i,k,inorth)=0.0_r8
!>            tl_Zw(k)=tl_z_w0+                                         &
!>   &                 BOUNDARY(ng)%tl_zeta_north(i)*                   &
!>   &                 (1.0_r8+z_w0*hinv)+                              &
!>   &                 BOUNDARY(ng)%zeta_north(i)*                      &
!>   &                 (tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=BOUNDARY(ng)%zeta_north(i)*ad_Zw(k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_Zw(k)
              ad_hinv=ad_hinv+z_w0*adfac
              BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%               &
     &                                        ad_zeta_north(i)+         &
     &                                      (1.0_r8+z_w0*hinv)*ad_Zw(k)
              ad_Zw(k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(isouth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
!>            tl_Hz_bry(i,k,isouth)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,isouth)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,isouth)
              ad_Hz_bry(i,k,isouth)=0.0_r8
!>            tl_Zw(k)=tl_z_w0+                                         &
!>   &                 BOUNDARY(ng)%tl_zeta_south(i)*                   &
!>   &                 (1.0_r8+z_w0*hinv)+                              &
!>   &                 BOUNDARY(ng)%zeta_south(i)*                      &
!>   &                 (tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=BOUNDARY(ng)%zeta_south(i)*ad_Zw(k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_Zw(k)
              ad_hinv=ad_hinv+z_w0*adfac
              BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%               &
     &                                        ad_zeta_south(i)+         &
     &                                      (1.0_r8+z_w0*hinv)*ad_Zw(k)
              ad_Zw(k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(ieast,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
!>            tl_Hz_bry(i,k,ieast)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,ieast)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,ieast)
              ad_Hz_bry(i,k,ieast)=0.0_r8
!>            tl_Zw(k)=tl_z_w0+                                         &
!>   &                 BOUNDARY(ng)%tl_zeta_east(j)*                    &
!>   &                 (1.0_r8+z_w0*hinv)+                              &
!>   &                 BOUNDARY(ng)%zeta_east(j)*                       &
!>   &                 (tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=BOUNDARY(ng)%zeta_east(j)*ad_Zw(k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_Zw(k)
              ad_hinv=ad_hinv+z_w0*adfac
              BOUNDARY(ng)%ad_zeta_east(j)=BOUNDARY(ng)%                &
     &                                       ad_zeta_east(j)+           &
     &                                     (1.0_r8+z_w0*hinv)*ad_Zw(k)
              ad_Zw(k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(iwest,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
!>            tl_Hz_bry(i,k,iwest)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,iwest)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,iwest)
              ad_Hz_bry(i,k,iwest)=0.0_r8
!>            tl_Zw(k)=tl_z_w0+                                         &
!>   &                 BOUNDARY(ng)%tl_zeta_west(j)*                    &
!>   &                 (1.0_r8+z_w0*hinv)+                              &
!>   &                 BOUNDARY(ng)%zeta_west(j)*                       &
!>   ^                 (tl_z_w0*hinv+z_w0*tl_hinv)
!>
              adfac=BOUNDARY(ng)%zeta_west(j)*ad_Zw(k)
              ad_z_w0=ad_z_w0+hinv*adfac+ad_Zw(k)
              ad_hinv=ad_hinv+z_w0*adfac
              BOUNDARY(ng)%ad_zeta_west(j)=BOUNDARY(ng)%                &
     &                                       ad_zeta_west(j)+           &
     &                                     (1.0_r8+z_w0*hinv)*ad_Zw(k)
              ad_Zw(k)=0.0_r8
!>            tl_z_w0=cff1_w*tl_hwater
!>
              ad_hwater=ad_hwater+cff1_w*ad_z_w0
              ad_z_w0=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  New formulation: Compute vertical depths (meters, negative) at
!                   RHO- and W-points, and vertical grid thicknesses.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform(ng).eq.2) THEN

        IF (ad_LBC(inorth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz_bry(i,k,inorth)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,inorth)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,inorth)
              ad_Hz_bry(i,k,inorth)=0.0_r8
!>            tl_Zw(k)=BOUNDARY(ng)%tl_zeta_north(i)+                   &
!>   &                 (BOUNDARY(ng)%tl_zeta_north(i)+                  &
!>   &                  tl_hwater)*cff2_w+                              &
!>   &                 (BOUNDARY(ng)%zeta_north(i)+                     &
!>   &                  hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_Zw(k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (BOUNDARY(ng)%zeta_north(i)+hwater)*ad_Zw(k)
              ad_hwater=ad_hwater+adfac
              BOUNDARY(ng)%ad_zeta_north(i)=BOUNDARY(ng)%               &
     &                                        ad_zeta_north(i)+         &
     &                                      ad_Zw(k)+adfac
              ad_Zw(k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(isouth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz_bry(i,k,isouth)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,isouth)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,isouth)
              ad_Hz_bry(i,k,isouth)=0.0_r8
!>            tl_Zw(k)=BOUNDARY(ng)%tl_zeta_south(i)+                   &
!>   &                 (BOUNDARY(ng)%tl_zeta_south(i)+                  &
!>   &                  tl_hwater)*cff2_w+                              &
!>   &                 (BOUNDARY(ng)%zeta_south(i)+                     &
!>   &                  hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_Zw(k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (BOUNDARY(ng)%zeta_south(i)+hwater)*ad_Zw(k)
              ad_hwater=ad_hwater+adfac
              BOUNDARY(ng)%ad_zeta_south(i)=BOUNDARY(ng)%               &
     &                                        ad_zeta_south(i)+         &
     &                                      ad_Zw(k)+adfac
              ad_Zw(k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(ieast,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz_bry(i,k,ieast)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,ieast)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,ieast)
              ad_Hz_bry(i,k,ieast)=0.0_r8
!>            tl_Zw(k)=BOUNDARY(ng)%tl_zeta_east(j)+                    &
!>   &                 (BOUNDARY(ng)%tl_zeta_east(j)+                   &
!>   &                  tl_hwater)*cff2_w+                              &
!>   &                 (BOUNDARY(ng)%zeta_east(j)+                      &
!>   &                  hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_Zw(k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (BOUNDARY(ng)%zeta_east(j)+hwater)*ad_Zw(k)
              ad_hwater=ad_hwater+adfac
              BOUNDARY(ng)%ad_zeta_east(j)=BOUNDARY(ng)%                &
     &                                       ad_zeta_east(j)+           &
     &                                     ad_Zw(k)+adfac
              ad_Zw(k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
!
        IF (ad_LBC(iwest,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
!>            tl_Hz_bry(i,k,iwest)=tl_Zw(k)-tl_Zw(k-1)
!>
              ad_Zw(k-1)=ad_Zw(k-1)-ad_Hz_bry(i,k,iwest)
              ad_Zw(k  )=ad_Zw(k  )+ad_Hz_bry(i,k,iwest)
              ad_Hz_bry(i,k,iwest)=0.0_r8
!>            tl_Zw(k)=BOUNDARY(ng)%tl_zeta_west(j)+                    &
!>   &                 (BOUNDARY(ng)%tl_zeta_west(j)+                   &
!>   &                  tl_hwater)*cff2_w+                              &
!>   &                 (BOUNDARY(ng)%zeta_west(j)+                      &
!>   &                  hwater)*tl_cff2_w
!>
              adfac=cff2_w*ad_Zw(k)
              ad_cff2_w=ad_cff2_w+                                      &
     &                  (BOUNDARY(ng)%zeta_west(j)+hwater)*ad_Zw(k)
              ad_hwater=ad_hwater+adfac
              BOUNDARY(ng)%ad_zeta_west(j)=BOUNDARY(ng)%                &
     &                                       ad_zeta_west(j)+           &
     &                                     ad_Zw(k)+adfac
              ad_Zw(k)=0.0_r8
!>            tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
!>   &                  (cff_w+cff1_w*hwater)*tl_hinv
!>
              ad_hinv=ad_hinv+                                          &
     &                (cff_w+cff1_w*hwater)*ad_cff2_w
              ad_hwater=ad_hwater+                                      &
     &                  cff1_w*hinv*ad_cff2_w
              ad_cff2_w=0.0_r8
            END DO
!>          tl_Zw(0)=-tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)-ad_Zw(0)
            ad_Zw(0)=0.0_r8
!>          tl_hinv=-hinv*hinv*tl_hwater
!>
            ad_hwater=ad_hwater-hinv*hinv*ad_hinv
            ad_hinv=0.0_r8
!>          tl_hwater=tl_h(i,j)
!>
            ad_h(i,j)=ad_h(i,j)+ad_hwater
            ad_hwater=0.0_r8
          END DO
        END IF
      END IF

      RETURN
      END SUBROUTINE ad_set_depth_bry_tile
# endif
#endif
      END MODULE ad_set_depth_mod
