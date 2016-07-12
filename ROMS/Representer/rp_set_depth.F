#include "cppdefs.h"
      MODULE rp_set_depth_mod
#ifdef TL_IOMS
!
!svn $Id: rp_set_depth.F 795 2016-05-11 01:42:43Z arango $
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
!  Independent Variables: tl_Hz, tl_z_r, tl_z_w                        !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
# ifdef SOLVE3D
      PUBLIC  :: rp_set_depth, rp_set_depth_tile
#  ifdef ADJUST_BOUNDARY
      PUBLIC  :: rp_set_depth_bry
#  endif
# endif
      PUBLIC  :: rp_bath, rp_bath_tile
!
      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE rp_set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
#  if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
#  endif
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
#  ifdef PROFILE
      CALL wclock_on (ng, iRPM, 12)
#  endif
      CALL rp_set_depth_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % h,                             &
     &                        GRID(ng) % tl_h,                          &
#  ifdef ICESHELF
     &                        GRID(ng) % zice,                          &
#  endif
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
     &                        SEDBED(ng) % tl_bed_thick,                &
#  endif
     &                        COUPLING(ng) % Zt_avg1,                   &
     &                        COUPLING(ng) % tl_Zt_avg1,                &
     &                        GRID(ng) % tl_Hz,                         &
     &                        GRID(ng) % tl_z_r,                        &
     &                        GRID(ng) % tl_z_w)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
# endif

      RETURN
      END SUBROUTINE rp_set_depth
!
!***********************************************************************
      SUBROUTINE rp_set_depth_tile (ng, tile,                           &
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
      USE exchange_2d_mod
      USE exchange_3d_mod
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
      real(r8), intent(in):: tl_bed_thick(LBi:,LBj:,:)
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
      real(r8), intent(in):: tl_bed_thick(LBi:UBi,LBj:UBi,2)
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
      integer :: i, j, k

      real(r8) :: cff, cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: tl_cff2_r, tl_cff2_w
      real(r8) :: tl_hinv, tl_hwater, tl_z_r0, tl_z_w0

#  ifdef WET_DRY
      real(r8), parameter :: eps = 1.0E-14_r8
#  endif

#  include "set_bounds.h"
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
          DO i=IstrT,IendT
#  if defined SEDIMENT && defined SED_MORPH
!>          h(i,j)=h(i,j)-bed_thick(i,j,nstp)+bed_thick(i,j,nnew)
!>
            tl_h(i,j)=tl_h(i,j)-                                        &
     &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
#  endif
!>          z_w(i,j,0)=-h(i,j)
!>
            tl_z_w(i,j,0)=-tl_h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
#  ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
#  endif
              tl_hwater=tl_h(i,j)
              hinv=1.0_r8/hwater
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_hinv=-hinv*hinv*tl_hwater+                             &
#   ifdef TL_IOMS
     &                2.0_r8*hinv
#   endif
#  endif
              z_w0=cff_w+cff1_w*hwater
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w0=cff1_w*tl_hwater+                                 &
#   ifdef TL_IOMS
     &                cff_w
#   endif
#  endif
!>            z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
!>
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w(i,j,k)=tl_z_w0+                                    &
     &                      tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+         &
     &                      Zt_avg1(i,j)*z_w0*tl_hinv+                  &
     &                      Zt_avg1(i,j)*tl_z_w0*hinv-                  &
#   ifdef TL_IOMS
     &                      2.0_r8*Zt_avg1(i,j)*z_w0*hinv
#   endif
#  else
              tl_z_w(i,j,k)=tl_Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)+         &
#   ifdef TL_IOMS
     &                      z_w0
#   endif
#  endif
              z_r0=cff_r+cff1_r*hwater
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_r0=cff1_r*tl_hwater+                                 &
#   ifdef TL_IOMS
     &                cff_r
#   endif
#  endif
!>            z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
!>
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_r(i,j,k)=tl_z_r0+                                    &
     &                      tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+         &
     &                      Zt_avg1(i,j)*z_r0*tl_hinv+                  &
     &                      Zt_avg1(i,j)*tl_z_r0*hinv-                  &
#   ifdef TL_IOMS
     &                      2.0_r8*Zt_avg1(i,j)*z_r0*hinv
#   endif
#  else
              tl_z_r(i,j,k)=tl_Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)+         &
#   ifdef TL_IOMS
     &                      z_r0
#   endif
#  endif
#  ifdef ICESHELF
!>            z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
!>            z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
#  endif
!>            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
!>
              tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
            END DO
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
          DO i=IstrT,IendT
#  if defined SEDIMENT && defined SED_MORPH
!>          h(i,j)=h(i,j)-bed_thick(i,j,nstp)+bed_thick(i,j,nnew)
!>
            tl_h(i,j)=tl_h(i,j)-                                        &
     &                tl_bed_thick(i,j,nstp)+tl_bed_thick(i,j,nnew)
#  endif
!>          z_w(i,j,0)=-h(i,j)
!>
            tl_z_w(i,j,0)=-tl_h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*SCALARS(ng)%sc_r(k)
            cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
#  ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
#  endif
              tl_hwater=tl_h(i,j)
              hinv=1.0_r8/(hc(ng)+hwater)
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_hinv=-hinv*hinv*tl_hwater+                             &
#   ifdef TL_IOMS
     &                2.0_r8*hinv
#   endif
#  endif

              cff2_w=(cff_w+cff1_w*hwater)*hinv
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
     &                  (cff_w+cff1_w*hwater)*tl_hinv-                  &
#   ifdef TL_IOMS
     &                  cff1_w*hwater*hinv
#   endif
#  endif

!>            z_w(i,j,k)=Zt_avg1(i,j)+                                  &
!>   &                   (Zt_avg1(i,j)+hwater)*cff2_w
!>
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w(i,j,k)=tl_Zt_avg1(i,j)+                            &
     &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_w+         &
     &                      (Zt_avg1(i,j)+hwater)*tl_cff2_w-            &
#   ifdef TL_IOMS
     &                      (Zt_avg1(i,j)+hwater)*cff2_w
#   endif
#  else
              tl_z_w(i,j,k)=tl_Zt_avg1(i,j)+                            &
     &                      tl_Zt_avg1(i,j)*cff2_w+                     &
#   ifdef TL_IOMS
     &                      hwater*cff2_w
#   endif
#  endif

              cff2_r=(cff_r+cff1_r*hwater)*hinv
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_r=cff1_r*tl_hwater*hinv+                          &
     &                  (cff_r+cff1_r*hwater)*tl_hinv-                  &
#   ifdef TL_IOMS
     &                  hwater*hinv
#   endif
#  endif

!>            z_r(i,j,k)=Zt_avg1(i,j)+                                  &
!>   &                   (Zt_avg1(i,j)+hwater)*cff2_r
!>
#  if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_r(i,j,k)=tl_Zt_avg1(i,j)+                            &
     &                      (tl_Zt_avg1(i,j)+tl_hwater)*cff2_r+         &
     &                      (Zt_avg1(i,j)+hwater)*tl_cff2_r-            &
#   ifdef TL_IOMS
     &                      (Zt_avg1(i,j)+hwater)*cff2_r
#   endif
#  else
              tl_z_r(i,j,k)=tl_Zt_avg1(i,j)+                            &
     &                      tl_Zt_avg1(i,j)*cff2_r+                     &
#   ifdef TL_IOMS
     &                      hwater*cff2_r
#   endif
#  endif
#  ifdef ICESHELF
!>            z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
!>            z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
#  endif
!>            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
!>
              tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!>      CALL exchange_r2d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          h)
!>
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_h)
!>      CALL exchange_w3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
!>   &                          z_w)
!>
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tl_z_w)
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          z_r)
!>
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_z_r)
!>      CALL exchange_r3d_tile (ng, tile,                               &
!>   &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
!>   &                          Hz)
!>
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          tl_Hz)
      END IF

#  ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    h)
!>
      CALL mp_exchange2d (ng, tile, iRPM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_h)
!>    CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
!>   &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    z_w)
!>
      CALL mp_exchange3d (ng, tile, iRPM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_z_w)
!>    CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    z_r, Hz)
!>
      CALL mp_exchange3d (ng, tile, iRPM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_z_r, tl_Hz)
#  endif

      RETURN
      END SUBROUTINE rp_set_depth_tile

#  ifdef ADJUST_BOUNDARY
!
!***********************************************************************
      SUBROUTINE rp_set_depth_bry (ng, tile)
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
#   include "tile.h"
!
#   ifdef PROFILE
      CALL wclock_on (ng, iRPM, 12)
#   endif
      CALL rp_set_depth_bry_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % tl_h,                      &
#   ifdef ICESHELF
     &                            GRID(ng) % zice,                      &
#   endif
     &                            GRID(ng) % tl_Hz_bry)
#   ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
#   endif

      RETURN
      END SUBROUTINE rp_set_depth_bry
!
!***********************************************************************
      SUBROUTINE rp_set_depth_bry_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj, LBij, UBij, &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  h, tl_h,                        &
#   ifdef ICESHELF
     &                                  zice,                           &
#   endif
     &                                  tl_Hz_bry)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
#   ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#   ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: tl_h(LBi:,LBj:)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#    endif
      real(r8), intent(out) :: tl_Hz_bry(LBij:,:,:)
#   else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_h(LBi:UBi,LBj:UBj)
#    ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(out) :: tl_Hz_bry(LBij:UBij,N(ng),4)
#   endif
!
!  Local variable declarations.
!
      integer :: i, ibry, j, k

      real(r8) :: cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_w0
      real(r8) :: tl_cff2_w, tl_hinv, tl_hwater, tl_z_w0

      real(r8), dimension(0:N(ng)) :: tl_Zw

#   include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  This routine is special and equivalent to "tl_set_bry_tile".  The
!  additional terms that are added to RPM transformation are not needed
!  here because of the way that the adjustments to the boundary are
!  computed.
!-----------------------------------------------------------------------
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

        IF (tl_LBC(iwest,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w0=cff1_w*tl_hwate
#   endif
!>            Zw(k)=z_w0+BOUNDARY(ng)%zeta_west(j)*(1.0_r8+z_w0*hinv)
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=tl_z_w0+                                         &
     &                 BOUNDARY(ng)%tl_zeta_west(j)*                    &
     &                 (1.0_r8+z_w0*hinv)+                              &
     &                 BOUNDARY(ng)%zeta_west(j)*z_w0*tl_hinv+          &
     &                 BOUNDARY(ng)%zeta_west(j)*tl_z_w0*hinv
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_west(j)*(1.0_r8+z_w0*hinv)
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(j,k,iwest)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(j,k,iwest)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(ieast,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w0=cff1_w*tl_hwater
#   endif
!>            Zw(k)=z_w0+BOUNDARY(ng)%zeta_east(j)*(1.0_r8+z_w0*hinv)
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=tl_z_w0+                                         &
     &                 BOUNDARY(ng)%tl_zeta_east(j)*                    &
     &                 (1.0_r8+z_w0*hinv)+                              &
     &                 BOUNDARY(ng)%zeta_east(j)*z_w0*tl_hinv+          &
     &                 BOUNDARY(ng)%zeta_east(j)*tl_z_w0*hinv
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_east(j)*(1.0_r8+z_w0*hinv)
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(j,k,ieast)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(j,k,ieast)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(isouth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w0=cff1_w*tl_hwater
#   endif
!>            Zw(k)=z_w0+BOUNDARY(ng)%zeta_south(i)*(1.0_r8+z_w0*hinv)
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=tl_z_w0+                                         &
     &                 BOUNDARY(ng)%tl_zeta_south(i)*                   &
     &                 (1.0_r8+z_w0*hinv)+                              &
     &                 BOUNDARY(ng)%zeta_south(i)*z_w0*tl_hinv+         &
     &                 BOUNDARY(ng)%zeta_south(i)*tl_z_w0*hinv
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_south(i)*(1.0_r8+z_w0*hinv)
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(i,k,isouth)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(i,k,isouth)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(inorth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_z_w0=cff1_w*tl_hwater
#   endif
!>            Zw(k)=z_w0+BOUNDARY(ng)%zeta_north(i)*(1.0_r8+z_w0*hinv)
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=tl_z_w0+                                         &
     &                 BOUNDARY(ng)%tl_zeta_north(i)*                   &
     &                 (1.0_r8+z_w0*hinv)+                              &
     &                 BOUNDARY(ng)%zeta_north(i)*z_w0*tl_hinv+         &
     &                 BOUNDARY(ng)%zeta_north(i)*tl_z_w0*hinv
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_north(i)*(1.0_r8+z_w0*hinv)
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(i,k,inorth)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(i,k,inorth)=tl_Zw(k)-tl_Zw(k-1)
            END DO
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

        IF (tl_LBC(iwest,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/(hc(ng)+hwater)
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
     &                  (cff_w+cff1_w*hwater)*tl_hinv
#   endif
!>            Zw(k)=BOUNDARY(ng)%zeta_west(j)+
!>   &              (BOUNDARY(ng)%zeta_west(j)+hwater)*cff2_w
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_west(j)+                    &
     &                 (BOUNDARY(ng)%tl_zeta_west(j)+                   &
     &                  tl_hwater)*cff2_w+                              &
     &                 (BOUNDARY(ng)%zeta_west(j)+                      &
     &                  hwater)*tl_cff2_w
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_west(j)+                    &
     &                 BOUNDARY(ng)%tl_zeta_west(j)*cff2_w
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(j,k,iwest)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(j,k,iwest)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(ieast,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/(hc(ng)+hwater)
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
     &                  (cff_w+cff1_w*hwater)*tl_hinv
#   endif
!>            Zw(k)=BOUNDARY(ng)%zeta_east(j)+                          &
!>   &              (BOUNDARY(ng)%zeta_east(j)+hwater)*cff2_w
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_east(j)+                    &
     &                 (BOUNDARY(ng)%tl_zeta_east(j)+                   &
     &                  tl_hwater)*cff2_w+                              &
     &                 (BOUNDARY(ng)%zeta_east(j)+                      &
     &                  hwater)*tl_cff2_w
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_east(j)+                    &
     &                 BOUNDARY(ng)%tl_zeta_east(j)*cff2_w
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(j,k,ieast)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(j,k,ieast)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(isouth,isFsur,ng)%acquire.and.                       &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/(hc(ng)+hwater)
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
     &                  (cff_w+cff1_w*hwater)*tl_hinv
#   endif
!>            Zw(k)=BOUNDARY(ng)%zeta_south(i)+                         &
!>   &              (BOUNDARY(ng)%zeta_south(i)+hwater)*cff2_w
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_south(i)+                   &
     &                 (BOUNDARY(ng)%tl_zeta_south(i)+                  &
     &                  tl_hwater)*cff2_w+                              &
     &                 (BOUNDARY(ng)%zeta_south(i)+                     &
     &                  hwater)*tl_cff2_w
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_south(i)+                   &
     &                 BOUNDARY(ng)%tl_zeta_south(i)*cff2_w
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(i,k,isouth)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(i,k,isouth)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF

        IF (tl_LBC(inorth,isFsur,ng)%acquire.and.                        &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#   ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#   endif
            tl_hwater=tl_h(i,j)
            hinv=1.0_r8/(hc(ng)+hwater)
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
            tl_hinv=-hinv*hinv*tl_hwater
#   endif
!>          Zw(0)=-h(i,j)
!>
            tl_Zw(0)=-tl_h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_cff2_w=cff1_w*tl_hwater*hinv+                          &
     &                  (cff_w+cff1_w*hwater)*tl_hinv
#   endif
!>            Zw(k)=BOUNDARY(ng)%zeta_north(i)+                         &
!>   &              (BOUNDARY(ng)%zeta_north(i)+hwater)*cff2_w
!>
#   if defined SEDIMENT_NOT_YET && defined SED_MORPH_NOT_YET
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_north(i)+                   &
     &                 (BOUNDARY(ng)%tl_zeta_north(i)+                  &
     &                  tl_hwater)*cff2_w+                              &
     &                 (BOUNDARY(ng)%zeta_north(i)+                     &
     &                  hwater)*tl_cff2_w
#   else
              tl_Zw(k)=BOUNDARY(ng)%tl_zeta_north(i)+                   &
     &                 BOUNDARY(ng)%tl_zeta_north(i)*cff2_w
#   endif
#   ifdef ICESHELF
!>            Zw(k)=Zw(k)-ABS(zice(i,j))
#   endif
!>            Hz_bry(i,k,inorth)=Zw(k)-Zw(k-1)
!>
              tl_Hz_bry(i,k,inorth)=tl_Zw(k)-tl_Zw(k-1)
            END DO
          END DO
        END IF
      END IF

#   ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      DO ibry=1,4
!>      CALL mp_exchange3d_bry (ng, tile, iNLM, 1, ibry,                &
!>   &                          LBij, UBij, 1, N(ng),                   &
!>   &                          NghostPoints,                           &
!>   &                          EWperiodic(ng), NSperiodic(ng),         &
!>   &                          Hz_bry(:,:,ibry))
!>
        CALL mp_exchange3d_bry (ng, tile, iTLM, 1, ibry,                &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          tl_Hz_bry(:,:,ibry))
      END DO
#   endif

      RETURN
      END SUBROUTINE rp_set_depth_bry_tile
#  endif
# endif
!
!***********************************************************************
      SUBROUTINE rp_bath (ng, tile)
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
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 12)
# endif
      CALL rp_bath_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % tl_h)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 12)
# endif

      RETURN
      END SUBROUTINE rp_bath
!
!***********************************************************************
      SUBROUTINE rp_bath_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         h, tl_h)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
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
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(out) :: tl_h(LBi:,LBj:)
# else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_h(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize tangent linear bathymetry tl_h(i,j) to h(i,j) so some of
!  the terms are cancelled in the barotropic pressure gradient.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          tl_h(i,j)=h(i,j)
        END DO
      END DO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_h)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iRPM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_h)
# endif

      RETURN
      END SUBROUTINE rp_bath_tile
#endif
      END MODULE rp_set_depth_mod
