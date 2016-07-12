#include "cppdefs.h"
      MODULE set_depth_mod
#ifdef SOLVE3D
!
!svn $Id: set_depth.F 798 2016-05-21 17:04:58Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!  Currently, two vertical coordinate transformations are available    !
!  with various possible vertical stretching, C(s), functions, (see    !
!  routine "set_scoord.F" for details).                                !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_depth0, set_depth0_tile
      PUBLIC  :: set_depth,  set_depth_tile
# ifdef ADJUST_BOUNDARY
      PUBLIC  :: set_depth_bry
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_depth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_grid
      USE mod_ocean
# if defined SEDIMENT && defined SED_MORPH
      USE mod_sedbed
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
      CALL wclock_on (ng, iNLM, 12)
# endif
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp(ng), nnew(ng),                          &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     SEDBED(ng) % bed_thick,                      &
# endif
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
# endif

      RETURN
      END SUBROUTINE set_depth
!
!***********************************************************************
      SUBROUTINE set_depth_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp, nnew,                            &
     &                           h,                                     &
# ifdef ICESHELF
     &                           zice,                                  &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                           bed_thick,                             &
# endif
     &                           Zt_avg1,                               &
     &                           Hz, z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
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
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(out) :: Hz(LBi:,LBj:,:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0
      real(r8) :: C2_r, C2_w, hh2, vert_n1, vert_a, vert_h0, vert_s0
# ifdef WET_DRY
      real(r8), parameter :: eps = 1.0E-14_r8
# endif

# include "set_bounds.h"
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
      IF (Vtransform(ng).eq.3) THEN
        vert_s0 = 90./120.
        vert_n1 = 2./3. !3./5.
        vert_h0 = 60./(vert_s0)**(vert_n1)
        vert_a = 1./(vert_s0-1)
      END IF

      IF (Vtransform(ng).eq.1) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
# if defined SEDIMENT && defined SED_MORPH
            h(i,j)=h(i,j)-bed_thick(i,j,nstp)+bed_thick(i,j,nnew)
# endif
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
              h(i,j)=eps
            END IF
# endif
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
            END DO
          END DO
        END DO

      ELSE IF (Vtransform(ng).eq.3) THEN
        DO j=JstrR,JendR
          DO k=1,N(ng)
            DO i=IstrR,IendR
              hwater=h(i,j)
              hinv=1.0_r8/hwater

              hh2=(min(vert_h0,hwater))*hinv
              IF (SCALARS(ng)%sc_w(k).gt.-vert_s0) THEN 
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1) 
                C2_r=-hh2*(-SCALARS(ng)%sc_r(k))**(vert_n1)
  
                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-C2_r)
                cff1_r=C2_r
  
                z_w0=cff_w+cff1_w*hwater
                z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
                z_r0=cff_r+cff1_r*hwater
                z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)
              ELSE
                C2_w=-hh2*(-SCALARS(ng)%sc_w(k))**(vert_n1)             &
     &               -(1-hh2)*(vert_a*(SCALARS(ng)%sc_w(k)+vert_s0))**2
                C2_r=-hh2*(-SCALARS(ng)%sc_r(k))**(vert_n1)             &
     &               -(1-hh2)*(vert_a*(SCALARS(ng)%sc_r(k)+vert_s0))**2

                cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-C2_w)
                cff1_w=C2_w
                cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-C2_r)
                cff1_r=C2_r

                z_w0=cff_w+cff1_w*hwater
                z_w(i,j,k)=z_w0+Zt_avg1(i,j)*(1.0_r8+z_w0*hinv)
                z_r0=cff_r+cff1_r*hwater
                z_r(i,j,k)=z_r0+Zt_avg1(i,j)*(1.0_r8+z_r0*hinv)

              END IF
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
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
# if defined SEDIMENT && defined SED_MORPH
            h(i,j)=h(i,j)-bed_thick(i,j,nstp)+bed_thick(i,j,nnew)
# endif
# if defined WET_DRY
            IF (h(i,j).eq.0.0_r8) THEN
              h(i,j)=eps
            END IF
# endif
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
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

              z_w(i,j,k)=Zt_avg1(i,j)+(Zt_avg1(i,j)+hwater)*cff2_w
              z_r(i,j,k)=Zt_avg1(i,j)+(Zt_avg1(i,j)+hwater)*cff2_r
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
              Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
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
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h)
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          z_w)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          z_r)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Hz)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    h)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    z_w)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    z_r, Hz)
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Compute level ticknesses at the open boundaries using the provided
!  free-surface values (zeta_west, zeta_east, zeta_south, zeta_north).
!-----------------------------------------------------------------------
!
      CALL set_depth_bry (ng, tile)
# endif

      RETURN
      END SUBROUTINE set_depth_tile

!
!***********************************************************************
      SUBROUTINE set_depth0 (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
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
      CALL wclock_on (ng, iNLM, 12)
# endif
      CALL set_depth0_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % h,                               &
# ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
# endif
     &                      GRID(ng) % z0_r,                            &
     &                      GRID(ng) % z0_w)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
# endif

      RETURN
      END SUBROUTINE set_depth0
!
!***********************************************************************
      SUBROUTINE set_depth0_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            h,                                    &
# ifdef ICESHELF
     &                            zice,                                 &
# endif
     &                            z_r, z_w)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(out) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: z_w(LBi:,LBj:,0:)
# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k

      real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_r0, z_w0

      real(r8), parameter :: zeta0 = 0.0_r8

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Original formulation: Compute time independent vertical depths
!                        (meters, negative) at RHO- and W-points.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = Zo_w + zeta(x,y,0) * [1.0 + Zo_w / h(x,y)]
!
!                 Zo_w = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!  where zeta(x,y,0) = 0 for time independent depths.
!
!-----------------------------------------------------------------------
!
      IF (Vtransform(ng).eq.1) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
            cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
            cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
            cff1_r=SCALARS(ng)%Cs_r(k)
            cff1_w=SCALARS(ng)%Cs_w(k)
            DO i=IstrT,IendT
              hwater=h(i,j)
# ifdef ICESHELF
              hwater=hwater-ABS(zice(i,j))
# endif
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_w(i,j,k)=z_w0+zeta0*(1.0_r8+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r(i,j,k)=z_r0+zeta0*(1.0_r8+z_r0*hinv)
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  New formulation: Compute time independent vertical depths
!                   (meters, negative) at RHO- and W-points.
!  Various stretching functions are possible.
!
!         z_w(x,y,s,t) = zeta(x,y,0) + [zeta(x,y,t)+ h(x,y)] * Zo_w
!
!                 Zo_w = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
!  where zeta(x,y,0) = 0 for time independent depths.
!
!-----------------------------------------------------------------------
!
      ELSE IF (Vtransform(ng).eq.2) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            z_w(i,j,0)=-h(i,j)
          END DO
          DO k=1,N(ng)
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

              z_w(i,j,k)=zeta0+(zeta0+hwater)*cff2_w
              z_r(i,j,k)=zeta0+(zeta0+hwater)*cff2_r
# ifdef ICESHELF
              z_w(i,j,k)=z_w(i,j,k)-ABS(zice(i,j))
              z_r(i,j,k)=z_r(i,j,k)-ABS(zice(i,j))
# endif
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
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          z_w)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          z_r)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    z_w)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    z_r)
# endif

      RETURN
      END SUBROUTINE set_depth0_tile

# ifdef ADJUST_BOUNDARY
!
!***********************************************************************
      SUBROUTINE set_depth_bry (ng, tile)
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
      CALL wclock_on (ng, iNLM, 12)
#  endif
      CALL set_depth_bry_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         GRID(ng) % h,                            &
#  ifdef ICESHELF
     &                         GRID(ng) % zice,                         &
#  endif
     &                         GRID(ng) % Hz_bry)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 12)
#  endif

      RETURN
      END SUBROUTINE set_depth_bry
!
!***********************************************************************
      SUBROUTINE set_depth_bry_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj, LBij, UBij,    &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               h,                                 &
#  ifdef ICESHELF
     &                               zice,                              &
#  endif
     &                               Hz_bry)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange3d_bry
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
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: Hz_bry(LBij:,:,:)
#  else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: Hz_bry(LBij:UBij,N(ng),4)
#  endif
!
!  Local variable declarations.
!
      integer :: i, ibry, j, k

      real(r8) :: cff_w, cff1_w, cff2_w
      real(r8) :: hinv, hwater, z_w0

      real(r8), dimension(0:N(ng)) :: Zw

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

        IF (LBC(iwest,isFsur,ng)%acquire.and.                           &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
              Zw(k)=z_w0+BOUNDARY(ng)%zeta_west(j)*(1.0_r8+z_w0*hinv)
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(j,k,iwest)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(ieast,isFsur,ng)%acquire.and.                           &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
              Zw(k)=z_w0+BOUNDARY(ng)%zeta_east(j)*(1.0_r8+z_w0*hinv)
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(j,k,ieast)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(isouth,isFsur,ng)%acquire.and.                          &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
              Zw(k)=z_w0+BOUNDARY(ng)%zeta_south(i)*(1.0_r8+z_w0*hinv)
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(i,k,isouth)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(inorth,isFsur,ng)%acquire.and.                          &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/hwater
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
              cff1_w=SCALARS(ng)%Cs_w(k)
              z_w0=cff_w+cff1_w*hwater
              Zw(k)=z_w0+BOUNDARY(ng)%zeta_north(i)*(1.0_r8+z_w0*hinv)
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(i,k,inorth)=Zw(k)-Zw(k-1)
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

        IF (LBC(iwest,isFsur,ng)%acquire.and.                           &
     &      DOMAIN(ng)%Western_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(iwest,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              Zw(k)=BOUNDARY(ng)%zeta_west(j)+                          &
     &              (BOUNDARY(ng)%zeta_west(j)+hwater)*cff2_w
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(j,k,iwest)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(ieast,isFsur,ng)%acquire.and.                           &
     &      DOMAIN(ng)%Eastern_Edge(tile)) THEN
          i=BOUNDS(ng)%edge(ieast,r2dvar)
          DO j=JstrT,JendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              Zw(k)=BOUNDARY(ng)%zeta_east(j)+                          &
     &              (BOUNDARY(ng)%zeta_east(j)+hwater)*cff2_w
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(j,k,ieast)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(isouth,isFsur,ng)%acquire.and.                          &
     &      DOMAIN(ng)%Southern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(isouth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              Zw(k)=BOUNDARY(ng)%zeta_south(i)+                         &
     &              (BOUNDARY(ng)%zeta_south(i)+hwater)*cff2_w
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(i,k,isouth)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF

        IF (LBC(inorth,isFsur,ng)%acquire.and.                          &
     &      DOMAIN(ng)%Northern_Edge(tile)) THEN
          j=BOUNDS(ng)%edge(inorth,r2dvar)
          DO i=IstrT,IendT
            hwater=h(i,j)
#  ifdef ICESHELF
            hwater=hwater-ABS(zice(i,j))
#  endif
            hinv=1.0_r8/(hc(ng)+hwater)
            Zw(0)=-h(i,j)
            DO k=1,N(ng)
              cff_w=hc(ng)*SCALARS(ng)%sc_w(k)
              cff1_w=SCALARS(ng)%Cs_w(k)
              cff2_w=(cff_w+cff1_w*hwater)*hinv
              Zw(k)=BOUNDARY(ng)%zeta_north(i)+                         &
     &              (BOUNDARY(ng)%zeta_north(i)+hwater)*cff2_w
#  ifdef ICESHELF
              Zw(k)=Zw(k)-ABS(zice(i,j))
#  endif
              Hz_bry(i,k,inorth)=Zw(k)-Zw(k-1)
            END DO
          END DO
        END IF
      END IF

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      DO ibry=1,4
        CALL mp_exchange3d_bry (ng, tile, iNLM, 1, ibry,                &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          Hz_bry(:,:,ibry))
      END DO
#  endif

      RETURN
      END SUBROUTINE set_depth_bry_tile
# endif
#endif
      END MODULE set_depth_mod
