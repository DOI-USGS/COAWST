#include "cppdefs.h"
      MODULE set_depth_mod
#ifdef SOLVE3D
!
!svn $Id: set_depth.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes the time evolving depths of the model grid    !
!  and its associated vertical transformation metric (thickness).      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_depth, set_depth_tile

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
      CALL set_depth_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     nstp(ng), nnew(ng),                          &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     GRID(ng) % bed_thick,                        &
# endif
     &                     COUPLING(ng) % Zt_avg1,                      &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
      RETURN
      END SUBROUTINE set_depth
!
!***********************************************************************
      SUBROUTINE set_depth_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
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
# if defined EW_PERIODIC || defined NS_PERIODIC
      USE exchange_2d_mod
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  if defined SEDIMENT && defined SED_MORPH
      real(r8), intent(in):: bed_thick(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: h(LBi:,LBj:)
      real(r8), intent(inout) :: Zt_avg1(LBi:,LBj:)
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
      real(r8), intent(inout) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
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
      integer :: i, j, k

      real(r8) :: cff_r, cff1_r, cff_w, cff1_w
      real(r8) :: hinv, hwater, z_r0, z_w0
# ifdef WET_DRY
      real(r8), parameter :: eps = 1.0E-14_r8
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute time evolving depths and vertical thicknesses.
!-----------------------------------------------------------------------
!
!  Compute vertical depths (meters, negative) at RHO- and W-points,
!  and vertical grid thicknesses.
!
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
          cff_w=hc(ng)*(SCALARS(ng)%sc_w(k)-SCALARS(ng)%Cs_w(k))
          cff1_w=SCALARS(ng)%Cs_w(k)
          cff_r=hc(ng)*(SCALARS(ng)%sc_r(k)-SCALARS(ng)%Cs_r(k))
          cff1_r=SCALARS(ng)%Cs_r(k)
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
# if defined EW_PERIODIC || defined NS_PERIODIC || defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
#  if defined EW_PERIODIC || defined NS_PERIODIC
      CALL exchange_r2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        h)
      CALL exchange_w3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        z_w)
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        z_r)
      CALL exchange_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Hz)
#  endif
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    h)
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    z_w)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    z_r, Hz)
#  endif
# endif
      RETURN
      END SUBROUTINE set_depth_tile
#endif
      END MODULE set_depth_mod
