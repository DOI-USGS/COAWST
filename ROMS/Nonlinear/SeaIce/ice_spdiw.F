#include "cppdefs.h"
      MODULE ice_spdiw_mod

#ifdef ICE_MODEL
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module computes the magnitude of the shear between the ice     !
!  and the surface water. In this case, the surface water is defined   !
!  as the water in a surface mixed layer, so that velocity must be     !
!  computed first.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ice
# ifdef LMD_SKPP
      USE mod_mixing
# endif
      USE mod_ocean
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_spdiw
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_spdiw (ng, tile, model)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 42, __LINE__, MyFile)
# endif
      CALL ice_spdiw_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef ICE_MODEL
     &                     liuol(ng),                                   &
# endif
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
# ifdef LMD_SKPP
     &                     MIXING(ng) % hsbl,                           &
# endif
     &                     ICE(ng) % Fi,                                &
     &                     ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_spdiw
!
!***********************************************************************
      SUBROUTINE ice_spdiw_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef ICE_MODEL
     &                           liuol,                                 &
# endif
     &                           z_r, z_w,                              &
     &                           u, v,                                  &
# ifdef LMD_SKPP
     &                           hsbl,                                  &
# endif
     &                           Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liuol
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: nlio, nbotu, nbotv, k
!
      integer,  dimension(IminS:ImaxS,JminS:JmaxS) :: nbot
!
      real(r8) :: dml, mlio, totml
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vw

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute magnitude of the shear between the ice and the surface water.
!-----------------------------------------------------------------------
!
!  Determine model level associated with the depth of the mixed layer.
!
      DO j=Jstrm2,Jendp2
        DO i=Istrm2,Iendp2
# ifdef LMD_SKPP
          mlio=MIN(-ABS(Hsbl(i,j)),-10.0_r8)
# else
          mlio=-10.0_r8
# endif
          nbot(i,j)=1
          DO k=N(ng),1,-1
            IF (z_r(i,j,k).lt.mlio) THEN
              nbot(i,j)=MIN(k,N(ng))
              nbot(i,j)=MAX(nbot(i,j),1)
              EXIT
            END IF
          END DO
        END DO
      END DO
!
!  Compute verticaly averaged U-velocity over mixed-layer thickness.
!
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend+1
          nlio=0
          nbotu=NINT(0.5_r8*(nbot(i-1,j)+nbot(i,j)))
          nbotu=MAX(MIN(nbotu,N(ng)),1)
          uw(i,j)=0.0_r8
          totml=0.0_r8
          DO k=N(ng),nbotu,-1
            nlio=nlio+1
            dml=0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+                    &
     &                  z_w(i  ,j,k)-z_w(i  ,j,k-1))
            uw(i,j)=uw(i,j)+u(i,j,k,nrhs)*dml
            totml=totml+dml
          END DO
          uw(i,j)=uw(i,j)/totml
        END DO
      END DO
!
!  Compute verticaly averaged U-velocity over mixed-layer thickness.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr,Iend
          nlio=0
          nbotv=NINT(0.5_r8*(nbot(i,j-1)+nbot(i,j)))
          nbotv=MAX(MIN(nbotv,N(ng)),1)
          vw(i,j)=0.0_r8
          totml=0.0_r8
          DO k=N(ng),nbotv,-1
            nlio=nlio+1
            dml=0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+                    &
     &                  z_w(i,j  ,k)-z_w(i,j  ,k-1))
            vw(i,j)=vw(i,j)+v(i,j,k,nrhs)*dml
            totml=totml+dml
          END DO
          vw(i,j)=vw(i,j)/totml
        END DO
      END DO
!
!  Compute magnitude of the shear between ice velocities and vertically
!  averaged mixed layer velocities.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Fi(i,j,icIOvs)=0.5_r8*                                        &
     &                   SQRT((uw(i  ,j)-Si(i  ,j,liuol,isUice)+        &
     &                         uw(i+1,j)-Si(i+1,j,liuol,isUice))**2+    &
     &                        (vw(i,j  )-Si(i,j  ,liuol,isVice)+        &
     &                         vw(i,j+1)-Si(i,j+1,liuol,isVice))**2)
        END DO
      END DO
!
!  Load vertically averaged mixed layer velocity components.
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          Fi(i,j,icUavg)=uw(i,j)
        END DO
      END DO
      DO j=JstrP,Jend
        DO i=Istr,Iend
          Fi(i,j,icVavg)=vw(i,j)
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fi(:,:,icIOvs))

        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fi(:,:,icUavg))

        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fi(:,:,icVavg))

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icIOvs),                               &
     &                    Fi(:,:,icUavg),                               &
     &                    Fi(:,:,icVavg))
# endif
!
      RETURN
      END SUBROUTINE ice_spdiw_tile
#endif
      END MODULE ice_spdiw_mod
