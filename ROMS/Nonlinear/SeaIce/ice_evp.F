#include "cppdefs.h"
      MODULE ice_evp_mod

#if defined ICE_MOMENTUM && defined ICE_EVP
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module computes fields for elastic-viscous-plastic rheology.   !
!  It load elastic ice velocities (uie,vie) and couputes ice pressure  !
!  (icPice), bulk viscosity (icBvis), and shear viscosity (icSvis).    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_scalars
!
      USE bc_2d_mod,       ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_evp
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_evp (ng, tile, model)
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
      CALL ice_evp_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   liold(ng), lieol(ng), liuol(ng),               &
# ifdef OUTFLOW_MASK
     &                   GRID(ng) % mask_outflow,                       &
# endif
     &                   GRID(ng) % pm,                                 &
     &                   GRID(ng) % pn,                                 &
     &                   ICE(ng) % Fi,                                  &
     &                   ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_evp
!
!***********************************************************************
      SUBROUTINE ice_evp_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         liold, lieol, liuol,                     &
# ifdef OUTFLOW_MASK
     &                         mask_outflow,                            &
# endif
     &                         pm, pn,                                  &
     &                         Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol, liuol
!
# ifdef ASSUMED_SHAPE
#  ifdef OUTFLOW_MASK
      real(r8), intent(in) :: mask_outflow(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
# else
#  ifdef OUTFLOW_MASK
      real(r8), intent(in) :: mask_outflow(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
# endif
!
! Local variable definitions
!
      integer :: i, j
!
      real(r8) :: delta, e2r, eone, epx, epy, etwos, zmax

      real(r8), parameter :: epso = 1.0E-12_r8
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_xx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_yy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_xy

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute elastic terms as function of the stress tensor and rates of
!  strain.
!-----------------------------------------------------------------------
!
!  Load ice velocity components into EVP velocity arrays.
!
      IF (iEVP(ng).eq.1) THEN
        DO j=Jstr,Jend
          DO i=IstrP,Iend
            Si(i,j,1,isUevp)=Si(i,j,liuol,isUice)
            Si(i,j,2,isUevp)=Si(i,j,liuol,isUice)
          END DO
        END DO
        DO j=JstrP,Jend
          DO i=Istr,Iend
            Si(i,j,1,isVevp)=Si(i,j,liuol,isVice)
            Si(i,j,2,isVevp)=Si(i,j,liuol,isVice)
          END DO
        END DO
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              Si(Istr,j,1,isUevp)=Si(Istr,j,liuol,isUice)
              Si(Istr,j,2,isUevp)=Si(Istr,j,liuol,isUice)
            END DO
            DO j=JstrP,Jend
              Si(Istr-1,j,1,isVevp)=Si(Istr-1,j,liuol,isVice)
              Si(Istr-1,j,2,isVevp)=Si(Istr-1,j,liuol,isVice)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              Si(Iend+1,j,1,isUevp)=Si(Iend+1,j,liuol,isUice)
              Si(Iend+1,j,2,isUevp)=Si(Iend+1,j,liuol,isUice)
            END DO
            DO j=JstrP,Jend
              Si(Iend+1,j,1,isVevp)=Si(Iend+1,j,liuol,isVice)
              Si(Iend+1,j,2,isVevp)=Si(Iend+1,j,liuol,isVice)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrP,Iend
              Si(i,Jstr-1,1,isUevp)=Si(i,Jstr-1,liuol,isUice)
              Si(i,Jstr-1,2,isUevp)=Si(i,Jstr-1,liuol,isUice)
            END DO
            DO i=Istr,Iend
              Si(i,Jstr,1,isVevp)=Si(i,Jstr,liuol,isVice)
              Si(i,Jstr,2,isVevp)=Si(i,Jstr,liuol,isVice)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrP,Iend
              Si(i,Jend+1,1,isUevp)=Si(i,Jend+1,liuol,isUice)
              Si(i,Jend+1,2,isUevp)=Si(i,Jend+1,liuol,isUice)
            END DO
            DO i=Istr,Iend
              Si(i,Jend+1,1,isVevp)=Si(i,Jend+1,liuol,isVice)
              Si(i,Jend+1,2,isVevp)=Si(i,Jend+1,liuol,isVice)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            Si(Istr  ,Jstr-1,1,isUevp)=Si(Istr  ,Jstr-1,liuol,isUice)
            Si(Istr  ,Jstr-1,2,isUevp)=Si(Istr  ,Jstr-1,liuol,isUice)
            Si(Istr-1,Jstr  ,1,isVevp)=Si(Istr-1,Jstr  ,liuol,isVice)
            Si(Istr-1,Jstr  ,2,isVevp)=Si(Istr-1,Jstr  ,liuol,isVice)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            Si(Iend+1,Jstr-1,1,isUevp)=Si(Iend+1,Jstr-1,liuol,isUice)
            Si(Iend+1,Jstr-1,2,isUevp)=Si(Iend+1,Jstr-1,liuol,isUice)
            Si(Iend+1,Jstr  ,1,isVevp)=Si(Iend+1,Jstr  ,liuol,isVice)
            Si(Iend+1,Jstr  ,2,isVevp)=Si(Iend+1,Jstr  ,liuol,isVice)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            Si(Istr  ,Jend+1,1,isUevp)=Si(Istr  ,Jend+1,liuol,isUice)
            Si(Istr  ,Jend+1,2,isUevp)=Si(Istr  ,Jend+1,liuol,isUice)
            Si(Istr-1,Jend+1,1,isVevp)=Si(Istr-1,Jend+1,liuol,isVice)
            Si(Istr-1,Jend+1,2,isVevp)=Si(Istr-1,Jend+1,liuol,isVice)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            Si(Iend+1,Jend+1,1,isUevp)=Si(Iend+1,Jend+1,liuol,isUice)
            Si(Iend+1,Jend+1,2,isUevp)=Si(Iend+1,Jend+1,liuol,isUice)
            Si(Iend+1,Jend+1,1,isVevp)=Si(Iend+1,Jend+1,liuol,isVice)
            Si(Iend+1,Jend+1,2,isVevp)=Si(Iend+1,Jend+1,liuol,isVice)
          END IF
        END IF
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,1,isUevp),                             &
     &                    Si(:,:,2,isUevp),                             &
     &                    Si(:,:,1,isVevp),                             &
     &                    Si(:,:,2,isVevp))
# endif
!
!  Compute strain rates. Then, compute ice pressure (strength; icPice),
!  bulk viscosity (icBvis), and shear viscosity (icSvis).
!
      e2r = 1.0_r8/(ellip_sq(ng))
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          eps_xx(i,j)=(Si(i+1,j,lieol,isUevp)-                          &
     &                Si(i  ,j,lieol,isUevp))*pm(i,j)
          eps_yy(i,j)=(Si(i,j+1,lieol,isVevp)-                          &
     &                Si(i,j  ,lieol,isVevp))*pn(i,j)
          epx=0.25_r8*(Si(i+1,j+1,lieol,isVevp)+                        &
     &                 Si(i+1,j  ,lieol,isVevp)-                        &
     &                 Si(i-1,j+1,lieol,isVevp)-                        &
     &                 Si(i-1,j  ,lieol,isVevp))*pm(i,j)
          epy=0.25_r8*(Si(i+1,j+1,lieol,isUevp)+                        &
     &                 Si(i  ,j+1,lieol,isUevp)-                        &
     &                 Si(i+1,j-1,lieol,isUevp)-                        &
     &                 Si(i  ,j-1,lieol,isUevp))*pn(i,j)
          eps_xy(i,j)=0.5_r8*(epx+epy)
!
          eone=eps_xx(i,j)+eps_yy(i,j)
          etwos=(eps_xx(i,j)-eps_yy(i,j))*                              &
                (eps_xx(i,j)-eps_yy(i,j))+                              &
     &          4.0_r8*eps_xy(i,j)*eps_xy(i,j)
!
          delta=ABS(eone**2+e2r*etwos)
          delta=MAX(SQRT(delta),epso)
# ifdef ICE_STRENGTH_QUAD
          Fi(i,j,icPice)=Fi(i,j,icPgrd)*                                &
                         Si(i,j,liold,isHice)*Si(i,j,liold,isHice)*     &
     &                   EXP(-Astrength(ng)*                            &
     &                       (1.0_r8-Si(i,j,liold,isAice)))
# else
          Fi(i,j,icPice)=Pstar(ng)*Si(i,j,liold,isHice)*                &
     &                   EXP(-Astrength(ng)*                            &
     &                       (1.0_r8-Si(i,j,liold,isAice)))
# endif
          Fi(i,j,icBvis)=Fi(i,j,icPice)/(2.0_r8*delta)
          zmax=2.5E+8_r8*Fi(i,j,icPice)
          Fi(i,j,icBvis)=MIN(Fi(i,j,icBvis), zetaMax(ng))
          Fi(i,j,icBvis)=MAX(Fi(i,j,icBvis), zetaMin(ng))
          Fi(i,j,icSvis)=e2r*Fi(i,j,icBvis)
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icPice))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icBvis))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icSvis))

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icPice),                               &
     &                    Fi(:,:,icBvis),                               &
     &                    Fi(:,:,icSvis))
# endif

# ifdef OUTFLOW_MASK
!
!  Apply outflow mask.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Fi(i,j,icPice)=Fi(i,j,icPice)*mask_outflow(i,j)
          Fi(i,j,icBvis)=Fi(i,j,icBvis)*mask_outflow(i,j)
          Fi(i,j,icSvis)=Fi(i,j,icSvis)*mask_outflow(i,j)
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE ice_evp_tile
#endif
      END MODULE ice_evp_mod
