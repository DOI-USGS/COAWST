#include "cppdefs.h"
      MODULE ice_evp_sig_mod

#if defined ICE_MOMENTUM && defined ICE_EVP
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license           W. Paul Budgell    !
!    See License_ROMS.md                            Katherine Hedstrom !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine timesteps the EVP stresses term.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ice
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE ice_bc2d_mod,    ONLY : ice_bc2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_evp_sig
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_evp_sig (ng, tile, model)
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
      CALL ice_evp_sig_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       liold(ng), lieol(ng), lienw(ng),           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       ICE(ng) % Fi,                              &
     &                       ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_evp_sig
!
!***********************************************************************
      SUBROUTINE ice_evp_sig_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             liold, lieol, lienw,                 &
     &                             pm, pn,                              &
     &                             Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, lieol, lienw
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
# else
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
# endif
!
!  Local variable definitions.
!
      integer :: i, j
!
      real(r8) :: alfa, beta, cff, f1, f2, f3, f4, gamma
      real(r8) :: E, E0, ee, ees, ep, epx, epy
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_xx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_xy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: eps_yy

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Timestep internal ice stress tensor term.
!-----------------------------------------------------------------------
!
!  Initial value for Young modulus, E.  It has values between 0 and 1,
!  and makes the elastic term small compared to other terms.
!
      E0=0.25_r8
!
!  Comput stress tensor.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          eps_xx(i,j)=(Si(i+1,j,lieol,isUevp)-                          &
     &                 Si(i  ,j,lieol,isUevp))*pm(i,j)
          eps_yy(i,j)=(Si(i,j+1,lieol,isVevp)-                          &
     &                 Si(i,j  ,lieol,isVevp))*pn(i,j)
          epx=0.25_r8*(Si(i+1,j+1,lieol,isVevp)+                        &
     &                 Si(i+1,j  ,lieol,isVevp)-                        &
     &                 Si(i-1,j+1,lieol,isVevp)-                        &
     &                 Si(i-1,j  ,lieol,isVevp))*pm(i,j)
          epy=0.25_r8*(Si(i+1,j+1,lieol,isUevp)+                        &
     &                 Si(i  ,j+1,lieol,isUevp)-                        &
     &                 Si(i+1,j-1,lieol,isUevp)-                        &
     &                 Si(i  ,j-1,lieol,isUevp))*pn(i,j)
          eps_xy(i,j)=0.5_r8*(epx+epy)
        END DO
      END DO
!
!  Timestep internal ice stress tensor term. A split explicit scheme is
!  used with shorter steps (dts) to resolve the elastic wave velocity.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (Si(i,j,liold,isHice).gt.0.01_r8) THEN
            E  =2.0_r8*E0*IceRho(ng)*Si(i,j,liold,isHice)/              &
     &          (pm(i,j)*dtevp(ng))**2
            ep =E*Fi(i,j,icPice)/(4.0_r8*Fi(i,j,icBvis)+1.0E-8_r8)
            ee =E/(2.0_r8*Fi(i,j,icSvis)+1.0E-8_r8)
            ees=E*(Fi(i,j,icSvis)-Fi(i,j,icBvis))/                      &
     &          (4.0_r8*Fi(i,j,icSvis)*Fi(i,j,icBvis)+1.0E-8_r8)
!
            cff  =1.0_r8/dtevp(ng)
            alfa =cff+ee+ees
            beta =ees
            gamma=1.0_r8/dtevp(ng)+ee
            f1=E*eps_xx(i,j)+cff*Si(i,j,lieol,isISxx)-ep
            f2=E*eps_yy(i,j)+cff*Si(i,j,lieol,isISyy)-ep
            f3=E*eps_xy(i,j)+cff*Si(i,j,lieol,isISxy)
            f4=1.0_r8/(alfa**2-beta**2)
            Si(i,j,lienw,isISxx)=f4*(alfa*f1-beta*f2)
            Si(i,j,lienw,isISxy)=f3/gamma
            Si(i,j,lienw,isISyy)=f4*(alfa*f2-beta*f1)
          ELSE
            Si(i,j,lienw,isISxx)=2.0_r8*Fi(i,j,icSvis)*eps_xx(i,j)+     &
     &                           (Fi(i,j,icBvis)-Fi(i,j,icSvis))*       &
     &                           (eps_xx(i,j)+eps_yy(i,j))-             &
     &                           0.5_r8*Fi(i,j,icPice)
            Si(i,j,lienw,isISyy)=2.0_r8*Fi(i,j,icSvis)*eps_yy(i,j)+     &
     &                           (Fi(i,j,icBvis)-Fi(i,j,icSvis))*       &
     &                           (eps_xx(i,j)+eps_yy(i,j))-             &
     &                           0.5_r8*Fi(i,j,icPice)
            Si(i,j,lienw,isISxy)=2.0_r8*Fi(i,j,icSvis)*eps_xy(i,j)
          END IF
        END DO
      END DO
!
!  Apply lateral boundary conditions.
!
      CALL ice_bc2d_tile (ng, tile, model, isISxx,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    lieol, lienw,                                 &
     &                    Si(:,:,:,isUevp),                             &
     &                    Si(:,:,:,isVevp),                             &
     &                    Si(:,:,:,isISxx),                             &
     &                    LBC(:,ibICE(isISxx),ng))

      CALL ice_bc2d_tile (ng, tile, model, isISxy,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    lieol, lienw,                                 &
     &                    Si(:,:,:,isUevp),                             &
     &                    Si(:,:,:,isVevp),                             &
     &                    Si(:,:,:,isISxy),                             &
     &                    LBC(:,ibICE(isISxy),ng))

      CALL ice_bc2d_tile (ng, tile, model, isISyy,                      &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    lieol, lienw,                                 &
     &                    Si(:,:,:,isUevp),                             &
     &                    Si(:,:,:,isVevp),                             &
     &                    Si(:,:,:,isISyy),                             &
     &                    LBC(:,ibICE(isISyy),ng))
!
     IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,lienw,isISxx))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,lienw,isISxy))

        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,lienw,isISyy))
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,lienw,isISxx),                         &
     &                    Si(:,:,lienw,isISxy),                         &
     &                    Si(:,:,lienw,isISyy))
# endif
!
      RETURN
      END SUBROUTINE ice_evp_sig_tile
#endif
      END MODULE ice_evp_sig_mod
