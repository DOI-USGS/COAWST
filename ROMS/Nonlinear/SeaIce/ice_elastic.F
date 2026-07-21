#include "cppdefs.h"
      MODULE ice_elastic_mod

#if defined ICE_MOMENTUM && defined ICE_EVP
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!=======================================================================
!                                                                      !
!  This routine timesteps the ice momentum equations.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef MICLM_NUDGING
      USE mod_clima
# endif
# if defined ICE_SHOREFAST || defined ICE_KEEL
      USE mod_coupling
# endif
# ifdef FASTICE_CLIMATOLOGY
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ice
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_u2d_tile, exchange_v2d_tile
      USE ice_uibc_mod,    ONLY : ice_uibc_tile
      USE ice_vibc_mod,    ONLY : ice_vibc_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_elastic
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_elastic (ng, tile, model)
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
      CALL ice_elastic_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       liold(ng), liuol(ng), liunw(ng),           &
     &                       lieol(ng), lienw(ng),                      &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef WET_DRY
     &                       GRID(ng) % rmask_wet,                      &
     &                       GRID(ng) % umask_wet,                      &
     &                       GRID(ng) % vmask_wet,                      &
# endif
# ifdef ICESHELF
     &                       GRID(ng) % zice,                           &
# endif
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_u,                           &
     &                       GRID(ng) % om_v,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % f,                              &
# if defined ICE_SHOREFAST || defined ICE_KEEL
     &                       GRID(ng) % h,                              &
     &                       COUPLING(ng) % Zt_avg1,                    &
# endif
# ifdef FASTICE_CLIMATOLOGY
     &                       FORCES(ng) % fastice_clm,                  &
# endif
# ifdef MICLM_NUDGING
     &                       CLIMA(ng) % uiclm,                         &
     &                       CLIMA(ng) % viclm,                         &
     &                       CLIMA(ng) % MInudgcof,                     &
# endif
     &                       ICE(ng) % Fi,                              &
     &                       ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_elastic
!
!***********************************************************************
      SUBROUTINE ice_elastic_tile (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             liold, liuol, liunw, lieol, lienw,   &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
# endif
# ifdef ICESHELF
     &                             zice,                                &
# endif
     &                             pm, pn, om_u, on_u, om_v, on_v,      &
     &                             f,                                   &
# if defined ICE_SHOREFAST || defined ICE_KEEL
     &                             h,                                   &
     &                             Zt_avg1,                             &
# endif
# ifdef FASTICE_CLIMATOLOGY
     &                             fastice_clm,                         &
# endif
# ifdef MICLM_NUDGING
     &                             uiclm, viclm, MInudgcof,             &
# endif
     &                             Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, liuol, liunw, lieol, lienw
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
# if defined ICE_SHOREFAST || defined ICE_KEEL
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
# endif
# ifdef FASTICE_CLIMATOLOGY
      real(r8), intent(in) :: fastice_clm(LBi:,LBj:)
# endif
# ifdef MICLM_NUDGING
      real(r8), intent(in) :: uiclm(LBi:,LBj:)
      real(r8), intent(in) :: viclm(LBi:,LBj:)
      real(r8), intent(in) :: MInudgcof(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Fi(LBi:,LBj:,:)
      real(r8), intent(inout) :: Si(LBi:,LBj:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
# if defined ICE_SHOREFAST || defined ICE_KEEL
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
# endif
# ifdef FASTICE_CLIMATOLOGY
      real(r8), intent(in) :: fastice_clm(LBi:UBi,LBj:UBj)
# endif
# ifdef MICLM_NUDGING
      real(r8), intent(in) :: uiclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: viclm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: MInudgcof(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Fi(LBi:UBi,LBj:UBj,nIceF)
      real(r8), intent(inout) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
#  endif
!
! Local variable definitions
!
      integer :: i, j
!
      real(r8) :: alfa, area, auf, avf, chux, chuy, clear, dsum, fakt
# ifdef ICE_KEEL
      real(r8) :: Cbu, Cbv, spd, thic, uvf, vuf
# endif
      real(r8) :: cff, f1, f2, s1, s2
# if defined ICE_SHOREFAST || defined ICE_KEEL
      real(r8) :: hu, hv
# endif
      real(r8) :: pmu, pnu, pmv, pnv
      real(r8) :: cosstang, sinstang
      real(r8) :: masu, mfu11, mfu21, mfu12, mfu22
      real(r8) :: masv, mfv11, mfv21, mfv12, mfv22
      real(r8) :: uforce, vforce

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Timestep ice U-momentum equations.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          uforce=0.0_r8
          area=om_u(i,j)*on_u(i,j)
          chux=0.5_r8*(Fi(i-1,j,icIOmt)+                                &
                       Fi(i  ,j,icIOmt))
          masu=0.5_r8*(Si(i  ,j,liold,isHice)+                          &
     &                 Si(i-1,j,liold,isHice))
          masu=MAX(masu, 0.1_r8)
          masu=masu*IceRho(ng)
          auf=MAX(0.5_r8*(Si(i-1,j,liold,isAice)+                       &
     &                    Si(i  ,j,liold,isAice)), 0.01_r8)
          pmu=0.5_r8*(pm(i,j)+pm(i-1,j))
          pnu=0.5_r8*(pn(i,j)+pn(i-1,j))
!
!  Forces from ice rheology x-direction.
!
          s1=(Si(i  ,j,lienw,isISxx)-                                   &
     &        Si(i-1,j,lienw,isISxx))*pmu
#ifdef MASKING
          IF ((umask(i,j  ).gt.0.0_r8).and.                             &
              (umask(i,j+1).lt.1.0_r8)) THEN
            f1=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# ifdef WET_DRY
          ELSE IF ((umask_wet(i,j  ).gt.0.0_r8).and.                    &
     &             (umask_wet(i,j+1).lt.1.0_r8)) THEN
            f1=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# endif
# ifdef ICESHELF
          ELSE IF ((zice(i-1,j).eq.0.0_r8).and.                         &
     &             (zice(i  ,j).eq.0.0_r8).and.                         &
     &             (zice(i-1,j+1)+zice(i,j+1)).ne.0.0_r8) THEN
            f1=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# endif
          ELSE
            f1=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                       &
     &                  Si(i  ,j+1,lienw,isISxy)+                       &
     &                  Si(i-1,j+1,lienw,isISxy)+                       &
     &                  Si(i-1,j  ,lienw,isISxy))
          END IF
#else
          f1=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                         &
     &                Si(i  ,j+1,lienw,isISxy)+                         &
     &                Si(i-1,j+1,lienw,isISxy)+                         &
     &                Si(i-1,j  ,lienw,isISxy))
#endif

#ifdef MASKING
          IF ((umask(i,j  ).gt.0.0_r8).and.                             &
     &        (umask(i,j-1).lt.1.0_r8)) THEN
            f2=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# ifdef WET_DRY
          ELSE IF ((umask_wet(i,j  ).gt.0.0_r8).and.                    &
     &             (umask_wet(i,j-1).lt.1.0_r8)) THEN
            f2=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# endif
# ifdef ICESHELF
          ELSE IF ((zice(i-1,j).eq.0.0_r8).and.                         &
     &             (zice(i  ,j).eq.0.0_r8).and.                         &
     &              zice(i-1,j-1)+zice(i,j-1).ne.0.0_r8) THEN
            f2=0.5_r8*(Si(i  ,j,lienw,isISxy)+                          &
     &                 Si(i-1,j,lienw,isISxy))
# endif
          ELSE
            f2=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                       &
     &                  Si(i-1,j  ,lienw,isISxy)+                       &
     &                  Si(i-1,j-1,lienw,isISxy)+                       &
     &                  Si(i  ,j-1,lienw,isISxy))
          END IF
#else
          f2=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                         &
     &                Si(i-1,j  ,lienw,isISxy)+                         &
     &                Si(i-1,j-1,lienw,isISxy)+                         &
     &                Si(i  ,j-1,lienw,isISxy))
#endif
          s2=(f1-f2)*pnu
          uforce=(s1+s2)*area
!
!  Add wind forcing.
!
          uforce=uforce+                                                &
     &           IceRho(ng)*Fi(i,j,icAIus)*area
!
!  Add pressure from tilting ocean surface.
!
          uforce=uforce-                                                &
     &           g*masu*on_u(i,j)*(Fi(i  ,j,icHsse)-                    &
     &                             Fi(i-1,j,icHsse))
!
!  Add Coriolis force.
!
          fakt=0.0_r8

#ifdef MASKING
# ifdef WET_DRY
          dsum=vmask(i-1,j  )*vmask_wet(i-1,j  )+                       &
     &         vmask(i  ,j  )*vmask_wet(i  ,j  )+                       &
     &         vmask(i-1,j+1)*vmask_wet(i-1,j+1)+                       &
     &         vmask(i  ,j+1)*vmask_wet(i  ,j+1)
          IF (dsum.gt.0.0_r8) fakt=1.0_r8/dsum

          mfv11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+            &
     &                  Si(i-1,j  ,liold,isHice)*f(i-1,j  ))*           &
     &          Si(i-1,j,lieol,isVevp)*                                 &
     &          vmask(i-1,j)*vmask_wet(i-1,j)
          mfv21=0.5_r8*(Si(i,j-1,liold,isHice)*f(i,j-1)+                &
     &                  Si(i,j  ,liold,isHice)*f(i,j  ))*               &
     &          Si(i,j,lieol,isVevp)*                                   &
     &          vmask(i,j)*vmask_wet(i,j)
          mfv12=0.5_r8*(Si(i-1,j  ,liold,isHice)*f(i-1,j  )+            &
     &                  Si(i-1,j+1,liold,isHice)*f(i-1,j+1))*           &
     &          Si(i-1,j+1,lieol,isVevp)*                               &
     &          vmask(i-1,j+1)*vmask_wet(i-1,j+1)
          mfv22=0.5_r8*(Si(i,j  ,liold,isHice)*f(i,j  )+                &
     &                  Si(i,j+1,liold,isHice)*f(i,j+1))*               &
     &          Si(i,j+1,lieol,isVevp)*                                 &
     &          vmask(i,j+1)*vmask_wet(i,j+1)
# else
          dsum=vmask(i-1,j  )+                                          &
     &         vmask(i  ,j  )+                                          &
     &         vmask(i-1,j+1)+                                          &
     &         vmask(i  ,j+1)
          IF (dsum.gt.0.0_r8) fakt=1.0_r8/dsum

          mfv11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+            &
     &                  Si(i-1,j  ,liold,isHice)*f(i-1,j  ))*           &
     &          Si(i-1,j,lieol,isVevp)*vmask(i-1,j)
          mfv21=0.5_r8*(Si(i,j-1,liold,isHice)*f(i,j-1)+                &
     &                  Si(i,j  ,liold,isHice)*f(i,j  ))*               &
     &          Si(i,j,lieol,isVevp)*vmask(i,j)
          mfv12=0.5_r8*(Si(i-1,j  ,liold,isHice)*f(i-1,j  )+            &
     &                  Si(i-1,j+1,liold,isHice)*f(i-1,j+1))*           &
     &          Si(i-1,j+1,lieol,isVevp)*vmask(i-1,j+1)
          mfv22=0.5_r8*(Si(i,j  ,liold,isHice)*f(i,j  )+                &
     &                  Si(i,j+1,liold,isHice)*f(i,j+1))*               &
     &          Si(i,j+1,lieol,isVevp)*vmask(i,j+1)
# endif
#else
           fakt=0.25_r8
           mfv11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+           &
     &                   Si(i-1,j  ,liold,isHice)*f(i-1,j  ))*          &
     &           Si(i-1,j,lieol,isVevp)
           mfv21=0.5_r8*(Si(i,j-1,liold,isHice)*f(i,j-1)+               &
     &                   Si(i,j  ,liold,isHice)*f(i,j  ))*              &
     &           Si(i,j,lieol,isVevp)
           mfv12=0.5_r8*(Si(i-1,j  ,liold,isHice)*f(i-1,j  )+           &
     &                   Si(i-1,j+1,liold,isHice)*f(i-1,j+1))*          &
     &           Si(i-1,j+1,lieol,isVevp)
           mfv22=0.5_r8*(Si(i,j  ,liold,isHice)*f(i,j  )+               &
     &                   Si(i,j+1,liold,isHice)*f(i,j+1))*              &
     &           Si(i,j+1,lieol,isVevp)
#endif
!
          uforce=uforce+                                                &
     &           fakt*IceRho(ng)*area*                                  &
     &           (mfv11+mfv21+mfv12+mfv22)
!
!  Add stress from ocean current.
!
          uforce=uforce/area +                                          &
     &           auf*chux*rho0*Fi(i,j,icUavg)
!
!  Compute basal stress as bottom drag to the deepest ice keels.
# ifdef ICE_KEEL
!  Use the ice keel basal stress parameterization of Lemieux (2005).
!
          hu=0.5_r8*(h(i-1,j)+MIN(Zt_avg1(i-1,j), 0.0_r8) +             &
     &               h(i  ,j)+MIN(Zt_avg1(i  ,j), 0.0_r8))
          thic=0.5_r8*(Si(i  ,j,liold,isHice)+                          &
     &                 Si(i-1,j,liold,isHice))/auf
          vuf=0.25*(Si(i  ,j  ,lieol,isVevp)+                           &
     &              Si(i  ,j+1,lieol,isVevp)+                           &
     &              Si(i-1,j  ,lieol,isVevp)+                           &
     &              Si(i-1,j+1,lieol,isVevp))
          spd=SQRT(Si(i,j,lieol,isUevp)*                                &
     &             Si(i,j,lieol,isUevp)+vuf*vuf)
!
!  Determine if the ice keel depth is greater than the water depth
!  (Lemieux k_1 parameter default value is 8.0)
!
          clear=MAX(thic-hu/8.0_r8, 0.0_r8)
!
!  Compute basal drag coefficient, C_bu. Use Lemieux k_2=20.0.
!
          Cbu=20.0_r8*clear*EXP(-20.0_r8*(1.0_r8-auf))/(spd+5.0E-5_r8)
!
!  For an implicit formulation, the basal stress just appears in the
!  RHS denominator.
!
          alfa=masu+dtevp(ng)*auf*rho0*chux +                           &
     &              dtevp(ng)*IceRho(ng)*Cbu
# else
          alfa=masu+dtevp(ng)*auf*rho0*chux
# endif
!
!  Solving the ice U-momentum equation.
!
          Si(i,j,lienw,isUevp)=(masu*Si(i,j,lieol,isUevp)+              &
    &                           dtevp(ng)*uforce)/alfa
!
# ifdef MICLM_NUDGING
          cff = 0.5*(MInudgcof(i,j)+MInudgcof(i-1,j))
          Si(i,j,lienw,isUevp)=Si(i,j,lienw,isUevp)+                    &
     &                         dtevp(ng)*cff*(uiclm(i,j)-               &
     &                                        Si(i,j,lienw,isUevp))
# endif
# ifdef ICE_SHOREFAST
          hu=0.5_r8*(h(i-1,j)+MIN(Zt_avg1(i-1,j), 0.0_r8)+              &
     &               h(i  ,j)+MIN(Zt_avg1(i  ,j), 0.0_r8))
          masu=0.5_r8*(Si(i  ,j,liold,isHice)+                          &
     &                 Si(i-1,j,liold,isHice))
          clear=hu-0.9_r8*masu
          clear=MAX(clear, 0.0_r8)
          IF (clear.lt.5.0_r8)                                          &
     &      Si(i,j,lienw,isUevp)=MAX(clear-1.0_r8, 0.0_r8)/             &
     &                           4.0_r8*Si(i,j,lienw,isUevp)
# endif
# ifdef FASTICE_CLIMATOLOGY
          Si(i,j,lienw,isUevp)=0.5_r8*(fastice_clm(i-1,j)+              &
     &                                 fastice_clm(i  ,j))*             &
     &                         Si(i,j,lienw,isUevp)
# endif
# ifdef MASKING
          Si(i,j,lienw,isUevp)=Si(i,j,lienw,isUevp)*umask(i,j)
# endif
# ifdef ICESHELF
          IF(zice(i-1,j)+zice(i,j).ne.0.0_r8) THEN
            Si(i,j,lienw,isUevp)=0.0_r8
          END IF
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Timestep ice V-momentum equations.
!-----------------------------------------------------------------------
!
      DO j=JstrP,Jend
        DO i=Istr,Iend
          vforce=0.0_r8
          area=om_v(i,j)*on_v(i,j)
          masv=0.5_r8*(Si(i,j  ,liold,isHice)+                          &
     &                 Si(i,j-1,liold,isHice))
          masv=MAX(masv, 0.1_r8)
          masv=masv*IceRho(ng)
          avf=MAX(0.5_r8*(Si(i,j-1,liold,isAice)+                       &
     &                    Si(i,j  ,liold,isAice)), 0.01_r8)
          chuy=0.5_r8*(Fi(i,j-1,icIOmt)+                                &
     &                 Fi(i,j  ,icIOmt))
          pmv=0.5_r8*(pm(i,j)+pm(i,j-1))
          pnv=0.5_r8*(pn(i,j)+pn(i,j-1))
!
!  Forces from ice rheology y-direction.
!
          s1=(Si(i,j  ,lienw,isISyy)-                                   &
     &        Si(i,j-1,lienw,isISyy))*pnv
#ifdef MASKING
          IF ((vmask(i  ,j).gt.0.0_r8).and.                             &
     &        (vmask(i+1,j).lt.1.0_r8)) THEN
            f1=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# ifdef WET_DRY
          ELSE IF ((vmask_wet(i  ,j).gt.0.0_r8).and.                    &
     &             (vmask_wet(i+1,j).lt.1.0_r8)) THEN
            f1=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# endif
# ifdef ICESHELF
          ELSE IF ((zice(i,j-1).eq.0.0_r8).and.                         &
     &             (zice(i,j  ).eq.0.0_r8).and.                         &
     &             (zice(i+1,j-1)+zice(i+1,j)).ne.0.0_r8) THEN
            f1=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# endif
          ELSE
            f1=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                       &
     &                  Si(i+1,j  ,lienw,isISxy)+                       &
     &                  Si(i+1,j-1,lienw,isISxy)+                       &
     &                  Si(i  ,j-1,lienw,isISxy))
          END IF
#else
          f1=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                         &
     &                Si(i+1,j  ,lienw,isISxy)+                         &
     &                Si(i+1,j-1,lienw,isISxy)+                         &
     &                Si(i  ,j-1,lienw,isISxy))
#endif

#ifdef MASKING
          IF ((vmask(i  ,j).gt.0.0_r8).and.                             &
     &        (vmask(i-1,j).lt.1.0_r8)) THEN
            f2=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# ifdef WET_DRY
          ELSE IF ((vmask_wet(i  ,j).gt.0.0_r8).and.                    &
     &             (vmask_wet(i-1,j).lt.1.0_r8)) THEN
            f2=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# endif
# ifdef ICESHELF
          ELSE IF ((zice(i,j-1).eq.0.0_r8).and.                         &
     &             (zice(i,j  ).eq.0.0_r8).and.                         &
     &             (zice(i-1,j-1)+zice(i-1,j)).ne.0.0_r8) THEN
            f2=0.5_r8*(Si(i,j  ,lienw,isISxy)+                          &
     &                 Si(i,j-1,lienw,isISxy))
# endif
          ELSE
            f2=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                       &
     &                  Si(i-1,j  ,lienw,isISxy)+                       &
     &                  Si(i-1,j-1,lienw,isISxy)+                       &
     &                  Si(i  ,j-1,lienw,isISxy))
          END IF
#else
          f2=0.25_r8*(Si(i  ,j  ,lienw,isISxy)+                         &
     &                Si(i-1,j  ,lienw,isISxy)+                         &
     &                Si(i-1,j-1,lienw,isISxy)+                         &
     &                Si(i  ,j-1,lienw,isISxy))
#endif
          s2=(f1-f2)*pmv
          vforce=(s1+s2)*area
!
!  Add wind forcing.
!
          vforce=vforce+                                                &
     &           IceRho(ng)*Fi(i,j,icAIvs)*area
!
!  Add pressure from tilting ocean surface.
!
          vforce=vforce-                                                &
     &           g*masv*om_v(i,j)*(Fi(i,j  ,icHsse)-                    &
     &                             Fi(i,j-1,icHsse))
!
!  Add Coriolis force.
!
          fakt=0.0_r8

#ifdef MASKING
# ifdef WET_DRY
          dsum=umask(i  ,j-1)*umask_wet(i  ,j-1)+                       &
     &         umask(i+1,j-1)*umask_wet(i+1,j-1)+                       &
     &         umask(i  ,j  )*umask_wet(i  ,j  )+                       &
     &         umask(i+1,j  )*umask_wet(i  ,j  )
          IF (dsum.gt.0.0_r8) fakt=1.0_r8/dsum

          mfu11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+            &
     &                  Si(i  ,j-1,liold,isHice)*f(i  ,j-1))*           &
     &          Si(i,j-1,lieol,isUevp)*                                 &
     &          umask(i,j-1)*umask_wet(i,j-1)
          mfu21=0.5_r8*(Si(i  ,j-1,liold,isHice)*f(i  ,j-1)+            &
     &                  Si(i+1,j-1,liold,isHice)*f(i+1,j-1))*           &
     &            Si(i+1,j-1,lieol,isUevp)*                             &
     &            umask(i+1,j-1)*umask_wet(i+1,j-1)
          mfu12=0.5_r8*(Si(i-1,j,liold,isHice)*f(i-1,j)+                &
     &                  Si(i  ,j,liold,isHice)*f(i  ,j))*               &
     &          Si(i,j,lieol,isUevp)*                                   &
     &          umask(i,j)*umask_wet(i,j)
          mfu22=0.5_r8*(Si(i  ,j,liold,isHice)*f(i  ,j)+                &
     &                  Si(i+1,j,liold,isHice)*f(i+1,j))*               &
     &          Si(i+1,j,lieol,isUevp)*                                 &
     &          umask(i+1,j)*umask_wet(i+1,j)
# else
          dsum=umask(i  ,j-1)+                                          &
     &         umask(i+1,j-1)+                                          &
     &         umask(i  ,j  )+                                          &
     &         umask(i+1,j  )
          IF (dsum.gt.0.0_r8) fakt=1.0_r8/dsum

          mfu11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+            &
     &                  Si(i  ,j-1,liold,isHice)*f(i  ,j-1))*           &
     &            Si(i,j-1,lieol,isUevp)*umask(i,j-1)
          mfu21=0.5_r8*(Si(i  ,j-1,liold,isHice)*f(i  ,j-1)+            &
     &                  Si(i+1,j-1,liold,isHice)*f(i+1,j-1))*           &
     &          Si(i+1,j-1,lieol,isUevp)*umask(i+1,j-1)
          mfu12=0.5_r8*(Si(i-1,j,liold,isHice)*f(i-1,j)+                &
     &                  Si(i  ,j,liold,isHice)*f(i  ,j))*               &
     &          Si(i,j,lieol,isUevp)*umask(i,j)
          mfu22=0.5_r8*(Si(i  ,j,liold,isHice)*f(i  ,j)+                &
     &                  Si(i+1,j,liold,isHice)*f(i+1,j))*               &
     &          Si(i+1,j,lieol,isUevp)*umask(i+1,j)
# endif
#else
          fakt=0.25_r8
          mfu11=0.5_r8*(Si(i-1,j-1,liold,isHice)*f(i-1,j-1)+            &
     &                  Si(i  ,j-1,liold,isHice)*f(i  ,j-1))*           &
     &          Si(i,j-1,lieol,isUevp)
          mfu21=0.5_r8*(Si(i  ,j-1,liold,isHice)*f(i  ,j-1)+            &
     &                  Si(i+1,j-1,liold,isHice)*f(i+1,j-1))*           &
     &          Si(i+1,j-1,lieol,isUevp)
          mfu12=0.5_r8*(Si(i-1,j,liold,isHice)*f(i-1,j)+                &
     &                  Si(i  ,j,liold,isHice)*f(i  ,j))*               &
     &          Si(i,j,lieol,isUevp)
          mfu22=0.5_r8*(Si(i  ,j,liold,isHice)*f(i  ,j)+                &
     &                  Si(i+1,j,liold,isHice)*f(i+1,j))*               &
     &          Si(i+1,j,lieol,isUevp)
#endif
!
          vforce=vforce-                                                &
     &           fakt*IceRho(ng)*area*                                  &
     &           (mfu11+mfu21+mfu12+mfu22)
!
!  Add stress from ocean current.
!
          vforce=vforce/area +                                          &
     &           avf*chuy*rho0*Fi(i,j,icVavg)
!
!  Compute basal stress as bottom drag to the deepest ice keels.
# ifdef ICE_KEEL
!  Use the ice keel basal stress parameterization of Lemieux (2005).
!
          hv=0.5_r8*(h(i,j-1)+MIN(Zt_avg1(i,j-1), 0.0_r8)+              &
     &               h(i,j  )+MIN(Zt_avg1(i,j  ), 0.0_r8))
          thic=0.5_r8*(Si(i,j  ,liold,isHice)+                          &
     &                 Si(i,j-1,liold,isHice))/avf
          uvf=0.25*(Si(i  ,j  ,lieol,isUevp)+                           &
     &              Si(i  ,j-1,lieol,isUevp)+                           &
     &              Si(i+1,j  ,lieol,isUevp)+
     &              Si(i+1,j-1,lieol,isUevp))
          spd=SQRT(Si(i,j,lieol,isVevp)*                                &
     &             Si(i,j,lieol,isVevp)+uvf*uvf)
!
!  Determine if the ice keel depth is greater than the water depth
!  (Lemieux k_1 parameter default value is 8.0)
!
          clear=MAX(thic-hv/8.0_r8, 0.0_r8)
!
!  Compute basal drag coefficient, C_bu. Use Lemieux k_2=20.0.
!
          Cbv=20.0_r8*clear*EXP(-20.0_r8*(1.0_r8-avf))/(spd+5.0E-5_r8)
!
!  For an implicit formulation, the basal stress just appears in the
!  RHS denominator.
!
          alfa=masv+dtevp(ng)*avf*rho0*chuy+                            &
     &              dtevp(ng)*IceRho(ng)*Cbv
# else
          alfa=masv+dtevp(ng)*avf*rho0*chuy
# endif
!
!  Solving the ice V-momentum equation.
!
          Si(i,j,lienw,isVevp)=(masv*Si(i,j,lieol,isVevp)+              &
     &                          dtevp(ng)*vforce)/alfa
!
# ifdef MICLM_NUDGING
          cff=0.5*(MInudgcof(i,j)+MInudgcof(i,j-1))
          Si(i,j,lienw,isVevp)=Si(i,j,lienw,isVevp)+                    &
     &                         dtevp(ng)*cff*(viclm(i,j)-               &
     &                                        Si(i,j,lienw,isVevp))
# endif
# ifdef ICE_SHOREFAST
          hv=0.5_r8*(h(i,j-1)+MIN(Zt_avg1(i,j-1), 0.0_r8) +             &
     &               h(i,j  )+MIN(Zt_avg1(i,j  ), 0.0_r8))
          masv=0.5_r8*(Si(i,j  ,liold,isHice)+                          &
     &                 Si(i,j-1,liold,isHice))
          clear=hv-0.9_r8*masv
          clear=MAX(clear, 0.0_r8)
          IF (clear.lt.5.0_r8)                                          &
     &       Si(i,j,lienw,isVevp)=MAX(clear-1.0_r8, 0.0_r8)/            &
     &                            4.0_r8*Si(i,j,lienw,isVevp)
# endif
# ifdef FASTICE_CLIMATOLOGY
          Si(i,j,lienw,isVevp)=0.5*(fastice_clm(i,j-1)+                 &
     &                              fastice_clm(i,j  ))*
     &                         Si(i,j,lienw,isVevp)
# endif
# ifdef MASKING
          Si(i,j,lienw,isVevp) = Si(i,j,lienw,isVevp)*vmask(i,j)
# endif
# ifdef ICESHELF
          IF (zice(i,j-1)+zice(i,j).ne.0.0_r8) THEN
            Si(i,j,lienw,isVevp)=0.0_r8
          END IF
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL ice_uibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    lieol, lienw,                                 &
     &                    Si(:,:,:,isUevp))

      CALL ice_vibc_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    lieol, lienw,                                 &
     &                    Si(:,:,:,isVevp))
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,lienw,isUevp))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,lienw,isVevp))
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,lienw,isUevp),                         &
     &                    Si(:,:,lienw,isVevp))
# endif

!
!-----------------------------------------------------------------------
!  Update ice velocity
!-----------------------------------------------------------------------
!
      IF(iEVP(ng).eq.nEVP(ng)) THEN
        DO j=Jstr,Jend
          DO i=IstrP,Iend
            Si(i,j,liunw,isUice)=Si(i,j,lienw,isUevp)
          END DO
        END DO
        DO j=JstrP,Jend
          DO i=Istr,Iend
            Si(i,j,liunw,isVice)=Si(i,j,lienw,isVevp)
          END DO
        END DO
!
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              Si(Istr,j,liunw,isUice)=Si(Istr,j,lienw,isUevp)
            END DO
            DO j=JstrP,Jend
              Si(Istr-1,j,liunw,isVice)=Si(Istr-1,j,lienw,isVevp)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              Si(Iend+1,j,liunw,isUice)=Si(Iend+1,j,lienw,isUevp)
            END DO
            DO j=JstrP,Jend
              Si(Iend+1,j,liunw,isVice)=Si(Iend+1,j,lienw,isVevp)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrP,Iend
              Si(i,Jstr-1,liunw,isUice)=Si(i,Jstr-1,lienw,isUevp)
            END DO
            DO i=Istr,Iend
              Si(i,Jstr,liunw,isVice)=Si(i,Jstr,lienw,isVevp)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrP,Iend
              Si(i,Jend+1,liunw,isUice)=Si(i,Jend+1,lienw,isUevp)
            END DO
            DO i=Istr,Iend
              Si(i,Jend+1,liunw,isVice)=Si(i,Jend+1,lienw,isVevp)
            END DO
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            Si(Istr  ,Jstr-1,liunw,isUice)=Si(Istr  ,Jstr-1,lienw,      &
     &                                        isUevp)
            Si(Istr-1,Jstr  ,liunw,isVice)=Si(Istr-1,Jstr  ,lienw,      &
     &                                        isVevp)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            Si(Iend+1,Jstr-1,liunw,isUice)=Si(Iend+1,Jstr-1,lienw,      &
     &                                        isUevp)
            Si(Iend+1,Jstr  ,liunw,isVice)=Si(Iend+1,Jstr  ,lienw,      &
     &                                        isVevp)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(iwest ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            Si(Istr  ,Jend+1,liunw,isUice)=Si(Istr  ,Jend+1,lienw,      &
     &                                        isUevp)
            Si(Istr-1,Jend+1,liunw,isVice)=Si(Istr-1,Jend+1,lienw,      &
     &                                        isVevp)
          END IF
        END IF
!
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng).or.        &
     &            CompositeGrid(ieast ,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            Si(Iend+1,Jend+1,liunw,isUice)=Si(Iend+1,Jend+1,lienw,      &
     &                                        isUevp)
            Si(Iend+1,Jend+1,liunw,isVice)=Si(Iend+1,Jend+1,lienw,      &
     &                                        isVevp)
          END IF
        END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,liunw,isUice))

        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Si(:,:,liunw,isVice))
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Si(:,:,liunw,isUice),                         &
     &                    Si(:,:,liunw,isVice))
# endif
      ENDIF
!
      RETURN
      END SUBROUTINE ice_elastic_tile
#endif
      END MODULE ice_elastic_mod
