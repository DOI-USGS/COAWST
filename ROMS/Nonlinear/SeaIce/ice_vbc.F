#include "cppdefs.h"
      MODULE ice_vbc_mod

#ifdef ICE_MODEL
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Paul Budgell       !
!    Licensed under a MIT/X style license           Katherine Hedstrom !
!    See License_ROMS.md                            Scott M. Durski    !
!=======================================================================
!                                                                      !
!  This routine sets the ice-water and ice-air stresses for the ice    !
!  momentum equations.                                                 !
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
      PUBLIC  :: ice_vbc
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_vbc (ng, tile, model)
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
      CALL ice_vbc_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   liold(ng), liuol(ng),                          &
# ifdef ICE_SHOREFAST
     &                   GRID(ng) % h,                                  &
# endif
# ifdef WET_DRY
     &                   GRID(ng) % umask_wet,                          &
     &                   GRID(ng) % vmask_wet,                          &
# endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
# ifdef ICESHELF
     &                   GRID(ng) % zice,                               &
# endif
     &                   COUPLING(ng) % Zt_avg1,                        &
     &                   OCEAN(ng) % rho,                               &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
# ifdef ICE_BULK_FLUXES
     &                   FORCES(ng) % sustr_ai,                         &
     &                   FORCES(ng) % svstr_ai,                         &
     &                   FORCES(ng) % sustr_ao,                         &
     &                   FORCES(ng) % svstr_ao,                         &
# endif
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   ICE(ng) % Fi,                                  &
     &                   ICE(ng) % Si)
# ifdef PROFILE
      CALL wclock_off (ng, model, 42, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_vbc
!
!***********************************************************************
      SUBROUTINE ice_vbc_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         liold, liuol,                            &
# ifdef ICE_SHOREFAST
     &                         h,                                       &
# endif
# ifdef WET_DRY
     &                         umask_wet, vmask_wet,                    &
# endif
     &                         z_r, z_w,                                &
# ifdef ICESHELF
     &                         zice,                                    &
# endif
     &                         Zt_avg1,                                 &
     &                         rho, u, v,                               &
# ifdef ICE_BULK_FLUXES
     &                         sustr_ai, svstr_ai,                      &
     &                         sustr_ao, svstr_ao,                      &
# endif
     &                         sustr, svstr,                            &
     &                         Fi, Si)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, liuol
!
# ifdef ASSUMED_SHAPE
#  ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef ICE_BULK_FLUXES
      real(r8), intent(in) :: sustr_ao(LBi:,LBj:)
      real(r8), intent(in) :: svstr_ao(LBi:,LBj:)
      real(r8), intent(in) :: sustr_ai(LBi:,LBj:)
      real(r8), intent(in) :: svstr_ai(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Si(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
      real(r8), intent(inout) :: Fi(LBi:,LBj:,:)
# else
#  ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef ICE_BULK_FLUXES
      real(r8), intent(in) :: sustr_ao(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr_ao(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sustr_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr_ai(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Si(LBi:UBi,LBj:UBj,2,nIceS)
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Fi(LBi:UBi,LBj:UBj,nIceF)
# endif
!
!  Local variable declarations.
!
      logical :: IceCavity
!
      integer :: i, j
!
      real(r8) :: aix, aiy, cff, chux, chuy, chuax, chuay, dztop
      real(r8) :: hix, hiy, rhoO, spd, thic, zdz0, z0
      real(r8) :: tauiwu, tauiwv
# ifdef ICE_SHOREFAST
      real(r8) :: clear, hh
# endif
!
      real(r8), parameter :: kappa = 0.4_r8
      real(r8), parameter :: z0ii = 0.01_r8
      real(r8), parameter :: eps = 1.0E-20
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: spdiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chuiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chuai
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utauiw
# ifndef ICE_MK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uwind
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vwind
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wind_speed
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: windu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: windv
# endif

# include "set_bounds.h"

# ifdef ICE_MK
!
!-----------------------------------------------------------------------
!  Compute ice-water and ice-air stresses according with the Mellor-
!  Kantha ice thermodynamics formulation.
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,Jend
        DO i=Istr-1,Iend
          rhoO=rho0+rho(i,j,N(ng))
          spd=Fi(i,j,icIOvs)
#  ifdef ICE_SHOREFAST
          spd=MAX(spd, 0.10_r8)
          clear=h(i,j)+                                                 &
     &          MIN(Zt_avg1(i,j), 0.0_r8)-0.9*Si(i,j,liold,isHice)
          clear=MAX(clear, 0.001_r8)
          IF (clear.lt.5.0_r8) THEN
            spd=0.2_r8*clear*spd                ! 0.2 = 1/5
          END IF
#  else
          spd=MAX(spd, 0.15_r8)
#  endif
#  ifdef ICE_MOM_BULK
          utauiw(i,j)=spd
          chuiw(i,j)=Cd_io(ng)*spd
#  else
          thic=Si(i,j,liold,isHice)/                                    &
     &         MAX(Si(i,j,liold,isAice), MAX(min_ai(ng), eps))
!
!  This parameterization for z0 seems to be designed to vary the
!  ice-ocean drag for ice between half a meter and 5 m thick (SMD).
!
          z0=MAX(z0ii*thic, 0.01_r8)
          z0=MIN(z0, 0.1_r8)
          dztop=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          zdz0=dztop/z0
!
!  This formulation for the ice-ocean momentum transfer coeficient,
!  chuiw, leads to extremely high values along the Alaskan coast where
!  the ice is particularly thick and the water is particularly shallow.
!  Therfore, the 'upper limit' is changed from 3 to 6 to cap the ice
!  water drag at a more 'reasonable' value in shallow water (SMD).
!
          IF (zdz0.lt.6.0_r8) zdz0=6.0_r8
          utauiw(i,j)=SQRT(Fi(i,j,icIOmt)*spd)
          utauiw(i,j)=MAX(utauiw(i,j), 1.0E-04_r8)
          chuiw(i,j)=kappa*utauiw(i,j)/LOG(zdz0)
!
!  Add an ice thickness dependent factor for the wind stress similar
!  to what is below for the othe formulation option (SMD).
!
          chuai(i,j)=1.0_r8-COS(1.0_r8*pi*MIN((thic+0.05_r8), 1.0_r8))
#   ifdef ICE_SHOREFAST
          hh=h(i,j)+MIN(Zt_avg1(i,j),0.0_r8)
          clear=hh-0.9_r8*Si(i,j,liold,isHice)
          clear=MAX(clear, 0.0_r8)
          IF (clear.lt.5.0_r8)                                          &
     &      chuiw(i,j)=(MAX(clear-1.0_r8, 0.0_r8)*0.25_r8)*chuiw(i,j)
#   endif
#  endif
        END DO
      END DO
!
!  Compute surface air-sea U-stress.
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          rhoO=1000.0_r8+0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          aix =0.5_r8*(Si(i,j,liold,isAice)+Si(i-1,j,liold,isAice))
          chux=0.5_r8*(chuiw(i,j)+chuiw(i-1,j))
!
!  Add thickness dependent wind ice drag (SMD).
!
#  ifdef ICE_MOM_BULK
          Fi(i,j,icAIus)=0.5_r8*aix*                                    &
     &                   (sustr_ai(i,j)+sustr_ai(i-1,j))/IceRho(ng)
#  else
          chuax=0.5_r8*(chuai(i,j)+chuai(i-1,j))
          Fi(i,j,icAIus)=0.5_r8*aix*chuax*                              &
     &                   (sustr_ai(i,j)+sustr_ai(i-1,j))/IceRho(ng)
#  endif
#  ifdef WET_DRY
          Fi(i,j,icAIus)=Fi(i,j,icAIus)*umask_wet(i,j)
#  endif
#  ifdef ICE_BULK_FLUXES
!
#   ifdef ICESHELF
          IceCavity=zice(i,j).ne.0.0_r8
#   else
          IceCavity=.FALSE.
#   endif
          IF (.not.IceCavity) THEN
            sustr(i,j)=aix*chux*(Si(i,j,liuol,isUice)-Fi(i,j,icUavg))+  &
     &                 (1.0_r8-aix)*sustr_ao(i,j)
#   ifdef WET_DRY
            sustr(i,j)=sustr(i,j)*umask_wet(i,j)
#   endif
          END IF
#  endif
        END DO
      END DO
!
!  Compute surface air-ice interface V-stress.
!
      DO j=JstrP,Jend
        DO i=Istr,Iend
          rhoO=1000.0_r8+0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          aiy =0.5_r8*(Si(i,j,liold,isAice)+Si(i,j-1,liold,isAice))
          chuy=0.5_r8*(chuiw(i,j)+chuiw(i,j-1))
!
!  Add thickness dependent wind ice drag (SMD).
!
#  ifdef ICE_MOM_BULK
          Fi(i,j,icAIvs)=0.5_r8*aiy*                                    &
     &                   (svstr_ai(i,j)+svstr_ai(i,j-1))/IceRho(ng)

#  else
          chuay=0.5_r8*(chuai(i,j)+chuai(i,j-1))
          Fi(i,j,icAIvs)=0.5_r8*aiy*chuay*                              &
     &                   (svstr_ai(i,j)+svstr_ai(i,j-1))/IceRho(ng)
#  endif
#  ifdef WET_DRY
          Fi(i,j,icAIvs)=Fi(i,j,icAIvs)*vmask_wet(i,j)
#  endif

#  ifdef ICE_BULK_FLUXES
!
#   ifdef ICESHELF
          IceCavity=zice(i,j).ne.0.0_r8
#   else
          IceCavity=.FALSE.
#   endif
          IF (.not.IceCavity) THEN
            svstr(i,j)=aiy*chuy*(Si(i,j,liuol,isVice)-Fi(i,j,icVavg))+          &
     &                 (1.0_r8-aiy)*svstr_ao(i,j)
#   ifdef WET_DRY
            svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#   endif
          END IF
#  endif
        END DO
      END DO
# else
!
!-----------------------------------------------------------------------
!  Alternate formulation to the Mellor-Kantha ice thermodynamics.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
           windu(i,j)=0.5_r8*(Uwind(i-1,j)+Uwind(i,j))
        END DO
      END DO
      DO j=JstrP,Jend
        DO i=Istr,Iend
           windv(i,j)=0.5_r8*(Vwind(i,j-1)+Vwind(i,j))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          cff=(Si(i,j,liold,isAice)+0.02)*rho0*Cd_io(ng)/IceRho(ng)
          spdiw(i,j)=cff*spd_iw(i,j)
          spdiw(i,j)=MAX(spdiw(i,j), cff*0.1_r8)
          wind_speed(i,j)=0.5*SQRT((windu(i,j)+windu(i+1,j))**2 +       &
     &                             (windv(i,j)+windv(i,j+1))**2)
#  ifdef ICE_SHOREFAST
          clear=h(i,j)+MIN(Zt_avg1(i,j), 0.0_r8)-0.9*hi(i,j,liold)
          clear=MAX(clear, 0.001_r8)
          IF (clear.lt.5.0_r8) THEN
            spd=0.2_r8*clear*spd
            spdiw(i,j)=0.2_r8*clear*spdiw(i,j)     ! 0.2 = 1/5
          END IF
#  endif
        END DO
      END DO
!
!  Compute surface air-ice interface U-stress.
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          rhoO=1000.0_r8+0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          aix=0.5_r8*(Si(i,j,liold,isAice)+Si(i-1,j,liold,isAice))
          hix=0.5_r8*(Si(i,j,liold,isHice)+Si(i-1,j,liold,isHice))
          spd=0.5_r8*(wind_speed(i,j)+wind_speed(i-1,j))
          Fi(i,j,icAIus)=aix*AirRho(ng)*                                &
     &                   (1.0_r8*Cd_ai(ng)*                             &
     &                    (1.0_r8-COS(1.0_r8*pi*                        &
     &                                MIN((hix/(aix+0.02_r8)+0.01_r8),  &
     &                                    1.0_r8))))*                   &
     &                   spd*windu(i,j)/IceRho(ng)

#  ifdef ICE_BULK_FLUXES
!
#   ifdef ICESHELF
          IceCavity=zice(i,j).ne.0.0_r8
#   else
          IceCavity=.FALSE.
#   endif
          IF (.not.IceCavity) THEN
            tauiwu=0.5_r8*(spdiw(i,j)+spdiw(i-1,j))*                    &
     &             (Si(i,j,liuol,IsUice)-Fi(i,j,icUavg))*IceRho(ng)/rhoO
            sustr(i,j)=tauiwu +                                         &
     &                 (1.0_r8-aix)*                                    &
     &                  0.5_r8*(sustr_ao(i,j)+sustr_ao(i-1,j))
#   ifdef WET_DRY
            sustr(i,j)=sustr(i,j)*umask_wet(i,j)
#   endif
          END IF
#  endif
        END DO
      END DO
!
!  Compute surface air-ice interface V-stress.
!
      DO j=JstrP,Jend
        DO i=Istr,Iend
          rhoO=1000.0_r8+0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          aiy=0.5_r8*(Si(i,j,liold,isAice)+Si(i,j-1,liold,isAice))
          hiy=0.5_r8*(Si(i,j,liold,isHice)+Si(i,j-1,liold,isHice))
          spd=0.5_r8*(wind_speed(i,j)+wind_speed(i,j-1))
          Fi(i,j,icAIvs)=aiy*AirRho(ng)*                                &
     &                   (0.5_r8*Cd_ai(ng)*                             &
     &                    (1.0_r8-COS(2.0_r8*pi*                        &
     &                                MIN((hiy/(aiy+0.02_r8)+0.1_r8),   &
     &                                    0.5_r8))))*                   &
     &                   spd*windv(i,j)/IceRho(ng)

#  ifdef ICE_BULK_FLUXES
!
#   ifdef ICESHELF
          IceCavity=zice(i,j).ne.0.0_r8
#   else
          IceCavity=.FALSE.
#   endif
          IF (.not.IceCavity) THEN
            tauiwv=0.5_r8*(spdiw(i,j)+spdiw(i,j-1))*                    &
     &             (Si(i,j,liuol,isVice)-Fi(i,j,icVavg))*IceRho(ng)/rho0
            svstr(i,j)=tauiwv+                                          &
     &                 (1.0_r8-aiy)*                                    &
     &                  0.5_r8*(svstr_ao(i,j)+svstr_ao(i,j-1))
#   ifdef WET_DRY
            svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
#   endif
          END IF
#  endif
        END DO
      END DO
# endif
!
# ifdef ICE_MK
!  Load sea surface elevation (isHsse), ice-ocean momentum transfer
!  coefficient (icIOmt), ice-ocean friction velocity (icIOfv).
# else
!  Load sea surface elevation (isHsse).
# endif
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Fi(i,j,icHsse)=Zt_avg1(i,j)
# ifdef ICE_MK
          Fi(i,j,icIOmt)=chuiw(i,j)
          Fi(i,j,icIOfv)=utauiw(i,j)
# endif
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icHsse))

      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icAIus))

      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icAIvs))

# ifdef ICE_MK
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icIOmt))

      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fi(:,:,icIOfv))
# endif

# ifdef ICE_BULK_FLUXES
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)

      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
# endif

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icHsse),                               &
     &                    Fi(:,:,icAIus),                               &
     &                    Fi(:,:,icAIvs))

#  ifdef ICE_MK
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    Fi(:,:,icIOmt),                               &
     &                    Fi(:,:,icIOfv))
#  endif

#  ifdef ICE_BULK_FLUXES
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sustr, svstr)
#  endif
# endif
!
      RETURN
      END SUBROUTINE ice_vbc_tile
!
#endif
      END MODULE ice_vbc_mod
