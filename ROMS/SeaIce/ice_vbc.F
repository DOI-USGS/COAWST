#include "cppdefs.h"
      MODULE ice_vbc_mod
#ifdef ICE_MODEL
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This module sets the ice-water and ice-air stresses for the
!  ice momentum equation.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC ice_vbc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_ice
      USE mod_coupling
#  ifdef LMD_SKPP
      USE mod_mixing
#  endif
      USE mod_stepping
!
      implicit none
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
# endif
      CALL ice_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   liold(ng), liuol(ng),                          &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   OCEAN(ng) % rho,                               &
     &                   COUPLING(ng) % Zt_avg1,                        &
     &                   ICE(ng) % ai,                                  &
     &                   ICE(ng) % hi,                                  &
     &                   ICE(ng) % ui,                                  &
     &                   ICE(ng) % vi,                                  &
     &                   ICE(ng) % tauaiu,                              &
     &                   ICE(ng) % tauaiv,                              &
     &                   ICE(ng) % uwater,                              &
     &                   ICE(ng) % vwater,                              &
     &                   ICE(ng) % sealev,                              &
# ifdef ICE_BULK_FLUXES
     &                   FORCES(ng) % sustr_aw,                         &
     &                   FORCES(ng) % svstr_aw,                         &
     &                   FORCES(ng) % tau_aix_n,                        &
     &                   FORCES(ng) % tau_aiy_n,                        &
# endif
# ifdef ICE_SHOREFAST
     &                   GRID(ng) % h,                                  &
# endif
     &                   ICE(ng) % utau_iw,                             &
     &                   ICE(ng) % chu_iw,                              &
     &                   ICE(ng) % spd_iw                               &
     &                   )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
# endif
      RETURN
      END SUBROUTINE ice_vbc
!
!***********************************************************************
      SUBROUTINE ice_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         liold, liuol,                            &
     &                         z_r, z_w,                                &
     &                         u, v,                                    &
     &                         sustr, svstr, rho,                       &
     &                         Zt_avg1,                                 &
     &                         ai, hi, ui, vi, tauaiu, tauaiv,          &
     &                         uwater, vwater, sealev,                  &
# ifdef ICE_BULK_FLUXES
     &                         sustr_aw, svstr_aw,                      &
     &                         tau_aix_n, tau_aiy_n,                    &
# endif
# ifdef ICE_SHOREFAST
     &                         h,                                       &
# endif
     &                         utau_iw, chu_iw, spd_iw)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      integer, intent(in) :: liold, liuol

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(out) :: tauaiu(LBi:,LBj:)
      real(r8), intent(out) :: tauaiv(LBi:,LBj:)
      real(r8), intent(in) :: uwater(LBi:,LBj:)
      real(r8), intent(in) :: vwater(LBi:,LBj:)
      real(r8), intent(out) :: sealev(LBi:,LBj:)
# ifdef ICE_BULK_FLUXES
      real(r8), intent(in) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(in) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(in) :: tau_aix_n(LBi:,LBj:)
      real(r8), intent(in) :: tau_aiy_n(LBi:,LBj:)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:,LBj:)
# endif
      real(r8), intent(out) :: utau_iw(LBi:,LBj:)
      real(r8), intent(inout) :: chu_iw(LBi:,LBj:)
      real(r8), intent(in) :: spd_iw(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: tauaiu(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tauaiv(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: uwater(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vwater(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: sealev(LBi:UBi,LBj:UBj)
# ifdef ICE_BULK_FLUXES
      real(r8), intent(in) :: sustr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tau_aix_n(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tau_aiy_n(LBi:UBi,LBj:UBj)
# endif
# ifdef ICE_SHOREFAST
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(out) :: utau_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: chu_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: spd_iw(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
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
# ifdef ICE_SHOREFAST
      real(r8) :: clear
      real(r8) :: masu
      real(r8) :: masv 
      real(r8) :: hu
      real(r8) :: hv
      real(r8) :: hh
# endif
      integer :: i, j

      integer :: k
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: spdiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: chuiw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: utauiw
# ifndef ICE_MK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uwind
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vwind
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wind_speed
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: windu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: windv
# endif
      real(r8) :: tauiwu
      real(r8) :: tauiwv
      real(r8) :: tauawu
      real(r8) :: tauawv
      real(r8) :: aix
      real(r8) :: aiy
      real(r8) :: spd
      real(r8) :: hix
      real(r8) :: hiy
      real(r8) :: chux
      real(r8) :: chuy
      real(r8) :: spdu
      real(r8) :: spdv

      real(r8) :: rhoO
      real(r8) :: dztop
      real(r8) :: thic
      real(r8) :: zdz0
      real(r8) :: z0

      real(r8), parameter :: kappa = 0.4_r8
      real(r8), parameter :: z0ii = 0.02_r8
      real(r8), parameter :: eps = 1.e-20

# include "set_bounds.h"

! *** Input from ocean model/data
# ifdef ICE_MK
      DO j=Jstr-1,Jend
        DO i=Istr-1,Iend
          rhoO = 1000._r8+rho(i,j,N(ng))
          spd = spd_iw(i,j)
#  ifdef ICE_SHOREFAST
          spd = max(spd,0.10_r8)
          clear = h(i,j)+Zt_avg1(i,j)-0.9*hi(i,j,liold)
          clear = MAX(clear,0.001_r8)
          IF(clear < 1.5_r8) THEN
            spd = (clear/1.5_r8)*spd
          ENDIF
#  else
          spd = max(spd,0.15_r8)
#  endif
#  ifdef ICE_MOM_BULK
          utauiw(i,j) = spd
          chuiw(i,j) = cdiw(ng)*spd
#  else
          thic = hi(i,j,liold)/max(ai(i,j,liold),max(min_a(ng),eps))
          z0 = max(z0ii*thic,0.01_r8)
          z0 = min(z0,0.1_r8)
          dztop=z_w(i,j,N(ng))-z_r(i,j,N(ng))
          zdz0 = dztop/z0
          if (zdz0 .lt. 3._r8) zdz0 = 3._r8
          utauiw(i,j) = sqrt(chu_iw(i,j)*spd)
          utauiw(i,j) = max(utauiw(i,j),1.E-04_r8)
          chuiw(i,j) = kappa*utauiw(i,j)/log(zdz0)
#   ifdef ICE_SHOREFAST
          hh = h(i,j)+Zt_avg1(i,j)
          clear = hh-0.9_r8*hi(i,j,liold)
          clear = MAX(clear,0.0_r8)
          IF(clear.lt.1.5_r8) chuiw(i,j) =                              &
     &       (MAX(clear-0.5_r8,0.0_r8)/1.0_r8)*chuiw(i,j)
#   endif
#  endif
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=IstrU,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          aix = 0.5_r8*(ai(i,j,liold)+ai(i-1,j,liold))
          hix = 0.5_r8*(hi(i,j,liold)+hi(i-1,j,liold))
          chux = 0.5_r8*(chuiw(i,j)+chuiw(i-1,j))
          tauaiu(i,j) = 0.5_r8*aix*(tau_aix_n(i,j)+tau_aix_n(i-1,j))    &
     &                        /rhoice(ng)
# ifdef ICE_BULK_FLUXES
          sustr(i,j) = aix*chux*(ui(i,j,liuol)-uwater(i,j))             &
     &                 + (1.0_r8-aix)*sustr_aw(i,j)
# endif
        END DO
      END DO

      DO j=JstrV,Jend
        DO i=Istr,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          aiy = 0.5_r8*(ai(i,j,liold)+ai(i,j-1,liold))
          hiy = 0.5_r8*(hi(i,j,liold)+hi(i,j-1,liold))
          chuy = 0.5_r8*(chuiw(i,j)+chuiw(i,j-1))
          tauaiv(i,j) = 0.5_r8*aiy*(tau_aiy_n(i,j)+tau_aiy_n(i,j-1))    &
     &                        /rhoice(ng)
# ifdef ICE_BULK_FLUXES
          svstr(i,j) = aiy*chuy*(vi(i,j,liuol)-vwater(i,j))             &
     &                 + (1.0_r8-aiy)*svstr_aw(i,j)
# endif
        END DO
      END DO
# else
      DO j=Jstr,Jend
        DO i=IstrU,Iend
           windu(i,j) = 0.5_r8*(Uwind(i-1,j)+Uwind(i,j))
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
           windv(i,j) = 0.5_r8*(Vwind(i,j-1)+Vwind(i,j))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          spd = spd_iw(i,j)
          spdiw(i,j) = (ai(i,j,liold)+0.02)*rho0*cdiw(ng)*spd           &
     &                  /rhoice(ng)
          spdiw(i,j) =                                                  &
     &       max(spdiw(i,j), (ai(i,j,liold)+0.02_r8)*cdiw(ng)           &
     &                                     *0.1_r8*rho0                 &
     &                        /rhoice(ng))
           wind_speed(i,j) = 0.5*sqrt((windu(i,j) +  windu(i+1,j))**2   &
     &                        +(windv(i,j) +  windv(i,j+1))**2)
# ifdef ICE_SHOREFAST
          clear = h(i,j)+Zt_avg1(i,j)-0.9*hi(i,j,liold)
          clear = MAX(clear,0.001_r8)
          IF(clear < 1.5_r8) THEN
            spdiw(i,j) = (clear/1.5_r8)*spdiw(i,j)
            wind(i,j) = (clear/1.5_r8)*wind(i,j)
          ENDIF
# endif
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=IstrU,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i-1,j,N(ng)))
          aix = 0.5_r8*(ai(i,j,liold)+ai(i-1,j,liold))
          hix = 0.5_r8*(hi(i,j,liold)+hi(i-1,j,liold))
          spd = 0.5_r8*(wind_speed(i,j)+wind_speed(i-1,j))
          tauaiu(i,j) = aix*rho_air(ng)*                                &
     &   (0.5_r8*cdai(ng)*(1.0_r8-COS(2.0_r8*pi*MIN(                    &
     &   (hix/(aix+0.02_r8)+0.1_r8),0.5_r8))))                          &
     &                *spd*windu(i,j)/rhoice(ng)
          tauiwu = 0.5_r8*(spdiw(i,j)+spdiw(i-1,j))                     &
     &               *(ui(i,j,liuol)-uwater(i,j))*rhoice(ng)/rhoO
# ifdef ICE_BULK_FLUXES
          sustr(i,j) = tauiwu + (1.0_r8-aix)*0.5_r8*                    &
     &                 (sustr_aw(i,j)+sustr_aw(i-1,j))/rhoO
# endif
        END DO
      END DO

      DO j=JstrV,Jend
        DO i=Istr,Iend
          rhoO = 1000._r8 + 0.5_r8*(rho(i,j,N(ng))+rho(i,j-1,N(ng)))
          aiy = 0.5_r8*(ai(i,j,liold)+ai(i,j-1,liold))
          hiy = 0.5_r8*(hi(i,j,liold)+hi(i,j-1,liold))
          spd = 0.5_r8*(wind_speed(i,j)+wind_speed(i,j-1))
          tauaiv(i,j) = aiy*rho_air(ng)*                                &
     &   (0.5_r8*cdai(ng)*(1.0_r8-COS(2.0_r8*pi*MIN(                    &
     &   (hiy/(aiy+0.02_r8)+0.1_r8),0.5_r8))))                          &
     &                *spd*windv(i,j)/rhoice(ng)
          tauiwv = 0.5_r8*(spdiw(i,j)+spdiw(i,j-1))                     &
     &               *(vi(i,j,liuol)-vwater(i,j))*rhoice(ng)/rho0
# ifdef ICE_BULK_FLUXES
          svstr(i,j) = tauiwv + (1.0_r8-aiy)*0.5_r8*                    &
     &                 (svstr_aw(i,j)+svstr_aw(i,j-1))/rhoO
# endif
        END DO
      END DO
# endif

      DO j=Jstr,Jend
        DO i=Istr,Iend
           sealev(i,j) = Zt_avg1(i,j)
           chu_iw(i,j) = chuiw(i,j)
           utau_iw(i,j) = utauiw(i,j)
        END DO
      END DO
!
!  Apply boundary conditions.
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sealev)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          utau_iw)
      CALL bc_r2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          chu_iw)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tauaiu)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tauaiv)
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sealev, utau_iw, chu_iw, sustr)
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    svstr, tauaiu, tauaiv)
#endif

      RETURN
      END SUBROUTINE ice_vbc_tile
#endif
      END MODULE ice_vbc_mod
